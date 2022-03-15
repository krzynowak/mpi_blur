

#include "MPIDataProcessor.h"
#include "Alloc.h"
#include <stdio.h>
#include <math.h> 
#include <iostream>
#include <iomanip>
#include <limits>
#include "mpi.h"

#define SET_ROW_CNT(x,y,z) ( z < ((x) % (y)) ) ? ( ((x) / (y)) + 1 ) : ((x) / (y))
#define ROOT_CORE 0
#define NO_TAG 0

double ** getArray2D(int row, int col)
{
	double **result;
	result = new double* [ row ]; // size rows
	for ( int i = 0; i < row; i++ )
		result[ i ] = new double[ col ]; // size cols
	return result;
}





void MPIDataProcessor::sendData()
{
	if(this->core_id)
	{

		for(int i=this->margin;i<this->assigned_rows + this->margin ;i++)
		{
			MPI_Recv(
			this->data[i],
			this->dataSize,
			MPI_DOUBLE,
			ROOT_CORE,
			NO_TAG,
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		}

		/* last core must receive it's margin now (0 is taken care of and the rest will have a routine called every iterration) */
		if(this->core_id == (this->core_count - 1))
		{
			for(int i = this->assigned_rows + this->margin; i < this->assigned_rows + this->margin * 2; i++ )
			{
				MPI_Recv(
				this->data[i],
				this->dataSize,
				MPI_DOUBLE,
				ROOT_CORE,
				NO_TAG,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
			}

		}
	}
	else
	{

		int lastRow = SET_ROW_CNT(this->dataSize - 2 * this->margin, this->core_count, ROOT_CORE);
		lastRow += this->margin;
		int destination_core_count = lastRow;

		//send data for each core
		for (int i=1;i<this->core_count;i++)
		{
			destination_core_count += SET_ROW_CNT(this->dataSize - 2 * this->margin, this->core_count, i);

			//send data for each row 
			for(int j=lastRow;j<destination_core_count;j++)
			{
				//send
				MPI_Send(
				this->data[j],
				this->dataSize,
				MPI_DOUBLE,
				i,
				NO_TAG,
				MPI_COMM_WORLD);
			}

			lastRow = destination_core_count;
		}

		for(int i = 0; i < this->margin; i++)
		{
			//send
			MPI_Send(
			this->data[lastRow + i],
			this->dataSize,
			MPI_DOUBLE,
			this->core_count - 1,
			NO_TAG,
			MPI_COMM_WORLD);
		}
	}
}

void MPIDataProcessor::syncMargins()
{
	int dst = (this->core_id+1)% this->core_count;
	int src = (this->core_id - 1 + this->core_count) % this->core_count;


	/****************************  
		FIRST PASS -> FORWARD
	****************************/

	/* Last and first core only perform one of the specified operrations */
	if( (this->core_id != ROOT_CORE) && (this->core_id != (this->core_count - 1)) )
	{
		for(int i=0; i<this->margin; i++)
		{
			MPI_Sendrecv(
			//send
			this->data[this->assigned_rows + i],
			this->dataSize,
			MPI_DOUBLE,
			dst,
			NO_TAG,
			//rec
			this->data[i],
			this->dataSize,
			MPI_DOUBLE,
			src,
			NO_TAG,
			//com
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE );
		}
	}
	else
	{
		if( this->core_id )
		{
			for(int i=0; i<this->margin; i++)
			{
				/* Last core only gets data */
				MPI_Recv(
				this->data[i],
				this->dataSize,
				MPI_DOUBLE,
				this->core_id - 1,
				NO_TAG,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
			}
		}
		else
		{
			for(int i=0; i<this->margin; i++)
			{	
				/* First core only sends data */
				MPI_Send(
				this->data[this->assigned_rows + i],
				this->dataSize,
				MPI_DOUBLE,
				this->core_id + 1,
				NO_TAG,
				MPI_COMM_WORLD);
			}
		}
	}

	/*****************************
		SECOND PASS -> BACKWARD
	*****************************/
	/* Last and first core only perform one of the specified operrations */
	if( (this->core_id != ROOT_CORE) && (this->core_id != (this->core_count - 1)) )
	{
		for(int i=0; i<this->margin; i++)
		{
			MPI_Sendrecv(
			//send
			this->data[this->margin + i],
			this->dataSize,
			MPI_DOUBLE,
			src,
			NO_TAG,
			//rec
			this->data[this->assigned_rows + this->margin + i],
			this->dataSize,
			MPI_DOUBLE,
			dst,
			NO_TAG,
			//com
			MPI_COMM_WORLD,
			MPI_STATUS_IGNORE );
		}
	}
	else
	{
		if( this->core_id )
		{
			for(int i=0; i<this->margin; i++)
			{	
				/* Last core only sends data */
				MPI_Send(
				this->data[this->margin + i],
				this->dataSize,
				MPI_DOUBLE,
				this->core_id - 1,
				NO_TAG,
				MPI_COMM_WORLD);
			}
		}
		else
		{
			
			for(int i=0; i<this->margin; i++)
			{
				/* First core only gets data */
				MPI_Recv(
				this->data[this->assigned_rows + this->margin + i],
				this->dataSize,
				MPI_DOUBLE,
				this->core_id + 1,
				NO_TAG,
				MPI_COMM_WORLD,
				MPI_STATUS_IGNORE);
			}
		}
	}
}


void MPIDataProcessor::updateNext()
{
	int lastRow = ( ROOT_CORE == this->core_id ) ? this->dataSize : (this->assigned_rows + 2 * this->margin);

	for(int i = 0; i < lastRow ; i++)
	{
		for(int j = 0; j < this->dataSize; j++)
		{
			this->nextData[i][j] = this->data[i][j];
		}
	}
}




void MPIDataProcessor::shareData()
{
	/* Obtain general information */
	MPI_Comm_size(MPI_COMM_WORLD, &this->core_count);
	MPI_Comm_rank(MPI_COMM_WORLD, &this->core_id);

	/* if more than 1 core in use -> prepare infor for parallel execution */
	if( 1 < this->core_count)
	{
		/* Share required configuration data with all cores */
		MPI_Bcast(
		&this->dataSize,
		1,
		MPI_INTEGER,
		NO_TAG,
		MPI_COMM_WORLD );

		/* Calculate how many rows each corre is responsible for */
		this->assigned_rows = SET_ROW_CNT(this->dataSize - 2 * this->margin, this->core_count, this->core_id);

		/* Allocate required memory */
		if(ROOT_CORE == this->core_id)
		{
			this->nextData = tableAlloc(this->dataSize);
		}
		else
		{
			this->data     = getArray2D( this->assigned_rows + 2 * this->margin, this->dataSize );
			this->nextData = getArray2D( this->assigned_rows + 2 * this->margin, this->dataSize );
		}
		
		/* Send each core, the data it is rresponsible for */
		this->sendData();

		MPI_Barrier(MPI_COMM_WORLD);
		
		/* Update margins based on neighbour cores */
		this->syncMargins();

		MPI_Barrier(MPI_COMM_WORLD);
	}
	else
	{
		/* If only 1 cores used, skip additional stuff */
		this->nextData = tableAlloc(this->dataSize);
		this->assigned_rows = this->dataSize - 2 * this->margin;
	}

	/* Copy data to buffer array */
	this->updateNext();
}

void MPIDataProcessor::singleExecution()
{
	double *buffer = new double[dataPortionSize];

	/* If more than 1 core is used -> update margins based on neighbors */
	if( 1 < this->core_count) this->syncMargins();

	for (int row = margin; row < (this->assigned_rows + margin); row++)
	{
		for (int col = margin; col < (dataSize - margin); col++)
		{
			createDataPortion(row, col, buffer);
			nextData[row][col] = function->calc(buffer);
		}
	}

	delete[] buffer;
	double **tmp = data;
	data = nextData;
	nextData = tmp;

}

void MPIDataProcessor::collectData()
{
	/* If more than 1 core is in use -> move all data back to root core */
	if( 1 < this->core_count)
	{
		/* Initialize data for root core data storage location */
		int lastRow = SET_ROW_CNT(this->dataSize - 2 * this->margin, this->core_count, ROOT_CORE);
		lastRow += this->margin;
		int destination_core_count = lastRow;


		if(ROOT_CORE == this->core_id)
		{
			for(int i = 1; i < this->core_count; i++)
			{
				destination_core_count += SET_ROW_CNT(this->dataSize - 2 * this->margin, this->core_count, i);

				for(int j=lastRow; j < destination_core_count; j++)
				{
					MPI_Recv(
					this->data[j],
					this->dataSize,
					MPI_DOUBLE,
					i,
					NO_TAG,
					MPI_COMM_WORLD,
					MPI_STATUS_IGNORE);
				}

				lastRow = destination_core_count;
			}
		}
		else
		{
			for(int j = this->margin; j < (this->assigned_rows + this->margin); j++)
			{
				MPI_Send(
				this->data[j],
				this->dataSize,
				MPI_DOUBLE,
				ROOT_CORE,
				NO_TAG,
				MPI_COMM_WORLD);
			}	
		}
	}
}


void MPIDataProcessor::createDataPortion(int row, int col, double *buffer)
{
	int counter = 0;

	for (int i = (row - margin); i <= (row + margin); i++)
    {
		for (int j = (col - margin); j <= (col + margin); j++)
        {
			buffer[counter++] = data[i][j];
        }
    }
}

