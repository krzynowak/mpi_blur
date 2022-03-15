#include "DataProcessor.h"


class MPIDataProcessor : public DataProcessor
{
private:
	void createDataPortion( int row, int col, double *buffer );

	int core_count;
	int core_id;
	int assigned_rows;

	void sendData();
	void syncMargins();
	void updateNext();

protected:
	void singleExecution();
	void collectData();
	void shareData();

public:
	double** getResult() {
		return data;
	}

};
