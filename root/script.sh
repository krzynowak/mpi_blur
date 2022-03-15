#!/bin/sh

# cd src
mpiCC src/*.cpp -o bin/prog.out && 
echo "Running for $CORES cores" &&
time mpirun --allow-run-as-root -np $CORES --oversubscribe ./bin/prog.out &&
rm bin/*
