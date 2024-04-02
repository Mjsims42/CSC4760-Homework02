#include<iostream>
#include<mpi.h>
#include <stdio.h>



int main(int argc, char* argv[])
{
	int size,rank;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	int local_value = rank + 1; // Each process has a different value
    	int result;
	int global_sum;
	int recv_data[size];
	printf("REDUCE ALL:\n");
	//BroadCast value among process
	MPI_Bcast(&local_value,1,MPI_INT,0,MPI_COMM_WORLD);
	//Perform a reduction
	MPI_Reduce(&local_value,&result,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	//BroadCast the result to all processes
	MPI_Bcast(&result,1,MPI_INT,0,MPI_COMM_WORLD);

	printf("Using Bcast and Reduce for reduce all. Rank %d: Result = %d\n", rank, result);
	
	printf("GATHER ALL:\n");
	MPI_Gather(&local_value, 1, MPI_INT, recv_data, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&recv_data,1,MPI_INT,0,MPI_COMM_WORLD);
	printf("Rank %d: Received data:\n", rank);
	for (int i = 0; i < size; i++) {
        	printf("Rank %d: %d\n", i, recv_data[i]);
	}
	
	printf("ALL TO ALL:\n");
	MPI_Alltoall(&local_value, 1, MPI_INT, recv_data, 1, MPI_INT, MPI_COMM_WORLD);

	printf("Rank %d: Received data:\n", rank);

    	for (int i = 0; i < size; i++) {
        	printf("Rank %d: Received from Rank %d: %d\n", rank, i, recv_data[i]);
    	}

	MPI_Finalize();

}
