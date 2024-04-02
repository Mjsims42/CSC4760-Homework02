
#include <iostream>
#include <stdio.h>
#include <mpi.h>



int main(int argc, char** argv) {
		int rank, size;
		int N = 5;

		MPI_Init(&argc,&argv);

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		int last_rank = (rank + size - 1) % size;
		int next_rank = (rank + 1) % size;

		int number = 0;
		for(int i = 0;i < N;i++)
		{
			if(rank == 0)
			{
			MPI_Send(&number,1,MPI_INT,next_rank,0, MPI_COMM_WORLD);
			}
			else{
			MPI_Recv(&number,1, MPI_INT,last_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			if(rank)
			{
			printf("process %d recieves number %d from process %d\n",rank,number,last_rank);
			}
			number++;
			MPI_Send(&number,1,MPI_INT,next_rank,0, MPI_COMM_WORLD);
			}
		MPI_Barrier(MPI_COMM_WORLD);
		}
		// do something in scope
		MPI_Finalize();
    	return 0;

}
