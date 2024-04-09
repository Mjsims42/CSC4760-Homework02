
#include <iostream>
#include <stdio.h>
#include <mpi.h>



int main(int argc, char** argv) {
		int rank, size;
		int N = 5;
		int flag;
		MPI_Status status;
		MPI_Init(&argc,&argv);

		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		int next_rank = (rank + 1) % size;
		int last_rank = (rank - 1 + size) % size;
		int number = 0;
		for(int i = 0;i < N; i++)
		{
			if(rank == 0)
			{
				MPI_Iprobe(last_rank, 0, MPI_COMM_WORLD, &flag, &status);
				if(flag)
				{
					MPI_Recv(&number,1, MPI_INT,3,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				}
				MPI_Send(&number,1,MPI_INT,next_rank,0, MPI_COMM_WORLD);	
                                number++;
			}
			else
			{
				MPI_Recv(&number,1, MPI_INT,last_rank,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				printf("process %d recieves number %d from process %d\n",rank,number,last_rank);
				if(rank == 3)
				{
					MPI_Send(&number,1,MPI_INT,0,0, MPI_COMM_WORLD);
				}
				else
				{
				MPI_Send(&number,1,MPI_INT,next_rank,0, MPI_COMM_WORLD);
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);	
		// do something in scope
		MPI_Finalize();
    	return 0;

}
