#include <stdio.h>
#include <iostream>
#include <string>
#include <mpi.h>

using namespace std;

class Domain
{
public:
  Domain(int _M, int _N, const char *_name="") : domain(new char[(_M+1)*(_N+1)]), M(_M), N(_N), name(_name)  {}
  virtual ~Domain() {delete[] domain;}
  char &operator()(int i, int j) {return domain[i*N+j];}
  char operator()(int i, int j)const {return domain[i*N+j];}

  int rows() const {return M;}
  int cols() const {return N;}

  const string & myname() const {return name;}

protected:
  char *domain;
  int M;
  int N;

  string name;
};

void parrallelize(int P,int Q,char** grid, int N, int iterations, int size, int myrank, MPI_Comm comm);
void zero_domain(Domain &domain);
void print_domain(Domain &domain);
void update_domain(Domain &new_domain, Domain &old_domain, int size, int myrank, MPI_Comm comm);

int main(int argc, char* argv[]) {
	
	int ndims = 2;
        int rank, size, reorder, my_cart_rank, ierr, nrows, ncols;
	int dims[ndims];
        int coords[ndims];
	int wrap_around;
	MPI_Comm comm2D;
	int N,P,Q;
	int iterations;
        MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	int array[4];
	if(rank == 0)
	{
		N = atoi(argv[1]);
        	P = atoi(argv[2]);
		Q = atoi(argv[3]);
		iterations = atoi(argv[4]);
		array[0] = N;
		array[1] = P;
		array[2] = Q;
		array[3] = iterations;
	}
	MPI_Bcast(array,4,MPI_INT,0,MPI_COMM_WORLD);
	if(rank != 0)
	{
		N = array[0];
		P = array[1];
		Q = array[2];
		iterations = array[3];
	}


	int local_N = N/P;
	int local_M = N/Q;
	
	char** local_grid = (char**)malloc(local_N * sizeof(char*));
       	for (int i = 0; i < local_N; i++) 
	{
        local_grid[i] = (char*)malloc(local_M * sizeof(char));
	}
	
	for (int i = 0; i < local_N; i++) {
        	for (int j = 0; j < local_M; j++) {
            		if ((i + j) % 2 == 0) {
                		local_grid[i][j] = '*';
            		} else {
                		local_grid[i][j] = ' ';
            		}
        	}
	}

	//2D Depomposition  of an NxN world of cells (let space be a dead space and ’*’ is a live space).
	//Represent cells as char	
	parrallelize(local_N,local_M,local_grid,N,iterations,size,rank,MPI_COMM_WORLD);
	for (int i = 0; i < local_N; i++) {

        free(local_grid[i]);

   	 }

    	free(local_grid);
	MPI_Finalize();

        return 0;

}

void parrallelize(int P, int Q,char** local_grid, int N, int iterations, int size, int rank, MPI_Comm COMM_2D)
{
	int m = P;//colomns
	int n = Q;//rows
	
	Domain even_domain(m,n,"even Domain");
  	Domain odd_domain(m,n,"odd Domain");
	
	zero_domain(even_domain);
  	zero_domain(odd_domain);	
	if((n >= 8) && (m >= 10))
 	 {
   		 even_domain(0,(n-1)) = 1;
   		 even_domain(0,0)     = 1;
   		 even_domain(0,1)     = 1;
    
   		 even_domain(3,5) = 1;
   		 even_domain(3,6) = 1;
   		 even_domain(3,7) = 1;

   		 even_domain(6,7) = 1;
   		 even_domain(7,7) = 1;
   		 even_domain(8,7) = 1;
   		 even_domain(9,7) = 1;
  	}
	 print_domain(even_domain);

  	Domain *odd, *even; // pointer swap magic
  	odd = &odd_domain;
  	even = &even_domain;
	for(int i = 0; i < iterations;i++){
		update_domain(*odd, *even, N, rank,COMM_2D);
		cout << "Iteration #" << i << endl;
		print_domain(*odd);
		
    		// swap pointers:
    		Domain *temp = odd;
    		odd  = even;
    		even = temp;

	}

}
void zero_domain(Domain &domain)
{
  for(int i = 0; i < domain.rows(); ++i)
    for(int j = 0; j < domain.cols(); ++j)
      domain(i,j) = 0;
}
void print_domain(Domain &domain)
{
  cout << domain.myname() << ":" <<endl;
  // this is naive; it doesn't understand big domains at all
  for(int i = 0; i < domain.rows(); ++i)
  {
    for(int j = 0; j < domain.cols(); ++j)
      cout << (domain(i,j) ? "*" : " ");
    cout << endl;
  }
}

void update_domain(Domain &new_domain, Domain &old_domain, int size, int myrank, MPI_Comm comm)
{
  MPI_Request request[4];

  int neighbor_count;
  int N = new_domain.cols();
  int m = new_domain.rows();

  char *top_row = new char[N];
  char *bottom_row = new char[N];

  char *top_halo = new char[N];
  char *bottom_halo = new char[N];

  // int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
  //             MPI_Comm comm, MPI_Request *request)
  for(int i = 0; i < 4; ++i) {
      request[i] = MPI_REQUEST_NULL;
  }
  // fill the top row
  for(int i = 0; i < N; ++i)
  {
      top_row[i] = old_domain(0,i);
  }
  MPI_Isend(top_row, N, MPI_CHAR, (myrank-1+size)%size, 0, comm, &request[0]);

  for(int i = 0; i < N; ++i)
  {
     bottom_row[i] = old_domain(m-1,i);
  }
  MPI_Isend(bottom_row, N, MPI_CHAR, (myrank+1)%size, 0, comm, &request[1]);

  // send my top row and bottom row to adjacent process
  // receive the halo from top and bottom process
  MPI_Irecv(top_halo, N, MPI_CHAR,    (myrank+size-1)%size, 0, comm, &request[2]);
  MPI_Irecv(bottom_halo, N, MPI_CHAR, (myrank+1)%size, 0, comm, &request[3]);

  MPI_Waitall(4, request, MPI_STATUSES_IGNORE); // complete all 4 transfers

  // at this point, I have halos from my neighbors

  // i=0:
  for(int j = 0; j < new_domain.cols(); ++j)
  {
      neighbor_count = 0;

      for(int delta_i = 0; delta_i <= 1; delta_i++)
      {
	for(int delta_j = -1; delta_j <= 1; delta_j++)
	{
	  if(delta_i == 0 && delta_j == 0) //skip self
	    continue;

	  // this first implementation is sequental and wraps the vertical
	  // and horizontal dimensions without dealing with halos (ghost cells)
	  if(old_domain((0+delta_i),
			(j+delta_j+old_domain.cols())%old_domain.cols()))
	     ++neighbor_count;
	}
      }
      for(int delta_j = -1; delta_j <= 1; delta_j++)
      {
        if(top_halo[(j+delta_j+old_domain.cols())%old_domain.cols()])
	  ++neighbor_count;
      }

      char mycell = old_domain(0,j);
      char newcell = 0;
      if(mycell == 0)
	newcell = (neighbor_count == 3) ? 1 : 0;
      else
	newcell = ((neighbor_count == 2)||(neighbor_count == 3)) ? 1 : 0;

      new_domain(0,j) = newcell;
  } // int j

  // i=m-1:
  for(int j = 0; j < new_domain.cols(); ++j)
  {
      neighbor_count = 0;

      for(int delta_i = -1; delta_i <= 0; delta_i++)
      {
	for(int delta_j = -1; delta_j <= 1; delta_j++)
	{
	  if(delta_i == 0 && delta_j == 0) //skip self
	    continue;

	  // this first implementation is sequental and wraps the vertical
	  // and horizontal dimensions without dealing with halos (ghost cells)
	  if(old_domain((m-1+delta_i),
			(j+delta_j+old_domain.cols())%old_domain.cols()))
	     ++neighbor_count;
	}
      }
      for(int delta_j = -1; delta_j <= 1; delta_j++)
      {
        if(bottom_halo[(j+delta_j+old_domain.cols())%old_domain.cols()])
	  ++neighbor_count;
      }

      char mycell = old_domain(m-1,j);
      char newcell = 0;
      if(mycell == 0)
	newcell = (neighbor_count == 3) ? 1 : 0;
      else
	newcell = ((neighbor_count == 2)||(neighbor_count == 3)) ? 1 : 0;

      new_domain(0,j) = newcell;
  } // int j

  // these update as in sequential case:
  for(int i = 1; i < (new_domain.rows()-1); ++i)
  {
    for(int j = 0; j < new_domain.cols(); ++j)
    {
      neighbor_count = 0;
      for(int delta_i = -1; delta_i <= 1; delta_i++)
      {
	for(int delta_j = -1; delta_j <= 1; delta_j++)
	{
	  if(delta_i == 0 && delta_j == 0) //skip self
	    continue;

	  // this first implementation is sequental and wraps the vertical
	  // and horizontal dimensions without dealing with halos (ghost cells)
	  if(old_domain((i+delta_i+old_domain.rows())%old_domain.rows(),
			(j+delta_j+old_domain.cols())%old_domain.cols()))
	     ++neighbor_count;

	}
      }
      char mycell = old_domain(i,j);
      char newcell = 0;
      if(mycell == 0)
	newcell = (neighbor_count == 3) ? 1 : 0;
      else
	newcell = ((neighbor_count == 2)||(neighbor_count == 3)) ? 1 : 0;

      new_domain(i,j) = newcell;
    } // int j
  } // int i


}

