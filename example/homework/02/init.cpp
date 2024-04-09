using namespace std;
#include <iostream>
#include <assert.h>
#include <string>

#include <mpi.h>

// forward declarations:

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

void zero_domain(Domain &domain);
void print_domain(Domain &domain);
void update_domain(Domain &new_domain, Domain &old_domain, int size, int myrank, MPI_Comm comm);
void parallel_code(int P, int Q, int N, int iterations, int size, int myrank, MPI_Comm comm);

int main(int argc, char **argv)
{
  int P, Q, N;
  int iterations;

  if(argc < 4)
  {
    cout << "usage: " << argv[0] << " M N iterations" << endl;
    exit(0);
  }

  int size, myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int array[3];
  if(myrank == 0)
  {
     P = atoi(argv[1]); Q = atoi(argv[2]);  N = atoi(argv[3]); iterations = atoi(argv[4]);

     array[0] = P;
     array[1] = Q;
     array[2] = N;
     array[3] = iterations;
    
  }
  MPI_Bcast(array, 4, MPI_INT, 0, MPI_COMM_WORLD);
  if(myrank != 0)
  {
    P = array[0];
    Q = array[1];
    N = array[2];
    iterations = array[3];
  }

  
  parallel_code(P, Q, N, iterations, size, myrank, MPI_COMM_WORLD);
  
  MPI_Finalize();
}

void parallel_code(int P, int Q, int N, int iterations, int size, int myrank, MPI_Comm comm)
{
    int rows_per_proc = N / P;
    int cols_per_proc = N / Q;
    int extra_rows = N % P;
    int extra_cols = N % Q;


    int n = (myrank / Q < extra_rows) ? (rows_per_proc + 1) : rows_per_proc;
    int m = (myrank % Q < extra_cols) ? (cols_per_proc + 1) : cols_per_proc;
  
  Domain even_domain(m,n,"even Domain");
  Domain odd_domain(m,n,"odd Domain");

  zero_domain(even_domain);
  zero_domain(odd_domain);

  // fill in even_domain with something meaningful (initial state)
  // this requires min size for default values to fit:
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

  // here is where I might print out my picture of the initial domain
  cout << "Initial:"<<endl; print_domain(even_domain);

  Domain *odd, *even; // pointer swap magic
  odd = &odd_domain;
  even = &even_domain;

  for(int i = 0; i < iterations; ++i)
  {
    update_domain(*odd, *even, size, myrank, comm);
    // here is where I might print out my picture of the updated domain
    cout << "Iteration #" << i << endl; print_domain(*odd);

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
  MPI_Request request[16];
  
  int neighbor_count;
  int n = new_domain.cols();
  int m = new_domain.rows();

  char *top_row = new char[n];
  char *right_row = new char[m];
  char *left_row = new char[m];
  char *bottom_row = new char[n];
  
  char *top_halo = new char[n];
  char *right_halo = new char[m];
  char *left_halo = new char[m];
  char *bottom_halo = new char[n];

  const int top = 0,topleft = 1,topright = 2,right = 3,left = 4, bottomleft = 5, bottomright = 6, bottom = 7;

  // int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
  //             MPI_Comm comm, MPI_Request *request)

  // fill the top row
  for(int i = 0; i < n; ++i)
  {
      top_row[i] = old_domain(0,i);
  }
  MPI_Isend(top_row, n, MPI_CHAR, (myrank-1+size)%size, bottom, comm, &request[0]);
  for(int i = 0; i < m; ++i)
  {
      right_row[i] = old_domain(i,m-1);
  }
  MPI_Isend(right_row, m, MPI_CHAR, (myrank-1+size)%size, left, comm, &request[1]);
  for(int i = 0; i < m; ++i)
  {
      left_row[i] = old_domain(i,0);
  }
  MPI_Isend(top_row, m, MPI_CHAR, (myrank-1+size)%size, right, comm, &request[2]);
  for(int i = 0; i < n; ++i)
  {
     bottom_row[i] = old_domain(m-1,i);
  }
  MPI_Isend(bottom_row, n, MPI_CHAR, (myrank+1)%size, top, comm, &request[3]);
  
  
  char *topleft_row = new char[1];
  char *topright_row  = new char[1];
  char *bottomleft_row  = new char[1];
  char *bottomright_row  = new char[1];

  topleft_row[0] = old_domain(1,1);
  topright_row[0] = old_domain(m,1);
  bottomleft_row[0] = old_domain(1,n);
  bottomright_row[0] = old_domain(m,n);
  
  MPI_Isend(topleft_row, 1, MPI_CHAR, (myrank+1)%size, topleft, comm, &request[4]);
  MPI_Isend(topright_row, 1, MPI_CHAR, (myrank+1)%size, topright, comm, &request[5]);
  MPI_Isend(bottomleft_row, 1, MPI_CHAR, (myrank+1)%size, topleft, comm, &request[6]);
  MPI_Isend(bottomright_row, 1, MPI_CHAR, (myrank+1)%size, topright, comm, &request[7]);

  char topleft_halo[1];
  char topright_halo[1];
  char bottomleft_halo[1];
  char bottomright_halo[1];

  // send my top row and bottom row to adjacent process
  // receive the halo from top and bottom process
  MPI_Irecv(top_halo, n, MPI_CHAR,    (myrank+size-1)%size, top, comm, &request[0]);
  MPI_Irecv(topleft_halo, 1, MPI_CHAR,    (myrank+size-1)%size, topleft, comm, &request[1]);
  MPI_Irecv(topright_halo, 1, MPI_CHAR,    (myrank+size-1)%size, topright, comm, &request[2]);
  MPI_Irecv(right_halo, m, MPI_CHAR,    (myrank+size-1)%size, right, comm, &request[3]);
  MPI_Irecv(left_halo, m, MPI_CHAR,    (myrank+size-1)%size, left, comm, &request[4]);
  MPI_Irecv(bottomleft_halo, 1, MPI_CHAR,    (myrank+size-1)%size, bottomleft, comm, &request[5]);
  MPI_Irecv(bottomright_halo, 1, MPI_CHAR,    (myrank+size-1)%size, bottomright, comm, &request[6]);
  MPI_Irecv(bottom_halo, n, MPI_CHAR, (myrank+1)%size, bottom, comm, &request[7]);

  MPI_Waitall(16, request, MPI_STATUSES_IGNORE); // complete all 4 transfers

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
  for(int j = 0; j < new_domain.rows(); ++j)
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
          if(old_domain(((j+delta_j+old_domain.rows())%old_domain.rows()),
                        (0+delta_i)))
             ++neighbor_count;
        }
      }
      for(int delta_j = -1; delta_j <= 1; delta_j++)
      {
        if(left_halo[(j+delta_j+old_domain.cols())%old_domain.cols()])
          ++neighbor_count;
      }
      char mycell = old_domain(0,j);
      char newcell = 0;
      if(mycell == 0)
        newcell = (neighbor_count == 3) ? 1 : 0;
      else
        newcell = ((neighbor_count == 2)||(neighbor_count == 3)) ? 1 : 0;
      new_domain(j,0) = newcell;
  }
  for(int j = 0; j < new_domain.rows(); ++j)
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
          if(old_domain(((j+delta_j+old_domain.rows())%old_domain.rows()),
                        (0+delta_i)))
             ++neighbor_count;
        }
      }
      for(int delta_j = -1; delta_j <= 1; delta_j++)
      {
        if(right_halo[(j+delta_j+old_domain.cols())%old_domain.cols()])
          ++neighbor_count;
      }
      char mycell = old_domain(0,j);
      char newcell = 0;
      if(mycell == 0)
        newcell = (neighbor_count == 3) ? 1 : 0;
      else
        newcell = ((neighbor_count == 2)||(neighbor_count == 3)) ? 1 : 0;
      new_domain(j,n-1) = newcell;
  }
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
      
      new_domain(m-1,j) = newcell;
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



