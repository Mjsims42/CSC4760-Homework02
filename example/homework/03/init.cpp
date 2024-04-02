#include <Kokkos_Core.hpp>
#include <cstdio>
#include <iostream>

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
	int n = 2000;
  // Make View of length n > 1000 and add values
  	Kokkos::View<int*> A("A",n); 
       	Kokkos::parallel_for("Loop1", n, KOKKOS_LAMBDA (const int i) {
                A(i) = i*i;
        });	
	Kokkos::fence();
  // create two additional views of same size and datatype
	Kokkos::View<int*> B("B",n);
	Kokkos::View<int*> C("C",n);
	Kokkos::Timer timer;
  // deep_copy
  	Kokkos::deep_copy(B,A);
	double time0 = timer.seconds();
	timer.reset();
  // user copy
  	Kokkos::parallel_for("usercopy", n, KOKKOS_LAMBDA (const int i) {
                B(i) = A(i);
        });
	double time1 = timer.seconds();
	std::cout << "deep_copy Time: " << time0 << std::endl;
  	std::cout << "User_copy Time: " << time1 << std::endl;
  // Output times
  /* 
   *From the looks like, deep_copy is faster than parrallel for on a cuda backend, this is probably due to cuda's optimization
   *
   */ 
  }
  Kokkos::finalize();
}
