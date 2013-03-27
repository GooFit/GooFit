// A small test program to debug a CUDA issue on Macs
#include <iostream>

__device__ double testDouble;
double hostDouble = 5.1;

int main (int argc, char** argv) {
 cudaError_t err = cudaMemcpyToSymbol(testDouble, (void*) &hostDouble, sizeof(double));
 std::cout << cudaGetErrorString(err) << std::endl;

 return 0;
}
