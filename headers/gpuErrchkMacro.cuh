#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// #define gpuCurandErrchk(x) do { if((x)!=CURAND_STATUS_SUCCESS)
// { \ printf("Error at %s:%d\n",__FILE__,__LINE__);\ return EXIT_FAILURE;}} while(0)

// The maximum value of float = 3.4028234664e+38
// The minimum value of float = 1.1754943508e-38