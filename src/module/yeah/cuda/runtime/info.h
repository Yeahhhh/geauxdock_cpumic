#ifndef _YEAH_CUDA_RUNTIME_INFO_H_
#define _YEAH_CUDA_RUNTIME_INFO_H_


#include <cuda_runtime.h>
#include <yeah/cuda/runtime/wrapper.h>


// reference:
// nvidia/cuda7/helper_cuda.h

int _ConvertSMVer2Cores_2 (int major, int minor)
{
    typedef struct {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] =
        {
            { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
            { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
            { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
            { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
            { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
            { 0x37, 192}, // Kepler Generation (SM 3.7) GK21x class
            { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
            { 0x52, 128}, // Maxwell Generation (SM 5.2) GM20x class
            {   -1, -1 }
        };

    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1) {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
            return nGpuArchCoresPerSM[index].Cores;
        index++;
    }
  
    // If did not find the values
    printf("MapSMtoCores for SM %d.%d is undefined.\n", major, minor);
    return 0;
  
    
    // an alternative implementation
    /*
      const int sp_per_mp =
      major == 1 ? 8 :
      major == 2 ? (minor == 0 ? 32 : 48) :
      major == 3 ? 192 :
      major == 5 ? 128 :
      0;
    */
}





//cudaDeviceProp prop;
//CUDA_ERR (cudaGetDeviceProperties (&prop, d));
//PrintCudaDeviceProp
void
PrintCudaDeviceProp (const cudaDeviceProp *prop, const int d)
{

    if (prop->major == 9999 && prop->minor == 9999) {
        printf ("device %d is invalid\n", d);
    }

    const char *NOorYES[] = {"NO", "YES"}; // 0-NO, 1-YES

    const int sp_per_mp = _ConvertSMVer2Cores_2 (prop->major, prop->minor);

    const double chip_bw_Bps =
        2 * prop->memoryClockRate * 1000.0 * (prop->memoryBusWidth >> 3);

    const double chip_sp_flops =
        1000.0 * sp_per_mp * prop->clockRate * prop->multiProcessorCount;


    printf ("GPU %d: %s with CC %d.%d\n", d, prop->name, prop->major, prop->minor);
    printf ("---------------------------------------\n");
    // core
    printf ("core clock rate \t\t%.2f GHz\n", prop->clockRate / 1e6);
    printf ("MP/GPU \t\t\t\t%d\n", prop->multiProcessorCount);
    printf ("SP/MP \t\t\t\t%d\n", sp_per_mp);
    printf ("SP FLOPS/GPU \t\t\t%.2f GFLOPS/s\n", chip_sp_flops * 1e-9);
    putchar ('\n');
    // memory
    printf ("memory clock rate \t\t%.2f GHz\n", (float) prop->memoryClockRate / 1e6);
    printf ("memory bus width \t\t%d b\n", prop->memoryBusWidth);
    printf ("memory bandwidth \t\t%.2f GB/s\n", chip_bw_Bps * 1e-9);
    printf ("size reg/block \t\t\t%.3f KiB\n", (float) prop->regsPerBlock / 1024.0f );
    printf ("size scratchpad/block \t\t%.3f KiB\n", (float) prop->sharedMemPerBlock / 1024.0f);
    printf ("size L2 cache/GPU \t\t%.2f KiB\n", (float) prop->l2CacheSize / 1024.0f);
    printf ("size constant/GPU \t\t%.2f KiB\n", (float) prop->totalConstMem / 1024.0f);
    printf ("size global/GPU \t\t%.2f MiB\n", (float) prop->totalGlobalMem / (1024.0f * 1024.0f));
    printf ("\n");
    // core / memory
    printf ("SP FLOPS/BW ratio \t\t%.2f\n", 4 * chip_sp_flops / chip_bw_Bps);
    putchar ('\n');
    // feature
    printf ("unifiedAddressing: \t\t%s\n", NOorYES[prop->unifiedAddressing]);
    printf ("canMapHostMemory: \t\t%s\n", NOorYES[prop->canMapHostMemory]); // support cudaHostAlloc(), cudaHostGetDevicePointer()
    printf ("concurrentKernels: \t\t%s\n", NOorYES[prop->concurrentKernels]); // simultaneously executing multiple kernels
    printf ("deviceOverlap: \t\t\t%s\n", NOorYES[prop->deviceOverlap]); // overlaping kernel and memcpy
    printf ("asyncEngineCount: \t\t%d\n", prop->asyncEngineCount); // overlaping kernel and memcpy (single/bi-direction)
    printf ("ECCEnabled:\t\t\t%s\n", NOorYES [prop->ECCEnabled]);
    putchar ('\n');
    // CUDA restrictions
    printf ("size of warp: \t\t\t%d\n", prop->warpSize);
    printf ("max threads/block: \t\t%d\n", prop->maxThreadsPerBlock);
    printf ("kernelExecTimeoutEnabled: \t%s\n", NOorYES[prop->kernelExecTimeoutEnabled]);
    putchar ('\n');
}






void PrintCudaFuncArributes (const cudaFuncAttributes * fa, const char * s)
{
    printf ("static cuda kernel information for \"%s\"\n", s);
    printf ("ptxVersion: \t\t\t%d\n", fa->ptxVersion);
    printf ("binaryVersion: \t\t\t%d\n", fa->binaryVersion);
    printf ("reg usage/T: \t\t\t%d\n", fa->numRegs);
    printf ("lmem usage/T: \t\t\t%.3f KiB\n", (float) fa->localSizeBytes / 1024.0f);
    printf ("smem usage/TB: \t\t\t%.3f KiB\n", (float) fa->sharedSizeBytes / 1024.0f);
    printf ("cmem usage/GPU: \t\t%.3f KiB\n", (float) fa->constSizeBytes / 1024.0f);
    printf ("max T per TB: \t\t\t%d\n", fa->maxThreadsPerBlock);
    putchar ('\n');
}



//GetPrintCudaFuncArributes ((void (*)) MyKernel, "MyKernel");
void GetPrintCudaFuncArributes (void (*func), const char * s)
{
    cudaFuncAttributes fa;
    CUDA_ERR (cudaFuncGetAttributes (&fa, func));
    PrintCudaFuncArributes (&fa, s);
}

#endif
