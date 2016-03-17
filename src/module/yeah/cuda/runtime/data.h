#ifndef _YEAH_CUDA_RUNTIME_DATA_H_
#define _YEAH_CUDA_RUNTIME_DATA_H_


// page locked memory
// pined memory
//
// Linux: mlock, munlock


// 2 API calls: what are the differences?
//cudaHostAlloc (pointer, size, flag)                       cudaHostAllocPortable: The memory returned by this call will be considered as pinned memory by all CUDA contexts, not just the one that performed the allocation
//cudaMallocHost


// zero copy, hostmapped memory



#include <assert.h>
#include <cuda_runtime.h>
#include <yeah/cuda/runtime/wrapper.h>


namespace yeah {

namespace cuda {

  template <class T0>
  class Data
  {
  public:
    T0 *dh;
    T0 *dd;
    size_t sz;
    int n;

    //Data (int num = 1) {n = num; sz = sizeof (T0) * n;}
    //~Data () {}
    void Init (int num = 1) {n = num; sz = sizeof (T0) * n;}
 
    // sync copy
    void H2Dsync () {CUDA_ERR (cudaMemcpy (dd, dh, sz, cudaMemcpyHostToDevice));}
    void D2Hsync () {CUDA_ERR (cudaMemcpy (dh, dd, sz, cudaMemcpyDeviceToHost));}
    // async copy
    void H2Dasync (cudaStream_t s = 0) {CUDA_ERR (cudaMemcpyAsync (dd, dh, sz, cudaMemcpyHostToDevice, s));}
    void D2Hasync (cudaStream_t s = 0) {CUDA_ERR (cudaMemcpyAsync (dh, dd, sz, cudaMemcpyDeviceToHost, s));}

  };


  template <class T0>
  class DataDummy: public Data <T0>
  {
  public:
    DataDummy (int num = 1) {this->n = num; this->sz = sizeof (T0) * this->n;}
    ~DataDummy () {}
  };


  template <class T0>
  class DataNonPinned: public Data <T0>
  {
  public:
    DataNonPinned (int num = 1) {this->n = num; this->sz = sizeof (T0) * this->n;}
    ~DataNonPinned () {}
    void Alloc () {
      this->dh = (T0 *) malloc (this->sz);
      assert (this->dh != NULL);
      CUDA_ERR (cudaMalloc ((void **) &this->dd, this->sz));
    }
    void Free () {
      CUDA_ERR (cudaFree (this->dd));
      free (this->dh);
    }
  };


  template <class T0>
  class DataPinned: public Data <T0>
  {
  public:
    DataPinned (int num = 1) {this->n = num; this->sz = sizeof (T0) * this->n;}
    ~DataPinned () {}
    void Alloc () {
      CUDA_ERR (cudaMallocHost ((void **) &this->dh, this->sz));
      assert (this->dh != NULL);
      CUDA_ERR (cudaMalloc ((void **) &this->dd, this->sz));
    }
    void Free () {
      CUDA_ERR (cudaFree (this->dd));
      CUDA_ERR (cudaFreeHost (this->dh));
    }
  };



  template <class T0>
  class DataHostMapped: public Data <T0>
  {
  public:
  };


  template <class T0>
  class DataAuto: public DataPinned <T0>
  {
  public:
      DataAuto (int num = 1) {this->n = num; this->sz = sizeof (T0) * this->n; this->Alloc ();}
      ~DataAuto () {this->Free ();}
  };


}
}


#endif
