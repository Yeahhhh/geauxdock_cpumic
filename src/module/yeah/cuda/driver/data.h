
#ifndef _YEAH_CUDA_DRIVER_DATA_H_
#define _YEAH_CUDA_DRIVER_DATA_H_

#include <assert.h>
#include <cuda.h>
#include <yeah/cuda/driver/wrapper.h>


namespace yeah {

namespace cu {

  template <class T0>
  class Data
  {
  public:
    T0 *dh;
    CUdeviceptr dd;
    size_t sz;
    int n;

    //Data (int num = 1) {n = num; sz = sizeof (T0) * n;}
    //~Data () {}
    void Init (int num = 1) {n = num; sz = sizeof (T0) * n;}
 
    // sync copy
    void H2Dsync () {CU_ERR (cuMemcpyHtoD (dd, dh, sz));}
    void D2Hsync () {CU_ERR (cuMemcpyDtoH (dh, dd, sz));}
    // async copy
    void H2Dasync (CUstream s = 0) {CU_ERR (cuMemcpyHtoDAsync (dd, dh, sz, s));}
    void D2Hasync (CUstream s = 0) {CU_ERR (cuMemcpyDtoHAsync (dh, dd, sz, s));}
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
      CU_ERR (cuMemAlloc (&this->dd, this->sz));
    ;}
    void Free () {
      CU_ERR (cuMemFree (this->dd));
      free (this->dh);
    }
  };


  template <class T0>
  class DataPinned: public Data <T0>
  {
  public:
  };



  template <class T0>
  class DataHostMapped: public Data <T0>
  {
  public:
  };


}
}


#endif
