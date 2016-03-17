#ifndef _YEAH_OPENCL_DATA_HPP_
#define _YEAH_OPENCL_DATA_HPP_


#include <assert.h>
#include <yeah/opencl/cl.h>


namespace yeah {

namespace opencl {

  template <class T0>
  class Data0
  {
  public:
    T0 *dh;
    cl_mem dd;
    size_t sz;
    int n;

    //Data (int num = 1) {n = num; sz = sizeof (T0) * n;}
    //~Data () {}
    void Init (int num = 1) {n = num; sz = sizeof (T0) * n; }
 
  };


  template <class T0>
  class Data: public Data0 <T0>
  {
  public:
    Data(int num = 1) {this->n = num; this->sz = sizeof (T0) * this->n;}
    ~Data() {}

    void AllocH () { this->dh = (T0 *) malloc (this->sz); assert (this->dh != NULL); }
    void FreeH () { free (this->dh); }
    //GE (cudaMalloc ((void **) &this->dd, this->sz))
    //clReleaseMemObject (this->dd);
  };


}
}

#endif
