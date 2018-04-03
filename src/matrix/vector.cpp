#include <cstdio>
#include <cstdlib>
#include <cstring> // for memset
#include <limits>
#include <iostream>
#include <vector>
#include <math.h>


class Vector {
 
  public:

    Vector() : size(0), data(NULL) {}
 
    Vector(int size_) : Vector() {
      size = size_;
      allocate(size_);
    }
 
    ~Vector() {
      deallocate();
    }
 
    double & operator() (int i) {
      return data[i]; }
    double  operator() (int i) const {
      return data[i]; }
 
    Vector& operator=(const Vector& source) {
 

      if (this != &source) { 
        if ( size != (source.size) ) {  
	        allocate(source.size);         
        }

        std::copy(source.data, source.data + source.size, data);
      }
      return *this;
    }
 
    void allocate(int size_) {
 
      deallocate();
 
      size = size_;
      data = new double[size_];
      memset(data,0,size_*sizeof(double));
 
    }
 
    void deallocate() {
 
      if (data)
        delete[] data;
 
      data = NULL;
 
    }    
 
    double norm() {
      double sum = 0;
      for (int i = 0; i < size; i++) sum += (*this)(i) * (*this)(i);
      return sqrt(sum);
    }
 
    void rescale(double factor) {
      for (int i = 0; i < size; i++) (*this)(i) /= factor;
    }
 
    void rescale_unit() {
      double factor = norm();
      rescale(factor);
    }
 
    int size;
 
  private:
    double* data;
 
};