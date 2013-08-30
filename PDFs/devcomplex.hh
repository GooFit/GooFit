#ifndef DEVCOMPLEX_HH
#define DEVCOMPLEX_HH

template <typename T> struct devcomplex {
  T real;
  T imag; 

  __device__ devcomplex<T> () : real(0), imag(0) {} 
  __device__ devcomplex<T> (T r, T i) : real(r), imag(i) {}

  __device__ devcomplex<T>& operator= (T other) {
    real = other;
    imag = 0; 
    return *this; 
  }

    __device__ devcomplex<T>& operator+= (const devcomplex<T>& other) {
    real += other.real;
    imag += other.imag;
    return *this;
  }

    __device__ devcomplex<T>& operator-= (const devcomplex<T>& other) {
    real -= other.real;
    imag -= other.imag;
    return *this;
  }

  __device__ devcomplex<T>& operator*= (T mult) {
    real *= mult;
    imag *= mult;
    return *this;
  }

  __device__ devcomplex<T>& operator/= (T mult) {
    real /= mult;
    imag /= mult;
    return *this;
  }

  __device__ devcomplex<T>& operator*= (const devcomplex<T>& other) {
    multiply(other.real, other.imag);
    return *this;
  }

  __device__ devcomplex<T>& multiply (const T other_real, const T other_imag) {
    T nreal = real * other_real - imag * other_imag;
    T nimag = real * other_imag + imag * other_real; 
    real = nreal;
    imag = nimag; 
    return *this; 
  }

  __device__ T abs2 () {return real*real + imag*imag;} 
};

template <typename T> __device__ devcomplex<T> operator+ (const devcomplex<T>& one, const devcomplex<T>& other) {
  return devcomplex<T>(one.real+other.real, one.imag+other.imag); 
}

template <typename T> __device__ devcomplex<T> operator- (const devcomplex<T>& one, const devcomplex<T>& other) {
  return devcomplex<T>(one.real-other.real, one.imag-other.imag); 
}

template <typename T> __device__ devcomplex<T> operator+ (const devcomplex<T>& one, T other) {
  return devcomplex<T>(one.real+other, one.imag); 
}

template <typename T> __device__ devcomplex<T> operator/ (const devcomplex<T>& one, const devcomplex<T>& other) {
  T inverse(1);
  inverse /= (other.real*other.real + other.imag*other.imag);
  return devcomplex<T>(inverse*(one.real*other.real+one.imag*other.imag), inverse*(one.imag*other.real - one.real*other.imag)); 
}

template <typename T> __device__ devcomplex<T> operator* (const devcomplex<T>& one, const devcomplex<T>& other) {
  return devcomplex<T>(one.real*other.real - one.imag*other.imag, one.real*other.imag + one.imag*other.real); 
}

template <typename T> __device__ devcomplex<T> operator* (const devcomplex<T>& one, const int& other) {
  return devcomplex<T>(one.real*other, one.imag*other); 
}

template<typename T> __device__ devcomplex<T> operator/ (const T& a, const devcomplex<T>& b) {
  T inverse(1);
  inverse /= (b.real*b.real + b.imag*b.imag);
  return devcomplex<T>(inverse*a*b.real, -inverse*a*b.imag); 
}

template<typename T> __device__ devcomplex<T> operator* (T a, const devcomplex<T>& b) {
  return devcomplex<T>(a*b.real, a*b.imag); 
}

template <typename T> __device__ T norm (const devcomplex<T>& z) {
  return SQRT(z.real*z.real + z.imag*z.imag); 
}

template <typename T> __device__ T norm2 (const devcomplex<T>& z) {
  return (z.real*z.real + z.imag*z.imag); 
}

template <typename T> __device__ devcomplex<T> exp (const devcomplex<T>& z) {
  T mult = EXP(z.real); 
  return devcomplex<T>(mult*COS(z.imag), mult*SIN(z.imag)); 
}

template <typename T> __device__ devcomplex<T> conj (const devcomplex<T>& z) {
  return devcomplex<T>(z.real, -(z.imag)); 
}

#endif
