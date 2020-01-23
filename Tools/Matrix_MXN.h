//!#####################################################################
//! \file Matrix_MxN.h
//!#####################################################################
// Class Matrix_MxN
//######################################################################
#ifndef __Matrix_MxN__
#define __Matrix_MxN__

namespace Nova{

template<class T>
class Matrix_MxN
{
public:
    int m,n; // size of the m by n matrix
    T *x; // pointer to the one dimensional data

    Matrix_MxN()
        :m(0),n(0),x(0)
    {}

    Matrix_MxN(const int m_input,const int n_input)
        :m(m_input),n(n_input)
    {
        x=new T[m*n];
        for(int k=0;k<m*n;k++) x[k]=0;
    }

    ~Matrix_MxN()
    {delete[] x;}

    int Rows() const
    {return m;}

    int Columns() const
    {return n;}

    void Resize(const int m_new,const int n_new)
    {T* x_new=new T[m_new*n_new];for(int t=0;t<m_new*n_new;t++) x_new[t]=(T)0;
    if(!x)delete[] x;x=x_new;m=m_new;n=n_new;}

    T& operator()(const int i,const int j)
    {assert(i>=0 && i<m);assert(j>=0 && j<n);return x[j*m+i];}

    const T& operator()(const int i,const int j) const
    {assert(i>=0 && i<m);assert(j>=0 && j<n);return x[j*m+i];}

    bool Valid_Index(const int i,const int j) const
    {return i>=0 && i<m && j>=0 && j<n;}


};
template<class T>
inline Matrix_MxN<T> operator*(const T a,const Matrix_MxN<T>& A)
{return A*a;}
}
#endif
