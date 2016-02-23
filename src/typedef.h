#ifndef _HARDROD_TYPEDEF_H_
#define _HARDROD_TYPEDEF_H_
#include "Array.h"
#include "fftw++.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>

//Complex number
typedef std::complex<double> Complex;

//fftw stuff
typedef Array::array3<Complex> Ac3;
typedef Array::array3<double> Ad3;
typedef Array::array3<int> Ai3;

typedef Array::array1<Complex> Ac1;
typedef Array::array1<double> Ad1;
typedef Array::array1<int> Ai1;

// ublas
namespace bnu = boost::numeric::ublas;
typedef boost::numeric::ublas::mapped_matrix<double> sprMatD;
typedef boost::numeric::ublas::mapped_matrix< std::complex <double> > sprMatC;
typedef boost::numeric::ublas::vector<double> vD;
typedef boost::numeric::ublas::vector< std::complex <double> > vC;
typedef boost::numeric::ublas::matrix<double> matD;
typedef boost::numeric::ublas::matrix< std::complex <double> > matC;


#endif
