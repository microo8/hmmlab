#include <string>
#include <sstream>
#include <gsl_math.h>
#include <gsl_vector.h>
#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gsl_linalg.h>
#include <gsl_eigen.h>

using namespace std;

#ifndef VMLIB_H
#define VMLIB_H


class Matrix;

class Vector
{
    gsl_vector* v;
public:
    Vector(unsigned int);
    Vector(unsigned int, double);
    ~Vector();
    unsigned int size();
    double operator[](unsigned int);
    void operator()(unsigned int, double);
    bool operator==(Vector&);
    void operator=(Vector&);
    Vector operator+(Vector&);
    Vector& operator+=(Vector&);
    Vector operator-(Vector&);
    Vector& operator-=(Vector&);
    double operator*(Vector&);
    Vector operator+(double);
    Vector& operator+=(double);
    Vector operator-(double);
    Vector& operator-=(double);
    Vector operator/(double);
    Vector& operator/=(double);
    Vector operator*(double);
    Vector& operator*=(double);
    double norm();
    double mean();
    string __repr__();

    friend class Matrix;
};

class Matrix
{
    gsl_matrix* m;
public:
    Matrix(unsigned int, unsigned int);
    Matrix(unsigned int, unsigned int, double);
    ~Matrix();
    unsigned int get_m();
    unsigned int get_n();
    bool operator==(Matrix&);
    void operator=(Matrix&);
    Matrix operator+(Matrix&);
    Matrix& operator+=(Matrix&);
    Matrix operator-(Matrix&);
    Matrix& operator-=(Matrix&);
    Matrix operator*(Matrix&);
    Matrix& operator*=(Matrix&);
    Vector operator*(Vector&);
    double operator()(unsigned int, unsigned int);
    void operator()(unsigned int, unsigned int, double);
    double det();
    double mean();
    void diagonalize();
    string __repr__();
};

#endif
