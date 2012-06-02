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
    Vector(int);
    Vector(int, double);
    ~Vector();
    int size();
    double operator[](int);
    void operator()(int, double);
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
    Matrix(int, int);
    Matrix(int, int, double);
    ~Matrix();
    int get_m();
    int get_n();
    bool operator==(Matrix&);
    void operator=(Matrix&);
    Matrix operator+(Matrix&);
    Matrix& operator+=(Matrix&);
    Matrix operator-(Matrix&);
    Matrix& operator-=(Matrix&);
    Matrix operator*(Matrix&);
    Matrix& operator*=(Matrix&);
    Vector operator*(Vector&);
    double operator()(int, int);
    void operator()(int, int, double);
    double det();
    double mean();
    void diagonalize();
    string __repr__();
};

#endif
