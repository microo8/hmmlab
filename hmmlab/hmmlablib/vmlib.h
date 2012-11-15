/*
    This file is part of HMMLab.

    HMMLab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HMMLab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HMMLab.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <assert.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <gsl_math.h>
#include <gsl_matrix.h>
#include <gsl_blas.h>
#include <gsl_linalg.h>
#include <gsl_eigen.h>
#include <gsl/gsl_statistics.h>

using namespace std;

#ifndef VMLIB_H
#define VMLIB_H

void gsl_vector_print(gsl_vector*);
void gsl_matrix_print(gsl_matrix*);
gsl_matrix* gsl_pca(const gsl_matrix*, unsigned int);

class Matrix;

class Vector
{
    gsl_vector* v;
public:
    Vector(gsl_vector*);
    Vector(unsigned int);
    Vector(unsigned int, double);
    ~Vector();
    unsigned int size();
    gsl_vector* get_vector();
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
    Matrix(gsl_matrix*);
    Matrix(unsigned int, unsigned int);
    Matrix(unsigned int, unsigned int, double);
    ~Matrix();
    unsigned int get_m();
    unsigned int get_n();
    gsl_matrix* get_matrix();
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
