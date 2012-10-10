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
#include "vmlib.h"

#ifndef VMLIB_CPP
#define VMLIB_CPP

Vector::Vector(unsigned int s)
{
    v = gsl_vector_calloc(s);
};

Vector::Vector(unsigned int s, double value)
{
    v = gsl_vector_alloc(s);
    gsl_vector_set_all(v, value);
};

Vector::~Vector()
{
    gsl_vector_free(v);
};

unsigned int Vector::size()
{
    return v->size;
};

double Vector::operator[](unsigned int i)
{
    return gsl_vector_get(v, i);
};

void Vector::operator()(unsigned int i, double d)
{
    gsl_vector_set(v, i, d);
};

bool Vector::operator==(Vector& vec)
{
    return gsl_vector_equal(v, vec.v);
};

void Vector::operator=(Vector& vec)
{
    gsl_vector_memcpy(v, vec.v);
};

Vector Vector::operator+(Vector& vec)
{
    Vector result(v->size);
    result = *this;
    gsl_vector_add(result.v, vec.v);
    return result;
};

Vector& Vector::operator+=(Vector& vec)
{
    gsl_vector_add(v, vec.v);
    return *this;
};

Vector Vector::operator-(Vector& vec)
{
    Vector result(v->size);
    result = *this;
    gsl_vector_sub(result.v, vec.v);
    return result;
};

Vector& Vector::operator-=(Vector& vec)
{
    gsl_vector_sub(v, vec.v);
    return *this;
};

double Vector::operator*(Vector& vec)
{
    double result;
    gsl_blas_ddot(v, vec.v, &result);
    return result;
};

Vector Vector::operator+(double d)
{
    Vector result(v->size);
    result = *this;
    gsl_vector_add_constant(result.v, d);
    return result;
};

Vector& Vector::operator+=(double d)
{
    gsl_vector_add_constant(v, d);
    return *this;
};

Vector Vector::operator-(double d)
{
    Vector result(v->size);
    result = *this;
    gsl_vector_add_constant(result.v, -d);
    return result;
};

Vector& Vector::operator-=(double d)
{
    gsl_vector_add_constant(v, -d);
    return *this;
};

Vector Vector::operator/(double d)
{
    Vector result(v->size);
    result = *this;
    gsl_vector_scale(result.v, 1.0 / d);
    return result;
};

Vector& Vector::operator/=(double d)
{
    gsl_vector_scale(v, 1.0 / d);
    return *this;
};

Vector Vector::operator*(double d)
{
    Vector result(v->size);
    result = *this;
    gsl_vector_scale(result.v, d);
    return result;
};

Vector& Vector::operator*=(double d)
{
    gsl_vector_scale(v, 1.0 / d);
    return *this;
};

double Vector::norm()
{
    return gsl_blas_dnrm2(v);
};

double Vector::mean()
{
    double result = 0;
    for(unsigned int i = 0; i < v->size; i++) {
        result += gsl_vector_get(v, i);
    }
    return result / v->size;
};

string Vector::__repr__()
{
    stringstream result;
    result << '[';
    for(unsigned int i = 0; i < v->size; i++) {
        result << gsl_vector_get(v, i);
        if(i < v->size - 1) {
            result << ", ";
        }
    }
    result << ']';
    return result.str();
};

Matrix::Matrix(unsigned int s1, unsigned int s2)
{
    m = gsl_matrix_calloc(s1, s2);
};

Matrix::Matrix(unsigned int s1, unsigned int s2 , double d)
{
    m = gsl_matrix_alloc(s1, s2);
    gsl_matrix_set_all(m, d);
};

Matrix::~Matrix()
{
    gsl_matrix_free(m);
};

unsigned int Matrix::get_m()
{
    return m->size1;
};

unsigned int Matrix::get_n()
{
    return m->size2;
};

bool Matrix::operator==(Matrix& mat)
{
    return gsl_matrix_equal(m, mat.m);
};

void Matrix::operator=(Matrix& mat)
{
    gsl_matrix_memcpy(m, mat.m);
};

Matrix Matrix::operator+(Matrix& mat)
{
    Matrix result(m->size1, m->size2);
    result = *this;
    gsl_matrix_add(result.m, mat.m);
    return result;
};

Matrix& Matrix::operator+=(Matrix& mat)
{
    gsl_matrix_add(m, mat.m);
    return *this;
};

Matrix Matrix::operator-(Matrix& mat)
{
    Matrix result(m->size1, m->size2);
    result = *this;
    gsl_matrix_sub(result.m, mat.m);
    return result;
};

Matrix& Matrix::operator-=(Matrix& mat)
{
    gsl_matrix_sub(m, mat.m);
    return *this;
};

Matrix Matrix::operator*(Matrix& mat)
{
    Matrix result(m->size1, mat.m->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m, mat.m, 0.0, result.m);
    return result;
};

Matrix& Matrix::operator*=(Matrix& mat)
{
    gsl_matrix* result = gsl_matrix_calloc(m->size1, mat.m->size2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m, mat.m, 0.0, result);
    gsl_matrix_free(m);
    m = result;
    return *this;
};

Vector Matrix::operator*(Vector& vec)
{
    Vector result(m->size1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, m, vec.v, 0.0, result.v);
    return result;
};

double Matrix::operator()(unsigned int i, unsigned int j)
{
    return gsl_matrix_get(m, i, j);
};

void Matrix::operator()(unsigned int i, unsigned int j, double d)
{
    gsl_matrix_set(m, i, j, d);
};

double Matrix::det()
{
    return gsl_linalg_LU_det(m, 0);
};

double Matrix::mean()
{
    double result = 0;
    for(unsigned int i = 0; i < m->size1; i++) {
        for(unsigned int j = 0; j < m->size2; j++) {
            result += gsl_matrix_get(m, i, j);
        }
    }
    return result / (m->size1 * m->size2);
};

void Matrix::diagonalize()
{
    gsl_vector* eval = gsl_vector_alloc(m->size1);
    gsl_matrix* evec = gsl_matrix_alloc(m->size1, m->size1);

    gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc(m->size1);
    gsl_eigen_symmv(m, eval, evec, w);
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(evec);

    for(unsigned int i = 0; i < m->size1; i++) {
        for(unsigned int j = 0; j < m->size2; j++) {
            gsl_matrix_set(m, i, j, i == j ? gsl_vector_get(eval, m->size1 - i - 1) : 0);
        }
    }
    gsl_vector_free(eval);
};

string Matrix::__repr__()
{
    stringstream result;
    result << '[';
    for(unsigned int i = 0; i < m->size1; i++) {
        if(i > 0) {
            result << ' ';
        }
        for(unsigned int j = 0; j < m->size2; j++) {
            result << gsl_matrix_get(m, i, j);
            if(j < m->size2 - 1) {
                result << ", ";
            }
        }
        if(i < m->size1 - 1) {
            result << '\n';
        }
    }
    result << ']';
    return result.str();
};

#endif
