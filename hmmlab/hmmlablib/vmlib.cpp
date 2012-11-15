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
#include <iostream>

#ifndef VMLIB_CPP
#define VMLIB_CPP

gsl_matrix* pca(const gsl_matrix* data, unsigned int L)
{
    /*
    @param data - matrix of data vectors, MxN matrix, each column is a data vector, M - dimension, N - data vector count
    @param L - dimension reduction
    */
    assert(data != NULL);
    assert(L > 0 && L < data->size2);
    unsigned int i;
    unsigned int rows = data->size1;
    unsigned int cols = data->size2;
    gsl_vector* mean = gsl_vector_alloc(rows);

    for(i = 0; i < rows; i++) {
        gsl_vector_set(mean, i, gsl_stats_mean(data->data + i * cols, 1, cols));
    }

    // Get mean-substracted data into matrix mean_substracted_data.
    cout << " Get mean-substracted data into matrix mean_substracted_data.\n";
    gsl_matrix* mean_substracted_data = gsl_matrix_alloc(rows, cols);
    gsl_matrix_memcpy(mean_substracted_data, data);
    for(i = 0; i < cols; i++) {
        gsl_vector_view mean_substracted_point_view = gsl_matrix_column(mean_substracted_data, i);
        gsl_vector_sub(&mean_substracted_point_view.vector, mean);
    }
    gsl_vector_free(mean);

    // Compute Covariance matrix
    cout << " Compute Covariance matrix\n";
    gsl_matrix* covariance_matrix = gsl_matrix_alloc(rows, rows);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0 / (double)(cols - 1), mean_substracted_data, mean_substracted_data, 0.0, covariance_matrix);
    gsl_matrix_free(mean_substracted_data);

    // Get eigenvectors, sort by eigenvalue.
    cout << " Get eigenvectors, sort by eigenvalue.\n";
    gsl_vector* eigenvalues = gsl_vector_alloc(rows);
    gsl_matrix* eigenvectors = gsl_matrix_alloc(rows, rows);
    gsl_eigen_symmv_workspace* workspace = gsl_eigen_symmv_alloc(rows);
    gsl_eigen_symmv(covariance_matrix, eigenvalues, eigenvectors, workspace);
    gsl_eigen_symmv_free(workspace);
    gsl_matrix_free(covariance_matrix);

    // Sort the eigenvectors
    cout << " Sort the eigenvectors\n";
    gsl_eigen_symmv_sort(eigenvalues, eigenvectors, GSL_EIGEN_SORT_ABS_DESC);
    gsl_vector_free(eigenvalues);

    // Project the original dataset
    cout << " Project the original dataset\n";
    gsl_matrix_view L_eigenvectors = gsl_matrix_submatrix(eigenvectors, 0, 0, rows, L);
    gsl_matrix* result = gsl_matrix_alloc(rows, L);
    gsl_matrix_memcpy(result, &L_eigenvectors.matrix);
    gsl_matrix_free(eigenvectors);

    // Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L
    cout << " Result is n LxN matrix, each column is the original data vector with reduced dimension from M to L\n";
    return result;
}

Vector::Vector(gsl_vector* vv): v(vv) {};

Vector::Vector(unsigned int s)
{
    v = gsl_vector_calloc(s);
};

Vector::Vector(unsigned int s, double value)
{
    if(value != 0) {
        v = gsl_vector_alloc(s);
        gsl_vector_set_all(v, value);
    } else {
        v = gsl_vector_calloc(s);
    }
};

Vector::~Vector()
{
    gsl_vector_free(v);
};

unsigned int Vector::size()
{
    return v->size;
};

gsl_vector* Vector::get_vector()
{
    gsl_vector* result = gsl_vector_alloc(v->size);
    gsl_vector_memcpy(result, v);
    return result;
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

Matrix::Matrix(gsl_matrix* mm): m(mm) {};

Matrix::Matrix(unsigned int s1, unsigned int s2)
{
    m = gsl_matrix_calloc(s1, s2);
};

Matrix::Matrix(unsigned int s1, unsigned int s2 , double d)
{
    if(d != 0) {
        m = gsl_matrix_alloc(s1, s2);
        gsl_matrix_set_all(m, d);
    } else {
        m = gsl_matrix_calloc(s1, s2);
    }
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

gsl_matrix* Matrix::get_matrix()
{
    gsl_matrix* result = gsl_matrix_alloc(m->size1, m->size2);
    gsl_matrix_memcpy(result, m);
    return result;
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
