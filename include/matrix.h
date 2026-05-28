// matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <stdexcept>
#include <complex>
#include <string>
#include "vector.h"

template <typename T>
class Matrix
{
private:
    int m, n;
    T* data;

public:
    Matrix(int im, int in): m(im), n(in)
    {
        data = new T[im * in]();
    }

    
    Matrix(const Matrix& other) : m(other.m), n(other.n) 
    {
        data = new T[m * n];
        for (int i = 0; i < m * n; i++) {
            data[i] = other.data[i];
        }
    }

    Matrix& operator=(const Matrix& other) 
    {
        if (this != &other) {
            delete[] data;
            m = other.m;
            n = other.n;
            data = new T[m * n];
            for (int i = 0; i < m * n; i++) {
                data[i] = other.data[i];
            }
        }
        return *this;
    }

    ~Matrix()
    {
        delete[] data;
    }

    // Move Constructor
    Matrix(Matrix&& other) noexcept : m(other.m), n(other.n), data(other.data) {
        other.m = 0;
        other.n = 0;
        other.data = nullptr; // Leave in a safe, destructible state
    }

    // Move Assignment Operator
    Matrix& operator=(Matrix&& other) noexcept {
        if (this != &other) {
            delete[] data;
            m = other.m;
            n = other.n;
            data = other.data;
            
            other.m = 0;
            other.n = 0;
            other.data = nullptr;
        }
        return *this;
    }

    // Column-major storing of data, for cache optimization
    T& at(int i, int j)
    {
        if(i < 0 || i >= m || j < 0 || j >= n)
            throw std::out_of_range("accessed element (" + std::to_string(i) + ", " + std::to_string(j) + 
                        ") of matrix of dimension (" + std::to_string(m) + ", " + std::to_string(n) + ")");
        
        return data[(j * m) + i];
    }
    const T& at(int i, int j) const
    {
        if(i < 0 || i >= m || j < 0 || j >= n)
            throw std::out_of_range("accessed element (" + std::to_string(i) + ", " + std::to_string(j) + 
                        ") of matrix of dimension (" + std::to_string(m) + ", " + std::to_string(n) + ")");
        
        return data[(j * m) + i];
    }

    // will be compiled on c++23 or above
    T& operator[](int i, int j) {
        return data[(j * m) + i];
    }

    const T& operator[](int i, int j) const {
        return data[(j * m) + i];
    }


    Vector<T> row(int i) const
    {
        Vector<T> vec = Vector<T>(n);

        for(int j = 0; j < n; j++)
        {
            vec[j] = (*this)[i, j];
        }

        return vec;
    }

    Vector<T> col(int j) const
    {
        Vector<T> vec = Vector<T>(m);

        for(int i = 0; i < m; i++)
        {
            vec[i] = (*this)[i, j];
        }

        return vec;
    }

    T col_inner_product(const int& j, const Vector<T>& vec) const
    {
        if(vec.length() != m)
            throw std::invalid_argument("invalid dimension in column multiplication");


        T val = 0;
        for(int i = 0; i < m; i++)
        {
            val += std::conj((*this)[i, j]) * vec[i];
        }

        return val;
    }

    static Matrix constant(int im, int in, T value)
    {
        Matrix mat(im, in);
        for(int i = 0; i < im * in; i++) {
            mat.data[i] = value;
        }
        return mat; 
    }

    static Matrix identity(int n)
    {
        Matrix mat(n, n);

        for(int i = 0; i < n; i++) {
            mat.data[i * (n + 1)] = 1;
        }
        
        return mat;
    }
};

#endif