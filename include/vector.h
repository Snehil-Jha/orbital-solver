// vector.h
#ifndef VECTOR_H
#define VECTOR_H

#include <stdexcept>
#include <complex>
#include <string>

template <typename T>
class Vector
{
private:
    int size;
    T *data;

public:
    Vector(int n): size(n)
    {
        data = new T[n]();
    }

    Vector(const Vector& other): size(other.size)
    {
        data = new T[size];
        for(int i = 0; i < size; i++)
        {
            data[i] = other.data[i];
        }
    }

    ~Vector()
    {
        delete[] data;
    }

    Vector& operator=(const Vector& other) 
    {
        if (this != &other) {
            delete[] data;
            size = other.size;
            data = new T[size];
            for (int i = 0; i < size; i++) {
                data[i] = other.data[i];
            }
        }
        return *this;
    }

    // Move Constructor
    Vector(Vector&& other) noexcept : size(other.size), data(other.data) {
        other.size = 0;
        other.data = nullptr; // Leave in a safe, destructible state
    }

    // Move Assignment Operator
    Vector& operator=(Vector&& other) noexcept {
        if (this != &other) {
            delete[] data;
            size = other.size;
            data = other.data;
            
            other.size = 0;
            other.data = nullptr;
        }
        return *this;
    }


    T& at(int i) { 
        if(i < 0 || i >= size) 
            throw std::out_of_range("accessed element " + std::to_string(i) + " of vector of size " + std::to_string(size));
        return this->data[i]; 
    }
    const T& at(int i) const {
        if(i < 0 || i >= size) 
            throw std::out_of_range("accessed element " + std::to_string(i) + " of vector of size " + std::to_string(size));
        return this->data[i];
    }

    T& operator[](int i) { return data[i]; }
    const T& operator[](int i) const { return data[i]; }

    int length () const { return this->size; }


    static T complex_inner(const Vector<T>& a, const Vector<T>& b)
    {
        if(a.size != b.size)
            throw std::invalid_argument("lengths of vectors do not match");
        
        T val = 0;
        for(int i = 0; i < a.size; i++)
        {
            val += std::conj(a[i]) * b[i];
        }

        return val;
    }
};

template <typename T, typename U>
Vector<T> operator*(U a, const Vector<T>& b)
{
    Vector<T> vec = b;

    for(int i = 0; i < vec.length(); i++)
    {
        vec[i] *= static_cast<T>(a);
    }

    return vec;
}

template <typename T, typename U>
Vector<T> operator*(const Vector<T>& b, U a)
{
    return a * b;
}

template <typename T, typename U>
Vector<T> operator/(const Vector<T>& b, U a)
{
    Vector<T> vec = b;
    for(int i = 0; i < vec.length(); i++)
    {
        vec[i] /= static_cast<T>(a);
    }
    return vec;
}

template <typename T>
Vector<T> operator+(const Vector<T>& a, const Vector<T>& b)
{
    if(a.length() != b.length())
        throw std::invalid_argument("lengths of vectors do not match");

    auto vec = a;
    for(int i = 0; i < b.length(); i++)
    {
        vec[i] += b[i];
    }

    return vec;
}

template <typename T>
Vector<T> operator-(const Vector<T>& a, const Vector<T>& b)
{
    if(a.length() != b.length())
        throw std::invalid_argument("lengths of vectors do not match");

    auto vec = a;
    for(int i = 0; i < b.length(); i++)
    {
        vec[i] -= b[i];
    }

    return vec;
}

#endif