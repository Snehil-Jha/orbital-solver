#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

#include <complex>
#include <functional>
#include "matrix.h"


/**
 * @brief computes the eigenvalues of a symmetric tridiagonal matrix
 * 
 * @param main main diagonal of the tridiagonal matrix (n)
 * @param sub side diagonal of the tridiagonal matrix (n - 1)
 * @param m1 index of the smallest returned eigenvalue
 * @param m2 index of the largest returned eigenvalue
 * @param result obtained eigenvalues of the system (m2 - m1 + 1)
 * @param max_iter maximum number of iterations performed per eigensvalue
 */
void bisect(
    const Vector<double>& main,
    const Vector<double>& sub,
    int m1, int m2,
    Vector<double>& result,
    const double epsilon = 1e-12,
    const int max_iter = 100
);

/**
 * @brief fills the tridiagonal matrix and recurrence vectors
 * 
 * @param H function which applies the hermitian matrix on the specified column
 * @param T stores the tridiagonal matrix (m × m). Must be all zeros at the start
 * @param Q stores the recurrence vectors (N × m). Must be all zeros at the start
 * @param epsilon numerical precision for stopping
 *
 * @return the number of iterations performed by the loop, typically m, but may be less
 */
int lanczos_get_tridiagonal(
    const std::function<void(
        const Matrix<std::complex<double>>& state, int col, Vector<std::complex<double>>& out
    )>& H,
    Vector<double>& T_main, Vector<double>& T_sub, Matrix<std::complex<double>>& Q,
    double epsilon = 1e-12
);

#endif