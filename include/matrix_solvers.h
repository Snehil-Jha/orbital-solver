#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

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
 * @param epsilon numerical precision for stopping midway
 *
 * @return the number of iterations performed by the loop, typically m, but may be less
 */
int lanczos(
    const std::function<void(
        const Matrix<double>& state, int col, Vector<double>& out
    )>& H,
    Vector<double>& T_main, Vector<double>& T_sub, Matrix<double>& Q,
    double epsilon = 1e-12
);

/**
 * @brief fills the tridiagonal matrix and recurrence vectors
 * 
 * @param H function which applies the hermitian matrix on the specified column
 * @param T stores the tridiagonal matrix (m × m). Must be all zeros at the start
 * @param Q stores the recurrence vectors (N × m). Must be all zeros at the start
 * @param shift the shift used for invert and shift method
 * @param epsilon numerical precision for stopping midway
 *
 * @return the number of iterations performed by the loop, typically m, but may be less
 */
int lanczos_invert_and_shift(
    const std::function<void(const Vector<double>& in, Vector<double>& out)>& H,
    Vector<double>& T_main, Vector<double>& T_sub, Matrix<double>& Q,
    const double shift, const double epsilon = 1e-12, const double tol = 1e-12, const int max_iter=1000
);


class Minres
{
    Vector<double> r;
    Vector<double> v;
    Vector<double> v_old;
    Vector<double> v_hat;
    Vector<double> w;
    Vector<double> wt;
    Vector<double> wtt;

    Vector<double> Av;

    const std::function<void(const Vector<double>& in, Vector<double>& out)> A;

    const int N;

    public:
    /**
     * @brief Construct a new Minres object
     * 
     * @param iA Function which applies A on the given vector
     * @param iN The dimensionality of the matrix system
     */
    Minres(
        const std::function<void(const Vector<double>& in, Vector<double>& out)>& iA,
        const int iN
    ): r(iN), v(iN), v_old(iN), v_hat(iN), w(iN), wt(iN), wtt(iN), Av(iN), A(iA), N(iN)
    {}

    /**
     * @brief Solves the system Ax = b
     * 
     * @param state Matrix which contains b
     * @param col_b the column in Q in which b is stored
     * @param x vector to store the final output
     * @param tol Tolerance of the solver
     * @param max_iter Maximum number of iterations the solver will go through
     *
     * @returns the number of iterations taken to reach tol
     */
    int solve(const Matrix<double>& state, const int col_b, Vector<double>& x, const double tol = 1e-12, const int max_iter = 1000);
};

/**
 * @brief Implementation of the Thomas algorithm, specifically for symmetric tridiagonal matrices Ax = b
 * 
 * @param main Main diagonal of the matrix A (n)
 * @param sub Sub diagonal of the matrix A (n-1)
 * @param b RHS of the equation (n)
 * @param x solution of the algorithm (n)
 * @param tmp a temporary vector needed for intermediate calculations (n - 1)
 * @param epsilon threshold for the denominator of the solver
 *
 * @returns true if the system encountered a zero denominator, false otherwise
 */
bool symmetric_thomas_solver(
    const Vector<double> &main,
    const Vector<double> &sub,
    const Vector<double> &b,
    Vector<double>& x,
    Vector<double>& tmp,
    const double epsilon = 1e-16
);

/**
 * @brief Uses the Rayleigh Quotient Iteration to find the eigenvector corresponding to the specified eigenvalue
 * 
 * @param eigenvalue approximate guess to the initial eigenvalue
 * @param main main diagonal of the matrix NOTE: This will be restored in the end, but might be modified during the function call
 * @param sub sub diagonal of the matrix
 * @param value better approximation to the eigenvalue
 * @param vector computed eigenvector
 *
 * @param tolerance margin for the residual
 * @param max_iter maximum number of iterations the program runs for
 *
 * @return the number of iterations taken by the algorithm
 */
int symmetric_tridag_rqi(
    const double eigenvalue,
    Vector<double> &main,
    const Vector<double> &sub,

    double& value,
    Vector<double>& vector,

    const double tolerance = 1e-12,
    const int max_iter = 100
);

#endif