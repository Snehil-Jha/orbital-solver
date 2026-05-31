#include <algorithm>
#include <cstdlib>
#include <random>
#include <functional>
#include <stdexcept>
#include <limits>
#include <cmath>

#include "matrix.h"
#include "vector.h"

#include "matrix_solvers.h"


void bisect(
    const Vector<double>& main,
    const Vector<double>& sub,
    int m1, int m2,
    Vector<double>& result,
    const double epsilon, const int max_iter
)
{
    // verification
    int n = main.length();
    if(sub.length() != n - 1 || result.length() != m2 - m1 + 1)
        throw std::invalid_argument("vectors of unexpected sizes obtained");


    // calculation of bounds of eigenvalues
    double xmin = main[0] - std::abs(sub[0]);
    double xmax = main[0] + std::abs(sub[0]);
    double h;
    for(int i = 1; i < n; i++)
    {
        h = std::abs(sub[i - 1]);
        if(i != n - 1)
        {
            h += std::abs(sub[i]);
        }

        xmax = std::max(xmax, main[i] + h);
        xmin = std::min(xmin, main[i] - h);
    }

    // start the bisection method

    // create array of lower bounds. (upper bounds go inside results)
    auto lower = Vector<double>(result.length(), xmin);
    for(int i = 0; i < result.length(); i++) result[i] = xmax;

    double x_upper = xmax;

    // loop for each eigenvalue
    for(int k = m2 - m1; k >= 0; k--)
    {
        double x_lower = xmin;
        for(int i = k; i >= 0; i--)
        {
            if(x_lower < lower[i])
            {
                x_lower = lower[i];
                break;
            }
        }

        double x_mid;
        for(int iter_count = 0; iter_count < max_iter; iter_count++)
        {
            x_mid = 0.5 * (x_upper + x_lower);
            if(x_upper - x_lower < epsilon * std::max(1., std::abs(x_mid))) break;

            // sturms sequence
            int a = 0;
            double q = 1;
            for(int i = 0; i < n; i++)
            {
                double sub_element = (i == 0) ? 0 : sub[i-1];

                if(q == 0) q = main[i] - x_mid - std::abs(sub_element) / std::numeric_limits<double>::epsilon();
                else q = main[i] - x_mid - sub_element * sub_element / q;

                if(q < 0) a++;
            }

            if(a <= k + m1)
            {
                if(a < m1)
                    x_lower = lower[0] = x_mid;
                else
                    x_lower = lower[a - m1] = x_mid;
            }
            else
            {
                x_upper = x_mid;
            }
        }

        result[k] = 0.5 * (x_lower + x_upper);
    }
}

int lanczos_get_tridiagonal(
    const std::function<void(
        const Matrix<double>& state, int col, Vector<double>& out
    )>& H,
    Vector<double>& T_main, Vector<double>& T_sub, Matrix<double>& Q,
    double epsilon
)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0); 

    int N = Q.row_count();
    int m = T_main.length();

    if(T_sub.length() != m - 1 || Q.col_count() != m)
        throw std::invalid_argument("Provided matrices not of appropriate shape.");

    // generate the random vector
    double sq_sum = 0;
    for(int j = 0; j < N; j++)
    {
        double val = dis(gen);

        Q[j, 0] = val;
        sq_sum += val * val;
    }

    double normalization = 1 / sqrt(sq_sum);
    for(int j = 0; j < N; j++)
    {
        Q[j, 0] *= normalization;
    }

    double beta_old = 0;

    auto w = Vector<double>(N);

    // iteration loop
    for(int j = 0; j < m; j++)
    {
        // inital vector
        H(Q, j, w);

        // primary orthogonalization, with raw loops for speed
        if(j != 0)
        {
            for(int i = 0; i < N; i++)
            {
                w[i] = w[i] - beta_old * Q[i, j-1];
            }
        }

        double alpha_j = Q.col_dot_product(j, w);

        for(int i = 0; i < N; i++)
        {
            w[i] = w[i] - alpha_j * Q[i, j];
        }

        // full reorthogonalization, again with raw loops
        for(int k = 0; k < j; k++)
        {
            double inner = Q.col_dot_product(k, w); 
            
            for(int i = 0; i < N; i++)
            {
                w[i] = w[i] - inner * Q[i, k];
            }
        }

        // normalize and store
        beta_old = sqrt(std::max(0., Vector<double>::dot_product(w, w)));
        
        T_main[j] = alpha_j; // main diagonal
        if(j < m-1)
        {
            T_sub[j] = beta_old;
        }

        if(beta_old < epsilon)
        {
            return j + 1;
        }

        if(j < m-1)
        {
            for(int i = 0; i < N; i++)
            {
                Q[i, j+1] = w[i] / beta_old;
            }
        }
    }

    return m;
}


int Minres::solve(const Matrix<double>& state, const int col_b, Vector<double>& x, const double tol, const int max_iter)
{
    // initilization
    double rho_sq = 0;
    for(int i = 0; i < N; i++)
    {
        x[i] = 0;
        r[i] = state[i, col_b];
        rho_sq += r[i] * r[i];
        v_old[i] = 0;
        w[i] = 0;
    }

    double rho = std::sqrt(rho_sq);
    double norm_factor = 1 / rho;
    for(int i = 0; i < N; i++)
    {
        v[i] = r[i] * norm_factor;
        wtt[i] = v[i];
    }

    double beta = 0, betat = 0, c = -1, s = 0;

    int iter = 0;

    while(std::abs(rho) > tol && iter < max_iter)
    {
        ++iter;
        // Lanczos
        A(v, Av);
        for(int i = 0; i < N; i++)
        {
            v_hat[i] = Av[i] - beta * v_old[i];
        }

        double alpha = Vector<double>::dot_product(v, v_hat);

        double beta_sq = 0;
        for(int i = 0; i < N; i++)
        {
            v_hat[i] = v_hat[i] - alpha * v[i];
            beta_sq += v_hat[i] * v_hat[i];
        }

        beta = std::sqrt(beta_sq);
        if(beta < tol)
        {
            break;
        }

        // QR-decomposition of the Lanczos matrix
        double l1= s * alpha - c * betat;
        double l2 = s * beta;
        double alphat = -s * betat - c * alpha;
        betat = c * beta;
        double l0 = std::sqrt(alphat * alphat + beta * beta);
        c = alphat / l0;
        s = beta / l0;

        for(int i = 0; i < N; i++)
        {
            v_old[i] = v[i];
            v[i] = v_hat[i] / beta;

            // The search vector 
            wt[i] = wtt[i] - l1 * w[i];
            wtt[i] = v[i] - l2 * w[i];
            w[i] = wt[i] / l0;

            // the approximate solution
            x[i] += rho * c * w[i];
        }

        rho = s * rho;
    }

    return iter;
}