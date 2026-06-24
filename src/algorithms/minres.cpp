#include "matrix_algo.h"
#include <cmath>

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
