#include <complex>
#include <random>
#include <functional>

#include "matrix.h"
#include "vector.h"



void lanczos(
    const std::function<void(
        const Matrix<std::complex<double>>& state, int col, Vector<std::complex<double>>& out
    )>& H,
    int N, int m, double epsilon=1e-12)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0); 

    auto Q = Matrix<std::complex<double>>(N, m+1);

    // generate the random vector
    double sq_sum = 0;
    for(int j = 0; j < N; j++)
    {
        double real = dis(gen);
        double imag = dis(gen);

        Q[j, 1] = std::complex<double>(real, imag);
        sq_sum += real * real + imag * imag;
    }

    double normalization = 1 / sqrt(sq_sum);
    for(int j = 0; j < N; j++)
    {
        Q[j, 1] *= normalization;
    }


    auto T = Matrix<double>(m, m);
    double beta_old = 0;

    auto w = Vector<std::complex<double>>(N);

    // iteration loop
    for(int j = 1; j <= m; j++)
    {
        // inital vector
        H(Q, j, w);

        // primary orthogonalization, with raw loops for
        for(int i = 0; i < N; i++)
        {
            w[i] = w[i] - beta_old * Q[i, j-1];
        }

        std::complex<double> alpha_j = Q.col_inner_product(j, w);

        for(int i = 0; i < N; i++)
        {
            w[i] = w[i] - alpha_j * Q[i, j];
        }

        // full reorthogonalization, again with raw loops
        for(int k = 1; k <= j; k++)
        {

            std::complex<double> inner = Q.col_inner_product(k, w); 
            
            for(int i = 0; i < N; i++)
            {
                w[i] = w[i] - inner * Q[i, k];
            }
        }

        // normalize and store
        beta_old = sqrt(std::real(Vector<std::complex<double>>::complex_inner(w, w)));
        
        T[j - 1, j - 1] = std::real(alpha_j); // main diagonal
        if(j != m)
        {
            T[j, j - 1] = beta_old;
            T[j - 1, j] = beta_old;
        }

        if(beta_old < epsilon)
            break;

        if(j < m)
        {
            for(int i = 0; i < N; i++)
            {
                Q[i, j+1] = w[i] / beta_old;
            }
        }
    }


}