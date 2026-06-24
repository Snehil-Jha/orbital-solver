#include "matrix_algo.h"
#include <random>

int lanczos(
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

int lanczos_invert_and_shift(
    const std::function<void(const Vector<double>& in, Vector<double>& out)>& H,
    Vector<double>& T_main, Vector<double>& T_sub, Matrix<double>& Q,
    const double shift, const double epsilon, const double tol, const int max_iter
)
{
    int N = Q.row_count();

    auto H_shift = [&shift, &H, &N](const Vector<double>& in, Vector<double>& out){
        H(in, out);

        for(int i = 0; i < N; i++)
        {
            out[i] -= shift * in[i];
        }
    };

    auto solver = Minres(H_shift, N);

    return lanczos([&solver, &tol, &max_iter](const Matrix<double>& state, int col, Vector<double>& out){

        solver.solve(state, col, out, tol, max_iter);

    }, T_main, T_sub, Q);
}
