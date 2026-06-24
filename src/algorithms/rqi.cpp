#include "matrix_algo.h"
#include <random>

int symmetric_tridag_rqi(
    const double eigenvalue,
    Vector<double> &main,
    const Vector<double> &sub,

    double& value,
    Vector<double>& vector,

    const double tolerance,
    const int max_iter
)
{
    const int n = main.length();
    if(sub.length() != n-1 || vector.length() != n)
        throw std::invalid_argument("array size mismatch");

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    value = eigenvalue;

    // generates a random guess eigenvector
    double vec_norm_sq = 0;
    for(int i = 0; i < n; i++)
    {
        vector[i] = dis(gen);
        vec_norm_sq += vector[i] * vector[i];
    }

    double vec_norm = sqrt(vec_norm_sq);

    for(int i = 0; i < n; i++)
    {
        vector[i] /= vec_norm;
    }


    // iteration
    auto tmp = Vector<double>(n-1);

    for(int iter_count = 0; iter_count < max_iter; iter_count++)
    {
        // shift the main diagonal
        for(int i = 0; i < n; i++)
            main[i] -= value;

        // solve the resulting system
        symmetric_thomas_solver(main, sub, vector, vector, tmp);

        // normalize the vector
        vec_norm_sq = 0;
        for(int i = 0; i < n; i++)
            vec_norm_sq += vector[i] * vector[i];

        vec_norm = sqrt(vec_norm_sq);
        for(int i = 0; i < n; i++)
            vector[i] /= vec_norm;        

        // compute the shift to the current eigenvalue and the current residual
        double term = main[0] * vector[0] + sub[0] * vector[1];
        double res_sq = term * term;
        double val_shift = vector[0] * term;
        for(int i = 1; i < n - 1; i++)
        {
            term = main[i] * vector[i] + sub[i] * vector[i+1] + sub[i-1] * vector[i-1];
            res_sq += term * term;
            val_shift += vector[i] * term;
        }
        term = main[n-1] * vector[n-1] + sub[n-2] * vector[n-2];
        res_sq += term * term;
        val_shift += vector[n-1] * term;
        
        // restores the main diagonal
        for(int i = 0; i < n; i++)
            main[i] += value;

        // update the eigenvalue
        value += val_shift;

        if(res_sq < tolerance * tolerance)
            return iter_count;

    }

    return max_iter;
}


int symmetric_tridag_rqi(
    const double eigenvalue,
    Vector<double> &main,
    const Vector<double> &sub,
    const Vector<double> &weight,

    double& value,
    Vector<double>& vector,

    const double tolerance,
    const int max_iter
)
{
    const int n = main.length();
    if(sub.length() != n-1 || vector.length() != n)
        throw std::invalid_argument("array size mismatch");

    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);

    value = eigenvalue;

    // generates a random guess eigenvector
    double vec_norm_sq = 0;
    for(int i = 0; i < n; i++)
    {
        vector[i] = dis(gen);
        vec_norm_sq +=  weight[i] * vector[i] * vector[i];
    }

    double vec_norm = sqrt(vec_norm_sq);

    for(int i = 0; i < n; i++)
    {
        vector[i] /= vec_norm;
    }


    // iteration
    auto tmp = Vector<double>(n-1);

    for(int iter_count = 0; iter_count < max_iter; iter_count++)
    {
        // shift the main diagonal
        for(int i = 0; i < n; i++)
        {
            main[i] -= value * weight[i];
            vector[i] *= weight[i];
        }

        // solve the resulting system
        symmetric_thomas_solver(main, sub, vector, vector, tmp);

        // normalize the vector
        vec_norm_sq = 0;
        for(int i = 0; i < n; i++)
            vec_norm_sq += weight[i] * vector[i] * vector[i];

        vec_norm = sqrt(vec_norm_sq);
        for(int i = 0; i < n; i++)
            vector[i] /= vec_norm;        

        // compute the shift to the current eigenvalue and the current residual
        double term = main[0] * vector[0] + sub[0] * vector[1];
        double res_sq = term * term;
        double val_shift = vector[0] * term;
        for(int i = 1; i < n - 1; i++)
        {
            term = main[i] * vector[i] + sub[i] * vector[i+1] + sub[i-1] * vector[i-1];
            res_sq += term * term;
            val_shift += vector[i] * term;
        }
        term = main[n-1] * vector[n-1] + sub[n-2] * vector[n-2];
        res_sq += term * term;
        val_shift += vector[n-1] * term;
        
        // restores the main diagonal
        for(int i = 0; i < n; i++)
            main[i] += value * weight[i];

        // update the eigenvalue
        value += val_shift;

        if(res_sq < tolerance * tolerance)
            return iter_count;

    }

    return max_iter;
}