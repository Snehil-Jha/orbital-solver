#include "system.h"

#include "vector.h"
#include "matrix.h"
#include "matrix_solvers.h"

void RadialSystem::update_potential(const std::function<double(double)>& updated_V, int updated_l)
{
    if(updated_l != -1)
        l = updated_l;


    for(int i = 0; i < N; i++)
    {
        ham_main[i] = kinetic_main[i] + ri_sq[i] * updated_V(ri[i]) + 0.5 * l * (l + 1);
    }
}

void RadialSystem::solve_energies(const double shift, Vector<double> &energies, const int m)
{
    const int k = energies.length();

    auto T_main = Vector<double>(m);
    auto T_sub = Vector<double>(m - 1);

    auto Q = Matrix<double>(N, m);

    auto shifted_main = Vector<double>(N);
    for(int i = 0; i < N; i++)
    {
        shifted_main[i] = ham_main[i] - shift * ri_sq[i];
    }

    auto result = Vector<double>(k);


    auto H_shift_and_invert = [this, &shifted_main](const Matrix<double>& state, const int col, Vector<double>& out){
        
        for(int i = 0; i < N; i++)
        {
            out[i] = ri[i] * state[i, col];
        }

        symmetric_thomas_solver(shifted_main, ham_sub, out, out, thomas_tmp);

        for(int i = 0; i < N; i++)
        {
            out[i] = ri[i] * out[i];
        }
    };

    lanczos(H_shift_and_invert, T_main, T_sub, Q);    

    bisect(T_main, T_sub, m-k, m-1, result);

    for(int i = 0; i < k; i++)
    {
        energies[i] = shift + (1. / result[k - 1 - i]);
    }
}
