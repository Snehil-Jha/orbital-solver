#include "system.h"

#include "vector.h"
#include "matrix.h"
#include "matrix_algo.h"
#include <cmath>
#include <stdexcept>


template <typename F>
void RadialSystem::update_potential(const F&& updated_V, int updated_l)
{
    if(updated_l != -1)
        l = updated_l;


    for(int i = 0; i < N; i++)
    {
        ham_main[i] = kinetic_main[i] + ri_sq[i] * updated_V(ri[i]) + 0.5 * l * (l + 1);
    }
}

void RadialSystem::solve_energy(Vector<double> &energies)
{
    const int k = energies.length();

    bisect(ham_main, ham_sub, ri_sq, 0, k-1, energies, 1e-12, 1000);
}


void RadialSystem::solve_wavefunction(Vector<double> &energies, Matrix<double> &wavefunctions)
{
    const int k_energy = energies.length();
    const int k_state = wavefunctions.row_count();

    if(wavefunctions.col_count() != N || k_energy < k_state)
        throw std::invalid_argument("invalid matrix size");

    bisect(ham_main, ham_sub, ri_sq, 0, k_energy - 1, energies, 1e-12, 1000);

    auto current_eigenvector = Vector<double>(N);

    for (int i = 0; i < k_state; i++)
    {
        // solve for the wavefunction
        double current_energy = energies[i];

        symmetric_tridag_rqi(
            current_energy,

            ham_main, ham_sub,
            ri_sq,
            
            current_energy, current_eigenvector,

            1e-12, 1000
        );
        
        // Gram-Schmidt to reduce numerical errors
        for(int lower = 0; lower < i; lower++)
        {
            double overlap = 0;
            for(int j = 0; j < N; j++)
            {
                overlap += wavefunctions[lower, j] * current_eigenvector[j] * ri_sq[j];
            }

            for(int j = 0; j < N; j++)
            {
                current_eigenvector[j] -= overlap * wavefunctions[lower, j];
            }
        }
        
        // Renormalize
        double norm_sq = 0;
        for(int j = 0; j < N; j++)
        {
            norm_sq += current_eigenvector[j] * current_eigenvector[j] * ri_sq[j];
        }
        double norm = sqrt(norm_sq);
        for(int j = 0; j < N; j++)
        {
            current_eigenvector[j] /= norm;
        }
        
        // Recompute energy again
        double final_energy = 0;
        for(int j = 0; j < N; j++)
        {
            double H_psi = ham_main[j] * current_eigenvector[j];
            if(j > 0) H_psi += ham_sub[j-1] * current_eigenvector[j-1];
            if(j < N-1) H_psi += ham_sub[j] * current_eigenvector[j+1];

            final_energy += current_eigenvector[j] * H_psi;
        }
        energies[i] = final_energy;

        for(int j = 0; j < N; j++)
        {
            wavefunctions[i, j] = current_eigenvector[j];
        }
    }
}

double RadialSystem::radial_coordinate(int index)
{
    return exp(x_min + index * dx);
}
