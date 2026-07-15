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
        pot_main[i] = ri_sq[i] * updated_V(ri[i]);
        ham_main[i] = kinetic_main[i] + pot_main[i] + 0.5 * l * (l + 1);
    }
}


void RadialSystem::update_azimuthal(int updated_l)
{
    l = updated_l;

    for(int i = 0; i < N; i++)
    {
        ham_main[i] = kinetic_main[i] + pot_main[i] + 0.5 * l * (l + 1);
    }
}


const void RadialSystem::solve_energy(Vector<double> &energies)
{
    const int k = energies.length();

    bisect(ham_main, ham_sub, ri_sq, 0, k-1, energies, 1e-12, 1000);
}

const void RadialSystem::solve_wavefunction(Vector<double> &energies, Matrix<double> &wavefunctions, const int index_start, int state_count, const bool warm_start)
{
    if(state_count == -1)
    {
        if(wavefunctions.row_count() != energies.length())
            throw std::invalid_argument("invalid matrix size");
        state_count = energies.length();
    }

    auto energies_bisect = Vector<double>(state_count);

    bisect(ham_main, ham_sub, ri_sq, 0, state_count - 1, energies_bisect, 1e-12, 1000);

    auto current_eigenvector = Vector<double>(N);

    for (int i = index_start; i < index_start + state_count; i++)
    {
        // solve for the wavefunction
        double current_energy = energies_bisect[i - index_start];
        
        if(warm_start)
        {
            for(int j = 0; j < N; j++)
            {
                current_eigenvector[j] = wavefunctions[i, j];
            }
        }

        symmetric_tridag_rqi(
            current_energy,

            ham_main, ham_sub,
            ri_sq,
            
            current_energy, current_eigenvector,

            1e-12, 1000,

            warm_start
        );
        
        // Gram-Schmidt to reduce numerical errors
        for(int lower = index_start; lower < i; lower++)
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
        double norm = sqrt(norm_sq * dx);
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
        energies[i] = final_energy * dx;

        for(int j = 0; j < N; j++)
        {
            wavefunctions[i, j] = current_eigenvector[j];
        }
    }
}

const double RadialSystem::radial_coordinate(int index)
{
    return exp(x_min + index * dx);
}
