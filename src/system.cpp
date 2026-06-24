#include "system.h"

#include "vector.h"
#include "matrix.h"
#include "matrix_algo.h"
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

}