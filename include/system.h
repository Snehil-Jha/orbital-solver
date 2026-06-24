#ifndef SYSTEM_H
#define SYSTEM_H

#include <cmath>
#include <functional>

#include "vector.h"
#include "matrix.h"

class RadialSystem {

    private:
    const double r_min;
    const double r_max;
    const int N;

    const double x_min;
    const double x_max;
    const double dx;
    const double exp_dx;

    double l;

    Vector<double> ri;
    Vector<double> ri_sq;

    Vector<double> kinetic_main;
    Vector<double> ham_sub;

    Vector<double> ham_main;


    public:
    template <typename F>
    RadialSystem(
        const double ir_min,
        const double ir_max,
        const int iN,
        const int il,
        const F&& V
    ): 
        r_min(ir_min), r_max(ir_max), N(iN),
        x_min(log(ir_min)), x_max(log(ir_max)), 
        dx( (x_max - x_min) / (N - 1) ), exp_dx(exp(dx)), l(il),
        ri(N), ri_sq(N), kinetic_main(N), ham_sub(N - 1), ham_main(N)
    {

        double xi;
        int i;
        for(i = 0; i < N - 1; i++)
        {
            xi = x_min + i * dx;
            ri[i] = exp(xi);
            ri_sq[i] = ri[i] * ri[i];

            kinetic_main[i] = 0.125 + (1 / (dx * dx));
            ham_main[i] = kinetic_main[i] + ri_sq[i] * V(ri[i]) + 0.5 * l * (l + 1);
            ham_sub[i] = -0.5 / (dx * dx);
        }

        xi = x_min + i * dx;
        ri[i] = exp(xi);
        ri_sq[i] = ri[i] * ri[i];

        kinetic_main[i] = 0.125 + (1 / (dx * dx));
        ham_main[i] = kinetic_main[i] + ri_sq[i] * V(ri[i]) + 0.5 * l * (l + 1);
    }


    /**
     * @brief Changes the potential for this solver, to avoid having to create new objects in a SCF loop
     * 
     * @param updated_V the new radial potential function
     * @param updated_l the new azimuthal number, use -1 for not changing
     */
    template <typename F>
    void update_potential(const F&& updated_V, int updated_l = -1);

    
    /**
     * @brief Computes the first few eigenvalues of the system
     * 
     * @param energies the energies calcualted by the solver
     */
    void solve_energy(Vector<double> &energies);
    
    /**
     * @brief Computes the first few eigenstates of the system, the eigenvalues will be computed as well
     * 
     * @param energies the energies calcualted by the solver
     * @param wavefunctions the wavefunctions calculated by the solver
     */
    void solve_wavefunction(Vector<double> &energies, Matrix<double> &wavefunctions);
};


#endif