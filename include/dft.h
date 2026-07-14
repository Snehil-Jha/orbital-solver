#ifndef DFT_H
#define DFT_H

#include "system.h"

#include<iostream>

class DFT {
    private:

    class OccupationState
    {
        public:
        int s_up;
        int s_down;
        int p_up;
        int p_down;

        Vector<double> energies_s_up;
        Vector<double> energies_s_down;
        Vector<double> energies_p_up;
        Vector<double> energies_p_down;

        Vector<double> energies;
        Matrix<double> wavefunctions;

        Vector<double> new_density_up;
        Vector<double> new_density_down;

        OccupationState(int iNe, int iN):
            energies_s_up(iNe), energies_s_down(iNe),
            energies_p_up(1 + ceil(iNe / 3.)), energies_p_down(1 + ceil(iNe / 3.)),
            energies(iNe), wavefunctions(iNe, iN),
            new_density_up(iN), new_density_down(iN)
        {}
    };

    const double xc_norm = - pow(1.5 / (std::numbers::pi * std::numbers::pi), 1. / 3.);

    const int N;
    const int Z;
    const int Ne;
    RadialSystem system;

    // variables for the poisson solver
    Vector<double> V_hatree;
    Vector<double> poisson_back;

    // densities
    Vector<double> density;
    Vector<double> density_up;
    Vector<double> density_down;
    
    // LDA
    Vector<double> V_xc_up;
    Vector<double> V_xc_down;
    

    // occupation states
    OccupationState occupation;

    /**
     * @brief Initializes the density to the independent hydrogenic case
     * 
     */
    void initialize_density();

    /**
     * @brief Solves the poisson equation for the current given density
     */
    void solve_poisson();

    /**
     * @brief Setups the exchange-correlation potential according to the current densities
     */
    void setup_xc_potential();

    /**
     * @brief Sets the internal potential for the system
     * 
     * @param l the azimuthal quantum number to be used
     * @param exchange_potential the exchange potential to be used 
     */
    void set_system_potential(int l, const Vector<double>& exchange_potential);

    /**
     * @brief Solves for the eigenvalues of the schrodinger equation
     */
    void solve_schrodinger();

    /**
     * @brief Computes the densities according to the solved energy levels 
     *
     * @param residual the difference norm between the old and new densities
     * @param mixing the mixing parameter between old and new densitites, where 1 means completely new
     */
    void compute_density(double& residual, const double mixing = 0.3);
    

    public:
    DFT(
        const int iZ,
        const int iNe,
        const double ir_min,
        const double ir_max,
        const int iN
    ): 
        Z(iZ), N(iN), Ne(iNe), system(ir_min, ir_max, iN, 0, [iZ](double r){return - iZ / r;}),
        V_hatree(N), poisson_back(N),
        density(N), density_up(N), density_down(N),
        V_xc_up(N), V_xc_down(N),
        occupation(iNe, iN)
    {}

    /**
     * @brief Performs the entire DFT algorithm to obtain the ground state energy
     * 
     * @param energy ground state energy as computed by the solver
     * @param max_iter the maximum number of iterations the solver will run for
     * @param residual_tol the tolerance factor for the difference norm between two iterations
     */
    void compute_ground_state(double &energy, const double mixing=0.3, const int max_iter = 8192, const double residual_tol = 1e-8);


    void print_status()
    {
        std::cout << occupation.s_up << std::endl;
        std::cout << occupation.s_down << std::endl;
        std::cout << occupation.p_up << std::endl;
        std::cout << occupation.p_down << std::endl;

    }
};

#endif