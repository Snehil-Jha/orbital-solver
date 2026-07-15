#include "dft.h"
#include <cmath>

#include <iostream>
#include <random>

void DFT::randomize_wavefunctions()
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    for(int j = 0; j < N; j++)
    {
        for(int i = 0; i < Ne; i++)
        {
            occupation.wavefunctions_s_up[i, j] = dis(gen);
            occupation.wavefunctions_s_down[i, j] = dis(gen);
        }
        for(int i = 0; i < 1 + ceil(Ne / 3.); i++)
        {
            occupation.wavefunctions_p_up[i, j] = dis(gen);
            occupation.wavefunctions_p_down[i, j] = dis(gen);
        }
    }
}

void DFT::initialize_density()
{
    // start with an initial guess of hydrogenic density, scaled properly for a logarithmic grid
    const double density_norm = 4 * Ne * Z * Z * Z;
    for(int i = 0; i < N; i++)
    {
        density[i] = density_norm * system.ri_sq[i] * system.ri[i] * exp(- 2 * Z * system.ri[i]);

        // slight bias to up spin, to induce asymmetry
        density_up[i] = 0.515625 * density[i];
        density_down[i] = 0.484375 * density[i];
    }
}

void DFT::solve_poisson()
{
    // solves the poisson equation by gauss law and applying Adams-Moulton to both integrals individually
    // implicit assumption that number of grid points is larger than 2

    constexpr double coeff = 1. / 24.;

    // forward pass, using trapezoidal rule to kickstart
    V_hatree[0] = 0;
    V_hatree[1] = 0.5 * system.dx * (density[0] + density[1]);
    V_hatree[2] = V_hatree[1] + 0.5 * system.dx * (density[1] + density[2]);

    // backward pass
    poisson_back[N - 1] = 0;
    poisson_back[N - 2] = 0.5 * system.dx * ((density[N - 1] / system.ri[N - 1]) + (density[N - 2] / system.ri[N - 2]));
    poisson_back[N - 3] = poisson_back[N - 2] + 0.5 * system.dx * ((density[N - 2] / system.ri[N - 2]) + (density[N - 3] / system.ri[N - 3]));

    for(int i = 3; i < N; i++)
    {
        V_hatree[i] = V_hatree[i - 1] + coeff * system.dx * (9 * density[i] + 19 * density[i - 1] - 5 * density[i - 2] + density[i - 3]);
        
        poisson_back[N - 1 - i] = poisson_back[N - i] + coeff * system.dx * (
            (density[N - 1 - i] / system.ri[N - 1 - i]) *   9 +
            (density[N     - i] / system.ri[N     - i]) *  19 +
            (density[N + 1 - i] / system.ri[N + 1 - i]) * (-5) +
            (density[N + 2 - i] / system.ri[N + 2 - i]) *   1
        );
    }

    // full sweep
    for(int i = 0; i < N; i++)
    {
        V_hatree[i] = (V_hatree[i] / system.ri[i]) + poisson_back[i];
    }
}

void DFT::setup_xc_potential()
{
    // note that this is not the actual LSDA potential, but rather a scaled version which will be fixed later
    for(int i = 0; i < N; i++)
    {
        if(density_up[i] < 0)
        {
            std::cout << "negative density_up[" << i << "]" << std::endl;
        }
        if(density_down[i] < 0)
        {
            std::cout << "negative density_down[" << i << "]" << std::endl;
        }

        V_xc_up[i] = xc_norm * pow(std::max(density_up[i],   0.0), 1. / 3.);
        V_xc_down[i] = xc_norm * pow(std::max(density_down[i],   0.0), 1. / 3.);
    }   
}

void DFT::set_system_potential(int l, const Vector<double>& exchange_potential)
{
    system.l = l;
    for(int i = 0; i < N; i++)
    {
        system.ham_main[i] = 
            system.kinetic_main[i] +
            system.pot_main[i] +
            0.5 * l * (l + 1) +
            system.ri_sq[i] * V_hatree[i] +
            system.ri[i] * exchange_potential[i];
    }
}

void DFT::solve_schrodinger()
{
    // computes all the possible energy states
    set_system_potential(0, V_xc_up);
    system.solve_energy(occupation.energies_s_up);

    set_system_potential(0, V_xc_down);
    system.solve_energy(occupation.energies_s_down);

    set_system_potential(1, V_xc_up);
    system.solve_energy(occupation.energies_p_up);
    
    set_system_potential(1, V_xc_down);
    system.solve_energy(occupation.energies_p_down);

    // resets the occupation
    occupation.s_up = occupation.s_down = occupation.p_up = occupation.p_down = 0;

    // finds the number of electrons in each set of orbitals
    while(
        occupation.s_up + occupation.s_down + occupation.p_up + occupation.p_down < Ne
    )
    {
        double min_energy = occupation.energies_s_up[occupation.s_up];
        int choice = 0;

        if(occupation.energies_s_down[occupation.s_down] < min_energy)
        {
            min_energy = occupation.energies_s_down[occupation.s_down];
            choice = 1;
        }

        // division by three due to p-orbital degeneracy
        if(occupation.energies_p_up[occupation.p_up / 3] < min_energy)
        {
            min_energy = occupation.energies_p_up[occupation.p_up / 3];
            choice = 2;
        }

        if(occupation.energies_p_down[occupation.p_down / 3] < min_energy)
        {
            min_energy = occupation.energies_p_down[occupation.p_down / 3];
            choice = 3;
        }

        if(choice == 0) occupation.s_up++;
        else if(choice == 1) occupation.s_down++;
        else if(choice == 2) occupation.p_up++;
        else if(choice == 3) occupation.p_down++;
    }
}

void DFT::compute_density(double& residual, const double mixing)
{
    auto add_from_state = [this](
        const Vector<double>& exchange_potential,
        Vector<double>& energies, Matrix<double>& wavefunctions, 
        int occupation_count, int degeneracy, int l,
        Vector<double>& updating_density
    ){

        if(occupation_count == 0) return;

        int maximum_state = ceil((1. * occupation_count) / degeneracy);

        set_system_potential(l, exchange_potential);
        system.solve_wavefunction(energies, wavefunctions, 0, maximum_state, true);

        for(int state = 0; state < maximum_state; state++)
        {
            for(int i = 0; i < N; i++)
            {
                int specific_occupation = (state != maximum_state - 1 || occupation_count % degeneracy == 0) ? degeneracy : (occupation_count % degeneracy);
                updating_density[i] += specific_occupation * system.ri_sq[i] * (wavefunctions[state, i] * wavefunctions[state, i]);
            }
        }
    };

    for(int i = 0; i < N; i++)
    {
        occupation.new_density_up[i] = occupation.new_density_down[i] = 0;
    }


    add_from_state(V_xc_up, occupation.energies_s_up, occupation.wavefunctions_s_up, occupation.s_up, 1, 0, occupation.new_density_up);
    add_from_state(V_xc_down, occupation.energies_s_down, occupation.wavefunctions_s_down, occupation.s_down, 1, 0, occupation.new_density_down);
    
    add_from_state(V_xc_up, occupation.energies_p_up, occupation.wavefunctions_p_up, occupation.p_up, 3, 1, occupation.new_density_up);
    add_from_state(V_xc_down, occupation.energies_p_down, occupation.wavefunctions_p_down, occupation.p_down, 3, 1, occupation.new_density_down);

    residual = 0;
    for(int i = 0; i < N; i++)
    {
        double new_density_up = occupation.new_density_up[i] * mixing + density_up[i] * (1 - mixing);
        double new_density_down = occupation.new_density_down[i] * mixing + density_down[i] * (1 - mixing);

        residual += (new_density_up - density_up[i]) * (new_density_up - density_up[i]);
        residual += (new_density_down - density_down[i]) * (new_density_down - density_down[i]);

        density_up[i] = new_density_up;
        density_down[i] = new_density_down;

        density[i] = density_up[i] + density_down[i];
    }

    residual = sqrt(residual);
}

void DFT::compute_ground_state(double& energy, const double mixing, const int max_iter, const double residual_tol)
{
    double residual = 0;

    // TODO:
    randomize_wavefunctions();

    initialize_density();
    for(int iter_number = 0; iter_number < max_iter; iter_number++)
    {
        solve_poisson();
        
        setup_xc_potential();

        solve_schrodinger();

        compute_density(residual, mixing);
        
        std::cout << "iteration: " << iter_number << "\t residual: " << residual
           << "\t occ: " << occupation.s_up << "," << occupation.s_down
           << "," << occupation.p_up << "," << occupation.p_down << std::endl;


        if(residual < residual_tol)
            break;
    }

    // compute the final energy
    
    double E_eig = 0;
    for(int i = 0; i < occupation.s_up; i++)
    {
        E_eig += occupation.energies_s_up[i];
    }
    for(int i = 0; i < occupation.s_down; i++)
    {
        E_eig += occupation.energies_s_down[i];
    }
    for(int i = 0; i < occupation.p_up; i++)
    {
        E_eig += occupation.energies_p_up[i / 3];
    }
    for(int i = 0; i < occupation.p_down; i++)
    {
        E_eig += occupation.energies_p_down[i / 3];
    }

    double E_hartree_correction = 0;
    double E_xc_correction = 0;

    // trapezoidal integration
    for(int i = 0; i < N; i++)
    {
        double weight = system.dx;
        if (i == 0 || i == N - 1) weight *= 0.5;

        E_hartree_correction -= 0.5 * V_hatree[i] * density[i] * weight;

        double V_xc_phys_up   = V_xc_up[i] / system.ri[i];
        double V_xc_phys_down = V_xc_down[i] / system.ri[i];

        E_xc_correction -= 0.25 * (V_xc_phys_up * density_up[i] + V_xc_phys_down * density_down[i]) * weight;
    }

    energy = E_eig + E_hartree_correction + E_xc_correction;
}