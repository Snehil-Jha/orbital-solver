#include <cstdlib>
#include <iostream>
#include <fstream>
#include "dft.h"
using namespace std;

int main()
{
    ofstream File("data.csv");

    File << "atomic_number, charge, r_min, r_max, grid_points, mixing, energy, iterations, residual, occupation_s_up, occupation_s_down, occupation_p_up, occupation_p_down" << endl;

    cout << "Starting Solver" << endl;
    for(int Z = 1; Z <= 10; Z++)
    {
        for(int charge = 0; charge <= 1; charge++)
        {

            int ne = Z - charge;
            if(ne <= 0 || ne > 10) continue;

            auto atom = DFT(Z, ne, 1e-10, 30, 1000);
        
            cout << "Initialized Z = " << Z << ", Charge = " << charge << endl;
        
            double energy = 0;
            auto status = atom.compute_ground_state(energy, 0.05);
            int s_up, s_down, p_up, p_down;
            atom.get_occupation(s_up, s_down, p_up, p_down);

            File << Z << ", " << charge << ", " << "1e-10, 30, 1000, 0.05, " << energy << ", " << status.first << ", " << status.second <<
                 ", " << s_up << ", " << s_down << ", " << p_up << ", " << p_down << endl;
            cout << "Computed Energy: " << energy << " Hatrees in " << status.first << " iterations and a residual of " << status.second << endl;
        }
    }

    File.close();

    system("pause");

    return 0;

    // RadialSystem hydrogen(1e-10, 150, 10000, 0, [](double r){return - 1.0/r;});

    // Vector<double> energies(5);
    // Matrix<double> wave(5, 10000);

    // cout << "starting solver" << endl;
    // // hydrogen.solve_energy(energies);
    // hydrogen.solve_wavefunction(energies, wave);
    // cout << "ended solver" << endl;

    // for(int i = 0; i < 5; i++)
    // {
    //     cout << i << '\t' << energies[i] << endl;
    // }
}