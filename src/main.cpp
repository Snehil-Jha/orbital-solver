#include <iostream>
#include "dft.h"
using namespace std;

int main()
{
    int Z = 6;
    auto atom = DFT(Z, Z, 1e-10, 30, 1000);

    cout << "Initialized Atom: Z = " << Z << endl;

    double energy = 0;
    atom.compute_ground_state(energy, 0.01);
    cout << "Computed Energy: " << energy << endl;


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