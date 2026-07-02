#include <iostream>
#include <fstream>
#include <ostream>

#include "vector.h"
#include "system.h"

using namespace std;

int main()
{
    int N = 10000;
    int m = 4;
    int k = 2;

    auto system = RadialSystem(
        1e-12, 150., N, 0, [](double r){
            return -1. / r;
        }
    );

    auto energies = Vector<double>(m);
    auto wavefunctions = Matrix<double>(m, N);
    system.solve_wavefunction(energies, wavefunctions);

    for(int i = 0; i < m; i++)
    {
        cout << "E = " << energies[i] << '\t';
        cout << 1. / (2 * (i + 1) * (i + 1)) << endl;
    }

    ofstream outputFile("data.csv");

    if (outputFile.is_open()) {
        outputFile << "r,psi\n";

        for (int i = 0; i < N; i++) {
            outputFile << system.radial_coordinate(i) << "," << wavefunctions[k, i] << "\n";
        }

        outputFile.close();
        cout << "Data successfully exported to data.csv" << endl;
    } else {
        cerr << "Failed to open file for writing." << endl;
    }
}