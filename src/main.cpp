#include <iostream>
#include "vector.h"
#include "system.h"

using namespace std;

// #define dx 0.125
// #define dx_inv 8.
// #define x0 -4.
// #define N 65

// #define m 50

// void H(const Vector<double>& in, Vector<double>& out)
// {
//     double x, pot, kin;
//     for(int i = 0; i < N; i++)
//     {
//         x = x0 + i * dx;

//         // harmonic potential
//         pot = 0.5 * x * x;

//         // kinetic
//         if(i == 0)
//             kin = -0.5 * dx_inv * dx_inv * (in[i + 1] - 2 * in[i]);
//         else if (i == N - 1)
//             kin = -0.5 * dx_inv * dx_inv * (- 2 * in[i] + in[i - 1]);
//         else
//             kin = -0.5 * dx_inv * dx_inv * (in[i + 1] - 2 * in[i] + in[i - 1]);

//         out[i] = pot * in[i] + kin;
//     }
// }

// int main()
// {
//     auto T_main = Vector<double>(m);
//     auto T_sub = Vector<double>(m - 1);

//     auto Q = Matrix<double>(N, m);

//     int k = 5;
//     auto result = Vector<double>(k);

//     double shift = 0;

//     cout << "Starting lanczos" << endl;
//     lanczos_invert_and_shift(H, T_main, T_sub, Q, shift);

//     cout << "Starting bisection" << endl;
//     bisect(T_main, T_sub, m-k, m-1, result);

//     cout << "Printing results" << endl;

//     for(int i = k-1; i >= 0; i--)
//     {
//         cout << "n = " << k-i-1 << "\t\tE=" << (shift + 1./result[i]) << endl;
//     }
    
//     cout << "Starting RQI" << endl;
//     auto eigenstate = Vector<double>(m);
//     double eigenvalue = 0;
//     cout << "Finished RQI in " << symmetric_tridag_rqi(result[2], T_main, T_sub, eigenvalue, eigenstate) << " steps" << endl;

//     cout << "Updated eigenvalue: " << (shift + 1./eigenvalue) << endl;

//     auto position_eigenstate = Vector<double>(N);
//     for(int i = 0; i < m; i++)
//     {
//         for(int j = 0; j < N; j++)
//         {
//             position_eigenstate[j] += eigenstate[i] * Q[j, i];
//         }
//     }

//     cout << "printing final wavefunction: " << endl;
//     cout << fixed;
//     for(int i = 0; i < N; i++)
//     {
//         cout << position_eigenstate[i] << ", ";
//     }
    

//     cout << endl << "Finished" << endl;
// }

int main()
{
    auto system = RadialSystem(
        1e-12, 120., 10000, 0, [](double r){
            return -1. / r;
        }
    );

    auto energies = Vector<double>(5);
    system.solve_energy(energies);

    for(int i = 0; i < 5; i++)
    {
        cout << "E = " << energies[i] << endl;
    }
}