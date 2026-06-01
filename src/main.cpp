#include <iostream>
#include "vector.h"
#include "matrix.h"
#include "matrix_solvers.h"


using namespace std;

#define dx 0.015625
#define dx_inv 64.
#define x0 -16.
#define N 2048

#define m 50

void H(const Vector<double>& in, Vector<double>& out)
{
    double x, pot, kin;
    for(int i = 0; i < N; i++)
    {
        x = x0 + i * dx;

        // harmonic potential
        pot = 0.5 * x * x;

        // kinetic
        if(i == 0)
            kin = -0.5 * dx_inv * dx_inv * (in[i + 1] - 2 * in[i]);
        else if (i == N - 1)
            kin = -0.5 * dx_inv * dx_inv * (- 2 * in[i] + in[i - 1]);
        else
            kin = -0.5 * dx_inv * dx_inv * (in[i + 1] - 2 * in[i] + in[i - 1]);

        out[i] = pot * in[i] + kin;
    }
}

int main()
{
    auto T_main = Vector<double>(m);
    auto T_sub = Vector<double>(m - 1);

    auto Q = Matrix<double>(N, m);

    int k = 10;
    auto result = Vector<double>(k);

    double shift = 0.49;

    cout << "Starting lanczos" << endl;
    lanczos_invert_and_shift(H, T_main, T_sub, Q, shift);

    cout << "Starting bisection" << endl;
    bisect(T_main, T_sub, m-k, m-1, result);

    cout << "Printing results" << endl;

    for(int i = k-1; i >= 0; i--)
    {
        cout << "n = " << k-i-1 << "\t\tE=" << (shift + 1./result[i]) << endl;
    }

    cout << "Finished" << endl;
}