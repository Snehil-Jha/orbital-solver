#include <cmath>
#include <iostream>
#include "vector.h"
#include "matrix.h"
#include "matrix_solvers.h"

using namespace std;

#define dx 0.0001220703125
#define dx_inv 8192.
#define x0 -4.
#define N 65536

#define m 50

void H(const Matrix<std::complex<double>>& state, int col, Vector<std::complex<double>>& out)
{
    double x;
    std::complex<double> potential, kinetic;

    for(int i = 0; i < N; i++)
    {
        x = x0 + i * dx;

        // square well
        potential = (abs(x) < 0.5) ? 0. : 100;

        if(i == 0)
        {
            kinetic = -0.5 * dx_inv * dx_inv * (state[i + 2, col] - 2. * state[i + 1, col] + state[i, col]);
        }
        else if(i == N - 1)
        {
            kinetic = -0.5 * dx_inv * dx_inv * (state[i, col] - 2. * state[i - 1, col] + state[i - 2, col]);
        }
        else
        {
            kinetic = -0.5 * dx_inv * dx_inv * (state[i + 1, col] + state[i - 1, col] - 2. * state[i, col]);
        }

        out[i] = potential + kinetic;
    }
}

int main()
{
    cout << "Starting sim" << endl;

    auto T_main = Vector<double>(m);
    auto T_sub = Vector<double>(m - 1);
    auto Q = Matrix<std::complex<double>>(N, m);

    cout << "Initialized values" << endl;

    lanczos_get_tridiagonal(H, T_main, T_sub, Q);

    cout << "Done lanczos." << endl;
    
    auto result = Vector<double>(m);

    bisect(T_main, T_sub, 0, m-1, result);

    for(int i = 0; i < m; i++)
    {
        cout << result[i] << endl;
    }
}