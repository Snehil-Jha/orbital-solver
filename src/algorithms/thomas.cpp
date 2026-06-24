#include "matrix_algo.h"

bool symmetric_thomas_solver(
    const Vector<double> &main,
    const Vector<double> &sub,
    const Vector<double> &b,
    Vector<double>& x,
    Vector<double>& tmp,
    const double epsilon
)
{
    bool blow_up = false;

    int n = main.length();
    if(sub.length() != n - 1 || b.length() != n || x.length() != n || tmp.length() != n - 1)
        throw std::invalid_argument("array size mismatch");


    double bet = main[0];
    if (std::abs(bet) < epsilon) { 
        bet = (bet < 0 ? -epsilon : epsilon);
        blow_up = true; 
    }

    tmp[0] = sub[0] / bet;
    x[0] = b[0] / bet;

    for(int j = 1; j < n - 1; j++)
    {
        bet = main[j] - sub[j - 1] * tmp[j - 1];
        if (std::abs(bet) < epsilon) { 
            bet = (bet < 0 ? -epsilon : epsilon);
            blow_up = true; 
        }

        tmp[j] = sub[j] / bet;
        x[j] = (b[j] - sub[j - 1] * x[j-1]) / bet;
    }

    bet = main[n-1] - sub[n-2] * tmp[n-2];
    if (std::abs(bet) < epsilon) { 
        bet = (bet < 0 ? -epsilon : epsilon); 
        blow_up = true; 
    }
    x[n-1] = (b[n-1] - sub[n-2] * x[n-2]) / bet;


    for(int j = n - 2; j >= 0; j--)
    {
        x[j] -= tmp[j] * x[j+1];
    }

    return blow_up;
}
