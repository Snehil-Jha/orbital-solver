#include "matrix_algo.h"
#include <limits>

void bisect(
    const Vector<double>& main,
    const Vector<double>& sub,
    int m1, int m2,
    Vector<double>& result,
    const double epsilon, const int max_iter
)
{
    // verification
    int n = main.length();
    if(sub.length() != n - 1 || result.length() != m2 - m1 + 1)
        throw std::invalid_argument("vectors of unexpected sizes obtained");


    // calculation of bounds of eigenvalues
    double xmin = main[0] - std::abs(sub[0]);
    double xmax = main[0] + std::abs(sub[0]);
    double h;
    for(int i = 1; i < n; i++)
    {
        h = std::abs(sub[i - 1]);
        if(i != n - 1)
        {
            h += std::abs(sub[i]);
        }

        xmax = std::max(xmax, main[i] + h);
        xmin = std::min(xmin, main[i] - h);
    }

    // start the bisection method

    // create array of lower bounds. (upper bounds go inside results)
    auto lower = Vector<double>(result.length(), xmin);
    for(int i = 0; i < result.length(); i++) result[i] = xmax;

    double x_upper = xmax;

    // loop for each eigenvalue
    for(int k = m2 - m1; k >= 0; k--)
    {
        double x_lower = xmin;
        for(int i = k; i >= 0; i--)
        {
            if(x_lower < lower[i])
            {
                x_lower = lower[i];
                break;
            }
        }

        double x_mid;
        for(int iter_count = 0; iter_count < max_iter; iter_count++)
        {
            x_mid = 0.5 * (x_upper + x_lower);
            if(x_upper - x_lower < epsilon * std::max(1., std::abs(x_mid))) break;

            // sturms sequence
            int a = 0;
            double q = 1;
            for(int i = 0; i < n; i++)
            {
                double sub_element = (i == 0) ? 0 : sub[i-1];

                if(q == 0) q = main[i] - x_mid - std::abs(sub_element) / std::numeric_limits<double>::epsilon();
                else q = main[i] - x_mid - sub_element * sub_element / q;

                if(q < 0) a++;
            }

            if(a <= k + m1)
            {
                if(a < m1)
                    x_lower = lower[0] = x_mid;
                else
                    x_lower = lower[a - m1] = x_mid;
            }
            else
            {
                x_upper = x_mid;
            }
        }

        result[k] = 0.5 * (x_lower + x_upper);
    }
}


void bisect(
    const Vector<double>& main,
    const Vector<double>& sub,
    const Vector<double>& weight,
    int m1, int m2,
    Vector<double>& result,
    const double epsilon, const int max_iter
)
{
    // verification
    int n = main.length();
    if(sub.length() != n - 1 || result.length() != m2 - m1 + 1)
        throw std::invalid_argument("vectors of unexpected sizes obtained");


    // calculation of bounds of eigenvalues
    double xmin = main[0] - std::abs(sub[0]);
    double xmax = main[0] + std::abs(sub[0]);
    double h;
    for(int i = 1; i < n; i++)
    {
        h = std::abs(sub[i - 1]);
        if(i != n - 1)
        {
            h += std::abs(sub[i]);
        }

        xmax = std::max(xmax, main[i] + h);
        xmin = std::min(xmin, main[i] - h);
    }

    // start the bisection method

    // create array of lower bounds. (upper bounds go inside results)
    auto lower = Vector<double>(result.length(), xmin);
    for(int i = 0; i < result.length(); i++) result[i] = xmax;

    double x_upper = xmax;

    // loop for each eigenvalue
    for(int k = m2 - m1; k >= 0; k--)
    {
        double x_lower = xmin;
        for(int i = k; i >= 0; i--)
        {
            if(x_lower < lower[i])
            {
                x_lower = lower[i];
                break;
            }
        }

        double x_mid;
        for(int iter_count = 0; iter_count < max_iter; iter_count++)
        {
            x_mid = 0.5 * (x_upper + x_lower);
            if(x_upper - x_lower < epsilon * std::max(1., std::abs(x_mid))) break;

            // sturms sequence
            int a = 0;
            double q = 1;
            for(int i = 0; i < n; i++)
            {
                double sub_element = (i == 0) ? 0 : sub[i-1];

                if(q == 0) q = main[i] - (x_mid * weight[i]) - std::abs(sub_element) / std::numeric_limits<double>::epsilon();
                else q = main[i] - (x_mid * weight[i]) - sub_element * sub_element / q;

                if(q < 0) a++;
            }

            if(a <= k + m1)
            {
                if(a < m1)
                    x_lower = lower[0] = x_mid;
                else
                    x_lower = lower[a - m1] = x_mid;
            }
            else
            {
                x_upper = x_mid;
            }
        }

        result[k] = 0.5 * (x_lower + x_upper);
    }
}