///  -----------------------------------------------------------------
///  Integral Class
///  Coded by JP Champeaux 2019
///
///  Implemented methods :
///  - rectangular
///  - Trapezium
///  - simpson
///  - simpson 3/8
///  - Bode
///  - Weddle
///  - HighOrder 7 8 9 10th order
///  - Romberg
///  - Simpson 2d (for data or continus fonction)
///  - nD Monte-Carlo integration
///  ----------------------------------------------------------------


#ifndef INTEGRATION_H_INCLUDED
#define INTEGRATION_H_INCLUDED

#define _USE_MATH_DEFINES

#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <vector>
#include <functional>
#include <algorithm>

/// ------------------------
/// STD integration function
/// ------------------------
template<typename Method, typename F, typename Float>
inline double Std_integrate(F f, Float a, Float b, size_t steps, const Method& m)
{
    double sum = 0.0;
    const double h = (b - a) / steps;
    for (size_t i = 0; i < steps; ++i)
        sum += m(f, a + h * i, h);
    return h * sum;
}

/// -------------------------------------------
/// IMPLEMENTED INTEGRATION METHODES
/// -------------------------------------------

class rectangular {
public:
    enum position_type { left, middle, right };
    explicit rectangular(position_type pos) : position(pos) {}

    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        switch (position) {
            case left: return f(x);
            case middle: return f(x + h / 2.0);
            case right: return f(x + h);
        }
        return 0.0; // fallback
    }
private:
    const position_type position;
};

class trapezium {
public:
    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        return (f(x) + f(x + h)) / 2.0;
    }
};

class simpson {
public:
    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        return (f(x) + 4.0 * f(x + h / 2.0) + f(x + h)) / 6.0;
    }
};

class simpson_3_8 {
public:
    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        return (f(x) + 3.0 * f(x + h / 3.0) + 3.0 * f(x + 2.0 * h / 3.0) + f(x + h)) / 8.0;
    }
};

class Bode {
public:
    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        return (14.0 * f(x) + 64.0 * f(x + h / 4.0) + 24.0 * f(x + h / 2.0) +
                64.0 * f(x + 3.0 * h / 4.0) + 14.0 * f(x + h)) / (4.0 * 45.0);
    }
};

class Weddle {
public:
    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        return (41.0 * f(x) + 216.0 * f(x + h / 6.0) + 27.0 * f(x + h / 3.0) +
                272.0 * f(x + h / 2.0) + 27.0 * f(x + 2.0 * h / 3.0) +
                216.0 * f(x + 5.0 * h / 6.0) + 41.0 * f(x + h)) / 840.0;
    }
};

class HighOrder {
public:
    enum position_type { _7th, _8th, _9th, _10th };
    explicit HighOrder(position_type pos) : position(pos) {}

    template<typename F, typename Float>
    inline double operator()(F f, Float x, Float h) const {
        switch (position) {
            case _7th:
                return (1.0 / 17280.0) * (751.0 * (f(x) + f(x + h)) +
                    3577.0 * (f(x + h / 7.0) + f(x + 6.0 * h / 7.0)) +
                    1323.0 * (f(x + 2.0 * h / 7.0) + f(x + 5.0 * h / 7.0)) +
                    2989.0 * (f(x + 3.0 * h / 7.0) + f(x + 4.0 * h / 7.0)));
            case _8th:
                return (1.0 / 28350.0) * (989.0 * (f(x) + f(x + h)) +
                    5888.0 * (f(x + h / 8.0) + f(x + 7.0 * h / 8.0)) -
                    928.0 * (f(x + 2.0 * h / 8.0) + f(x + 6.0 * h / 8.0)) +
                    10496.0 * (f(x + 3.0 * h / 8.0) + f(x + 5.0 * h / 8.0)) -
                    4540.0 * f(x + 4.0 * h / 8.0));
            case _9th:
                return (1.0 / 89600.0) * (2857.0 * (f(x) + f(x + h)) +
                    15741.0 * (f(x + h / 9.0) + f(x + 8.0 * h / 9.0)) +
                    1080.0 * (f(x + 2.0 * h / 9.0) + f(x + 7.0 * h / 9.0)) +
                    19344.0 * (f(x + 3.0 * h / 9.0) + f(x + 6.0 * h / 9.0)) +
                    5778.0 * (f(x + 4.0 * h / 9.0) + f(x + 5.0 * h / 9.0)));
            case _10th:
                return (1.0 / 598752.0) * (16067.0 * (f(x) + f(x + h)) +
                    106300.0 * (f(x + h / 10.0) + f(x + 9.0 * h / 10.0)) -
                    48525.0 * (f(x + 2.0 * h / 10.0) + f(x + 8.0 * h / 10.0)) +
                    272400.0 * (f(x + 3.0 * h / 10.0) + f(x + 7.0 * h / 10.0)) -
                    260550.0 * (f(x + 4.0 * h / 10.0) + f(x + 6.0 * h / 10.0)) +
                    427368.0 * f(x + 5.0 * h / 10.0));
        }
        return 0.0;
    }
private:
    const position_type position;
};

/// ---------------------------------
/// Romberg integration
/// ---------------------------------
template<typename F, typename Float>
double Romberg_integrate(F f, Float a, Float b, size_t max_steps, double acc)
{
    std::vector<double> R1(max_steps), R2(max_steps);
    double* Rp = R1.data();
    double* Rc = R2.data();
    double h = b - a;
    Rp[0] = 0.5 * h * (f(a) + f(b));

    for (size_t i = 1; i < max_steps; ++i)
    {
        h *= 0.5;
        double sum = 0.0;
        size_t ep = 1 << (i - 1);
        for (size_t j = 1; j <= ep; ++j)
            sum += f(a + (2 * j - 1) * h);

        Rc[0] = h * sum + 0.5 * Rp[0];
        for (size_t j = 1; j <= i; ++j)
        {
            double n_k = std::pow(4.0, j);
            Rc[j] = (n_k * Rc[j - 1] - Rp[j - 1]) / (n_k - 1.0);
        }

        if (i > 1 && std::abs(Rp[i - 1] - Rc[i]) < acc)
            return Rc[i];

        std::swap(Rp, Rc);
    }

    return Rp[max_steps - 1];
}

/// ----------------------------
/// Monte Carlo integration nD
/// ----------------------------
double MC_integrate(const std::function<double(const std::vector<double>&)>& f,
                    const std::vector<double>& a, const std::vector<double>& b,
                    double /*prec*/, size_t max_steps)
{
    std::default_random_engine gen;
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    const size_t dim = a.size();
    std::vector<double> X(dim);
    double volume = 1.0;

    for (size_t i = 0; i < dim; ++i)
        volume *= (b[i] - a[i]);

    double sum = 0.0;
    for (size_t m = 0; m < max_steps; ++m)
    {
        for (size_t i = 0; i < dim; ++i)
            X[i] = a[i] + (b[i] - a[i]) * dist(gen);
        sum += f(X);
    }

    return (sum / max_steps) * volume;
}

#endif // INTEGRATION_H_INCLUDED

