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
#include <bits/stdc++.h>

using namespace std;


/// ------------------------
/// STD integration function
/// ------------------------
template<typename Method, typename F, typename Float>
double Std_integrate(F f, Float a, Float b, size_t steps,Method m)
{
    double s = 0;
    double h = (b-a)/steps;
    for (int i = 0; i < steps; ++i) s +=m(f, a + h*i, h);
    return h*s;
};

/// -------------------------------------------
/// IMPLEMENTED INTEGRATION METHODES
/// 1er -> 10th order
/// -------------------------------------------

class rectangular
{
public:
    enum position_type { left, middle, right };
    rectangular(position_type pos): position(pos) {}
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        switch(position)
        {
        case left:
            return f( x );
        case middle:
            return f( x+h/2.);
        case right:
            return f( x+h );
        }
    }
private:
    const position_type position;
};

class trapezium
{
public:
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        return (f( x ) + f( x+h ))/2.;
    }
};

class simpson
{
public:
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        return (f(x) + 4.*f(x+h/2) + f(x+h))/6.;
    }
};

class simpson_3_8
{
public:
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        return (f( x ) + 3.*f( x + h/3. ) + 3*f( x + 2.*h/3. ) + f( x + h ))/8.;
    }

};

class Bode
{
public:
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        return ( (14.)*f(x) + (64.)*f(x+h/4.) + (24.)*f(x+2.*h/4.) + (64.)*f(x+3*h/4.) + (14.)*f(x+h) )/(4.*45.);
    }

};

class Weddle
{
public:
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        return ( (41.)*f(x) + (216.)*f(x+h/6.) + (27.)*f(x+2.*h/6.) + (272.)*f(x+3*h/6.) + (27.)*f(x+4.*h/6.) + (216)*f(x+5.*h/6.) + 41.*f(x+h) )/(840.);
    }

};

class HighOrder
{
public:
    enum position_type { _7th, _8th, _9th, _10th };
    HighOrder(position_type pos): position(pos) {}
    template<typename F, typename Float>
    double operator()(F f, Float x, Float h) const
    {
        switch(position)
        {
        case _7th:
            return (1./17280.)*(751.*(f(x)+f(x+h)) + 3577.*(f(x+h/7.)+f(x+6.*h/7.)) + 1323.*(f(x+2.*h/7.)+f(x+5.*h/7.)) + 2989.*(f(x+3*h/7.)+f(x+4.*h/7)) );
        case _8th:
            return (1./28350.)*(989.*(f(x)+f(x+h)) + 5888.*(f(x+h/8.)+f(x+7.*h/8.)) - 928.*(f(x+2.*h/8.)+f(x+6.*h/8.)) + 10496.*(f(x+3.*h/8.)+f(x+5.*h/8.)) - 4540*f(x+4.*h/8.));
        case _9th:
            return (1./89600.)*(2857.*(f(x)+f(x+h)) + 15741*(f(x+h/9.)+f(x+8.*h/9.)) + 1080*(f(x+2.*h/9.)+f(x+7.*h/9.)) + 19344.*(f(x+3.*h/9.)+f(x+6.*h/9.)) + 5778*((f(x+4.*h/9.)+f(x+5*h/9)))  );
        case _10th:
            return (1./598752.)*(16067.*(f(x)+f(x+h)) + 106300*(f(x+h/10.)+f(x+9.*h/10)) - 48525.*(f(x+2*h/10.)+f(x+8*h/10)) + 272400.*(f(x+3.*h/10.)+f(x+7.*h/10.)) - 260550.*(f(x+4.*h/10.)+f(x+6.*h/10.)) + 427368*f(x+5.*h/10.)  )    ;
        }
    }
private:
    const position_type position;
};

/// ---------------------------------
/// Romberg integration
/// ---------------------------------
template<typename F, typename Float>
double Romberg_integrate (F f, Float a, Float b,size_t max_steps, float acc)
{
    double R1[max_steps], R2[max_steps];
    double *Rp = &R1[0];
    double *Rc = &R2[0];
    double h = (b-a);
    Rp[0] = (f(a) + f(b))*h*.5;
    for(size_t i = 1; i < max_steps; ++i)
    {
        h /= 2.;
        double c = 0;
        size_t ep = 1 << (i-1);
        for(size_t j = 1; j <= ep; ++j)
        {
            c += f(a+(2*j-1)*h);
        }
        Rc[0] = h*c + .5*Rp[0];

        for(size_t j = 1; j <= i; ++j)
        {
            double n_k = pow(4, j);
            Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1);
        }

        if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc)
        {
            return Rc[i-1];
        }
        double *rt = Rp;
        Rp = Rc;
        Rc = rt;
    }
    return Rp[max_steps-1];
};

/// ----------------------------
/// MontheCarlo integration nD;
/// ----------------------------
double MC_integrate(double(*f) (vector<double>), vector<double> a, vector<double> b, double prec,size_t max_steps)
{

    default_random_engine gene;
    uniform_real_distribution<double> distri(0,1);

    double I=0.0;
    double p=1.0;
    int m=1;

    vector<double> X(a.size());
    for(unsigned int i=0; i<a.size(); i++) p=p*(b[i]-a[i]);

    while(1)
    {
        m++;
        for(int i=0; i<a.size(); i++)
        {
            X[i]=a[i]+(b[i]-a[i])*distri(gene);   // génére les x,y,z dans les intervales a,b
        }
        I+=f(X)*p;
        if(m>=max_steps) break;
    }
    return I*p/m;
};

///------------------------------------------------------------
/// Double integration by Simpson 1/3 rule for DATA I_i(xi,yj)
///------------------------------------------------------------
double Simpson_2D(vector< vector<double> > Datas, double step_x, double step_y)
{
    int nx=Datas.size();
    int ny=Datas[0].size();
    //  cout<<nx<<" "<<ny<<endl;
    double h=step_x ;
    double k=step_y ;
    vector<double> ax ;
    ax.resize(nx);

    for (int i = 0; i < nx; ++i)
    {
        ax[i] = 0;
        for (int j = 0; j < ny; ++j)
        {
            if (j == 0 || j == ny - 1)
                ax[i] += Datas[i][j];
            else if (j % 2 == 0)
                ax[i] += 2. * Datas[i][j];
            else
                ax[i] += 4. * Datas[i][j];
        }
        ax[i] *= (k / 3.);
    }

    double answer = 0;

    for (int i = 0; i < nx; ++i)
    {
        if (i == 0 || i == nx - 1)
            answer += ax[i];
        else if (i % 2 == 0)
            answer += 2. * ax[i];
        else
            answer += 4. * ax[i];
    }
    answer *= (h / 3.);

    return answer;

}

/// --------------------------------------------------------------
/// Double integration by Simpson 1/3 rule for continuous f(x,y)
/// --------------------------------------------------------------

double Simpson_2D(double(*f)(double x, double y), double h, double k, double lx, double ux, double ly, double uy)
{
    int nx, ny;
    double answer;

    // number of points

    nx=(ux-lx)/h+1;
    ny=(uy-ly)/k+1;

    // Generate Grid of values

    vector< vector<double> > z(nx, vector<double>(ny));
    vector<double> ax ;
    ax.resize(nx);

    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            z[i][j]=( f(lx + i * h, ly + j * k) );
        }
    }

    for (int i = 0; i < nx; ++i)
    {
        ax[i] = 0;
        for (int j = 0; j < ny; ++j)
        {
            if (j == 0 || j == ny - 1)
                ax[i] += z[i][j];
            else if (j % 2 == 0)
                ax[i] += 2. * z[i][j];
            else
                ax[i] += 4. * z[i][j];
        }
        ax[i] *= (k / 3.);
    }

    answer = 0;

    for (int i = 0; i < nx; ++i)
    {
        if (i == 0 || i == nx - 1)
            answer += ax[i];
        else if (i % 2 == 0)
            answer += 2. * ax[i];
        else
            answer += 4. * ax[i];
    }
    answer *= (h / 3.);

    return answer;

};


#endif // INTEGRATION_H_INCLUDED
