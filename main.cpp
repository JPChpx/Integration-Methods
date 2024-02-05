///  ------------------------------------------------------
///  Integration Class
///  Coded by J.P Champeaux 2019
///  stand alone integration class : Integration.h
///  For educational purpose
///  Main : Exemples of use for Pi calculation
///  ------------------------------------------------------

#include <iostream>
#include <iostream>
#include "version.h"
#include "math.h"
#include "Integration.h"

using namespace std;

// Example of Function to integrate
double f(double x)
{
    return 4./(1.+x*x);   // fonction 1D
}


// Example of Data function to integrate
double fc(vector<double> X) // fonciton nD
{
    double y=0;
    for (int j = 0; j < X.size(); j = j+1)
    {
        y=y+4.*(1./(1.+X[j]*X[j]));
    }
    return y;
}

// Example on 2d function to integrate
double fxy(double x, double y)
{
    return exp(-(x*x+y*y));
}


int main()
{
    cout<<" -----------------------------------------------------------"<<endl;
    cout<<" Integration Class by JPC  "<<endl;
    cout<<" Vers :"<<AutoVersion::FULLVERSION_STRING<<endl;
    cout<<" Rev  :"<<AutoVersion::REVISION<<"."<<AutoVersion::BUILD<<endl;
    cout<<" Date :"<<AutoVersion::DATE<<"/"<<AutoVersion::MONTH<<"/"<<AutoVersion::YEAR<<endl;
    cout<<" -----------------------------------------------------------"<<endl<<endl;

    cout.precision(5);
    double pi=3.141592653589793;
    cout<<"PI ="<<pi<<endl;

    // rectangular
    double v=Std_integrate (f, 0., 1., 1000,rectangular(rectangular::left));
    double err=fabs(pi-v);
    cout<<"PI rectangular left     ="<<v<<" ("<<err<<")"<<  endl;
    // rectangular
    v=   Std_integrate (f, 0., 1., 1000,rectangular(rectangular::right));
    err=fabs(pi-v);
    cout<<"PI rectangular right    ="<<v<<" ("<<err<<")"<<  endl;
    // rectangular
    v=   Std_integrate (f, 0., 1., 1000,rectangular(rectangular::middle));
    err=fabs(pi-v);
    cout<<"PI rectangular center   ="<<v<<" ("<<err<<")"<<  endl;

    // Trapezium
    v=   Std_integrate (f, 0., 1., 1000,trapezium());
    err=fabs(pi-v);
    cout<<"PI trapeze              ="<<v<<" ("<<err<<")"<<  endl;

    // Simpson
    v=Std_integrate (f, 0., 1., 1000,simpson());
    err=fabs(pi-v);
    cout<<"PI simpson              ="<<v<<" ("<<err<<")"<<  endl;

    // Simpson 3/8
    v=Std_integrate (f, 0., 1., 1000,simpson_3_8());
    err=fabs(pi-v);
    cout<<"PI simpson 1/3          ="<<v<<" ("<<err<<")"<<  endl;

    // Bode
    v=Std_integrate (f, 0., 1., 1000,Bode());
    err=fabs(pi-v);
    cout<<"PI Bode                 ="<<v<<" ("<<err<<")"<<  endl;

    // Weddle
    v=Std_integrate (f, 0., 1., 1000,Weddle());
    err=fabs(pi-v);
    cout<<"PI Weddle               ="<<v<<" ("<<err<<")"<<  endl;

    // 7th Order
    v=Std_integrate (f, 0., 1., 1000,HighOrder(HighOrder::_7th));
    err=fabs(pi-v);
    cout<<"PI _7th                 ="<<v<<" ("<<err<<")"<<  endl;

    // 8th Order
    v=Std_integrate (f, 0., 1., 1000,HighOrder(HighOrder::_8th));
    err=fabs(pi-v);
    cout<<"PI _8th                 ="<<v<<" ("<<err<<")"<<  endl;

    // 9th order
    v=Std_integrate (f, 0., 1., 1000,HighOrder(HighOrder::_9th));
    err=fabs(pi-v);
    cout<<"PI _9th                 ="<<v<<" ("<<err<<")"<<  endl;

    // 10th order
    v=Std_integrate (f, 0., 1., 1000,HighOrder(HighOrder::_10th));
    err=fabs(pi-v);
    cout<<"PI _10th                ="<<v<<" ("<<err<<")"<<  endl;

    // Romberg
    v=Romberg_integrate(f,0.,1.,1000,1.e-6);
    err=fabs(pi-v);
    cout<<"PI Romberg              ="<<v<<" ("<<err<<")"<<  endl;

    // Monte-carlo integration
    vector<double> A = {0.0};//, 0.0, 0.0, 0.0, 0.0, 0.0}; /* left end-points */
    vector<double> B = {1.0};// 1.0, 1.0, 1.0, 1.0, 1.0}; /* right end-points */
    v=MC_integrate(fc,A,B,1.e-8,1000);
    err=fabs(pi-v);
    cout<<"PI MC                   ="<<v<<" ("<<err<<")"<<  endl;

    // Integration 2D simpson 1/3

    float h, k, lx, ux, ly, uy;

    lx = -100., ux = 100., ly = -100,
    uy = 100, h = 0.1, k = 0.2;
    v=Simpson_2D(fxy, h, k, lx, ux, ly, uy);
    err=fabs(pi-v);
    cout<<"PI 2dSimpson 1/3        ="<<v<<" ("<<err<<")"<<  endl;


    // Integration 2D Datas smpson 1/3
    int nx=100;
    int ny=100;
    h=0.1;
    k=0.2;
    // building 2d datas
    vector< vector<double> > Datas(nx, vector<double>(ny));
    double x,y;
    for(int i=0; i<nx; i++)
    {
        x=i*h ;
        for(int j=0; j<ny; j++)
        {
            y=j*k;
            Datas[i][j]=fxy(x,y);
        }

    }
    // integration of 2d datas
    v=4.*Simpson_2D(Datas,h,k); // 4 car intégré sur x=0->10 y=0->5;
    err=fabs(pi-v);
    cout<<"2d GIRD 100x100 DATAs   ="<<v<<" ("<<err<<")"<<  endl;

    return 0;
}
