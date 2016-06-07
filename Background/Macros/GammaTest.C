#include <iostream>
#include <iomanip>
#include <TMath.h>

using namespace std;

void GammaTest() {

    double n = 2;
    double x = 0;

    for (unsigned i=1;i<5;i++) {

//        cout << setw(6) << i << setw(12) << inc_gamma(1+i,-x) << endl;
//        cout << setw(6) << i << setw(12) << integral(i,-x) << endl;
        double temp = integral(i,1);
        cout << setw(6) << i << setw(12) << temp << endl; 

    }

    cout << fact(25)*pow(-1,25) << endl;
    
}

double fact(unsigned &n) {

    double out(1);
    for (unsigned i(0);i<n;i++) out *= n-i;
    return out;

}

double expoSum(unsigned n, double x) {
    double out;
    for (unsigned k=0;k<=n;k++) {
        out += pow(x,k)/fact(k);
    }
    return out;
}

double inc_gamma(unsigned n, double x) {

    cout << "inc_gamma input:  " << setw(12) << n << setw(12) << x << endl;
    cout << "inc_gamma output: " << setw(12) << fact(n-1)*exp(-x)*expoSum(n-1,x) << endl;

    return fact(n-1)*exp(-x)*expoSum(n-1,x);

}

double integral(unsigned n, double x) {

    double incGamma = inc_gamma(n+1,-log(x));
    double factor = pow(-1.0,n);

    cout << setw(24) << factor << setw(12) << incGamma << endl;

    return factor*incGamma;

}


