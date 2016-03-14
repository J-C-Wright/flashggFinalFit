
#include <vector>
#include <iostream>
#include <math>

using namespace std;

void HyperGeoTest() {

    gSystem->Load("libMathMore");

    vector<double> parameters;
    parameters.push_back(1);
    parameters.push_back(1);
    parameters.push_back(1);
    parameters.push_back(1);
    parameters.push_back(1);

    double x(0.9);

    double polyLog(0.0);
    for (unsigned i(1);i<parameters.size();i++) {

        if (i == 1) {
            polyLog += parameters[i];
        }else if (i == 2) {
            polyLog += parameters[i]*log(x);
        }else{
            polyLog += parameters[i]*pow(log(x),i);
        }
    }

    float integral = (pow(x,1+polyLog) * ROOT::Math::hyperg(-parameters[0],3*(1+polyLog),3*(1+polyLog)+1,pow(x,1.0/3.0)))/(1+polyLog);

    cout << "Integral is " << integral << endl;

}

