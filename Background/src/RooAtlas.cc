
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"

#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "../interface/RooAtlas.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooArgList.h"

#include "Math/SpecFuncMathMore.h"

using namespace std;

RooAtlas::RooAtlas(const char *name, const char *title, RooAbsReal &x, const RooArgList &paramList, unsigned N):
    RooAbsPdf(name,title),
    x_("x","Dependent",this,x),
    paramList_("paramList","List of parameters",this) {

        paramList_.add(paramList);
        gSystem->Load("libMathMore");
        N_ = N;

    }

RooAtlas::RooAtlas(const RooAtlas & other, const char* name):
    RooAbsPdf(other,name),
    x_("x",this,other.x_),
    paramList_("paramList",this,other.paramList_),
    N_(other.N_) {

        gSystem->Load("libMathMore");

    }

Double_t RooAtlas::evaluate() const {

    double x = x_;

    //Part 1: (1-x^1/3)^b
    float b = static_cast<RooAbsReal*>(paramList_.at(0))->getVal();
    float part1 = pow((1-pow(x,1.0/3.0)),b); 

    //Part 2: inside the power of x (the sum of powers of logs)
    float part2(0.0);
    for (unsigned order(0); order <= N_; order++) {
        if (order == 0) {
            part2 += static_cast<RooAbsReal*>(paramList_.at(order+1))->getVal();
        }else if (order == 1) {
            part2 += static_cast<RooAbsReal*>(paramList_.at(order+1))->getVal()*log(x);
        }else{
            part2 += static_cast<RooAbsReal*>(paramList_.at(order+1))->getVal()*pow(log(x),order);
        }
    }
    
    //Part 3: x raised to the power of the logs part (2)
    float part3 = pow(x,part2);

    return part1*part3;
}

Int_t RooAtlas::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {

    if (rangeName && strlen(rangeName)) {
        return 0;
    }  
    if (matchArgs(allVars,analVars,x_)) return 1 ;
    return 0 ;
}

Double_t RooAtlas::analyticalIntegral(Int_t code, const char* rangeName) const {

    switch(code) {
        case 1:
            {
                double xMax = x_.max(rangeName);
                double xMin = x_.min(rangeName);
                return indefiniteAtlasIntegral(xMax) - indefiniteAtlasIntegral(xMin);
            }
    }

    assert(0);
    return 0;
}

Double_t RooAtlas::indefiniteAtlasIntegral(double xVal) const {

    double polyLog(0.0);
    for (unsigned order(0);order<=N_;order++) {
        if (order == 0) {
            polyLog += static_cast<RooAbsReal*>(paramList_.at(order+1))->getVal();
        }else if (order == 1) {
            polyLog += static_cast<RooAbsReal*>(paramList_.at(order+1))->getVal()*log(xVal);
        }else{
            polyLog += static_cast<RooAbsReal*>(paramList_.at(order+1))->getVal()*pow(log(xVal),order);
        }
    }
   
    Double_t ret =  (pow(xVal,1+polyLog) * ROOT::Math::hyperg( static_cast<RooAbsReal*>(paramList_.at(0))->getVal()*-1.0,
                                                               3*(1+polyLog),
                                                               3*(1+polyLog)+1,
                                                               pow(xVal,1.0/3.0) ))/(1+polyLog);
    return ret;
}




