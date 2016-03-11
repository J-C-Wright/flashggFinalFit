


#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "../interface/RooAtlas.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooArgList.h"

using namespace std;

ClassImp(RooAtlas);

RooAtlas::RooAtlas(const char *name, const char *title, RooAbsReal &x, const RooArgList &paramList):
    RooAbsPdf(name,title),
    x_("x","Dependent",this,x),
    paramList_("paramList","List of parameters",this) {
        paramList_.add(paramList);
    }

RooAtlas::RooAtlas(const RooAtlas & other, const char* name):
    RooAbsPdf(other,name),
    x_("x",this,other.x_),
    paramList_("paramList",this,other.paramList_) {
    }

RooAtlas::evaluate() const {

    //Part 1: (1-x^1/3)&^b
    float b = paramList_.at(0)->getVal();
    float part1 = pow((1-pow(x_,1.0/3.0)),b); 

    //Part 2: inside the power of x (the sum of powers of logs)
    float part2(0.0);
    for (unsigned order(0); order <= N; order++) {
        if (order == 0) {
            part2 += paramList_.at(order+1)->getVal();
        }else if (order == 1) {
            part2 += paramList_.at(order+1)->getVal()*log(x_);
        }else{
            part2 += paramList_.at(order+1)->getVal()*pow(log(x_),order);
        }
    }
    
    //Part 3: x raised to the power of the logs part (2)
    float part3 = pow(x_,part2);

    return part1*part3;
}

Int_t RooPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName = 0) const {
    if (rangeName && strlen(rangeName)) {
        return 0;
    }  
    if (matchArgs(allVars,analVars,x)) return 1 ;
    return 0 ;
}

Double_t RooPowerLaw::analyticalIntegral(Int_t code, const char* rangeName) const 
{






