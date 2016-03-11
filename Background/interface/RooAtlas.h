#ifndef ROO_ATLAS
#define ROO_ATLAS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"


template <unsigned N> class RooAtlas : public RooAbsPdf {

    public:
        RooAtlas(){};
        RooAtlas(const char *name, const char *title, RooAbsReal &x, const RooArgList &paramList);
        RooAtlas(const RooAtlas & other, const char* name = 0);
        virtual TObject* clone(const char* newname) const { return new RooAtlas(*this,newname); }
        virtual ~RooBernsteinFast(){};

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
        Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

    protected:
        RooRealProxy x_;
        RooListProxy paramList_;

        Double_t evaluate() const;

        ClassDef(RooAtlas,1)

};

#endif
