#ifndef ROO_ATLAS
#define ROO_ATLAS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooAtlas : public RooAbsPdf {

    public:
        RooAtlas(){};
        RooAtlas(const char *name, const char *title, RooAbsReal &x, const RooArgList &paramList,unsigned N);
        RooAtlas(const RooAtlas & other, const char* name = 0);
        virtual TObject* clone(const char* newname) const { return new RooAtlas(*this,newname); }
        virtual ~RooAtlas(){};

        Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
        Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

        Double_t indefiniteAtlasIntegral(double xVal);

    protected:
        RooRealProxy x_;
        RooListProxy paramList_;
        unsigned N_;

        Double_t evaluate() const;

        ClassDef(RooAtlas,1)

};

#endif
