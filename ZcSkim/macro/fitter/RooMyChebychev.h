#ifndef ROO_MYCHEBYCHEV
#define ROO_MYCHEBYCHEV

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList ;

class RooMyChebychev : public RooAbsPdf {
 public:

  RooMyChebychev() ;
  RooMyChebychev(const char *name, const char *title,
		 RooAbsReal& _x, const RooArgList& _coefList) ;

  RooMyChebychev(const RooMyChebychev& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { return new RooMyChebychev(*this, newname); }
  inline virtual ~RooMyChebychev() { }

 private:

  RooRealProxy _x;
  RooListProxy _coefList ;

  //double fA; 
  //double fB; 
  //std::vector<double> fT; // polynomial
  //std::vector<double> fC; // coefficients

  Double_t evaluate() const;


  ClassDef(RooMyChebychev,1) // Chebychev polynomial PDF
    };

#endif
