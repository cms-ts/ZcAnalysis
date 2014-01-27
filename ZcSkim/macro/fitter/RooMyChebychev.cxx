#include <cmath>
#include <iostream>

#include "RooFit.h"

#include "Riostream.h"

#include "RooMyChebychev.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooArgList.h"

#if defined(__my_func__)
#undef __my_func__
#endif
#if defined(WIN32)
#define __my_func__ __FUNCTION__
#else
#define __my_func__ __func__
#endif


ClassImp(RooMyChebychev)
  ;


//_____________________________________________________________________________
RooMyChebychev::RooMyChebychev()
{
}


//_____________________________________________________________________________
RooMyChebychev::RooMyChebychev(const char* name, const char* title, 
			       RooAbsReal& x, const RooArgList& coefList): 
  RooAbsPdf(name, title),
  _x("x", "Dependent", this, x),
  _coefList("coefficients","List of coefficients",this)
{
  // Constructor

  TIterator* coefIter = coefList.createIterator() ;
  RooAbsArg* coef ;
  while((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      std::cerr << "RooMyChebychev::ctor(" << GetName() <<
	") ERROR: coefficient " << coef->GetName() <<
	" is not of type RooAbsReal" << std::endl ;
      assert(0) ;
    }
    _coefList.add(*coef) ;
  }
  delete coefIter ;
}



//_____________________________________________________________________________
RooMyChebychev::RooMyChebychev(const RooMyChebychev& other, const char* name) :
  RooAbsPdf(other, name), 
  _x("x", this, other._x), 
  _coefList("coefList",this,other._coefList)
{
}


//_____________________________________________________________________________
Double_t RooMyChebychev::evaluate() const 
{

  Double_t x = (_x-_x.min()-_x.max())/(_x.max()-_x.min());

  const int order = _coefList.getSize()-1;

  std::vector<double> fT(order); // polynomial
  std::vector<double> fC(order); // coefficients
  
  if (order == 0) return ((RooAbsReal&)_coefList[0]).getVal(); 
  if (order == 1) return ((RooAbsReal&)_coefList[0]).getVal() + x*((RooAbsReal&)_coefList[1]).getVal(); 
  
  // build the polynomials
  fT[0] = 1;
  fT[1] = x; 
  for (int i = 1; i< order; ++i) { 
    fT[i+1] =  2 *x * fT[i] - fT[i-1]; 
  }
  Double_t sum = ((RooAbsReal&)_coefList[0]).getVal()*fT[0]; 
  for (int i = 1; i<= order; ++i) { 
    sum += ((RooAbsReal&)_coefList[i]).getVal() * fT[i]; 
  }

  return sum;

}

#undef __my_func__
