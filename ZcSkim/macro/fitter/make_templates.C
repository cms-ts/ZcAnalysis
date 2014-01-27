#include <vector>

TFile * f[8];
TH1F * h_b[2][12];
TH1F * h_c[2][12];
TH1F * h_l[2][12];

TH1F * b_template;
TH1F * c_template;
TH1F * l_template;

TF1 * f_b;
TF1 * f_c;
TF1 * f_l;


const double Lumi2012_ele  = 19789.0;
const double Lumi2012_muon = 19751.0;

double w[] = { //31200./57709905.,
	       //54.838/10000431.,
	       //33.21/10000283.,
	       //8.059/9799908.,
	       //225.197/6923750.,
	       3503.71/30459503. };


class Chebyshev {
public: 
   Chebyshev(int n, double xmin, double xmax) : 
      fA(xmin), fB(xmax),
      fT(std::vector<double>(n) )  {}

   double operator() (const double * xx, const double *p) { 
      double x = (xx[0] - fA -fB)/(fB-fA);
      int order = fT.size(); 
      if (order == 1) return p[0]; 
      if (order == 2) return p[0] + x*p[1]; 
      // build the polynomials
      fT[0] = 1;
      fT[1] = x; 
      for (int i = 1; i< order; ++i) { 
         fT[i+1] =  2 *x * fT[i] - fT[i-1]; 
      }
      double sum = p[0]*fT[0]; 
      for (int i = 1; i<= order; ++i) { 
         sum += p[i] * fT[i]; 
      }
      return sum; 
   }

private: 
   double fA; 
   double fB; 
   std::vector<double> fT; // polynomial
   std::vector<double> fC; // coefficients
};





void make_templates(){

  const TString h_name = "hc_BJP";


  const TString version = "v01";

  const TString path = "/gpfs/cms/users/casarsa/analysis/Zc/work/output/" + version;


  vector <TString> file_names;
  //file_names.push_back("Wj.root");
  //file_names.push_back("WW.root");
  //file_names.push_back("WZ.root");
  //file_names.push_back("ZZ.root");
  //file_names.push_back("TTbar.root");
  file_names.push_back("DYJetsToLL.root");



  for (int ifile=0; ifile<file_names.size(); ++ifile){

    f[ifile] = new TFile(path + "/" + file_names[ifile]);

    h_b[0][ifile] =  (TH1F*) f[ifile]->Get("anaEle/" + h_name + "_ee_b")->Clone(); 
    h_b[1][ifile] =  (TH1F*) f[ifile]->Get("anaMuo/" + h_name + "_mm_b")->Clone(); 
    h_b[0][ifile]->Scale(Lumi2012_ele*w[ifile]); 
    h_b[1][ifile]->Scale(Lumi2012_muon*w[ifile]);

    h_c[0][ifile] =  (TH1F*) f[ifile]->Get("anaEle/" + h_name + "_ee_c")->Clone(); 
    h_c[1][ifile] =  (TH1F*) f[ifile]->Get("anaMuo/" + h_name + "_mm_c")->Clone(); 
    h_c[0][ifile]->Scale(Lumi2012_ele*w[ifile]); 
    h_c[1][ifile]->Scale(Lumi2012_muon*w[ifile]);

    h_l[0][ifile] =  (TH1F*) f[ifile]->Get("anaEle/" + h_name + "_ee_l")->Clone(); 
    h_l[1][ifile] =  (TH1F*) f[ifile]->Get("anaMuo/" + h_name + "_mm_l")->Clone(); 
    h_l[0][ifile]->Scale(Lumi2012_ele*w[ifile]); 
    h_l[1][ifile]->Scale(Lumi2012_muon*w[ifile]);

  }

  b_template = (TH1F*) h_b[0][0]->Clone("b_template");
  for (int i=1;i<file_names.size();++i) b_template->Add(h_b[0][i]);
  for (int i=0;i<file_names.size();++i) b_template->Add(h_b[1][i]);

  c_template = (TH1F*) h_c[0][0]->Clone("c_template");
  for (int i=1;i<file_names.size();++i) c_template->Add(h_c[0][i]);
  for (int i=0;i<file_names.size();++i) c_template->Add(h_c[1][i]);

  l_template = (TH1F*) h_l[0][0]->Clone("l_template");
  for (int i=1;i<file_names.size();++i) l_template->Add(h_l[0][i]);
  for (int i=0;i<file_names.size();++i) l_template->Add(h_l[1][i]);

  b_template->Draw("EP");


  f_b = new TF1("f_b","[0]*TMath::Exp(-[1]*x)*(1.+TMath::Erf((x-[2])/[3]))",0.,10.);
  f_b->SetParameters(1500.,1.,1.7,1.44);
  //b_template->Fit(f_b);
  
  f_c = new TF1("f_c","[0]*TMath::Exp(-[1]*x)*(1.+TMath::Erf((x-[2])/[3]))",0.,10.);
  f_c->SetParameters(1500.,1.,1.7,1.44);
  //c_template->Fit(f_c);
  
  f_l = new TF1("f_l","[0]*TMath::Exp(-[1]*x)*(1.+TMath::Erf((x-[2])/[3]))",0.,10.);
  f_l->SetParameters(1500.,1.,1.7,1.44);
  //l_template->Fit(f_l);


  // --- Use the Chebishev polynomials:

  int n = 12; 
  double xmin =  0.;
  double xmax = 10.;

  Chebyshev * cheb = new Chebyshev(n,xmin,xmax);
  TF1 * f1 = new TF1("f1",cheb,xmin,xmax,n+1,"Chebyshev");
  for (int i = 0; i <=n; ++i) f1->SetParameter(i,1);  


  b_template->SetMarkerStyle(20);
  b_template->SetMarkerColor(4);
  //b_template->Rebin();
  b_template->Scale(1./b_template->Integral());

  //b_template->Fit("f1");

  c_template->SetMarkerStyle(20);
  c_template->SetMarkerColor(4);
  //c_template->Rebin();
  c_template->Scale(1./c_template->Integral());

  c_template->Fit("f1");

  l_template->SetMarkerStyle(20);
  l_template->SetMarkerColor(4);
  //l_template->Rebin();
  l_template->Scale(1./l_template->Integral());

  //l_template->Fit("f1");


//  Double_t xbins[83];
//
//  //for (int i=1; i<=b_template->GetNbinsX(); ++i)
//  for (int i=1; i<=81; ++i)
//    xbins[i-1] = b_template->GetBinLowEdge(i);
//
//  //xbins[81] = 8.;
//  //xbins[82] = 10.;
//   
//  TH1F * h_rebinned = b_template->Rebin(82,"h_rebinned",xbins);
//  h_rebinned->Sumw2();
//  h_rebinned->Scale(1.,"width");
//  
//  
//  h_rebinned->Fit(f1,"I");


}

