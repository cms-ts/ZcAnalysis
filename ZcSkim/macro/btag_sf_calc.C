void btag_sf_calc(TString flavour="b", TString tagger="CSVL"){

 
  float ptmin[16] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
  float ptmax[16] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};



  // ==========================================================================================
  //  b-tag SF for heavy flavours

  if (  flavour == "b" ) {

    // from: https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFb-pt_NOttbar_payload_EPS13.txt
  

    float SFb_error_CSVL[16] = {
      0.033408,
      0.015446,
      0.0146992,
      0.0183964,
      0.0185363,
      0.0145547,
      0.0176743,
      0.0203609,
      0.0143342,
      0.0148771,
      0.0157936,
      0.0176496,
      0.0209156,
      0.0278529,
      0.0346877,
      0.0350101 };
  

    float SFb_error_CSVT[16] = {
      0.0511028,
      0.0306671,
      0.0317498,
      0.032779,
      0.0291528,
      0.0249308,
      0.0301118,
      0.032047,
      0.0348072,
      0.0357745,
      0.0378756,
      0.0412608,
      0.0777516,
      0.0860741,
      0.0942209,
      0.104106 };
  
 
    float SFb[16] = {0.};
    float SFb_err[16] = {0.};


    for (int ipt=0; ipt<16; ++ipt){

      double x   = 0.5*(ptmax[ipt]+ptmin[ipt]);
      float etamin = -2.4;
      float etamax =  2.4;


      if ( tagger == "CSVL" ){
	SFb[ipt] = 1.00572*((1.+(0.013676*x))/(1.+(0.0143279*x)));
	SFb_err[ipt] = SFb_error_CSVL[ipt];
      }
      else if ( tagger == "CSVT" ){
	SFb[ipt] = (0.9203+(-3.32421e-05*x))+(-7.74664e-08*(x*x));
	SFb_err[ipt] = SFb_error_CSVT[ipt];
      }
    
      cout << ptmin[ipt]   << "\t" 
	   << ptmax[ipt]   << "\t"
	   << etamin     << "\t" 
	   << etamax     << "\t"
	   << SFb[ipt]     << "\t"
	   << SFb_err[ipt] << "\t"
	   << SFb_err[ipt] << endl;


    }

  }


  // ==========================================================================================
  //  b-tag SF for light flavours

  if (  flavour == "l" ) {

    // from: https://twiki.cern.ch/twiki/pub/CMS/BtagPOG/SFlightFuncs_EPS2013.C 
    //
    // if( Atagger == "CSVL" && sEtamin == "0.0" && sEtamax == "0.5")
    // {
    // if( meanminmax == "mean" ) tmpSFl = new TF1("SFlight","((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "min" ) tmpSFl = new TF1("SFlightMin","((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "max" ) tmpSFl = new TF1("SFlightMax","((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)))", 20.,ptmax);
    // }
    // if( Atagger == "CSVL" && sEtamin == "0.5" && sEtamax == "1.0")
    // {
    // if( meanminmax == "mean" ) tmpSFl = new TF1("SFlight","((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "min" ) tmpSFl = new TF1("SFlightMin","((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "max" ) tmpSFl = new TF1("SFlightMax","((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)))", 20.,ptmax);
    // }
    // if( Atagger == "CSVL" && sEtamin == "1.0" && sEtamax == "1.5")
    // {
    // if( meanminmax == "mean" ) tmpSFl = new TF1("SFlight","((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "min" ) tmpSFl = new TF1("SFlightMin","((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "max" ) tmpSFl = new TF1("SFlightMax","((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)))", 20.,ptmax);
    // }
    // if( Atagger == "CSVL" && sEtamin == "1.5" && sEtamax == "2.4")
    // {
    // if( meanminmax == "mean" ) tmpSFl = new TF1("SFlight","((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "min" ) tmpSFl = new TF1("SFlightMin","((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)))", 20.,ptmax);
    // if( meanminmax == "max" ) tmpSFl = new TF1("SFlightMax","((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)))", 20.,ptmax);
    // }


    float etamin[4] = {0.0, 0.5, 1.0, 1.5};
    float etamax[4] = {0.5, 1.0, 1.5, 2.4};


    for (int ipt=0; ipt<16; ++ipt){

      double x = 0.5*(ptmax[ipt]+ptmin[ipt]);

      double SFl     = 1.;    
      double SFl_min = 1.;
      double SFl_max = 1.;
	

      if ( tagger == "CSVL" ) {


	// eta 0.0-0.5

	SFl     = ((1.01177+(0.0023066*x))+(-4.56052e-06*(x*x)))+(2.57917e-09*(x*(x*x)));
	SFl_min = ((0.977761+(0.00170704*x))+(-3.2197e-06*(x*x)))+(1.78139e-09*(x*(x*x)));
	SFl_max = ((1.04582+(0.00290226*x))+(-5.89124e-06*(x*x)))+(3.37128e-09*(x*(x*x)));

	cout << ptmin[ipt] << " \t" 
	     << ptmax[ipt] << " \t"
	     << etamin[0]  << " \t" 
	     << etamax[0]  << " \t"
	     << SFl        << " \t"
	     << SFl_min    << " \t"
	     << SFl_max    << endl;
      

	// eta 0.5-1.0

	SFl     = ((0.975966+(0.00196354*x))+(-3.83768e-06*(x*x)))+(2.17466e-09*(x*(x*x)));
	SFl_min = ((0.945135+(0.00146006*x))+(-2.70048e-06*(x*x)))+(1.4883e-09*(x*(x*x)));
	SFl_max = ((1.00683+(0.00246404*x))+(-4.96729e-06*(x*x)))+(2.85697e-09*(x*(x*x)));

	cout << ptmin[ipt] << " \t" 
	     << ptmax[ipt] << " \t"
	     << etamin[1]  << " \t" 
	     << etamax[1]  << " \t"
	     << SFl        << " \t"
	     << SFl_min    << " \t"
	     << SFl_max    << endl;


	// eta 1.0-1.5

	SFl     = ((0.93821+(0.00180935*x))+(-3.86937e-06*(x*x)))+(2.43222e-09*(x*(x*x)));
	SFl_min = ((0.911657+(0.00142008*x))+(-2.87569e-06*(x*x)))+(1.76619e-09*(x*(x*x)));
	SFl_max = ((0.964787+(0.00219574*x))+(-4.85552e-06*(x*x)))+(3.09457e-09*(x*(x*x)));

	cout << ptmin[ipt] << " \t" 
	     << ptmax[ipt] << " \t"
	     << etamin[2]  << " \t" 
	     << etamax[2]  << " \t"
	     << SFl        << " \t"
	     << SFl_min    << " \t"
	     << SFl_max    << endl;


	// eta 1.5-2.4

	SFl     = ((1.00022+(0.0010998*x))+(-3.10672e-06*(x*x)))+(2.35006e-09*(x*(x*x)));
	SFl_min = ((0.970045+(0.000862284*x))+(-2.31714e-06*(x*x)))+(1.68866e-09*(x*(x*x)));
	SFl_max = ((1.03039+(0.0013358*x))+(-3.89284e-06*(x*x)))+(3.01155e-09*(x*(x*x)));

	cout << ptmin[ipt] << " \t" 
	     << ptmax[ipt] << " \t"
	     << etamin[3]  << " \t" 
	     << etamax[3]  << " \t"
	     << SFl        << " \t"
	     << SFl_min    << " \t"
	     << SFl_max    << endl;


      }
      else if ( tagger == "CSVT" ){

	SFl     = ((1.00462+(0.00325971*x))+(-7.79184e-06*(x*x)))+(5.22506e-09*(x*(x*x)));
	SFl_min = ((0.845757+(0.00186422*x))+(-4.6133e-06*(x*x)))+(3.21723e-09*(x*(x*x)));
	SFl_max = ((1.16361+(0.00464695*x))+(-1.09467e-05*(x*x)))+(7.21896e-09*(x*(x*x)));

	cout << ptmin[ipt] << " \t" 
	     << ptmax[ipt] << " \t"
	     << -etamax[3] << " \t" 
	     << etamax[3]  << " \t"
	     << SFl        << " \t"
	     << SFl_min    << " \t"
	     << SFl_max    << endl;

      }


    }



//    gSystem->CompileMacro("SFlightFuncs_EPS2013.C");
//
//    TF1* SFlight[4]  = { GetSFlmean("CSV","L",0.0, 0.5, "ABCD"),
//			 GetSFlmean("CSV","L",0.5, 1.0, "ABCD"),
//			 GetSFlmean("CSV","L",1.0, 1.5, "ABCD"),
//			 GetSFlmean("CSV","L",1.5, 2.4, "ABCD") };
//
//    TF1* SFlightmin[4]  = { GetSFlmin("CSV","L",0.0, 0.5, "ABCD"),
//			    GetSFlmin("CSV","L",0.5, 1.0, "ABCD"),
//			    GetSFlmin("CSV","L",1.0, 1.5, "ABCD"),
//			    GetSFlmin("CSV","L",1.5, 2.4, "ABCD") };
//
//    TF1* SFlightmax[4]  = { GetSFlmax("CSV","L",0.0, 0.5, "ABCD"),
//			    GetSFlmax("CSV","L",0.5, 1.0, "ABCD"),
//			    GetSFlmax("CSV","L",1.0, 1.5, "ABCD"),
//			    GetSFlmax("CSV","L",1.5, 2.4, "ABCD") };
//
//    float etamin[4] = {0.0, 0.5, 1.0, 1.5};
//    float etamax[4] = {0.5, 1.0, 1.5, 2.4};
//
//
//    for (int ipt=0; ipt<16; ++ipt){
//
//      double x = 0.5*(ptmax[ipt]+ptmin[ipt]);
//
//      for ( int ieta=0; ieta<4; ++ieta ) {
//
//	cout << ptmin[ipt]   << "\t" 
//	     << ptmax[ipt]   << "\t"
//	     << etamin[ieta] << "\t" 
//	     << etamax[ieta] << "\t"
//	     << SFlight[ieta]->Eval(x) << "\t"
//	     << SFlightmin[ieta]->Eval(x) << "\t"
//	     << SFlightmax[ieta]->Eval(x) <<  endl;
//
//      }
//
//    }
//
//
//  }





  }



}

