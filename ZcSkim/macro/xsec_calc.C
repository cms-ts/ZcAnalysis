

void xsec_calc(){

  const double xsec_ee = 0.0534965;
  const double xsec_mm = 0.0606789;


  double N_tag_obs[2][2] = { 23513., sqrt(23513.),   // ee
  			     32719., sqrt(32719.)};  // mm
  double N_tag_bkg[2][2] = { 1266., 28.,
  			     1715., 31.};

  double N_inc_obs[2][2] = { 808161.,  sqrt(808161.),
  			    1139560., sqrt(1.13956e+06)};
  double N_inc_bkg[2][2] = {  9815., 69.,
			     13080., 71.};

  //N_tag_obs[0][0] = 23513;
  //N_tag_obs[1][0] = 32719;

  //N_tag_bkg[0][0] += -0.2*1266.;
  //N_tag_bkg[1][0] += -0.2*1715.;

  //N_inc_obs[0][0] = 808161;
  //N_inc_obs[1][0] = 1.13956e+06;

  //N_inc_bkg[0][0] += -0.2*9815.;
  //N_inc_bkg[1][0] += -0.2*13080.;


  const double f_c[2][2] = { 0.2424, 0.0181,
			     0.2726, 0.0158 }; 

  //const double eps_tag[2][2] = { 0.1351, 0.0022,   // ee
  //				 0.1319, 0.0018 };   // mm
  //const double eps_svx[2][2] = { 1., 0.,   // ee
  //  				 1., 0. };   // mm

  const double eps_tag[2][2] = { 0.40518, 0.00206,   // ee
  				 0.41344, 0.00176 }; // mm: 
  const double eps_svx[2][2] = { 0.33350, 0.00311,   
  				 0.31906, 0.00260 }; 

  const double A_c[2][2] = { 0.7481, 0.0043,
			     0.7013, 0.0038 };
  const double A_inc[2][2] = { 0.8006,  0.0011,
			       0.74818, 0.00097 };


  double xsec_ratio[2][3] = {0.,0.,0.,0.,0.,0.};

  for(int i=0; i<2; ++i) {

    xsec_ratio[i][0] = (N_tag_obs[i][0]-N_tag_bkg[i][0])/(N_inc_obs[i][0]-N_inc_bkg[i][0])*f_c[i][0]/
      (eps_tag[i][0]*eps_svx[i][0])*A_inc[i][0]/A_c[i][0];

    double w_NtagObs =  A_inc[i][0]*f_c[i][0]/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0])); 
    double w_NtagBkg = -A_inc[i][0]*f_c[i][0]/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0]));
    double w_NincObs = -A_inc[i][0]*f_c[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*pow(N_inc_obs[i][0]-N_inc_bkg[i][0],2));
    double w_NincBkg =  A_inc[i][0]*f_c[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*pow(N_inc_obs[i][0]-N_inc_bkg[i][0],2));
    double w_fc      =  A_inc[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0]));
    double w_epsc    = -A_inc[i][0]*f_c[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*eps_tag[i][0]*eps_tag[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0]));
    double w_epssvx  = -A_inc[i][0]*f_c[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0]));
    double w_Ac      = -A_inc[i][0]*f_c[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0]));
    double w_Ainc    =  f_c[i][0]*(N_tag_obs[i][0]-N_tag_bkg[i][0])/(A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*(N_inc_obs[i][0]-N_inc_bkg[i][0]));

    // statistical error:
    xsec_ratio[i][1] = sqrt( pow(w_NtagObs*N_tag_obs[i][1],2) + pow(w_NincObs*N_inc_obs[i][1],2) +
			     pow(w_NtagBkg*N_tag_bkg[i][1],2) + pow(w_NincBkg*N_inc_bkg[i][1],2) );

    //cout <<  xsec_ratio[i][1]  << " ";
    //xsec_ratio[i][1] = sqrt( pow(-A_inc[i][0]*f_c[i][0]*(N_inc_bkg[i][0]-N_tag_bkg[i][0]-N_inc_obs[i][0]+N_tag_obs[i][0])/
    //				 (A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*pow(N_inc_obs[i][0]-N_inc_bkg[i][0],2)),2)*N_tag_obs[i][0] +
    //			     pow( A_inc[i][0]*f_c[i][0]*(N_inc_bkg[i][0]-N_tag_bkg[i][0]-N_inc_obs[i][0]+N_tag_obs[i][0])/
    //				  (A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*pow(N_inc_obs[i][0]-N_inc_bkg[i][0],2)),2)*N_tag_bkg[i][0] +
    //			     pow( A_inc[i][0]*f_c[i][0]*(N_tag_bkg[i][0]-N_tag_obs[i][0])/
    //				  (A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*pow(N_inc_bkg[i][0]-N_inc_obs[i][0],2)),2)*(N_inc_obs[i][0]-N_tag_obs[i][0]) +
    //			     pow(-A_inc[i][0]*f_c[i][0]*(N_tag_bkg[i][0]-N_tag_obs[i][0])/
    //				 (A_c[i][0]*eps_tag[i][0]*eps_svx[i][0]*pow(N_inc_bkg[i][0]-N_inc_obs[i][0],2)),2)*(N_inc_bkg[i][0]-N_tag_bkg[i][0]));
    //cout << xsec_ratio[i][1]  << endl;


    // systematic error:
    xsec_ratio[i][2] = sqrt( pow(w_fc*f_c[i][1], 2) +   
			     pow(w_epsc*eps_tag[i][1], 2) + 
			     pow(w_epssvx*eps_svx[i][1], 2) + 
			     pow(w_Ac*A_c[i][1], 2) +   
			     pow(w_Ainc*A_inc[i][1], 2) ); 



    //xsec_ratio[i][1] = 
    //  pow(f_c[i][0]/(N_inc_obs[i][0]-N_inc_bkg[i][0])*A_inc[i][0]/A_c[i][0]*N_tag_obs[i][1]/eps_tag[i][0],2) +
    //  pow(f_c[i][0]/(N_inc_obs[i][0]-N_inc_bkg[i][0])*A_inc[i][0]/A_c[i][0]*N_tag_bkg[i][1]/eps_tag[i][0],2) +
    //  pow((N_tag_obs[i][0]-N_tag_bkg[i][0])/(N_inc_obs[i][0]-N_inc_bkg[i][0])*A_inc[i][0]/A_c[i][0]*f_c[i][1]/eps_tag[i][0],2) +
    //  pow((N_tag_obs[i][0]-N_tag_bkg[i][0])/pow(N_inc_obs[i][0]-N_inc_bkg[i][0],2)*f_c[i][0]/eps_tag[i][0]*A_inc[i][0]/A_c[i][0]*N_inc_obs[i][1],2) +
    //  pow((N_tag_obs[i][0]-N_tag_bkg[i][0])/pow(N_inc_obs[i][0]-N_inc_bkg[i][0],2)*f_c[i][0]/eps_tag[i][0]*A_inc[i][0]/A_c[i][0]*N_inc_bkg[i][1],2) +
    //  pow((N_tag_obs[i][0]-N_tag_bkg[i][0])/(N_inc_obs[i][0]-N_inc_bkg[i][0])*f_c[i][0]/pow(eps_tag[i][0],2)*A_inc[i][0]/A_c[i][0]*eps_tag[i][1],2) +
    //  pow((N_tag_obs[i][0]-N_tag_bkg[i][0])/(N_inc_obs[i][0]-N_inc_bkg[i][0])*f_c[i][0]/eps_tag[i][0]*A_inc[i][1]/A_c[i][0],2) +
    //  pow((N_tag_obs[i][0]-N_tag_bkg[i][0])/(N_inc_obs[i][0]-N_inc_bkg[i][0])*f_c[i][0]/eps_tag[i][0]*A_inc[i][0]/pow(A_c[i][0],2)*A_c[i][1],2);
    //
    //xsec_ratio[i][1] = sqrt(xsec_ratio[i][1]);
      
  }


  cout << xsec_ratio[0][0] << " +- " << xsec_ratio[0][1] << " +- " << xsec_ratio[0][2] << endl;
  cout << xsec_ratio[1][0] << " +- " << xsec_ratio[1][1] << " +- " << xsec_ratio[1][2] << endl;

  cout << endl;
  cout <<  xsec_ee- xsec_ratio[0][0] << " " <<  xsec_mm- xsec_ratio[1][0] << endl;

}
