
#ifndef _COMPTON_H
#define _COMPTON_H

namespace Compton
{
  //input: photon incident energy in GeV and 
  //       photon scattering angle in CM frame in radian 
  //return RCS XS in nb/GeV^2/sr
  double GetRCSXS(double Ei, double Theta_cm_g);
  double GetRCSKLL(double Ei, double Theta_cm_g);
  double GetRCSALL(double Ei, double Theta_cm_g);
	  
  double GetTheta_cm(double Ei, double Theta_lab_g);

  double GetRCSXS_lab(double Ei, double Theta_lab_g);
  double GetRCSKLL_lab(double Ei, double Theta_lab_g);
  double GetRCSALL_lab(double Ei, double Theta_lab_g);

  //root fcn
  double XSVsE_FixedThetalab_deg(double *x, double *par);
  double XSVsE_FixedThetacm_deg(double *x, double *par);

  //x is cos(Theta_cm)
  double RCSXSVsCosTh_cm(double *x, double *par);
  double RCSKLLVsCosTh_cm(double *x, double *par);
  double RCSALLVsCosTh_cm(double *x, double *par);

  //x is Theta_cm
  double RCSXSVsTh_cm(double *x, double *par);
  double RCSKLLVsTh_cm(double *x, double *par);
  double RCSALLVsTh_cm(double *x, double *par);

}

#endif  //_COMPTON_H

