//2009 version of P. Bosted's electron XS model
//
//Fits to inclusive inelastic electron scattering (see M.E. Christy and P.E. Bosted, 
//arXiv:0711.0159 (submitted to PRC), for proton fit, and P.E. Bosted and M.E. Christy, 
//arXiv:0711.0159 [also Phys. Rev. C 77, 065206 (2008)], for deuteron and neutron fit.
// This 2007 code also includes the quasi-elastic and inelastic model for nuclei shown 
//in Phys. Rev. C. 78, 015202 (2008). (arXiv:0712.2438) (Ratios of 15N/12C and 4He/12C 
//inclusive electroproduction cross sections in the nucleon resonance region, 
//P.E. Bosted, R. Fersch {\it et al.}). 
//The 2009 fit is the same for proton and deuteron as the 2007 fit, but for A>2 the 
//fit has been greatly improved by the addition of new data from JLab, and an additional 
//25 parameters for A-dependence. This code works well up to copper, but may have problems 
//for heavier nuclei. It works well for 4He, but the model for 3He still needs tweaking. 
//The reference for the 2009 $A>2$ fit is: http://arxiv.org/abs/1203.2262. 
//

namespace PBosted 
{
	extern "C" 
	{
		void bosted_(double* Z, double* A, double* Ei, double* Ep, double* ang, double* xs, double *Tb, double *Ta);
	}

	//input
	//       Z,N: proton and neutron number of the nucleus.	;
	//Ei, Ef: incoming and outgoing electron energy in GeV;
	//Theta: scattering angle for outgoing particle in radian;
	//Tb and Ta will be used for radiated XS only if they are both positive
	//Tb: material thickness in unit of radiation length before scattering;
	//Ta: material thickness in unit of radiation length after scattering;
	double GetXS(int Z, int N, double Ei, double Ef, double theta, double Tb, double Ta)
	{
		double NZ, NA;
		NZ = Z;
		NA = Z+N;
		double XS;
#ifdef WIN32
		//this fortran routine does not work in windows, do not know why
		return 1.0;
#else
		bosted_(&NZ, &NA, &Ei, &Ef, &theta, &XS, &Tb, &Ta);
#endif
		return XS;
	}
}

