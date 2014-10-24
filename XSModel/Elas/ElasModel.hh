//return elas xs in micro barn
namespace  ElasModel
{
	//input
	//		Z,N: proton and neutron number of the nucleus.	;
	//		Ei: incoming electron energy in GeV;
	//		Theta: scattering angle for outgoing particle in radian;
	//		Mtg_GeV: target mass if it is given
	//		iUseLarry: in case for C12, use larry's model
	double GetXS(int Z, int N, double Ei, double Theta, double Mtg_GeV=0.0, int iUseLarry=0);

	//John Arrington's Form factor
	double XS_H(double Ei, double theta);
	//Xiaohui's Form factor
	double XS_H1(double Ei, double theta);
	double XS_He3(double Ei, double theta);
	double XS_He4(double Ei, double theta);
	double XS_N14(double Ei, double theta);
	double XS_Ta184(double Ei, double theta);
	//wrapper for Stansfield's elas C12 XS
	double XS_C12(double Ei, double theta);
	//wrapper for L. Cardman's elas C12 XS
	double Carbon(double Ei, double theta);
	double GetXS_All(double Ei, double theta, int iZ, double Mtg_GeV);
}

