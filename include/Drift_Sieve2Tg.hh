#ifndef _DRIFTINFIELD_H_
#define _DRIFTINFIELD_H_

namespace DriftSieve2Tg
{
	void SetUseSeptumInDrift(int val);
	
	//This routine is used to drift particle in the field, will stop once z_tr across z_tr_limit 
	//or track_length>tracklengthlimit 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	int Drift2Ztr_Fast(TVector3& VX, TVector3& VP, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget);

	//This routine is used to drift particle in the field, will stop once z_tr across z_tr_limit 
	//or track_length>tracklengthlimit 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	//if steplength<1.0E-4, it will be used while approching the limit, otherwise
	//it will be used to control print interval:  print_interval = 1.0E-3/steplength
	int Drift2Ztr(TVector3& VX, TVector3& VP, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget);


	int Drift2Ztr_EqualStep(TVector3& VX, TVector3& VP, double angle_rad, double m_gev, double q,
		double z_tr_limit, double tracklengthlimit, double steplength, double ztarget);

	
	int Drift2Ztr(const double *V5In, double *V5Out, double& z_tr, double P0_gev, double angle_rad, 
		double m_gev, double q, double z_tr_limit, double tracklengthlimit, double steplength, double ztarget);

	//This routine is used to drift particle in the field, will stop once z across zlimit 
	//or track_length>tracklengthlimit 
	//definition of backward and forward: If Pz<0 go backward, otherwise forward
	//if steplength<1.0E-4, it will be used while approching the limit, otherwise
	//it will be used to control print interval:  print_interval = 1.0E-3/steplength
	int Drift2Z(TVector3& VX, TVector3& VP, double m_gev, double q, double zlimit,
		double tracklengthlimit, double steplength);

	//The following were developed to test if we can match 2 trajectories to find the 
	//most close trajectory for reconstruction
	//It turns out that matching trajectory does not work

	//will return the trajectary
	int DriftPath(TVector3& VX, TVector3& VP, double m_gev, double q, double zlimit,
		double tracklengthlimit, double steplength, vector < vector <double> > &traj);
	
	//compare 2 trjectories
	bool Cmp2Traj( vector < vector <double> > &trajrec, vector < vector <double> > &trajbpm, 
		size_t &indexI, size_t &indexJ, double &dist);
}

#endif //_DRIFTINFIELD_H_
