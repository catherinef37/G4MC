#ifndef HRS_TRANSBASE_H
#define HRS_TRANSBASE_H

class HRSTransBase
{
public:
  HRSTransBase();
  virtual ~HRSTransBase();
  
  virtual bool TransLeftHRS(double* v) = 0;
  virtual bool TransRightHRS(double* v) = 0;
  virtual bool TransRightHRS_C(double* v) = 0;
  virtual void ReconLeftHRS(double* v) = 0;
  virtual void ReconRightHRS(double* v) = 0;
  virtual void Acceptance(double* v, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool) = 0;
  virtual void Acceptance_C(double* v, double* x_check, double* theta_check, double* y_check, double* phi_check, int* acc_bool) = 0;

  virtual double GetAngle() = 0;
  
  virtual void FPCorrLeft(const double* V5tg, double* V5fp);
  virtual void FPCorrRight(const double* V5tg, double* V5fp);
  
  // Rotate on X axis in transport coordinate
  // Notice the positive direction is anti-clockwise
  // Z is assumed to be 0
  void CoordsCorrection(double angle, double* v);
};

#endif
