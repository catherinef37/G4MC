C Forward transfer functions for right hrs with septum based on prex_Rn_dir.dat
c  setup for (x|theta)=(y|phi)=0 143 cm downstream of the 1st VDC
c HRS + PREX room temperature septum (right side)
c                     -JJL 3/29/2010

c Tune "B"

c typical call: answer = function(x,5)
c INPUTS: x = 5 or more element array 
c              x(1)=x0  (meters)
c              x(2)=theta0 (really tan(theta0))
c              x(3)=y0   (meters)
c              x(4)=phi0 (really tan(phi0))
c              x(5)=delta (fractional value NOT percent)
c         M=5
c
c OUTPUT: units are the same as inputs
c 
c NOMENCLATURE: function name = prefix + _sp_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means path length
c     
c           suffixes:     fp means target to focus
c                         q1ex means target to Q1 exit
c                         q2ex means target to Q2 exit
c                         den  means target to Dipole entrance
c                         dex  means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3m  means target to Q3 middle
c                         q3ex means target to Q3 exit
c                         sen  means septum entrance
c                         sm   means halfway through the septum
c                         sex  means septum exit
c                         col  means target to Q1 collimator
c                         cq1x means collimator to q1 exit
c                         cden means collimator to dipole entrance
c                         cdex means collimator to dipole exit
c                         cq3e means collimator to q3 entrance
c                         cq3x means collimator to q3 exit
c                         cfp  means collimator to focus
c                               i.e. The plane perpendicular to the optic axis
c                                    that crosses the 1st VDC
c
c          _sp_ is for septum PREX
c
c APERTURES:
c    Coordinate systems are different in the spectrometer with septum model compared to the spectrometer
c     without septum. So some apertures even in the body of the spectrometer appear to be different. Numbers
c    given here supercede any numbers given in regard to transfer functions for the spectrometers alone.
c
c     sen: +0.088 < x < +0.382
c          -0.120 < y < 0.120
c      sm: +0.088 < x < +0.382
c          -0.120 < y < 0.120
c     sex: +0.088 < x < +0.382
c          -0.120 < y < 0.120
c     q1ex is a circle of radius 0.1492 m
c     q2ex is a circle of radius 0.30 m
c     den is a trapazoid:
c                                   -5.22008 < x < -4.98099
c             -(-0.192436784*x -0.192436784) < y < -0.192436784*x -0.192436784  
c     dex is also a trapazoid: 
c                                   -0.46188 < x < 0.46188
c                   -(-0.01610808*x + 0.125) < y < -0.01610808*x + 0.125
c     q3en is a circle of radius 0.3 m
c     q3m  is a circle of radius 0.3 m
c     q3ex is a circle of radius 0.3 m
c
      function x_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2230503E-02/
      data xmin/
     1 -0.49991E-02,-0.51642E-01,-0.19979E-01,-0.27954E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16881146E-02, 0.12166227E+00, 0.46531553E-02,-0.45058406E-02,
     +  0.16646635E-02, 0.13803524E-02,-0.23503867E-02, 0.20201562E-02,
     +  0.32496466E-02,-0.76183904E-03, 0.22430087E-02,-0.39648658E-03,
     +  0.22263671E-03, 0.22040548E-03, 0.11533081E-02,-0.22741142E-02,
     + -0.21509791E-02, 0.61018619E-03,-0.73291356E-03,-0.33171859E-03,
     + -0.62142807E-03, 0.90805785E-03,-0.14014813E-02, 0.34017154E-04,
     + -0.15810892E-03,-0.48574744E-04, 0.19936594E-03, 0.23102765E-03,
     + -0.22580511E-03, 0.38600256E-03,-0.67539362E-03, 0.20528927E-02,
     +  0.22254861E-02, 0.14769036E-03,-0.34486738E-03,-0.33123375E-03,
     +  0.12915561E-03,-0.76602330E-04, 0.17838451E-03, 0.13374116E-02,
     +  0.13854977E-02, 0.34569330E-04, 0.59218367E-03,-0.19015036E-03,
     + -0.25176993E-03, 0.29583040E-03,-0.77898676E-05,-0.84275857E-03,
     + -0.34330951E-03,-0.24747330E-03,-0.66306151E-04, 0.27772678E-04,
     +  0.12913611E-04, 0.50088613E-04,-0.70766575E-04,-0.80182865E-04,
     + -0.59350248E-03,-0.48804528E-03, 0.62724166E-05, 0.13219015E-03,
     +  0.24593057E-04,-0.10526277E-03,-0.70986178E-04, 0.27334271E-03,
     +  0.20675342E-03,-0.38659917E-04, 0.54054850E-04, 0.28668207E-03,
     + -0.10026641E-03,-0.81974395E-04, 0.17081213E-03,-0.70855411E-03,
     +  0.13296583E-03,-0.60886534E-04,-0.10244772E-03, 0.78106450E-03,
     +  0.64975850E-03,-0.11675007E-03,-0.66184206E-04, 0.61581602E-04,
     +  0.83694402E-04,-0.87552813E-04,-0.53503365E-04,-0.14775005E-03,
     + -0.36074438E-04,-0.44224864E-04,-0.11281381E-03,-0.65346831E-04,
     +  0.48996637E-04, 0.24571348E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x23*x31        
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x23            
      x_sp_col    =x_sp_col    
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)        *x31        
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x42    
      x_sp_col    =x_sp_col    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x24            
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21*x32        
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)                *x51
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)    *x22        *x51
      x_sp_col    =x_sp_col    
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x22*x31        
     4  +coeff( 31)    *x24*x31        
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x23    *x42    
     7  +coeff( 34)*x11*x21*x33        
     8  +coeff( 35)        *x31*x41    
      x_sp_col    =x_sp_col    
     9  +coeff( 36)            *x42    
     1  +coeff( 37)    *x21*x31    *x51
     2  +coeff( 38)    *x23        *x51
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)    *x22*x31*x41    
     5  +coeff( 41)    *x22    *x42    
     6  +coeff( 42)*x12*x21            
     7  +coeff( 43)    *x23*x32        
     8  +coeff( 44)*x11*x24            
      x_sp_col    =x_sp_col    
     9  +coeff( 45)*x11*x23    *x41    
     1  +coeff( 46)    *x21*x31*x43    
     2  +coeff( 47)*x12    *x32        
     3  +coeff( 48)    *x24    *x42    
     4  +coeff( 49)*x11*x24    *x41    
     5  +coeff( 50)*x11*x23*x33        
     6  +coeff( 51)        *x32        
     7  +coeff( 52)            *x41*x51
     8  +coeff( 53)*x11            *x51
      x_sp_col    =x_sp_col    
     9  +coeff( 54)*x11*x21*x31        
     1  +coeff( 55)    *x22    *x41*x51
     2  +coeff( 56)    *x21*x32*x41    
     3  +coeff( 57)    *x21*x31*x42    
     4  +coeff( 58)    *x21    *x43    
     5  +coeff( 59)*x12    *x31        
     6  +coeff( 60)*x11*x21*x31*x41    
     7  +coeff( 61)*x11*x21    *x42    
     8  +coeff( 62)    *x23*x31    *x51
      x_sp_col    =x_sp_col    
     9  +coeff( 63)    *x23    *x41*x51
     1  +coeff( 64)    *x22*x31*x42    
     2  +coeff( 65)    *x22    *x43    
     3  +coeff( 66)*x11    *x31    *x52
     4  +coeff( 67)*x12*x21    *x41    
     5  +coeff( 68)    *x21*x32*x42    
     6  +coeff( 69)    *x21*x31*x42*x51
     7  +coeff( 70)    *x21    *x43*x51
     8  +coeff( 71)*x11*x22*x31*x41    
      x_sp_col    =x_sp_col    
     9  +coeff( 72)    *x24*x31*x41    
     1  +coeff( 73)*x11*x22    *x42    
     2  +coeff( 74)        *x33*x42    
     3  +coeff( 75)    *x24*x31    *x51
     4  +coeff( 76)    *x23*x31*x42    
     5  +coeff( 77)    *x23    *x43    
     6  +coeff( 78)*x11*x24*x31        
     7  +coeff( 79)*x11    *x33*x41    
     8  +coeff( 80)*x11    *x31*x43    
      x_sp_col    =x_sp_col    
     9  +coeff( 81)    *x22*x33    *x51
     1  +coeff( 82)*x11    *x31*x41*x52
     2  +coeff( 83)*x11        *x42*x52
     3  +coeff( 84)*x11*x23*x31*x41    
     4  +coeff( 85)*x11*x22*x32*x41    
     5  +coeff( 86)*x12        *x43    
     6  +coeff( 87)*x11*x22    *x43    
     7  +coeff( 88)*x11*x22    *x42*x51
     8  +coeff( 89)*x12*x23*x31        
      x_sp_col    =x_sp_col    
     9  +coeff( 90)    *x23*x33*x41    
c
      return
      end
      function t_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 36)
      data ncoeff/ 35/
      data avdat/ -0.1197655E-02/
      data xmin/
     1 -0.49991E-02,-0.51642E-01,-0.19979E-01,-0.27954E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.44580273E-03, 0.45200732E-01,-0.54941308E-02, 0.86960726E-03,
     + -0.28106126E-02, 0.19448526E-02,-0.39623649E-03, 0.37828812E-02,
     + -0.63787313E-03, 0.29004968E-03, 0.11645001E-04, 0.11453175E-02,
     +  0.44331534E-03, 0.19118312E-02, 0.88220941E-05,-0.32418815E-03,
     +  0.10698733E-03, 0.54767699E-03,-0.18840941E-02,-0.17891770E-02,
     +  0.36126268E-03, 0.29017687E-04,-0.31841121E-03, 0.15123433E-03,
     +  0.13196243E-03,-0.17242179E-04, 0.60401385E-03, 0.18503099E-02,
     +  0.18127330E-02, 0.46188634E-05,-0.13539138E-03,-0.44792853E-03,
     +  0.61015424E-04, 0.66311492E-04, 0.23600231E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x23            
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x23    *x41    
      t_sp_col    =t_sp_col    
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x22*x33        
     7  +coeff( 16)        *x31        
     8  +coeff( 17)*x11*x21            
      t_sp_col    =t_sp_col    
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x23*x32        
     4  +coeff( 22)                *x51
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11    *x33        
      t_sp_col    =t_sp_col    
     9  +coeff( 27)*x11*x22    *x41    
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)*x11*x22*x33        
     4  +coeff( 31)*x11    *x31        
     5  +coeff( 32)    *x21*x32        
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)*x11*x21*x31        
     8  +coeff( 35)*x11*x22*x31        
c
      return
      end
      function y_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.1003641E-01/
      data xmin/
     1 -0.49991E-02,-0.51642E-01,-0.19979E-01,-0.27954E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.52737882E-02, 0.20815948E-01, 0.65020688E-01, 0.46378723E-02,
     + -0.13932830E-02,-0.40121782E-02,-0.44864428E-03,-0.25729337E-02,
     +  0.97783294E-03,-0.12216633E-03,-0.24410311E-03,-0.54334465E-03,
     + -0.21397347E-04,-0.12370078E-02,-0.44584883E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22    *x41    
      y_sp_col    =y_sp_col    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22*x33        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)        *x33        
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x21*x33        
c
      return
      end
      function p_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4542711E-02/
      data xmin/
     1 -0.49991E-02,-0.51642E-01,-0.19979E-01,-0.27954E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.31281598E-02,-0.16825011E-02, 0.12383327E-02, 0.28635845E-01,
     +  0.65813903E-02,-0.45697456E-02,-0.53972984E-03, 0.13653857E-02,
     + -0.28262297E-02,-0.80694335E-04,-0.12674473E-02,-0.34787238E-03,
     + -0.14469721E-02, 0.51751465E-03,-0.62588387E-03, 0.10926946E-02,
     +  0.10021763E-02, 0.26649630E-03,-0.44176122E-03, 0.12179316E-02,
     + -0.22346398E-02,-0.24406484E-02,-0.91977490E-04, 0.24574660E-03,
     +  0.43269600E-04,-0.51246148E-04,-0.38899280E-03,-0.63527358E-03,
     + -0.54383301E-03,-0.31583695E-03, 0.21482425E-03,-0.12546682E-03,
     + -0.54974429E-03, 0.43945012E-03, 0.54427917E-03,-0.35859953E-03,
     + -0.46768422E-04,-0.31299824E-04,-0.70987924E-04,-0.18506647E-03,
     +  0.74425465E-04, 0.52409800E-04,-0.48919374E-04, 0.42283817E-04,
     +  0.15641954E-03,-0.10566545E-03,-0.19171697E-03, 0.57039990E-04,
     + -0.13645674E-03, 0.12814770E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x23            
      p_sp_col    =p_sp_col    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x22*x33        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)            *x42    
      p_sp_col    =p_sp_col    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11                
     6  +coeff( 24)        *x32        
     7  +coeff( 25)    *x21        *x51
     8  +coeff( 26)            *x41*x51
      p_sp_col    =p_sp_col    
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x21    *x42    
     2  +coeff( 29)        *x31*x42    
     3  +coeff( 30)            *x43    
     4  +coeff( 31)*x11*x22            
     5  +coeff( 32)*x11*x21*x31        
     6  +coeff( 33)    *x22*x32        
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)*x11*x23    *x43    
      p_sp_col    =p_sp_col    
     9  +coeff( 36)*x12*x22*x32    *x51
     1  +coeff( 37)        *x31    *x51
     2  +coeff( 38)*x12                
     3  +coeff( 39)    *x21*x32        
     4  +coeff( 40)        *x32*x41    
     5  +coeff( 41)    *x21    *x41*x51
     6  +coeff( 42)    *x21        *x52
     7  +coeff( 43)*x11    *x32        
     8  +coeff( 44)*x11*x21        *x51
      p_sp_col    =p_sp_col    
     9  +coeff( 45)    *x21*x32*x41    
     1  +coeff( 46)    *x21    *x43    
     2  +coeff( 47)        *x31*x43    
     3  +coeff( 48)    *x22*x31    *x51
     4  +coeff( 49)    *x22    *x41*x51
     5  +coeff( 50)    *x21    *x42*x51
c
      return
      end
      function l_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3490818E-03/
      data xmin/
     1 -0.49991E-02,-0.51642E-01,-0.19979E-01,-0.27954E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.79817994E-03,-0.26572333E-02,-0.58851512E-02,-0.30263294E-02,
     + -0.55855606E-04,-0.89377048E-03, 0.12682959E-04,-0.13026781E-03,
     + -0.12646907E-03, 0.33594429E-03, 0.62282925E-04, 0.12035081E-03,
     +  0.12458320E-03, 0.20915733E-04,-0.70904010E-04,-0.13070719E-03,
     +  0.41082407E-04, 0.39301314E-04,-0.48375838E-04, 0.16690238E-03,
     + -0.21350266E-04,-0.37381981E-05,-0.11655752E-04, 0.13568913E-03,
     + -0.49826831E-04, 0.14140662E-05, 0.55767123E-04, 0.99453682E-05,
     +  0.10790685E-04, 0.66120215E-05, 0.44141448E-05, 0.13810441E-04,
     + -0.22704004E-04,-0.22417056E-04, 0.39715360E-05, 0.13940682E-04,
     + -0.92134396E-05,-0.93464014E-05,-0.39437632E-05,-0.17440661E-04,
     + -0.61751060E-04, 0.16231979E-04, 0.25418898E-04,-0.10791795E-04,
     +  0.17759892E-04, 0.25984233E-04, 0.21603120E-04,-0.27420525E-04,
     + -0.16236303E-04, 0.28733166E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22        *x51
     8  +coeff(  8)                *x51
      l_sp_col    =l_sp_col    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21            
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x31        
      l_sp_col    =l_sp_col    
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x22*x33*x41    
     4  +coeff( 22)                *x52
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)    *x22*x31*x41    
     7  +coeff( 25)*x11*x23    *x43    
     8  +coeff( 26)*x11                
      l_sp_col    =l_sp_col    
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)        *x31*x42    
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)*x11    *x32        
     5  +coeff( 32)    *x22*x32        
     6  +coeff( 33)    *x21*x32*x41    
     7  +coeff( 34)        *x32*x42    
     8  +coeff( 35)    *x23        *x51
      l_sp_col    =l_sp_col    
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)*x11*x22    *x41    
     2  +coeff( 38)*x11    *x31*x41*x51
     3  +coeff( 39)    *x23*x32        
     4  +coeff( 40)    *x23*x31*x41    
     5  +coeff( 41)    *x23    *x42    
     6  +coeff( 42)    *x21*x31*x43    
     7  +coeff( 43)    *x21*x31*x41*x52
     8  +coeff( 44)*x11    *x31*x42*x51
      l_sp_col    =l_sp_col    
     9  +coeff( 45)*x12*x21*x31*x41    
     1  +coeff( 46)    *x21*x32*x43    
     2  +coeff( 47)    *x21*x31*x42*x52
     3  +coeff( 48)*x11*x23    *x42    
     4  +coeff( 49)*x11*x22*x31    *x52
     5  +coeff( 50)*x11*x21    *x42*x52
c
      return
      end
      function x_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3243843E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35046928E-02, 0.13820717E+00, 0.21235705E-02,-0.10606142E-01,
     +  0.34711678E-02, 0.35224396E-02, 0.26484472E-02,-0.54140864E-02,
     +  0.35459842E-02,-0.17646390E-02, 0.22064387E-04, 0.70822774E-02,
     + -0.92774438E-03, 0.47117021E-03, 0.53146742E-02, 0.13098646E-02,
     + -0.16762203E-02, 0.10409401E-03, 0.26683994E-02,-0.70732972E-03,
     + -0.52500917E-02,-0.52432134E-02,-0.11920842E-02, 0.17146377E-02,
     +  0.13960109E-02,-0.22268046E-05,-0.74056850E-03,-0.72990055E-03,
     +  0.13076246E-03,-0.13256518E-02,-0.15625122E-03, 0.61060576E-03,
     +  0.23306895E-02, 0.29354915E-02,-0.47221460E-03, 0.73501829E-03,
     + -0.34601018E-02, 0.63010734E-02, 0.60825166E-02,-0.63876942E-03,
     +  0.41855710E-04, 0.31078144E-03,-0.34170644E-04,-0.18637627E-03,
     + -0.34944245E-03,-0.12651231E-04, 0.20968294E-03, 0.59067964E-03,
     +  0.83787912E-04, 0.16964776E-03, 0.57608949E-03,-0.81158103E-03,
     +  0.10163980E-04,-0.68672482E-04, 0.13967622E-04,-0.39445418E-04,
     + -0.11889340E-02,-0.19764120E-03, 0.46745929E-04, 0.82616454E-04,
     +  0.22626195E-04, 0.41680768E-03, 0.18111030E-03,-0.56259532E-03,
     +  0.51449111E-03,-0.10328821E-02, 0.42032171E-02, 0.38268734E-03,
     + -0.30964628E-03,-0.13040742E-02, 0.16168610E-03,-0.61785424E-03,
     + -0.66962332E-03,-0.39734092E-03, 0.89981710E-04, 0.20007046E-03,
     +  0.37591203E-03, 0.27817930E-02, 0.10792819E-02, 0.16222079E-03,
     +  0.32944838E-03, 0.22416514E-03, 0.20585577E-03, 0.26865592E-02,
     + -0.22150436E-02, 0.12300088E-03, 0.39742165E-03,-0.68522466E-03,
     +  0.44488319E-03,-0.16703029E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x41    
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)        *x31        
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x24*x31        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 18)                *x51
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)    *x21    *x42    
     5  +coeff( 23)    *x24            
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)    *x23*x32        
     8  +coeff( 26)*x11    *x33        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 27)        *x31*x41    
     1  +coeff( 28)            *x42    
     2  +coeff( 29)*x11            *x51
     3  +coeff( 30)    *x21*x32        
     4  +coeff( 31)    *x21        *x52
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)    *x22    *x42    
     8  +coeff( 35)*x11*x23            
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 36)*x11*x22*x31        
     1  +coeff( 37)    *x24    *x41    
     2  +coeff( 38)    *x23*x31*x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)*x11*x24            
     5  +coeff( 41)    *x24*x32        
     6  +coeff( 42)*x11*x21*x33        
     7  +coeff( 43)*x11*x22*x33        
     8  +coeff( 44)        *x32        
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 45)*x11    *x31        
     1  +coeff( 46)    *x21    *x41*x51
     2  +coeff( 47)*x11*x21*x31        
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x23*x31    *x51
     6  +coeff( 51)    *x23    *x41*x51
     7  +coeff( 52)*x11*x23*x33        
     8  +coeff( 53)            *x41*x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 54)    *x22        *x51
     1  +coeff( 55)*x12                
     2  +coeff( 56)    *x23        *x51
     3  +coeff( 57)    *x21    *x43    
     4  +coeff( 58)*x11*x22        *x51
     5  +coeff( 59)        *x33    *x51
     6  +coeff( 60)        *x32*x41*x51
     7  +coeff( 61)        *x31*x41*x52
     8  +coeff( 62)*x11*x21*x31*x41    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 63)*x12*x21    *x41    
     1  +coeff( 64)*x11*x23    *x41    
     2  +coeff( 65)    *x21*x31*x43    
     3  +coeff( 66)    *x24    *x42    
     4  +coeff( 67)    *x23    *x43    
     5  +coeff( 68)*x11*x24        *x51
     6  +coeff( 69)*x11*x23*x31*x41    
     7  +coeff( 70)    *x21*x33*x42    
     8  +coeff( 71)*x12*x21*x31    *x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 72)*x11*x23    *x41*x51
     1  +coeff( 73)    *x21*x31*x42*x52
     2  +coeff( 74)    *x21    *x43*x52
     3  +coeff( 75)*x12    *x32    *x51
     4  +coeff( 76)*x12*x23*x31        
     5  +coeff( 77)*x11*x21*x32*x42    
     6  +coeff( 78)    *x23*x32*x42    
     7  +coeff( 79)    *x23*x31*x43    
     8  +coeff( 80)*x11*x21*x32*x41*x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 81)*x11*x21    *x43*x51
     1  +coeff( 82)    *x22*x31*x42*x52
     2  +coeff( 83)*x12*x24*x31        
     3  +coeff( 84)    *x23*x33*x42    
     4  +coeff( 85)    *x23*x32*x43    
     5  +coeff( 86)*x11*x21*x32    *x53
     6  +coeff( 87)    *x23    *x42*x53
     7  +coeff( 88)*x11*x24*x32*x41    
     8  +coeff( 89)*x12*x22*x31*x41*x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 90)*x11    *x33*x41*x52
c
      return
      end
      function t_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.9797944E-04/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.40009449E-03,-0.11608735E-01,-0.20338397E-02,-0.23080267E-02,
     +  0.27984958E-02, 0.83215121E-03, 0.28170069E-03,-0.11615301E-02,
     +  0.44864279E-03,-0.29523356E-03, 0.89663008E-04, 0.26198491E-05,
     +  0.19150944E-02, 0.62392191E-05,-0.14769388E-03, 0.33012457E-04,
     +  0.39820799E-04, 0.60589134E-03,-0.11115954E-02,-0.20476615E-03,
     + -0.14961119E-03, 0.12159404E-03, 0.81593298E-05, 0.31266178E-03,
     + -0.45774264E-04, 0.18358114E-03,-0.12029395E-03, 0.27182867E-03,
     + -0.12972092E-03,-0.78953546E-03, 0.11224839E-03,-0.80590362E-05,
     +  0.29877588E-03,-0.53885211E-04, 0.36459074E-04,-0.51752191E-04,
     +  0.70518450E-04, 0.43465127E-03, 0.11481462E-03, 0.21840294E-05,
     +  0.33808534E-03,-0.40093175E-04, 0.16406924E-04,-0.89121415E-04,
     + -0.10783842E-03, 0.88532706E-05, 0.28014034E-04,-0.23477829E-04,
     +  0.19806463E-04, 0.32025063E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x41    
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)        *x33        
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)    *x21    *x42*x52
     6  +coeff( 15)        *x31        
     7  +coeff( 16)                *x51
     8  +coeff( 17)*x11*x21            
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)    *x22*x33        
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)    *x21*x33*x41    
     8  +coeff( 26)    *x23*x33*x41    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)    *x21*x31*x41    
     4  +coeff( 31)*x11*x22*x31        
     5  +coeff( 32)*x11    *x33        
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)*x11    *x31        
     8  +coeff( 35)    *x22        *x51
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)*x11*x21    *x41    
     2  +coeff( 38)    *x22    *x42    
     3  +coeff( 39)    *x23        *x51
     4  +coeff( 40)            *x42*x52
     5  +coeff( 41)    *x23    *x42    
     6  +coeff( 42)    *x22*x33*x41    
     7  +coeff( 43)        *x33*x41*x52
     8  +coeff( 44)        *x31*x41    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 45)            *x42    
     1  +coeff( 46)*x12                
     2  +coeff( 47)*x11*x21*x31        
     3  +coeff( 48)*x11        *x42    
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x22*x31*x41    
c
      return
      end
      function y_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/ -0.1149580E-01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16589741E-02, 0.32782622E-01, 0.15225032E+00, 0.20286979E-01,
     + -0.59838700E-02, 0.42002425E-02,-0.14641795E-01,-0.18742003E-02,
     + -0.18368568E-02,-0.12008469E-02,-0.32101199E-02,-0.10179684E-01,
     +  0.45564543E-02, 0.65029995E-03, 0.37975307E-03, 0.33783622E-02,
     + -0.19715137E-02, 0.29082878E-04, 0.49081608E-03,-0.51829917E-02,
     + -0.12964390E-01,-0.30582938E-02,-0.29139865E-01, 0.24116378E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)*x11*x21            
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x22*x33        
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)    *x21*x31        
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff( 18)        *x33        
     1  +coeff( 19)        *x34        
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x22*x31*x42    
     5  +coeff( 23)    *x22*x33*x41    
     6  +coeff( 24)    *x22*x35*x41    
c
      return
      end
      function p_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4893997E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.98830741E-03,-0.28498601E-02, 0.10413361E-01, 0.64244770E-01,
     +  0.10496110E-01,-0.71069663E-02,-0.18439719E-02,-0.11740138E-02,
     + -0.36603806E-02,-0.24039480E-02,-0.76332199E-03, 0.20071079E-02,
     + -0.20787118E-02, 0.14187178E-02, 0.27823778E-05,-0.13041527E-02,
     +  0.22998715E-02, 0.23600121E-02,-0.50522725E-03, 0.57531771E-03,
     + -0.91710925E-03,-0.14570614E-02,-0.87013459E-02,-0.67322440E-02,
     +  0.34830711E-03,-0.19559430E-03, 0.52467152E-03, 0.15816269E-03,
     +  0.37313451E-03,-0.26036642E-03, 0.25602879E-02,-0.70006126E-02,
     + -0.13064522E-02, 0.60269440E-03,-0.44398829E-02,-0.66051311E-04,
     + -0.10293686E-02,-0.17186401E-04,-0.17188201E-02, 0.29935312E-03,
     +  0.98449258E-04,-0.25855942E-03,-0.30312777E-03,-0.47559585E-03,
     + -0.28410321E-02, 0.15723536E-02, 0.46326409E-03, 0.54622366E-03,
     +  0.31804823E-03,-0.28398263E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)*x11*x21            
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)        *x33    *x51
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)        *x31*x41    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 18)            *x42    
     1  +coeff( 19)        *x31    *x51
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x22*x32        
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x22*x31*x42    
     7  +coeff( 25)    *x22*x33*x41    
     8  +coeff( 26)*x11                
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 27)        *x32        
     1  +coeff( 28)    *x21        *x51
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x22*x31*x41    
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)*x11*x23            
     8  +coeff( 35)    *x22    *x43    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 36)*x11        *x41    
     1  +coeff( 37)    *x21*x31*x41    
     2  +coeff( 38)        *x32*x41    
     3  +coeff( 39)    *x21    *x42    
     4  +coeff( 40)        *x31*x42    
     5  +coeff( 41)    *x21        *x52
     6  +coeff( 42)        *x31*x43    
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)    *x23*x32        
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 45)    *x22*x32*x41    
     1  +coeff( 46)    *x23    *x42    
     2  +coeff( 47)        *x32*x43    
     3  +coeff( 48)*x11*x21*x32*x41    
     4  +coeff( 49)*x11*x22*x33        
     5  +coeff( 50)*x12*x23*x31        
c
      return
      end
      function l_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1093464E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15361847E-02,-0.26322296E-02,-0.56418604E-02,-0.35600334E-02,
     + -0.45516681E-04,-0.72036864E-03,-0.32229065E-02,-0.11585718E-03,
     + -0.94768056E-03, 0.10188132E-02,-0.26597252E-04,-0.87978115E-04,
     +  0.37278078E-03,-0.82791739E-04, 0.24148171E-03,-0.78931553E-04,
     +  0.13242166E-03,-0.97694603E-04, 0.75646589E-03, 0.29982322E-04,
     +  0.35679666E-04, 0.11398840E-03, 0.80716192E-04, 0.10416959E-03,
     + -0.33576231E-03, 0.45343090E-03, 0.90637324E-04, 0.27373183E-03,
     +  0.40557232E-04,-0.16564991E-04, 0.21861833E-04, 0.40160205E-04,
     +  0.11267841E-03,-0.21844424E-03,-0.99935442E-05,-0.35491955E-03,
     +  0.46055960E-04,-0.24750937E-05, 0.15115266E-04, 0.16026416E-03,
     + -0.24733370E-04, 0.11410602E-04, 0.12315980E-04, 0.18795827E-04,
     +  0.89923333E-05, 0.37505808E-04, 0.42650023E-04,-0.31150343E-04,
     + -0.16325655E-04,-0.35999517E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21            
     3  +coeff( 12)                *x51
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22        *x51
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x21        *x51
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x22*x31*x41    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 27)    *x21*x31        
     1  +coeff( 28)    *x21    *x42    
     2  +coeff( 29)        *x31*x41*x51
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)        *x31*x43    
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x21*x33*x41    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x21*x31*x41*x52
     2  +coeff( 38)*x11    *x31        
     3  +coeff( 39)*x11        *x41    
     4  +coeff( 40)    *x21*x31*x41    
     5  +coeff( 41)            *x43    
     6  +coeff( 42)        *x31    *x52
     7  +coeff( 43)                *x53
     8  +coeff( 44)*x11*x21*x31        
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 45)*x11    *x32        
     1  +coeff( 46)    *x22*x32        
     2  +coeff( 47)    *x22*x31    *x51
     3  +coeff( 48)*x11*x23            
     4  +coeff( 49)*x11*x21*x32        
     5  +coeff( 50)*x11*x22    *x41    
c
      return
      end
      function x_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5811551E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.43507488E-02, 0.18971764E+00,-0.40765321E-02,-0.26228964E-01,
     +  0.11775252E-01, 0.87674903E-02, 0.60311928E-02,-0.13326469E-01,
     +  0.73981616E-02,-0.42698686E-02, 0.74378506E-04, 0.18584013E-01,
     + -0.21921780E-02, 0.13039014E-01,-0.14129227E-01, 0.29321981E-02,
     + -0.38474507E-02, 0.14107072E-01, 0.31267849E-03, 0.10889688E-02,
     +  0.63126986E-02,-0.17044984E-02, 0.56004414E-03,-0.12914772E-01,
     + -0.84679894E-03,-0.27003135E-02, 0.40978701E-02, 0.30937344E-02,
     +  0.13418450E-01,-0.19659830E-04,-0.31406030E-02, 0.14766543E-02,
     + -0.10892390E-02, 0.18693525E-02,-0.83690286E-02, 0.64435400E-04,
     + -0.16456248E-02,-0.16882342E-02,-0.86346088E-03, 0.21630267E-03,
     + -0.26037954E-03, 0.50659146E-03, 0.13591120E-02, 0.52987230E-02,
     +  0.65221381E-02, 0.18037928E-03, 0.58870490E-04,-0.67150882E-04,
     + -0.15467937E-02,-0.18120933E-02,-0.39999586E-03, 0.53807591E-04,
     +  0.41522591E-04,-0.34260610E-03, 0.48111731E-03,-0.40150134E-03,
     +  0.37143531E-03,-0.13380445E-02, 0.59880194E-03,-0.15333528E-02,
     +  0.64886769E-03, 0.68590245E-02,-0.45582812E-04, 0.87347097E-03,
     + -0.12604342E-02, 0.39688282E-03, 0.91499649E-02, 0.77834157E-02,
     + -0.81980403E-03,-0.14778236E-02, 0.58365357E-02, 0.12287355E-02,
     +  0.37264227E-03, 0.19784255E-02,-0.16266836E-02, 0.15291426E-03,
     +  0.10591101E-03, 0.11920339E-03,-0.23899855E-04, 0.26951058E-03,
     + -0.21636298E-03,-0.47041420E-04,-0.53946528E-03,-0.29838469E-02,
     + -0.31029633E-02,-0.26534830E-03, 0.15766179E-03,-0.12009691E-03,
     + -0.12576311E-03, 0.29365878E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x41    
     2  +coeff( 11)        *x33        
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)        *x31        
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x24*x31        
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)                *x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)    *x24            
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 27)*x11*x22    *x41    
     1  +coeff( 28)    *x23*x32        
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)*x11    *x33        
     4  +coeff( 31)    *x21*x32        
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)*x11*x22*x31        
     8  +coeff( 35)    *x24    *x41    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 36)*x11*x22*x33        
     1  +coeff( 37)        *x31*x41    
     2  +coeff( 38)            *x42    
     3  +coeff( 39)*x11    *x31        
     4  +coeff( 40)    *x21*x31    *x51
     5  +coeff( 41)    *x21    *x41*x51
     6  +coeff( 42)*x11*x21*x31        
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)    *x22*x31*x41    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 45)    *x22    *x42    
     1  +coeff( 46)*x12*x21            
     2  +coeff( 47)        *x32*x42    
     3  +coeff( 48)        *x32    *x52
     4  +coeff( 49)*x11*x24            
     5  +coeff( 50)*x11*x23*x33        
     6  +coeff( 51)        *x32        
     7  +coeff( 52)    *x22        *x51
     8  +coeff( 53)*x12                
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 54)*x11        *x42    
     1  +coeff( 55)    *x21    *x42*x51
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)*x12*x21    *x41    
     4  +coeff( 58)*x11*x23    *x41    
     5  +coeff( 59)    *x21*x31*x43    
     6  +coeff( 60)    *x24    *x42    
     7  +coeff( 61)*x11*x21*x33        
     8  +coeff( 62)    *x23    *x43    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 63)    *x23*x31    *x52
     1  +coeff( 64)*x11*x24        *x51
     2  +coeff( 65)*x11*x23    *x41*x51
     3  +coeff( 66)*x12*x23*x31        
     4  +coeff( 67)    *x23*x32*x42    
     5  +coeff( 68)    *x23*x31*x43    
     6  +coeff( 69)*x11*x21*x31*x42*x51
     7  +coeff( 70)    *x24    *x42*x52
     8  +coeff( 71)    *x23*x33*x42    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 72)    *x23*x32*x41*x52
     1  +coeff( 73)*x11*x21*x32    *x53
     2  +coeff( 74)    *x23*x33*x43    
     3  +coeff( 75)*x12*x24*x32*x41    
     4  +coeff( 76)        *x31    *x51
     5  +coeff( 77)            *x41*x51
     6  +coeff( 78)        *x32*x41    
     7  +coeff( 79)        *x31    *x52
     8  +coeff( 80)    *x23        *x51
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 81)*x11    *x31*x41    
     1  +coeff( 82)*x11            *x52
     2  +coeff( 83)    *x21*x32*x41    
     3  +coeff( 84)    *x21*x31*x42    
     4  +coeff( 85)    *x21    *x43    
     5  +coeff( 86)    *x24        *x51
     6  +coeff( 87)        *x32*x41*x51
     7  +coeff( 88)        *x31    *x53
     8  +coeff( 89)            *x41*x53
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 90)*x11*x21*x31*x41    
c
      return
      end
      function t_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1714351E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10224550E-02, 0.47820229E-01,-0.23766253E-02,-0.90404125E-02,
     +  0.31167238E-02, 0.32885966E-02, 0.13326252E-02,-0.46246108E-02,
     +  0.23313819E-02, 0.77936340E-04,-0.65245468E-03,-0.12966688E-02,
     +  0.71793380E-02, 0.19102865E-03, 0.21860279E-03, 0.12386366E-02,
     +  0.25317855E-02,-0.50783572E-02,-0.38105497E-03,-0.66174762E-05,
     +  0.64052321E-03, 0.22652901E-02,-0.58255700E-03, 0.17597625E-02,
     +  0.11852824E-03,-0.52833837E-03,-0.93936757E-03,-0.36690328E-02,
     + -0.18455727E-03, 0.48317033E-03,-0.42085714E-04, 0.12173961E-02,
     +  0.36997537E-02,-0.33890823E-03,-0.46517904E-03,-0.22954337E-03,
     +  0.23141199E-03, 0.31937956E-03, 0.95014388E-04, 0.18457911E-02,
     +  0.82093332E-03,-0.14326711E-03, 0.14590102E-02, 0.17669071E-02,
     + -0.22361509E-02, 0.40775842E-04, 0.33767203E-04, 0.16869952E-03,
     +  0.50485192E-04, 0.19661586E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x43    
     2  +coeff( 11)        *x31        
     3  +coeff( 12)            *x41    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x22    *x41    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)                *x53
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)    *x21*x33*x41    
     6  +coeff( 24)    *x23*x33*x41    
     7  +coeff( 25)                *x51
     8  +coeff( 26)*x11        *x41    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)    *x21    *x41*x51
     3  +coeff( 30)*x11*x22*x31        
     4  +coeff( 31)*x11    *x33        
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)    *x23    *x42    
     7  +coeff( 34)        *x31*x41    
     8  +coeff( 35)            *x42    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 36)*x11    *x31        
     1  +coeff( 37)        *x31*x42    
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)*x12*x21            
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)    *x23*x32        
     6  +coeff( 42)*x11*x21*x33        
     7  +coeff( 43)    *x22*x33*x41    
     8  +coeff( 44)    *x22*x31*x41*x52
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 45)    *x22*x33*x41*x52
     1  +coeff( 46)            *x41*x51
     2  +coeff( 47)*x12                
     3  +coeff( 48)        *x32*x41    
     4  +coeff( 49)    *x21*x31    *x51
     5  +coeff( 50)*x11*x21*x31        
c
      return
      end
      function y_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 23)
      data ncoeff/ 22/
      data avdat/ -0.1542621E-01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26910913E-02, 0.36815662E-01, 0.20537122E+00, 0.31306189E-01,
     + -0.92907827E-02, 0.73181405E-02, 0.87267370E-03,-0.21872345E-01,
     + -0.28610476E-02,-0.14282908E-01, 0.73219403E-02,-0.12186230E-02,
     + -0.46563642E-02,-0.71645752E-02,-0.44527254E-03,-0.25407584E-01,
     +  0.65569910E-02,-0.22043148E-02,-0.52305206E-03, 0.17102008E-02,
     + -0.19902330E-01,-0.45669232E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22            
      y_sp_q2ex   =y_sp_q2ex   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x21*x33        
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)        *x31*x41    
      y_sp_q2ex   =y_sp_q2ex   
     9  +coeff( 18)    *x21*x31        
     1  +coeff( 19)        *x31*x42    
     2  +coeff( 20)        *x34        
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22*x34        
c
      return
      end
      function p_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2299192E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.28446689E-03,-0.77269683E-02,-0.28859049E-01,-0.31830885E-02,
     +  0.22089083E-02, 0.78678248E-03, 0.40223910E-02, 0.76257647E-03,
     + -0.67707751E-03, 0.91540680E-03,-0.63519739E-03, 0.73706475E-03,
     + -0.22616735E-03, 0.36548171E-03,-0.50348478E-06, 0.69359934E-03,
     +  0.84768084E-03, 0.72416085E-04, 0.28196040E-02, 0.63799183E-04,
     +  0.38350944E-03,-0.61767123E-03, 0.70956441E-04,-0.19813728E-03,
     + -0.22341189E-03,-0.13287488E-03, 0.23523069E-03,-0.60593092E-03,
     +  0.21576867E-02, 0.62330430E-04, 0.81335958E-04, 0.28497435E-03,
     +  0.18842531E-03,-0.13653377E-03,-0.10252730E-04, 0.17894083E-04,
     + -0.99094221E-04,-0.14918791E-03,-0.72186877E-03, 0.22180173E-04,
     +  0.10349078E-03,-0.82031642E-04,-0.13679519E-03,-0.42233380E-04,
     +  0.11298771E-03,-0.62721439E-04,-0.26271207E-03, 0.40415881E-03,
     + -0.90375761E-04,-0.19025756E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31    *x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)                *x52
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21            
     2  +coeff( 11)    *x22        *x51
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x21*x33        
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x22    *x41    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 18)        *x31*x42    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)*x11                
     3  +coeff( 21)    *x21*x31        
     4  +coeff( 22)            *x42    
     5  +coeff( 23)        *x32*x41    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)*x11*x22            
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)        *x33*x41    
     4  +coeff( 31)        *x32*x42    
     5  +coeff( 32)        *x31*x43    
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)            *x42*x53
     8  +coeff( 35)*x11*x23*x31        
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 36)    *x23*x33        
     1  +coeff( 37)        *x33*x43    
     2  +coeff( 38)        *x32        
     3  +coeff( 39)        *x31*x41    
     4  +coeff( 40)*x11        *x41    
     5  +coeff( 41)    *x21    *x42    
     6  +coeff( 42)    *x21*x31    *x51
     7  +coeff( 43)    *x21    *x41*x51
     8  +coeff( 44)                *x53
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 45)*x11*x21*x31        
     1  +coeff( 46)*x11*x21        *x51
     2  +coeff( 47)    *x23*x31        
     3  +coeff( 48)    *x22*x32        
     4  +coeff( 49)    *x22*x31    *x51
     5  +coeff( 50)*x11*x23            
c
      return
      end
      function l_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1788704E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.51759E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22303602E-02,-0.96968806E-05,-0.25890521E-02,-0.53465017E-02,
     + -0.48884712E-02,-0.95362659E-04,-0.14685178E-02,-0.55307881E-02,
     + -0.23379177E-03,-0.17075357E-02,-0.15131859E-03, 0.16635378E-02,
     +  0.31626852E-04, 0.60905691E-03, 0.12008159E-03,-0.12097808E-03,
     +  0.44698155E-03, 0.29042963E-03, 0.12228634E-02, 0.15059886E-03,
     +  0.48858863E-04, 0.16630710E-03, 0.14139136E-03,-0.56086713E-03,
     +  0.74078864E-03, 0.44194845E-03, 0.56614459E-04,-0.10852234E-04,
     +  0.11057209E-03,-0.20229530E-04, 0.43022425E-04,-0.17511616E-03,
     +  0.19567547E-03,-0.35242952E-03, 0.72275529E-04,-0.41052066E-04,
     +  0.23175668E-04, 0.26857932E-03,-0.13220111E-03,-0.20057332E-04,
     +  0.27248103E-04, 0.21489222E-04, 0.43433891E-04,-0.41319548E-04,
     +  0.61953528E-04,-0.46932735E-04,-0.32103690E-04,-0.61843166E-03,
     + -0.28408629E-04, 0.48085686E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21    *x43    
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22*x31        
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)    *x21    *x42    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)        *x32    *x51
     2  +coeff( 29)        *x31*x41*x51
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)    *x23*x31        
     6  +coeff( 33)    *x22    *x41*x51
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x22*x32    *x51
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 36)                *x51
     1  +coeff( 37)*x11        *x41    
     2  +coeff( 38)    *x21*x31*x41    
     3  +coeff( 39)            *x43    
     4  +coeff( 40)    *x21        *x52
     5  +coeff( 41)        *x31    *x52
     6  +coeff( 42)                *x53
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)        *x32*x42    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 45)    *x22*x31    *x51
     1  +coeff( 46)*x11*x22    *x41    
     2  +coeff( 47)*x11    *x31*x41*x51
     3  +coeff( 48)    *x23    *x42    
     4  +coeff( 49)        *x33*x42    
     5  +coeff( 50)    *x22    *x43    
c
      return
      end
      function x_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5095685E+01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18589617E-02,-0.11099883E+00, 0.39965459E-02, 0.18558145E-01,
     + -0.70644896E-02,-0.61902008E-02, 0.93237972E-02,-0.48386175E-02,
     +  0.26913974E-02,-0.91350824E-02,-0.13517370E-01, 0.13813346E-02,
     + -0.92623092E-03,-0.44709807E-02,-0.39949163E-03, 0.99946149E-02,
     + -0.19554135E-02,-0.72594522E-02, 0.16068591E-05, 0.23751540E-03,
     + -0.36286428E-02, 0.11963994E-02, 0.20079112E-02, 0.90966877E-02,
     +  0.21246050E-02,-0.41875592E-02,-0.27915682E-02,-0.71739922E-02,
     +  0.54922908E-04,-0.18342678E-03,-0.10563035E-02, 0.10147621E-02,
     +  0.92251122E-03, 0.49305073E-03,-0.18077231E-04, 0.56215271E-03,
     + -0.94626227E-03,-0.33418906E-02, 0.70969138E-03,-0.11156291E-02,
     +  0.59451750E-02,-0.17376698E-02, 0.99769386E-03, 0.43373773E-03,
     + -0.21828088E-03, 0.26361190E-03,-0.10258287E-03, 0.57125645E-03,
     + -0.20574182E-03,-0.57382561E-03,-0.84669743E-03,-0.11225059E-03,
     +  0.26989048E-02, 0.89098717E-03,-0.76337885E-04, 0.22893641E-03,
     + -0.17741128E-03,-0.10034630E-03,-0.15362604E-03,-0.38083471E-03,
     + -0.10353556E-02, 0.76719699E-03, 0.12414458E-02,-0.29445591E-02,
     + -0.14930581E-03,-0.52581494E-04, 0.26140068E-03, 0.10895098E-03,
     +  0.18962148E-03, 0.32002272E-03, 0.85844443E-03, 0.19099672E-02,
     +  0.13195849E-02,-0.12953515E-03, 0.13178622E-03,-0.18049635E-03,
     + -0.10880282E-01,-0.82227979E-02,-0.53205161E-03,-0.38436253E-03,
     + -0.53196744E-03, 0.21074634E-03, 0.44428010E-03, 0.30591784E-03,
     + -0.21910612E-03,-0.79637265E-03, 0.22787465E-03,-0.17284683E-02,
     + -0.11355936E-02, 0.65236335E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23*x31        
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x23            
      x_sp_den  =x_sp_den  
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)        *x31        
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      x_sp_den  =x_sp_den  
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)    *x21*x33*x41    
     2  +coeff( 20)    *x21*x32    *x52
     3  +coeff( 21)    *x23*x33*x41    
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x32        
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x24            
     8  +coeff( 26)    *x22    *x42    
      x_sp_den  =x_sp_den  
     9  +coeff( 27)*x11*x22    *x41    
     1  +coeff( 28)    *x23    *x42    
     2  +coeff( 29)*x11    *x33        
     3  +coeff( 30)                *x51
     4  +coeff( 31)    *x22            
     5  +coeff( 32)        *x31*x41    
     6  +coeff( 33)            *x42    
     7  +coeff( 34)    *x22        *x51
     8  +coeff( 35)    *x21    *x41*x51
      x_sp_den  =x_sp_den  
     9  +coeff( 36)    *x21        *x52
     1  +coeff( 37)*x11*x21    *x41    
     2  +coeff( 38)    *x22*x31*x41    
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)*x11*x22*x31        
     5  +coeff( 41)    *x24    *x41    
     6  +coeff( 42)    *x23*x32        
     7  +coeff( 43)*x11*x24            
     8  +coeff( 44)*x11*x23*x31        
      x_sp_den  =x_sp_den  
     9  +coeff( 45)*x11*x22*x33        
     1  +coeff( 46)        *x32        
     2  +coeff( 47)            *x41*x51
     3  +coeff( 48)*x11    *x31        
     4  +coeff( 49)    *x21*x31    *x51
     5  +coeff( 50)*x11*x21*x31        
     6  +coeff( 51)    *x22*x32        
     7  +coeff( 52)*x12*x21            
     8  +coeff( 53)    *x24*x31        
      x_sp_den  =x_sp_den  
     9  +coeff( 54)*x11*x23    *x41    
     1  +coeff( 55)        *x31    *x51
     2  +coeff( 56)*x11        *x42    
     3  +coeff( 57)    *x21    *x42*x51
     4  +coeff( 58)*x11*x22        *x51
     5  +coeff( 59)*x11*x21*x31*x41    
     6  +coeff( 60)*x12*x21    *x41    
     7  +coeff( 61)    *x21*x31*x43    
     8  +coeff( 62)    *x21    *x43*x51
      x_sp_den  =x_sp_den  
     9  +coeff( 63)    *x24    *x42    
     1  +coeff( 64)    *x23    *x43    
     2  +coeff( 65)*x11*x21*x32    *x51
     3  +coeff( 66)    *x23*x31    *x52
     4  +coeff( 67)*x11    *x33*x41    
     5  +coeff( 68)*x11    *x32*x42    
     6  +coeff( 69)*x11    *x32    *x52
     7  +coeff( 70)*x11    *x31*x41*x52
     8  +coeff( 71)*x11*x23    *x41*x51
      x_sp_den  =x_sp_den  
     9  +coeff( 72)    *x21*x31*x42*x52
     1  +coeff( 73)    *x21    *x43*x52
     2  +coeff( 74)*x12    *x32    *x51
     3  +coeff( 75)*x12        *x41*x52
     4  +coeff( 76)*x12*x23*x31        
     5  +coeff( 77)    *x23*x32*x42    
     6  +coeff( 78)    *x23*x31*x43    
     7  +coeff( 79)*x11*x21*x32*x41*x51
     8  +coeff( 80)*x11*x21*x31*x41*x52
      x_sp_den  =x_sp_den  
     9  +coeff( 81)*x11*x24*x31*x41    
     1  +coeff( 82)    *x22*x32*x41*x52
     2  +coeff( 83)*x11*x23*x33        
     3  +coeff( 84)*x12*x21*x32*x41    
     4  +coeff( 85)*x12    *x32*x41*x51
     5  +coeff( 86)*x11*x22*x31*x42*x51
     6  +coeff( 87)    *x24*x31    *x53
     7  +coeff( 88)    *x23*x33*x42    
     8  +coeff( 89)    *x23*x32*x42*x51
      x_sp_den  =x_sp_den  
     9  +coeff( 90)*x12*x22    *x42*x51
c
      return
      end
      function t_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1297449E+01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.62028418E-03, 0.48603322E-01, 0.63390721E-03,-0.24074048E-02,
     +  0.47908239E-02,-0.94101662E-02, 0.30957304E-02, 0.33037567E-02,
     + -0.13120968E-02,-0.47748489E-02, 0.24772251E-02,-0.62058686E-03,
     + -0.35765287E-03, 0.50697435E-03, 0.46204170E-03, 0.72999038E-02,
     +  0.22855074E-03,-0.55691809E-02, 0.62510290E-03, 0.28886441E-02,
     + -0.57799515E-03,-0.23066960E-03, 0.17094451E-02,-0.54949062E-03,
     +  0.75093209E-03,-0.93162811E-03, 0.13246220E-02,-0.39479025E-02,
     + -0.23447030E-03,-0.31584763E-03, 0.48228569E-03,-0.23382703E-04,
     +  0.13287198E-02, 0.43698200E-02,-0.24631372E-03, 0.69184578E-04,
     +  0.32490189E-03, 0.25738205E-03, 0.14236138E-02, 0.18490168E-02,
     + -0.60603638E-05, 0.95287530E-03,-0.33347544E-04, 0.61606370E-04,
     +  0.33190569E-04, 0.78975492E-04, 0.55700308E-04, 0.76037315E-04,
     +  0.99164557E-04,-0.65869259E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23*x31        
      t_sp_den  =t_sp_den  
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)        *x31        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)*x11            *x51
      t_sp_den  =t_sp_den  
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x23*x31*x41    
     3  +coeff( 21)    *x21*x33*x41    
     4  +coeff( 22)    *x21*x32    *x52
     5  +coeff( 23)    *x23*x33*x41    
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)    *x21*x32        
      t_sp_den  =t_sp_den  
     9  +coeff( 27)    *x22    *x41    
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)    *x21    *x41*x51
     3  +coeff( 30)    *x21        *x52
     4  +coeff( 31)*x11*x22*x31        
     5  +coeff( 32)*x11    *x33        
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)*x11    *x31        
      t_sp_den  =t_sp_den  
     9  +coeff( 36)        *x31*x41*x51
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)        *x33    *x51
     6  +coeff( 42)    *x23*x32        
     7  +coeff( 43)        *x32        
     8  +coeff( 44)        *x31    *x51
      t_sp_den  =t_sp_den  
     9  +coeff( 45)*x12                
     1  +coeff( 46)    *x22        *x51
     2  +coeff( 47)    *x21*x31    *x51
     3  +coeff( 48)            *x41*x52
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)*x11        *x42    
c
      return
      end
      function y_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.6570819E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16595600E-02, 0.69077327E-02, 0.94333746E-01, 0.18806910E-01,
     + -0.58873501E-02, 0.26444090E-02, 0.16874226E-01, 0.17119162E-02,
     +  0.29989404E-02,-0.13161405E-01,-0.16757286E-02,-0.10871584E-01,
     +  0.41972282E-02, 0.72069321E-03,-0.48344494E-02,-0.18553897E-02,
     +  0.41415165E-02, 0.47828467E-02, 0.24731306E-03,-0.11067399E-02,
     + -0.10591220E-02, 0.10483731E-02,-0.84497058E-03,-0.12017983E-01,
     + -0.15143179E-01,-0.29185358E-02,-0.27698071E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x31    *x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)                *x52
      y_sp_den  =y_sp_den  
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x22            
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22*x33        
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)        *x31*x41    
      y_sp_den  =y_sp_den  
     9  +coeff( 18)            *x42    
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)            *x43    
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)        *x34        
     5  +coeff( 23)        *x33*x42    
     6  +coeff( 24)    *x22*x31*x41    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x22    *x41*x51
      y_sp_den  =y_sp_den  
     9  +coeff( 27)    *x22*x34        
c
      return
      end
      function p_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.7766830E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12794066E-02, 0.40573338E-02,-0.27289793E-01,-0.91905147E-01,
     + -0.87851789E-02, 0.85438620E-02,-0.42855702E-02,-0.23029933E-01,
     + -0.40574502E-02, 0.31461937E-02, 0.16199362E-01, 0.28816788E-02,
     + -0.98301750E-02,-0.43525440E-02, 0.10211241E-02, 0.11740278E-02,
     + -0.80453546E-03, 0.28270742E-03, 0.60361177E-02, 0.14493561E-03,
     + -0.11967461E-02, 0.13590476E-02, 0.85864076E-02,-0.55042664E-02,
     + -0.16144492E-02,-0.98099734E-03,-0.55181835E-03, 0.40000677E-02,
     + -0.93642075E-03, 0.12226967E-03, 0.91557042E-03, 0.27718124E-03,
     +  0.48772726E-03, 0.54005799E-02, 0.56552974E-03, 0.42204550E-03,
     + -0.30890580E-02, 0.59850798E-02,-0.69880957E-03,-0.42957524E-02,
     + -0.45456752E-03,-0.10547175E-02, 0.12301210E-03,-0.10528703E-02,
     +  0.18044144E-03,-0.21099372E-03,-0.16876523E-03,-0.26199265E-03,
     + -0.14075715E-03,-0.17132951E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sp_den  =p_sp_den  
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)            *x42    
      p_sp_den  =p_sp_den  
     9  +coeff( 18)*x11    *x31        
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)        *x31*x42    
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)    *x21*x33*x41    
     8  +coeff( 26)            *x42*x53
      p_sp_den  =p_sp_den  
     9  +coeff( 27)    *x22*x33*x41    
     1  +coeff( 28)    *x23*x33*x41    
     2  +coeff( 29)        *x31*x41    
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)    *x23            
     5  +coeff( 32)    *x22*x31        
     6  +coeff( 33)    *x21*x32        
     7  +coeff( 34)    *x21*x31*x41    
     8  +coeff( 35)        *x32*x41    
      p_sp_den  =p_sp_den  
     9  +coeff( 36)*x11*x21*x31        
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)    *x22*x31*x41    
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)    *x23    *x42    
     5  +coeff( 41)        *x32*x43    
     6  +coeff( 42)    *x23*x33        
     7  +coeff( 43)*x11                
     8  +coeff( 44)            *x43    
      p_sp_den  =p_sp_den  
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)        *x31    *x52
     2  +coeff( 47)                *x53
     3  +coeff( 48)*x11*x22            
     4  +coeff( 49)*x11*x21        *x51
     5  +coeff( 50)*x11        *x41*x51
c
      return
      end
      function l_sp_den  (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5965943E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73767728E-02, 0.21498601E+00,-0.95878253E-02,-0.77038156E-02,
     + -0.10467454E-01,-0.34655921E-01, 0.13820356E-01, 0.11870707E-01,
     + -0.47718203E-02,-0.17862283E-01, 0.91417413E-02,-0.38664946E-02,
     + -0.89136614E-02, 0.16523420E-02, 0.12338130E-01, 0.25655890E-01,
     +  0.54849114E-03,-0.16445735E-02, 0.81866549E-03, 0.52495701E-02,
     + -0.18483235E-01, 0.23897022E-02, 0.68097999E-02,-0.30288526E-02,
     + -0.11274141E-02, 0.89252386E-02,-0.20825362E-02,-0.36284563E-02,
     + -0.13450448E-01, 0.94260363E-03, 0.11067445E-01,-0.51334838E-03,
     +  0.16885042E-02,-0.18978967E-03, 0.46940604E-02, 0.12625728E-01,
     + -0.14912477E-02,-0.89440739E-03,-0.11727209E-02,-0.75704511E-03,
     +  0.46808855E-03,-0.11114210E-02, 0.36610534E-04, 0.90367009E-03,
     +  0.77356640E-02,-0.27563854E-03, 0.37083279E-02, 0.30426754E-03,
     + -0.21705800E-03, 0.11322309E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_den  =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23*x31        
      l_sp_den  =l_sp_den  
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x22*x33        
      l_sp_den  =l_sp_den  
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)    *x21*x33*x41    
     7  +coeff( 25)    *x21*x32    *x52
     8  +coeff( 26)    *x23*x33*x41    
      l_sp_den  =l_sp_den  
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)    *x21*x32        
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)    *x22    *x42    
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)*x11*x22*x31        
     7  +coeff( 34)*x11    *x33        
     8  +coeff( 35)*x11*x22    *x41    
      l_sp_den  =l_sp_den  
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)    *x22*x33*x41    
     2  +coeff( 38)*x11    *x31        
     3  +coeff( 39)    *x22        *x51
     4  +coeff( 40)    *x21    *x41*x51
     5  +coeff( 41)        *x31*x41*x51
     6  +coeff( 42)    *x21        *x52
     7  +coeff( 43)                *x53
     8  +coeff( 44)*x11*x21    *x41    
      l_sp_den  =l_sp_den  
     9  +coeff( 45)    *x22*x31*x41    
     1  +coeff( 46)        *x32    *x52
     2  +coeff( 47)    *x23*x32        
     3  +coeff( 48)                *x51
     4  +coeff( 49)        *x32        
     5  +coeff( 50)*x12                
c
      return
      end
      function x_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.8640270E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12057847E-02, 0.29819360E+00,-0.80028130E-02, 0.13309248E+00,
     + -0.14829338E-01, 0.21874694E-01,-0.56826949E-01, 0.41628830E-01,
     +  0.19791514E-01,-0.28857352E-01, 0.17183969E-01,-0.41758069E-02,
     + -0.35204873E-02, 0.43116439E-01, 0.23849830E-01,-0.34452047E-01,
     +  0.57188058E-02,-0.61099562E-02,-0.65695629E-02, 0.20626325E-01,
     + -0.52925432E-03, 0.11809538E-01,-0.56739845E-02, 0.11008143E-01,
     + -0.35576087E-02,-0.63975942E-02,-0.28247597E-01,-0.43303305E-02,
     +  0.85318433E-02,-0.67456713E-04, 0.23601539E-01,-0.16677384E-03,
     + -0.97688171E-03,-0.39122989E-02,-0.78251272E-04,-0.10199836E-02,
     +  0.14598056E-02, 0.41491116E-03,-0.14411261E-02,-0.11008361E-02,
     +  0.29407640E-02, 0.24342644E-02,-0.16716202E-02, 0.35606474E-02,
     + -0.15330155E-01, 0.55279741E-02,-0.30371665E-02,-0.24300742E-03,
     + -0.20555772E-02, 0.37063452E-03, 0.23340720E-02,-0.16227366E-02,
     +  0.15108566E-02,-0.30142270E-03, 0.17332432E-02, 0.25023243E-02,
     +  0.84179286E-02, 0.97542908E-02, 0.37247650E-03,-0.22296949E-03,
     + -0.24615466E-02, 0.76501799E-03,-0.26873311E-02, 0.85103144E-04,
     + -0.37600764E-03, 0.10344200E-02,-0.15751362E-02, 0.18764050E-02,
     +  0.11978856E-02, 0.12505713E-02, 0.84329542E-03, 0.18188575E-02,
     + -0.28980374E-02, 0.13326597E-02,-0.39351391E-03, 0.63964096E-02,
     + -0.83790597E-03,-0.80814952E-03, 0.10729609E-02, 0.59827801E-03,
     + -0.21387399E-02,-0.38398907E-03, 0.47026671E-03, 0.63335826E-03,
     +  0.61943894E-03,-0.19821757E-02, 0.28418144E-01, 0.25610121E-01,
     +  0.12987597E-02, 0.20275228E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)        *x31        
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      x_sp_dex    =x_sp_dex    
     9  +coeff( 18)    *x24            
     1  +coeff( 19)    *x24*x31        
     2  +coeff( 20)    *x23*x31*x41    
     3  +coeff( 21)    *x21*x33*x41    
     4  +coeff( 22)    *x23*x33*x41    
     5  +coeff( 23)            *x42    
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)    *x21*x32        
      x_sp_dex    =x_sp_dex    
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)        *x33*x41    
     4  +coeff( 31)    *x23    *x42    
     5  +coeff( 32)*x11    *x33        
     6  +coeff( 33)        *x32        
     7  +coeff( 34)        *x31*x41    
     8  +coeff( 35)        *x31    *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 36)            *x41*x51
     1  +coeff( 37)*x11*x21            
     2  +coeff( 38)*x11            *x51
     3  +coeff( 39)    *x21*x31    *x51
     4  +coeff( 40)    *x21        *x52
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x23        *x51
     7  +coeff( 43)*x11*x23            
     8  +coeff( 44)*x11*x22*x31        
      x_sp_dex    =x_sp_dex    
     9  +coeff( 45)    *x24    *x41    
     1  +coeff( 46)    *x23*x32        
     2  +coeff( 47)*x11*x24            
     3  +coeff( 48)*x11    *x31    *x52
     4  +coeff( 49)*x11*x23*x31        
     5  +coeff( 50)*x11*x22*x33        
     6  +coeff( 51)    *x24*x33*x41    
     7  +coeff( 52)*x11    *x31        
     8  +coeff( 53)    *x22        *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 54)            *x42*x51
     1  +coeff( 55)*x11*x21*x31        
     2  +coeff( 56)    *x22*x32        
     3  +coeff( 57)    *x22*x31*x41    
     4  +coeff( 58)    *x22    *x42    
     5  +coeff( 59)*x12*x21            
     6  +coeff( 60)*x11*x22        *x51
     7  +coeff( 61)*x11*x23    *x41    
     8  +coeff( 62)*x11*x21*x32    *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 63)    *x24    *x42*x52
     1  +coeff( 64)*x12                
     2  +coeff( 65)*x11        *x42    
     3  +coeff( 66)    *x21    *x42*x51
     4  +coeff( 67)    *x24        *x51
     5  +coeff( 68)*x11*x21*x31*x41    
     6  +coeff( 69)    *x23*x31    *x51
     7  +coeff( 70)    *x23    *x41*x51
     8  +coeff( 71)*x12*x21    *x41    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 72)    *x21*x31*x43    
     1  +coeff( 73)    *x21    *x43*x51
     2  +coeff( 74)    *x24    *x41*x51
     3  +coeff( 75)            *x43*x52
     4  +coeff( 76)    *x23    *x43    
     5  +coeff( 77)*x11    *x33*x41    
     6  +coeff( 78)*x11    *x32*x42    
     7  +coeff( 79)*x11*x24        *x51
     8  +coeff( 80)    *x22*x33    *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 81)*x11*x23    *x41*x51
     1  +coeff( 82)    *x21*x33    *x52
     2  +coeff( 83)*x12    *x32    *x51
     3  +coeff( 84)*x12    *x31*x41*x51
     4  +coeff( 85)*x12*x23*x31        
     5  +coeff( 86)*x11*x21*x33*x41    
     6  +coeff( 87)    *x23*x32*x42    
     7  +coeff( 88)    *x23*x31*x43    
     8  +coeff( 89)*x11*x21*x32*x41*x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 90)    *x23*x32*x41*x51
c
      return
      end
      function t_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.5357819E+00/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.38800418E-03,-0.69248095E-01, 0.15631110E-02, 0.22752339E-01,
     +  0.24440526E-02,-0.45277979E-02, 0.11711914E-01,-0.39880374E-02,
     + -0.40043122E-02, 0.59551452E-02,-0.12764016E-02,-0.38125301E-02,
     +  0.79605461E-03,-0.88258311E-02,-0.26316044E-03,-0.21800292E-02,
     +  0.59958855E-02,-0.76403789E-03,-0.13115729E-03,-0.45212107E-02,
     + -0.15676640E-03, 0.66016120E-03,-0.10376581E-02, 0.45636785E-02,
     + -0.17247087E-03, 0.66726370E-03, 0.46871896E-03, 0.52535161E-03,
     +  0.53168704E-04,-0.15357530E-02,-0.99620630E-03,-0.51623276E-02,
     + -0.15761238E-03, 0.10849949E-02,-0.24935379E-03, 0.16364678E-03,
     +  0.12574838E-03,-0.36497886E-03,-0.17616566E-02,-0.22486818E-02,
     + -0.59874519E-03, 0.88000852E-04, 0.44817178E-04, 0.28455327E-03,
     + -0.45480429E-04,-0.79029865E-04, 0.72130475E-04,-0.14956371E-03,
     + -0.11307062E-03,-0.29584154E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sp_dex    =t_sp_dex    
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23            
     4  +coeff( 13)        *x31        
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x42    
      t_sp_dex    =t_sp_dex    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x22*x33        
     2  +coeff( 20)    *x23*x31*x41    
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x22*x31        
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)        *x31*x41*x51
      t_sp_dex    =t_sp_dex    
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)    *x21        *x52
     2  +coeff( 29)*x11    *x33        
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)        *x32        
     7  +coeff( 34)    *x21*x32        
     8  +coeff( 35)    *x22        *x51
      t_sp_dex    =t_sp_dex    
     9  +coeff( 36)    *x21    *x41*x51
     1  +coeff( 37)        *x31    *x52
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)*x11*x22*x31        
     6  +coeff( 42)            *x42    
     7  +coeff( 43)        *x31    *x51
     8  +coeff( 44)*x11    *x31        
      t_sp_dex    =t_sp_dex    
     9  +coeff( 45)*x12                
     1  +coeff( 46)        *x32*x41    
     2  +coeff( 47)        *x32    *x51
     3  +coeff( 48)*x11*x21*x31        
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x22*x32        
c
      return
      end
      function y_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 46)
      data ncoeff/ 45/
      data avdat/  0.3980987E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.60068240E-03,-0.40456116E-01,-0.21679915E-01, 0.14515941E-01,
     + -0.12667259E-02, 0.72157261E-03, 0.40345932E-02, 0.55447798E-02,
     +  0.10021523E-01, 0.54854408E-01, 0.75811120E-02,-0.58896998E-02,
     + -0.40160015E-01,-0.84052207E-02,-0.64446381E-02, 0.13483377E-02,
     + -0.41273418E-02,-0.52642575E-02,-0.12266361E-02,-0.52940734E-02,
     + -0.28612126E-01,-0.94887037E-02, 0.75940676E-02,-0.34106028E-03,
     + -0.81328917E-02,-0.28223807E-02,-0.47831470E-02,-0.17049694E-02,
     +  0.40763738E-02,-0.97239629E-03,-0.85726375E-03,-0.26887399E-02,
     +  0.29268593E-03,-0.98351936E-03,-0.40968551E-03, 0.24096428E-03,
     +  0.19471204E-02, 0.16235754E-03,-0.79399487E-02,-0.96365521E-02,
     + -0.48513690E-03, 0.22712124E-02, 0.26495710E-02, 0.12453014E-02,
     + -0.17793487E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_dex    =y_sp_dex    
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)    *x22            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)            *x41*x52
      y_sp_dex    =y_sp_dex    
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x22    *x41    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x23            
     6  +coeff( 24)            *x45    
     7  +coeff( 25)    *x22    *x41*x51
     8  +coeff( 26)        *x31*x42    
      y_sp_dex    =y_sp_dex    
     9  +coeff( 27)            *x43    
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)    *x21    *x42    
     3  +coeff( 30)                *x53
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)    *x21*x33*x41    
     7  +coeff( 34)        *x31    *x54
     8  +coeff( 35)        *x32*x41    
      y_sp_dex    =y_sp_dex    
     9  +coeff( 36)*x11    *x31        
     1  +coeff( 37)    *x21*x31*x41    
     2  +coeff( 38)*x11            *x51
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)    *x23*x31        
     7  +coeff( 43)    *x23    *x41    
     8  +coeff( 44)    *x23        *x51
      y_sp_dex    =y_sp_dex    
     9  +coeff( 45)    *x22*x34        
c
      return
      end
      function p_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2147559E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13986450E-03,-0.11566449E-01,-0.21188619E-01, 0.60353783E-03,
     +  0.38521350E-03,-0.69356040E-03,-0.82351100E-02, 0.57101162E-03,
     + -0.20393103E-02, 0.17085631E-02, 0.11060227E-01, 0.17669873E-02,
     +  0.12402703E-02,-0.19432130E-03,-0.47368784E-02,-0.21394601E-02,
     + -0.14992320E-02, 0.36323743E-03, 0.28868581E-03, 0.34201110E-03,
     +  0.10482410E-02,-0.10291014E-02,-0.17283059E-03,-0.16013352E-02,
     +  0.29728489E-03,-0.15026292E-04, 0.54768729E-03,-0.53262163E-03,
     + -0.30223074E-03,-0.14414453E-03,-0.11997944E-03, 0.54841419E-03,
     + -0.45998412E-03,-0.41035557E-04, 0.45598215E-04, 0.38391780E-04,
     +  0.53238728E-04,-0.25402059E-03, 0.65260989E-04, 0.11128523E-03,
     +  0.52107527E-03,-0.18180315E-03,-0.23372870E-03,-0.98790822E-03,
     + -0.80679765E-03, 0.30763476E-04,-0.90067158E-04, 0.73969502E-04,
     +  0.17334615E-03, 0.20458191E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      p_sp_dex    =p_sp_dex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21    *x41*x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 18)    *x21            
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)    *x22    *x41*x51
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x23        *x52
      p_sp_dex    =p_sp_dex    
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)        *x31*x42    
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)                *x53
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)    *x22*x31    *x51
     7  +coeff( 34)*x11*x21            
     8  +coeff( 35)*x11    *x31        
      p_sp_dex    =p_sp_dex    
     9  +coeff( 36)*x11            *x51
     1  +coeff( 37)    *x21*x32        
     2  +coeff( 38)    *x21        *x52
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)*x11*x21    *x41    
     5  +coeff( 41)    *x23    *x41    
     6  +coeff( 42)    *x21    *x41*x52
     7  +coeff( 43)*x11*x22    *x41    
     8  +coeff( 44)    *x22    *x42*x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 45)    *x22*x33*x41*x51
     1  +coeff( 46)        *x32        
     2  +coeff( 47)        *x32*x41    
     3  +coeff( 48)        *x31*x41*x51
     4  +coeff( 49)    *x23*x31        
     5  +coeff( 50)    *x21    *x43    
c
      return
      end
      function l_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.5122880E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.36668594E-02,-0.42765871E+00, 0.50848066E-02,-0.98709986E-01,
     +  0.18899549E-01,-0.38787622E-01, 0.78461275E-01,-0.42277347E-01,
     + -0.27764939E-01, 0.40084288E-01,-0.24846904E-01,-0.60868479E-01,
     +  0.27592250E-02,-0.20664269E-02,-0.19641060E-02, 0.17410739E-02,
     + -0.13711909E-01, 0.44994127E-01,-0.51300339E-02,-0.59492368E-03,
     + -0.23417721E-01, 0.38609738E-02,-0.11154121E-01, 0.44849589E-02,
     + -0.12000247E-02,-0.65240981E-02, 0.83055384E-02, 0.33004571E-01,
     + -0.20069627E-02, 0.38812354E-02, 0.21615173E-02, 0.23953825E-03,
     + -0.10314750E-01,-0.32497905E-01,-0.53025532E-03, 0.10327203E-02,
     + -0.22848451E-02,-0.98924991E-02,-0.13592545E-01,-0.26785883E-02,
     + -0.42363084E-02, 0.42775754E-03,-0.75302618E-02, 0.11427928E-02,
     +  0.19426820E-02,-0.30382277E-03, 0.54590555E-03, 0.52029587E-03,
     + -0.90583129E-03,-0.75486250E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_sp_dex    =l_sp_dex    
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)        *x31        
     5  +coeff( 14)            *x42    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)                *x52
     8  +coeff( 17)    *x22    *x41    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x22*x33        
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)    *x21*x33*x41    
     5  +coeff( 23)    *x23*x33*x41    
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x22*x31        
      l_sp_dex    =l_sp_dex    
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)    *x22        *x51
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)    *x21        *x52
     5  +coeff( 32)*x11    *x33        
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)    *x23    *x42    
     8  +coeff( 35)        *x31    *x51
      l_sp_dex    =l_sp_dex    
     9  +coeff( 36)            *x42*x51
     1  +coeff( 37)*x11*x21    *x41    
     2  +coeff( 38)    *x22*x31*x41    
     3  +coeff( 39)    *x22    *x42    
     4  +coeff( 40)    *x23        *x51
     5  +coeff( 41)*x11*x22*x31        
     6  +coeff( 42)*x11    *x31    *x52
     7  +coeff( 43)    *x23*x32        
     8  +coeff( 44)        *x31*x41    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 45)*x11    *x31        
     1  +coeff( 46)*x12                
     2  +coeff( 47)    *x21*x31    *x51
     3  +coeff( 48)            *x41*x52
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)*x12*x21            
c
      return
      end
      function x_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5441129E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.18263733E-02, 0.18718468E+00,-0.52033416E-02, 0.13837682E+00,
     + -0.10355782E-01, 0.23248855E-01,-0.37469260E-01, 0.32334350E-01,
     +  0.13393950E-01,-0.19215807E-01,-0.59893942E-02, 0.12435450E-01,
     + -0.27093554E-02, 0.29110657E-01,-0.45670709E-02, 0.28997448E-02,
     + -0.23401381E-01, 0.36969520E-02,-0.70955786E-02, 0.14794544E-01,
     +  0.42735835E-03,-0.33503285E-03, 0.78655267E-02,-0.28907296E-02,
     + -0.10089526E-02, 0.54140673E-02,-0.22459344E-02, 0.12753326E-01,
     + -0.44422341E-02,-0.19330794E-01,-0.41364087E-02,-0.33556856E-02,
     +  0.25199063E-02, 0.53872205E-02, 0.14971131E-01, 0.28076602E-04,
     + -0.66236878E-03,-0.15737793E-02,-0.81671448E-03, 0.18285272E-02,
     + -0.91839151E-03, 0.24915924E-02, 0.44932100E-03, 0.39248862E-02,
     + -0.18327587E-02, 0.70092210E-03,-0.43098000E-03,-0.57900539E-04,
     +  0.32967571E-03,-0.11301547E-02, 0.37913810E-03, 0.74387679E-03,
     +  0.52788961E-02, 0.73972847E-02, 0.29374988E-03,-0.28637075E-02,
     +  0.75731601E-03, 0.26769959E-03,-0.39175446E-02,-0.23298268E-02,
     +  0.80346745E-02, 0.95332049E-04, 0.80978731E-04, 0.28043694E-03,
     +  0.11834399E-02, 0.66204101E-03, 0.13215985E-04, 0.10344903E-02,
     +  0.17445317E-02, 0.22646438E-02, 0.89006306E-03, 0.49070281E-03,
     + -0.12262848E-02, 0.18633036E-02,-0.19542240E-02, 0.56355278E-03,
     +  0.38957342E-02, 0.34227472E-03,-0.88613003E-03,-0.11128535E-02,
     +  0.24750145E-03, 0.13526175E-02, 0.20552967E-01, 0.16212558E-01,
     +  0.86854666E-03, 0.14813585E-02, 0.17773765E-02,-0.12626362E-02,
     +  0.46120756E-03,-0.17092389E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23            
     4  +coeff( 13)        *x31        
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21    *x42    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x24    *x41    
     2  +coeff( 20)    *x23*x31*x41    
     3  +coeff( 21)    *x22*x33        
     4  +coeff( 22)    *x21*x33*x41    
     5  +coeff( 23)    *x23*x33*x41    
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)    *x22*x31        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)    *x21*x31*x41    
     4  +coeff( 31)    *x21    *x41*x51
     5  +coeff( 32)    *x24            
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)*x11*x22    *x41    
     8  +coeff( 35)    *x23    *x42    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 36)*x11    *x33        
     1  +coeff( 37)        *x32        
     2  +coeff( 38)    *x21*x31    *x51
     3  +coeff( 39)    *x21        *x52
     4  +coeff( 40)*x11*x21    *x41    
     5  +coeff( 41)*x11*x23            
     6  +coeff( 42)*x11*x22*x31        
     7  +coeff( 43)*x11*x22        *x51
     8  +coeff( 44)    *x23*x32        
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 45)*x11*x24            
     1  +coeff( 46)*x11*x21*x33        
     2  +coeff( 47)*x11*x22*x33        
     3  +coeff( 48)        *x31    *x51
     4  +coeff( 49)*x11*x21            
     5  +coeff( 50)*x11    *x31        
     6  +coeff( 51)                *x53
     7  +coeff( 52)*x11*x21*x31        
     8  +coeff( 53)    *x22*x31*x41    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 54)    *x22    *x42    
     1  +coeff( 55)*x12*x21            
     2  +coeff( 56)    *x24*x31        
     3  +coeff( 57)    *x22*x32*x41    
     4  +coeff( 58)    *x24*x32        
     5  +coeff( 59)    *x24    *x42    
     6  +coeff( 60)*x11*x23*x33        
     7  +coeff( 61)*x12*x23    *x43*x52
     8  +coeff( 62)*x11            *x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 63)*x12                
     1  +coeff( 64)        *x31*x41*x51
     2  +coeff( 65)    *x22*x32        
     3  +coeff( 66)    *x21    *x42*x51
     4  +coeff( 67)        *x31*x41*x52
     5  +coeff( 68)*x11*x21*x31*x41    
     6  +coeff( 69)    *x23*x31    *x51
     7  +coeff( 70)    *x23    *x41*x51
     8  +coeff( 71)    *x23        *x52
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 72)    *x22*x32    *x51
     1  +coeff( 73)*x11*x23    *x41    
     2  +coeff( 74)    *x21*x31*x43    
     3  +coeff( 75)    *x21    *x43*x51
     4  +coeff( 76)    *x24        *x52
     5  +coeff( 77)    *x23    *x43    
     6  +coeff( 78)*x11*x21*x32    *x51
     7  +coeff( 79)    *x22    *x42*x52
     8  +coeff( 80)*x11*x23    *x41*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 81)*x12*x23*x31        
     1  +coeff( 82)*x11*x21*x32*x42    
     2  +coeff( 83)    *x23*x32*x42    
     3  +coeff( 84)    *x23*x31*x43    
     4  +coeff( 85)*x11*x21*x32*x41*x51
     5  +coeff( 86)    *x23*x32*x41*x51
     6  +coeff( 87)*x11*x24*x31*x41    
     7  +coeff( 88)    *x22*x32*x41*x52
     8  +coeff( 89)*x12    *x32*x41*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 90)*x11*x22*x31*x41*x52
c
      return
      end
      function t_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.7466303E-03/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22261303E-03,-0.69679700E-01, 0.15967498E-02, 0.22706876E-01,
     +  0.24695715E-02,-0.58396826E-02, 0.11721873E-01,-0.49429825E-02,
     + -0.42956783E-02, 0.60272873E-02,-0.14645333E-02,-0.37433193E-02,
     +  0.77386911E-03,-0.91379378E-02, 0.64591258E-02,-0.79311727E-03,
     + -0.38382402E-02, 0.56814845E-03,-0.15970926E-02, 0.68731967E-03,
     + -0.21269887E-03, 0.12179749E-02,-0.23527993E-02, 0.49750376E-02,
     + -0.46562785E-03, 0.42951771E-03, 0.61182422E-04,-0.15556683E-02,
     + -0.43495354E-04,-0.51194024E-02, 0.41100068E-03, 0.13895734E-03,
     + -0.10483894E-02, 0.28988410E-03,-0.39246021E-03, 0.39244225E-03,
     + -0.48758846E-03, 0.22034009E-03,-0.62131876E-03,-0.11659012E-02,
     +  0.39131584E-03, 0.13459165E-03, 0.14910683E-03,-0.79011086E-04,
     + -0.57699701E-04, 0.29378527E-03,-0.45658315E-04,-0.14502440E-03,
     +  0.11328267E-03,-0.15945997E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23            
     4  +coeff( 13)        *x31        
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x23*x31*x41    
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 18)    *x21*x33*x41    
     1  +coeff( 19)    *x23*x33*x41    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x21*x32        
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)    *x21        *x52
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 27)*x11    *x33        
     1  +coeff( 28)*x11*x22    *x41    
     2  +coeff( 29)    *x22*x33        
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)        *x31*x41    
     5  +coeff( 32)*x11*x21            
     6  +coeff( 33)    *x22*x31        
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)*x11*x21    *x41    
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 36)        *x32*x42    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)            *x42*x52
     3  +coeff( 39)*x11*x22*x31        
     4  +coeff( 40)    *x23*x32        
     5  +coeff( 41)        *x33*x43    
     6  +coeff( 42)        *x32        
     7  +coeff( 43)            *x42    
     8  +coeff( 44)        *x31    *x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 45)            *x41*x51
     1  +coeff( 46)*x11    *x31        
     2  +coeff( 47)*x12                
     3  +coeff( 48)        *x31*x41*x51
     4  +coeff( 49)                *x53
     5  +coeff( 50)*x11*x21*x31        
c
      return
      end
      function y_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/  0.5992274E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.60555164E-03,-0.52393198E-01,-0.41159671E-01, 0.15595970E-01,
     + -0.14471516E-02, 0.82288031E-03, 0.44011907E-02, 0.60403524E-02,
     +  0.12787987E-01, 0.68366982E-01, 0.93741994E-02,-0.44806423E-02,
     + -0.45636129E-01,-0.10048630E-01,-0.67933202E-02, 0.13959680E-02,
     + -0.53402483E-02,-0.71916054E-02,-0.13193161E-02,-0.56421943E-02,
     + -0.33116285E-01,-0.12208416E-01, 0.97408658E-02,-0.10124039E-01,
     + -0.62632151E-02,-0.20266371E-02, 0.54989513E-02,-0.12898825E-02,
     + -0.80108660E-03,-0.40434678E-02, 0.26167254E-03,-0.12864423E-02,
     +  0.18836326E-03,-0.26902193E-02, 0.28432508E-02,-0.10875728E-02,
     +  0.22359550E-03,-0.47880152E-03,-0.82278019E-02,-0.88139074E-02,
     + -0.77944098E-03,-0.19991179E-02,-0.23101200E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_q3en   =y_sp_q3en   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)    *x22            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)            *x41*x52
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x22    *x41    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x23            
     6  +coeff( 24)    *x22    *x41*x51
     7  +coeff( 25)            *x43    
     8  +coeff( 26)            *x42*x51
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)                *x53
     2  +coeff( 29)        *x33*x42    
     3  +coeff( 30)    *x22*x31    *x51
     4  +coeff( 31)    *x21*x33*x41    
     5  +coeff( 32)        *x31    *x54
     6  +coeff( 33)        *x33        
     7  +coeff( 34)        *x31*x42    
     8  +coeff( 35)    *x21*x31*x41    
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)    *x21        *x52
     2  +coeff( 38)        *x34*x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)    *x22*x34        
     7  +coeff( 43)    *x23        *x52
c
      return
      end
      function p_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2082444E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12009752E-03,-0.12030553E-01,-0.19446593E-01, 0.12406372E-02,
     +  0.10447614E-03,-0.14487929E-02,-0.99438634E-02, 0.73233491E-03,
     + -0.20716570E-02, 0.18877115E-02, 0.12397034E-01, 0.20057440E-02,
     +  0.14274076E-02,-0.59861731E-03,-0.61766924E-02,-0.24598034E-02,
     +  0.35630952E-03, 0.36420178E-03, 0.13587947E-02,-0.13874469E-02,
     + -0.21627052E-02, 0.40859741E-03, 0.90726302E-03,-0.11979523E-02,
     + -0.82878780E-03,-0.98365053E-04, 0.78975972E-04,-0.67053968E-03,
     +  0.12659325E-03,-0.18293780E-03,-0.17521573E-03,-0.27312647E-03,
     + -0.15307886E-03,-0.12460812E-03, 0.36357786E-03, 0.51366624E-04,
     +  0.10525883E-03,-0.13960565E-03,-0.74990647E-04,-0.31360054E-04,
     +  0.10608682E-03,-0.12328963E-03, 0.51445637E-04, 0.45757426E-03,
     + -0.96485292E-03,-0.68124844E-03, 0.38952665E-04, 0.52452491E-04,
     + -0.36962061E-04,-0.40602270E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21            
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)        *x31*x41    
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x22*x31    *x51
     8  +coeff( 26)*x11*x21            
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 27)*x11    *x31        
     1  +coeff( 28)        *x31*x42    
     2  +coeff( 29)        *x31*x41*x51
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)    *x21        *x52
     5  +coeff( 32)            *x41*x52
     6  +coeff( 33)                *x53
     7  +coeff( 34)*x11*x21        *x51
     8  +coeff( 35)    *x23        *x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 36)        *x32        
     1  +coeff( 37)    *x21*x32        
     2  +coeff( 38)        *x32*x41    
     3  +coeff( 39)    *x21*x31    *x51
     4  +coeff( 40)*x11*x21*x31        
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x22        *x52
     7  +coeff( 43)*x11            *x53
     8  +coeff( 44)    *x23    *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 45)    *x22    *x42*x51
     1  +coeff( 46)    *x22*x33*x41*x51
     2  +coeff( 47)        *x32    *x51
     3  +coeff( 48)*x11*x22            
     4  +coeff( 49)*x11    *x32        
     5  +coeff( 50)*x11        *x41*x51
c
      return
      end
      function l_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3503092E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28317224E-02,-0.27826852E+00,-0.32169659E-01, 0.11541987E-01,
     + -0.33092663E-01, 0.50581560E-01,-0.19590307E-01,-0.18187188E-01,
     +  0.25751619E-01,-0.16607232E-01,-0.62726103E-02,-0.39572030E-01,
     + -0.16778590E-02,-0.22291050E-02,-0.10161637E-02, 0.27733974E-01,
     + -0.33222437E-02,-0.15763484E-01, 0.19328902E-02,-0.56712450E-02,
     + -0.48067240E-03, 0.27725284E-02, 0.51600737E-02, 0.20716716E-01,
     + -0.22290214E-02, 0.14052034E-02, 0.11642862E-02, 0.10338690E-02,
     + -0.61423727E-02,-0.21388087E-01,-0.33586381E-02, 0.16136508E-02,
     + -0.79495439E-04,-0.35943717E-03, 0.53957780E-03,-0.39210063E-02,
     + -0.49936298E-04,-0.77138818E-02, 0.18869230E-02,-0.14964947E-02,
     + -0.19948492E-02,-0.45940615E-02,-0.12440157E-03, 0.97012974E-03,
     +  0.57359127E-03,-0.16664922E-03,-0.55619672E-04, 0.90283266E-03,
     +  0.46633123E-03,-0.61945728E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23*x31        
      l_sp_q3en   =l_sp_q3en   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)    *x21*x33*x41    
     2  +coeff( 20)    *x23*x33*x41    
     3  +coeff( 21)        *x31    *x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x32        
     6  +coeff( 24)    *x21*x31*x41    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)    *x21    *x41*x51
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 27)    *x21        *x52
     1  +coeff( 28)*x11    *x33        
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)*x11*x22*x33        
     5  +coeff( 32)            *x41    
     6  +coeff( 33)        *x32        
     7  +coeff( 34)                *x52
     8  +coeff( 35)*x11*x21            
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 36)    *x22*x31        
     1  +coeff( 37)        *x33        
     2  +coeff( 38)    *x22    *x41    
     3  +coeff( 39)            *x42*x51
     4  +coeff( 40)*x11*x21    *x41    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x23*x32        
     7  +coeff( 43)    *x22*x33        
     8  +coeff( 44)        *x31        
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 45)*x11    *x31        
     1  +coeff( 46)*x12                
     2  +coeff( 47)        *x32    *x51
     3  +coeff( 48)        *x31*x41*x51
     4  +coeff( 49)            *x41*x52
     5  +coeff( 50)*x11*x21*x31        
c
      return
      end
      function x_sp_q3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5276050E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26690743E-02, 0.99666663E-01,-0.17005422E-02,-0.31252822E-02,
     +  0.19226252E+00,-0.76683396E-02, 0.16653227E-01,-0.23468202E-01,
     +  0.27976913E-01,-0.99408692E-02,-0.11805682E-01, 0.80005582E-02,
     + -0.41312408E-02,-0.16294403E-01,-0.44048685E-02, 0.18631617E-01,
     +  0.11913989E-02, 0.52892067E-03,-0.12876577E-02, 0.87444847E-02,
     +  0.30980676E-02,-0.28697962E-02,-0.13528015E-01,-0.16113600E-02,
     +  0.23275116E-02,-0.32337387E-02, 0.80012307E-02, 0.15014684E-02,
     + -0.17888485E-02,-0.66723231E-04,-0.20360209E-02, 0.36572390E-02,
     + -0.14292356E-02,-0.37849706E-03, 0.67475141E-03, 0.37263341E-02,
     + -0.25788162E-04, 0.62230349E-04,-0.35031553E-03,-0.27024813E-03,
     + -0.68867794E-03, 0.11236360E-02, 0.20417127E-03, 0.14022002E-02,
     + -0.48110448E-02, 0.11855580E-01, 0.90667652E-02,-0.12611533E-02,
     + -0.29384240E-03, 0.64649631E-03, 0.11434071E-03,-0.15802676E-03,
     +  0.59645005E-04, 0.70392241E-03, 0.35282355E-02, 0.42837180E-03,
     + -0.12844678E-02, 0.73366813E-04, 0.22520409E-02, 0.16389057E-02,
     +  0.29707586E-02,-0.71858848E-03, 0.11592550E-02,-0.96611306E-03,
     +  0.55735378E-03, 0.12897847E-02, 0.15228703E-02,-0.74301541E-04,
     +  0.43097781E-02, 0.29955080E-03, 0.74737583E-03, 0.78503806E-04,
     +  0.10855510E-01, 0.79545416E-02,-0.84228301E-03, 0.36369823E-03,
     +  0.16457829E-03, 0.36754415E-02, 0.19034880E-02,-0.79466566E-03,
     +  0.27480046E-03,-0.39534055E-03,-0.38680565E-03, 0.70986530E-03,
     +  0.21342556E-02, 0.35699242E-03,-0.23458352E-03,-0.43577020E-03,
     + -0.32781149E-03, 0.32607536E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      x_sp_q3m   =x_sp_q3m   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x23            
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x33*x41    
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 18)    *x23*x33        
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)    *x21*x32        
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)    *x24            
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)    *x24*x31        
     3  +coeff( 30)        *x33*x41    
     4  +coeff( 31)        *x31*x41    
     5  +coeff( 32)    *x22*x31        
     6  +coeff( 33)*x11        *x41    
     7  +coeff( 34)    *x21        *x52
     8  +coeff( 35)                *x53
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)*x11    *x33        
     2  +coeff( 38)*x11*x24*x31        
     3  +coeff( 39)        *x32        
     4  +coeff( 40)        *x31    *x51
     5  +coeff( 41)*x11    *x31        
     6  +coeff( 42)*x11*x21    *x41    
     7  +coeff( 43)*x12*x21            
     8  +coeff( 44)*x11*x22*x31        
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 45)    *x24    *x41    
     1  +coeff( 46)    *x23*x31*x41    
     2  +coeff( 47)    *x23    *x42    
     3  +coeff( 48)*x11*x24            
     4  +coeff( 49)*x11*x23*x31        
     5  +coeff( 50)    *x22*x31*x41*x52
     6  +coeff( 51)    *x22    *x42*x52
     7  +coeff( 52)*x11            *x51
     8  +coeff( 53)*x12                
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 54)*x11*x21*x31        
     1  +coeff( 55)    *x22    *x42    
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)    *x24        *x51
     4  +coeff( 58)        *x31*x41*x52
     5  +coeff( 59)    *x23*x32        
     6  +coeff( 60)    *x23*x31    *x51
     7  +coeff( 61)    *x23    *x41*x51
     8  +coeff( 62)*x11*x23    *x41    
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 63)    *x21*x31*x43    
     1  +coeff( 64)    *x21    *x43*x51
     2  +coeff( 65)    *x24*x32        
     3  +coeff( 66)    *x24*x31*x41    
     4  +coeff( 67)    *x24    *x41*x51
     5  +coeff( 68)            *x42*x53
     6  +coeff( 69)    *x23    *x43    
     7  +coeff( 70)*x11*x21*x32    *x51
     8  +coeff( 71)    *x22*x33*x41    
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 72)    *x22*x33    *x51
     1  +coeff( 73)    *x23*x32*x42    
     2  +coeff( 74)    *x23*x31*x43    
     3  +coeff( 75)*x11*x23*x33        
     4  +coeff( 76)*x12    *x32*x41*x51
     5  +coeff( 77)        *x33*x41*x53
     6  +coeff( 78)    *x23*x33*x42    
     7  +coeff( 79)*x11*x24*x31*x42    
     8  +coeff( 80)*x12*x22*x32    *x52
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 81)*x11*x21            
     1  +coeff( 82)        *x31*x41*x51
     2  +coeff( 83)            *x42*x51
     3  +coeff( 84)    *x22*x32        
     4  +coeff( 85)    *x22*x31*x41    
     5  +coeff( 86)    *x22*x31    *x51
     6  +coeff( 87)    *x22        *x52
     7  +coeff( 88)*x11*x23            
     8  +coeff( 89)    *x21*x33        
      x_sp_q3m   =x_sp_q3m   
     9  +coeff( 90)    *x21*x32*x41    
c
      return
      end
      function t_sp_q3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.1313364E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15330833E-02,-0.37301030E-01, 0.73331239E-03, 0.67290090E-01,
     + -0.57545104E-02,-0.15913027E-02, 0.47594151E-02, 0.96563733E-03,
     + -0.17937506E-02, 0.30981074E-03, 0.33045185E-03, 0.24046688E-02,
     + -0.59370196E-03,-0.15089466E-02,-0.36999513E-02,-0.29248910E-03,
     +  0.55356801E-03,-0.34323963E-03, 0.10022512E-03, 0.14783227E-04,
     +  0.26202970E-03,-0.10261423E-03, 0.10184891E-02, 0.92221989E-03,
     + -0.38871789E-03, 0.62543542E-04,-0.99412129E-04, 0.11243291E-03,
     + -0.34114614E-03, 0.16770299E-03,-0.65378321E-03,-0.18763520E-03,
     + -0.33115898E-03,-0.29484334E-03, 0.23577601E-03,-0.27466929E-03,
     +  0.22699813E-03,-0.54417481E-03,-0.31718516E-03,-0.20080607E-03,
     + -0.15921501E-03,-0.32440352E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_sp_q3m   =t_sp_q3m   
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)        *x31        
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)                *x53
      t_sp_q3m   =t_sp_q3m   
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)        *x32        
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)*x11        *x41    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)        *x31*x41    
      t_sp_q3m   =t_sp_q3m   
     9  +coeff( 27)        *x31    *x51
     1  +coeff( 28)*x11    *x31        
     2  +coeff( 29)    *x22*x31        
     3  +coeff( 30)    *x21*x32        
     4  +coeff( 31)    *x22    *x41    
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)        *x31*x41*x51
     7  +coeff( 34)    *x22*x31*x41    
     8  +coeff( 35)    *x21*x31*x41*x51
      t_sp_q3m   =t_sp_q3m   
     9  +coeff( 36)    *x22        *x52
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)        *x32*x42*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x21*x31    *x53
     6  +coeff( 42)*x11*x22*x33        
c
      return
      end
      function y_sp_q3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/  0.8192663E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30633478E-03,-0.62903054E-01,-0.63371375E-01, 0.15538334E-01,
     +  0.10382978E-02, 0.48455964E-02, 0.66091148E-02, 0.13484588E-01,
     +  0.76953731E-01, 0.11045639E-01,-0.58974731E-02,-0.53817075E-01,
     + -0.11396248E-01,-0.58865412E-02, 0.17640688E-02,-0.48952006E-04,
     +  0.59227399E-02,-0.84244972E-02,-0.12525142E-02,-0.12566219E-02,
     + -0.54675606E-02,-0.37091400E-01,-0.14293799E-01, 0.93676588E-02,
     + -0.22344126E-02,-0.12593025E-01,-0.22687064E-02, 0.34010231E-02,
     + -0.13053796E-02,-0.17018528E-02,-0.48262225E-02,-0.59563448E-02,
     + -0.42072185E-02,-0.11319834E-02,-0.26278335E-02,-0.88327443E-02,
     + -0.92203189E-02,-0.89926866E-03,-0.10079789E-02,-0.26149792E-02,
     +  0.14116827E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sp_q3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x32        
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sp_q3m   =y_sp_q3m   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)    *x22            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x42    
      y_sp_q3m   =y_sp_q3m   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)    *x22    *x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x23            
     7  +coeff( 25)            *x43*x52
     8  +coeff( 26)    *x22    *x41*x51
      y_sp_q3m   =y_sp_q3m   
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)    *x21*x31*x41    
     2  +coeff( 29)                *x53
     3  +coeff( 30)        *x31*x42*x52
     4  +coeff( 31)    *x22*x31    *x51
     5  +coeff( 32)            *x43    
     6  +coeff( 33)            *x41*x52
     7  +coeff( 34)    *x21*x31    *x51
     8  +coeff( 35)        *x31*x44    
      y_sp_q3m   =y_sp_q3m   
     9  +coeff( 36)    *x22*x31*x41    
     1  +coeff( 37)    *x22    *x42    
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)        *x31    *x54
     4  +coeff( 40)    *x22*x34        
     5  +coeff( 41)*x11*x22*x33        
c
      return
      end
      function p_sp_q3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.9917168E-04/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.32212277E-03, 0.55244658E-03, 0.43015797E-02,-0.40101134E-02,
     + -0.28970882E-02, 0.17787850E-02, 0.28183140E-03, 0.45639807E-02,
     + -0.94604882E-03, 0.11475001E-02,-0.26377139E-02,-0.88733053E-02,
     + -0.75078366E-03,-0.12095548E-02, 0.41801841E-02, 0.10407848E-02,
     +  0.21847303E-02, 0.56405264E-04,-0.68124803E-03, 0.26686021E-03,
     +  0.12459691E-02, 0.84371597E-03, 0.49930147E-03, 0.33396817E-03,
     + -0.85743464E-04, 0.79908832E-05,-0.12517210E-03,-0.11924611E-03,
     +  0.58816595E-03, 0.24448370E-03, 0.39865225E-03, 0.38732300E-03,
     +  0.18303966E-02,-0.47372014E-03,-0.44327072E-03,-0.28467071E-04,
     +  0.16775086E-03, 0.11788191E-03,-0.10226573E-03, 0.39395843E-04,
     +  0.68204594E-04,-0.43836949E-03,-0.60848560E-03, 0.15568013E-02,
     +  0.23941603E-03, 0.19528788E-03,-0.25201286E-03,-0.27890952E-03,
     + -0.21122373E-03,-0.81516610E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sp_q3m   =p_sp_q3m   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)            *x41*x52
      p_sp_q3m   =p_sp_q3m   
     9  +coeff( 18)    *x22*x33        
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)            *x43    
     5  +coeff( 23)        *x31    *x52
     6  +coeff( 24)                *x53
     7  +coeff( 25)        *x32        
     8  +coeff( 26)*x11    *x31        
      p_sp_q3m   =p_sp_q3m   
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)    *x21    *x42    
     2  +coeff( 29)        *x31*x42    
     3  +coeff( 30)    *x21*x31    *x51
     4  +coeff( 31)    *x21    *x41*x51
     5  +coeff( 32)            *x42*x51
     6  +coeff( 33)    *x22    *x42    
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x22*x33*x41    
      p_sp_q3m   =p_sp_q3m   
     9  +coeff( 36)*x11            *x51
     1  +coeff( 37)        *x32*x41    
     2  +coeff( 38)        *x31*x41*x51
     3  +coeff( 39)    *x21        *x52
     4  +coeff( 40)*x11    *x32        
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)    *x23*x31        
     7  +coeff( 43)    *x23    *x41    
     8  +coeff( 44)    *x22*x31*x41    
      p_sp_q3m   =p_sp_q3m   
     9  +coeff( 45)    *x22*x31    *x51
     1  +coeff( 46)    *x22    *x41*x51
     2  +coeff( 47)    *x22        *x52
     3  +coeff( 48)    *x21    *x41*x52
     4  +coeff( 49)            *x41*x53
     5  +coeff( 50)    *x23    *x42    
c
      return
      end
      function l_sp_q3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4672457E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40067271E-02,-0.27834362E+00,-0.32095458E-01, 0.11551093E-01,
     + -0.35978667E-01, 0.50659467E-01,-0.16766712E-01,-0.17490057E-01,
     +  0.25602091E-01,-0.16765334E-01,-0.56796684E-02,-0.39290633E-01,
     + -0.10969775E-02,-0.24651620E-02,-0.15180449E-02,-0.10878083E-02,
     +  0.27892783E-01,-0.33322573E-02,-0.16556611E-01, 0.18930712E-02,
     + -0.54639373E-02,-0.61931857E-03, 0.59940823E-03, 0.27020634E-02,
     +  0.54090433E-02, 0.20898072E-01,-0.20424866E-02, 0.11050794E-02,
     +  0.11733299E-02, 0.18348638E-02,-0.60947663E-02,-0.22153869E-01,
     + -0.33198441E-02, 0.16016199E-02,-0.23055765E-03,-0.54577487E-02,
     +  0.14015023E-02,-0.96712419E-03,-0.57576792E-02,-0.63231657E-02,
     + -0.15603992E-02,-0.53067440E-02,-0.27185585E-02, 0.32618616E-02,
     +  0.13672257E-02, 0.18478910E-02, 0.23969610E-02,-0.26669665E-03,
     + -0.34492931E-02,-0.44165673E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_q3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23*x31        
      l_sp_q3m   =l_sp_q3m   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21    *x42    
      l_sp_q3m   =l_sp_q3m   
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x23*x31*x41    
     2  +coeff( 20)    *x21*x33*x41    
     3  +coeff( 21)    *x23*x33*x41    
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)*x11*x21            
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x21*x32        
     8  +coeff( 26)    *x21*x31*x41    
      l_sp_q3m   =l_sp_q3m   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)    *x21        *x52
     3  +coeff( 30)*x11    *x33        
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)*x11*x22*x33        
     7  +coeff( 34)            *x41    
     8  +coeff( 35)        *x32        
      l_sp_q3m   =l_sp_q3m   
     9  +coeff( 36)    *x22    *x41    
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x22*x31*x41    
     4  +coeff( 40)    *x22    *x42    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x23*x32        
     7  +coeff( 43)    *x22*x33        
     8  +coeff( 44)        *x33*x42    
      l_sp_q3m   =l_sp_q3m   
     9  +coeff( 45)        *x32*x43    
     1  +coeff( 46)    *x22*x31*x41*x51
     2  +coeff( 47)        *x33    *x52
     3  +coeff( 48)        *x32*x41*x52
     4  +coeff( 49)    *x22*x33    *x52
     5  +coeff( 50)        *x33*x42*x52
c
      return
      end
      function x_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1436917E-01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.77625280E-02, 0.70348993E-01,-0.29707579E-02, 0.32324407E+00,
     + -0.92662312E-02,-0.22930797E-01, 0.33500828E-01,-0.20956283E-01,
     + -0.49533471E-02, 0.19566106E-01,-0.11639781E-01, 0.80488026E-02,
     + -0.16966362E-02,-0.65931045E-02,-0.18157788E-01,-0.19770160E-02,
     +  0.92657702E-02, 0.35684137E-02,-0.13024800E-01,-0.61278101E-02,
     + -0.93623169E-03, 0.19368052E-02, 0.24982223E-02, 0.19286532E-01,
     + -0.15279453E-02, 0.58728434E-04,-0.61003509E-04,-0.28345590E-02,
     +  0.46253923E-03, 0.37298473E-02,-0.14149972E-02,-0.31330977E-02,
     + -0.21525659E-02, 0.80100261E-02, 0.13337355E-02, 0.36652200E-02,
     + -0.78630364E-05,-0.66400040E-03,-0.26014377E-03, 0.12216776E-02,
     +  0.15020695E-02,-0.44720341E-02, 0.10168410E-01, 0.95234551E-02,
     + -0.29173761E-03,-0.43156842E-03,-0.85960247E-03,-0.54787809E-03,
     +  0.42987906E-03, 0.61023519E-02,-0.16418936E-02, 0.17919428E-03,
     + -0.21536453E-03, 0.56960701E-03,-0.19236455E-02, 0.27545695E-02,
     +  0.23185555E-02, 0.43019517E-02,-0.88741782E-03,-0.10834764E-02,
     +  0.21841261E-02,-0.13289881E-03,-0.18082196E-02,-0.40144581E-03,
     + -0.11839271E-03,-0.20530495E-03, 0.13770237E-03, 0.15261942E-03,
     +  0.99528220E-03, 0.40993858E-02,-0.30423279E-03,-0.11359083E-02,
     +  0.29756516E-03, 0.10125490E-03, 0.50676009E-03, 0.56239363E-03,
     +  0.28952741E-03,-0.16997634E-03,-0.11024334E-02, 0.92523341E-03,
     + -0.60323207E-03,-0.14821842E-03,-0.31477390E-02, 0.12297248E-02,
     +  0.34048580E-03, 0.15985488E-02, 0.51132531E-03, 0.42626099E-02,
     +  0.30243492E-02,-0.48266616E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff(  9)    *x24            
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x23            
     4  +coeff( 13)        *x31        
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)    *x22    *x41    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)                *x53
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x24*x31        
     8  +coeff( 26)        *x33*x41    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 27)    *x23*x33        
     1  +coeff( 28)        *x31*x41    
     2  +coeff( 29)*x11*x21            
     3  +coeff( 30)    *x22*x31        
     4  +coeff( 31)*x11        *x41    
     5  +coeff( 32)    *x21*x32        
     6  +coeff( 33)    *x21*x31    *x51
     7  +coeff( 34)    *x23*x31        
     8  +coeff( 35)    *x23        *x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)        *x31    *x53
     2  +coeff( 38)*x11    *x31        
     3  +coeff( 39)*x11            *x51
     4  +coeff( 40)*x11*x21    *x41    
     5  +coeff( 41)*x11*x22*x31        
     6  +coeff( 42)    *x24    *x41    
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)    *x23    *x42    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 45)        *x32        
     1  +coeff( 46)        *x31    *x51
     2  +coeff( 47)        *x31*x41*x51
     3  +coeff( 48)            *x42*x51
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)    *x22    *x42    
     6  +coeff( 51)    *x22        *x52
     7  +coeff( 52)*x12*x21            
     8  +coeff( 53)*x11*x23            
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 54)*x11*x22        *x51
     1  +coeff( 55)    *x24        *x51
     2  +coeff( 56)    *x23*x32        
     3  +coeff( 57)    *x23*x31    *x51
     4  +coeff( 58)    *x23    *x41*x51
     5  +coeff( 59)*x11*x24            
     6  +coeff( 60)    *x24*x31*x41    
     7  +coeff( 61)    *x24    *x41*x51
     8  +coeff( 62)        *x32    *x53
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 63)*x11*x23*x33        
     1  +coeff( 64)*x12*x21    *x42*x51
     2  +coeff( 65)    *x24*x32    *x52
     3  +coeff( 66)        *x32    *x51
     4  +coeff( 67)            *x41*x52
     5  +coeff( 68)*x11*x21        *x51
     6  +coeff( 69)    *x22*x32        
     7  +coeff( 70)    *x22*x31*x41    
     8  +coeff( 71)    *x22*x31    *x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 72)    *x21*x31*x42    
     1  +coeff( 73)    *x21        *x53
     2  +coeff( 74)        *x32*x41*x51
     3  +coeff( 75)        *x31*x41*x52
     4  +coeff( 76)            *x42*x52
     5  +coeff( 77)*x11*x21*x31*x41    
     6  +coeff( 78)*x11    *x31    *x52
     7  +coeff( 79)*x11*x23    *x41    
     8  +coeff( 80)    *x21*x31*x43    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 81)    *x21    *x43*x51
     1  +coeff( 82)*x11*x22*x32        
     2  +coeff( 83)    *x24    *x42    
     3  +coeff( 84)    *x24*x31    *x51
     4  +coeff( 85)*x11*x22        *x52
     5  +coeff( 86)    *x24        *x52
     6  +coeff( 87)*x11*x21*x33        
     7  +coeff( 88)    *x23*x31*x42    
     8  +coeff( 89)    *x23    *x43    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 90)    *x23*x32    *x51
c
      return
      end
      function t_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 36)
      data ncoeff/ 35/
      data avdat/  0.5915408E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27999158E-02,-0.14091840E-01, 0.61565668E-04, 0.11757703E+00,
     +  0.69183009E-02,-0.10434740E-01,-0.17062097E-02, 0.21725507E-02,
     + -0.18139103E-02, 0.69243379E-03,-0.10479802E-02,-0.25084049E-02,
     + -0.12910351E-02,-0.88192354E-03, 0.33257683E-03, 0.26215034E-03,
     + -0.60672796E-03, 0.80732760E-04,-0.52122772E-03, 0.57321158E-04,
     + -0.11600377E-02, 0.55807777E-03, 0.14978177E-02,-0.98045298E-03,
     + -0.42435317E-03,-0.12582314E-03,-0.29953758E-03, 0.44475595E-03,
     + -0.24911307E-03, 0.11017892E-02, 0.57057920E-03, 0.31880112E-03,
     +  0.39793612E-03,-0.33829070E-03, 0.30106425E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)        *x32*x41*x51
     8  +coeff( 17)    *x21*x31        
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 18)        *x32        
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)                *x53
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x22*x31*x41*x51
     7  +coeff( 25)        *x31*x41    
     8  +coeff( 26)        *x31    *x51
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)        *x31*x42*x51
     7  +coeff( 34)    *x22        *x52
     8  +coeff( 35)*x11*x23            
c
      return
      end
      function y_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/  0.6077809E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51365979E-03,-0.40026993E-01,-0.49192239E-01, 0.78222537E-02,
     +  0.32090698E-03, 0.26648878E-02, 0.52082879E-02, 0.40089905E-01,
     +  0.63174176E-02,-0.41761557E-02,-0.34239508E-01,-0.65371767E-02,
     + -0.25721467E-02, 0.12215396E-02, 0.11396506E-02, 0.85580323E-04,
     +  0.39144433E-02,-0.53854357E-02,-0.64966135E-03,-0.28836014E-03,
     + -0.28670826E-02,-0.23121845E-01,-0.87125609E-02, 0.56618997E-02,
     + -0.42514415E-02,-0.76664658E-02, 0.17310581E-02,-0.10072044E-02,
     +  0.10347315E-02, 0.25779142E-02,-0.66971676E-02,-0.26783843E-02,
     + -0.22034249E-02, 0.21638337E-03,-0.24572797E-02,-0.63457107E-03,
     + -0.14903315E-02,-0.18360305E-02,-0.46224361E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x32        
     6  +coeff(  6)            *x42    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)    *x22            
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)    *x21*x32        
     8  +coeff( 17)    *x21    *x42    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)    *x22    *x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x23            
     7  +coeff( 25)            *x45    
     8  +coeff( 26)    *x22    *x41*x51
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 27)        *x31*x41    
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)        *x31    *x52
     3  +coeff( 30)    *x21*x31*x41    
     4  +coeff( 31)        *x31*x44    
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)            *x43    
     7  +coeff( 34)*x11    *x31        
     8  +coeff( 35)        *x32*x43    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 36)*x11*x21        *x51
     1  +coeff( 37)    *x22        *x52
     2  +coeff( 38)    *x23        *x52
     3  +coeff( 39)    *x22    *x42*x53
c
      return
      end
      function p_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2157484E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.58519212E-03, 0.83471817E-03, 0.19094884E-01, 0.10380418E-01,
     + -0.63186157E-02, 0.30966003E-02, 0.20771793E-02,-0.40703616E-03,
     +  0.18394530E-01,-0.25050731E-02, 0.36338030E-02,-0.59086415E-02,
     + -0.25375636E-01,-0.29958833E-02,-0.35112319E-02, 0.26019479E-02,
     +  0.13352458E-01, 0.41364171E-02, 0.41788500E-02,-0.18972713E-02,
     +  0.58098766E-03,-0.57365972E-03,-0.17904476E-02, 0.25273249E-02,
     +  0.19108355E-02, 0.13954056E-02,-0.97940199E-03, 0.14733317E-02,
     +  0.51242631E-03, 0.84396661E-03, 0.48239724E-03, 0.29854905E-02,
     + -0.10040066E-02, 0.12651550E-03, 0.27267719E-03,-0.28662058E-03,
     +  0.57112102E-04, 0.23530367E-03,-0.10066436E-02,-0.17037687E-02,
     +  0.38306590E-02, 0.44451314E-02, 0.12533482E-02, 0.53217018E-03,
     + -0.55348891E-03, 0.80422440E-04,-0.45797580E-04,-0.10814783E-03,
     +  0.82120253E-03,-0.42104814E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x22    *x41    
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)        *x31    *x52
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)        *x31*x42    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)                *x53
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)    *x22*x33    *x51
     8  +coeff( 35)        *x32*x41    
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)*x11    *x32        
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x22*x31*x41    
     6  +coeff( 42)    *x22    *x42    
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)    *x22*x32    *x52
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 45)*x11*x22*x33        
     1  +coeff( 46)*x11                
     2  +coeff( 47)*x11            *x51
     3  +coeff( 48)        *x32    *x51
     4  +coeff( 49)    *x22*x32        
     5  +coeff( 50)    *x21    *x43    
c
      return
      end
      function l_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4740702E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16873382E-02,-0.27926224E+00,-0.31160258E-01, 0.11454606E-01,
     + -0.36277626E-01, 0.50801482E-01,-0.17839808E-01, 0.25639333E-01,
     + -0.12038626E-01,-0.88500138E-02,-0.16579079E-01,-0.57715848E-02,
     + -0.39442889E-01,-0.16638924E-02,-0.21885149E-02, 0.27587388E-01,
     + -0.32223330E-02,-0.15871011E-01, 0.17862869E-02,-0.55003716E-02,
     + -0.45790407E-03, 0.62404532E-03, 0.27440106E-02,-0.79626503E-03,
     +  0.51316149E-02,-0.65713576E-02, 0.20432789E-01, 0.19424834E-03,
     + -0.20063827E-02, 0.12257225E-02, 0.23067722E-02, 0.14318291E-02,
     +  0.12572483E-02,-0.64840280E-02,-0.22806555E-03,-0.22422934E-01,
     +  0.26130190E-03, 0.10126340E-02, 0.15974847E-02,-0.35564634E-02,
     +  0.14276435E-02,-0.10852116E-02,-0.47587766E-02,-0.21934153E-02,
     + -0.19936296E-02, 0.11230313E-02,-0.49471683E-02,-0.86441805E-03,
     + -0.42847809E-02,-0.37136883E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23*x31        
     8  +coeff(  8)    *x21*x31        
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x41*x51
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)    *x21*x33*x41    
     2  +coeff( 20)    *x23*x33*x41    
     3  +coeff( 21)        *x32        
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x21*x32        
     8  +coeff( 26)    *x22    *x41    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)            *x43    
     2  +coeff( 29)    *x22        *x51
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)                *x53
     6  +coeff( 33)*x11    *x33        
     7  +coeff( 34)*x11*x22    *x41    
     8  +coeff( 35)    *x22*x33        
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 36)    *x23    *x42    
     1  +coeff( 37)        *x33*x42    
     2  +coeff( 38)        *x31        
     3  +coeff( 39)            *x41    
     4  +coeff( 40)    *x22*x31        
     5  +coeff( 41)        *x31*x41*x51
     6  +coeff( 42)*x11*x21    *x41    
     7  +coeff( 43)    *x22    *x42    
     8  +coeff( 44)    *x23        *x51
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 45)*x11*x22*x31        
     1  +coeff( 46)*x11    *x31    *x52
     2  +coeff( 47)    *x23*x32        
     3  +coeff( 48)*x12*x23            
     4  +coeff( 49)    *x22*x33*x41    
     5  +coeff( 50)        *x31    *x51
c
      return
      end
      function x_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 14)
      data ncoeff/ 13/
      data avdat/  0.1199169E+00/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.61399904E-02,-0.20119749E-01,-0.34063052E-01,-0.20369707E-04,
     + -0.44569519E-04, 0.69311274E-04,-0.14518814E-04, 0.32548476E-04,
     +  0.29814788E-04, 0.19118565E-04, 0.74164936E-05,-0.52011101E-05,
     + -0.76125880E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      x_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)        *x31*x41    
      x_sp_sen    =x_sp_sen    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x24*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)        *x32        
     4  +coeff( 13)*x11*x21            
c
      return
      end
      function t_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/  0.9713093E-01/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.34517851E-02,-0.20394729E-04,-0.61169069E-03,-0.27432231E-01,
     + -0.24466933E-03,-0.46682992E-03,-0.28729939E-03,-0.23643364E-03,
     +  0.54845418E-03, 0.14793029E-05,-0.87229499E-04, 0.24362477E-03,
     + -0.27637541E-04,-0.18442812E-03, 0.32987425E-03,-0.79926867E-04,
     + -0.23493242E-05, 0.14127679E-03, 0.32737880E-04, 0.60072729E-04,
     + -0.16622203E-03, 0.23572923E-03, 0.21346210E-03, 0.20691441E-04,
     +  0.11335528E-04, 0.59894523E-04, 0.49216436E-04, 0.24118122E-04,
     + -0.32582255E-04, 0.33051823E-04, 0.71351220E-04,-0.65407810E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      t_sp_sen    =t_sp_sen    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x22*x33        
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)*x11                
      t_sp_sen    =t_sp_sen    
     9  +coeff( 18)    *x21*x31        
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)        *x31    *x51
     7  +coeff( 25)                *x52
     8  +coeff( 26)    *x21*x31*x41    
      t_sp_sen    =t_sp_sen    
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)    *x22*x32        
     5  +coeff( 32)*x11*x23            
c
      return
      end
      function y_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1154357E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11073431E-02,-0.26216899E-05,-0.44813441E-05, 0.60342641E-05,
     +  0.66440858E-01, 0.49854005E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)*x11                
c
      return
      end
      function p_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.1003938E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.71856362E-03, 0.50688326E-01,-0.12312130E-02,-0.63126319E-03,
     + -0.16688246E-03,-0.12469997E-03, 0.20115524E-03, 0.40474802E-03,
     + -0.91656235E-04, 0.57605757E-04, 0.15482117E-03,-0.22819002E-03,
     + -0.18309843E-03, 0.29355532E-03, 0.22359542E-04,-0.53987867E-04,
     +  0.89238216E-04,-0.66314526E-04, 0.78969686E-04, 0.17487146E-03,
     + -0.45379005E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      p_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)            *x41    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      p_sp_sen    =p_sp_sen    
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x22*x31        
      p_sp_sen    =p_sp_sen    
     9  +coeff( 18)    *x21*x32        
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)*x11    *x33        
c
      return
      end
      function l_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1428167E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12576432E-02, 0.17548436E-02, 0.30322545E-02,-0.17193785E-02,
     + -0.42714091E-05,-0.45570897E-03, 0.36880338E-05, 0.18665110E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21            
     8  +coeff(  8)                *x51
c
      return
      end
      function x_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1811226E+00/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.82209995E-02,-0.20643258E-01,-0.47868252E-01,-0.85000711E-03,
     +  0.10102852E-02, 0.29975115E-03, 0.15242514E-02,-0.50048926E-03,
     +  0.43495325E-03, 0.14767141E-03,-0.24220918E-03,-0.45727831E-03,
     +  0.13058543E-02, 0.10229647E-05, 0.43833817E-04, 0.93276345E-03,
     + -0.38987841E-03,-0.77274442E-03, 0.24115854E-03,-0.12262755E-03,
     + -0.36636347E-03, 0.52145992E-04, 0.29901389E-03, 0.29076840E-03,
     + -0.21390388E-03, 0.22540640E-03, 0.15036118E-02,-0.11923812E-03,
     +  0.78962614E-04,-0.31559821E-03,-0.92723159E-04,-0.22327098E-03,
     +  0.11288867E-03,-0.20185624E-03, 0.28249980E-04, 0.77792007E-04,
     +  0.11535065E-03,-0.29474779E-03, 0.74185699E-03, 0.48976904E-03,
     + -0.21781854E-03,-0.91156154E-03,-0.79308625E-03,-0.20033980E-03,
     +  0.28165035E-04,-0.15899714E-04,-0.31763888E-04,-0.11783482E-03,
     + -0.35817397E-04,-0.31082549E-04, 0.39629595E-03, 0.12448702E-03,
     +  0.80730104E-04, 0.86484055E-04,-0.80915917E-04, 0.98199635E-04,
     +  0.11543326E-03, 0.88596258E-04, 0.38989814E-03,-0.10934273E-03,
     + -0.13693364E-03, 0.30533316E-04, 0.26290965E-04,-0.24541454E-04,
     +  0.18361196E-04, 0.12269578E-04,-0.78722014E-05, 0.16277196E-04,
     +  0.36499389E-05,-0.16417765E-04,-0.34507935E-04,-0.19812436E-04,
     + -0.12331548E-03, 0.11498803E-04, 0.24009902E-04, 0.78480634E-05,
     + -0.20296102E-04, 0.24397250E-04,-0.12032590E-03,-0.11188458E-03,
     + -0.13090619E-04, 0.86549269E-04,-0.97436759E-05, 0.22458789E-04,
     + -0.79682519E-04,-0.43334891E-03,-0.32969032E-03, 0.32023207E-04,
     + -0.13361868E-04,-0.25881911E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x24*x31        
      x_sp_sm     =x_sp_sm     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x21*x33        
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x23    *x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 18)    *x24    *x41    
     1  +coeff( 19)    *x21*x31        
     2  +coeff( 20)        *x32        
     3  +coeff( 21)        *x31*x41    
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x23*x31        
     8  +coeff( 26)*x11*x21    *x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)        *x32*x42    
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)*x11*x23*x31        
     5  +coeff( 32)    *x24*x32        
     6  +coeff( 33)    *x24*x31*x43    
     7  +coeff( 34)            *x42    
     8  +coeff( 35)        *x31    *x51
      x_sp_sm     =x_sp_sm     
     9  +coeff( 36)    *x21*x32        
     1  +coeff( 37)*x11*x21*x31        
     2  +coeff( 38)    *x23*x31*x41    
     3  +coeff( 39)    *x22*x31*x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)*x11*x23    *x41    
     6  +coeff( 42)    *x24*x31*x41    
     7  +coeff( 43)    *x24    *x42    
     8  +coeff( 44)    *x24*x32*x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 45)*x11                
     1  +coeff( 46)    *x21        *x51
     2  +coeff( 47)    *x22        *x51
     3  +coeff( 48)*x11*x22            
     4  +coeff( 49)        *x31*x42    
     5  +coeff( 50)            *x43    
     6  +coeff( 51)    *x22*x32        
     7  +coeff( 52)    *x21*x31*x42    
     8  +coeff( 53)    *x21    *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 54)        *x31*x43    
     1  +coeff( 55)    *x23*x32        
     2  +coeff( 56)*x11*x21*x31*x41    
     3  +coeff( 57)*x11*x21    *x42    
     4  +coeff( 58)*x11*x24            
     5  +coeff( 59)    *x22*x32*x41    
     6  +coeff( 60)    *x22*x31*x43    
     7  +coeff( 61)*x11*x23    *x42    
     8  +coeff( 62)    *x21*x32*x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 63)*x11        *x41    
     1  +coeff( 64)        *x32*x41    
     2  +coeff( 65)        *x31*x41*x51
     3  +coeff( 66)            *x42*x51
     4  +coeff( 67)*x11*x21        *x51
     5  +coeff( 68)    *x23        *x51
     6  +coeff( 69)*x11    *x32        
     7  +coeff( 70)    *x22*x31    *x51
     8  +coeff( 71)    *x22    *x41*x51
      x_sp_sm     =x_sp_sm     
     9  +coeff( 72)*x12*x21            
     1  +coeff( 73)*x11*x22    *x41    
     2  +coeff( 74)*x11*x21*x32        
     3  +coeff( 75)    *x22*x33        
     4  +coeff( 76)    *x22*x32    *x51
     5  +coeff( 77)*x11*x22    *x42    
     6  +coeff( 78)*x12*x23            
     7  +coeff( 79)    *x23*x31*x42    
     8  +coeff( 80)    *x23    *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 81)*x11*x24*x31        
     1  +coeff( 82)*x11*x24    *x41    
     2  +coeff( 83)*x11    *x33*x41    
     3  +coeff( 84)    *x22    *x43*x51
     4  +coeff( 85)*x11*x23*x31*x41    
     5  +coeff( 86)    *x24*x31*x42    
     6  +coeff( 87)    *x24    *x43    
     7  +coeff( 88)        *x33*x43    
     8  +coeff( 89)*x11*x21    *x43*x51
      x_sp_sm     =x_sp_sm     
     9  +coeff( 90)    *x22*x31*x41*x53
c
      return
      end
      function t_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1612890E+00/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.46212212E-02,-0.11843542E-02,-0.28738189E-01,-0.34045817E-02,
     +  0.24837183E-02,-0.78865839E-03, 0.24102766E-02, 0.10692427E-02,
     +  0.13653038E-02, 0.13058621E-02, 0.56593434E-03,-0.56236377E-03,
     +  0.66152972E-03,-0.78998401E-03,-0.64013363E-03, 0.17318946E-03,
     +  0.54620480E-03,-0.12476469E-02, 0.17757345E-02, 0.18358649E-02,
     + -0.20701485E-03, 0.10628992E-03, 0.63707161E-03, 0.67252695E-03,
     + -0.19360195E-03, 0.11327402E-03, 0.98484074E-04, 0.51843858E-03,
     + -0.44662401E-03,-0.69869735E-03, 0.84928912E-03, 0.79924089E-03,
     + -0.54401095E-03,-0.23910565E-03,-0.86484433E-04, 0.35785615E-04,
     +  0.86057262E-04,-0.50165130E-04, 0.15602646E-03, 0.10149336E-03,
     + -0.10380759E-03,-0.47937172E-03, 0.79077092E-03,-0.10760709E-03,
     +  0.12251665E-03, 0.10424407E-03,-0.54987019E-03,-0.12746346E-03,
     + -0.65005245E-03,-0.75100292E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21            
      t_sp_sm     =t_sp_sm     
     9  +coeff(  9)    *x22*x31        
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23*x31        
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)                *x52
     8  +coeff( 17)*x11*x21    *x41    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x22*x31*x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)        *x32        
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)*x11*x21*x31        
      t_sp_sm     =t_sp_sm     
     9  +coeff( 27)*x11    *x32        
     1  +coeff( 28)    *x22*x32        
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)    *x23    *x42    
     4  +coeff( 31)    *x22*x31*x42    
     5  +coeff( 32)    *x22    *x43    
     6  +coeff( 33)*x11*x23    *x41    
     7  +coeff( 34)    *x23*x31*x43    
     8  +coeff( 35)    *x21        *x51
      t_sp_sm     =t_sp_sm     
     9  +coeff( 36)        *x31    *x51
     1  +coeff( 37)    *x21*x32        
     2  +coeff( 38)*x11*x22            
     3  +coeff( 39)        *x31*x43    
     4  +coeff( 40)    *x23        *x51
     5  +coeff( 41)    *x22    *x41*x51
     6  +coeff( 42)    *x23*x31*x41    
     7  +coeff( 43)    *x22*x32*x41    
     8  +coeff( 44)        *x32*x43    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 45)    *x22*x32    *x51
     1  +coeff( 46)*x11*x23        *x52
     2  +coeff( 47)    *x22*x32*x43    
     3  +coeff( 48)*x12*x23*x32        
     4  +coeff( 49)*x12*x21*x32*x41*x52
     5  +coeff( 50)*x12*x21*x32*x42*x52
c
      return
      end
      function y_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  9)
      data ncoeff/  8/
      data avdat/ -0.1638224E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14781052E-02,-0.46219102E-04,-0.96269599E-04, 0.44634617E-05,
     +  0.90653002E-01, 0.49458798E-02,-0.10970745E-02,-0.50142780E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21*x31        
c
      return
      end
      function p_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/ -0.1034744E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.47822454E-03, 0.49160589E-01,-0.45956848E-02,-0.22179100E-02,
     +  0.20576830E-02, 0.83483994E-03, 0.30482931E-02,-0.50366018E-03,
     +  0.93237881E-03, 0.15970382E-02,-0.22223483E-03,-0.24502695E-03,
     + -0.26137917E-03, 0.43278045E-03,-0.17096923E-02,-0.16838197E-02,
     +  0.39843959E-03, 0.14158576E-02, 0.15423821E-02, 0.10550923E-03,
     +  0.17889842E-03,-0.23350159E-03,-0.45029196E-03, 0.44925077E-03,
     + -0.95655625E-04, 0.78427271E-04, 0.78135832E-04,-0.18624462E-05,
     +  0.21463561E-03,-0.11407232E-03,-0.15663095E-03, 0.14444000E-03,
     +  0.36040071E-03,-0.10835053E-03, 0.12189511E-03, 0.30155547E-03,
     + -0.25247564E-03,-0.19887484E-03,-0.10626677E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)            *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23*x31        
     2  +coeff( 11)    *x21*x31*x42    
     3  +coeff( 12)        *x31        
     4  +coeff( 13)*x11                
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      p_sp_sm     =p_sp_sm     
     9  +coeff( 18)    *x23*x31*x41    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)    *x21        *x51
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x32        
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)*x11    *x31        
     8  +coeff( 26)    *x21    *x41*x51
      p_sp_sm     =p_sp_sm     
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x21*x32*x41    
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)        *x32*x42    
     4  +coeff( 31)*x11*x23            
     5  +coeff( 32)*x11*x22*x31        
     6  +coeff( 33)    *x23*x32        
     7  +coeff( 34)    *x22*x32*x41    
     8  +coeff( 35)*x11*x22    *x42    
      p_sp_sm     =p_sp_sm     
     9  +coeff( 36)    *x22*x33*x41    
     1  +coeff( 37)    *x21*x32*x43    
     2  +coeff( 38)        *x33*x43    
     3  +coeff( 39)    *x23        *x53
c
      return
      end
      function l_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/ -0.2042680E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17765037E-02, 0.18157605E-02, 0.47169193E-02, 0.11344279E-03,
     + -0.24323021E-02, 0.18977624E-04,-0.62240777E-03,-0.33386874E-04,
     + -0.28928720E-04,-0.33463803E-04,-0.15666714E-04,-0.44999219E-04,
     + -0.56297391E-04, 0.13776565E-04,-0.22812401E-04, 0.10718090E-04,
     + -0.68723257E-05, 0.24695297E-04,-0.82467332E-05, 0.22688248E-04,
     + -0.45932211E-04,-0.33507768E-04,-0.29305343E-05,-0.28852430E-05,
     + -0.10894526E-04,-0.76481447E-05,-0.54334223E-05, 0.42369138E-05,
     + -0.55054215E-05,-0.14737119E-04, 0.87650660E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21            
      l_sp_sm     =l_sp_sm     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x32        
     8  +coeff( 17)                *x52
      l_sp_sm     =l_sp_sm     
     9  +coeff( 18)    *x23            
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x22*x31*x41    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11                
     6  +coeff( 24)        *x31    *x51
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)    *x21    *x42    
      l_sp_sm     =l_sp_sm     
     9  +coeff( 27)        *x31*x42    
     1  +coeff( 28)*x11*x22            
     2  +coeff( 29)*x11*x21*x31        
     3  +coeff( 30)    *x22*x32        
     4  +coeff( 31)*x11*x23            
c
      return
      end
      function x_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2763546E+00/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10956788E-01,-0.21334158E-01,-0.62454261E-01,-0.34068886E-02,
     +  0.39249989E-02, 0.10032493E-02, 0.38789201E-02, 0.10505909E-02,
     +  0.49492036E-03,-0.67716127E-03, 0.22725244E-02,-0.13651216E-02,
     +  0.57962036E-03,-0.95892424E-03, 0.18007576E-03, 0.54791238E-03,
     + -0.95266773E-03, 0.36048032E-02, 0.32306649E-02,-0.19551557E-02,
     +  0.21414034E-03,-0.60019526E-03, 0.14303140E-03,-0.13327527E-03,
     +  0.75031939E-03, 0.75184391E-03, 0.27666768E-03,-0.51259110E-03,
     + -0.32908807E-03,-0.12138417E-02,-0.79359184E-03,-0.52229391E-03,
     + -0.18702135E-02, 0.31654249E-03, 0.74823714E-04,-0.26855810E-03,
     +  0.56133063E-04, 0.18218133E-03,-0.26853607E-03, 0.87255903E-03,
     +  0.28731365E-03, 0.24614960E-03,-0.72834757E-03, 0.29808623E-03,
     +  0.14714815E-02, 0.10633927E-02,-0.20831588E-02,-0.35883876E-03,
     + -0.59006368E-04, 0.29724842E-03, 0.22370630E-03,-0.17646464E-03,
     +  0.18456741E-03, 0.76441705E-03,-0.22491840E-03,-0.46329148E-03,
     + -0.35204456E-03,-0.33584944E-03, 0.78909536E-04, 0.31134678E-04,
     +  0.28360382E-04,-0.27078291E-04, 0.48917511E-04, 0.12930260E-04,
     + -0.77048382E-04,-0.47396177E-04, 0.53905140E-04,-0.32430605E-03,
     +  0.24647583E-04, 0.75437194E-04,-0.93637187E-04, 0.60632818E-04,
     + -0.27534147E-03,-0.26458115E-03, 0.10147189E-03, 0.21732374E-03,
     +  0.61299797E-04,-0.19894862E-03,-0.75508462E-03,-0.67226408E-03,
     +  0.30316798E-04,-0.40738440E-04, 0.38748789E-04,-0.88150837E-05,
     + -0.42324064E-04,-0.90331714E-05, 0.32736385E-04,-0.25036910E-04,
     + -0.15052101E-03, 0.39920396E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21    *x41    
      x_sp_sex    =x_sp_sex    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11*x21    *x41    
     8  +coeff( 17)    *x23    *x41    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 18)    *x22*x31*x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x24    *x41    
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)            *x42    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)    *x21    *x42    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 27)*x11*x21*x31        
     1  +coeff( 28)    *x23*x31        
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)    *x24*x31        
     4  +coeff( 31)    *x23    *x42    
     5  +coeff( 32)*x11*x23    *x41    
     6  +coeff( 33)    *x24    *x42    
     7  +coeff( 34)    *x24*x31*x43    
     8  +coeff( 35)*x11                
      x_sp_sex    =x_sp_sex    
     9  +coeff( 36)        *x32        
     1  +coeff( 37)        *x31    *x51
     2  +coeff( 38)    *x21*x32        
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)    *x22*x32        
     5  +coeff( 41)        *x31*x43    
     6  +coeff( 42)*x11*x21*x31*x41    
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)*x11*x21    *x42    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 45)    *x22*x31*x42    
     1  +coeff( 46)    *x22    *x43    
     2  +coeff( 47)    *x24*x31*x41    
     3  +coeff( 48)    *x24*x32*x41    
     4  +coeff( 49)    *x21        *x51
     5  +coeff( 50)    *x21*x31*x42    
     6  +coeff( 51)    *x21    *x43    
     7  +coeff( 52)    *x23*x32        
     8  +coeff( 53)*x11*x24            
      x_sp_sex    =x_sp_sex    
     9  +coeff( 54)    *x22*x32*x41    
     1  +coeff( 55)*x11*x23*x31        
     2  +coeff( 56)    *x24*x32        
     3  +coeff( 57)    *x22*x31*x43    
     4  +coeff( 58)*x11*x23    *x42    
     5  +coeff( 59)*x11        *x41    
     6  +coeff( 60)        *x31*x41*x51
     7  +coeff( 61)            *x42*x51
     8  +coeff( 62)*x11*x21        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 63)    *x23        *x51
     1  +coeff( 64)*x11    *x32        
     2  +coeff( 65)    *x22    *x41*x51
     3  +coeff( 66)*x12*x21            
     4  +coeff( 67)    *x21*x32*x41    
     5  +coeff( 68)*x11*x22    *x41    
     6  +coeff( 69)*x11*x21*x32        
     7  +coeff( 70)    *x22*x33        
     8  +coeff( 71)*x11*x22    *x42    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 72)*x12*x23            
     1  +coeff( 73)    *x23*x31*x42    
     2  +coeff( 74)    *x23    *x43    
     3  +coeff( 75)*x11*x24*x31        
     4  +coeff( 76)*x11*x24    *x41    
     5  +coeff( 77)    *x22    *x43*x51
     6  +coeff( 78)*x11*x23*x31*x41    
     7  +coeff( 79)    *x24*x31*x42    
     8  +coeff( 80)    *x24    *x43    
      x_sp_sex    =x_sp_sex    
     9  +coeff( 81)        *x33*x43    
     1  +coeff( 82)*x11*x21    *x43*x51
     2  +coeff( 83)*x11    *x31        
     3  +coeff( 84)    *x21    *x41*x51
     4  +coeff( 85)            *x43    
     5  +coeff( 86)                *x53
     6  +coeff( 87)*x11        *x42    
     7  +coeff( 88)    *x22*x31    *x51
     8  +coeff( 89)*x11*x22*x31        
      x_sp_sex    =x_sp_sex    
     9  +coeff( 90)        *x33*x41    
c
      return
      end
      function t_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2224501E+00/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.58545186E-02, 0.19266934E-02,-0.12634640E-02,-0.29238068E-01,
     + -0.64482912E-02, 0.59437780E-02, 0.28867109E-02, 0.10029743E-02,
     + -0.14286792E-02, 0.14480886E-02, 0.12701112E-02, 0.31737273E-03,
     +  0.29952378E-04,-0.51325321E-03, 0.54500636E-03, 0.63615135E-03,
     + -0.93544478E-03,-0.95542637E-03, 0.78771182E-03, 0.80197223E-03,
     + -0.36403324E-03, 0.53284079E-03,-0.12091510E-02, 0.21944577E-02,
     + -0.71364373E-03, 0.28856238E-03, 0.11950361E-03,-0.24376602E-03,
     +  0.11681785E-03,-0.23758350E-03, 0.46559871E-03, 0.17294112E-02,
     + -0.14993975E-03,-0.71762438E-03,-0.72316185E-03, 0.48821414E-03,
     + -0.49799687E-03,-0.55968808E-03, 0.32263368E-04, 0.85236527E-04,
     + -0.10314193E-03, 0.84652384E-04,-0.58884118E-04, 0.11519585E-03,
     +  0.12137384E-03, 0.42486325E-04,-0.10574621E-03,-0.12818727E-03,
     +  0.10243048E-03,-0.89009198E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)*x11*x21            
      t_sp_sex    =t_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22*x31        
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)        *x31*x42    
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)    *x22*x31*x42    
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)        *x31*x41    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 18)            *x42    
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)*x11*x21    *x41    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)*x11*x23            
     8  +coeff( 26)    *x22*x33*x41    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 27)*x11                
     1  +coeff( 28)        *x32        
     2  +coeff( 29)            *x41*x51
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)    *x22*x32        
     5  +coeff( 32)    *x22*x31*x41    
     6  +coeff( 33)    *x21        *x53
     7  +coeff( 34)    *x23*x31*x41    
     8  +coeff( 35)    *x23    *x42    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)*x11*x22*x31*x41*x52
     3  +coeff( 39)*x12                
     4  +coeff( 40)    *x21*x32        
     5  +coeff( 41)    *x21    *x41*x51
     6  +coeff( 42)*x11*x21*x31        
     7  +coeff( 43)*x11*x21        *x51
     8  +coeff( 44)    *x21    *x43    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 45)    *x23        *x51
     1  +coeff( 46)        *x33    *x51
     2  +coeff( 47)    *x22    *x41*x51
     3  +coeff( 48)    *x21    *x42*x51
     4  +coeff( 49)    *x21*x31    *x52
     5  +coeff( 50)*x11*x22*x31        
c
      return
      end
      function y_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.2158495E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.16213538E-02,-0.89662040E-04,-0.20222079E-03, 0.17569961E-04,
     +  0.11355332E+00, 0.48785456E-02,-0.24439837E-02,-0.10454732E-02,
     +  0.65206323E-03, 0.17029969E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21*x31        
      y_sp_sex    =y_sp_sex    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x23            
c
      return
      end
      function p_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 46)
      data ncoeff/ 45/
      data avdat/ -0.1245795E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.42967338E-03, 0.46958644E-01,-0.56358888E-02,-0.27071775E-02,
     +  0.20434251E-02,-0.38603129E-03, 0.87196229E-03, 0.38345088E-02,
     + -0.59096422E-03, 0.11016635E-02,-0.22587262E-02,-0.18402447E-02,
     +  0.40435526E-03, 0.18677744E-02,-0.28640375E-03, 0.21177075E-03,
     +  0.17468713E-03, 0.51210303E-03, 0.21571119E-02, 0.18222861E-02,
     +  0.43152504E-04,-0.27404583E-03,-0.60350122E-03, 0.10744778E-03,
     + -0.12811313E-03, 0.55157125E-03, 0.12976509E-02, 0.98600925E-04,
     + -0.45799468E-04, 0.55195556E-04, 0.83627070E-04, 0.20792433E-03,
     + -0.99718200E-04,-0.13634261E-03,-0.10256811E-04,-0.71039642E-04,
     +  0.60381461E-03, 0.88574382E-03, 0.62673091E-04, 0.11657290E-03,
     +  0.30161694E-03,-0.18025801E-03, 0.26453685E-03,-0.11407798E-02,
     + -0.14799309E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21*x31        
     5  +coeff(  5)    *x23            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23    *x41    
      p_sp_sex    =p_sp_sex    
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23*x31        
     6  +coeff( 15)        *x31        
     7  +coeff( 16)    *x21        *x51
     8  +coeff( 17)*x11*x21            
      p_sp_sex    =p_sp_sex    
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x23*x31*x41    
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)                *x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21*x32        
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)*x11    *x33        
     8  +coeff( 26)*x11*x22    *x41    
      p_sp_sex    =p_sp_sex    
     9  +coeff( 27)    *x21*x31*x43    
     1  +coeff( 28)    *x21*x32*x41*x51
     2  +coeff( 29)    *x22        *x51
     3  +coeff( 30)    *x21*x31    *x51
     4  +coeff( 31)    *x21    *x41*x51
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)        *x32*x42    
     7  +coeff( 34)*x11*x23            
     8  +coeff( 35)*x11*x22*x31        
      p_sp_sex    =p_sp_sex    
     9  +coeff( 36)*x11    *x31    *x52
     1  +coeff( 37)    *x23*x32        
     2  +coeff( 38)    *x21*x32*x42    
     3  +coeff( 39)*x11*x21*x33        
     4  +coeff( 40)*x11*x22    *x42    
     5  +coeff( 41)    *x22*x33*x41    
     6  +coeff( 42)        *x33*x43    
     7  +coeff( 43)*x11*x22*x33        
     8  +coeff( 44)    *x23*x32*x42    
      p_sp_sex    =p_sp_sex    
     9  +coeff( 45)    *x23*x31*x43    
c
      return
      end
      function l_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/ -0.2982771E-02/
      data xmin/
     1 -0.49991E-02,-0.52045E-01,-0.19982E-01,-0.27969E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.52006E-01, 0.19991E-01, 0.24958E-01, 0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25527810E-02, 0.19364314E-02, 0.74662627E-02, 0.60481677E-03,
     + -0.34597879E-02, 0.88863897E-04,-0.76061540E-03,-0.15957170E-03,
     + -0.12011458E-03,-0.83260442E-04,-0.20755499E-03,-0.13357858E-03,
     + -0.39779181E-04, 0.95968608E-04,-0.15873586E-03, 0.64255750E-04,
     + -0.86620181E-04, 0.33449098E-04, 0.40017774E-04,-0.51655712E-04,
     +  0.11528462E-03,-0.20537566E-03,-0.18506426E-03,-0.11891045E-04,
     +  0.16351765E-04,-0.83644572E-05,-0.83214465E-04,-0.75721808E-04,
     + -0.19248891E-05, 0.16532149E-04,-0.15224146E-04,-0.54405831E-04,
     +  0.46132780E-04,-0.11947868E-04, 0.29450973E-05,-0.14335408E-04,
     +  0.14882879E-04, 0.71817754E-04, 0.72902039E-04,-0.66532200E-04,
     + -0.48960621E-04, 0.45924819E-04, 0.10023146E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21            
      l_sp_sex    =l_sp_sex    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)    *x21*x31        
      l_sp_sex    =l_sp_sex    
     9  +coeff( 18)        *x32        
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x22*x31*x41    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)*x11                
     7  +coeff( 25)    *x21        *x51
     8  +coeff( 26)        *x31    *x51
      l_sp_sex    =l_sp_sex    
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x21    *x42    
     2  +coeff( 29)        *x31*x42    
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)*x11*x21*x31        
     5  +coeff( 32)    *x22*x32        
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)    *x21*x32        
     8  +coeff( 35)            *x43    
      l_sp_sex    =l_sp_sex    
     9  +coeff( 36)    *x23        *x51
     1  +coeff( 37)    *x22    *x41*x51
     2  +coeff( 38)    *x23*x31*x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)    *x22*x31*x42    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)*x11*x23    *x41    
     7  +coeff( 43)*x11*x21        *x53
c
      return
      end
      function x_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3991251E-01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19601865E-01, 0.92578921E-02,-0.30905807E-02, 0.82994294E+00,
     + -0.16693905E-01, 0.64388566E-01,-0.66113457E-01, 0.33020210E-01,
     + -0.27743572E-01, 0.11370937E-01, 0.95327832E-02,-0.14278428E-01,
     + -0.14902183E-01,-0.29068302E-01,-0.39336178E-02,-0.20558783E-02,
     + -0.45656767E-02,-0.12715875E-01,-0.10498507E-01, 0.72072842E-02,
     +  0.24611676E-01,-0.50753588E-03,-0.34169329E-03,-0.55780686E-02,
     +  0.92148141E-03, 0.13755094E-01, 0.71460255E-02,-0.46594217E-02,
     + -0.18734135E-01,-0.39154543E-02, 0.19063387E-02,-0.65367267E-03,
     + -0.94154995E-03,-0.99010649E-03, 0.45704986E-02,-0.17554401E-02,
     + -0.33231343E-02, 0.26995223E-03,-0.26906142E-02, 0.45839287E-02,
     + -0.84804534E-03,-0.56908361E-03,-0.11739754E-02,-0.20378181E-02,
     +  0.85052347E-03, 0.53173867E-02, 0.80720372E-02, 0.11280141E-02,
     +  0.58148947E-03, 0.18470719E-02,-0.59056031E-02, 0.13216719E-02,
     + -0.69798515E-02, 0.13188892E-01, 0.12184985E-01, 0.88023758E-02,
     + -0.28227542E-02, 0.56530456E-02, 0.18086417E-02, 0.47013286E-03,
     +  0.55723201E-03, 0.33373645E-03, 0.15648264E-02, 0.26897175E-03,
     +  0.15999225E-03,-0.90776238E-03, 0.28517866E-02, 0.27490037E-02,
     +  0.33485489E-02, 0.44369330E-02, 0.20498405E-02, 0.28427339E-02,
     +  0.18946287E-02,-0.42820055E-03,-0.12959576E-02, 0.18226462E-02,
     + -0.25134729E-03,-0.11936065E-02,-0.84058277E-03, 0.28408854E-03,
     + -0.17819535E-03,-0.19123794E-03, 0.42452902E-03,-0.53443230E-03,
     +  0.48054202E-03, 0.44889955E-03,-0.71862742E-03, 0.94189250E-03,
     +  0.47075617E-03, 0.97359688E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22            
      x_sp_fp     =x_sp_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23*x31        
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)        *x31        
     8  +coeff( 17)            *x41*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x24            
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x21*x33*x41    
     5  +coeff( 23)        *x32        
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x22    *x41    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x21*x32        
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)    *x21*x31    *x51
     4  +coeff( 31)*x11*x22            
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)    *x24*x31        
     7  +coeff( 34)        *x31    *x51
     8  +coeff( 35)    *x22*x31        
      x_sp_fp     =x_sp_fp     
     9  +coeff( 36)*x11        *x41    
     1  +coeff( 37)        *x31*x41*x51
     2  +coeff( 38)    *x23        *x51
     3  +coeff( 39)    *x22        *x52
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)*x11    *x31        
     6  +coeff( 42)*x11            *x51
     7  +coeff( 43)        *x32    *x51
     8  +coeff( 44)            *x42*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 45)*x11*x21    *x41    
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)    *x22    *x42    
     3  +coeff( 48)    *x21    *x43    
     4  +coeff( 49)    *x21        *x53
     5  +coeff( 50)*x11*x22*x31        
     6  +coeff( 51)    *x24    *x41    
     7  +coeff( 52)*x11*x22        *x51
     8  +coeff( 53)    *x24        *x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 54)    *x23*x31*x41    
     1  +coeff( 55)    *x23    *x42    
     2  +coeff( 56)    *x23    *x41*x51
     3  +coeff( 57)    *x23        *x52
     4  +coeff( 58)    *x24    *x41*x51
     5  +coeff( 59)    *x23*x31*x41*x51
     6  +coeff( 60)            *x41*x52
     7  +coeff( 61)*x11*x21*x31        
     8  +coeff( 62)*x11*x21        *x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 63)    *x22*x32        
     1  +coeff( 64)*x12*x21            
     2  +coeff( 65)    *x21*x31*x42    
     3  +coeff( 66)    *x21    *x42*x51
     4  +coeff( 67)        *x31*x41*x52
     5  +coeff( 68)            *x42*x52
     6  +coeff( 69)    *x23*x32        
     7  +coeff( 70)    *x23*x31    *x51
     8  +coeff( 71)    *x21*x31*x43    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 72)    *x24*x31    *x51
     1  +coeff( 73)    *x23*x31    *x52
     2  +coeff( 74)*x11    *x33    *x52
     3  +coeff( 75)*x11*x23*x33        
     4  +coeff( 76)*x11*x22    *x42*x52
     5  +coeff( 77)        *x32*x41    
     6  +coeff( 78)        *x31*x42    
     7  +coeff( 79)            *x43    
     8  +coeff( 80)        *x31    *x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 81)*x11    *x31    *x51
     1  +coeff( 82)*x11        *x41*x51
     2  +coeff( 83)*x11*x23            
     3  +coeff( 84)    *x21*x31    *x52
     4  +coeff( 85)        *x32    *x52
     5  +coeff( 86)*x11*x21*x31*x41    
     6  +coeff( 87)*x11*x24            
     7  +coeff( 88)    *x22*x31*x42    
     8  +coeff( 89)    *x22*x32    *x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 90)    *x22    *x41*x52
c
c      print *, "We are at the end of x_sp_fp"
      return
      end
      function t_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 36)
      data ncoeff/ 35/
      data avdat/  0.5915433E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27998511E-02,-0.14091826E-01, 0.61458173E-04, 0.11757769E+00,
     +  0.69178492E-02,-0.10434853E-01,-0.17061528E-02, 0.21724768E-02,
     + -0.18134452E-02, 0.69245417E-03,-0.10476375E-02,-0.25082750E-02,
     + -0.12909487E-02,-0.88333141E-03, 0.33270079E-03, 0.26354019E-03,
     + -0.60670526E-03, 0.80785365E-04,-0.52120636E-03, 0.57297220E-04,
     + -0.11602406E-02, 0.55497797E-03, 0.14981178E-02,-0.98205439E-03,
     + -0.42383585E-03,-0.12556506E-03,-0.29955711E-03, 0.44504530E-03,
     + -0.24917058E-03, 0.11017938E-02, 0.57079946E-03, 0.31870900E-03,
     +  0.39949658E-03,-0.33808558E-03, 0.30126882E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11                
     8  +coeff(  8)    *x22            
      t_sp_fp     =t_sp_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x21    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)    *x23*x31        
     7  +coeff( 16)        *x32*x41*x51
     8  +coeff( 17)    *x21*x31        
      t_sp_fp     =t_sp_fp     
     9  +coeff( 18)        *x32        
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)                *x53
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x22*x31*x41*x51
     7  +coeff( 25)        *x31*x41    
     8  +coeff( 26)        *x31    *x51
      t_sp_fp     =t_sp_fp     
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)        *x31*x42*x51
     7  +coeff( 34)    *x22        *x52
     8  +coeff( 35)*x11*x23            
c
      return
      end
      function y_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.3219719E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22484763E-02, 0.42413820E-01,-0.46000322E-02,-0.19319044E-01,
     +  0.35614420E-02,-0.60792564E-03,-0.69519510E-02,-0.19698650E-01,
     + -0.68153054E-01,-0.68500065E-02, 0.43873130E-02, 0.43401133E-01,
     +  0.83600879E-02, 0.29125803E-02, 0.11777391E-01, 0.19177670E-01,
     +  0.27608691E-03, 0.19087141E-02,-0.12085072E-02, 0.89300107E-02,
     +  0.35449859E-01, 0.92968112E-02,-0.95948838E-02, 0.21415099E-03,
     +  0.76763652E-03,-0.18750588E-03, 0.19959926E-02,-0.49751960E-02,
     +  0.16340335E-02, 0.52588312E-02, 0.70851436E-02, 0.75302261E-03,
     + -0.12505609E-02, 0.49914937E-02,-0.36677895E-02, 0.14564143E-02,
     +  0.21847598E-02, 0.25210858E-02,-0.64716004E-02,-0.37798507E-03,
     +  0.20978970E-02,-0.13143525E-02,-0.13246762E-02, 0.90240892E-02,
     +  0.12398713E-01, 0.61131100E-03, 0.29683465E-02, 0.44775042E-02,
     + -0.34279563E-02,-0.49050250E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      y_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sp_fp     =y_sp_fp     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)            *x42*x51
     6  +coeff( 15)    *x22            
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)    *x21*x32        
      y_sp_fp     =y_sp_fp     
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x22    *x41    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x23            
     6  +coeff( 24)            *x45    
     7  +coeff( 25)        *x33    *x52
     8  +coeff( 26)        *x35*x41    
      y_sp_fp     =y_sp_fp     
     9  +coeff( 27)        *x31    *x54
     1  +coeff( 28)        *x31*x41    
     2  +coeff( 29)        *x32*x41    
     3  +coeff( 30)        *x31*x42    
     4  +coeff( 31)            *x43    
     5  +coeff( 32)        *x31*x41*x51
     6  +coeff( 33)*x11        *x41    
     7  +coeff( 34)        *x31    *x52
     8  +coeff( 35)    *x21    *x42    
      y_sp_fp     =y_sp_fp     
     9  +coeff( 36)                *x53
     1  +coeff( 37)    *x21*x31    *x51
     2  +coeff( 38)    *x21    *x41*x51
     3  +coeff( 39)            *x41*x53
     4  +coeff( 40)    *x21*x33*x41    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x21*x31*x41    
     7  +coeff( 43)        *x31    *x53
     8  +coeff( 44)    *x22*x31*x41    
      y_sp_fp     =y_sp_fp     
     9  +coeff( 45)    *x22    *x42    
     1  +coeff( 46)*x11*x21        *x51
     2  +coeff( 47)    *x22*x31    *x51
     3  +coeff( 48)    *x22    *x41*x51
     4  +coeff( 49)    *x23*x31        
     5  +coeff( 50)    *x23    *x41    
c
      return
      end
      function p_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2157776E-02/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.58462261E-03, 0.83458406E-03, 0.19093899E-01, 0.10378784E-01,
     + -0.63183457E-02, 0.30966867E-02, 0.20780840E-02,-0.40713054E-03,
     +  0.18395679E-01,-0.25049581E-02, 0.36333755E-02,-0.59084920E-02,
     + -0.25372114E-01,-0.29984615E-02,-0.35115990E-02, 0.26020079E-02,
     +  0.13353568E-01, 0.41372804E-02, 0.41898522E-02,-0.18972631E-02,
     +  0.58101007E-03,-0.57333242E-03,-0.17911128E-02, 0.25283762E-02,
     +  0.19123912E-02, 0.14037542E-02,-0.97985764E-03, 0.14744105E-02,
     +  0.51438110E-03, 0.84397529E-03, 0.48062860E-03, 0.29859750E-02,
     + -0.10166095E-02, 0.12456752E-03, 0.27286337E-03,-0.28672314E-03,
     +  0.57132867E-04, 0.23537231E-03,-0.10067968E-02,-0.17037931E-02,
     +  0.38304094E-02, 0.44443049E-02, 0.12541673E-02, 0.53499726E-03,
     + -0.55400009E-03, 0.80456492E-04,-0.45655921E-04,-0.10845815E-03,
     +  0.82119048E-03,-0.42074459E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_fp     =p_sp_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x22    *x41    
      p_sp_fp     =p_sp_fp     
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)            *x43    
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)        *x31    *x52
      p_sp_fp     =p_sp_fp     
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)        *x31*x42    
     2  +coeff( 29)    *x21*x31    *x51
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)                *x53
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)    *x22*x33    *x51
     8  +coeff( 35)        *x32*x41    
      p_sp_fp     =p_sp_fp     
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)*x11    *x32        
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x22*x31*x41    
     6  +coeff( 42)    *x22    *x42    
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)    *x22*x32    *x52
      p_sp_fp     =p_sp_fp     
     9  +coeff( 45)*x11*x22*x33        
     1  +coeff( 46)*x11                
     2  +coeff( 47)*x11            *x51
     3  +coeff( 48)        *x32    *x51
     4  +coeff( 49)    *x22*x32        
     5  +coeff( 50)    *x21    *x43    
c
      return
      end
      function l_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1415866E-01/
      data xmin/
     1 -0.49991E-02,-0.50379E-01,-0.19979E-01,-0.27919E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49960E-02, 0.50623E-01, 0.19991E-01, 0.24090E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11090697E-01,-0.27801663E+00,-0.35796918E-01, 0.11526894E-01,
     + -0.36965717E-01, 0.51164970E-01,-0.45446530E-02,-0.39209262E-01,
     + -0.17423632E-01, 0.25637401E-01,-0.18416414E-01,-0.62714978E-02,
     + -0.41701761E-02, 0.66068862E-02,-0.38943958E-01,-0.26777566E-02,
     +  0.58122085E-04,-0.23753908E-02, 0.25197402E-01,-0.34126400E-02,
     + -0.29618275E-04,-0.14107388E-01,-0.13889724E-02, 0.23567851E-02,
     +  0.27928581E-02,-0.63005928E-02, 0.17589867E-01,-0.24778456E-02,
     +  0.16066256E-02, 0.49050557E-02,-0.29785824E-02, 0.99098566E-03,
     + -0.47015096E-02,-0.94505656E-03,-0.19631160E-01,-0.11339951E-03,
     +  0.87909569E-03, 0.80988952E-03, 0.58056810E-03,-0.15918326E-02,
     +  0.12012862E-02, 0.16907010E-02, 0.45894069E-03,-0.31722183E-02,
     +  0.75726205E-03, 0.38352946E-02,-0.95414341E-03,-0.58906963E-02,
     + -0.74858610E-02, 0.21501521E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      l_sp_fp     =l_sp_fp     
     9  +coeff(  9)    *x23*x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23            
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)        *x31    *x51
      l_sp_fp     =l_sp_fp     
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)        *x32*x42    
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)        *x32        
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x21*x32        
     8  +coeff( 26)    *x22    *x41    
      l_sp_fp     =l_sp_fp     
     9  +coeff( 27)    *x21*x31*x41    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)    *x21    *x41*x51
     3  +coeff( 30)            *x42*x51
     4  +coeff( 31)    *x23        *x51
     5  +coeff( 32)*x11    *x33        
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)    *x22*x33        
     8  +coeff( 35)    *x23    *x42    
      l_sp_fp     =l_sp_fp     
     9  +coeff( 36)    *x22*x31*x41*x51
     1  +coeff( 37)        *x33*x41*x51
     2  +coeff( 38)        *x33    *x52
     3  +coeff( 39)            *x43*x52
     4  +coeff( 40)    *x22    *x42*x52
     5  +coeff( 41)        *x31        
     6  +coeff( 42)            *x41    
     7  +coeff( 43)*x11*x21            
     8  +coeff( 44)    *x22*x31        
      l_sp_fp     =l_sp_fp     
     9  +coeff( 45)        *x32    *x51
     1  +coeff( 46)        *x31*x41*x51
     2  +coeff( 47)*x11*x21    *x41    
     3  +coeff( 48)    *x22*x31*x41    
     4  +coeff( 49)    *x22    *x42    
     5  +coeff( 50)    *x21    *x42*x51
c
      return
      end
      function x_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.3243843E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11988E+00, 0.47532E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50856019E-02, 0.69771208E-01, 0.70643641E-01, 0.56444341E-03,
     +  0.17562961E-04, 0.23289346E-02, 0.21732172E-03,-0.53536336E-04,
     + -0.14095762E-03, 0.36204183E-04,-0.11946642E-04, 0.62376307E-05,
     + -0.49464247E-03,-0.15855368E-03,-0.15448334E-03,-0.39811334E-06,
     + -0.90156609E-05,-0.47960320E-05,-0.64237538E-04,-0.95820127E-04,
     +  0.82912673E-04, 0.59935162E-04,-0.12024258E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11            *x51
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x21*x31*x41    
      x_sp_cq1x   =x_sp_cq1x   
     9  +coeff(  9)    *x21        *x52
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x23*x32        
     6  +coeff( 15)    *x23    *x42    
     7  +coeff( 16)        *x32        
     8  +coeff( 17)        *x31*x41    
      x_sp_cq1x   =x_sp_cq1x   
     9  +coeff( 18)                *x52
     1  +coeff( 19)    *x21*x32        
     2  +coeff( 20)*x11        *x42    
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)*x11*x22*x31        
     5  +coeff( 23)    *x21*x33*x41    
c
      return
      end
      function t_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/ -0.9797944E-04/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11988E+00, 0.47532E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.91266364E-03, 0.34572173E-01, 0.28139013E-05,-0.44768609E-01,
     + -0.86370001E-05, 0.54101337E-03, 0.21204236E-02,-0.12683906E-03,
     +  0.13340091E-03,-0.12778303E-03,-0.90867157E-04, 0.96263108E-03,
     +  0.17132499E-04,-0.91675729E-05, 0.27760027E-04,-0.16171447E-02,
     + -0.11447556E-04,-0.31642511E-03, 0.67264823E-05, 0.15420950E-03,
     +  0.12738224E-03,-0.61133396E-04, 0.13395013E-03,-0.20745193E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)*x11            *x51
     8  +coeff(  8)    *x21    *x42    
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21        *x52
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x23*x31*x41    
     4  +coeff( 13)                *x51
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)*x11*x22*x31*x41    
     8  +coeff( 17)        *x31*x41    
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)    *x21*x32*x41    
     4  +coeff( 22)    *x21*x33*x41    
     5  +coeff( 23)    *x23    *x42*x51
     6  +coeff( 24)    *x21*x33*x43    
c
      return
      end
      function y_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.1149580E-01/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11988E+00, 0.47532E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.40565506E-02, 0.73072977E-01, 0.62544182E-01,-0.14997720E-02,
     +  0.28508162E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x31 = x3
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)        *x31    *x51
     5  +coeff(  5)                *x51
c
      return
      end
      function p_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.4893997E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11988E+00, 0.47532E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.33020772E-03, 0.17150077E-05, 0.21382164E-01, 0.39617952E-01,
     + -0.11207436E-02, 0.22462581E-03,-0.46737044E-03,-0.30070113E-03,
     +  0.57680609E-04, 0.12847416E-03, 0.72211515E-04,-0.37919913E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)        *x31    *x51
     6  +coeff(  6)                *x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22*x31        
      p_sp_cq1x   =p_sp_cq1x   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)        *x31    *x52
     3  +coeff( 12)*x11*x23    *x41    
c
      return
      end
      function l_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.5206256E-03/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11988E+00, 0.47532E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50773507E-03,-0.68734204E-04, 0.10355274E-03, 0.23126008E-03,
     +  0.10027197E-05,-0.13935735E-02,-0.20593379E-03, 0.78056019E-07,
     + -0.86204271E-03,-0.11267633E-02, 0.22541952E-04,-0.10314824E-04,
     +  0.18687186E-04, 0.17935133E-02,-0.93208707E-03, 0.57830515E-04,
     +  0.43552198E-04,-0.10484012E-04,-0.40525657E-06,-0.52888146E-04,
     +  0.51057315E-04,-0.25741217E-05,-0.44886815E-05,-0.45331890E-05,
     +  0.13609899E-05, 0.13755039E-05,-0.32064943E-05,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)        *x32    *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x42*x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x12                
     7  +coeff( 16)*x12            *x51
     8  +coeff( 17)*x11                
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)                *x52
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)        *x31*x41*x51
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)        *x31*x43    
     6  +coeff( 24)        *x32    *x52
     7  +coeff( 25)        *x32*x41    
     8  +coeff( 26)        *x31    *x52
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 27)    *x22        *x52
c
      return
      end
      function x_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 73)
      data ncoeff/ 72/
      data avdat/ -0.5095685E+01/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.53157490E-02,-0.17699793E+00, 0.61717346E-01, 0.25856295E-02,
     +  0.85912999E-02,-0.83881151E-02,-0.69519784E-03,-0.33694621E-04,
     +  0.10931958E-03,-0.67694862E-02,-0.13366359E-02, 0.68751011E-04,
     + -0.10428562E-03,-0.95070027E-04, 0.21488452E-02,-0.10614399E-02,
     + -0.18465749E-03, 0.46928684E-03, 0.38327032E-03, 0.80097211E-03,
     +  0.26460187E-03,-0.11472932E-03,-0.18021783E-03, 0.42086822E-03,
     + -0.22579555E-03, 0.14603960E-03, 0.10248648E-02, 0.61061501E-03,
     +  0.12071130E-04, 0.22188060E-02, 0.74118425E-05,-0.68246709E-06,
     + -0.12202212E-03, 0.13851989E-04,-0.19357567E-03, 0.33710308E-04,
     +  0.25856029E-03, 0.13994334E-03,-0.32875643E-04,-0.10523475E-03,
     + -0.91217553E-04,-0.38387843E-04,-0.12738371E-02, 0.24120396E-02,
     + -0.24386002E-03, 0.64976804E-03,-0.93355734E-03, 0.12781710E-04,
     +  0.44646708E-05, 0.12793271E-02, 0.65522909E-04, 0.36131203E-05,
     + -0.16944472E-04,-0.10494928E-04,-0.28943553E-03,-0.13715358E-03,
     + -0.20739744E-04,-0.58110588E-03, 0.31632037E-04,-0.22616675E-03,
     +  0.10449634E-03,-0.89493464E-04,-0.11572105E-04, 0.95379364E-04,
     +  0.20049678E-03,-0.38556720E-03,-0.22159982E-02, 0.43099839E-03,
     +  0.87823038E-03, 0.26506230E-02, 0.20536401E-02,-0.16193165E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)*x11            *x51
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21*x31        
      x_sp_cden   =x_sp_cden   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)    *x21        *x52
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x23*x31*x41    
     8  +coeff( 17)    *x21*x31*x43    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 18)    *x23    *x43    
     1  +coeff( 19)    *x21*x32        
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x21*x32*x41    
     5  +coeff( 23)    *x21    *x43    
     6  +coeff( 24)    *x21*x33*x41    
     7  +coeff( 25)*x11*x22    *x42*x51
     8  +coeff( 26)*x11    *x32        
      x_sp_cden   =x_sp_cden   
     9  +coeff( 27)*x11    *x31*x41    
     1  +coeff( 28)*x11            *x52
     2  +coeff( 29)    *x21*x32    *x53
     3  +coeff( 30)*x11*x24    *x42    
     4  +coeff( 31)            *x41*x51
     5  +coeff( 32)                *x52
     6  +coeff( 33)*x11    *x31        
     7  +coeff( 34)    *x22*x31        
     8  +coeff( 35)*x11        *x41    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 36)    *x21*x31    *x51
     1  +coeff( 37)*x11*x21        *x51
     2  +coeff( 38)    *x23        *x51
     3  +coeff( 39)    *x22*x32        
     4  +coeff( 40)    *x21*x33        
     5  +coeff( 41)    *x21*x32    *x51
     6  +coeff( 42)    *x21        *x53
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)*x11*x22*x31*x41    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 45)    *x21*x32*x42*x51
     1  +coeff( 46)    *x23*x32*x42    
     2  +coeff( 47)*x12*x21*x31*x42    
     3  +coeff( 48)            *x41    
     4  +coeff( 49)    *x21    *x41*x51
     5  +coeff( 50)*x12                
     6  +coeff( 51)    *x24            
     7  +coeff( 52)        *x32*x41    
     8  +coeff( 53)            *x43    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 54)        *x32    *x51
     1  +coeff( 55)    *x23*x31        
     2  +coeff( 56)    *x22*x31*x41    
     3  +coeff( 57)    *x22        *x52
     4  +coeff( 58)*x11*x22    *x41    
     5  +coeff( 59)    *x24    *x41    
     6  +coeff( 60)*x11*x22        *x51
     7  +coeff( 61)    *x23*x32        
     8  +coeff( 62)*x12*x22            
      x_sp_cden   =x_sp_cden   
     9  +coeff( 63)*x11*x24            
     1  +coeff( 64)    *x21*x31*x42*x51
     2  +coeff( 65)    *x23*x33        
     3  +coeff( 66)    *x21*x33*x42    
     4  +coeff( 67)    *x23*x33*x41    
     5  +coeff( 68)*x11*x24*x32        
     6  +coeff( 69)    *x21*x33*x43    
     7  +coeff( 70)*x11*x22*x33*x41    
     8  +coeff( 71)*x11*x22*x31*x43    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 72)*x11*x24*x32*x43    
c
      return
      end
      function t_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/  0.1297449E+01/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11268437E-02, 0.92829995E-01,-0.41797195E-01, 0.71699307E-02,
     + -0.20215723E-02, 0.54125063E-03,-0.43333304E-04, 0.33168049E-03,
     +  0.44817985E-02,-0.25478119E-03,-0.10816216E-03,-0.47801345E-03,
     +  0.27545315E-03, 0.18318990E-06,-0.98171367E-04,-0.32246411E-02,
     +  0.15195942E-03,-0.27935364E-03,-0.62437257E-03, 0.19590181E-03,
     + -0.53685377E-04, 0.89664175E-03,-0.49141141E-04,-0.75578736E-03,
     +  0.14872427E-03,-0.27590550E-04, 0.18695995E-03,-0.34891255E-03,
     +  0.43111481E-03,-0.23978205E-03,-0.60312275E-03,-0.77893503E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x51
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      t_sp_cden   =t_sp_cden   
     9  +coeff(  9)*x11            *x51
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)        *x31        
     6  +coeff( 15)            *x41    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)        *x31*x41*x51
      t_sp_cden   =t_sp_cden   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)*x11*x23            
     2  +coeff( 20)    *x21*x31        
     3  +coeff( 21)    *x21    *x41    
     4  +coeff( 22)            *x42    
     5  +coeff( 23)                *x52
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)        *x31    *x51
      t_sp_cden   =t_sp_cden   
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x23*x32        
     4  +coeff( 31)    *x21*x33*x41    
     5  +coeff( 32)*x11*x22*x31*x41    
c
      return
      end
      function y_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.6570819E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20113306E-02, 0.93240775E-02, 0.89804955E-01,-0.14299111E-02,
     + -0.74729283E-03, 0.42882096E-03, 0.13987158E-03, 0.59260661E-02,
     +  0.99322973E-02, 0.37348121E-02,-0.32807377E-02, 0.10178301E-02,
     + -0.12681177E-02, 0.71486842E-03,-0.82977826E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)        *x31    *x51
      y_sp_cden   =y_sp_cden   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x22*x31        
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)    *x22            
     6  +coeff( 15)        *x31    *x52
c
      return
      end
      function p_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.7766830E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.45502596E-02, 0.50308402E-02,-0.64852782E-01,-0.92247231E-02,
     + -0.22529955E-02,-0.23275472E-01,-0.19483123E-01, 0.34616017E-03,
     +  0.12888947E-03, 0.68534552E-02, 0.80581252E-02,-0.11658904E-01,
     + -0.15056717E-02, 0.27607856E-02,-0.14471369E-01, 0.10398210E-01,
     + -0.69812217E-04,-0.10999387E-03,-0.19980168E-03,-0.51922066E-03,
     +  0.65101716E-02,-0.67177281E-03, 0.10530554E-03,-0.15482797E-02,
     +  0.75285435E-02,-0.12645497E-02,-0.44716729E-03,-0.36638524E-03,
     +  0.13247560E-02,-0.78947241E-05, 0.30255411E-02, 0.56995233E-02,
     + -0.30122134E-02, 0.22140985E-03, 0.14252309E-03, 0.70356741E-03,
     + -0.30413677E-03, 0.79359830E-03,-0.33702259E-03,-0.38186883E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      p_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cden   =p_sp_cden   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)                *x51
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)        *x32*x41*x52
      p_sp_cden   =p_sp_cden   
     9  +coeff( 18)        *x33        
     1  +coeff( 19)        *x32*x41    
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)*x11*x21*x31        
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)                *x52
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)        *x31*x42    
      p_sp_cden   =p_sp_cden   
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)    *x22*x31*x42    
     2  +coeff( 29)*x12*x22    *x41    
     3  +coeff( 30)    *x21        *x51
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)*x11    *x31    *x51
     7  +coeff( 34)*x11        *x41*x51
     8  +coeff( 35)    *x21*x33        
      p_sp_cden   =p_sp_cden   
     9  +coeff( 36)    *x21*x31*x42    
     1  +coeff( 37)    *x22*x31    *x51
     2  +coeff( 38)*x11*x22    *x41    
     3  +coeff( 39)*x12*x21*x31        
     4  +coeff( 40)    *x21*x33*x41    
c
      return
      end
      function l_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.5398299E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12852056E-01, 0.34343985E+00,-0.12069261E+00,-0.29953165E-01,
     + -0.51303455E-02, 0.48684442E-03,-0.10100976E-02,-0.23860317E-02,
     +  0.17616607E-01, 0.16528741E-01,-0.42844918E-02, 0.43604634E-03,
     +  0.72204845E-03, 0.23107023E-02,-0.12137282E-02, 0.70625934E-03,
     + -0.97848615E-03, 0.46380705E-03,-0.18179751E-02,-0.36553065E-02,
     +  0.17232706E-02,-0.35214319E-02,-0.17945534E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)        *x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)        *x31*x41    
      l_sp_cden   =l_sp_cden   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)*x11            *x51
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)    *x21        *x52
      l_sp_cden   =l_sp_cden   
     9  +coeff( 18)            *x41    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)*x11*x22*x32        
c
      return
      end
      function x_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 60)
      data ncoeff/ 59/
      data avdat/ -0.5441129E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.90839993E-02, 0.37612125E+00, 0.26658412E-04, 0.13861294E+00,
     + -0.17901815E+00, 0.33319883E-01, 0.13005653E-03, 0.25534121E-01,
     + -0.58785533E-02,-0.84913522E-03,-0.87943382E-03, 0.22876661E-02,
     + -0.48270330E-03, 0.63144782E-03,-0.16337829E-01, 0.56246235E-02,
     + -0.55377604E-03,-0.31450619E-02,-0.21934947E-02,-0.25329392E-02,
     +  0.35236977E-03, 0.66752117E-02, 0.19129536E-03,-0.11968149E-03,
     + -0.25700885E-03, 0.43296436E-03,-0.61984238E-03, 0.29802840E-03,
     +  0.24564838E-03, 0.40945929E-03,-0.71733125E-03, 0.55215828E-03,
     +  0.20815208E-03, 0.88290387E-03,-0.37308810E-02,-0.11390964E-01,
     + -0.25073544E-02,-0.51925093E-04, 0.44876738E-02, 0.46494129E-03,
     +  0.20491319E-03,-0.17151261E-02,-0.41948509E-03,-0.33970024E-02,
     + -0.45497471E-03,-0.78074879E-03,-0.96681342E-03, 0.25901553E-03,
     +  0.22821453E-03, 0.16834980E-02, 0.10421813E-02, 0.27415324E-02,
     + -0.11591144E-02,-0.20306707E-03, 0.20940333E-03,-0.61785511E-04,
     + -0.11459909E-03,-0.51289279E-03,-0.12198400E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21        *x52
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 18)    *x24            
     1  +coeff( 19)            *x42    
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x23*x31*x41    
     5  +coeff( 23)*x11*x24            
     6  +coeff( 24)*x11    *x33        
     7  +coeff( 25)    *x21*x33*x41    
     8  +coeff( 26)*x11*x22    *x42    
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 27)*x12*x23            
     1  +coeff( 28)        *x31        
     2  +coeff( 29)    *x21*x31        
     3  +coeff( 30)        *x31*x41*x51
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)    *x23*x31        
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)    *x21*x32*x41    
     8  +coeff( 35)    *x21*x32*x42    
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 36)*x11*x22*x31*x41    
     1  +coeff( 37)*x11*x24    *x42    
     2  +coeff( 38)        *x31    *x51
     3  +coeff( 39)    *x23            
     4  +coeff( 40)*x11    *x31        
     5  +coeff( 41)    *x22*x31        
     6  +coeff( 42)    *x23    *x41    
     7  +coeff( 43)    *x22*x32        
     8  +coeff( 44)*x12*x21            
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 45)    *x24        *x51
     1  +coeff( 46)    *x24    *x42*x52
     2  +coeff( 47)*x11    *x32        
     3  +coeff( 48)    *x22    *x42    
     4  +coeff( 49)*x11            *x52
     5  +coeff( 50)*x11*x23            
     6  +coeff( 51)    *x21*x31*x42    
     7  +coeff( 52)*x11*x22    *x41    
     8  +coeff( 53)    *x23*x32        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 54)    *x23        *x52
     1  +coeff( 55)    *x22*x32    *x51
     2  +coeff( 56)*x11            *x53
     3  +coeff( 57)*x12            *x52
     4  +coeff( 58)*x11*x23*x31*x41    
     5  +coeff( 59)        *x31*x43*x52
c
      return
      end
      function t_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/  0.7466303E-03/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.19662841E-02,-0.11150031E+00,-0.38696522E-04, 0.22771409E-01,
     +  0.38742457E-01,-0.11337404E-01,-0.33052528E-03,-0.13957653E-02,
     +  0.72410521E-02, 0.41699252E-03,-0.22508288E-03,-0.40183226E-02,
     +  0.28499888E-03,-0.26679953E-03, 0.45845067E-03,-0.17770245E-02,
     +  0.25485453E-03, 0.69547634E-04, 0.10931690E-03, 0.18867494E-02,
     + -0.42420684E-03,-0.14825322E-03, 0.16696581E-02, 0.87141403E-03,
     + -0.19978052E-03,-0.68357587E-03,-0.35084807E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)            *x43    
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)*x11    *x32*x41    
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)    *x21*x33*x41    
     7  +coeff( 25)        *x33*x42    
     8  +coeff( 26)    *x23    *x42*x51
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 27)*x12*x23    *x41    
c
      return
      end
      function y_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/  0.5992274E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15038303E-01,-0.13945477E+00, 0.14999770E+00,-0.61103529E-02,
     +  0.46869656E-02, 0.89728140E-03,-0.58881595E-03, 0.16462017E-02,
     +  0.27401041E-01, 0.34805752E-01,-0.19882831E-02, 0.47577644E-03,
     + -0.17637987E-01,-0.57209574E-01, 0.32387509E-02, 0.89100283E-02,
     + -0.21421441E-02,-0.16009185E-02,-0.55170734E-02,-0.93468576E-03,
     + -0.90956315E-02,-0.19098578E-01,-0.47344118E-03,-0.76159829E-03,
     + -0.45774556E-02, 0.17261840E-01,-0.21147339E-02, 0.38097030E-03,
     + -0.18967206E-02, 0.18156674E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x22            
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)        *x31    *x52
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff( 18)    *x21*x31    *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)    *x22    *x41    
     5  +coeff( 23)        *x35        
     6  +coeff( 24)        *x32*x41    
     7  +coeff( 25)        *x31*x42    
     8  +coeff( 26)*x11        *x41    
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)*x11            *x51
     2  +coeff( 29)    *x22*x31    *x51
     3  +coeff( 30)    *x23*x31        
c
      return
      end
      function p_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2082444E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29622621E-02,-0.31053862E-01, 0.22163894E-01,-0.94526884E-03,
     + -0.66825259E-03,-0.65080528E-02,-0.99228239E-02, 0.12584173E-03,
     +  0.92345203E-04, 0.36995648E-02, 0.81853000E-02,-0.24903473E-03,
     + -0.95217852E-02,-0.15095869E-04,-0.60955186E-04,-0.47272115E-03,
     +  0.38955096E-03, 0.33216055E-02,-0.76539058E-03, 0.14865078E-02,
     + -0.79697842E-03, 0.19729836E-03, 0.68819107E-04,-0.38573264E-04,
     +  0.94818433E-05,-0.64216543E-03, 0.46986034E-02,-0.14594734E-03,
     +  0.91498927E-03, 0.34824901E-03,-0.18673980E-04, 0.25866975E-02,
     +  0.35931794E-05, 0.38203914E-04,-0.16277978E-03,-0.19039820E-03,
     + -0.17452994E-03,-0.80012594E-03,-0.32895015E-03,-0.29922015E-03,
     + -0.87136046E-04,-0.18395153E-03, 0.11956114E-03,-0.40093658E-03,
     +  0.79429432E-04, 0.40813506E-04, 0.20559978E-03,-0.44562086E-03,
     +  0.14534251E-03, 0.17212494E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x22            
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 18)*x11    *x31        
     1  +coeff( 19)    *x22*x31    *x51
     2  +coeff( 20)    *x21            
     3  +coeff( 21)        *x31*x42    
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)            *x42*x51
     7  +coeff( 25)        *x31    *x52
     8  +coeff( 26)*x11*x21*x31        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x23*x31        
     2  +coeff( 29)    *x23    *x41    
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21        *x51
     5  +coeff( 32)*x11        *x41    
     6  +coeff( 33)    *x23            
     7  +coeff( 34)    *x21*x32        
     8  +coeff( 35)        *x33        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 36)    *x21*x31*x41    
     1  +coeff( 37)            *x41*x52
     2  +coeff( 38)*x11        *x41*x51
     3  +coeff( 39)    *x22    *x41*x51
     4  +coeff( 40)    *x21*x31    *x52
     5  +coeff( 41)        *x31    *x53
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)    *x23        *x52
     8  +coeff( 44)*x11*x23*x31        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 45)        *x32        
     1  +coeff( 46)        *x33    *x51
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)    *x22*x32*x41    
     4  +coeff( 49)*x12        *x42    
     5  +coeff( 50)    *x23    *x41*x51
c
      return
      end
      function l_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/ -0.2946497E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.76944320E-02,-0.48103100E+00, 0.21497597E-03,-0.31803742E-01,
     +  0.18957625E+00,-0.49882296E-01, 0.26138867E-02,-0.27198677E-02,
     +  0.19537401E-02, 0.10590277E-02,-0.43902805E-03, 0.27063169E-01,
     + -0.19740658E-01,-0.17024436E-02,-0.11109697E-02,-0.42947344E-03,
     +  0.16004356E-02, 0.10881806E-02, 0.10017959E-02,-0.19934038E-02,
     + -0.98994588E-02, 0.95261103E-02, 0.22333090E-02, 0.58314172E-02,
     +  0.47419369E-02,-0.55231187E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x32        
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)        *x31        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)    *x21*x32        
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff( 18)        *x32    *x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)            *x42    
     3  +coeff( 21)    *x23            
     4  +coeff( 22)*x11*x22            
     5  +coeff( 23)    *x23*x32        
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)    *x21*x32*x42    
     8  +coeff( 26)    *x23*x31*x42    
c
      return
      end
      function x_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5276050E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.74291918E-02, 0.24634120E+00, 0.33618970E-03, 0.19258599E+00,
     + -0.14033425E+00, 0.20213157E-01,-0.62707339E-04, 0.27125912E-01,
     + -0.99104941E-02,-0.15948399E-03,-0.24730801E-02,-0.15809933E-02,
     +  0.17825816E-02, 0.17820795E-02,-0.74857948E-02,-0.97401493E-03,
     + -0.36974230E-02, 0.70095289E-03,-0.30906880E-02,-0.34129314E-03,
     + -0.11450694E-02, 0.60921692E-03,-0.11409811E-02,-0.71300799E-03,
     + -0.14599601E-02, 0.65559754E-04, 0.12917044E-02,-0.49345044E-03,
     +  0.17401675E-02, 0.78915857E-03,-0.23999659E-02, 0.72018662E-03,
     + -0.20088824E-02, 0.12580790E-03, 0.57178666E-02, 0.12157867E-03,
     +  0.45529706E-03,-0.46297754E-02,-0.96704735E-03,-0.24237682E-03,
     +  0.29334790E-03, 0.66688779E-03,-0.58770570E-03,-0.22688118E-03,
     +  0.25149409E-03,-0.12048014E-02,-0.43748547E-04,-0.80842059E-03,
     + -0.52931067E-03, 0.15779086E-04,-0.37137957E-03,-0.59822210E-04,
     +  0.11534806E-02, 0.67508721E-04,-0.73910994E-03, 0.13640863E-02,
     +  0.21582270E-03, 0.87961578E-03,-0.87604352E-03,-0.32502892E-04,
     +  0.59051043E-03, 0.15440468E-02, 0.29022447E-03,-0.13623354E-02,
     +  0.10817289E-03,-0.23817943E-03,-0.29158959E-03,-0.49479953E-04,
     + -0.38708461E-03,-0.67804212E-04,-0.11972132E-03,-0.69238944E-03,
     +  0.35499275E-03,-0.26485864E-02,-0.17752328E-02, 0.15029732E-03,
     + -0.14119121E-03, 0.22831179E-03,-0.14309406E-03, 0.27841052E-05,
     + -0.59005711E-03,-0.11297813E-03, 0.33008637E-04, 0.45341521E-03,
     +  0.93806174E-03, 0.71785651E-03, 0.41832704E-04, 0.98980789E-03,
     + -0.49343356E-03, 0.43586144E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff(  9)                *x52
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)    *x24            
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 18)                *x53
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x22        *x52
     3  +coeff( 21)    *x24        *x51
     4  +coeff( 22)    *x21*x31        
     5  +coeff( 23)            *x42    
     6  +coeff( 24)        *x32    *x51
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)        *x33    *x51
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 27)    *x23*x32        
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)*x12*x22            
     3  +coeff( 30)    *x22*x31*x41*x51
     4  +coeff( 31)    *x21*x32*x42    
     5  +coeff( 32)    *x23*x32*x41    
     6  +coeff( 33)*x12*x21*x32        
     7  +coeff( 34)        *x31        
     8  +coeff( 35)    *x23            
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 36)*x11    *x31        
     1  +coeff( 37)*x11            *x51
     2  +coeff( 38)*x11*x22            
     3  +coeff( 39)    *x23    *x41    
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)    *x21*x33        
     6  +coeff( 42)    *x21*x32*x41    
     7  +coeff( 43)    *x21*x31*x42    
     8  +coeff( 44)    *x21        *x53
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 45)        *x32    *x52
     1  +coeff( 46)    *x23    *x42    
     2  +coeff( 47)    *x23*x31    *x51
     3  +coeff( 48)    *x23        *x52
     4  +coeff( 49)*x11*x24            
     5  +coeff( 50)        *x33    *x52
     6  +coeff( 51)*x11*x21*x32    *x51
     7  +coeff( 52)        *x33*x42*x51
     8  +coeff( 53)        *x31*x41*x51
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 54)    *x22*x32        
     1  +coeff( 55)*x11    *x31*x41    
     2  +coeff( 56)    *x22    *x42    
     3  +coeff( 57)    *x22*x31    *x51
     4  +coeff( 58)*x11            *x52
     5  +coeff( 59)    *x21    *x42*x51
     6  +coeff( 60)*x12    *x31        
     7  +coeff( 61)*x11*x22*x31        
     8  +coeff( 62)*x11*x22    *x41    
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 63)        *x31*x41*x52
     1  +coeff( 64)*x11*x21*x31*x41    
     2  +coeff( 65)    *x22    *x43    
     3  +coeff( 66)*x11    *x32    *x51
     4  +coeff( 67)    *x22        *x53
     5  +coeff( 68)    *x21*x33    *x51
     6  +coeff( 69)    *x24        *x52
     7  +coeff( 70)        *x32*x41*x52
     8  +coeff( 71)            *x42*x53
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 72)    *x23*x33        
     1  +coeff( 73)    *x23*x33    *x51
     2  +coeff( 74)*x11*x24*x31*x41    
     3  +coeff( 75)*x11*x22*x32*x42    
     4  +coeff( 76)        *x31    *x51
     5  +coeff( 77)            *x41*x51
     6  +coeff( 78)    *x22    *x41    
     7  +coeff( 79)    *x21*x31    *x51
     8  +coeff( 80)        *x33        
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 81)            *x42*x51
     1  +coeff( 82)        *x31    *x52
     2  +coeff( 83)            *x41*x52
     3  +coeff( 84)*x11*x21        *x51
     4  +coeff( 85)    *x21    *x43    
     5  +coeff( 86)    *x21*x32    *x51
     6  +coeff( 87)    *x21*x31    *x52
     7  +coeff( 88)*x11*x22        *x51
     8  +coeff( 89)            *x42*x52
      x_sp_cq3m   =x_sp_cq3m   
     9  +coeff( 90)        *x31    *x53
c
      return
      end
      function t_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 25)
      data ncoeff/ 24/
      data avdat/ -0.1313364E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73739380E-03,-0.38552616E-01, 0.10628402E-03, 0.67402683E-01,
     + -0.57477537E-02,-0.53188200E-02, 0.29684673E-02, 0.23047339E-03,
     + -0.96901122E-03, 0.49376856E-02, 0.52453176E-03,-0.15497782E-02,
     +  0.19436186E-03,-0.65546471E-03,-0.42133810E-03, 0.12754189E-03,
     + -0.12150341E-03,-0.31239673E-03,-0.15245473E-04, 0.41830100E-03,
     + -0.32409391E-03,-0.37400669E-03,-0.27597349E-03,-0.23704558E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x32        
      t_sp_cq3m   =t_sp_cq3m   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)                *x53
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x23            
     5  +coeff( 14)*x11*x23            
     6  +coeff( 15)*x11*x23        *x51
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)    *x21    *x42    
      t_sp_cq3m   =t_sp_cq3m   
     9  +coeff( 18)        *x32    *x51
     1  +coeff( 19)        *x31*x41*x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x22        *x52
     4  +coeff( 22)*x11*x22        *x51
     5  +coeff( 23)*x11    *x31*x41*x51
     6  +coeff( 24)    *x23        *x52
c
      return
      end
      function y_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/  0.8192663E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17258970E-01,-0.16505122E+00, 0.16188601E+00,-0.65800059E-02,
     +  0.84502306E-02, 0.46011287E-03, 0.82980993E-03, 0.92268351E-03,
     +  0.28535157E-01, 0.42435959E-01,-0.48225936E-02, 0.60059316E-03,
     + -0.36295719E-01,-0.45154061E-01, 0.38080816E-02, 0.24414891E-01,
     + -0.97113766E-03,-0.74443654E-02,-0.12567249E-02,-0.10011457E-01,
     + -0.23215771E-01, 0.33960256E-03,-0.23477513E-02,-0.33788143E-04,
     + -0.31094321E-02,-0.32955238E-02,-0.49505439E-02, 0.27545020E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cq3m   =y_sp_cq3m   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)    *x22            
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)    *x21*x31    *x51
      y_sp_cq3m   =y_sp_cq3m   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x22    *x41    
     4  +coeff( 22)        *x35        
     5  +coeff( 23)            *x45    
     6  +coeff( 24)        *x33    *x52
     7  +coeff( 25)        *x31    *x52
     8  +coeff( 26)    *x22*x31    *x51
      y_sp_cq3m   =y_sp_cq3m   
     9  +coeff( 27)        *x32*x41    
     1  +coeff( 28)    *x23    *x41    
c
      return
      end
      function p_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.9917168E-04/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16994345E-02,-0.97293872E-04, 0.12344286E-01,-0.21521583E-01,
     +  0.11514194E-02,-0.98020309E-05,-0.35655688E-03,-0.11608900E-02,
     + -0.66384237E-04, 0.91781104E-02,-0.53022999E-04,-0.98533863E-04,
     + -0.60662087E-02,-0.88748644E-03, 0.27481893E-02, 0.89287008E-04,
     +  0.10584109E-02,-0.22534547E-03, 0.90431131E-04, 0.54550666E-03,
     +  0.58139849E-03, 0.77922683E-03, 0.86731964E-03, 0.18985884E-03,
     + -0.35550148E-03,-0.10861146E-03,-0.38367559E-02, 0.66359571E-04,
     + -0.44324668E-03,-0.27070288E-03,-0.18751006E-03, 0.32122040E-03,
     +  0.54166601E-04, 0.40543229E-04,-0.10936907E-03,-0.28940261E-03,
     +  0.37544763E-04,-0.12852738E-02,-0.25330836E-03,-0.21122221E-03,
     +  0.33490680E-03, 0.23615458E-03, 0.50812925E-03, 0.17844184E-04,
     + -0.10013715E-02,-0.43479195E-04, 0.10013642E-02,-0.15409663E-03,
     +  0.11648051E-03, 0.31132123E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)        *x31    *x52
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 18)                *x52
     1  +coeff( 19)        *x33        
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)*x11    *x31        
     6  +coeff( 24)        *x32*x41    
     7  +coeff( 25)    *x22*x31    *x51
     8  +coeff( 26)    *x21        *x51
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)    *x21    *x41*x51
     3  +coeff( 30)    *x21*x31    *x52
     4  +coeff( 31)        *x31    *x53
     5  +coeff( 32)        *x31*x42    
     6  +coeff( 33)    *x21        *x52
     7  +coeff( 34)                *x53
     8  +coeff( 35)    *x23*x31        
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 36)    *x21*x31*x42    
     1  +coeff( 37)        *x33    *x51
     2  +coeff( 38)    *x22    *x41*x51
     3  +coeff( 39)        *x32*x41*x51
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)    *x23*x31*x41    
     6  +coeff( 42)*x11*x23*x31        
     7  +coeff( 43)*x11*x21    *x43    
     8  +coeff( 44)*x11*x22            
      p_sp_cq3m   =p_sp_cq3m   
     9  +coeff( 45)*x11*x21*x31        
     1  +coeff( 46)*x11*x23            
     2  +coeff( 47)*x11*x21    *x41*x51
     3  +coeff( 48)    *x23*x32        
     4  +coeff( 49)    *x23*x31    *x51
     5  +coeff( 50)*x11*x23    *x41    
c
      return
      end
      function l_sp_cq3m   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.4115843E-02/
      data xmin/
     1 -0.11959E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.49964E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49943E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.65204906E-02,-0.48415804E+00, 0.15537156E-03,-0.31637784E-01,
     +  0.19225679E+00,-0.55372443E-01, 0.70479559E-02,-0.30104467E-02,
     +  0.44369358E-02, 0.10953059E-02,-0.11640308E-03,-0.13003998E-02,
     +  0.29800849E-01,-0.21186830E-01,-0.17646606E-02, 0.10976662E-02,
     + -0.11715513E-02,-0.38601438E-03, 0.22874693E-02, 0.91860478E-03,
     + -0.21451542E-02, 0.98191365E-03, 0.63268712E-03, 0.12831580E-01,
     + -0.71252775E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_cq3m   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x32        
      l_sp_cq3m   =l_sp_cq3m   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)        *x31        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)                *x52
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)        *x31*x41*x51
     8  +coeff( 17)    *x21*x31        
      l_sp_cq3m   =l_sp_cq3m   
     9  +coeff( 18)        *x31    *x51
     1  +coeff( 19)    *x21*x32        
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)            *x42    
     4  +coeff( 22)    *x23*x32        
     5  +coeff( 23)        *x32    *x53
     6  +coeff( 24)*x12*x21*x31*x43    
     7  +coeff( 25)    *x23*x32*x43    
c
      return
      end
      function x_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1436917E-01/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11109993E-01, 0.25605854E+00, 0.46934270E-04, 0.82582649E-03,
     +  0.32360753E+00,-0.17808634E+00, 0.34468137E-01,-0.20908363E-01,
     +  0.25928705E-02,-0.38658930E-02, 0.72892900E-02,-0.41654445E-02,
     +  0.48347847E-02, 0.25449437E-02, 0.13244073E-02, 0.28362705E-02,
     + -0.68625133E-02, 0.16450228E-01,-0.22934210E-01, 0.14353174E-03,
     + -0.17405875E-03,-0.29147787E-02,-0.12840595E-02, 0.21347108E-03,
     + -0.51612775E-02,-0.14241824E-02, 0.25362228E-02,-0.50610225E-02,
     +  0.33082019E-02,-0.22937316E-02, 0.29741541E-01,-0.29055043E-02,
     + -0.31288741E-02,-0.26092408E-02,-0.31934825E-02,-0.88905671E-03,
     + -0.34271166E-03, 0.60138712E-03,-0.15619940E-03,-0.72806506E-04,
     + -0.28750877E-03, 0.90251124E-03, 0.96896978E-03,-0.10321315E-02,
     + -0.20028327E-02, 0.53129165E-03,-0.15137564E-04,-0.65011547E-04,
     + -0.16414974E-03,-0.28158029E-03,-0.11031688E-02,-0.23257913E-03,
     +  0.12604546E-03, 0.80259371E-03,-0.19815203E-03, 0.41406988E-02,
     +  0.57915695E-04,-0.10742890E-02, 0.29948039E-03, 0.25162851E-02,
     +  0.34787107E-03,-0.49191335E-03,-0.87313273E-03, 0.77524933E-03,
     +  0.33398057E-03,-0.10720262E-02,-0.32663860E-02, 0.11312770E-03,
     + -0.22590798E-02, 0.21044776E-03, 0.50410512E-03,-0.48654381E-03,
     +  0.15527906E-02,-0.10606380E-01, 0.57767046E-03, 0.10480012E-02,
     + -0.52745559E-03, 0.55378349E-03,-0.61012118E-03,-0.26243905E-03,
     +  0.21783527E-03, 0.14944744E-03, 0.79672120E-03,-0.47457943E-03,
     +  0.18418397E-03,-0.77563664E-03, 0.87363593E-03,-0.26602234E-03,
     + -0.12940544E-02, 0.84801577E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x21*x31        
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21    *x42    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 18)    *x22            
     1  +coeff( 19)    *x24            
     2  +coeff( 20)        *x32        
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22*x31*x41    
     5  +coeff( 23)            *x42    
     6  +coeff( 24)        *x31    *x51
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)        *x32    *x51
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 27)        *x31*x41*x51
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)    *x22    *x42    
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)*x11*x23            
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)    *x24        *x51
     7  +coeff( 34)    *x23        *x52
     8  +coeff( 35)    *x21*x32*x42    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 36)*x12    *x31*x41    
     1  +coeff( 37)    *x21*x33*x42    
     2  +coeff( 38)*x12*x23    *x41    
     3  +coeff( 39)            *x41*x51
     4  +coeff( 40)    *x22*x31        
     5  +coeff( 41)*x11        *x41    
     6  +coeff( 42)    *x21*x32    *x51
     7  +coeff( 43)        *x32    *x52
     8  +coeff( 44)*x11*x24            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 45)    *x22*x31*x41*x51
     1  +coeff( 46)*x12*x21*x31        
     2  +coeff( 47)    *x21*x33    *x51
     3  +coeff( 48)        *x33    *x52
     4  +coeff( 49)            *x43*x53
     5  +coeff( 50)    *x21*x31    *x51
     6  +coeff( 51)            *x42*x51
     7  +coeff( 52)        *x31    *x52
     8  +coeff( 53)    *x22    *x41*x51
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 54)    *x21*x32*x41    
     1  +coeff( 55)    *x21        *x53
     2  +coeff( 56)*x11*x22        *x51
     3  +coeff( 57)        *x33    *x51
     4  +coeff( 58)        *x31*x41*x52
     5  +coeff( 59)            *x42*x52
     6  +coeff( 60)    *x22    *x42*x51
     7  +coeff( 61)*x11            *x53
     8  +coeff( 62)    *x22        *x53
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 63)*x11*x22*x32        
     1  +coeff( 64)    *x23    *x43    
     2  +coeff( 65)    *x22    *x41    
     3  +coeff( 66)*x11            *x51
     4  +coeff( 67)*x12                
     5  +coeff( 68)            *x41*x52
     6  +coeff( 69)    *x22*x31*x41    
     7  +coeff( 70)    *x21    *x41*x52
     8  +coeff( 71)        *x31    *x53
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 72)            *x41*x53
     1  +coeff( 73)*x11*x21        *x52
     2  +coeff( 74)*x12*x22            
     3  +coeff( 75)    *x22*x32    *x51
     4  +coeff( 76)*x11        *x42*x51
     5  +coeff( 77)    *x21*x32    *x52
     6  +coeff( 78)*x12    *x32        
     7  +coeff( 79)    *x24        *x52
     8  +coeff( 80)        *x32    *x53
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 81)    *x23*x32    *x51
     1  +coeff( 82)    *x23*x31    *x52
     2  +coeff( 83)*x12*x22        *x51
     3  +coeff( 84)*x11*x24        *x51
     4  +coeff( 85)*x11        *x42*x52
     5  +coeff( 86)*x12*x21    *x42    
     6  +coeff( 87)*x12*x21        *x52
     7  +coeff( 88)*x11*x22        *x53
     8  +coeff( 89)*x11*x24*x31*x41    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 90)*x11    *x32*x43    
c
      return
      end
      function t_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.5915408E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28916742E-02, 0.24830580E-01,-0.13281094E-04, 0.32389790E-03,
     +  0.11763829E+00,-0.38108818E-01, 0.71238615E-02,-0.10592095E-01,
     + -0.14403452E-02, 0.23281589E-03, 0.15045023E-02, 0.37980286E-03,
     + -0.19059366E-02,-0.71351806E-03, 0.58956415E-03,-0.11876902E-02,
     +  0.43050470E-02,-0.13630139E-03, 0.31038493E-03,-0.24731096E-02,
     + -0.11892287E-02, 0.35829277E-03,-0.28912630E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)                *x53
     7  +coeff( 16)*x11*x23            
     8  +coeff( 17)*x11*x21            
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)*x11    *x32        
     4  +coeff( 22)        *x32    *x52
     5  +coeff( 23)*x11*x22        *x51
c
      return
      end
      function y_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.6077809E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10383888E-01,-0.10462663E+00, 0.92847377E-01,-0.27359899E-02,
     +  0.53940164E-02, 0.17479743E-03, 0.63541997E-03, 0.43673278E-03,
     +  0.10080937E-01, 0.29576633E-01,-0.29036163E-02,-0.24952043E-01,
     + -0.27022405E-01, 0.19228695E-02, 0.15727755E-01, 0.20475432E-02,
     + -0.61392756E-02,-0.64982497E-03,-0.42432938E-02,-0.16891273E-01,
     + -0.36782082E-03,-0.37765069E-03,-0.31031466E-02,-0.43032249E-03,
     + -0.31119050E-02,-0.12817305E-02, 0.57683745E-03,-0.97515108E-03,
     +  0.21587631E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22            
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x21    *x41*x51
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)        *x35        
     4  +coeff( 22)        *x32*x41    
     5  +coeff( 23)    *x22*x31    *x51
     6  +coeff( 24)                *x52
     7  +coeff( 25)        *x31*x42    
     8  +coeff( 26)            *x41*x52
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)        *x31    *x53
     2  +coeff( 29)    *x23    *x41    
c
      return
      end
      function p_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2157484E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.54493658E-02,-0.17638362E-02, 0.51485654E-01,-0.60605183E-01,
     +  0.26393028E-02, 0.71243319E-03, 0.46280203E-02,-0.33108963E-03,
     +  0.25744438E-01, 0.49242903E-04,-0.38895448E-03,-0.12545122E-01,
     + -0.93179429E-02, 0.11262241E-01, 0.10570735E-01, 0.33843440E-02,
     + -0.12569571E-02,-0.62792207E-03,-0.18009989E-02, 0.65600802E-03,
     +  0.13953102E-02,-0.27575524E-03,-0.51167689E-03,-0.93144306E-03,
     + -0.90475939E-02, 0.21241312E-02,-0.33981851E-03,-0.38592017E-03,
     +  0.13720535E-03,-0.13654655E-01,-0.10427622E-02,-0.41398173E-03,
     + -0.25009603E-03, 0.29220665E-03, 0.59627509E-05, 0.26199492E-03,
     + -0.22595801E-03,-0.44239813E-03, 0.14283776E-03, 0.83547653E-04,
     + -0.44345944E-02, 0.15588156E-02, 0.32677187E-03, 0.83845528E-03,
     + -0.29749196E-03, 0.34404738E-03,-0.69892936E-04, 0.61822580E-02,
     + -0.11940367E-03,-0.40434269E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x22            
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 18)                *x52
     1  +coeff( 19)*x11    *x31        
     2  +coeff( 20)        *x33        
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)        *x32*x41    
     6  +coeff( 24)        *x31    *x53
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)        *x31*x42    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)*x12*x23    *x43    
     2  +coeff( 29)        *x32    *x51
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x21*x31*x42    
     6  +coeff( 33)        *x33    *x51
     7  +coeff( 34)*x11*x21            
     8  +coeff( 35)    *x21*x31*x41    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 36)    *x21    *x42    
     1  +coeff( 37)            *x43    
     2  +coeff( 38)            *x41*x52
     3  +coeff( 39)                *x53
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x23*x31        
     8  +coeff( 44)*x11*x23    *x41    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 45)        *x31*x42*x53
     1  +coeff( 46)*x11*x23*x31    *x51
     2  +coeff( 47)    *x22        *x51
     3  +coeff( 48)*x12    *x31        
     4  +coeff( 49)    *x22*x31*x41    
     5  +coeff( 50)    *x21*x32*x41    
c
      return
      end
      function l_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.4147748E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.71523525E-02,-0.48243490E+00,-0.30866485E-01, 0.18833755E+00,
     + -0.55622347E-01, 0.61144559E-02,-0.87363068E-02, 0.12758211E-02,
     + -0.39007699E-02, 0.62910854E-02, 0.29169332E-01, 0.19037785E-02,
     + -0.15611746E-01,-0.17290056E-02,-0.19423025E-03, 0.14691585E-02,
     + -0.32506408E-02,-0.11725666E-02,-0.58187527E-03, 0.16048868E-02,
     +  0.15183842E-02,-0.77438861E-03, 0.23044914E-02,-0.15387042E-01,
     +  0.14466957E-01,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)        *x31        
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)            *x42*x51
     7  +coeff( 16)                *x53
     8  +coeff( 17)            *x42    
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff( 18)    *x21*x31        
     1  +coeff( 19)        *x31    *x51
     2  +coeff( 20)    *x21*x32        
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)    *x23*x32        
     6  +coeff( 24)    *x23        *x52
     7  +coeff( 25)*x11*x22        *x52
c
      return
      end
      function x_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3991251E-01/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23454748E-01, 0.36369902E+00,-0.32189535E-03, 0.26533396E-02,
     +  0.83042765E+00,-0.34329548E+00, 0.66942655E-01,-0.65779984E-01,
     +  0.81306195E-03, 0.16177766E-02,-0.10491073E-01,-0.48821350E-02,
     +  0.31273791E-02,-0.13072439E-01, 0.16160473E-01,-0.40287443E-01,
     +  0.76029073E-02,-0.16373228E-02,-0.19002808E-01, 0.48465870E-01,
     +  0.86297449E-02, 0.51289448E-02, 0.35005347E-02, 0.92079549E-03,
     + -0.69527916E-03, 0.78516481E-02,-0.49741506E-02, 0.86309630E-02,
     + -0.88687229E-03,-0.15767843E-01,-0.60430127E-02, 0.11012269E-01,
     + -0.68207891E-02, 0.39710198E-01,-0.29135949E-02,-0.15143335E-01,
     + -0.45929803E-02, 0.45046315E-03,-0.29237913E-02, 0.31430179E-02,
     + -0.21752284E-02,-0.67789771E-03,-0.75617344E-02,-0.14248281E-03,
     + -0.11909435E-01,-0.36525328E-02, 0.63630575E-02,-0.14405443E-02,
     + -0.36502073E-02,-0.43449251E-03,-0.20948256E-03,-0.63319219E-03,
     + -0.64385793E-03,-0.20440742E-01, 0.11885255E-02, 0.57513290E-03,
     +  0.68717964E-04, 0.17582638E-02, 0.11138505E-02, 0.11125242E-01,
     +  0.24478321E-02,-0.24106419E-02, 0.26962282E-02,-0.86083892E-03,
     + -0.21382363E-02, 0.16187652E-02,-0.11136075E-01,-0.27294341E-02,
     +  0.17784000E-02,-0.12862264E-02, 0.17442698E-03,-0.15397640E-02,
     +  0.73344522E-03, 0.38653992E-02,-0.91949869E-02, 0.72574738E-03,
     +  0.44939201E-03, 0.87796040E-02,-0.31046625E-02, 0.29005988E-02,
     +  0.36029413E-03,-0.33967721E-02,-0.15919014E-03,-0.10074320E-02,
     + -0.31685871E-02, 0.48457598E-02,-0.24556289E-02, 0.92009371E-02,
     +  0.97518187E-03,-0.39158734E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x24            
     8  +coeff( 17)                *x53
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)        *x32    *x52
     6  +coeff( 24)        *x31    *x51
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)*x11*x22            
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 27)        *x32    *x51
     1  +coeff( 28)        *x31*x41*x51
     2  +coeff( 29)        *x31    *x52
     3  +coeff( 30)    *x23        *x51
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)    *x22        *x52
     7  +coeff( 34)*x11*x23            
     8  +coeff( 35)    *x21    *x42*x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 36)    *x24        *x51
     1  +coeff( 37)    *x23        *x52
     2  +coeff( 38)    *x22*x31        
     3  +coeff( 39)            *x42*x51
     4  +coeff( 40)    *x21*x32    *x51
     5  +coeff( 41)        *x31*x41*x52
     6  +coeff( 42)            *x42*x52
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)    *x23*x31    *x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 45)*x11*x24            
     1  +coeff( 46)    *x22*x31*x41*x51
     2  +coeff( 47)    *x22    *x42*x51
     3  +coeff( 48)    *x21*x33*x41    
     4  +coeff( 49)    *x21*x32*x42    
     5  +coeff( 50)    *x23*x32*x41    
     6  +coeff( 51)    *x23        *x53
     7  +coeff( 52)        *x33*x41*x53
     8  +coeff( 53)    *x21*x31    *x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 54)*x12                
     1  +coeff( 55)    *x23*x31        
     2  +coeff( 56)    *x22*x31    *x51
     3  +coeff( 57)    *x21*x33        
     4  +coeff( 58)    *x21*x32*x41    
     5  +coeff( 59)*x11*x22    *x41    
     6  +coeff( 60)*x11*x22        *x51
     7  +coeff( 61)        *x31    *x53
     8  +coeff( 62)            *x41*x53
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 63)    *x22        *x53
     1  +coeff( 64)            *x42    
     2  +coeff( 65)*x11            *x51
     3  +coeff( 66)*x11            *x52
     4  +coeff( 67)*x12*x21            
     5  +coeff( 68)    *x21*x31*x41*x51
     6  +coeff( 69)    *x21*x31    *x52
     7  +coeff( 70)    *x21    *x41*x52
     8  +coeff( 71)        *x33    *x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 72)    *x23*x32        
     1  +coeff( 73)*x11*x21*x31*x41    
     2  +coeff( 74)*x11*x21        *x52
     3  +coeff( 75)*x12*x22            
     4  +coeff( 76)*x11        *x43    
     5  +coeff( 77)*x11            *x53
     6  +coeff( 78)*x11*x23        *x51
     7  +coeff( 79)    *x21*x32    *x52
     8  +coeff( 80)    *x21*x31*x41*x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 81)*x12    *x32        
     1  +coeff( 82)*x12        *x42    
     2  +coeff( 83)        *x33*x41*x51
     3  +coeff( 84)    *x24        *x52
     4  +coeff( 85)        *x32    *x53
     5  +coeff( 86)        *x31*x41*x53
     6  +coeff( 87)            *x42*x53
     7  +coeff( 88)*x12*x23            
     8  +coeff( 89)    *x23*x32    *x51
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 90)*x11*x21        *x53
c
      return
      end
      function t_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.5915433E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28916688E-02, 0.24829695E-01,-0.13184985E-04, 0.32360086E-03,
     +  0.11763906E+00,-0.38108028E-01, 0.71233232E-02,-0.10592667E-01,
     + -0.14422776E-02, 0.23317224E-03, 0.15022459E-02, 0.37953301E-03,
     + -0.19066088E-02,-0.71512797E-03, 0.58635714E-03,-0.11877090E-02,
     +  0.43068859E-02,-0.13569075E-03, 0.31053610E-03,-0.24710917E-02,
     + -0.11883281E-02, 0.36000853E-03,-0.28935706E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sp_cfp    =t_sp_cfp    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)    *x21        *x52
     6  +coeff( 15)                *x53
     7  +coeff( 16)*x11*x23            
     8  +coeff( 17)*x11*x21            
      t_sp_cfp    =t_sp_cfp    
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x21    *x42    
     3  +coeff( 21)*x11    *x32        
     4  +coeff( 22)        *x32    *x52
     5  +coeff( 23)*x11*x22        *x51
c
      return
      end
      function y_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/ -0.3219719E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13156720E-01, 0.11787751E+00,-0.16958174E+00, 0.87750880E-02,
     + -0.32876355E-02,-0.39696961E-03,-0.10389967E-02,-0.65317110E-03,
     + -0.44799112E-01,-0.10411357E-01, 0.13819630E-02,-0.29573604E-02,
     +  0.25430266E-02, 0.66736288E-01, 0.15928912E-02,-0.97760733E-03,
     + -0.23678061E-01, 0.14413179E-01, 0.12610476E-01, 0.13826165E-01,
     + -0.52768947E-02, 0.22995404E-02,-0.54203096E-03, 0.57342383E-02,
     +  0.45383018E-02, 0.32047290E-03,-0.40555266E-02,-0.51810043E-02,
     + -0.13001270E-02,-0.78968247E-02,-0.81134140E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cfp    =y_sp_cfp    
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)        *x31*x41*x51
     7  +coeff( 16)            *x42*x51
     8  +coeff( 17)*x11        *x41    
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 18)        *x31    *x52
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)        *x31    *x53
     4  +coeff( 22)*x12                
     5  +coeff( 23)        *x35        
     6  +coeff( 24)        *x32*x41    
     7  +coeff( 25)    *x21*x31    *x51
     8  +coeff( 26)                *x53
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)*x11*x21            
     2  +coeff( 29)    *x21        *x53
     3  +coeff( 30)    *x22*x31*x42*x51
     4  +coeff( 31)    *x23    *x43    
c
      return
      end
      function p_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.2157776E-02/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.54466403E-02,-0.16401999E-02, 0.51500559E-01,-0.60619857E-01,
     +  0.26463962E-02, 0.60437690E-03, 0.45856587E-02,-0.31027466E-03,
     +  0.25654128E-01,-0.29684103E-04,-0.33595675E-03,-0.12551730E-01,
     + -0.93114031E-02, 0.50137779E-02, 0.10445638E-01, 0.33953718E-02,
     + -0.14541141E-02,-0.63989288E-03,-0.17485732E-02, 0.60784351E-03,
     +  0.12917693E-02,-0.26662904E-03,-0.27696809E-03,-0.92773320E-03,
     + -0.89947032E-02, 0.17257580E-02,-0.31114923E-03,-0.57876389E-03,
     +  0.13381855E-03,-0.44860400E-03, 0.12864127E-03,-0.12554487E-02,
     + -0.24356938E-03, 0.47615019E-03, 0.77583914E-04, 0.46717640E-03,
     + -0.42652367E-02, 0.21683151E-03,-0.10487661E-02,-0.98610530E-03,
     +  0.15294147E-02, 0.35557905E-03, 0.77850552E-03,-0.31565651E-03,
     +  0.33590139E-03,-0.19707900E-03,-0.81433871E-04, 0.15031845E-03,
     + -0.61681005E-03, 0.55607670E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_cfp    =p_sp_cfp    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x22            
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 18)                *x52
     1  +coeff( 19)*x11    *x31        
     2  +coeff( 20)        *x33        
     3  +coeff( 21)    *x21*x31    *x51
     4  +coeff( 22)    *x21        *x51
     5  +coeff( 23)        *x32*x41    
     6  +coeff( 24)        *x31    *x53
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)        *x31*x42    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)*x12*x23    *x43    
     2  +coeff( 29)        *x32    *x51
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x11*x21*x31        
     6  +coeff( 33)        *x33    *x51
     7  +coeff( 34)*x11*x21            
     8  +coeff( 35)    *x23            
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 36)    *x21*x31*x41    
     1  +coeff( 37)*x11*x21    *x41    
     2  +coeff( 38)    *x21*x33        
     3  +coeff( 39)    *x23    *x41    
     4  +coeff( 40)    *x21*x32*x41    
     5  +coeff( 41)    *x22    *x43    
     6  +coeff( 42)*x11*x23*x31        
     7  +coeff( 43)*x11*x23    *x41    
     8  +coeff( 44)        *x31*x42*x53
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 45)*x11*x23*x31    *x51
     1  +coeff( 46)    *x21*x32        
     2  +coeff( 47)    *x22        *x51
     3  +coeff( 48)    *x21    *x41*x51
     4  +coeff( 49)    *x21        *x52
     5  +coeff( 50)*x11            *x52
c
      return
      end
      function l_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/ -0.1356569E-01/
      data xmin/
     1 -0.11756E+00,-0.45200E-01,-0.59619E-01,-0.31961E-01,-0.42669E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11915E+00, 0.45626E-01, 0.40973E-01, 0.29997E-01, 0.49790E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19518697E-02,-0.48356077E+00,-0.22044443E-02,-0.35340935E-01,
     +  0.19003467E+00,-0.52822310E-01,-0.30486840E-02,-0.38584523E-01,
     + -0.10049467E-01, 0.80718109E-02,-0.44982694E-02, 0.69872080E-02,
     +  0.32681446E-02, 0.45529054E-04, 0.26010524E-01, 0.38697619E-02,
     +  0.49086902E-02, 0.97609486E-03,-0.11750294E-02, 0.24927191E-02,
     + -0.18753174E-02, 0.15153291E-02,-0.61311631E-03, 0.15796501E-01,
     + -0.10870544E-01,-0.13872015E-02, 0.89127549E-04,-0.31682742E-02,
     + -0.69938842E-02,-0.40986211E-03,-0.92613924E-03, 0.18911157E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      l_sp_cfp    =l_sp_cfp    
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)                *x53
     4  +coeff( 13)        *x31        
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)    *x22        *x51
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 18)        *x31*x41*x51
     1  +coeff( 19)    *x23        *x51
     2  +coeff( 20)    *x21*x31*x41*x51
     3  +coeff( 21)            *x42*x52
     4  +coeff( 22)*x11            *x53
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)        *x31*x41    
     7  +coeff( 25)            *x42    
     8  +coeff( 26)        *x31    *x51
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 27)    *x21*x32        
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)*x11    *x31    *x51
     4  +coeff( 31)    *x22*x32        
     5  +coeff( 32)    *x23*x32        
c
      return
      end
      function x_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 45)
      data ncoeff/ 44/
      data avdat/  0.2088940E-02/
      data xmin/
     1 -0.11901E+00,-0.44663E-01,-0.59223E-01,-0.30909E-01,-0.49804E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11976E+00, 0.46730E-01, 0.40960E-01, 0.29710E-01, 0.49988E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10627676E-01, 0.56822175E+00, 0.36194775E-03,-0.29398878E-04,
     +  0.13358927E+00,-0.25170255E+00, 0.19946054E-01, 0.35761023E-03,
     +  0.27826775E-01,-0.34589400E-02,-0.31766144E-03,-0.97595033E-03,
     + -0.22474190E-02, 0.64616936E-03, 0.11393138E-01,-0.87292138E-02,
     +  0.82248123E-02,-0.92627894E-03,-0.32668919E-02,-0.23485634E-02,
     + -0.15953569E-02,-0.31177569E-02,-0.54216823E-02, 0.72385767E-03,
     + -0.27034255E-02,-0.79093287E-02, 0.97950615E-05,-0.56014903E-03,
     +  0.32726486E-03, 0.10792756E-02, 0.63259056E-03, 0.21235947E-03,
     +  0.12321063E-03,-0.32004486E-02, 0.24725595E-02,-0.13447503E-02,
     + -0.16376661E-02, 0.12726040E-02, 0.44992487E-02,-0.96302749E-04,
     + -0.49995695E-03,-0.71667065E-03,-0.26153564E-02,-0.19797441E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      x_sp_cdex   =x_sp_cdex   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)        *x32        
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)            *x42    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x21*x31*x43    
     5  +coeff( 23)*x12*x21*x31*x41    
     6  +coeff( 24)    *x21*x31        
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)*x11*x22            
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 27)            *x43    
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)*x11*x22*x31        
     5  +coeff( 32)        *x32    *x51
     6  +coeff( 33)                *x53
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)*x11*x21        *x51
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 36)*x11    *x32        
     1  +coeff( 37)*x11        *x42    
     2  +coeff( 38)    *x21*x31*x42    
     3  +coeff( 39)*x11*x22    *x41    
     4  +coeff( 40)        *x33    *x51
     5  +coeff( 41)*x11*x21*x31*x41    
     6  +coeff( 42)*x11*x24            
     7  +coeff( 43)    *x21*x33*x41    
     8  +coeff( 44)*x11*x22*x32        
c
      return
      end
      function t_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 33)
      data ncoeff/ 32/
      data avdat/  0.5338570E+00/
      data xmin/
     1 -0.11901E+00,-0.44663E-01,-0.59223E-01,-0.30909E-01,-0.49804E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11976E+00, 0.46730E-01, 0.40960E-01, 0.29710E-01, 0.49988E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16497379E-02,-0.11206154E+00,-0.54484850E-03, 0.22837935E-01,
     +  0.38889069E-01,-0.84792208E-02, 0.20349722E-02,-0.12001591E-02,
     + -0.26480897E-03, 0.67273766E-03, 0.51212643E-03,-0.15791805E-02,
     +  0.53992369E-02,-0.52718082E-02, 0.12431269E-02, 0.62350533E-03,
     +  0.39740998E-03, 0.33668515E-02,-0.19075408E-02,-0.23909924E-03,
     + -0.31068013E-03, 0.28602837E-04, 0.51619201E-04,-0.21612676E-03,
     + -0.78444410E-03, 0.82043363E-04, 0.99968740E-04, 0.30148751E-03,
     +  0.85135034E-04,-0.23977409E-03, 0.52540249E-03, 0.60438825E-03,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      t_sp_cdex   =t_sp_cdex   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)    *x21*x31*x41    
     2  +coeff( 11)        *x31        
     3  +coeff( 12)        *x32        
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)    *x21        *x52
      t_sp_cdex   =t_sp_cdex   
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)        *x31    *x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)        *x31*x41*x51
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)        *x31    *x52
      t_sp_cdex   =t_sp_cdex   
     9  +coeff( 27)                *x53
     1  +coeff( 28)    *x21*x31*x41*x51
     2  +coeff( 29)    *x22        *x52
     3  +coeff( 30)        *x32    *x52
     4  +coeff( 31)    *x23*x32        
     5  +coeff( 32)    *x22*x31*x41*x51
c
      return
      end
      function y_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 34)
      data ncoeff/ 33/
      data avdat/  0.5636352E-02/
      data xmin/
     1 -0.11901E+00,-0.44663E-01,-0.59223E-01,-0.30909E-01,-0.49804E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11976E+00, 0.46730E-01, 0.40960E-01, 0.29710E-01, 0.49988E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11721020E-01,-0.10880479E+00, 0.12513180E+00,-0.43857554E-02,
     +  0.42893179E-02, 0.22166217E-03, 0.30110264E-03, 0.89893019E-03,
     +  0.21375520E-01, 0.27872415E-01,-0.17918543E-02,-0.16918959E-01,
     + -0.48483260E-01, 0.41685919E-02, 0.67202915E-02,-0.17850354E-02,
     + -0.34929599E-03,-0.19080024E-02,-0.31193220E-02,-0.23513082E-02,
     + -0.19310676E-01,-0.15035187E-01, 0.15687810E-03,-0.15667756E-02,
     +  0.35685897E-03, 0.43531469E-03, 0.33843925E-03,-0.36627976E-02,
     +  0.17352311E-01,-0.15243378E-02, 0.10483290E-01,-0.13883281E-02,
     +  0.48379144E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x34 = x33*x3
      x35 = x34*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x21            
     6  +coeff(  6)        *x32        
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      y_sp_cdex   =y_sp_cdex   
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)    *x22            
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x21    *x42    
      y_sp_cdex   =y_sp_cdex   
     9  +coeff( 18)    *x21*x31    *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)    *x22    *x41    
     5  +coeff( 23)        *x35        
     6  +coeff( 24)        *x31*x44    
     7  +coeff( 25)                *x52
     8  +coeff( 26)    *x21        *x51
      y_sp_cdex   =y_sp_cdex   
     9  +coeff( 27)        *x33        
     1  +coeff( 28)        *x32*x41    
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)*x11*x21*x31        
     5  +coeff( 32)    *x22*x31    *x51
     6  +coeff( 33)    *x23*x31*x42    
c
      return
      end
      function p_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 49)
      data ncoeff/ 48/
      data avdat/  0.2779033E-02/
      data xmin/
     1 -0.11901E+00,-0.44663E-01,-0.59223E-01,-0.30909E-01,-0.49804E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11976E+00, 0.46730E-01, 0.40960E-01, 0.29710E-01, 0.49988E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22316561E-02,-0.29419145E-01, 0.17607536E-01,-0.76241401E-03,
     + -0.35265024E-03,-0.29474264E-02,-0.11867854E-01,-0.70568203E-04,
     +  0.19998956E-03, 0.34400548E-02, 0.68773464E-02, 0.46483363E-03,
     + -0.80893189E-02,-0.22029453E-02,-0.65916823E-03, 0.16467790E-02,
     +  0.75249677E-03, 0.15551146E-03, 0.33720091E-04, 0.11299707E-03,
     + -0.30531391E-03,-0.52356627E-03,-0.10297283E-02, 0.37175475E-02,
     + -0.17183153E-03,-0.24786868E-03, 0.38430083E-02, 0.32566820E-03,
     + -0.33261000E-04,-0.54346408E-04, 0.81224985E-04, 0.80300157E-03,
     + -0.13185979E-03,-0.10020250E-03,-0.73177449E-04,-0.64105139E-03,
     + -0.22042764E-03, 0.64612884E-03,-0.35114141E-03, 0.11573351E-03,
     + -0.37541497E-04, 0.26013457E-04,-0.83727879E-03,-0.46313257E-03,
     +  0.20378656E-02, 0.13010879E-03, 0.44334665E-04, 0.59253311E-04,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31*x41    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)        *x31    *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)*x11    *x31        
     8  +coeff( 17)    *x21            
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 18)    *x22            
     1  +coeff( 19)        *x31    *x52
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x22*x31    *x51
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)        *x33        
     8  +coeff( 26)            *x41*x52
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)    *x21*x31*x41    
     3  +coeff( 30)    *x21    *x42    
     4  +coeff( 31)    *x22        *x51
     5  +coeff( 32)    *x21*x31    *x51
     6  +coeff( 33)    *x21*x31    *x52
     7  +coeff( 34)    *x21    *x41*x52
     8  +coeff( 35)        *x31    *x53
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 36)    *x22*x31*x42    
     1  +coeff( 37)    *x23*x31    *x51
     2  +coeff( 38)    *x23    *x41*x51
     3  +coeff( 39)*x11*x23*x31        
     4  +coeff( 40)        *x32        
     5  +coeff( 41)    *x23            
     6  +coeff( 42)    *x21        *x52
     7  +coeff( 43)*x11*x21*x31        
     8  +coeff( 44)*x11    *x31    *x51
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 45)    *x23    *x41    
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)        *x33    *x51
     3  +coeff( 48)*x11*x23            
c
      return
      end
      function l_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 27)
      data ncoeff/ 26/
      data avdat/ -0.1889313E-01/
      data xmin/
     1 -0.11901E+00,-0.44663E-01,-0.59223E-01,-0.30909E-01,-0.49804E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11976E+00, 0.46730E-01, 0.40960E-01, 0.29710E-01, 0.49988E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.10331889E-01,-0.77322590E+00, 0.40179942E-03,-0.98721527E-01,
     +  0.31917179E+00,-0.53658985E-01,-0.14226877E-01, 0.72340732E-02,
     + -0.13941515E-02, 0.18556184E-02, 0.56499720E-03, 0.30857901E-03,
     + -0.19041725E-02, 0.28310958E-01,-0.29626375E-02,-0.24071135E-01,
     + -0.12294577E-02,-0.11159690E-03, 0.18941102E-02, 0.13565740E-02,
     +  0.11868177E-01,-0.24482215E-03,-0.13640049E-02, 0.15073100E-02,
     +  0.54848962E-02,-0.64786705E-02,
     +      0.      /
      data ientry/0/
c
      if (ientry.ne.0) go to 10
      ientry=1
      do 5 i=1,m
      if (xmin(i).eq.xmax(i)) go to 5
      scale(i)=2./(xmax(i)-xmin(i))
   5  continue
c
  10  continue
c      normalize variables between -1 and +1
      x1 =1.+(x(  1)-xmax(  1))*scale(  1)
      x2 =1.+(x(  2)-xmax(  2))*scale(  2)
      x3 =1.+(x(  3)-xmax(  3))*scale(  3)
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x31 = x3
      x32 = x31*x3
      x33 = x32*x3
      x41 = x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x21*x31*x41    
      l_sp_cdex   =l_sp_cdex   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)                *x52
     2  +coeff( 11)        *x31        
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x22        *x51
      l_sp_cdex   =l_sp_cdex   
     9  +coeff( 18)        *x32    *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x23*x31*x41    
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)    *x21*x31    *x51
     6  +coeff( 24)        *x31*x41*x51
     7  +coeff( 25)    *x21*x33*x41    
     8  +coeff( 26)    *x23*x32*x41    
c
      return
      end
