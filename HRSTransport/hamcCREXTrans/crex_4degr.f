C Forward transfer functions for hrs with septum based on prex2_c_dir.dat
C Uses Paul Brindza's 2-coil septum magnet for PREX II at 4 degrees
c  prex2_4deg_dir.dat   (right side) (+y is towards the beam)
c                     -JJL 4/19/12
c
c
c typical call: answer = function(x,M)
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
c NOMENCLATURE: function name = prefix + _s_ +suffix
c           prefixes:     x means xfinal
c                         t means thetafinal
c                         y means yfinal
c                         p means phifinal
c                         l means path length relative to the central trajectory
c     
c           suffixes:     fp means target to focus
c                         sext means septum exit
c                         q1en means Q1 entrance
c                         q1ex means target to Q1 exit
c                         den  means target to dipole entrance
c                         dmd  means dipole middle
c                         dex  means target to dipole exit
c                         q3en means target to Q3 entrance
c                         q3ex means target to Q3 exit
c
c          _s4_ is for hrs with septum at 4 degrees
c
      function x_s4_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2358617E-01/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37154972E-02, 0.15773648E-01,-0.32738692E-03, 0.46166650E-03,
     +  0.78074908E+00,-0.37668891E-01, 0.51882338E-01,-0.58657780E-01,
     + -0.15506763E-01,-0.39856243E-02, 0.17269088E-01, 0.30271782E-01,
     + -0.11437974E-01,-0.39490648E-01, 0.36367249E-01,-0.47608274E-02,
     + -0.23124317E-01,-0.57290015E-02,-0.16579219E-02,-0.70326370E-02,
     +  0.58147558E-02, 0.20169357E-01, 0.47512953E-02,-0.26532365E-02,
     +  0.36972535E-02, 0.25748855E-02, 0.55560526E-02, 0.46855561E-02,
     +  0.14914843E-01,-0.52419226E-02, 0.36498622E-03, 0.12518354E-02,
     + -0.99167333E-03,-0.25794907E-02,-0.26795676E-03, 0.70245273E-03,
     + -0.28431122E-02,-0.22777736E-02,-0.79905689E-02, 0.10094310E-01,
     +  0.38204042E-02, 0.15284722E-02, 0.52890312E-02,-0.10866139E-03,
     + -0.27811122E-02,-0.75064210E-03, 0.28097471E-02, 0.76128443E-03,
     + -0.66100266E-02, 0.63230055E-02,-0.59584570E-02,-0.48982548E-02,
     +  0.12602021E-01, 0.11255069E-02, 0.11274848E-01,-0.41426495E-02,
     +  0.13439402E-02, 0.55713789E-03, 0.40592169E-02, 0.39704377E-03,
     +  0.12178464E-01, 0.77623860E-02,-0.41692096E-03,-0.78479303E-02,
     + -0.24689147E-02, 0.53886964E-03,-0.25592647E-01, 0.10414559E-02,
     +  0.25390915E-02,-0.24803381E-02,-0.10063264E-02,-0.65524323E-03,
     +  0.23525326E-01,-0.17845955E-01,-0.13811835E-01, 0.52106689E-03,
     +  0.26638045E-02, 0.12423628E-01, 0.32814863E-03, 0.14004315E-02,
     + -0.40565428E-03,-0.42631021E-02,-0.48556039E-03,-0.20999499E-02,
     + -0.12057929E-02,-0.34482493E-02,-0.25866760E-02, 0.40815253E-03,
     + -0.75907679E-03,-0.10348629E-02,
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
      x33 = x32*x3
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
      x_s4_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_s4_fp     =x_s4_fp     
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22            
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)    *x21    *x41*x51
      x_s4_fp     =x_s4_fp     
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)                *x53
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)            *x41*x51
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)    *x23*x31        
      x_s4_fp     =x_s4_fp     
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)    *x23    *x42    
     3  +coeff( 30)*x11*x24    *x41    
     4  +coeff( 31)*x11    *x32*x43    
     5  +coeff( 32)    *x22*x31        
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)        *x33    *x51
      x_s4_fp     =x_s4_fp     
     9  +coeff( 36)*x11*x21            
     1  +coeff( 37)*x11        *x41    
     2  +coeff( 38)            *x41*x54
     3  +coeff( 39)    *x24    *x42    
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)*x11*x23    *x41*x51
     6  +coeff( 42)*x11*x23*x31*x42    
     7  +coeff( 43)    *x22    *x41    
     8  +coeff( 44)        *x31*x42    
      x_s4_fp     =x_s4_fp     
     9  +coeff( 45)            *x43    
     1  +coeff( 46)    *x21*x31    *x51
     2  +coeff( 47)    *x21    *x41*x52
     3  +coeff( 48)    *x21        *x53
     4  +coeff( 49)    *x24    *x41    
     5  +coeff( 50)    *x23*x31*x41    
     6  +coeff( 51)    *x21    *x44    
     7  +coeff( 52)    *x24        *x51
     8  +coeff( 53)    *x23    *x41*x51
      x_s4_fp     =x_s4_fp     
     9  +coeff( 54)    *x22    *x42*x51
     1  +coeff( 55)    *x21    *x43*x51
     2  +coeff( 56)            *x44*x51
     3  +coeff( 57)    *x21*x31*x41*x52
     4  +coeff( 58)*x11*x21    *x41    
     5  +coeff( 59)    *x23*x31*x42    
     6  +coeff( 60)*x11*x21        *x51
     7  +coeff( 61)    *x21    *x44*x51
     8  +coeff( 62)    *x21    *x42*x53
      x_s4_fp     =x_s4_fp     
     9  +coeff( 63)*x11    *x33        
     1  +coeff( 64)    *x24*x31*x42    
     2  +coeff( 65)    *x22*x32*x43    
     3  +coeff( 66)*x11*x22        *x51
     4  +coeff( 67)    *x23    *x43*x51
     5  +coeff( 68)        *x33*x43*x51
     6  +coeff( 69)    *x24        *x53
     7  +coeff( 70)*x11*x24            
     8  +coeff( 71)    *x21*x33*x44    
      x_s4_fp     =x_s4_fp     
     9  +coeff( 72)    *x24*x31*x42*x51
     1  +coeff( 73)    *x24    *x43*x51
     2  +coeff( 74)    *x23    *x44*x51
     3  +coeff( 75)    *x22    *x43*x53
     4  +coeff( 76)    *x23*x33*x43    
     5  +coeff( 77)*x11*x23    *x41*x52
     6  +coeff( 78)*x11*x24*x31*x42*x52
     7  +coeff( 79)        *x31*x41*x51
     8  +coeff( 80)    *x22*x31*x41    
      x_s4_fp     =x_s4_fp     
     9  +coeff( 81)        *x33*x41    
     1  +coeff( 82)    *x21*x31*x42    
     2  +coeff( 83)        *x32*x42    
     3  +coeff( 84)            *x44    
     4  +coeff( 85)    *x21*x31*x41*x51
     5  +coeff( 86)    *x21    *x42*x51
     6  +coeff( 87)            *x43*x51
     7  +coeff( 88)    *x21*x31    *x52
     8  +coeff( 89)            *x42*x52
      x_s4_fp     =x_s4_fp     
     9  +coeff( 90)                *x54
c
      return
      end
      function t_s4_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3687019E-02/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87420939E-03,-0.10274219E-01,-0.12135938E-04, 0.39078106E-03,
     +  0.12714450E+00, 0.71680923E-02,-0.10717158E-01,-0.20630085E-02,
     +  0.39863000E-02,-0.26759952E-02, 0.15655167E-02,-0.48108953E-02,
     + -0.16461181E-02,-0.34133473E-02,-0.93193533E-03, 0.34184931E-02,
     + -0.42820684E-03,-0.20838126E-03, 0.10443778E-02,-0.55523135E-03,
     + -0.64183783E-03, 0.23628962E-02, 0.56352979E-03,-0.19920980E-03,
     +  0.78640075E-03, 0.12127740E-02,-0.23969197E-03,-0.74841996E-03,
     + -0.64281041E-04,-0.36636982E-03, 0.87100023E-04, 0.82307619E-04,
     + -0.19959873E-03, 0.22960540E-03, 0.27737158E-03, 0.54204895E-03,
     +  0.10685308E-02,-0.11082349E-02, 0.20486594E-03,-0.49144524E-03,
     + -0.31919451E-03,-0.10752243E-03, 0.10290339E-03,-0.52443676E-03,
     + -0.33889126E-03, 0.10210630E-03, 0.78391022E-04, 0.14304498E-02,
     +  0.20352798E-03, 0.31151030E-04, 0.13365029E-03, 0.24687932E-03,
     + -0.20999195E-03, 0.11784067E-03, 0.16700395E-02, 0.11721385E-02,
     +  0.78027963E-03,-0.72499533E-03, 0.75457629E-03,-0.29281850E-02,
     +  0.20372945E-03,-0.25808820E-03,-0.24246515E-02,-0.26551745E-04,
     +  0.25532319E-03,-0.56532343E-04,-0.41290701E-03,-0.57785186E-04,
     + -0.30869580E-03,-0.16977327E-03,-0.10859642E-03,-0.10028867E-03,
     +  0.21336726E-04, 0.42400497E-03,-0.52356202E-03,-0.41953928E-04,
     + -0.39737599E-03, 0.24945894E-03,-0.19399337E-03,-0.40665534E-03,
     + -0.59846472E-04, 0.39987073E-04,-0.25155637E-03, 0.60320942E-03,
     +  0.13824779E-03,-0.10322456E-03, 0.33531879E-03,-0.56550308E-03,
     + -0.37253564E-03,-0.26944058E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_s4_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_s4_fp     =t_s4_fp     
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x31        
      t_s4_fp     =t_s4_fp     
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)            *x41*x53
     7  +coeff( 25)*x11*x22    *x41    
     8  +coeff( 26)    *x22    *x41*x52
      t_s4_fp     =t_s4_fp     
     9  +coeff( 27)    *x24    *x41*x51
     1  +coeff( 28)            *x44*x53
     2  +coeff( 29)        *x31    *x51
     3  +coeff( 30)            *x41*x51
     4  +coeff( 31)*x11*x21            
     5  +coeff( 32)    *x22*x31        
     6  +coeff( 33)*x11        *x41    
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x23*x31        
      t_s4_fp     =t_s4_fp     
     9  +coeff( 36)    *x21    *x41*x52
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x21    *x44    
     3  +coeff( 39)*x11*x22        *x51
     4  +coeff( 40)    *x22    *x41    
     5  +coeff( 41)            *x43    
     6  +coeff( 42)    *x21*x31    *x51
     7  +coeff( 43)        *x31*x41*x51
     8  +coeff( 44)            *x42*x51
      t_s4_fp     =t_s4_fp     
     9  +coeff( 45)            *x41*x52
     1  +coeff( 46)*x11        *x42    
     2  +coeff( 47)*x11*x21        *x51
     3  +coeff( 48)    *x22    *x41*x51
     4  +coeff( 49)    *x21    *x42*x51
     5  +coeff( 50)*x11            *x52
     6  +coeff( 51)    *x21        *x53
     7  +coeff( 52)    *x22    *x43    
     8  +coeff( 53)    *x21*x31*x43    
      t_s4_fp     =t_s4_fp     
     9  +coeff( 54)*x11*x21    *x41*x51
     1  +coeff( 55)    *x23    *x41*x51
     2  +coeff( 56)    *x22    *x42*x51
     3  +coeff( 57)    *x21    *x43*x51
     4  +coeff( 58)    *x24    *x42    
     5  +coeff( 59)    *x22    *x42*x52
     6  +coeff( 60)    *x23    *x43*x51
     7  +coeff( 61)*x11*x23        *x52
     8  +coeff( 62)    *x21*x33*x41*x53
      t_s4_fp     =t_s4_fp     
     9  +coeff( 63)    *x24    *x43*x52
     1  +coeff( 64)        *x31*x42    
     2  +coeff( 65)    *x22        *x51
     3  +coeff( 66)        *x32*x42    
     4  +coeff( 67)    *x21    *x43    
     5  +coeff( 68)        *x31*x43    
     6  +coeff( 69)            *x43*x51
     7  +coeff( 70)    *x22        *x52
     8  +coeff( 71)            *x42*x52
      t_s4_fp     =t_s4_fp     
     9  +coeff( 72)                *x54
     1  +coeff( 73)*x12*x21            
     2  +coeff( 74)    *x23*x31*x41    
     3  +coeff( 75)    *x24        *x51
     4  +coeff( 76)*x11    *x32    *x51
     5  +coeff( 77)    *x23        *x52
     6  +coeff( 78)    *x21*x31*x41*x52
     7  +coeff( 79)    *x21        *x54
     8  +coeff( 80)            *x41*x54
      t_s4_fp     =t_s4_fp     
     9  +coeff( 81)*x12*x22            
     1  +coeff( 82)    *x23*x31*x42    
     2  +coeff( 83)    *x21*x33*x42    
     3  +coeff( 84)    *x23    *x43    
     4  +coeff( 85)*x11*x21    *x41*x52
     5  +coeff( 86)        *x33*x41*x52
     6  +coeff( 87)    *x21    *x42*x53
     7  +coeff( 88)            *x42*x54
     8  +coeff( 89)*x11*x24    *x41    
      t_s4_fp     =t_s4_fp     
     9  +coeff( 90)    *x24*x31*x42    
c
      return
      end
      function y_s4_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/  0.1341614E-01/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.56707440E-02, 0.19492542E-02,-0.18871230E+00, 0.19805687E-02,
     + -0.54769773E-01,-0.12499028E-01, 0.60710765E-01,-0.25553803E-02,
     + -0.21492960E-01, 0.79160824E-03, 0.62504925E-01, 0.11632221E-01,
     +  0.17643200E-02,-0.12708101E-01, 0.71799964E-01, 0.81550190E-02,
     + -0.14880341E-01, 0.14410715E-02,-0.34446508E-01, 0.11940046E-01,
     + -0.25858847E-01,-0.59741330E-06, 0.43200798E-03,-0.15629136E-02,
     + -0.28619883E-02, 0.28311743E-02,-0.65574730E-02,-0.52847220E-02,
     + -0.17517271E-01, 0.71101529E-02,-0.71692624E-03,-0.20335946E-01,
     + -0.55997842E-02,-0.16437288E-03,-0.71326445E-04,-0.10800847E-01,
     +  0.61956164E-03, 0.67737810E-02,-0.12951860E-02, 0.84961660E-03,
     +  0.25501046E-02,-0.17428319E-02,-0.24637317E-02,-0.78317849E-02,
     + -0.21126557E-02,-0.18253505E-01, 0.82256142E-02,-0.67775667E-03,
     + -0.12643996E-02, 0.67013747E-03,-0.25147391E-02, 0.17163332E-02,
     + -0.15859805E-02,-0.31433448E-02,-0.25491898E-02,-0.97163022E-02,
     +  0.24536038E-02, 0.77203182E-02, 0.21777493E-03,-0.51710120E-03,
     +  0.16619721E-01,-0.66013178E-02, 0.25476664E-02, 0.52885298E-03,
     + -0.61503274E-03, 0.44151209E-03, 0.70142024E-03, 0.13282893E-02,
     +  0.20839633E-03,-0.22408545E-02, 0.26819936E-03, 0.11389952E-02,
     +  0.39258557E-02, 0.20496456E-02,-0.12848276E-01, 0.19005421E-02,
     +  0.33532865E-02,-0.24259984E-02, 0.33927799E-03,-0.61705811E-02,
     + -0.45758295E-02,-0.18376749E-02, 0.71883514E-02,-0.72765863E-02,
     +  0.24897710E-02,-0.31605964E-02,-0.37143105E-02, 0.17725291E-01,
     + -0.17434841E-02, 0.71472754E-02,-0.12191181E-01, 0.40285019E-02,
     +  0.70025520E-02,-0.60650432E-02,-0.30577057E-02,
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
      y_s4_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_s4_fp     =y_s4_fp     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11                
     2  +coeff( 11)    *x22            
     3  +coeff( 12)    *x21        *x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      y_s4_fp     =y_s4_fp     
     9  +coeff( 18)        *x33*x41    
     1  +coeff( 19)            *x44    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)        *x33*x42    
     5  +coeff( 23)    *x22*x33        
     6  +coeff( 24)        *x35*x41    
     7  +coeff( 25)        *x31*x41    
     8  +coeff( 26)    *x21*x31        
      y_s4_fp     =y_s4_fp     
     9  +coeff( 27)        *x31*x42    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)        *x31*x43    
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x24    *x41    
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)        *x32        
     8  +coeff( 35)        *x33        
      y_s4_fp     =y_s4_fp     
     9  +coeff( 36)    *x21    *x42    
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)    *x22*x31        
     3  +coeff( 39)*x11        *x41    
     4  +coeff( 40)            *x41*x52
     5  +coeff( 41)    *x22        *x51
     6  +coeff( 42)    *x21        *x52
     7  +coeff( 43)                *x53
     8  +coeff( 44)    *x23    *x41    
      y_s4_fp     =y_s4_fp     
     9  +coeff( 45)*x11*x22            
     1  +coeff( 46)            *x45    
     2  +coeff( 47)    *x23    *x42    
     3  +coeff( 48)        *x32*x41    
     4  +coeff( 49)    *x21*x31*x41    
     5  +coeff( 50)        *x31    *x52
     6  +coeff( 51)        *x32*x42    
     7  +coeff( 52)    *x22*x31*x41    
     8  +coeff( 53)    *x23*x31        
      y_s4_fp     =y_s4_fp     
     9  +coeff( 54)    *x22    *x41*x51
     1  +coeff( 55)            *x41*x53
     2  +coeff( 56)        *x31*x44    
     3  +coeff( 57)    *x23        *x51
     4  +coeff( 58)    *x21    *x44    
     5  +coeff( 59)    *x21*x31*x44    
     6  +coeff( 60)*x11            *x53
     7  +coeff( 61)    *x22    *x44    
     8  +coeff( 62)*x11*x23    *x41    
      y_s4_fp     =y_s4_fp     
     9  +coeff( 63)*x11*x21    *x44    
     1  +coeff( 64)        *x31*x41*x51
     2  +coeff( 65)*x11        *x42    
     3  +coeff( 66)*x11*x21*x31        
     4  +coeff( 67)*x11        *x41*x51
     5  +coeff( 68)    *x21    *x41*x52
     6  +coeff( 69)*x12                
     7  +coeff( 70)        *x32*x43    
     8  +coeff( 71)*x11            *x52
      y_s4_fp     =y_s4_fp     
     9  +coeff( 72)    *x22        *x52
     1  +coeff( 73)    *x21*x31*x43    
     2  +coeff( 74)            *x44*x51
     3  +coeff( 75)    *x22    *x43    
     4  +coeff( 76)    *x23*x31*x41    
     5  +coeff( 77)*x11*x21    *x42    
     6  +coeff( 78)    *x24*x31        
     7  +coeff( 79)*x12        *x41    
     8  +coeff( 80)    *x21    *x41*x53
      y_s4_fp     =y_s4_fp     
     9  +coeff( 81)            *x41*x54
     1  +coeff( 82)    *x23        *x52
     2  +coeff( 83)    *x22*x31*x43    
     3  +coeff( 84)        *x32*x42*x52
     4  +coeff( 85)    *x23    *x43    
     5  +coeff( 86)*x11*x22    *x42    
     6  +coeff( 87)*x11*x21    *x42*x51
     7  +coeff( 88)    *x22    *x45    
     8  +coeff( 89)*x11*x21*x32*x42    
      y_s4_fp     =y_s4_fp     
     9  +coeff( 90)        *x34*x42*x52
     1  +coeff( 91)    *x24    *x44    
     2  +coeff( 92)*x11*x21*x32*x42*x51
     3  +coeff( 93)*x11*x24    *x42    
     4  +coeff( 94)*x11*x23    *x44    
     5  +coeff( 95)*x11*x24*x31*x42    
c
      return
      end
      function p_s4_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3995606E-02/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21738787E-02, 0.41830304E-03, 0.13333930E-02,-0.70307948E-01,
     + -0.21881592E-01, 0.32924445E-03, 0.24877395E-01, 0.26719294E-01,
     + -0.55843396E-02, 0.46618795E-02,-0.14387258E-01, 0.32212755E-02,
     +  0.30365473E-01,-0.10367173E-01,-0.10677462E-01,-0.68700254E-04,
     +  0.14032943E-02,-0.59132581E-02, 0.26920235E-02,-0.32666868E-02,
     +  0.33390028E-02,-0.12127994E-01,-0.13573027E-02,-0.96225663E-03,
     + -0.79051638E-03,-0.39295079E-02, 0.29329592E-02, 0.55495771E-02,
     + -0.57977214E-02,-0.22632764E-02,-0.88885268E-02,-0.11022699E-03,
     + -0.50908784E-04,-0.51651255E-03,-0.62842574E-03,-0.82945509E-03,
     + -0.51093353E-02, 0.14792371E-02, 0.28193719E-02, 0.35545533E-02,
     + -0.45937384E-03,-0.24671908E-03,-0.83441148E-03, 0.22219191E-03,
     +  0.22149386E-02,-0.23959328E-02,-0.31302269E-02, 0.17087877E-02,
     + -0.12809202E-02,-0.57078787E-03,-0.11193435E-03,-0.12018988E-02,
     + -0.46683366E-04, 0.19621864E-03,-0.43166272E-03,-0.44401496E-03,
     +  0.12602446E-03, 0.17607099E-03,-0.10842204E-02,-0.21470556E-03,
     + -0.23522950E-02,-0.32149374E-02, 0.59469132E-03, 0.21709672E-03,
     +  0.15072597E-02,-0.69231138E-03,-0.74184238E-03, 0.18452842E-02,
     + -0.14518815E-02, 0.37624256E-02, 0.68496685E-02,-0.12832133E-02,
     + -0.75550220E-03, 0.38167724E-03, 0.31056788E-02,-0.19622145E-02,
     + -0.13231096E-02,-0.14294048E-02,-0.15405403E-02, 0.21692391E-02,
     +  0.10749449E-03, 0.15307445E-03, 0.19920959E-03, 0.10124376E-03,
     +  0.13479315E-03, 0.14357663E-03,-0.21492585E-03,-0.24448865E-03,
     + -0.35035456E-03, 0.13446881E-03,
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
      p_s4_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      p_s4_fp     =p_s4_fp     
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x24            
     7  +coeff( 16)        *x33*x41    
     8  +coeff( 17)    *x21*x31        
      p_s4_fp     =p_s4_fp     
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)        *x31*x42    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)            *x44    
     5  +coeff( 23)        *x31*x41    
     6  +coeff( 24)        *x31    *x51
     7  +coeff( 25)                *x52
     8  +coeff( 26)    *x21    *x42    
      p_s4_fp     =p_s4_fp     
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)        *x31*x43    
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)    *x24    *x41    
     5  +coeff( 32)*x11        *x43    
     6  +coeff( 33)        *x32        
     7  +coeff( 34)        *x32*x41    
     8  +coeff( 35)    *x21    *x41*x51
      p_s4_fp     =p_s4_fp     
     9  +coeff( 36)*x11*x22            
     1  +coeff( 37)    *x23    *x41    
     2  +coeff( 38)    *x22*x31*x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)            *x44*x51
     5  +coeff( 41)*x11        *x41    
     6  +coeff( 42)        *x32*x42    
     7  +coeff( 43)    *x24*x31        
     8  +coeff( 44)    *x23*x31*x41    
      p_s4_fp     =p_s4_fp     
     9  +coeff( 45)    *x21    *x44    
     1  +coeff( 46)        *x31*x44    
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)    *x21*x31*x44    
     4  +coeff( 49)            *x43*x53
     5  +coeff( 50)*x11*x23    *x42    
     6  +coeff( 51)*x11            *x51
     7  +coeff( 52)            *x42*x51
     8  +coeff( 53)    *x21        *x52
      p_s4_fp     =p_s4_fp     
     9  +coeff( 54)        *x31    *x52
     1  +coeff( 55)            *x41*x52
     2  +coeff( 56)                *x53
     3  +coeff( 57)*x12                
     4  +coeff( 58)*x11*x21*x31        
     5  +coeff( 59)    *x23*x31        
     6  +coeff( 60)*x11        *x42    
     7  +coeff( 61)    *x21*x31*x42    
     8  +coeff( 62)    *x21    *x43    
      p_s4_fp     =p_s4_fp     
     9  +coeff( 63)    *x23        *x51
     1  +coeff( 64)*x12        *x41    
     2  +coeff( 65)*x11*x21    *x42    
     3  +coeff( 66)    *x24        *x51
     4  +coeff( 67)    *x23        *x52
     5  +coeff( 68)            *x43*x52
     6  +coeff( 69)*x11*x22    *x42    
     7  +coeff( 70)    *x23*x31*x42    
     8  +coeff( 71)    *x23    *x43    
      p_s4_fp     =p_s4_fp     
     9  +coeff( 72)        *x32*x44    
     1  +coeff( 73)    *x23    *x42*x51
     2  +coeff( 74)*x12    *x31*x42    
     3  +coeff( 75)*x11*x24    *x42    
     4  +coeff( 76)    *x22    *x43*x53
     5  +coeff( 77)*x12*x21    *x41*x53
     6  +coeff( 78)        *x32*x43*x54
     7  +coeff( 79)*x12*x24*x34        
     8  +coeff( 80)*x12*x24*x32    *x52
      p_s4_fp     =p_s4_fp     
     9  +coeff( 81)    *x21*x31    *x51
     1  +coeff( 82)    *x21*x33        
     2  +coeff( 83)    *x21*x32*x41    
     3  +coeff( 84)    *x22*x31    *x51
     4  +coeff( 85)    *x21*x31*x41*x51
     5  +coeff( 86)        *x32*x41*x51
     6  +coeff( 87)        *x32    *x52
     7  +coeff( 88)            *x42*x52
     8  +coeff( 89)*x11*x22    *x41    
      p_s4_fp     =p_s4_fp     
     9  +coeff( 90)*x11*x21*x31*x41    
c
      return
      end
      function l_s4_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.1810104E-01/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17167531E-01,-0.23128204E+00,-0.30718939E-01, 0.84986668E-02,
     + -0.26304059E-01, 0.60626134E-01,-0.11366935E-02,-0.40177032E-01,
     + -0.24216248E-01, 0.84609678E-02,-0.14059028E-01,-0.78985570E-02,
     +  0.37122238E-01, 0.63410630E-02,-0.52903378E-02,-0.46469934E-01,
     + -0.26378650E-02, 0.50705704E-02,-0.86014403E-03, 0.39994842E-02,
     +  0.44516167E-02, 0.35669911E-02,-0.39725979E-02,-0.61216732E-02,
     + -0.85372571E-02,-0.17719749E-01, 0.47492018E-03,-0.19141027E-02,
     +  0.50720428E-02,
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
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_s4_fp     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      l_s4_fp     =l_s4_fp     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)            *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)                *x53
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)            *x41    
      l_s4_fp     =l_s4_fp     
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x22    *x41    
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)    *x23*x31        
     7  +coeff( 25)*x11*x22    *x41    
     8  +coeff( 26)    *x23    *x42    
      l_s4_fp     =l_s4_fp     
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)            *x42*x51
     2  +coeff( 29)    *x21    *x43    
c
      return
      end
      function x_s4_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1039176E-01/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.63384214E-03, 0.52294407E-01,-0.32991517E-03,-0.40235865E-03,
     +  0.32579792E+00,-0.27047092E-01, 0.26606919E-01,-0.20427575E-01,
     + -0.82284845E-02,-0.16180557E-02, 0.12072363E-01, 0.16052833E-01,
     + -0.32701278E-02,-0.60876911E-02,-0.20934008E-01, 0.22583757E-01,
     + -0.10498925E-01,-0.16240937E-02, 0.39879270E-02,-0.11864039E-02,
     + -0.12452110E-02, 0.37569946E-02,-0.45234221E-02, 0.89549337E-03,
     +  0.23141028E-02, 0.15538823E-02, 0.13652147E-01, 0.83581044E-03,
     +  0.25738759E-02, 0.12194772E-02,-0.18303404E-02, 0.86600138E-02,
     + -0.38101750E-02,-0.41764957E-04, 0.41818065E-02,-0.21260398E-01,
     + -0.19266647E-02,-0.40271613E-03,-0.65639557E-03,-0.46200482E-02,
     +  0.33900281E-03,-0.29873699E-02, 0.36740450E-02, 0.28826282E-02,
     +  0.70436590E-03, 0.57607698E-02,-0.79615007E-03, 0.21071799E-03,
     + -0.55061001E-03,-0.21579806E-01, 0.13718233E-01, 0.23155459E-02,
     + -0.22720178E-02, 0.12693086E-02, 0.86612366E-02,-0.11751410E-02,
     +  0.87353336E-02,-0.37122918E-02,-0.15460063E-02, 0.28890339E-02,
     +  0.21798916E-02, 0.44612079E-02,-0.64639578E-04,-0.10276624E-02,
     +  0.48671197E-02,-0.37814935E-02, 0.27117049E-03,-0.15332982E-02,
     + -0.23293165E-02,-0.41962642E-03,-0.18100580E-03,-0.51886257E-03,
     +  0.38904374E-03, 0.65519073E-03, 0.64281337E-02,-0.34455347E-03,
     +  0.80431709E-02,-0.19366357E-02, 0.14437888E-02, 0.69126376E-03,
     + -0.48906035E-02,-0.81841201E-02,-0.21848609E-02, 0.36990712E-02,
     + -0.32346894E-02, 0.92831213E-03, 0.50169607E-02,-0.14505964E-02,
     +  0.61512948E-03,-0.11766677E-02,
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
      x33 = x32*x3
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
      x_s4_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21    *x41*x51
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)*x11*x22            
     2  +coeff( 20)        *x31*x41    
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)    *x22    *x41    
     5  +coeff( 23)    *x21*x31*x41    
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)                *x53
     8  +coeff( 26)    *x23*x31        
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)    *x22    *x41*x51
     4  +coeff( 31)*x11        *x41    
     5  +coeff( 32)    *x23    *x42    
     6  +coeff( 33)    *x21    *x44    
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)*x11*x22    *x41    
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 36)    *x24    *x43*x52
     1  +coeff( 37)            *x43    
     2  +coeff( 38)    *x21*x31    *x51
     3  +coeff( 39)            *x42*x51
     4  +coeff( 40)    *x21    *x43    
     5  +coeff( 41)*x11*x21            
     6  +coeff( 42)    *x24    *x41    
     7  +coeff( 43)    *x23*x31*x41    
     8  +coeff( 44)    *x22    *x42*x51
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 45)        *x31*x43*x51
     1  +coeff( 46)    *x23*x31*x42    
     2  +coeff( 47)    *x21*x31*x44    
     3  +coeff( 48)*x11*x21        *x51
     4  +coeff( 49)*x11*x21    *x41*x51
     5  +coeff( 50)    *x23    *x43*x51
     6  +coeff( 51)    *x22    *x43*x52
     7  +coeff( 52)    *x23    *x41*x53
     8  +coeff( 53)*x11*x24            
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 54)*x11*x21*x32*x41    
     1  +coeff( 55)    *x23    *x43*x52
     2  +coeff( 56)    *x21*x32*x42*x53
     3  +coeff( 57)    *x21    *x44*x53
     4  +coeff( 58)    *x24*x33*x42    
     5  +coeff( 59)    *x24*x32*x43    
     6  +coeff( 60)*x11*x23    *x41*x51
     7  +coeff( 61)*x11*x23*x31*x42    
     8  +coeff( 62)    *x21*x33*x44*x52
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 63)        *x31    *x51
     1  +coeff( 64)            *x41*x52
     2  +coeff( 65)    *x22*x31*x41    
     3  +coeff( 66)    *x21*x31*x42    
     4  +coeff( 67)        *x31*x43    
     5  +coeff( 68)            *x44    
     6  +coeff( 69)            *x43*x51
     7  +coeff( 70)                *x54
     8  +coeff( 71)*x11    *x31        
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 72)    *x24*x31        
     1  +coeff( 73)    *x22*x33        
     2  +coeff( 74)    *x21*x32*x42    
     3  +coeff( 75)    *x23    *x41*x51
     4  +coeff( 76)    *x21*x31*x42*x51
     5  +coeff( 77)    *x21    *x43*x51
     6  +coeff( 78)            *x44*x51
     7  +coeff( 79)    *x22    *x41*x52
     8  +coeff( 80)    *x21*x31*x41*x52
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 81)    *x24*x31*x41    
     1  +coeff( 82)    *x24    *x42    
     2  +coeff( 83)    *x22*x32*x42    
     3  +coeff( 84)    *x23    *x43    
     4  +coeff( 85)    *x22*x31*x43    
     5  +coeff( 86)    *x21*x32*x43    
     6  +coeff( 87)    *x22    *x43*x51
     7  +coeff( 88)    *x21*x31*x43*x51
     8  +coeff( 89)    *x23*x31    *x52
      x_s4_q3ex   =x_s4_q3ex   
     9  +coeff( 90)    *x22*x31*x41*x52
c
      return
      end
      function t_s4_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3687125E-02/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.87402645E-03,-0.10274568E-01,-0.12168965E-04, 0.39023106E-03,
     +  0.12714262E+00, 0.71683414E-02,-0.10716294E-01,-0.20629768E-02,
     +  0.39862972E-02,-0.26760502E-02, 0.15655776E-02,-0.48104185E-02,
     + -0.16459684E-02,-0.34115489E-02,-0.92957437E-03, 0.34186931E-02,
     + -0.42825850E-03,-0.20831694E-03, 0.10503867E-02,-0.55533554E-03,
     + -0.64179150E-03, 0.23624464E-02, 0.56411617E-03,-0.19769048E-03,
     +  0.78592834E-03, 0.12067966E-02,-0.23503682E-03,-0.73809747E-03,
     + -0.64133397E-04,-0.36635649E-03, 0.87169370E-04, 0.82352744E-04,
     + -0.19958746E-03, 0.22971435E-03, 0.27750179E-03, 0.53974713E-03,
     +  0.10692453E-02,-0.11085777E-02, 0.20493644E-03,-0.49046357E-03,
     + -0.31807760E-03,-0.10757001E-03, 0.10266571E-03,-0.52111055E-03,
     + -0.33679354E-03, 0.10200653E-03, 0.78223129E-04, 0.14243264E-02,
     +  0.20178234E-03, 0.30624218E-04, 0.13264309E-03, 0.24535129E-03,
     + -0.20937134E-03, 0.11781816E-03, 0.16700648E-02, 0.11688194E-02,
     +  0.77938865E-03,-0.72503171E-03, 0.75102231E-03,-0.29290305E-02,
     +  0.20300003E-03,-0.25811087E-03,-0.24154359E-02,-0.26538179E-04,
     +  0.25558859E-03,-0.56363806E-04,-0.41240623E-03,-0.57952191E-04,
     + -0.30656991E-03,-0.16995637E-03,-0.10753183E-03,-0.10179763E-03,
     +  0.21410895E-04, 0.42403501E-03,-0.52337279E-03,-0.42043441E-04,
     + -0.39754290E-03, 0.24894191E-03,-0.19319767E-03,-0.40093545E-03,
     + -0.59821992E-04, 0.39463572E-04,-0.25125867E-03, 0.60175394E-03,
     +  0.13855506E-03,-0.10294318E-03, 0.33525511E-03,-0.55949978E-03,
     + -0.37201605E-03,-0.26949047E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_s4_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)*x11                
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21*x31        
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 18)        *x31*x41    
     1  +coeff( 19)                *x53
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)            *x41*x53
     7  +coeff( 25)*x11*x22    *x41    
     8  +coeff( 26)    *x22    *x41*x52
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 27)    *x24    *x41*x51
     1  +coeff( 28)            *x44*x53
     2  +coeff( 29)        *x31    *x51
     3  +coeff( 30)            *x41*x51
     4  +coeff( 31)*x11*x21            
     5  +coeff( 32)    *x22*x31        
     6  +coeff( 33)*x11        *x41    
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)    *x23*x31        
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 36)    *x21    *x41*x52
     1  +coeff( 37)    *x23    *x42    
     2  +coeff( 38)    *x21    *x44    
     3  +coeff( 39)*x11*x22        *x51
     4  +coeff( 40)    *x22    *x41    
     5  +coeff( 41)            *x43    
     6  +coeff( 42)    *x21*x31    *x51
     7  +coeff( 43)        *x31*x41*x51
     8  +coeff( 44)            *x42*x51
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 45)            *x41*x52
     1  +coeff( 46)*x11        *x42    
     2  +coeff( 47)*x11*x21        *x51
     3  +coeff( 48)    *x22    *x41*x51
     4  +coeff( 49)    *x21    *x42*x51
     5  +coeff( 50)*x11            *x52
     6  +coeff( 51)    *x21        *x53
     7  +coeff( 52)    *x22    *x43    
     8  +coeff( 53)    *x21*x31*x43    
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 54)*x11*x21    *x41*x51
     1  +coeff( 55)    *x23    *x41*x51
     2  +coeff( 56)    *x22    *x42*x51
     3  +coeff( 57)    *x21    *x43*x51
     4  +coeff( 58)    *x24    *x42    
     5  +coeff( 59)    *x22    *x42*x52
     6  +coeff( 60)    *x23    *x43*x51
     7  +coeff( 61)*x11*x23        *x52
     8  +coeff( 62)    *x21*x33*x41*x53
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 63)    *x24    *x43*x52
     1  +coeff( 64)        *x31*x42    
     2  +coeff( 65)    *x22        *x51
     3  +coeff( 66)        *x32*x42    
     4  +coeff( 67)    *x21    *x43    
     5  +coeff( 68)        *x31*x43    
     6  +coeff( 69)            *x43*x51
     7  +coeff( 70)    *x22        *x52
     8  +coeff( 71)            *x42*x52
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 72)                *x54
     1  +coeff( 73)*x12*x21            
     2  +coeff( 74)    *x23*x31*x41    
     3  +coeff( 75)    *x24        *x51
     4  +coeff( 76)*x11    *x32    *x51
     5  +coeff( 77)    *x23        *x52
     6  +coeff( 78)    *x21*x31*x41*x52
     7  +coeff( 79)    *x21        *x54
     8  +coeff( 80)            *x41*x54
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 81)*x12*x22            
     1  +coeff( 82)    *x23*x31*x42    
     2  +coeff( 83)    *x21*x33*x42    
     3  +coeff( 84)    *x23    *x43    
     4  +coeff( 85)*x11*x21    *x41*x52
     5  +coeff( 86)        *x33*x41*x52
     6  +coeff( 87)    *x21    *x42*x53
     7  +coeff( 88)            *x42*x54
     8  +coeff( 89)*x11*x24    *x41    
      t_s4_q3ex   =t_s4_q3ex   
     9  +coeff( 90)    *x24*x31*x42    
c
      return
      end
      function y_s4_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.8832722E-03/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20242664E-02,-0.27513376E-02, 0.65571219E-01, 0.30601348E-03,
     +  0.23742631E-01, 0.55016908E-02,-0.32399006E-01, 0.29112877E-01,
     + -0.27870838E-03,-0.26411202E-01,-0.56442004E-02, 0.49282913E-02,
     +  0.60343801E-03, 0.58568758E-02, 0.29107076E-03,-0.36431432E-01,
     + -0.52247290E-02,-0.12546183E-02,-0.35598774E-02, 0.66621881E-02,
     + -0.89992676E-02,-0.47259041E-03, 0.14078603E-01,-0.53758551E-04,
     +  0.57256250E-02, 0.12190686E-01, 0.10287588E-04,-0.13613347E-03,
     +  0.44645130E-03, 0.13047529E-02,-0.17788154E-02, 0.87498396E-03,
     +  0.28748680E-02, 0.44432604E-02,-0.29104403E-02, 0.58778789E-03,
     +  0.71147275E-02,-0.25316058E-02,-0.34312536E-02, 0.76555839E-03,
     +  0.82066962E-02, 0.25569182E-02,-0.43632830E-02, 0.11651267E-03,
     + -0.12280566E-02,-0.10348952E-02, 0.46629831E-02, 0.46218775E-04,
     + -0.10275522E-02, 0.54482808E-02,-0.31601009E-03, 0.18200768E-02,
     + -0.37011527E-02, 0.23736649E-02, 0.21644121E-02,-0.59609530E-02,
     +  0.63834297E-04,-0.98965422E-04, 0.11019749E-02,-0.22754583E-02,
     +  0.12174863E-02, 0.36962032E-02, 0.12713009E-02, 0.60397829E-03,
     + -0.22596512E-02,-0.39455670E-03,-0.12239151E-02, 0.10120629E-02,
     + -0.20180396E-02,-0.23845565E-02,-0.37226412E-02, 0.24576671E-02,
     + -0.16399261E-02, 0.10374006E-02,-0.12984261E-03, 0.11561172E-03,
     +  0.13926541E-02, 0.22860916E-03, 0.27418518E-02,-0.17315301E-03,
     + -0.51234744E-03, 0.18102191E-03,-0.89460635E-04, 0.20455013E-02,
     + -0.67782070E-03,-0.14881366E-02,-0.26482739E-03,-0.90883882E-03,
     + -0.29901267E-03,-0.17276853E-03, 0.12519016E-02,-0.80501591E-03,
     + -0.25220609E-02, 0.29276386E-02,-0.69536013E-03,
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
      y_s4_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x41*x51
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)        *x32*x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21    *x41*x51
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)    *x23            
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)        *x33*x41    
     5  +coeff( 23)            *x44    
     6  +coeff( 24)    *x21*x33        
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x24            
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 27)        *x33*x42    
     1  +coeff( 28)    *x22*x33        
     2  +coeff( 29)        *x35*x41    
     3  +coeff( 30)        *x31*x41    
     4  +coeff( 31)    *x21*x31        
     5  +coeff( 32)        *x31    *x51
     6  +coeff( 33)        *x31*x42    
     7  +coeff( 34)    *x21    *x42    
     8  +coeff( 35)    *x22*x31        
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 36)*x11        *x41    
     1  +coeff( 37)        *x31*x43    
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)    *x22    *x41*x51
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)    *x24    *x41    
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)    *x22    *x44    
     8  +coeff( 44)        *x32        
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)                *x53
     2  +coeff( 47)    *x21    *x43    
     3  +coeff( 48)            *x43*x51
     4  +coeff( 49)    *x22*x31*x41    
     5  +coeff( 50)            *x45    
     6  +coeff( 51)*x11*x21        *x51
     7  +coeff( 52)            *x44*x51
     8  +coeff( 53)    *x23    *x42    
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 54)    *x24        *x51
     1  +coeff( 55)        *x31*x44*x51
     2  +coeff( 56)    *x23    *x43    
     3  +coeff( 57)        *x31*x41*x51
     4  +coeff( 58)    *x21*x31    *x51
     5  +coeff( 59)        *x32*x42    
     6  +coeff( 60)    *x22    *x42    
     7  +coeff( 61)    *x23*x31        
     8  +coeff( 62)        *x31*x44    
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 63)    *x23        *x51
     1  +coeff( 64)    *x22        *x52
     2  +coeff( 65)    *x21    *x44    
     3  +coeff( 66)    *x23*x31*x41    
     4  +coeff( 67)    *x22    *x42*x51
     5  +coeff( 68)    *x24*x31        
     6  +coeff( 69)    *x22*x31*x43    
     7  +coeff( 70)    *x21    *x44*x51
     8  +coeff( 71)    *x23*x31*x42    
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 72)*x11*x23    *x41    
     1  +coeff( 73)    *x21*x31*x45    
     2  +coeff( 74)    *x21    *x41*x54
     3  +coeff( 75)*x11*x21    *x44    
     4  +coeff( 76)*x11            *x51
     5  +coeff( 77)    *x21*x31*x42    
     6  +coeff( 78)*x11        *x42    
     7  +coeff( 79)    *x21    *x42*x51
     8  +coeff( 80)*x11*x21*x31        
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 81)    *x22*x31    *x51
     1  +coeff( 82)*x11        *x41*x51
     2  +coeff( 83)*x12                
     3  +coeff( 84)        *x31*x43*x51
     4  +coeff( 85)    *x21*x31*x42*x51
     5  +coeff( 86)            *x43*x52
     6  +coeff( 87)*x11*x21*x31*x41    
     7  +coeff( 88)*x11*x21    *x42    
     8  +coeff( 89)        *x31*x41*x53
      y_s4_q3ex   =y_s4_q3ex   
     9  +coeff( 90)*x12        *x41    
     1  +coeff( 91)    *x22    *x41*x52
     2  +coeff( 92)            *x41*x54
     3  +coeff( 93)    *x21    *x45    
     4  +coeff( 94)            *x45*x51
     5  +coeff( 95)        *x32*x42*x52
c
      return
      end
      function p_s4_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3995553E-02/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21280209E-02, 0.42160146E-03, 0.13373302E-02,-0.70310146E-01,
     + -0.21870857E-01, 0.28539298E-03, 0.24819840E-01, 0.26541643E-01,
     + -0.56538642E-02, 0.46847337E-02,-0.14431493E-01, 0.32550711E-02,
     +  0.30343294E-01,-0.10384324E-01,-0.10652846E-01,-0.27347971E-04,
     +  0.14913734E-02,-0.59055374E-02, 0.26453559E-02,-0.33905189E-02,
     +  0.33723095E-02,-0.12100256E-01,-0.13515224E-02,-0.90279296E-03,
     + -0.89477043E-03,-0.37751736E-02, 0.28308528E-02, 0.57174503E-02,
     + -0.58626141E-02,-0.22943323E-02,-0.88420399E-02, 0.11768008E-03,
     + -0.15992005E-03,-0.23782659E-03,-0.69553405E-03,-0.47076191E-02,
     +  0.13936546E-02, 0.63086831E-04, 0.24803439E-02, 0.33875008E-02,
     + -0.60226623E-03,-0.86749053E-04,-0.11909283E-02,-0.18818327E-03,
     + -0.91477431E-03, 0.21045173E-02,-0.28963205E-02, 0.40985076E-02,
     + -0.15412435E-02,-0.26346129E-03,-0.65211923E-03,-0.12113504E-03,
     + -0.12091085E-02,-0.51885290E-04, 0.21325606E-03, 0.23544501E-03,
     + -0.44297121E-03, 0.77431170E-04, 0.18814929E-03,-0.15539217E-02,
     + -0.29043395E-02, 0.55885682E-03,-0.28691269E-03, 0.21144973E-03,
     +  0.12212682E-02, 0.61358942E-03,-0.21407150E-02,-0.81563205E-03,
     + -0.66936255E-03, 0.15405116E-02, 0.27460564E-03,-0.12451060E-02,
     + -0.10346151E-02,-0.23113337E-03, 0.62958272E-02,-0.12693639E-02,
     + -0.10634362E-02, 0.22410929E-03, 0.31136992E-03,-0.70170808E-03,
     + -0.80993649E-03, 0.13041956E-03, 0.61251951E-04, 0.20177114E-03,
     + -0.13006870E-03, 0.12527961E-03, 0.19303366E-03,-0.34724656E-03,
     +  0.15756195E-03, 0.12493457E-03,
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
      p_s4_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x24            
     7  +coeff( 16)        *x33*x41    
     8  +coeff( 17)    *x21*x31        
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x22*x31        
     2  +coeff( 20)        *x31*x42    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)            *x44    
     5  +coeff( 23)        *x31*x41    
     6  +coeff( 24)        *x31    *x51
     7  +coeff( 25)                *x52
     8  +coeff( 26)    *x21    *x42    
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)        *x31*x43    
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)    *x24    *x41    
     5  +coeff( 32)*x11        *x43    
     6  +coeff( 33)        *x32        
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)*x11*x22            
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 36)    *x23    *x41    
     1  +coeff( 37)    *x22*x31*x41    
     2  +coeff( 38)        *x34*x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)            *x44*x51
     5  +coeff( 41)*x11        *x41    
     6  +coeff( 42)    *x21*x31*x41    
     7  +coeff( 43)    *x23*x31        
     8  +coeff( 44)        *x32*x42    
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 45)    *x24*x31        
     1  +coeff( 46)    *x21    *x44    
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)    *x23*x31*x42    
     4  +coeff( 49)            *x43*x53
     5  +coeff( 50)*x11*x23    *x42    
     6  +coeff( 51)        *x32*x41    
     7  +coeff( 52)*x11            *x51
     8  +coeff( 53)            *x42*x51
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 54)    *x21        *x52
     1  +coeff( 55)        *x31    *x52
     2  +coeff( 56)            *x41*x52
     3  +coeff( 57)                *x53
     4  +coeff( 58)*x12                
     5  +coeff( 59)*x11*x21*x31        
     6  +coeff( 60)    *x21*x31*x42    
     7  +coeff( 61)    *x21    *x43    
     8  +coeff( 62)    *x23        *x51
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 63)        *x31*x42*x51
     1  +coeff( 64)*x12        *x41    
     2  +coeff( 65)*x11*x21    *x42    
     3  +coeff( 66)    *x21*x31*x43    
     4  +coeff( 67)        *x31*x44    
     5  +coeff( 68)    *x24        *x51
     6  +coeff( 69)    *x23        *x52
     7  +coeff( 70)            *x43*x52
     8  +coeff( 71)    *x21*x31    *x53
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 72)    *x21    *x41*x53
     1  +coeff( 73)            *x41*x54
     2  +coeff( 74)*x11    *x32*x42    
     3  +coeff( 75)    *x23    *x43    
     4  +coeff( 76)        *x32*x44    
     5  +coeff( 77)    *x23    *x42*x51
     6  +coeff( 78)*x12            *x52
     7  +coeff( 79)*x12    *x31*x42    
     8  +coeff( 80)    *x23*x34    *x52
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 81)*x12*x24*x32        
     1  +coeff( 82)    *x22*x32        
     2  +coeff( 83)*x11    *x31*x41    
     3  +coeff( 84)    *x21*x32*x41    
     4  +coeff( 85)*x11        *x42    
     5  +coeff( 86)    *x22*x31    *x51
     6  +coeff( 87)        *x32*x41*x51
     7  +coeff( 88)            *x42*x52
     8  +coeff( 89)    *x22*x33        
      p_s4_q3ex   =p_s4_q3ex   
     9  +coeff( 90)*x11*x21*x31*x41    
c
      return
      end
      function l_s4_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.8507318E-02/
      data xmin/
     1 -0.49937E-02,-0.38003E-01,-0.49948E-02,-0.26471E-01,-0.41514E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.75396304E-02,-0.23166338E+00,-0.28702687E-01, 0.84415935E-02,
     + -0.28156254E-01, 0.60285501E-01,-0.24004098E-01, 0.81185596E-02,
     + -0.10262465E-01,-0.52404338E-02,-0.43840423E-01, 0.10815352E-02,
     + -0.20460140E-01,-0.20123261E-02,-0.51880986E-02,-0.63207555E-02,
     + -0.31104086E-02, 0.32299448E-01, 0.38417305E-02, 0.76838071E-02,
     + -0.52561616E-02, 0.66296379E-02,-0.78919977E-02,-0.46289933E-03,
     + -0.71624550E-03, 0.13390484E-02, 0.25114494E-02,-0.55222558E-02,
     + -0.21430517E-02,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_s4_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x21*x31        
      l_s4_q3ex   =l_s4_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x22            
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x23        *x51
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)            *x41    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)    *x21        *x51
     8  +coeff( 17)            *x41*x51
      l_s4_q3ex   =l_s4_q3ex   
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)                *x53
      l_s4_q3ex   =l_s4_q3ex   
     9  +coeff( 27)    *x24            
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)*x11        *x43*x51
c
      return
      end
      function x_s4_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1068199E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.44784476E-02, 0.15006925E+00,-0.22107263E-02, 0.12458183E+00,
     +  0.16406594E-01,-0.45238774E-01, 0.20693509E-01, 0.18605074E-01,
     + -0.78255059E-02,-0.52023605E-02,-0.43478170E-02,-0.23253897E-01,
     +  0.33970587E-01,-0.47448440E-02, 0.62526925E-02,-0.17559148E-02,
     +  0.98623661E-02,-0.68479208E-02, 0.15244351E-02,-0.34814994E-02,
     +  0.31550745E-02,-0.30732693E-02, 0.14711015E-01, 0.70233946E-02,
     + -0.58665010E-03,-0.93703653E-03,-0.14461400E-02,-0.39435350E-02,
     +  0.20200077E-02, 0.29598409E-02,-0.44220922E-03, 0.72905822E-02,
     + -0.14410425E-02,-0.89172460E-02, 0.72526699E-02,-0.12662162E-02,
     +  0.42271432E-02, 0.30949973E-02, 0.31233365E-02,-0.18303922E-02,
     +  0.76245791E-02, 0.20544257E-02,-0.96117361E-02, 0.12999977E-02,
     + -0.34401864E-02,-0.22692401E-01, 0.92696836E-02,-0.35595803E-02,
     +  0.38847302E-02, 0.59596240E-02,-0.81571359E-02,-0.17520860E-02,
     +  0.14333106E-01, 0.23302981E-02,-0.20478664E-02, 0.13768942E-02,
     + -0.10370808E-04, 0.44178241E-02,-0.45957495E-02,-0.24452126E-02,
     + -0.21779658E-02, 0.11148334E-02,-0.39285197E-03, 0.67866064E-03,
     + -0.60862425E-03,-0.22194661E-03,-0.22104189E-02,-0.19880519E-02,
     +  0.73968014E-03,-0.29211917E-02, 0.12727352E-03,-0.16253344E-02,
     + -0.21434990E-02, 0.31462761E-02,-0.34246589E-02, 0.16239756E-02,
     +  0.70465903E-03,-0.45684201E-03, 0.66059030E-03,-0.29075216E-03,
     + -0.48787436E-02, 0.13721394E-02,-0.66485646E-03,-0.53282003E-02,
     +  0.42646988E-02, 0.77354256E-02, 0.36891794E-02,-0.19192852E-02,
     +  0.24728500E-02, 0.22476506E-02,
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
      x_s4_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      x_s4_q3en   =x_s4_q3en   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)    *x22    *x41    
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)        *x31        
     8  +coeff( 26)        *x31*x41    
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 27)    *x24            
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11*x23    *x41*x51
     4  +coeff( 31)    *x21*x31    *x51
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)            *x42*x52
     7  +coeff( 34)    *x24    *x41    
     8  +coeff( 35)    *x23*x31*x41    
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 36)    *x22*x32*x41    
     1  +coeff( 37)    *x23*x31*x42    
     2  +coeff( 38)    *x24    *x41*x51
     3  +coeff( 39)    *x23    *x42*x51
     4  +coeff( 40)    *x24*x32    *x51
     5  +coeff( 41)    *x23*x31*x42*x51
     6  +coeff( 42)    *x22*x32    *x53
     7  +coeff( 43)    *x22*x31*x41*x53
     8  +coeff( 44)            *x44*x53
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 45)*x11*x24            
     1  +coeff( 46)    *x23    *x44*x51
     2  +coeff( 47)    *x22*x33*x41*x53
     3  +coeff( 48)*x11*x24    *x41*x51
     4  +coeff( 49)    *x24*x34*x41*x51
     5  +coeff( 50)*x11*x23    *x42*x51
     6  +coeff( 51)    *x24*x32*x42*x52
     7  +coeff( 52)    *x24*x33*x44    
     8  +coeff( 53)*x11*x24*x31*x43*x52
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 54)    *x22*x31        
     1  +coeff( 55)            *x43    
     2  +coeff( 56)        *x31*x41*x51
     3  +coeff( 57)            *x42*x51
     4  +coeff( 58)    *x22*x31*x41    
     5  +coeff( 59)    *x21*x31*x42    
     6  +coeff( 60)            *x44    
     7  +coeff( 61)    *x21*x31*x41*x51
     8  +coeff( 62)    *x21    *x42*x51
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 63)    *x21    *x41*x52
     1  +coeff( 64)    *x21        *x53
     2  +coeff( 65)                *x54
     3  +coeff( 66)*x11    *x31        
     4  +coeff( 67)    *x24*x31        
     5  +coeff( 68)    *x22*x31*x42    
     6  +coeff( 69)    *x21*x32*x42    
     7  +coeff( 70)    *x21    *x44    
     8  +coeff( 71)*x11            *x51
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 72)        *x33*x41*x51
     1  +coeff( 73)    *x21*x31*x42*x51
     2  +coeff( 74)    *x21    *x43*x51
     3  +coeff( 75)            *x44*x51
     4  +coeff( 76)    *x21    *x42*x52
     5  +coeff( 77)        *x31*x42*x52
     6  +coeff( 78)    *x21        *x54
     7  +coeff( 79)*x11*x21    *x41    
     8  +coeff( 80)*x11    *x31*x41    
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 81)    *x24*x31*x41    
     1  +coeff( 82)    *x23*x32*x41    
     2  +coeff( 83)*x11        *x42    
     3  +coeff( 84)    *x24    *x42    
     4  +coeff( 85)    *x22    *x44    
     5  +coeff( 86)    *x21    *x44*x51
     6  +coeff( 87)            *x44*x52
     7  +coeff( 88)    *x22    *x41*x53
     8  +coeff( 89)    *x21*x31*x41*x53
      x_s4_q3en   =x_s4_q3en   
     9  +coeff( 90)    *x21    *x41*x54
c
      return
      end
      function t_s4_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.9392374E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.19097147E-02,-0.59374485E-01, 0.10134414E-02, 0.20410536E-01,
     +  0.17865680E-02,-0.38874156E-02, 0.14808417E-01,-0.26799331E-02,
     + -0.60088034E-02, 0.19162217E-02,-0.10608340E-02,-0.10004588E-01,
     + -0.38373880E-02, 0.50499784E-02,-0.17274793E-02, 0.15228019E-02,
     +  0.40615810E-03, 0.10556196E-02,-0.34029458E-02, 0.20945591E-02,
     + -0.13918593E-02,-0.22818255E-02, 0.26173114E-02, 0.18479997E-03,
     + -0.14793074E-03,-0.12069282E-02, 0.15215100E-03, 0.23076760E-02,
     + -0.21123949E-02, 0.16553627E-02,-0.10831722E-02,-0.60074643E-04,
     +  0.47804770E-03, 0.18343957E-03,-0.87405229E-03,-0.34158679E-02,
     +  0.64838573E-03,-0.91336129E-04,-0.17424590E-03,-0.16595189E-03,
     +  0.62247860E-03,-0.90420846E-03, 0.19281509E-02, 0.53877890E-03,
     +  0.17884648E-02, 0.62856416E-03, 0.15830786E-03, 0.13078478E-02,
     +  0.30792765E-02,-0.28272858E-02, 0.57667232E-03, 0.68436330E-03,
     + -0.43395418E-03,-0.37224818E-03,-0.48711221E-03,-0.75030304E-03,
     + -0.53445797E-03, 0.11574519E-02,-0.11616414E-02, 0.72724622E-04,
     + -0.14208327E-02,-0.14741274E-02, 0.10441428E-02,-0.13164891E-02,
     +  0.56918262E-03, 0.72829111E-03,-0.35137124E-02,-0.60596084E-02,
     + -0.23121176E-04, 0.68135669E-05, 0.90878282E-04,-0.71071193E-03,
     + -0.66964283E-04, 0.51275193E-03, 0.79618039E-04, 0.59598595E-04,
     +  0.10561554E-03,-0.17844689E-03,-0.14436800E-04, 0.15494980E-03,
     +  0.17354035E-03, 0.57622841E-04, 0.39246396E-03,-0.43180349E-03,
     +  0.24823754E-03,-0.14981272E-03, 0.19115627E-03, 0.67250844E-03,
     + -0.13237864E-03,-0.11225967E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_s4_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_s4_q3en   =t_s4_q3en   
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)            *x42    
     8  +coeff( 17)            *x41*x51
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)        *x31        
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x21    *x41*x51
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 27)    *x24            
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)    *x23*x31*x41    
     3  +coeff( 30)*x11*x23    *x42    
     4  +coeff( 31)*x11*x23    *x41*x51
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)    *x21        *x52
     8  +coeff( 35)*x11*x21    *x41    
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 36)    *x22    *x42    
     1  +coeff( 37)    *x21*x31*x42    
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)            *x42*x52
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)    *x22*x32*x41    
     6  +coeff( 42)*x11*x21    *x42    
     7  +coeff( 43)    *x21    *x44    
     8  +coeff( 44)    *x22*x32    *x51
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)*x11*x24            
     2  +coeff( 47)*x11*x23*x31        
     3  +coeff( 48)*x11*x23    *x41    
     4  +coeff( 49)    *x24    *x42    
     5  +coeff( 50)    *x23    *x43    
     6  +coeff( 51)*x11        *x44    
     7  +coeff( 52)*x11*x22    *x41*x51
     8  +coeff( 53)    *x24    *x41*x51
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 54)*x11        *x43*x51
     1  +coeff( 55)    *x21*x31*x41*x53
     2  +coeff( 56)    *x21    *x42*x53
     3  +coeff( 57)    *x21    *x41*x54
     4  +coeff( 58)    *x24*x31*x42    
     5  +coeff( 59)    *x24*x32    *x51
     6  +coeff( 60)    *x24*x31    *x52
     7  +coeff( 61)    *x22*x31*x42*x52
     8  +coeff( 62)*x11*x23*x31*x42    
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 63)*x11*x22    *x42*x52
     1  +coeff( 64)    *x21*x31*x43*x53
     2  +coeff( 65)*x12    *x33*x42    
     3  +coeff( 66)*x11    *x33*x44    
     4  +coeff( 67)*x11*x22*x31*x42*x52
     5  +coeff( 68)*x11*x22    *x44*x52
     6  +coeff( 69)        *x31    *x51
     7  +coeff( 70)*x11*x21            
     8  +coeff( 71)*x11    *x31        
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 72)    *x22*x31        
     1  +coeff( 73)        *x32*x41    
     2  +coeff( 74)            *x43    
     3  +coeff( 75)                *x53
     4  +coeff( 76)*x11    *x32        
     5  +coeff( 77)*x11    *x31*x41    
     6  +coeff( 78)    *x21*x32*x41    
     7  +coeff( 79)        *x33*x41    
     8  +coeff( 80)        *x32*x42    
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 81)        *x31*x43    
     1  +coeff( 82)*x11    *x31    *x51
     2  +coeff( 83)    *x21*x31*x41*x51
     3  +coeff( 84)    *x21    *x42*x51
     4  +coeff( 85)            *x43*x51
     5  +coeff( 86)    *x21*x31    *x52
     6  +coeff( 87)*x11*x23            
     7  +coeff( 88)    *x24*x31        
     8  +coeff( 89)*x11*x21*x32        
      t_s4_q3en   =t_s4_q3en   
     9  +coeff( 90)    *x23*x32        
c
      return
      end
      function y_s4_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.9196191E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.40974650E-02,-0.31673741E-02, 0.14912929E+00, 0.60042640E-03,
     +  0.48706338E-01,-0.56705721E-01, 0.19012865E-02, 0.39682124E-01,
     + -0.61579474E-03,-0.54554302E-01,-0.10170667E-01, 0.45943460E-02,
     +  0.12007087E-01,-0.68933710E-01,-0.30237171E-02,-0.71549476E-02,
     +  0.12540860E-01,-0.91111157E-02, 0.80722850E-04, 0.20604711E-02,
     +  0.14044907E-01, 0.28794589E-01,-0.22510785E-03, 0.23995381E-01,
     + -0.24562798E-03, 0.78085321E-02,-0.24527137E-03, 0.10475091E-01,
     +  0.87649142E-02,-0.30327106E-02,-0.73089013E-02,-0.57972623E-02,
     +  0.16222858E-02, 0.11831084E-01, 0.18988460E-01, 0.51760823E-02,
     +  0.24636136E-02,-0.26003236E-02, 0.13231694E-02, 0.58553931E-02,
     + -0.57365042E-02, 0.12166932E-02, 0.51790895E-02,-0.26392918E-02,
     +  0.94228042E-02,-0.17879168E-02,-0.71847974E-02,-0.10059370E-01,
     +  0.58852509E-02,-0.91280113E-03,-0.72206130E-04,-0.85500500E-03,
     + -0.25593459E-02, 0.19453852E-02, 0.24083843E-02,-0.82118288E-02,
     + -0.67105484E-02, 0.17045358E-03,-0.51617343E-03,-0.88659021E-04,
     +  0.26069517E-03,-0.44992351E-03, 0.24875684E-02, 0.55350078E-03,
     + -0.62478933E-03,-0.39959382E-03, 0.17626804E-02, 0.62610477E-03,
     + -0.28258836E-03,-0.20391413E-03, 0.12450615E-02,-0.49451198E-02,
     +  0.28990086E-02, 0.10806918E-01,-0.24432039E-02,-0.22619427E-02,
     + -0.38589860E-03, 0.10928344E-02, 0.41529333E-03, 0.40162019E-02,
     +  0.85878989E-03,-0.71743382E-02, 0.16543329E-02,-0.55846403E-03,
     + -0.12488099E-01,-0.54125552E-03, 0.11068549E-02,-0.57716295E-02,
     +  0.13835547E-02, 0.43235342E-02, 0.10965524E-02,-0.43591973E-02,
     +  0.23629421E-02,-0.15570270E-02, 0.30496870E-02,
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
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_s4_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)        *x31    *x51
     8  +coeff(  8)            *x41*x51
      y_s4_q3en   =y_s4_q3en   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)                *x52
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x23            
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)        *x33*x41    
     2  +coeff( 20)        *x32*x42    
     3  +coeff( 21)        *x31*x43    
     4  +coeff( 22)            *x44    
     5  +coeff( 23)    *x21*x33        
     6  +coeff( 24)    *x24            
     7  +coeff( 25)        *x33*x42    
     8  +coeff( 26)        *x31*x44    
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 27)    *x22*x33        
     1  +coeff( 28)            *x42    
     2  +coeff( 29)    *x21    *x42    
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)    *x22    *x42    
     5  +coeff( 32)*x11*x21    *x41    
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)            *x45    
     8  +coeff( 35)    *x24    *x41    
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 36)*x11*x23            
     1  +coeff( 37)        *x31*x41    
     2  +coeff( 38)    *x21*x31        
     3  +coeff( 39)        *x32*x41    
     4  +coeff( 40)        *x31*x42    
     5  +coeff( 41)    *x22*x31        
     6  +coeff( 42)*x11        *x41    
     7  +coeff( 43)    *x21    *x43    
     8  +coeff( 44)    *x22*x31*x41    
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 45)    *x23    *x41    
     1  +coeff( 46)    *x22    *x41*x51
     2  +coeff( 47)    *x23    *x42    
     3  +coeff( 48)    *x23    *x43    
     4  +coeff( 49)*x11*x23    *x41    
     5  +coeff( 50)                *x53
     6  +coeff( 51)            *x43*x51
     7  +coeff( 52)    *x21*x31*x43    
     8  +coeff( 53)*x11*x21    *x42    
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 54)    *x24*x31        
     1  +coeff( 55)        *x31*x44*x51
     2  +coeff( 56)    *x22    *x44    
     3  +coeff( 57)    *x22*x31*x45    
     4  +coeff( 58)        *x32        
     5  +coeff( 59)        *x31*x41*x51
     6  +coeff( 60)    *x21*x31    *x51
     7  +coeff( 61)*x11            *x51
     8  +coeff( 62)    *x21        *x52
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 63)    *x21*x31*x42    
     1  +coeff( 64)*x11        *x42    
     2  +coeff( 65)    *x21*x31*x41*x51
     3  +coeff( 66)*x11*x21*x31        
     4  +coeff( 67)    *x23*x31        
     5  +coeff( 68)    *x21    *x41*x52
     6  +coeff( 69)*x12                
     7  +coeff( 70)*x11*x21        *x51
     8  +coeff( 71)    *x22        *x52
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 72)    *x21    *x44    
     1  +coeff( 73)        *x31*x43*x51
     2  +coeff( 74)    *x22    *x43    
     3  +coeff( 75)            *x43*x52
     4  +coeff( 76)    *x22    *x42*x51
     5  +coeff( 77)*x12        *x41    
     6  +coeff( 78)    *x22    *x41*x52
     7  +coeff( 79)    *x24        *x51
     8  +coeff( 80)            *x45*x51
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 81)    *x21*x32*x42*x51
     1  +coeff( 82)    *x23*x31*x42    
     2  +coeff( 83)*x11*x22    *x42    
     3  +coeff( 84)    *x24*x31    *x51
     4  +coeff( 85)    *x22    *x45    
     5  +coeff( 86)*x12    *x31*x42    
     6  +coeff( 87)*x12*x21    *x41*x51
     7  +coeff( 88)    *x21*x31*x44*x52
     8  +coeff( 89)    *x24        *x53
      y_s4_q3en   =y_s4_q3en   
     9  +coeff( 90)    *x23*x31*x42*x52
     1  +coeff( 91)*x12*x22*x32        
     2  +coeff( 92)*x11*x24    *x42    
     3  +coeff( 93)*x12*x21    *x42*x51
     4  +coeff( 94)    *x24*x32    *x52
     5  +coeff( 95)*x11*x23    *x42*x52
c
      return
      end
      function p_s4_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1234570E-04/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.66980475E-03, 0.47018286E-03,-0.11155268E-02, 0.17519567E-01,
     +  0.74682790E-02,-0.93143288E-04,-0.80893347E-02,-0.10507341E-01,
     +  0.37767718E-03, 0.18303730E-02,-0.17054270E-02, 0.87755909E-02,
     +  0.13857298E-02,-0.11089657E-02, 0.19772451E-02,-0.12319405E-01,
     +  0.30686420E-02,-0.28034076E-02,-0.12495609E-02, 0.39460440E-02,
     + -0.54656563E-03, 0.34232534E-03,-0.98545139E-03, 0.14689231E-02,
     +  0.10158803E-02,-0.42689205E-03, 0.40958021E-02, 0.19748668E-03,
     +  0.26718210E-03,-0.73851971E-03, 0.18360832E-02, 0.20643410E-02,
     +  0.81453269E-03, 0.30953828E-02,-0.12303195E-02, 0.91886177E-05,
     +  0.17352062E-03,-0.32635216E-03,-0.12274635E-02, 0.34418155E-03,
     + -0.60905822E-04,-0.60694566E-03, 0.60433347E-03, 0.26421619E-03,
     + -0.68687246E-03, 0.49948873E-03, 0.65693847E-03, 0.81886316E-03,
     + -0.28583402E-03, 0.29237749E-03, 0.55350720E-04,-0.19691611E-03,
     +  0.29275662E-03,-0.71675371E-04, 0.16355896E-03, 0.17380615E-03,
     +  0.35826716E-03,-0.33344037E-03, 0.43895215E-03, 0.76525606E-03,
     +  0.52765664E-03,-0.63822226E-03,-0.72804955E-03,-0.18412161E-03,
     +  0.72935031E-03,-0.67647279E-03,-0.12420814E-02,-0.25647591E-03,
     +  0.39889445E-03, 0.27425087E-04,-0.50207269E-04,-0.55520526E-04,
     +  0.34425073E-03, 0.10279134E-03, 0.10337594E-02,-0.88270943E-04,
     + -0.52581603E-04,-0.38646067E-04,-0.93181974E-04,-0.98804398E-04,
     + -0.13338722E-03, 0.28059841E-03, 0.20609573E-03, 0.98011456E-04,
     +  0.20669804E-03, 0.54813434E-04,-0.84235551E-04,-0.19226108E-03,
     + -0.10263455E-02,-0.19501671E-02,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_s4_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      p_s4_q3en   =p_s4_q3en   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x23            
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)            *x43    
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)    *x24            
     3  +coeff( 21)    *x21*x31        
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)    *x22*x31        
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)        *x31*x42    
     8  +coeff( 26)            *x41*x52
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 27)            *x44    
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)        *x31*x43    
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)    *x24    *x41    
     8  +coeff( 35)    *x23    *x42    
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 36)        *x32        
     1  +coeff( 37)        *x32*x41    
     2  +coeff( 38)    *x21        *x52
     3  +coeff( 39)    *x22    *x42    
     4  +coeff( 40)        *x32*x42    
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)    *x22    *x41*x51
     7  +coeff( 43)            *x43*x51
     8  +coeff( 44)    *x21    *x41*x52
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 45)    *x21    *x44    
     1  +coeff( 46)    *x24        *x51
     2  +coeff( 47)    *x22    *x41*x52
     3  +coeff( 48)*x11*x23    *x41    
     4  +coeff( 49)    *x24*x31*x41    
     5  +coeff( 50)    *x23        *x53
     6  +coeff( 51)        *x31*x41*x51
     7  +coeff( 52)                *x53
     8  +coeff( 53)    *x21*x31*x42    
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 54)        *x31*x42*x51
     1  +coeff( 55)    *x22        *x52
     2  +coeff( 56)            *x42*x52
     3  +coeff( 57)    *x24*x31        
     4  +coeff( 58)*x11*x21    *x42    
     5  +coeff( 59)    *x22    *x43    
     6  +coeff( 60)        *x31*x44    
     7  +coeff( 61)    *x23    *x41*x51
     8  +coeff( 62)            *x43*x52
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 63)    *x22*x31*x43    
     1  +coeff( 64)    *x23*x31*x41*x51
     2  +coeff( 65)        *x31*x44*x51
     3  +coeff( 66)            *x44*x52
     4  +coeff( 67)    *x24    *x43    
     5  +coeff( 68)        *x31*x44*x52
     6  +coeff( 69)*x12*x21    *x42*x51
     7  +coeff( 70)*x11            *x51
     8  +coeff( 71)*x12                
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 72)*x11*x21*x31        
     1  +coeff( 73)    *x23*x31        
     2  +coeff( 74)*x11        *x42    
     3  +coeff( 75)    *x21    *x43    
     4  +coeff( 76)    *x22*x31    *x51
     5  +coeff( 77)*x12        *x41    
     6  +coeff( 78)*x11*x22    *x41    
     7  +coeff( 79)    *x21*x31*x43    
     8  +coeff( 80)*x11*x21    *x41*x51
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 81)    *x22*x31*x41*x51
     1  +coeff( 82)        *x31*x43*x51
     2  +coeff( 83)    *x23        *x52
     3  +coeff( 84)    *x21*x31*x41*x52
     4  +coeff( 85)    *x22        *x53
     5  +coeff( 86)*x12    *x32        
     6  +coeff( 87)    *x23*x33        
     7  +coeff( 88)*x11*x22    *x42    
     8  +coeff( 89)    *x23*x31*x42    
      p_s4_q3en   =p_s4_q3en   
     9  +coeff( 90)    *x23    *x43    
c
      return
      end
      function l_s4_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.8236276E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10469259E-01,-0.23309463E+00,-0.28983921E-01, 0.85280510E-02,
     + -0.25510063E-01, 0.60960684E-01,-0.11059019E-01,-0.24662714E-01,
     + -0.21840178E-02, 0.80571156E-02,-0.44715941E-01,-0.20845875E-01,
     +  0.30929128E-01,-0.55066254E-02,-0.18839854E-02,-0.15742261E-02,
     +  0.41933670E-02, 0.74219750E-02,-0.52914890E-02, 0.69759325E-02,
     + -0.94658425E-02,-0.27996852E-03,-0.64966083E-03, 0.26845322E-02,
     + -0.18342952E-02,-0.48989160E-02,-0.10250163E-01, 0.34372339E-02,
     +  0.94474582E-02,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_s4_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      l_s4_q3en   =l_s4_q3en   
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)            *x42    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)*x11        *x41    
      l_s4_q3en   =l_s4_q3en   
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x23*x31        
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)        *x31        
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x24            
     7  +coeff( 25)            *x42*x52
     8  +coeff( 26)    *x23*x31*x41    
      l_s4_q3en   =l_s4_q3en   
     9  +coeff( 27)    *x22    *x43    
     1  +coeff( 28)*x11*x23    *x42    
     2  +coeff( 29)    *x24    *x43    
c
      return
      end
      function x_s4_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1451097E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.55977432E-02, 0.24403606E+00,-0.31536976E-02, 0.11976951E+00,
     +  0.14674178E-01,-0.69479637E-01, 0.26258392E-01, 0.26517402E-01,
     + -0.11132443E-01,-0.87028239E-02,-0.35505444E-01, 0.52033599E-01,
     + -0.29277080E-02, 0.87151388E-02,-0.81933625E-02, 0.13999910E-01,
     + -0.10263134E-01, 0.15887178E-02,-0.33253958E-02, 0.61777020E-02,
     + -0.13992612E-02,-0.47581471E-02, 0.26012333E-01, 0.10172082E-01,
     + -0.67906163E-03,-0.70967479E-02,-0.11008090E-01, 0.10275433E-01,
     + -0.21658430E-01, 0.11176680E-01,-0.63956836E-02,-0.14855781E-02,
     + -0.20404204E-02, 0.19595448E-02, 0.22457264E-01,-0.43844800E-02,
     +  0.20931922E-02, 0.32429325E-02, 0.16449547E-02,-0.96003601E-03,
     + -0.77277038E-03,-0.55185058E-02, 0.30398532E-03,-0.16146150E-02,
     + -0.11663295E-01, 0.75943395E-03,-0.19403201E-02,-0.79810256E-02,
     + -0.28915456E-02, 0.10908956E-01,-0.15130949E-02, 0.72233314E-02,
     + -0.41139456E-02, 0.15376909E-02,-0.45994110E-03, 0.96917397E-03,
     +  0.17923415E-02,-0.30372492E-02,-0.44610887E-02, 0.64381240E-02,
     +  0.59637125E-02, 0.84706210E-02, 0.43457672E-02, 0.93773156E-02,
     + -0.40101898E-02, 0.85019879E-02, 0.37025873E-01,-0.40606977E-02,
     + -0.10550403E-02, 0.27139367E-04,-0.21878369E-02, 0.79869470E-02,
     +  0.31778968E-02, 0.42764275E-03, 0.84333058E-03,-0.41181102E-03,
     +  0.56996229E-02,-0.74714096E-03, 0.71325083E-02,-0.10005988E-02,
     +  0.14769239E-02,-0.29025022E-02, 0.22245108E-02,-0.31516133E-02,
     +  0.69876539E-03,-0.69552637E-03,-0.96825613E-02, 0.14737174E-02,
     + -0.96413313E-03,-0.41856412E-02,
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
      x_s4_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      x_s4_dex    =x_s4_dex    
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)                *x52
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)            *x42    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21*x31*x41    
      x_s4_dex    =x_s4_dex    
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x24            
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)        *x31        
     8  +coeff( 26)    *x21    *x43    
      x_s4_dex    =x_s4_dex    
     9  +coeff( 27)    *x24    *x41    
     1  +coeff( 28)    *x23*x31*x41    
     2  +coeff( 29)    *x24    *x42    
     3  +coeff( 30)*x11*x23    *x41*x51
     4  +coeff( 31)*x11*x24*x32*x41*x51
     5  +coeff( 32)        *x31*x41    
     6  +coeff( 33)            *x41*x51
     7  +coeff( 34)    *x22*x31        
     8  +coeff( 35)    *x22    *x42    
      x_s4_dex    =x_s4_dex    
     9  +coeff( 36)    *x21*x31*x42    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)    *x22    *x41*x51
     3  +coeff( 39)*x11*x21            
     4  +coeff( 40)    *x24*x31        
     5  +coeff( 41)    *x22*x32*x41    
     6  +coeff( 42)    *x21    *x44    
     7  +coeff( 43)*x11            *x51
     8  +coeff( 44)    *x22*x32    *x51
      x_s4_dex    =x_s4_dex    
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)    *x22*x31*x41*x52
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)    *x24*x31*x42    
     4  +coeff( 49)*x11*x21    *x41*x51
     5  +coeff( 50)    *x23    *x41*x53
     6  +coeff( 51)            *x44*x53
     7  +coeff( 52)    *x21    *x42*x54
     8  +coeff( 53)*x11*x24            
      x_s4_dex    =x_s4_dex    
     9  +coeff( 54)*x11*x21*x32*x41    
     1  +coeff( 55)    *x24*x32*x42    
     2  +coeff( 56)*x11*x23        *x51
     3  +coeff( 57)*x11        *x43*x51
     4  +coeff( 58)*x11*x23    *x42    
     5  +coeff( 59)    *x24*x33*x42    
     6  +coeff( 60)    *x21*x34*x44    
     7  +coeff( 61)    *x24*x32    *x53
     8  +coeff( 62)    *x22*x31*x42*x54
      x_s4_dex    =x_s4_dex    
     9  +coeff( 63)*x11*x23*x31*x42    
     1  +coeff( 64)    *x23*x33*x44    
     2  +coeff( 65)*x11*x24    *x41*x51
     3  +coeff( 66)*x11*x22    *x42*x52
     4  +coeff( 67)*x11*x24*x31*x43*x54
     5  +coeff( 68)            *x43    
     6  +coeff( 69)    *x21*x31    *x51
     7  +coeff( 70)            *x42*x51
     8  +coeff( 71)    *x21        *x52
      x_s4_dex    =x_s4_dex    
     9  +coeff( 72)    *x22*x31*x41    
     1  +coeff( 73)    *x21    *x42*x51
     2  +coeff( 74)    *x21*x31    *x52
     3  +coeff( 75)    *x21    *x41*x52
     4  +coeff( 76)*x11    *x31        
     5  +coeff( 77)    *x22*x31*x42    
     6  +coeff( 78)        *x33*x42    
     7  +coeff( 79)    *x22    *x43    
     8  +coeff( 80)        *x31*x44    
      x_s4_dex    =x_s4_dex    
     9  +coeff( 81)    *x23*x31    *x51
     1  +coeff( 82)            *x44*x51
     2  +coeff( 83)    *x23        *x52
     3  +coeff( 84)    *x21    *x41*x53
     4  +coeff( 85)*x11*x21    *x41    
     5  +coeff( 86)*x11    *x31*x41    
     6  +coeff( 87)    *x24*x31*x41    
     7  +coeff( 88)    *x23*x32*x41    
     8  +coeff( 89)*x11        *x42    
      x_s4_dex    =x_s4_dex    
     9  +coeff( 90)    *x22*x32*x42    
c
      return
      end
      function t_s4_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.5920748E+00/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.22349658E-02,-0.80073439E-01, 0.15768836E-02, 0.27741617E-01,
     +  0.23626864E-02, 0.19308940E-01,-0.46316725E-02,-0.76448978E-02,
     +  0.24863563E-02,-0.26737277E-02, 0.10336210E-01,-0.14239429E-01,
     + -0.18307650E-02,-0.20203406E-02,-0.12055918E-02,-0.23460418E-02,
     +  0.13520006E-02, 0.27290811E-02, 0.13451842E-02,-0.16823307E-02,
     + -0.19012425E-02,-0.27967258E-02,-0.67872759E-02,-0.12805064E-03,
     +  0.26034392E-03,-0.67440083E-03, 0.10020057E-02, 0.30311963E-02,
     + -0.26292994E-02,-0.17026106E-02,-0.31231786E-03, 0.34066881E-02,
     + -0.24863121E-02, 0.11648388E-02, 0.19705221E-02,-0.13836939E-01,
     + -0.46703403E-03, 0.64717693E-03,-0.31424276E-03,-0.28068025E-04,
     +  0.42470612E-03, 0.11990516E-02, 0.26254624E-02, 0.63933624E-03,
     + -0.68989568E-04, 0.19703335E-02, 0.31753026E-02,-0.73895475E-03,
     +  0.15086997E-02, 0.32603036E-03, 0.10635020E-02,-0.65788929E-03,
     +  0.12552848E-02,-0.26783238E-02, 0.37584058E-03,-0.22261975E-03,
     + -0.30451585E-03,-0.28484864E-02,-0.20878664E-02,-0.24773295E-02,
     + -0.29524590E-02,-0.11656798E-02, 0.24924427E-02, 0.23086995E-02,
     + -0.48282193E-02,-0.15600384E-02,-0.14884037E-02,-0.10652646E-01,
     +  0.94207497E-02,-0.37070110E-01, 0.23914251E-01, 0.31462133E-01,
     +  0.22416566E-03, 0.19467155E-03,-0.10169834E-02,-0.19358300E-02,
     + -0.18427988E-03,-0.46761450E-03,-0.20294345E-02, 0.29742488E-03,
     + -0.38637946E-03,-0.63496525E-03,-0.26796167E-03,-0.14839854E-03,
     +  0.90536487E-03, 0.19498731E-02,-0.72648778E-03, 0.82322775E-03,
     + -0.26574323E-03,-0.97496522E-03,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_s4_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      t_s4_dex    =t_s4_dex    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x22            
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)                *x52
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11        *x41    
      t_s4_dex    =t_s4_dex    
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x21    *x41*x51
     8  +coeff( 26)            *x41*x52
      t_s4_dex    =t_s4_dex    
     9  +coeff( 27)    *x24            
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)    *x23*x31*x41    
     3  +coeff( 30)    *x22    *x43    
     4  +coeff( 31)    *x23        *x52
     5  +coeff( 32)    *x24    *x42    
     6  +coeff( 33)    *x24    *x41*x51
     7  +coeff( 34)    *x24    *x43    
     8  +coeff( 35)*x11*x23    *x43    
      t_s4_dex    =t_s4_dex    
     9  +coeff( 36)*x11*x24    *x44*x54
     1  +coeff( 37)*x11*x21            
     2  +coeff( 38)    *x21        *x52
     3  +coeff( 39)*x11*x21    *x41    
     4  +coeff( 40)        *x33*x41    
     5  +coeff( 41)*x11        *x42    
     6  +coeff( 42)    *x21*x31*x42    
     7  +coeff( 43)    *x22    *x41*x51
     8  +coeff( 44)*x11*x23            
      t_s4_dex    =t_s4_dex    
     9  +coeff( 45)*x11    *x33        
     1  +coeff( 46)    *x21    *x44    
     2  +coeff( 47)    *x23    *x41*x51
     3  +coeff( 48)        *x31*x43*x51
     4  +coeff( 49)    *x21    *x42*x52
     5  +coeff( 50)        *x31*x41*x53
     6  +coeff( 51)*x11*x24            
     7  +coeff( 52)*x11*x22    *x42    
     8  +coeff( 53)    *x22*x32*x42    
      t_s4_dex    =t_s4_dex    
     9  +coeff( 54)    *x23    *x43    
     1  +coeff( 55)*x11        *x42*x52
     2  +coeff( 56)*x11*x24*x31        
     3  +coeff( 57)    *x24*x32    *x51
     4  +coeff( 58)*x11*x23    *x41*x51
     5  +coeff( 59)    *x22*x31*x42*x52
     6  +coeff( 60)    *x23    *x41*x53
     7  +coeff( 61)    *x21    *x42*x54
     8  +coeff( 62)*x11*x22*x31*x42*x52
      t_s4_dex    =t_s4_dex    
     9  +coeff( 63)    *x24*x33*x43    
     1  +coeff( 64)*x11*x24*x32*x41*x51
     2  +coeff( 65)    *x23    *x43*x54
     3  +coeff( 66)*x12*x23*x32*x42    
     4  +coeff( 67)*x12*x22*x33*x41*x51
     5  +coeff( 68)*x12*x21    *x44*x53
     6  +coeff( 69)*x11*x24*x33*x44    
     7  +coeff( 70)*x11*x24*x31*x44*x52
     8  +coeff( 71)*x12*x23    *x44*x53
      t_s4_dex    =t_s4_dex    
     9  +coeff( 72)*x11*x24*x33*x44*x54
     1  +coeff( 73)        *x31        
     2  +coeff( 74)*x11    *x31        
     3  +coeff( 75)    *x22*x31        
     4  +coeff( 76)    *x22    *x41    
     5  +coeff( 77)        *x31*x42    
     6  +coeff( 78)            *x43    
     7  +coeff( 79)    *x22    *x42    
     8  +coeff( 80)        *x31*x43    
      t_s4_dex    =t_s4_dex    
     9  +coeff( 81)    *x23        *x51
     1  +coeff( 82)    *x21    *x42*x51
     2  +coeff( 83)            *x42*x52
     3  +coeff( 84)*x12*x21            
     4  +coeff( 85)    *x24*x31        
     5  +coeff( 86)    *x24    *x41    
     6  +coeff( 87)    *x24        *x51
     7  +coeff( 88)*x11*x21    *x41*x51
     8  +coeff( 89)*x11        *x42*x51
      t_s4_dex    =t_s4_dex    
     9  +coeff( 90)    *x21    *x43*x51
c
      return
      end
      function y_s4_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 96)
      data ncoeff/ 95/
      data avdat/ -0.9122630E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.32477973E-02,-0.20321016E-02, 0.13240196E+00, 0.12931012E-03,
     +  0.41493841E-01, 0.19410821E-02, 0.88306908E-02,-0.44322040E-01,
     +  0.14786102E-02, 0.32268353E-01,-0.52829669E-03,-0.46486974E-01,
     + -0.81016570E-02, 0.37251906E-02, 0.10036688E-01,-0.57941895E-01,
     + -0.18111564E-02,-0.61052707E-02, 0.10036280E-01,-0.70998049E-02,
     +  0.24462773E-01,-0.60878356E-04, 0.20212423E-01,-0.13738135E-04,
     + -0.20264040E-03, 0.45063971E-02, 0.72322907E-02,-0.49741161E-02,
     +  0.12737432E-01,-0.65110987E-02,-0.50916532E-02, 0.13155339E-02,
     +  0.16413935E-01, 0.43899785E-02, 0.16598008E-03,-0.20147590E-02,
     +  0.88072295E-03, 0.95674669E-03,-0.15884230E-02, 0.67416220E-02,
     + -0.20217595E-02, 0.10242884E-01,-0.60220286E-02, 0.19786322E-03,
     +  0.52979644E-02, 0.19072894E-02,-0.26680093E-03,-0.21092277E-02,
     +  0.73855245E-02,-0.44146287E-02,-0.22495864E-02, 0.16988346E-02,
     + -0.66918862E-03, 0.27517357E-02,-0.63807019E-02,-0.83818780E-02,
     + -0.33841177E-03, 0.20329740E-03,-0.50945766E-03, 0.39396510E-02,
     +  0.43021090E-03,-0.46092650E-03,-0.33595128E-03, 0.89387072E-03,
     + -0.13672261E-03,-0.16183370E-03, 0.59136353E-03,-0.59109903E-03,
     +  0.27242429E-02, 0.10626054E-01,-0.23414872E-02,-0.21032616E-02,
     +  0.89852198E-03, 0.34044783E-02, 0.17388639E-03,-0.71919952E-02,
     + -0.15930063E-02, 0.19403710E-02, 0.16319117E-02,-0.11874816E-01,
     + -0.30799946E-02, 0.10739617E-02,-0.44326033E-02, 0.19734667E-02,
     + -0.32839098E-02, 0.41306242E-02,-0.78526838E-02, 0.33724812E-03,
     + -0.10013775E-03, 0.82903355E-03, 0.91511471E-03, 0.45516412E-03,
     + -0.28697998E-03,-0.57279249E-03,-0.22779255E-03,
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
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      y_s4_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)        *x31*x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21    *x41    
      y_s4_dex    =y_s4_dex    
     9  +coeff(  9)        *x31    *x51
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11                
     3  +coeff( 12)    *x22            
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)            *x41*x52
      y_s4_dex    =y_s4_dex    
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)            *x44    
     4  +coeff( 22)    *x21*x33        
     5  +coeff( 23)    *x24            
     6  +coeff( 24)        *x33*x42    
     7  +coeff( 25)    *x22*x33        
     8  +coeff( 26)        *x31*x42    
      y_s4_dex    =y_s4_dex    
     9  +coeff( 27)    *x21    *x42    
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)        *x31*x43    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)    *x24    *x41    
     7  +coeff( 34)*x11*x23            
     8  +coeff( 35)        *x32        
      y_s4_dex    =y_s4_dex    
     9  +coeff( 36)    *x21*x31        
     1  +coeff( 37)        *x32*x41    
     2  +coeff( 38)*x11        *x41    
     3  +coeff( 39)    *x21    *x41*x51
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)    *x22    *x41*x51
     6  +coeff( 42)            *x45    
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)    *x22*x33*x41    
      y_s4_dex    =y_s4_dex    
     9  +coeff( 45)*x11*x23    *x41    
     1  +coeff( 46)        *x32*x42    
     2  +coeff( 47)            *x43*x51
     3  +coeff( 48)    *x22*x31*x41    
     4  +coeff( 49)        *x31*x44    
     5  +coeff( 50)    *x21    *x44    
     6  +coeff( 51)*x11*x21    *x42    
     7  +coeff( 52)    *x24*x31        
     8  +coeff( 53)    *x21*x33*x42    
      y_s4_dex    =y_s4_dex    
     9  +coeff( 54)        *x31*x44*x51
     1  +coeff( 55)    *x22    *x44    
     2  +coeff( 56)    *x22*x31*x45    
     3  +coeff( 57)        *x31*x41*x51
     4  +coeff( 58)*x11            *x51
     5  +coeff( 59)                *x53
     6  +coeff( 60)    *x21    *x43    
     7  +coeff( 61)*x11        *x42    
     8  +coeff( 62)    *x21*x31*x41*x51
      y_s4_dex    =y_s4_dex    
     9  +coeff( 63)*x11*x21*x31        
     1  +coeff( 64)    *x23*x31        
     2  +coeff( 65)*x12                
     3  +coeff( 66)*x11*x21        *x51
     4  +coeff( 67)    *x22        *x52
     5  +coeff( 68)    *x21*x31*x43    
     6  +coeff( 69)        *x31*x43*x51
     7  +coeff( 70)    *x22    *x43    
     8  +coeff( 71)            *x43*x52
      y_s4_dex    =y_s4_dex    
     9  +coeff( 72)    *x22    *x42*x51
     1  +coeff( 73)    *x24        *x51
     2  +coeff( 74)            *x45*x51
     3  +coeff( 75)    *x21*x32*x42*x51
     4  +coeff( 76)    *x23    *x43    
     5  +coeff( 77)    *x22*x31*x42*x51
     6  +coeff( 78)    *x21    *x43*x52
     7  +coeff( 79)*x11*x22    *x42    
     8  +coeff( 80)    *x22    *x45    
      y_s4_dex    =y_s4_dex    
     9  +coeff( 81)    *x23*x31*x43    
     1  +coeff( 82)*x12*x21    *x41*x51
     2  +coeff( 83)*x11*x24    *x42    
     3  +coeff( 84)*x12*x21    *x42*x51
     4  +coeff( 85)    *x24*x32*x43    
     5  +coeff( 86)*x11*x23    *x42*x52
     6  +coeff( 87)*x12*x23*x31*x44*x52
     7  +coeff( 88)    *x21*x31*x41    
     8  +coeff( 89)    *x21*x31    *x51
      y_s4_dex    =y_s4_dex    
     9  +coeff( 90)    *x21    *x42*x51
     1  +coeff( 91)        *x32*x43    
     2  +coeff( 92)    *x21        *x53
     3  +coeff( 93)*x11*x21*x31*x41    
     4  +coeff( 94)    *x22*x31*x41*x51
     5  +coeff( 95)*x12        *x41    
c
      return
      end
      function p_s4_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.4489905E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.64276182E-03, 0.53321628E-03,-0.11769668E-02, 0.12637868E-01,
     +  0.62110140E-02,-0.67506342E-04,-0.67436732E-02,-0.10346607E-01,
     +  0.37950766E-03, 0.16078211E-02,-0.19439793E-02, 0.29412122E-03,
     +  0.76828785E-02, 0.12418512E-02,-0.90058701E-03, 0.19743093E-02,
     + -0.10127721E-01, 0.14399877E-02, 0.29189300E-02,-0.22841368E-02,
     + -0.16411579E-02, 0.33568249E-02,-0.44835507E-03,-0.76859986E-03,
     +  0.91516029E-03,-0.43420109E-03,-0.31760090E-03, 0.34253886E-02,
     +  0.22615251E-03, 0.22241699E-03,-0.63690299E-03, 0.18636980E-02,
     +  0.78192231E-03, 0.16452221E-02,-0.10972257E-02, 0.98994968E-03,
     +  0.63191610E-03, 0.22743435E-02,-0.11067368E-02,-0.84471150E-03,
     +  0.18057859E-04, 0.15061536E-03, 0.39920845E-03,-0.13472706E-03,
     + -0.37391463E-03,-0.11520440E-02,-0.56149260E-04, 0.25459402E-03,
     +  0.58978101E-04,-0.65474084E-03, 0.43661121E-03,-0.77292870E-03,
     +  0.23143459E-03, 0.73322857E-03,-0.97176136E-03,-0.19942389E-02,
     +  0.50434464E-04, 0.34140452E-04,-0.45420169E-04, 0.30950308E-03,
     + -0.11496087E-03,-0.69249836E-05, 0.27724149E-03,-0.47816171E-04,
     + -0.12069695E-03,-0.40231098E-03,-0.24517821E-03, 0.50949567E-03,
     +  0.46218949E-03,-0.20007558E-03, 0.37000354E-03, 0.24800265E-03,
     + -0.50639838E-03, 0.55588163E-04, 0.42534879E-03, 0.45010846E-03,
     +  0.36138252E-03, 0.76368119E-03,-0.20278366E-03,-0.53004434E-04,
     +  0.67646382E-04, 0.20771350E-03, 0.26688993E-03, 0.38745169E-04,
     +  0.68800720E-04, 0.57261128E-04,-0.10016609E-03, 0.22762213E-04,
     + -0.49580303E-04,-0.70694237E-04,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_s4_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21    *x41    
      p_s4_dex    =p_s4_dex    
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22    *x41    
      p_s4_dex    =p_s4_dex    
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)            *x43    
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x24            
     5  +coeff( 23)    *x21*x31        
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)        *x31*x42    
     8  +coeff( 26)    *x21        *x52
      p_s4_dex    =p_s4_dex    
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)            *x44    
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)    *x23    *x41    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)        *x31*x43    
     8  +coeff( 35)    *x22    *x41*x51
      p_s4_dex    =p_s4_dex    
     9  +coeff( 36)            *x43*x51
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)    *x24    *x41    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)        *x34*x42    
     6  +coeff( 42)        *x32*x41    
     7  +coeff( 43)            *x42*x51
     8  +coeff( 44)                *x53
      p_s4_dex    =p_s4_dex    
     9  +coeff( 45)    *x22*x31*x41    
     1  +coeff( 46)    *x22    *x42    
     2  +coeff( 47)*x11*x21        *x51
     3  +coeff( 48)    *x23        *x51
     4  +coeff( 49)        *x31*x42*x51
     5  +coeff( 50)    *x21    *x44    
     6  +coeff( 51)    *x24        *x51
     7  +coeff( 52)    *x22    *x42*x51
     8  +coeff( 53)    *x22    *x41*x52
      p_s4_dex    =p_s4_dex    
     9  +coeff( 54)*x11*x23    *x41    
     1  +coeff( 55)    *x23*x31*x42    
     2  +coeff( 56)    *x23    *x43    
     3  +coeff( 57)*x11            *x51
     4  +coeff( 58)        *x31*x41*x51
     5  +coeff( 59)*x12                
     6  +coeff( 60)    *x23*x31        
     7  +coeff( 61)    *x21*x31*x41*x51
     8  +coeff( 62)    *x21    *x42*x51
      p_s4_dex    =p_s4_dex    
     9  +coeff( 63)    *x24*x31        
     1  +coeff( 64)*x12        *x41    
     2  +coeff( 65)*x11*x22    *x41    
     3  +coeff( 66)*x11*x21    *x42    
     4  +coeff( 67)    *x22*x31*x42    
     5  +coeff( 68)        *x31*x44    
     6  +coeff( 69)    *x23    *x41*x51
     7  +coeff( 70)    *x22*x31*x41*x51
     8  +coeff( 71)        *x31*x43*x51
      p_s4_dex    =p_s4_dex    
     9  +coeff( 72)    *x23        *x52
     1  +coeff( 73)            *x43*x52
     2  +coeff( 74)*x12    *x32        
     3  +coeff( 75)    *x24    *x41*x51
     4  +coeff( 76)        *x31*x44*x51
     5  +coeff( 77)*x12*x21    *x42*x51
     6  +coeff( 78)*x11*x23    *x42*x52
     7  +coeff( 79)*x12*x24*x31    *x51
     8  +coeff( 80)*x11*x21*x31        
      p_s4_dex    =p_s4_dex    
     9  +coeff( 81)*x11        *x42    
     1  +coeff( 82)    *x21*x31*x42    
     2  +coeff( 83)        *x32*x42    
     3  +coeff( 84)*x11        *x41*x51
     4  +coeff( 85)    *x22        *x52
     5  +coeff( 86)    *x21        *x53
     6  +coeff( 87)    *x21*x31*x43    
     7  +coeff( 88)        *x34    *x51
     8  +coeff( 89)*x11*x21    *x41*x51
      p_s4_dex    =p_s4_dex    
     9  +coeff( 90)*x11        *x42*x51
c
      return
      end
      function l_s4_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.1035815E-01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13342048E-01,-0.35492027E+00, 0.35175442E-03,-0.88429898E-01,
     +  0.14054799E-01,-0.33146918E-01, 0.94665863E-01,-0.26158322E-01,
     + -0.38427453E-01, 0.12796983E-01,-0.69221266E-01,-0.35362318E-01,
     +  0.49214393E-01,-0.82433047E-02, 0.63858451E-02,-0.15340879E-01,
     +  0.12206378E-01,-0.86900946E-02, 0.11476615E-01,-0.13637267E-01,
     +  0.13682231E-02,-0.79550111E-03, 0.67717195E-02, 0.11959628E-01,
     + -0.90008490E-02,-0.75488249E-02,-0.10082085E-01, 0.97725224E-02,
     + -0.75290905E-03,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_s4_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_s4_dex    =l_s4_dex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x22    *x41    
     8  +coeff( 17)    *x21*x31*x41    
      l_s4_dex    =l_s4_dex    
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)*x11*x22    *x41    
     3  +coeff( 21)                *x52
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x24            
     6  +coeff( 24)    *x24    *x41    
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)*x11*x23    *x41*x51
      l_s4_dex    =l_s4_dex    
     9  +coeff( 27)*x11*x23*x31*x43    
     1  +coeff( 28)*x11*x24    *x41*x53
     2  +coeff( 29)    *x22        *x51
c
      return
      end
      function x_s4_dmd    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5098173E+01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.74123306E-03,-0.94743609E-01, 0.22878356E-01,-0.36909720E-02,
     + -0.84635578E-02, 0.29598838E-02,-0.16519180E-01, 0.29441509E-02,
     +  0.11398600E-01, 0.14480717E-02,-0.61673843E-02, 0.34760828E-02,
     + -0.34671125E-02,-0.55952364E-03, 0.16835068E-02,-0.20860161E-02,
     +  0.15979338E-02,-0.97830584E-02,-0.33263704E-02, 0.35858722E-02,
     + -0.74556295E-03, 0.47829715E-02,-0.34776065E-02,-0.27974634E-03,
     + -0.21429146E-02, 0.51661226E-03, 0.49234269E-03,-0.38291251E-02,
     + -0.13489093E-02,-0.10810905E-01, 0.30865913E-03,-0.18313103E-02,
     + -0.47383588E-02, 0.13640245E-02, 0.23937646E-03, 0.39132452E-03,
     +  0.22132503E-03, 0.17698571E-02, 0.22199900E-03, 0.54854508E-02,
     + -0.11523253E-02,-0.63564745E-03,-0.16761553E-02,-0.79321995E-03,
     +  0.72372577E-03,-0.16152419E-03, 0.12512812E-02, 0.21861091E-02,
     + -0.17350540E-02, 0.55538593E-02, 0.90673557E-05,-0.59342292E-04,
     +  0.20681105E-02,-0.11823300E-03,-0.85941894E-03,-0.64828113E-03,
     + -0.53498517E-02, 0.15034906E-02,-0.86252019E-03, 0.21333571E-02,
     +  0.50485046E-02,-0.28042684E-02,-0.17817829E-01, 0.14647614E-01,
     +  0.40996043E-03,-0.76319498E-04,-0.49847698E-04, 0.76890667E-03,
     +  0.20127946E-03,-0.77561002E-04,-0.24946895E-03, 0.27990193E-03,
     +  0.15359772E-03,-0.90305955E-03,-0.21501764E-03, 0.10982765E-02,
     +  0.77986077E-03,-0.14338442E-02, 0.14242028E-02, 0.12983238E-03,
     + -0.51007955E-03,-0.69310435E-03, 0.10855313E-03,-0.36265232E-03,
     + -0.43449603E-03,-0.93621638E-03, 0.30782179E-03, 0.32600926E-03,
     + -0.17552474E-03,-0.20399684E-03,
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
      x_s4_dmd    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x23            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x21*x31        
      x_s4_dmd    =x_s4_dmd    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21*x31*x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x22            
     6  +coeff( 15)    *x24            
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)*x11        *x41    
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)*x11*x21            
     4  +coeff( 22)    *x24    *x41    
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x21    *x43*x51
     8  +coeff( 26)*x11*x23            
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 27)    *x23*x31*x43*x51
     1  +coeff( 28)*x11*x23    *x41*x51
     2  +coeff( 29)*x11*x23*x31*x43    
     3  +coeff( 30)*x11*x24    *x44*x52
     4  +coeff( 31)        *x31        
     5  +coeff( 32)    *x22*x31        
     6  +coeff( 33)    *x22    *x42    
     7  +coeff( 34)    *x21*x31*x42    
     8  +coeff( 35)        *x32*x42    
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 36)            *x44    
     1  +coeff( 37)*x11    *x31        
     2  +coeff( 38)    *x24*x31        
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)    *x24    *x42    
     5  +coeff( 41)    *x22    *x44    
     6  +coeff( 42)        *x31*x44*x51
     7  +coeff( 43)    *x23    *x41*x52
     8  +coeff( 44)    *x21    *x42*x53
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 45)        *x31*x42*x53
     1  +coeff( 46)*x11*x22*x31        
     2  +coeff( 47)*x11*x21    *x41*x51
     3  +coeff( 48)    *x24*x31*x41*x51
     4  +coeff( 49)    *x22*x33*x41*x51
     5  +coeff( 50)    *x23    *x43*x51
     6  +coeff( 51)        *x32*x42*x53
     7  +coeff( 52)    *x23        *x54
     8  +coeff( 53)*x11*x24            
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 54)*x11*x22    *x41*x51
     1  +coeff( 55)    *x24*x32*x41*x51
     2  +coeff( 56)*x11        *x43*x51
     3  +coeff( 57)    *x23    *x43*x52
     4  +coeff( 58)*x11*x23    *x42    
     5  +coeff( 59)*x11*x21*x31*x42*x51
     6  +coeff( 60)*x11*x24    *x41*x51
     7  +coeff( 61)    *x24*x32*x42*x52
     8  +coeff( 62)*x11*x22*x31*x42*x52
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 63)    *x23*x31*x44*x53
     1  +coeff( 64)    *x23*x33*x44*x53
     2  +coeff( 65)            *x42    
     3  +coeff( 66)        *x31    *x51
     4  +coeff( 67)                *x52
     5  +coeff( 68)            *x43    
     6  +coeff( 69)    *x22        *x51
     7  +coeff( 70)        *x32    *x51
     8  +coeff( 71)            *x42*x51
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 72)    *x21        *x52
     1  +coeff( 73)        *x31*x43    
     2  +coeff( 74)    *x21    *x42*x51
     3  +coeff( 75)    *x21*x31    *x52
     4  +coeff( 76)    *x21    *x41*x52
     5  +coeff( 77)    *x22*x31*x42    
     6  +coeff( 78)    *x22    *x43    
     7  +coeff( 79)    *x21    *x44    
     8  +coeff( 80)    *x23*x31    *x51
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 81)    *x23    *x41*x51
     1  +coeff( 82)        *x31*x43*x51
     2  +coeff( 83)            *x44*x51
     3  +coeff( 84)    *x21*x31*x41*x52
     4  +coeff( 85)        *x31*x42*x52
     5  +coeff( 86)            *x43*x52
     6  +coeff( 87)        *x31*x41*x53
     7  +coeff( 88)            *x41*x54
     8  +coeff( 89)    *x22*x34        
      x_s4_dmd    =x_s4_dmd    
     9  +coeff( 90)*x11*x21    *x41    
c
      return
      end
      function t_s4_dmd    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3606076E+01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.47797952E-01, 0.53976226E+00,-0.18068455E-01, 0.59466134E-02,
     + -0.25113015E-01, 0.13723616E+00,-0.15868871E+00, 0.22089923E-01,
     +  0.88368110E-01, 0.10120524E+00,-0.19090891E-01, 0.26507992E-01,
     +  0.12771260E-01,-0.31451475E-01,-0.78273945E-01,-0.49862736E-02,
     + -0.20473106E-01, 0.19692117E-01,-0.19581979E-02,-0.11485100E-01,
     +  0.90795709E-02, 0.11823181E-01, 0.34739893E-01, 0.68618864E-01,
     +  0.65552630E-02, 0.18322530E-02, 0.69565359E-02,-0.14451728E-01,
     + -0.32718349E-01,-0.52225980E-03, 0.12662000E-01, 0.21964513E-01,
     +  0.11179719E-01, 0.27516901E-01, 0.13044018E-01, 0.20171542E-01,
     +  0.25342101E-01,-0.21841249E-02, 0.51820895E-03, 0.68390905E-03,
     + -0.13590634E-02, 0.31647317E-02, 0.55047777E-02,-0.14589083E-02,
     +  0.42077721E-03,-0.21903443E-02,-0.13397035E-01,-0.48386301E-02,
     +  0.32158997E-02, 0.18063388E-02,-0.96525857E-02,-0.85612610E-02,
     +  0.41394257E-02, 0.30270710E-02, 0.44908398E-02,-0.10280288E-01,
     + -0.44142399E-02,-0.55442578E-02, 0.33158142E-01, 0.39703036E-02,
     + -0.33040920E-02, 0.13245905E-01, 0.39327708E-02, 0.22678024E-02,
     +  0.13753094E-02, 0.57143299E-02,-0.13333386E-02,-0.18176431E-01,
     + -0.53072725E-02,-0.13448779E-01,-0.18931592E-01,-0.35810929E-02,
     +  0.35214222E-02, 0.39261300E-02,-0.11314531E-01, 0.41962173E-02,
     +  0.36274928E-02,-0.91567980E-02, 0.57045473E-02,-0.50128568E-01,
     +  0.12483315E-01, 0.36776751E-01,-0.67074358E-01, 0.39933939E-01,
     +  0.51472891E-01, 0.31766649E-01,-0.37206423E-01, 0.78935490E-03,
     + -0.10537406E-02, 0.12293495E-02,
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
      t_s4_dmd    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      t_s4_dmd    =t_s4_dmd    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)            *x42    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x21*x31*x41    
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 18)*x11*x22            
     1  +coeff( 19)        *x31        
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x24            
     4  +coeff( 22)    *x23*x31        
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)    *x24        *x51
     8  +coeff( 26)*x11            *x51
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 27)    *x21    *x41*x51
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x24*x31        
     4  +coeff( 31)    *x24    *x41    
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)    *x23*x31*x42    
     7  +coeff( 34)*x11*x23    *x41*x51
     8  +coeff( 35)*x11*x23*x31*x42    
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 36)*x11*x24*x31*x41*x52
     1  +coeff( 37)*x11*x22    *x42*x54
     2  +coeff( 38)*x11    *x32*x42*x54
     3  +coeff( 39)        *x31    *x51
     4  +coeff( 40)                *x52
     5  +coeff( 41)*x11    *x31        
     6  +coeff( 42)            *x43    
     7  +coeff( 43)            *x42*x51
     8  +coeff( 44)*x11    *x31*x41    
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 45)        *x33*x41    
     1  +coeff( 46)*x11        *x42    
     2  +coeff( 47)    *x21*x31*x42    
     3  +coeff( 48)    *x22    *x41*x51
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)*x11*x22*x31        
     6  +coeff( 51)*x11*x21    *x41*x51
     7  +coeff( 52)    *x23    *x41*x51
     8  +coeff( 53)        *x31*x43*x51
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 54)        *x31*x42*x52
     1  +coeff( 55)            *x43*x52
     2  +coeff( 56)*x11*x24            
     3  +coeff( 57)    *x24*x31*x41    
     4  +coeff( 58)    *x22*x32*x42    
     5  +coeff( 59)    *x23    *x43    
     6  +coeff( 60)*x11        *x43*x51
     7  +coeff( 61)*x11*x22        *x52
     8  +coeff( 62)    *x23    *x41*x52
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 63)    *x22*x31*x41*x52
     1  +coeff( 64)    *x23        *x53
     2  +coeff( 65)    *x21*x31*x41*x53
     3  +coeff( 66)    *x21    *x42*x53
     4  +coeff( 67)*x12*x22*x31        
     5  +coeff( 68)*x11*x24    *x41    
     6  +coeff( 69)*x12*x21*x31*x41    
     7  +coeff( 70)    *x24*x31*x42    
     8  +coeff( 71)    *x23    *x44    
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 72)*x12*x21        *x52
     1  +coeff( 73)    *x22*x32    *x53
     2  +coeff( 74)*x12        *x44    
     3  +coeff( 75)*x11*x24    *x41*x51
     4  +coeff( 76)        *x33*x44*x51
     5  +coeff( 77)*x11*x23        *x53
     6  +coeff( 78)*x12*x24    *x41    
     7  +coeff( 79)*x12    *x32*x43    
     8  +coeff( 80)    *x24*x31*x43*x51
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 81)*x11*x22*x31*x42*x52
     1  +coeff( 82)    *x23*x31*x42*x53
     2  +coeff( 83)    *x24*x32*x42*x52
     3  +coeff( 84)    *x22*x34*x42*x52
     4  +coeff( 85)    *x24*x33*x43*x51
     5  +coeff( 86)*x11*x24    *x43*x52
     6  +coeff( 87)    *x23*x33*x42*x53
     7  +coeff( 88)        *x31*x41    
     8  +coeff( 89)        *x32*x41    
      t_s4_dmd    =t_s4_dmd    
     9  +coeff( 90)        *x31*x42    
c
      return
      end
      function y_s4_dmd    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1116656E-01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.46724724E-02, 0.38935591E-02, 0.12910350E+00,-0.40855794E-02,
     +  0.26633138E-01, 0.56279278E-02, 0.97462079E-02,-0.30934138E-01,
     +  0.54931617E-02,-0.31999469E-01,-0.38752903E-02,-0.24304619E-03,
     +  0.11136875E-01, 0.13723940E-02, 0.35412549E-02,-0.32529021E-02,
     +  0.20140046E-02, 0.15519449E-01,-0.33352452E-02, 0.67288638E-03,
     +  0.26788346E-02, 0.41618635E-03, 0.63498132E-03,-0.79025148E-03,
     +  0.75542978E-02,-0.12870451E-01, 0.26646990E-02, 0.15502926E-01,
     +  0.15407300E-03,-0.52985805E-03,-0.30005840E-03, 0.39440690E-03,
     + -0.11970360E-02,-0.63295162E-03,-0.12998244E-02, 0.88399202E-02,
     +  0.53464584E-02, 0.29553534E-02, 0.13938037E-03,-0.89862896E-03,
     +  0.99554926E-03, 0.37679162E-02, 0.24100860E-03,-0.18137078E-02,
     +  0.14507935E-02, 0.89156562E-02,-0.48621055E-02, 0.79022822E-04,
     + -0.30095331E-03,-0.37221867E-03, 0.11439205E-02,-0.11101312E-03,
     +  0.78991416E-03,-0.40674410E-03, 0.83834275E-04,-0.12308224E-03,
     +  0.10231070E-02, 0.27180419E-03,-0.90413971E-03, 0.27453415E-01,
     +  0.13675352E-02,-0.27384327E-02, 0.20021447E-02, 0.59756078E-03,
     +  0.73952338E-03,-0.10144874E-02, 0.43328886E-03, 0.84730040E-03,
     + -0.50419639E-03,-0.72047213E-03, 0.80025819E-03,-0.26112605E-01,
     + -0.76603465E-03,-0.23645952E-01, 0.15845489E-02,-0.48015628E-03,
     +  0.15294006E-02,-0.16668937E-02, 0.68192580E-03, 0.11651227E-02,
     + -0.30534226E-03, 0.10294025E-01,-0.57760649E-03, 0.47279688E-03,
     + -0.25998574E-03, 0.67053008E-03,-0.24182550E-02,-0.72803855E-03,
     +  0.16063316E-02,-0.17381692E-02,
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
      y_s4_dmd    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22            
      y_s4_dmd    =y_s4_dmd    
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)        *x33*x41    
     4  +coeff( 13)    *x24            
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)    *x22*x31        
     8  +coeff( 17)    *x23            
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 18)            *x44    
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)        *x33*x43    
     3  +coeff( 21)*x11*x23            
     4  +coeff( 22)        *x31    *x51
     5  +coeff( 23)                *x52
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)        *x31*x43    
     8  +coeff( 26)    *x22    *x42    
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x24    *x41    
     2  +coeff( 29)        *x34*x42    
     3  +coeff( 30)    *x21    *x41    
     4  +coeff( 31)*x11                
     5  +coeff( 32)        *x32*x41    
     6  +coeff( 33)    *x21    *x42    
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)    *x22*x31*x41    
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 36)            *x45    
     1  +coeff( 37)    *x24    *x42    
     2  +coeff( 38)*x11*x23    *x41    
     3  +coeff( 39)            *x42*x51
     4  +coeff( 40)    *x21    *x41*x51
     5  +coeff( 41)*x11*x22            
     6  +coeff( 42)        *x31*x44    
     7  +coeff( 43)    *x23        *x51
     8  +coeff( 44)*x11*x21    *x42    
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 45)    *x24*x31        
     1  +coeff( 46)    *x24    *x45    
     2  +coeff( 47)    *x24*x31*x45    
     3  +coeff( 48)        *x32        
     4  +coeff( 49)    *x21*x31*x41    
     5  +coeff( 50)*x11        *x41    
     6  +coeff( 51)        *x32*x42    
     7  +coeff( 52)        *x31*x42*x51
     8  +coeff( 53)            *x43*x51
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 54)*x11*x21*x31        
     1  +coeff( 55)*x11        *x41*x51
     2  +coeff( 56)*x12                
     3  +coeff( 57)        *x32*x43    
     4  +coeff( 58)    *x22        *x52
     5  +coeff( 59)    *x21    *x44    
     6  +coeff( 60)    *x22    *x43    
     7  +coeff( 61)    *x23    *x42    
     8  +coeff( 62)    *x22    *x42*x51
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 63)*x11*x22    *x41    
     1  +coeff( 64)    *x23    *x41*x51
     2  +coeff( 65)        *x31*x44*x51
     3  +coeff( 66)    *x23*x31*x41*x51
     4  +coeff( 67)*x11*x23*x31        
     5  +coeff( 68)    *x23*x31    *x52
     6  +coeff( 69)    *x21*x31    *x54
     7  +coeff( 70)*x11*x24            
     8  +coeff( 71)    *x22*x31*x44    
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 72)    *x22    *x45    
     1  +coeff( 73)*x11*x22    *x43    
     2  +coeff( 74)    *x24    *x43    
     3  +coeff( 75)*x11*x23    *x42    
     4  +coeff( 76)    *x24*x31*x41*x51
     5  +coeff( 77)    *x24    *x42*x51
     6  +coeff( 78)*x11*x24    *x41    
     7  +coeff( 79)*x12*x21    *x41*x51
     8  +coeff( 80)    *x21*x31*x45*x51
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 81)*x12        *x41*x52
     1  +coeff( 82)    *x24    *x44    
     2  +coeff( 83)*x12*x21*x31*x42    
     3  +coeff( 84)*x12*x22*x32        
     4  +coeff( 85)*x12*x21    *x41*x52
     5  +coeff( 86)*x11*x24*x31*x42    
     6  +coeff( 87)    *x22*x31*x45*x52
     7  +coeff( 88)    *x24*x34    *x52
     8  +coeff( 89)*x12*x23*x31*x41*x51
      y_s4_dmd    =y_s4_dmd    
     9  +coeff( 90)*x12*x22*x31*x44    
c
      return
      end
      function p_s4_dmd    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.8277779E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.20174594E-02,-0.43563033E-02,-0.39666165E-01,-0.11358267E-02,
     +  0.27572329E-02,-0.17674569E-01, 0.85394137E-03,-0.27644185E-02,
     +  0.98872837E-02, 0.16438370E-02, 0.22501645E-02,-0.13910227E-05,
     + -0.92962943E-02, 0.35200503E-02,-0.35351389E-02, 0.16411743E-04,
     +  0.32067371E-02,-0.95708389E-03, 0.50267775E-03, 0.64986193E-03,
     + -0.21123255E-02,-0.30166132E-02,-0.35378663E-03, 0.69412106E-03,
     + -0.76903915E-03, 0.56849961E-03, 0.29544842E-02,-0.12278573E-02,
     + -0.10884358E-02,-0.15285665E-02, 0.90032692E-04, 0.52280031E-03,
     + -0.48470273E-03, 0.17709474E-03,-0.15720061E-03, 0.66951656E-03,
     + -0.51327399E-03,-0.11625943E-03,-0.11766203E-03, 0.85067347E-03,
     + -0.49313578E-04,-0.36160392E-03, 0.52222994E-03,-0.13268791E-03,
     +  0.70144466E-04, 0.60982406E-04, 0.23412644E-04,-0.24139167E-05,
     +  0.66678028E-03, 0.21591198E-03, 0.33497968E-03, 0.29426374E-03,
     +  0.20478200E-03,-0.44691609E-03,-0.17066338E-03,-0.71651698E-03,
     +  0.11563865E-03, 0.11588043E-02,-0.11678705E-02,-0.86947526E-04,
     +  0.42036030E-04,-0.25747568E-03,-0.41137755E-04, 0.62493586E-04,
     + -0.54874277E-03,-0.50765386E-04, 0.92216906E-04, 0.17331314E-02,
     + -0.36493533E-04, 0.24694830E-03,-0.15946310E-03, 0.41113648E-03,
     +  0.20821231E-03,-0.40880444E-04,-0.30916207E-03,-0.11640477E-03,
     + -0.56612148E-03,-0.92840160E-03,-0.32675110E-02, 0.62129716E-03,
     +  0.24364592E-03,-0.45551200E-03,-0.22228339E-03,-0.22204532E-03,
     + -0.87382097E-03, 0.21862103E-04, 0.38831131E-04, 0.16777145E-03,
     + -0.83548242E-04,-0.82963837E-04,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_s4_dmd    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_s4_dmd    =p_s4_dmd    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23*x32        
     8  +coeff( 17)    *x21            
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 18)    *x21*x31        
     1  +coeff( 19)        *x31    *x51
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22*x31*x42    
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)*x11*x21    *x41    
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)            *x44    
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)    *x24    *x43    
     4  +coeff( 31)*x11*x21            
     5  +coeff( 32)    *x21*x31*x41    
     6  +coeff( 33)        *x31*x42    
     7  +coeff( 34)    *x21        *x52
     8  +coeff( 35)                *x53
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 36)    *x24            
     1  +coeff( 37)        *x31*x43    
     2  +coeff( 38)*x11*x21        *x51
     3  +coeff( 39)*x11        *x41*x51
     4  +coeff( 40)            *x43*x51
     5  +coeff( 41)    *x24*x31        
     6  +coeff( 42)*x11*x22    *x41    
     7  +coeff( 43)    *x24        *x51
     8  +coeff( 44)    *x22*x31*x41*x51
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 45)        *x31*x41    
     1  +coeff( 46)*x11            *x51
     2  +coeff( 47)        *x31*x41*x51
     3  +coeff( 48)        *x34        
     4  +coeff( 49)    *x23    *x41    
     5  +coeff( 50)    *x22*x31*x41    
     6  +coeff( 51)*x11        *x42    
     7  +coeff( 52)        *x31*x42*x51
     8  +coeff( 53)    *x22        *x52
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 54)    *x23*x31*x41    
     1  +coeff( 55)*x11*x21    *x41*x51
     2  +coeff( 56)    *x22    *x42*x51
     3  +coeff( 57)        *x32*x42*x51
     4  +coeff( 58)            *x44*x51
     5  +coeff( 59)    *x24    *x42    
     6  +coeff( 60)*x11                
     7  +coeff( 61)*x11    *x31        
     8  +coeff( 62)    *x22*x31        
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 63)        *x32*x41    
     1  +coeff( 64)    *x21*x31    *x51
     2  +coeff( 65)            *x42*x51
     3  +coeff( 66)        *x31    *x52
     4  +coeff( 67)*x11*x22            
     5  +coeff( 68)    *x21    *x43    
     6  +coeff( 69)*x12        *x41    
     7  +coeff( 70)    *x24    *x41    
     8  +coeff( 71)    *x21    *x44    
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 72)        *x31*x43*x51
     1  +coeff( 73)*x11*x24            
     2  +coeff( 74)*x12*x21    *x41    
     3  +coeff( 75)*x11*x23    *x41    
     4  +coeff( 76)    *x21*x34*x41    
     5  +coeff( 77)*x11*x22    *x42    
     6  +coeff( 78)    *x23*x31*x42    
     7  +coeff( 79)    *x23    *x43    
     8  +coeff( 80)    *x24    *x41*x51
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 81)*x12*x21*x31*x42    
     1  +coeff( 82)*x11*x22    *x42*x52
     2  +coeff( 83)    *x21*x31*x43*x53
     3  +coeff( 84)*x12*x24*x31    *x51
     4  +coeff( 85)*x11*x24*x31*x42*x52
     5  +coeff( 86)*x11*x21*x31        
     6  +coeff( 87)*x11    *x31*x41    
     7  +coeff( 88)    *x21*x31*x42    
     8  +coeff( 89)        *x32*x42    
      p_s4_dmd    =p_s4_dmd    
     9  +coeff( 90)    *x23        *x51
c
      return
      end
      function l_s4_dmd    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.8050958E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15713450E-02, 0.18315667E+00,-0.95342956E-02,-0.56483117E-02,
     + -0.44909813E-01, 0.76008891E-02, 0.16764561E-01, 0.33006281E-01,
     + -0.34339346E-02,-0.61855824E-02,-0.21856293E-01,-0.54210164E-02,
     +  0.17791145E-01,-0.38192922E-02, 0.41908980E-02,-0.60852608E-02,
     + -0.10699409E-02,-0.12072694E-02, 0.19397010E-02,-0.29935893E-02,
     +  0.45777387E-02, 0.66788807E-02,-0.12097901E-01, 0.14577692E-01,
     +  0.10202106E-02, 0.54734736E-03, 0.43826676E-02,-0.49476987E-02,
     + -0.22469780E-02,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      l_s4_dmd    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x23    *x41    
      l_s4_dmd    =l_s4_dmd    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x24            
     8  +coeff( 17)        *x31        
      l_s4_dmd    =l_s4_dmd    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11*x21            
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)    *x22*x31        
     8  +coeff( 26)*x11            *x51
      l_s4_dmd    =l_s4_dmd    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)*x11*x23            
c
      return
      end
      function x_s4_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.6067465E+01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.91782666E-03, 0.71866475E-01,-0.16665613E-01, 0.28103194E-02,
     +  0.62461430E-02,-0.20728183E-02, 0.11193485E-01, 0.17403966E-02,
     + -0.21678347E-02,-0.77981274E-02,-0.11781803E-02, 0.46195895E-02,
     + -0.23911572E-02, 0.23808496E-02,-0.80462225E-03, 0.15634656E-02,
     + -0.12399531E-02, 0.65941112E-02, 0.22017581E-02,-0.73962903E-03,
     + -0.39559952E-02,-0.35061853E-02, 0.25488455E-02, 0.18874153E-03,
     +  0.27661299E-03,-0.44369856E-02, 0.39457674E-02, 0.21930856E-02,
     + -0.25153169E-03, 0.88588189E-03, 0.44107600E-02,-0.14899696E-02,
     + -0.16839517E-03, 0.47986594E-03,-0.13909367E-03,-0.73271123E-03,
     + -0.30513052E-03,-0.20698535E-02,-0.41226836E-03,-0.50879759E-03,
     +  0.38968613E-02,-0.50360762E-03, 0.96639368E-03,-0.12620520E-02,
     + -0.52460324E-03,-0.30968009E-03,-0.12890298E-02, 0.29174975E-03,
     +  0.59262728E-02, 0.62481657E-03, 0.23685669E-03,-0.95745397E-03,
     + -0.23544843E-02, 0.37414485E-03,-0.36029916E-03, 0.22745030E-02,
     +  0.33074885E-02, 0.33416448E-02,-0.12986085E-02, 0.23863800E-02,
     + -0.64853020E-02, 0.44677872E-02, 0.47141411E-02, 0.31433406E-02,
     +  0.18917676E-02,-0.54909857E-02, 0.20826139E-01, 0.40556279E-04,
     + -0.29325162E-03, 0.90342619E-04, 0.91293899E-04,-0.27375031E-03,
     + -0.20825134E-03, 0.17671483E-02,-0.38397138E-03, 0.93295304E-04,
     + -0.76448763E-04, 0.48144144E-03, 0.15597159E-03,-0.31540619E-03,
     +  0.16114383E-03, 0.32538941E-03,-0.15273466E-03,-0.89567940E-03,
     + -0.20822756E-03, 0.51680178E-03, 0.23024406E-03,-0.14537146E-03,
     + -0.21402370E-02, 0.10848428E-02,
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
      x_s4_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)    *x23            
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x22            
      x_s4_den    =x_s4_den    
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)            *x41    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)    *x24            
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)*x11        *x41    
      x_s4_den    =x_s4_den    
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)            *x42    
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)    *x24    *x41    
     5  +coeff( 23)    *x23*x31*x41    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x21    *x43*x51
     8  +coeff( 26)    *x24    *x42    
      x_s4_den    =x_s4_den    
     9  +coeff( 27)*x11*x23    *x42*x51
     1  +coeff( 28)*x11*x23*x31*x43    
     2  +coeff( 29)        *x31        
     3  +coeff( 30)    *x22*x31        
     4  +coeff( 31)    *x22    *x42    
     5  +coeff( 32)    *x21*x31*x42    
     6  +coeff( 33)            *x43*x51
     7  +coeff( 34)*x11*x21            
     8  +coeff( 35)*x11    *x31        
      x_s4_den    =x_s4_den    
     9  +coeff( 36)    *x24*x31        
     1  +coeff( 37)    *x22*x32*x41    
     2  +coeff( 38)    *x21    *x44    
     3  +coeff( 39)    *x22*x32    *x51
     4  +coeff( 40)    *x21*x31*x42*x51
     5  +coeff( 41)    *x23    *x43    
     6  +coeff( 42)*x11*x23            
     7  +coeff( 43)    *x22*x32    *x53
     8  +coeff( 44)*x11*x24            
      x_s4_den    =x_s4_den    
     9  +coeff( 45)*x11        *x44    
     1  +coeff( 46)*x11*x22    *x41*x51
     2  +coeff( 47)*x11*x21    *x42*x51
     3  +coeff( 48)*x11        *x43*x51
     4  +coeff( 49)    *x23    *x43*x52
     5  +coeff( 50)    *x23    *x42*x53
     6  +coeff( 51)*x11*x24*x31        
     7  +coeff( 52)*x11*x23    *x42    
     8  +coeff( 53)    *x24*x33*x42    
      x_s4_den    =x_s4_den    
     9  +coeff( 54)*x11*x22    *x43    
     1  +coeff( 55)*x11*x21*x31*x43    
     2  +coeff( 56)    *x24*x31*x42*x52
     3  +coeff( 57)    *x23*x31*x42*x53
     4  +coeff( 58)    *x23*x33*x44    
     5  +coeff( 59)*x11*x22*x31*x42*x52
     6  +coeff( 60)*x11*x22    *x43*x52
     7  +coeff( 61)    *x24*x32*x43*x52
     8  +coeff( 62)    *x22*x32*x43*x54
      x_s4_den    =x_s4_den    
     9  +coeff( 63)*x11*x22    *x44*x52
     1  +coeff( 64)    *x23*x33*x44*x52
     2  +coeff( 65)*x11*x24*x31*x42*x52
     3  +coeff( 66)*x11*x24*x33*x44    
     4  +coeff( 67)*x11*x24*x31*x44*x52
     5  +coeff( 68)                *x51
     6  +coeff( 69)        *x31*x41    
     7  +coeff( 70)            *x41*x51
     8  +coeff( 71)        *x33        
      x_s4_den    =x_s4_den    
     9  +coeff( 72)            *x43    
     1  +coeff( 73)    *x21        *x52
     2  +coeff( 74)    *x22*x31*x41    
     3  +coeff( 75)    *x21*x32*x41    
     4  +coeff( 76)        *x33*x41    
     5  +coeff( 77)    *x23        *x51
     6  +coeff( 78)    *x21    *x42*x51
     7  +coeff( 79)    *x21*x31    *x52
     8  +coeff( 80)    *x21    *x41*x52
      x_s4_den    =x_s4_den    
     9  +coeff( 81)    *x21        *x53
     1  +coeff( 82)    *x23*x32        
     2  +coeff( 83)    *x21*x34        
     3  +coeff( 84)    *x21*x31*x43    
     4  +coeff( 85)        *x32*x42*x51
     5  +coeff( 86)    *x21*x31*x41*x52
     6  +coeff( 87)*x11*x21    *x41    
     7  +coeff( 88)*x11    *x31*x41    
     8  +coeff( 89)    *x24*x31*x41    
      x_s4_den    =x_s4_den    
     9  +coeff( 90)    *x23*x32*x41    
c
      return
      end
      function t_s4_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.3746166E+01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.31861275E-01, 0.56916499E+00,-0.26545517E-01, 0.10378453E+00,
     + -0.16221641E+00, 0.22046197E-01, 0.78277491E-01,-0.25704401E-02,
     +  0.10797761E+00,-0.11334241E-01,-0.21846170E-01,-0.72014697E-01,
     + -0.27735047E-02,-0.10728246E-02,-0.26359979E-01, 0.66493633E-02,
     +  0.22279624E-01,-0.17519746E-02,-0.63576549E-02,-0.11284018E-01,
     +  0.16362049E-01, 0.25669834E-01, 0.53444166E-01, 0.20947054E-02,
     + -0.33441193E-01, 0.47008949E-02, 0.31290825E-01,-0.14805477E-01,
     +  0.29926093E-01, 0.12650945E-01,-0.38393594E-01,-0.63652024E-01,
     + -0.13628684E-02,-0.15157580E-02, 0.63308780E-02,-0.33866125E-02,
     +  0.63637933E-02, 0.10058939E-02, 0.79788994E-02,-0.37417803E-02,
     + -0.23098681E-02, 0.64646220E-03,-0.55694296E-02, 0.16838489E-02,
     + -0.22448131E-02, 0.22366582E-02,-0.10656317E-02,-0.23408964E-01,
     + -0.91639301E-02,-0.12239128E-01,-0.17611731E-01,-0.26225499E-03,
     + -0.23900634E-02,-0.11685275E-01, 0.16926872E-02, 0.85449768E-02,
     +  0.39756256E-02, 0.37766255E-01, 0.20780582E-01,-0.11010682E-01,
     +  0.17391561E-01,-0.21364151E-02, 0.17708153E-01,-0.24615711E-03,
     +  0.36674082E-01,-0.19232854E-01, 0.39325887E-02,-0.91867186E-02,
     +  0.47578271E-02, 0.24467517E-01,-0.50713164E-02,-0.29134631E-01,
     + -0.29129856E-02, 0.78663444E-02,-0.25273846E-01, 0.14608954E-01,
     + -0.19715672E-01,-0.12298083E-01, 0.14541417E-01, 0.17181164E-01,
     +  0.37577309E-01, 0.40460262E-01, 0.12436664E-01, 0.48761502E-01,
     +  0.20549947E-01, 0.46410296E-01,-0.30297387E-01, 0.11203002E-01,
     +  0.19517038E-01,-0.28632049E-01,
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
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
      x54 = x53*x5
c
c                  function
c
      t_s4_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      t_s4_den    =t_s4_den    
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21*x31*x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x22            
      t_s4_den    =t_s4_den    
     9  +coeff( 18)        *x31        
     1  +coeff( 19)            *x42    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x23*x31        
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11            *x51
     7  +coeff( 25)    *x21    *x43    
     8  +coeff( 26)    *x21    *x42*x51
      t_s4_den    =t_s4_den    
     9  +coeff( 27)    *x23*x31*x41    
     1  +coeff( 28)*x11*x24            
     2  +coeff( 29)*x11*x23    *x41*x51
     3  +coeff( 30)*x11*x23*x31*x42    
     4  +coeff( 31)*x11*x24    *x41*x51
     5  +coeff( 32)    *x24*x33*x42    
     6  +coeff( 33)        *x31*x41    
     7  +coeff( 34)*x11    *x31        
     8  +coeff( 35)    *x21    *x41*x51
      t_s4_den    =t_s4_den    
     9  +coeff( 36)    *x21        *x52
     1  +coeff( 37)    *x24            
     2  +coeff( 38)    *x22*x32        
     3  +coeff( 39)*x11*x21    *x41    
     4  +coeff( 40)*x11        *x42    
     5  +coeff( 41)    *x21*x31*x42    
     6  +coeff( 42)*x11        *x41*x51
     7  +coeff( 43)    *x21    *x41*x52
     8  +coeff( 44)*x12*x21            
      t_s4_den    =t_s4_den    
     9  +coeff( 45)*x11*x23            
     1  +coeff( 46)*x11*x22*x31        
     2  +coeff( 47)    *x24    *x41    
     3  +coeff( 48)    *x21    *x44    
     4  +coeff( 49)*x11*x21    *x41*x51
     5  +coeff( 50)    *x23    *x41*x51
     6  +coeff( 51)    *x21*x31*x42*x51
     7  +coeff( 52)            *x44*x51
     8  +coeff( 53)*x12*x22            
      t_s4_den    =t_s4_den    
     9  +coeff( 54)*x11*x23    *x41    
     1  +coeff( 55)*x12        *x42    
     2  +coeff( 56)*x11*x22    *x42    
     3  +coeff( 57)    *x22*x32*x42    
     4  +coeff( 58)    *x23    *x43    
     5  +coeff( 59)    *x22    *x44    
     6  +coeff( 60)    *x21*x31*x44    
     7  +coeff( 61)*x11*x22    *x41*x51
     8  +coeff( 62)*x11    *x32*x41*x51
      t_s4_den    =t_s4_den    
     9  +coeff( 63)    *x23    *x41*x52
     1  +coeff( 64)*x12*x21*x31*x41    
     2  +coeff( 65)    *x22*x33*x42    
     3  +coeff( 66)    *x24    *x43    
     4  +coeff( 67)*x11*x21*x31*x42*x51
     5  +coeff( 68)    *x22    *x44*x51
     6  +coeff( 69)    *x22*x32    *x53
     7  +coeff( 70)    *x21*x31*x42*x53
     8  +coeff( 71)    *x21    *x42*x54
      t_s4_den    =t_s4_den    
     9  +coeff( 72)    *x24*x32*x42    
     1  +coeff( 73)*x11*x23*x31*x41*x51
     2  +coeff( 74)    *x24*x32*x41*x51
     3  +coeff( 75)    *x24    *x42*x52
     4  +coeff( 76)    *x23    *x42*x53
     5  +coeff( 77)*x12*x24    *x41    
     6  +coeff( 78)*x12*x22*x31*x42    
     7  +coeff( 79)*x12*x22    *x43    
     8  +coeff( 80)    *x21*x33*x43*x52
      t_s4_den    =t_s4_den    
     9  +coeff( 81)    *x23*x31*x42*x53
     1  +coeff( 82)    *x23    *x42*x54
     2  +coeff( 83)*x12*x23*x32*x41    
     3  +coeff( 84)*x11*x23*x31*x43*x51
     4  +coeff( 85)*x11*x23    *x44*x51
     5  +coeff( 86)*x11*x24*x31*x41*x52
     6  +coeff( 87)*x11*x22*x33*x41*x52
     7  +coeff( 88)    *x24*x31*x41*x54
     8  +coeff( 89)    *x22    *x44*x54
      t_s4_den    =t_s4_den    
     9  +coeff( 90)*x12*x23*x31*x43    
c
      return
      end
      function y_s4_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 74)
      data ncoeff/ 73/
      data avdat/ -0.1403265E-01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.51905699E-02, 0.51915473E-02, 0.14354537E+00,-0.36918032E-02,
     +  0.27985945E-01, 0.58497018E-02, 0.74135922E-02,-0.32273382E-01,
     +  0.99416887E-02,-0.27627906E-01,-0.39687795E-02, 0.35215219E-05,
     +  0.11528943E-01, 0.11651122E-02,-0.31665217E-02, 0.38420539E-02,
     + -0.30816945E-02, 0.23181532E-02, 0.16075205E-01,-0.35077501E-02,
     + -0.14970704E-02, 0.92793265E-02,-0.11006475E-01, 0.27732521E-02,
     +  0.11075906E-01, 0.25897976E-02, 0.87305598E-04, 0.27696046E-03,
     + -0.29787744E-03, 0.76031330E-03,-0.15785398E-02, 0.56099342E-02,
     +  0.62822462E-02, 0.29745162E-02,-0.40490751E-03, 0.18131694E-03,
     + -0.33813016E-03,-0.12717147E-03,-0.51343866E-03, 0.13465298E-02,
     +  0.84096642E-03, 0.37751973E-02,-0.20642030E-03,-0.20557097E-02,
     +  0.11176587E-02,-0.35021672E-03,-0.31934836E-03,-0.29974975E-03,
     + -0.24932652E-03, 0.41582051E-03,-0.87295244E-04,-0.16932054E-02,
     +  0.64073550E-03, 0.15508028E-02,-0.44061939E-03,-0.26537643E-02,
     +  0.25493203E-03,-0.47991934E-03,-0.31645554E-02, 0.22844321E-02,
     + -0.10178190E-02,-0.95419853E-03,-0.67852827E-03, 0.63650068E-02,
     + -0.12466279E-02,-0.24185267E-02,-0.16039455E-02, 0.24814678E-02,
     +  0.26858349E-02,-0.16572472E-02, 0.28808173E-02, 0.18147902E-02,
     +  0.12701049E-02,
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
      y_s4_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)                *x51
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22            
      y_s4_den    =y_s4_den    
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)        *x33*x41    
     4  +coeff( 13)    *x24            
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)    *x22*x31        
      y_s4_den    =y_s4_den    
     9  +coeff( 18)    *x23            
     1  +coeff( 19)            *x44    
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)        *x31*x45    
     4  +coeff( 22)        *x31*x43    
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x24    *x41    
     8  +coeff( 26)*x11*x23            
      y_s4_den    =y_s4_den    
     9  +coeff( 27)        *x32        
     1  +coeff( 28)        *x31    *x51
     2  +coeff( 29)*x11                
     3  +coeff( 30)        *x32*x41    
     4  +coeff( 31)    *x22*x31*x41    
     5  +coeff( 32)            *x45    
     6  +coeff( 33)    *x24    *x42    
     7  +coeff( 34)*x11*x23    *x41    
     8  +coeff( 35)    *x21*x31        
      y_s4_den    =y_s4_den    
     9  +coeff( 36)                *x52
     1  +coeff( 37)    *x21    *x42    
     2  +coeff( 38)            *x41*x52
     3  +coeff( 39)    *x22        *x51
     4  +coeff( 40)        *x32*x42    
     5  +coeff( 41)*x11*x22            
     6  +coeff( 42)        *x31*x44    
     7  +coeff( 43)        *x32*x42*x51
     8  +coeff( 44)*x11*x21    *x42    
      y_s4_den    =y_s4_den    
     9  +coeff( 45)    *x24*x31        
     1  +coeff( 46)    *x21*x31*x41    
     2  +coeff( 47)            *x42*x51
     3  +coeff( 48)*x11        *x41    
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)    *x23*x31        
     6  +coeff( 51)*x12                
     7  +coeff( 52)    *x21    *x44    
     8  +coeff( 53)    *x23    *x42    
      y_s4_den    =y_s4_den    
     9  +coeff( 54)*x11*x22    *x41    
     1  +coeff( 55)            *x41*x54
     2  +coeff( 56)    *x22    *x44    
     3  +coeff( 57)    *x21*x34*x42    
     4  +coeff( 58)*x11*x24            
     5  +coeff( 59)*x11*x22*x31*x42    
     6  +coeff( 60)*x11*x23    *x42    
     7  +coeff( 61)*x11*x24    *x41    
     8  +coeff( 62)    *x22*x31*x45    
      y_s4_den    =y_s4_den    
     9  +coeff( 63)*x11*x24*x31*x41    
     1  +coeff( 64)*x11*x24*x31*x42    
     2  +coeff( 65)*x11*x24    *x43    
     3  +coeff( 66)    *x22*x31*x45*x52
     4  +coeff( 67)*x12*x21*x31*x44    
     5  +coeff( 68)    *x21    *x45*x54
     6  +coeff( 69)*x11*x24    *x43*x51
     7  +coeff( 70)*x12*x21    *x43*x52
     8  +coeff( 71)*x11*x24*x31*x41*x52
      y_s4_den    =y_s4_den    
     9  +coeff( 72)*x11*x22*x31*x45*x51
     1  +coeff( 73)    *x23*x35*x41*x52
c
      return
      end
      function p_s4_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1037616E-01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.29139279E-02, 0.25948817E-02,-0.51492346E-02,-0.61680462E-01,
     + -0.55274707E-02, 0.62119188E-02,-0.86185504E-02,-0.54583844E-03,
     +  0.96781962E-02, 0.18456483E-02,-0.96234027E-03,-0.37033497E-02,
     + -0.23215236E-02,-0.67923736E-03,-0.22824240E-03, 0.49305981E-03,
     +  0.74984320E-03,-0.17876636E-02, 0.45366943E-03, 0.29959397E-02,
     +  0.25268414E-02,-0.12056913E-02, 0.10763606E-02,-0.60036156E-03,
     +  0.82755642E-03,-0.12454201E-02, 0.26372939E-02,-0.34118507E-02,
     +  0.81811735E-03,-0.31032409E-02, 0.49824180E-03, 0.23940155E-03,
     +  0.33182590E-03,-0.15895908E-02,-0.26994268E-03,-0.11208708E-02,
     + -0.54242898E-03, 0.15564201E-02, 0.76566372E-04,-0.14339766E-03,
     +  0.52274857E-03, 0.41779509E-03,-0.16565139E-03,-0.16738183E-03,
     + -0.23450225E-03,-0.11929744E-03, 0.32923985E-03,-0.31059713E-03,
     + -0.20343412E-02, 0.59665710E-03, 0.49076870E-03,-0.72391448E-03,
     +  0.33302189E-04, 0.34511948E-04,-0.26063135E-03,-0.56379973E-04,
     +  0.57979196E-04,-0.16236483E-03, 0.48642331E-04, 0.69770153E-03,
     +  0.66184584E-04, 0.20985110E-03,-0.24855541E-03,-0.39949251E-03,
     +  0.16051666E-03,-0.58746000E-03,-0.16289280E-03,-0.65069989E-03,
     +  0.26370506E-03,-0.62870595E-03, 0.56224572E-03, 0.24953636E-03,
     + -0.31168649E-04, 0.17286347E-04, 0.77843259E-04, 0.12419086E-04,
     + -0.19922234E-03,-0.11088662E-03,-0.46219368E-03,-0.95500218E-04,
     + -0.14206683E-03,-0.23723362E-03,-0.15038424E-03,-0.22428059E-03,
     + -0.70300681E-03, 0.14469052E-03,-0.24039682E-03, 0.14165351E-03,
     +  0.52694423E-03,-0.84091786E-04,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_s4_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      p_s4_den    =p_s4_den    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)        *x31*x41    
     7  +coeff( 16)        *x31    *x51
     8  +coeff( 17)*x11*x21            
      p_s4_den    =p_s4_den    
     9  +coeff( 18)    *x24            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x22    *x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)        *x31*x42    
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)            *x41*x52
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)    *x23    *x41    
      p_s4_den    =p_s4_den    
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)            *x44    
     2  +coeff( 29)            *x43*x51
     3  +coeff( 30)    *x24    *x41    
     4  +coeff( 31)    *x21*x31*x41    
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)        *x31*x43    
     8  +coeff( 35)*x11*x21        *x51
      p_s4_den    =p_s4_den    
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)            *x44*x51
     3  +coeff( 39)        *x31*x41*x53
     4  +coeff( 40)        *x34*x42    
     5  +coeff( 41)    *x23            
     6  +coeff( 42)    *x22*x31        
     7  +coeff( 43)        *x32*x41    
     8  +coeff( 44)                *x53
      p_s4_den    =p_s4_den    
     9  +coeff( 45)*x11*x22            
     1  +coeff( 46)*x11        *x41*x51
     2  +coeff( 47)        *x31*x42*x51
     3  +coeff( 48)*x11*x22    *x41    
     4  +coeff( 49)    *x23    *x42    
     5  +coeff( 50)    *x24        *x51
     6  +coeff( 51)        *x31*x43*x51
     7  +coeff( 52)*x11*x23    *x41    
     8  +coeff( 53)*x11    *x31        
      p_s4_den    =p_s4_den    
     9  +coeff( 54)        *x32    *x51
     1  +coeff( 55)            *x42*x51
     2  +coeff( 56)        *x31    *x52
     3  +coeff( 57)*x11*x21*x31        
     4  +coeff( 58)    *x23*x31        
     5  +coeff( 59)*x11        *x42    
     6  +coeff( 60)    *x21    *x43    
     7  +coeff( 61)        *x32*x41*x51
     8  +coeff( 62)    *x22        *x52
      p_s4_den    =p_s4_den    
     9  +coeff( 63)    *x24*x31        
     1  +coeff( 64)    *x23*x31*x41    
     2  +coeff( 65)*x11*x21    *x42    
     3  +coeff( 66)    *x22    *x43    
     4  +coeff( 67)*x11*x21    *x41*x51
     5  +coeff( 68)    *x22    *x42*x51
     6  +coeff( 69)*x11*x24            
     7  +coeff( 70)    *x24    *x42    
     8  +coeff( 71)    *x24    *x41*x51
      p_s4_den    =p_s4_den    
     9  +coeff( 72)*x12*x21*x31*x42    
     1  +coeff( 73)        *x32        
     2  +coeff( 74)*x11            *x51
     3  +coeff( 75)    *x21*x31    *x51
     4  +coeff( 76)*x12                
     5  +coeff( 77)    *x23        *x51
     6  +coeff( 78)    *x22*x31    *x51
     7  +coeff( 79)    *x21    *x42*x51
     8  +coeff( 80)    *x21    *x41*x52
      p_s4_den    =p_s4_den    
     9  +coeff( 81)            *x42*x52
     1  +coeff( 82)        *x31*x44    
     2  +coeff( 83)    *x21*x31*x42*x51
     3  +coeff( 84)*x11*x22    *x42    
     4  +coeff( 85)    *x23    *x43    
     5  +coeff( 86)*x11        *x44    
     6  +coeff( 87)        *x32*x44    
     7  +coeff( 88)*x11*x23        *x51
     8  +coeff( 89)    *x23    *x42*x51
      p_s4_den    =p_s4_den    
     9  +coeff( 90)*x11    *x31*x42*x51
c
      return
      end
      function l_s4_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.1384550E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.28256120E-02, 0.25745868E-03,-0.68531028E-03,-0.58807982E-02,
     +  0.48783819E-04,-0.50428319E-02,-0.44958625E-03,-0.42160647E-02,
     + -0.14104067E-02, 0.23401075E-03, 0.39739460E-02,-0.85966906E-03,
     + -0.87464380E-03,-0.14582823E-03, 0.17844724E-03, 0.15685629E-02,
     + -0.22224446E-02, 0.23535251E-03,-0.49308699E-04,-0.15841769E-03,
     + -0.29489840E-03, 0.14563510E-03, 0.30110707E-03,-0.61037875E-03,
     +  0.27361530E-03,-0.14044087E-04, 0.31103998E-04, 0.90275928E-04,
     + -0.22243477E-03,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_s4_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_s4_den    =l_s4_den    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x24            
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x22*x31        
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)    *x24    *x41    
      l_s4_den    =l_s4_den    
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)        *x31    *x51
     2  +coeff( 20)    *x23            
     3  +coeff( 21)        *x31*x42    
     4  +coeff( 22)            *x41*x52
     5  +coeff( 23)    *x22*x31*x41    
     6  +coeff( 24)            *x44    
     7  +coeff( 25)            *x44*x51
     8  +coeff( 26)*x11                
      l_s4_den    =l_s4_den    
     9  +coeff( 27)    *x21        *x51
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)    *x23    *x41    
c
      return
      end
      function x_s4_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.6826038E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.13915075E-02, 0.11959220E+00,-0.14303363E-01, 0.56590619E-02,
     +  0.23021754E-02, 0.23261392E-02, 0.97869663E-02, 0.15816615E-02,
     + -0.19065118E-02,-0.75259712E-02,-0.95433008E-03, 0.33634303E-02,
     + -0.21489454E-02, 0.19103440E-02,-0.77896641E-03, 0.13503793E-02,
     + -0.11281620E-02, 0.71862461E-02, 0.32539652E-02,-0.67847926E-03,
     +  0.23266920E-03,-0.24035790E-02,-0.27457674E-02, 0.22197554E-02,
     +  0.11799701E-03,-0.39165337E-02,-0.71633421E-03, 0.20879102E-02,
     +  0.83452919E-02,-0.16189588E-03, 0.56959479E-03, 0.38968001E-02,
     + -0.51351689E-03, 0.47077534E-04,-0.23258362E-03, 0.39078889E-03,
     + -0.10354992E-03,-0.26648634E-03,-0.22264055E-03,-0.12639441E-03,
     + -0.24270837E-03, 0.46559167E-03,-0.71674323E-04, 0.16007416E-02,
     +  0.28125738E-04, 0.12594467E-04,-0.38380173E-03,-0.42761844E-02,
     +  0.94240700E-03,-0.60939434E-03, 0.70866774E-03, 0.55203139E-03,
     +  0.30260456E-02, 0.20053282E-02,-0.17576780E-05,-0.84138475E-03,
     + -0.42934209E-03,-0.18341235E-03, 0.19088358E-03, 0.99397403E-04,
     + -0.19133806E-02,-0.71054453E-03, 0.19812630E-02,-0.21239503E-02,
     +  0.48789801E-03, 0.84925303E-03, 0.13631805E-02,-0.57315445E-02,
     + -0.85557008E-03, 0.29566016E-02,-0.66367636E-03,-0.86594670E-03,
     +  0.13656694E-02, 0.24434736E-01,-0.17430581E-01,-0.13221762E-03,
     + -0.40639780E-03,-0.22431101E-03, 0.51781033E-04, 0.40474979E-03,
     +  0.16210468E-03, 0.11141059E-03, 0.14081306E-03,-0.73996891E-03,
     +  0.15143264E-02, 0.10516907E-02,-0.25048494E-03,-0.39124759E-03,
     +  0.12826927E-03, 0.61167270E-03,
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
      x_s4_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)    *x22            
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff(  9)    *x21*x31        
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)            *x41    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)    *x24            
     7  +coeff( 16)    *x23*x31        
     8  +coeff( 17)*x11        *x41    
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)            *x42    
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x21    *x43    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)*x11            *x51
     8  +coeff( 26)    *x24    *x42    
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 27)    *x21*x31*x44    
     1  +coeff( 28)*x11*x24*x31*x41*x52
     2  +coeff( 29)*x11*x22    *x44*x52
     3  +coeff( 30)        *x31        
     4  +coeff( 31)    *x22*x31        
     5  +coeff( 32)    *x22    *x42    
     6  +coeff( 33)    *x21*x31*x42    
     7  +coeff( 34)        *x33    *x51
     8  +coeff( 35)        *x31*x42*x51
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 36)*x11*x21            
     1  +coeff( 37)*x11    *x31        
     2  +coeff( 38)    *x24*x31        
     3  +coeff( 39)    *x22*x32*x41    
     4  +coeff( 40)*x11    *x31*x41    
     5  +coeff( 41)*x11        *x42    
     6  +coeff( 42)        *x31*x44*x51
     7  +coeff( 43)*x11            *x52
     8  +coeff( 44)    *x23    *x41*x52
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 45)    *x22*x31*x41*x52
     1  +coeff( 46)    *x21    *x42*x53
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)    *x24*x31*x42    
     4  +coeff( 49)    *x22*x31*x44    
     5  +coeff( 50)*x11*x21    *x41*x51
     6  +coeff( 51)    *x23*x31*x42*x51
     7  +coeff( 52)    *x21*x33*x41*x52
     8  +coeff( 53)    *x23    *x42*x52
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 54)    *x22*x31*x42*x52
     1  +coeff( 55)    *x23        *x54
     2  +coeff( 56)*x11*x24            
     3  +coeff( 57)*x11*x22*x31*x41    
     4  +coeff( 58)*x11*x21*x31*x42    
     5  +coeff( 59)*x11*x23        *x51
     6  +coeff( 60)*x11*x24*x31        
     7  +coeff( 61)*x11*x24    *x41    
     8  +coeff( 62)*x11*x23    *x42    
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 63)*x11*x23    *x41*x51
     1  +coeff( 64)    *x22*x33*x42*x52
     2  +coeff( 65)*x11*x21*x34*x41    
     3  +coeff( 66)*x11*x23*x31*x42    
     4  +coeff( 67)    *x23*x33*x44    
     5  +coeff( 68)    *x24*x32*x42*x52
     6  +coeff( 69)*x11        *x44*x52
     7  +coeff( 70)    *x22*x32*x42*x54
     8  +coeff( 71)*x11*x24*x32    *x51
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 72)*x11*x22*x31*x42*x52
     1  +coeff( 73)*x11*x23*x31*x41*x53
     2  +coeff( 74)*x11*x24*x31*x44*x52
     3  +coeff( 75)*x11*x24*x33*x44*x52
     4  +coeff( 76)        *x31*x41    
     5  +coeff( 77)            *x43    
     6  +coeff( 78)    *x21        *x52
     7  +coeff( 79)            *x41*x52
     8  +coeff( 80)    *x22*x31*x41    
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 81)    *x21*x32*x41    
     1  +coeff( 82)    *x22    *x41*x51
     2  +coeff( 83)    *x21*x31    *x52
     3  +coeff( 84)    *x21    *x41*x52
     4  +coeff( 85)    *x22*x31*x42    
     5  +coeff( 86)    *x22    *x43    
     6  +coeff( 87)        *x31*x44    
     7  +coeff( 88)    *x22*x31*x41*x51
     8  +coeff( 89)        *x33*x41*x51
      x_s4_q1ex   =x_s4_q1ex   
     9  +coeff( 90)    *x21    *x43*x51
c
      return
      end
      function t_s4_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2436897E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10891698E-04,-0.54013273E-02,-0.21767995E-03,-0.17209012E-02,
     +  0.37979687E-03,-0.34086583E-02, 0.17336833E-02, 0.11467919E-02,
     +  0.89035463E-03,-0.13251943E-02, 0.37855926E-03, 0.28353825E-03,
     +  0.23638925E-02,-0.41561422E-03,-0.24306415E-03,-0.53790520E-03,
     +  0.68551331E-03, 0.64788313E-04,-0.20522303E-03,-0.70340390E-03,
     +  0.87825628E-03, 0.87610638E-03,-0.17060642E-04,-0.12508513E-03,
     + -0.48851136E-04,-0.87667155E-04,-0.87667584E-04, 0.68826182E-03,
     + -0.70835126E-03, 0.60113704E-04, 0.25159563E-03, 0.48846216E-03,
     + -0.17906133E-03, 0.37252254E-03, 0.30024821E-03,-0.27082031E-03,
     +  0.93178480E-03, 0.13382408E-03,-0.13319590E-02, 0.84791944E-04,
     + -0.29315215E-04,-0.73243624E-04, 0.41456664E-04,-0.94333904E-04,
     + -0.32825366E-03, 0.40658429E-04,-0.97626107E-04, 0.45387271E-04,
     +  0.81594240E-04,-0.43383774E-04,-0.51472004E-03,-0.38649603E-04,
     +  0.37252474E-04,-0.10404320E-03,-0.44519456E-04, 0.24713893E-03,
     + -0.48323467E-03,-0.42969983E-04, 0.23304460E-03,-0.70931361E-03,
     + -0.46271427E-04, 0.21526520E-03, 0.93019371E-04,-0.26250337E-03,
     +  0.88028813E-04, 0.18078972E-03, 0.44348359E-04, 0.19444774E-03,
     + -0.31560706E-03,-0.32051938E-03,-0.13434939E-03,-0.44388285E-04,
     + -0.49764139E-03, 0.18346377E-03,-0.33114393E-04, 0.30579683E-03,
     + -0.16045049E-03,-0.10092080E-03, 0.23115649E-03,-0.31263055E-03,
     + -0.18666114E-03, 0.10060283E-03, 0.64963073E-03, 0.27055919E-03,
     +  0.16155076E-03, 0.27738485E-03, 0.33495444E-03, 0.10154086E-03,
     + -0.28046384E-03,-0.80745311E-04,
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
      t_s4_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)*x11*x22            
     3  +coeff( 12)    *x23*x31        
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)*x11*x22    *x41    
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)    *x24            
     2  +coeff( 20)    *x24    *x41    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x23    *x43    
     5  +coeff( 23)                *x51
     6  +coeff( 24)            *x42    
     7  +coeff( 25)        *x33        
     8  +coeff( 26)    *x21    *x41*x51
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 27)    *x21        *x52
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x23        *x51
     4  +coeff( 31)    *x22*x33        
     5  +coeff( 32)    *x23*x31*x41    
     6  +coeff( 33)*x11*x24            
     7  +coeff( 34)    *x23*x31*x42    
     8  +coeff( 35)*x11*x23*x31*x42    
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 36)*x11*x24    *x41*x51
     1  +coeff( 37)*x11*x23    *x42*x51
     2  +coeff( 38)*x12    *x33*x42*x51
     3  +coeff( 39)    *x23*x31*x43*x54
     4  +coeff( 40)*x11*x21            
     5  +coeff( 41)*x11    *x31        
     6  +coeff( 42)            *x43    
     7  +coeff( 43)            *x42*x51
     8  +coeff( 44)*x11        *x42    
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 45)    *x21*x31*x42    
     1  +coeff( 46)*x12*x21            
     2  +coeff( 47)*x11*x23            
     3  +coeff( 48)*x11*x22*x31        
     4  +coeff( 49)    *x24*x31        
     5  +coeff( 50)    *x22*x32*x41    
     6  +coeff( 51)    *x21    *x44    
     7  +coeff( 52)    *x22*x32    *x51
     8  +coeff( 53)*x11        *x42*x51
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 54)            *x44*x51
     1  +coeff( 55)    *x22*x31    *x52
     2  +coeff( 56)    *x21*x31*x41*x52
     3  +coeff( 57)    *x21    *x42*x52
     4  +coeff( 58)*x12*x22            
     5  +coeff( 59)*x11*x22    *x42    
     6  +coeff( 60)    *x24    *x42    
     7  +coeff( 61)*x11    *x31*x43    
     8  +coeff( 62)    *x22    *x44    
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 63)*x11*x21*x31*x41*x51
     1  +coeff( 64)*x11*x21    *x42*x51
     2  +coeff( 65)*x11*x21    *x41*x52
     3  +coeff( 66)    *x23    *x41*x52
     4  +coeff( 67)    *x23        *x53
     5  +coeff( 68)    *x21    *x42*x53
     6  +coeff( 69)    *x24*x33        
     7  +coeff( 70)*x11*x24    *x41    
     8  +coeff( 71)*x11*x23    *x42    
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 72)*x12    *x31*x42    
     1  +coeff( 73)    *x24*x31*x42    
     2  +coeff( 74)    *x22*x33*x42    
     3  +coeff( 75)*x11*x24        *x51
     4  +coeff( 76)*x11*x23    *x41*x51
     5  +coeff( 77)    *x24*x31*x41*x51
     6  +coeff( 78)*x11*x21*x32*x41*x51
     7  +coeff( 79)    *x23*x31*x42*x51
     8  +coeff( 80)    *x21*x33*x42*x51
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 81)*x11*x21    *x43*x51
     1  +coeff( 82)*x11*x22    *x41*x52
     2  +coeff( 83)    *x23    *x42*x52
     3  +coeff( 84)    *x22*x31*x42*x52
     4  +coeff( 85)    *x22*x32    *x53
     5  +coeff( 86)    *x21*x31*x42*x53
     6  +coeff( 87)    *x21    *x42*x54
     7  +coeff( 88)*x11*x21*x34*x41    
     8  +coeff( 89)    *x24*x32*x42    
      t_s4_q1ex   =t_s4_q1ex   
     9  +coeff( 90)*x11*x24*x31    *x51
c
      return
      end
      function y_s4_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 60)
      data ncoeff/ 59/
      data avdat/ -0.1752748E-01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.57854517E-02, 0.13665304E+00,-0.26707994E-02, 0.19575985E-01,
     +  0.10535419E-02, 0.37507145E-02,-0.22830626E-02,-0.13948121E-03,
     + -0.22728911E-01, 0.21324008E-03, 0.31420547E-02, 0.83063003E-02,
     + -0.19270107E-01, 0.48592672E-04,-0.29309974E-02, 0.12027493E-01,
     +  0.80835931E-02,-0.22177689E-03, 0.78095277E-02,-0.16511063E-02,
     + -0.10219439E-02,-0.22202262E-02, 0.17787585E-02, 0.56572100E-02,
     + -0.75077643E-02,-0.26643558E-02, 0.88139120E-02, 0.20764491E-02,
     +  0.74965115E-04, 0.62373321E-03, 0.49311868E-04, 0.84633171E-03,
     + -0.13434393E-02, 0.19567993E-02,-0.33021939E-03,-0.20016407E-03,
     +  0.79067591E-04, 0.98261086E-03, 0.38380179E-03, 0.35807050E-02,
     +  0.88570919E-03, 0.32141036E-02, 0.25941522E-02, 0.91638445E-03,
     +  0.11058980E-03,-0.89567489E-04,-0.15339669E-03,-0.15228291E-03,
     + -0.17152041E-03, 0.33529568E-03, 0.10268719E-03, 0.23454556E-02,
     +  0.17419810E-03,-0.12114352E-02,-0.95213024E-03,-0.14381802E-03,
     + -0.71302650E-03,-0.14425882E-02,-0.80224930E-03,
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
c
c                  function
c
      y_s4_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)    *x21            
     4  +coeff(  4)                *x51
     5  +coeff(  5)        *x31*x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)        *x31    *x51
      y_s4_q1ex   =y_s4_q1ex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x33        
     2  +coeff( 11)        *x31*x42    
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)        *x31    *x52
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)            *x44    
     8  +coeff( 17)    *x24            
      y_s4_q1ex   =y_s4_q1ex   
     9  +coeff( 18)        *x35        
     1  +coeff( 19)        *x31        
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)                *x52
     4  +coeff( 22)    *x22*x31        
     5  +coeff( 23)    *x23            
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)*x11*x21    *x41    
      y_s4_q1ex   =y_s4_q1ex   
     9  +coeff( 27)    *x24    *x41    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)        *x32        
     3  +coeff( 30)        *x32*x41    
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)    *x22        *x51
     6  +coeff( 33)    *x22*x31*x41    
     7  +coeff( 34)    *x23    *x41    
     8  +coeff( 35)    *x21*x31        
      y_s4_q1ex   =y_s4_q1ex   
     9  +coeff( 36)*x11                
     1  +coeff( 37)    *x21    *x42    
     2  +coeff( 38)        *x32*x42    
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)            *x45    
     5  +coeff( 41)    *x24*x31        
     6  +coeff( 42)    *x24    *x42    
     7  +coeff( 43)*x11*x23    *x41    
     8  +coeff( 44)    *x22    *x41*x53
      y_s4_q1ex   =y_s4_q1ex   
     9  +coeff( 45)    *x21        *x51
     1  +coeff( 46)    *x21*x31*x41    
     2  +coeff( 47)        *x31*x41*x51
     3  +coeff( 48)        *x32*x41*x51
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)    *x23*x31        
     6  +coeff( 51)*x11        *x41*x51
     7  +coeff( 52)        *x31*x44    
     8  +coeff( 53)*x11*x21        *x51
      y_s4_q1ex   =y_s4_q1ex   
     9  +coeff( 54)    *x21    *x44    
     1  +coeff( 55)            *x44*x51
     2  +coeff( 56)*x11        *x43    
     3  +coeff( 57)*x11*x21    *x42    
     4  +coeff( 58)    *x22    *x44    
     5  +coeff( 59)    *x23*x31*x43    
c
      return
      end
      function p_s4_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5518497E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18198612E-02,-0.12072953E-02, 0.21354002E-02, 0.50981138E-01,
     +  0.97710257E-02,-0.10812701E-01, 0.18164448E-02,-0.88155707E-02,
     +  0.52591838E-03,-0.14254043E-02,-0.13269585E-02, 0.52378452E-02,
     +  0.38170372E-02,-0.11539401E-02,-0.60828112E-03, 0.73655846E-03,
     + -0.11201821E-02, 0.16132960E-02,-0.12527466E-02,-0.41291844E-02,
     +  0.53448086E-02, 0.39533852E-02, 0.46661883E-03,-0.68123796E-03,
     +  0.25241924E-02, 0.89456735E-03, 0.20943831E-03,-0.10120485E-03,
     + -0.10718891E-03, 0.26547862E-03, 0.57664415E-03, 0.10024967E-02,
     +  0.11511592E-02,-0.14895377E-03, 0.64113847E-04,-0.27087779E-03,
     + -0.33807987E-05, 0.33652113E-03, 0.85935746E-04, 0.33856367E-03,
     +  0.51158067E-03,-0.62017824E-03, 0.42204544E-03,-0.15999884E-02,
     +  0.16691608E-02, 0.36616795E-04,-0.29912288E-04, 0.41199371E-04,
     + -0.15594493E-03, 0.13097234E-03, 0.42911190E-04, 0.15652095E-03,
     + -0.57149102E-03,-0.41610023E-03, 0.61790802E-03,-0.11218817E-03,
     + -0.26688233E-03, 0.45952524E-03,-0.48521374E-05, 0.51720178E-03,
     +  0.21824666E-03,-0.25815296E-03, 0.10423301E-02, 0.32439639E-03,
     + -0.28372469E-03,-0.10748227E-02,-0.70027076E-04, 0.24021425E-04,
     +  0.24824360E-03,-0.69414775E-04, 0.32440765E-03,-0.39511728E-04,
     + -0.43373511E-04,-0.45733821E-04, 0.10023937E-03,-0.43036154E-04,
     +  0.18406937E-03, 0.23276715E-03,-0.13803999E-03,-0.17351916E-03,
     +  0.95895142E-04,-0.60245808E-03,-0.16341166E-03, 0.10190275E-03,
     + -0.24556177E-03, 0.16185101E-03, 0.52780171E-04,-0.41547406E-03,
     + -0.44219963E-04, 0.39970691E-03,
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
      p_s4_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x22    *x41    
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x41*x51
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)            *x43    
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x23            
     8  +coeff( 17)    *x22*x31        
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 18)        *x31*x42    
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)            *x44    
     4  +coeff( 22)    *x24    *x41    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x22*x31*x41    
     7  +coeff( 25)        *x31*x43    
     8  +coeff( 26)*x11*x23            
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 27)        *x34*x42    
     1  +coeff( 28)*x11                
     2  +coeff( 29)        *x31    *x51
     3  +coeff( 30)        *x32*x41    
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)    *x23    *x41    
     6  +coeff( 33)*x11*x23    *x41    
     7  +coeff( 34)    *x21*x31        
     8  +coeff( 35)    *x21        *x51
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 36)    *x21    *x42    
     1  +coeff( 37)        *x31*x41*x51
     2  +coeff( 38)*x11*x22            
     3  +coeff( 39)*x11*x21        *x51
     4  +coeff( 40)    *x22    *x41*x51
     5  +coeff( 41)    *x24*x31        
     6  +coeff( 42)*x11*x21    *x42    
     7  +coeff( 43)    *x23    *x42    
     8  +coeff( 44)            *x44*x51
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 45)    *x24    *x42    
     1  +coeff( 46)        *x32        
     2  +coeff( 47)    *x21*x31*x41    
     3  +coeff( 48)                *x53
     4  +coeff( 49)*x11*x21*x31        
     5  +coeff( 50)    *x23*x31        
     6  +coeff( 51)*x11        *x41*x51
     7  +coeff( 52)            *x42*x52
     8  +coeff( 53)    *x22    *x43    
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 54)    *x21    *x44    
     1  +coeff( 55)        *x31*x44    
     2  +coeff( 56)    *x24        *x51
     3  +coeff( 57)        *x31*x43*x51
     4  +coeff( 58)        *x32*x44    
     5  +coeff( 59)*x12            *x52
     6  +coeff( 60)*x11*x23    *x42    
     7  +coeff( 61)*x12*x21    *x41*x51
     8  +coeff( 62)    *x21*x32*x43*x51
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 63)            *x44*x53
     1  +coeff( 64)*x11    *x33*x42*x51
     2  +coeff( 65)*x12*x21    *x41*x52
     3  +coeff( 66)*x12*x24*x31*x44    
     4  +coeff( 67)*x11        *x41    
     5  +coeff( 68)        *x31    *x52
     6  +coeff( 69)            *x41*x52
     7  +coeff( 70)*x12                
     8  +coeff( 71)    *x21    *x43    
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 72)*x11    *x31    *x51
     1  +coeff( 73)    *x21*x31*x41*x51
     2  +coeff( 74)        *x32*x41*x51
     3  +coeff( 75)    *x21    *x41*x52
     4  +coeff( 76)*x12        *x41    
     5  +coeff( 77)*x11*x22    *x41    
     6  +coeff( 78)    *x22*x31*x42    
     7  +coeff( 79)    *x21*x31*x43    
     8  +coeff( 80)            *x43*x52
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 81)    *x22        *x53
     1  +coeff( 82)            *x42*x53
     2  +coeff( 83)            *x41*x54
     3  +coeff( 84)*x12*x22            
     4  +coeff( 85)*x11*x24            
     5  +coeff( 86)*x11*x23*x31        
     6  +coeff( 87)*x12    *x32        
     7  +coeff( 88)    *x23    *x43    
     8  +coeff( 89)*x11        *x44    
      p_s4_q1ex   =p_s4_q1ex   
     9  +coeff( 90)    *x22    *x44    
c
      return
      end
      function l_s4_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.9136269E-04/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26471E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15655009E-02, 0.10261835E-03,-0.73652988E-03,-0.64483904E-02,
     + -0.52198484E-04,-0.19881302E-02,-0.15309962E-03,-0.23604480E-02,
     + -0.74680458E-03, 0.16201988E-02,-0.60801016E-03,-0.27225754E-03,
     +  0.72656461E-03, 0.93399503E-04,-0.76837299E-04,-0.18415334E-03,
     +  0.10885035E-03, 0.18741662E-03,-0.46284869E-03,-0.91029040E-03,
     + -0.25823425E-04, 0.32433472E-04,-0.82173043E-04, 0.61969629E-04,
     +  0.10532138E-03,-0.16280863E-03,-0.21451332E-03, 0.15988297E-04,
     +  0.51143041E-04,
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
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_s4_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)        *x31*x41    
     8  +coeff(  8)            *x42    
      l_s4_q1ex   =l_s4_q1ex   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)                *x52
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)    *x22        *x51
      l_s4_q1ex   =l_s4_q1ex   
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)            *x44    
     2  +coeff( 20)    *x24    *x41    
     3  +coeff( 21)        *x31    *x51
     4  +coeff( 22)*x11*x21            
     5  +coeff( 23)    *x23            
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)        *x31*x43    
      l_s4_q1ex   =l_s4_q1ex   
     9  +coeff( 27)*x11*x23    *x41    
     1  +coeff( 28)    *x21        *x51
     2  +coeff( 29)            *x41*x52
c
      return
      end
      function x_s4_q1en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.6654511E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.64385845E-03, 0.10437776E+00,-0.61345608E-02, 0.45105079E-02,
     +  0.27945044E-02, 0.40881201E-02, 0.80403598E-03,-0.74192532E-03,
     + -0.29111353E-02,-0.45623124E-03, 0.21199112E-03, 0.16188340E-02,
     + -0.10147927E-02, 0.27655661E-02, 0.94726682E-03, 0.25936362E-03,
     + -0.43565029E-03, 0.48864045E-03,-0.12578227E-02,-0.51355449E-03,
     + -0.13233336E-02, 0.14559763E-02,-0.29016301E-03, 0.13386711E-02,
     +  0.10538530E-02,-0.79450067E-04, 0.38830293E-03,-0.50910929E-03,
     +  0.28132228E-03, 0.19115498E-03,-0.33649319E-03,-0.52732852E-03,
     +  0.12214492E-03,-0.11472726E-02, 0.11451374E-02, 0.18493166E-03,
     + -0.20614348E-03, 0.74215136E-05,-0.67958790E-04,-0.45847063E-03,
     +  0.41712882E-03,-0.32504118E-03,-0.20180049E-03,-0.27880944E-04,
     + -0.79300837E-03, 0.22787614E-03,-0.17135178E-02, 0.60717924E-03,
     +  0.56107907E-03, 0.45427025E-03, 0.39968826E-03, 0.42408047E-03,
     +  0.33120552E-03,-0.67190267E-04, 0.42749057E-03,-0.15475113E-02,
     +  0.10364270E-02, 0.40231892E-02,-0.12391590E-02,-0.48173158E-03,
     +  0.19345838E-04,-0.44271834E-04,-0.88579845E-05,-0.17855281E-05,
     + -0.13890630E-03,-0.23852393E-04, 0.31893684E-04,-0.17881284E-04,
     + -0.28952851E-04, 0.34838540E-04,-0.40883289E-04, 0.13941841E-03,
     + -0.47138623E-04,-0.19570712E-03, 0.40988667E-04, 0.25991234E-03,
     + -0.65169545E-04,-0.16290722E-03, 0.12720680E-03,-0.90560825E-04,
     +  0.26208608E-03, 0.11689476E-03,-0.26723255E-04,-0.59765414E-04,
     +  0.31551328E-04,-0.36781883E-04, 0.49099949E-03,-0.26027614E-04,
     +  0.86800617E-04,-0.23997398E-03,
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
      x33 = x32*x3
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      x_s4_q1en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x23    *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21*x31        
      x_s4_q1en   =x_s4_q1en   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x21*x31*x41    
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x24            
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 18)    *x23*x31        
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)    *x24    *x41    
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)            *x42    
     6  +coeff( 24)    *x22    *x42    
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)        *x31        
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 27)    *x22*x31        
     1  +coeff( 28)    *x21*x31*x42    
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)*x11*x21            
     4  +coeff( 31)    *x24*x31        
     5  +coeff( 32)    *x21    *x44    
     6  +coeff( 33)    *x21*x31*x42*x51
     7  +coeff( 34)    *x24    *x42    
     8  +coeff( 35)    *x23    *x43    
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 36)    *x22*x31*x41*x52
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11    *x33        
     3  +coeff( 39)*x11*x21    *x42    
     4  +coeff( 40)*x11*x24            
     5  +coeff( 41)*x11*x22    *x42    
     6  +coeff( 42)*x11    *x32*x42    
     7  +coeff( 43)*x11        *x44    
     8  +coeff( 44)*x11*x24*x31        
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 45)*x11*x24    *x41    
     1  +coeff( 46)*x11*x23    *x41*x51
     2  +coeff( 47)    *x23*x31*x43*x52
     3  +coeff( 48)    *x21*x33*x43*x52
     4  +coeff( 49)    *x23    *x44*x52
     5  +coeff( 50)    *x22*x31*x44*x52
     6  +coeff( 51)*x11*x23*x31*x42    
     7  +coeff( 52)*x11    *x32*x44    
     8  +coeff( 53)*x11*x23    *x42*x51
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 54)*x11*x23*x33*x41    
     1  +coeff( 55)*x11*x23*x31*x43    
     2  +coeff( 56)*x11*x22*x31*x42*x52
     3  +coeff( 57)*x11*x24*x31*x41*x52
     4  +coeff( 58)*x11*x24*x31*x42*x52
     5  +coeff( 59)    *x24*x34*x43*x52
     6  +coeff( 60)*x11*x24*x33*x43    
     7  +coeff( 61)                *x51
     8  +coeff( 62)        *x31*x41    
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 63)        *x31    *x51
     1  +coeff( 64)        *x31*x42    
     2  +coeff( 65)            *x43    
     3  +coeff( 66)    *x22        *x51
     4  +coeff( 67)            *x42*x51
     5  +coeff( 68)            *x41*x52
     6  +coeff( 69)    *x23        *x51
     7  +coeff( 70)        *x33    *x51
     8  +coeff( 71)    *x22    *x41*x51
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 72)    *x21    *x41*x52
     1  +coeff( 73)*x11    *x31        
     2  +coeff( 74)    *x22*x31*x42    
     3  +coeff( 75)    *x21*x32*x42    
     4  +coeff( 76)    *x22    *x43    
     5  +coeff( 77)    *x21*x31*x43    
     6  +coeff( 78)    *x23    *x41*x51
     7  +coeff( 79)    *x21    *x43*x51
     8  +coeff( 80)            *x44*x51
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 81)    *x21*x31*x41*x52
     1  +coeff( 82)            *x43*x52
     2  +coeff( 83)*x11*x21*x31        
     3  +coeff( 84)    *x24*x32        
     4  +coeff( 85)*x11*x21    *x41    
     5  +coeff( 86)*x11    *x31*x41    
     6  +coeff( 87)    *x23*x31*x42    
     7  +coeff( 88)*x11    *x31    *x51
     8  +coeff( 89)    *x23*x31*x41*x51
      x_s4_q1en   =x_s4_q1en   
     9  +coeff( 90)    *x23    *x42*x51
c
      return
      end
      function t_s4_q1en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5212933E-03/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23470957E-03, 0.30876638E-01,-0.64637461E-02, 0.24554217E-02,
     +  0.41301986E-02,-0.50936994E-03,-0.27124316E-02,-0.48224152E-04,
     + -0.43327131E-03, 0.78657817E-03, 0.29884526E-03, 0.16083079E-02,
     +  0.79849683E-03,-0.77388278E-03,-0.48899592E-03,-0.86786639E-03,
     +  0.32616736E-03,-0.47723638E-03,-0.14948609E-02, 0.13707876E-02,
     +  0.25630817E-02, 0.40322906E-04,-0.12786796E-02, 0.81479579E-03,
     + -0.70634866E-04,-0.27921426E-03, 0.15270947E-03, 0.33347585E-03,
     +  0.54016296E-03, 0.13692896E-02,-0.50871947E-03, 0.61164001E-05,
     + -0.16891969E-03,-0.37382200E-03,-0.35992014E-03, 0.44858220E-03,
     + -0.11711584E-02, 0.13116407E-02, 0.74868556E-04,-0.64850331E-03,
     + -0.81319828E-04, 0.19711972E-03, 0.10013370E-02, 0.13939926E-02,
     + -0.10155708E-02,-0.61333226E-03,-0.89132256E-03, 0.78744983E-03,
     +  0.87448733E-03, 0.15077284E-04,-0.36200774E-04,-0.63329862E-04,
     + -0.15635873E-03, 0.30280171E-04,-0.20610548E-03,-0.32865701E-04,
     +  0.14948489E-04, 0.72015231E-04, 0.16502051E-03,-0.12860420E-03,
     +  0.11237734E-03,-0.11997651E-03, 0.32356690E-03,-0.52144588E-03,
     + -0.29533982E-03,-0.17628948E-03,-0.37902202E-04, 0.11597706E-03,
     +  0.10621958E-03, 0.15936884E-03,-0.70181784E-04,-0.96007796E-04,
     +  0.78524419E-04, 0.68901616E-04, 0.39507230E-03,-0.10366549E-03,
     +  0.11411088E-03,-0.53952442E-03,-0.20555615E-03, 0.40781425E-03,
     +  0.13694144E-03, 0.12029154E-03,-0.19281305E-03, 0.39407113E-03,
     + -0.11053778E-03,-0.55235618E-04, 0.88420493E-04,-0.10526926E-03,
     +  0.97063996E-04,-0.14660617E-03,
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
      t_s4_q1en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x23    *x41    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)    *x21*x33        
      t_s4_q1en   =t_s4_q1en   
     9  +coeff(  9)            *x41    
     1  +coeff( 10)    *x22            
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)    *x21    *x41*x51
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 18)    *x24            
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)*x11*x22    *x41    
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x23*x33        
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23*x31*x41    
     7  +coeff( 25)        *x31        
     8  +coeff( 26)            *x42    
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 27)*x11*x21            
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)    *x23*x31        
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)    *x21*x31*x42    
     5  +coeff( 32)*x12*x21            
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)    *x24*x31        
     8  +coeff( 35)*x11*x24            
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)    *x24    *x42    
     2  +coeff( 38)    *x23    *x43    
     3  +coeff( 39)*x11        *x44    
     4  +coeff( 40)*x11*x24    *x41    
     5  +coeff( 41)*x11*x23    *x42    
     6  +coeff( 42)    *x23    *x44    
     7  +coeff( 43)*x11*x23    *x41*x51
     8  +coeff( 44)*x11*x23    *x42*x51
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 45)*x11*x24    *x41*x53
     1  +coeff( 46)*x11*x23*x34*x41*x51
     2  +coeff( 47)*x11*x24    *x42*x53
     3  +coeff( 48)*x11*x23*x33*x44    
     4  +coeff( 49)*x12*x22*x33*x43*x51
     5  +coeff( 50)                *x51
     6  +coeff( 51)        *x31*x41    
     7  +coeff( 52)*x11    *x31        
     8  +coeff( 53)            *x43    
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 54)*x11            *x51
     1  +coeff( 55)*x11        *x42    
     2  +coeff( 56)*x11    *x31    *x51
     3  +coeff( 57)        *x33    *x51
     4  +coeff( 58)*x11        *x41*x51
     5  +coeff( 59)    *x21    *x42*x51
     6  +coeff( 60)    *x21    *x41*x52
     7  +coeff( 61)*x11*x22*x31        
     8  +coeff( 62)*x11*x21*x31*x41    
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 63)    *x22    *x43    
     1  +coeff( 64)    *x21    *x44    
     2  +coeff( 65)*x11*x21    *x41*x51
     3  +coeff( 66)    *x23    *x41*x51
     4  +coeff( 67)            *x44*x51
     5  +coeff( 68)    *x23        *x52
     6  +coeff( 69)    *x22*x31    *x52
     7  +coeff( 70)    *x21    *x42*x52
     8  +coeff( 71)    *x21        *x54
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 72)*x12*x22            
     1  +coeff( 73)*x12*x21    *x41    
     2  +coeff( 74)*x11*x21*x32*x41    
     3  +coeff( 75)    *x23*x31*x42    
     4  +coeff( 76)    *x22*x32*x42    
     5  +coeff( 77)    *x23*x31*x41*x51
     6  +coeff( 78)*x11*x21    *x42*x51
     7  +coeff( 79)    *x23    *x42*x51
     8  +coeff( 80)    *x23    *x41*x52
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 81)    *x22*x31*x41*x52
     1  +coeff( 82)*x12*x23            
     2  +coeff( 83)*x12*x22    *x41    
     3  +coeff( 84)*x11*x23*x31*x41    
     4  +coeff( 85)    *x24*x32*x41    
     5  +coeff( 86)*x12    *x31*x42    
     6  +coeff( 87)*x12        *x43    
     7  +coeff( 88)*x11*x22*x32    *x51
     8  +coeff( 89)    *x24*x32    *x51
      t_s4_q1en   =t_s4_q1en   
     9  +coeff( 90)*x11*x22*x31*x41*x51
c
      return
      end
      function y_s4_q1en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.1816674E-01/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35381860E-02, 0.53161834E-02, 0.66067085E-01, 0.46949200E-02,
     +  0.13243441E-02,-0.62235030E-02, 0.24482096E-02,-0.59206239E-02,
     + -0.76449208E-03, 0.25636220E-03,-0.86385274E-03, 0.22209506E-02,
     + -0.70232432E-03, 0.95490721E-03, 0.34817180E-02,-0.25500827E-02,
     + -0.78152103E-03, 0.26970622E-02,-0.21598759E-03,-0.69239776E-03,
     +  0.53337665E-03, 0.20656097E-02, 0.56129252E-03, 0.63855125E-03,
     +  0.26455678E-04, 0.21147150E-03, 0.20296207E-03,-0.36276501E-03,
     +  0.13463151E-02,-0.38617899E-03, 0.12173373E-02, 0.67854422E-03,
     + -0.50301307E-04,-0.71520000E-04,-0.54846063E-04,-0.14647104E-03,
     + -0.90218295E-04, 0.25242366E-03, 0.88660337E-04, 0.78989012E-03,
     + -0.14638850E-03, 0.29928150E-03,
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
      x43 = x42*x4
      x44 = x43*x4
      x45 = x44*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_s4_q1en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x43    
     8  +coeff(  8)    *x22    *x41    
      y_s4_q1en   =y_s4_q1en   
     9  +coeff(  9)    *x21            
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)            *x44    
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)*x11*x21    *x41    
      y_s4_q1en   =y_s4_q1en   
     9  +coeff( 18)    *x24    *x41    
     1  +coeff( 19)                *x52
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x23            
     4  +coeff( 22)        *x31*x43    
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)        *x32        
     8  +coeff( 26)        *x32*x41    
      y_s4_q1en   =y_s4_q1en   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x22*x31*x41    
     2  +coeff( 29)            *x45    
     3  +coeff( 30)        *x31*x45    
     4  +coeff( 31)    *x24    *x42    
     5  +coeff( 32)*x11*x23    *x41    
     6  +coeff( 33)    *x21*x31        
     7  +coeff( 34)            *x41*x51
     8  +coeff( 35)*x11                
      y_s4_q1en   =y_s4_q1en   
     9  +coeff( 36)    *x21    *x42    
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)        *x32*x42    
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)        *x31*x44    
     5  +coeff( 41)*x11*x21    *x42    
     6  +coeff( 42)    *x24*x31        
c
      return
      end
      function p_s4_q1en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.4871928E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.21039415E-02,-0.79420861E-03, 0.25212241E-01, 0.66573243E-02,
     + -0.70900484E-02, 0.15110486E-02,-0.87661477E-03,-0.57513793E-02,
     +  0.34044432E-02, 0.35304294E-03, 0.42852430E-03, 0.23935253E-02,
     + -0.33002070E-02, 0.30234365E-02, 0.56339224E-03,-0.69650362E-03,
     + -0.27870623E-03, 0.48571709E-03, 0.95500553E-03,-0.79244433E-03,
     +  0.13802858E-02, 0.26173827E-02,-0.54361641E-04,-0.62396977E-03,
     +  0.23013930E-03, 0.55916270E-03,-0.60934195E-03, 0.57937612E-03,
     + -0.70775161E-04, 0.31707354E-04, 0.69811322E-04,-0.19928094E-03,
     + -0.30357314E-04, 0.22307619E-03, 0.24935592E-03,-0.46696974E-03,
     +  0.71921374E-03, 0.16125778E-02,-0.87836161E-04, 0.28818111E-04,
     + -0.10197550E-03,-0.56311699E-04, 0.10484069E-03,-0.14924281E-03,
     +  0.23053013E-03, 0.21049002E-03, 0.22253675E-03, 0.19238693E-03,
     + -0.98546770E-05,-0.46381101E-04,-0.60209986E-04, 0.20995969E-04,
     +  0.21866490E-04, 0.71522823E-04,-0.44244989E-04, 0.16705817E-03,
     +  0.85611216E-04,-0.12520469E-03, 0.79334772E-04,-0.17205399E-03,
     + -0.14447772E-03, 0.19988457E-03,-0.13785267E-04,-0.71793664E-04,
     +  0.38554677E-03, 0.11858132E-03, 0.12581385E-03, 0.18525205E-03,
     +  0.13691623E-03, 0.10556947E-03,-0.82219827E-04, 0.15177955E-03,
     + -0.25035677E-03, 0.18929873E-03,-0.44763705E-03, 0.16901425E-03,
     + -0.36213103E-04,-0.26127165E-04,-0.19985600E-04,-0.35615125E-04,
     +  0.33108103E-04, 0.15008866E-04,-0.37808841E-04,-0.39613380E-04,
     + -0.36370471E-04, 0.45329813E-04,-0.24825402E-04, 0.47638350E-04,
     + -0.47023357E-04, 0.43663225E-04,
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
      p_s4_q1en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22    *x41    
      p_s4_q1en   =p_s4_q1en   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)        *x31        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)            *x44    
     6  +coeff( 15)        *x31*x44    
     7  +coeff( 16)    *x21    *x41    
     8  +coeff( 17)                *x52
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 18)    *x23            
     1  +coeff( 19)        *x31*x42    
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)        *x31*x43    
     4  +coeff( 22)    *x24    *x41    
     5  +coeff( 23)            *x41*x51
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)    *x22        *x51
     8  +coeff( 26)    *x23    *x41    
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 27)    *x22*x31*x41    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)*x11                
     3  +coeff( 30)        *x32        
     4  +coeff( 31)        *x32*x41    
     5  +coeff( 32)    *x21    *x42    
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)        *x32*x42    
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)    *x24    *x42    
     3  +coeff( 39)    *x21*x31        
     4  +coeff( 40)    *x21        *x51
     5  +coeff( 41)    *x21*x31*x41    
     6  +coeff( 42)*x11*x21*x31        
     7  +coeff( 43)*x11*x21        *x51
     8  +coeff( 44)            *x43*x51
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 45)    *x24*x31        
     1  +coeff( 46)    *x23    *x42    
     2  +coeff( 47)        *x32*x43    
     3  +coeff( 48)    *x22    *x41*x53
     4  +coeff( 49)        *x31    *x51
     5  +coeff( 50)*x11        *x41    
     6  +coeff( 51)        *x31*x41*x51
     7  +coeff( 52)        *x31    *x52
     8  +coeff( 53)                *x53
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 54)    *x23*x31        
     1  +coeff( 55)        *x32*x41*x51
     2  +coeff( 56)*x11*x22    *x41    
     3  +coeff( 57)    *x23*x31*x41    
     4  +coeff( 58)    *x21    *x44    
     5  +coeff( 59)    *x22    *x42*x51
     6  +coeff( 60)            *x44*x51
     7  +coeff( 61)*x11*x24            
     8  +coeff( 62)    *x24*x31*x41    
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 63)*x12            *x52
     1  +coeff( 64)*x11*x21        *x53
     2  +coeff( 65)*x11*x23    *x42    
     3  +coeff( 66)*x11*x21*x32*x42    
     4  +coeff( 67)*x12*x21    *x41*x51
     5  +coeff( 68)*x11*x23    *x41*x51
     6  +coeff( 69)*x11    *x32*x42*x51
     7  +coeff( 70)*x11*x24*x31*x41    
     8  +coeff( 71)*x11*x21*x34    *x51
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 72)*x11*x24*x31*x42    
     1  +coeff( 73)*x11*x24    *x43    
     2  +coeff( 74)*x11*x23*x31*x41*x52
     3  +coeff( 75)*x11    *x34*x44*x51
     4  +coeff( 76)*x12    *x34*x44    
     5  +coeff( 77)*x11    *x31        
     6  +coeff( 78)    *x21    *x41*x51
     7  +coeff( 79)*x12                
     8  +coeff( 80)    *x22*x32        
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 81)    *x21    *x43    
     1  +coeff( 82)*x11        *x41*x51
     2  +coeff( 83)        *x31*x42*x51
     3  +coeff( 84)            *x41*x53
     4  +coeff( 85)                *x54
     5  +coeff( 86)*x11    *x33        
     6  +coeff( 87)*x12        *x41    
     7  +coeff( 88)    *x22*x32*x41    
     8  +coeff( 89)*x11*x21    *x41*x51
      p_s4_q1en   =p_s4_q1en   
     9  +coeff( 90)        *x33*x41*x51
c
      return
      end
      function l_s4_q1en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/  0.1074375E-02/
      data xmin/
     1 -0.49937E-02,-0.38763E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.38139E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50006429E-03, 0.81611361E-04,-0.75616973E-03,-0.69485572E-02,
     + -0.14816542E-03,-0.16675002E-02,-0.84901042E-03, 0.66324422E-03,
     + -0.11069790E-03,-0.26860493E-03,-0.20208802E-03,-0.33057921E-04,
     +  0.30836076E-03,-0.21229808E-03,-0.31316976E-03, 0.73731685E-05,
     +  0.57675501E-04, 0.56403576E-04,-0.42089421E-04, 0.39053848E-04,
     +  0.68334981E-06,-0.95717383E-04, 0.69687871E-04,-0.77654280E-04,
     + -0.51867453E-04, 0.29567365E-04,-0.55901837E-04,-0.67712375E-04,
     + -0.15969886E-03,
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
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
c
c                  function
c
      l_s4_q1en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x22    *x41    
      l_s4_q1en   =l_s4_q1en   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x24            
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)            *x44    
     6  +coeff( 15)    *x24    *x41    
     7  +coeff( 16)        *x33*x42    
     8  +coeff( 17)    *x21    *x41    
      l_s4_q1en   =l_s4_q1en   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)        *x33        
     4  +coeff( 22)        *x31*x42    
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x22*x31*x41    
      l_s4_q1en   =l_s4_q1en   
     9  +coeff( 27)*x11*x23            
     1  +coeff( 28)*x11*x23    *x41    
     2  +coeff( 29)    *x24    *x42    
c
      return
      end
      function x_s4_sext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2737899E+00/
      data xmin/
     1 -0.49937E-02,-0.40148E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.40930E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.32803116E-02, 0.76576328E-03,-0.53263293E-02,-0.61393525E-01,
     + -0.33223655E-02, 0.51414650E-02,-0.88237255E-03, 0.51078033E-02,
     + -0.28667129E-02,-0.19903900E-02,-0.12280627E-04, 0.83759637E-03,
     + -0.53859333E-03,-0.73499640E-03,-0.11732860E-02,-0.24921200E-02,
     +  0.64482098E-03,-0.24716058E-02, 0.73346891E-04, 0.13554898E-03,
     +  0.68935775E-03, 0.16642485E-02, 0.66726474E-03,-0.48083780E-03,
     + -0.20432284E-04,-0.28982112E-03,-0.21061509E-03,-0.76542015E-03,
     +  0.37379752E-03,-0.34562454E-05,-0.63636038E-03, 0.11431411E-03,
     +  0.18910784E-03,-0.44398453E-05, 0.54637752E-04,-0.14415945E-03,
     + -0.51206491E-04,-0.35011463E-03,-0.53249986E-03,-0.12163774E-03,
     + -0.87605388E-03,-0.36797606E-03,-0.35945322E-04, 0.13099874E-04,
     +  0.64146909E-04,-0.10294589E-03, 0.48868973E-04, 0.17339519E-03,
     +  0.11618751E-03,-0.19482557E-03,-0.39470900E-03,-0.14203566E-03,
     +  0.30027984E-04, 0.10426672E-03,-0.29749481E-04, 0.29900562E-03,
     +  0.59619033E-04,-0.14189350E-03,-0.14045319E-03,-0.32895507E-03,
     +  0.50849654E-03, 0.73394527E-04,-0.37424576E-04,-0.12239575E-04,
     +  0.24173254E-04,-0.44591976E-04,-0.19922362E-03, 0.27710317E-04,
     +  0.36771926E-04,-0.63520536E-04,-0.79970203E-04, 0.39047110E-04,
     +  0.62572020E-04, 0.12793478E-04, 0.14735550E-03, 0.26774769E-04,
     + -0.96890377E-04,-0.42897291E-05, 0.38626789E-04,-0.61020866E-04,
     +  0.44974917E-04,-0.16671109E-03,-0.51239520E-04, 0.25866716E-03,
     +  0.16021584E-03,-0.13692745E-03, 0.41101375E-05, 0.82270621E-04,
     + -0.15842036E-03,-0.18091507E-04,
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
      x_s4_sext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x22    *x41    
      x_s4_sext   =x_s4_sext   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)    *x24            
     2  +coeff( 11)        *x33*x41    
     3  +coeff( 12)    *x21    *x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)        *x31*x43    
     7  +coeff( 16)            *x44    
     8  +coeff( 17)*x11*x21            
      x_s4_sext   =x_s4_sext   
     9  +coeff( 18)    *x24    *x41    
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)                *x52
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)        *x32        
     8  +coeff( 26)        *x31*x41    
      x_s4_sext   =x_s4_sext   
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)        *x34*x41    
     4  +coeff( 31)*x11*x23    *x41    
     5  +coeff( 32)    *x21*x31        
     6  +coeff( 33)    *x21    *x42    
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)*x11                
      x_s4_sext   =x_s4_sext   
     9  +coeff( 36)        *x32*x42    
     1  +coeff( 37)    *x22    *x41*x51
     2  +coeff( 38)    *x24*x31        
     3  +coeff( 39)        *x31*x44    
     4  +coeff( 40)*x11*x22            
     5  +coeff( 41)    *x24    *x42    
     6  +coeff( 42)    *x22*x32*x43    
     7  +coeff( 43)    *x21        *x51
     8  +coeff( 44)        *x31    *x51
      x_s4_sext   =x_s4_sext   
     9  +coeff( 45)    *x21*x31*x41    
     1  +coeff( 46)    *x23*x31        
     2  +coeff( 47)        *x32*x41*x51
     3  +coeff( 48)            *x43*x51
     4  +coeff( 49)    *x22*x32*x41    
     5  +coeff( 50)    *x23    *x42    
     6  +coeff( 51)    *x22*x31*x42    
     7  +coeff( 52)        *x32*x43    
     8  +coeff( 53)    *x22*x31*x41*x51
      x_s4_sext   =x_s4_sext   
     9  +coeff( 54)*x11*x21*x31        
     1  +coeff( 55)*x11*x21        *x51
     2  +coeff( 56)*x11*x21    *x42    
     3  +coeff( 57)*x11        *x43    
     4  +coeff( 58)*x11*x23*x31        
     5  +coeff( 59)*x11*x24    *x41    
     6  +coeff( 60)*x11*x23    *x42    
     7  +coeff( 61)    *x24*x31*x44    
     8  +coeff( 62)*x11    *x32*x42*x52
      x_s4_sext   =x_s4_sext   
     9  +coeff( 63)        *x32*x41    
     1  +coeff( 64)    *x21    *x41*x51
     2  +coeff( 65)        *x31    *x52
     3  +coeff( 66)    *x21*x31*x42    
     4  +coeff( 67)    *x21    *x43    
     5  +coeff( 68)    *x23        *x51
     6  +coeff( 69)                *x54
     7  +coeff( 70)    *x23*x31*x41    
     8  +coeff( 71)    *x22    *x43    
      x_s4_sext   =x_s4_sext   
     9  +coeff( 72)    *x21    *x44    
     1  +coeff( 73)    *x24        *x51
     2  +coeff( 74)        *x34    *x51
     3  +coeff( 75)            *x44*x51
     4  +coeff( 76)        *x33    *x52
     5  +coeff( 77)        *x32*x41*x52
     6  +coeff( 78)            *x43*x52
     7  +coeff( 79)        *x31*x41*x53
     8  +coeff( 80)        *x31    *x54
      x_s4_sext   =x_s4_sext   
     9  +coeff( 81)            *x41*x54
     1  +coeff( 82)    *x24*x31*x41    
     2  +coeff( 83)    *x23*x32*x41    
     3  +coeff( 84)    *x23    *x43    
     4  +coeff( 85)    *x22    *x44    
     5  +coeff( 86)        *x32*x44    
     6  +coeff( 87)    *x24*x31    *x51
     7  +coeff( 88)    *x22*x31*x42*x51
     8  +coeff( 89)    *x22    *x43*x51
      x_s4_sext   =x_s4_sext   
     9  +coeff( 90)*x11*x22*x31        
c
      return
      end
      function t_s4_sext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2215329E+00/
      data xmin/
     1 -0.49937E-02,-0.40148E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.40930E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.28304323E-02, 0.12616999E-02,-0.26433457E-01,-0.67255595E-02,
     +  0.92966100E-02,-0.11433696E-02, 0.10643867E-02, 0.65937913E-02,
     + -0.35217495E-02,-0.37474363E-03,-0.38970844E-03,-0.27269884E-02,
     +  0.29816893E-02,-0.33749128E-02, 0.92519540E-03, 0.31621009E-03,
     + -0.70450886E-03, 0.74533612E-03,-0.97955053E-03, 0.82678435E-03,
     + -0.31616264E-02, 0.15238018E-03,-0.32949142E-03,-0.81426575E-03,
     +  0.58885559E-03,-0.14516134E-02,-0.64594758E-03, 0.83215025E-04,
     + -0.30300494E-04,-0.13283830E-04, 0.43332297E-03, 0.47322555E-04,
     + -0.15825726E-03, 0.21163614E-04,-0.76491805E-03,-0.16041065E-02,
     +  0.12821391E-03,-0.76061282E-04, 0.13477827E-03,-0.13558834E-03,
     + -0.19320597E-03,-0.53869506E-04,-0.53478903E-04, 0.20220804E-03,
     +  0.37239097E-04,-0.34716423E-03, 0.43097709E-03,-0.45037686E-03,
     + -0.61525474E-03,-0.23844132E-03, 0.65145854E-04,-0.50023050E-05,
     +  0.13783926E-03, 0.62611485E-04, 0.50594321E-04, 0.85049702E-04,
     + -0.25766736E-03,-0.13815657E-03,-0.25919368E-03,-0.22934547E-03,
     + -0.79111211E-04, 0.13826138E-03,-0.24791567E-04,-0.18839492E-03,
     + -0.26500036E-03, 0.31543066E-03,-0.52066061E-04, 0.49194728E-04,
     + -0.17233942E-03,-0.41210442E-03, 0.17990806E-03, 0.49770327E-03,
     + -0.72975454E-05,-0.22462441E-03, 0.11783613E-03, 0.11056875E-03,
     +  0.15330910E-03,-0.19377023E-03, 0.11631901E-03,-0.26633125E-03,
     + -0.27362246E-03, 0.24326303E-03,-0.24501927E-03,-0.14563675E-03,
     +  0.25114181E-03,-0.18297668E-03, 0.32092421E-03, 0.32789301E-03,
     +  0.15438947E-04, 0.37025231E-04,
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
      t_s4_sext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22    *x41    
      t_s4_sext   =t_s4_sext   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)        *x31        
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)    *x24            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x24    *x41    
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)                *x52
     8  +coeff( 17)    *x23            
      t_s4_sext   =t_s4_sext   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)        *x31*x42    
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)            *x44    
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x23    *x41    
     7  +coeff( 25)    *x22*x31*x41    
     8  +coeff( 26)        *x31*x43    
      t_s4_sext   =t_s4_sext   
     9  +coeff( 27)*x11*x23            
     1  +coeff( 28)*x11                
     2  +coeff( 29)        *x32        
     3  +coeff( 30)        *x32*x41    
     4  +coeff( 31)    *x21    *x42    
     5  +coeff( 32)            *x42*x51
     6  +coeff( 33)*x11*x22            
     7  +coeff( 34)    *x21*x31    *x52
     8  +coeff( 35)*x11*x23    *x41    
      t_s4_sext   =t_s4_sext   
     9  +coeff( 36)    *x24    *x42    
     1  +coeff( 37)    *x21*x31        
     2  +coeff( 38)    *x21        *x51
     3  +coeff( 39)    *x21*x31*x41    
     4  +coeff( 40)    *x23*x31        
     5  +coeff( 41)        *x32*x42    
     6  +coeff( 42)*x11*x21        *x51
     7  +coeff( 43)    *x22    *x41*x51
     8  +coeff( 44)            *x43*x51
      t_s4_sext   =t_s4_sext   
     9  +coeff( 45)        *x31    *x53
     1  +coeff( 46)    *x24*x31        
     2  +coeff( 47)*x11*x21    *x42    
     3  +coeff( 48)    *x23    *x42    
     4  +coeff( 49)        *x31*x44    
     5  +coeff( 50)    *x22*x32*x43    
     6  +coeff( 51)*x11        *x41    
     7  +coeff( 52)        *x31*x41*x51
     8  +coeff( 53)*x11*x21*x31        
      t_s4_sext   =t_s4_sext   
     9  +coeff( 54)    *x23        *x51
     1  +coeff( 55)        *x32*x41*x51
     2  +coeff( 56)*x11*x22*x31        
     3  +coeff( 57)*x11*x22    *x41    
     4  +coeff( 58)    *x23*x31*x41    
     5  +coeff( 59)    *x22    *x43    
     6  +coeff( 60)        *x32*x43    
     7  +coeff( 61)    *x22    *x42*x51
     8  +coeff( 62)            *x44*x51
      t_s4_sext   =t_s4_sext   
     9  +coeff( 63)        *x31    *x54
     1  +coeff( 64)*x11*x23*x31        
     2  +coeff( 65)    *x24*x31*x41    
     3  +coeff( 66)    *x22    *x44    
     4  +coeff( 67)*x11        *x43*x51
     5  +coeff( 68)*x12            *x52
     6  +coeff( 69)*x11*x24*x31        
     7  +coeff( 70)*x11*x23    *x42    
     8  +coeff( 71)*x11*x22    *x43    
      t_s4_sext   =t_s4_sext   
     9  +coeff( 72)    *x24    *x43    
     1  +coeff( 73)*x12*x21    *x41*x51
     2  +coeff( 74)        *x34*x41*x52
     3  +coeff( 75)        *x31*x43*x53
     4  +coeff( 76)        *x32*x41*x54
     5  +coeff( 77)*x12*x22*x31*x41    
     6  +coeff( 78)        *x34*x44    
     7  +coeff( 79)    *x24*x31*x42*x51
     8  +coeff( 80)    *x24    *x43*x51
      t_s4_sext   =t_s4_sext   
     9  +coeff( 81)*x11*x23*x32*x41*x51
     1  +coeff( 82)    *x23    *x44*x52
     2  +coeff( 83)*x12*x21    *x41*x53
     3  +coeff( 84)*x11    *x33*x44*x51
     4  +coeff( 85)    *x24*x32    *x54
     5  +coeff( 86)*x12*x24        *x54
     6  +coeff( 87)*x11    *x32*x44*x54
     7  +coeff( 88)*x12    *x32*x43*x54
     8  +coeff( 89)    *x21*x31    *x51
      t_s4_sext   =t_s4_sext   
     9  +coeff( 90)    *x21    *x41*x51
c
      return
      end
      function y_s4_sext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.6274225E-03/
      data xmin/
     1 -0.49937E-02,-0.40148E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.40930E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12514188E-02,-0.43560158E-04,-0.48404798E-03, 0.10307962E+00,
     + -0.54187793E-02, 0.46245093E-02, 0.25369795E-02, 0.37207981E-02,
     +  0.10152955E-02,-0.23224610E-02,-0.42119698E-03, 0.14953353E-02,
     +  0.60877681E-03,-0.41842589E-03, 0.15636861E-03,-0.59218117E-03,
     + -0.32586692E-03, 0.18692658E-03,-0.61650254E-03, 0.19242110E-02,
     +  0.66047953E-03,-0.94352686E-03, 0.71148219E-03,-0.18273527E-03,
     +  0.59707243E-04,-0.35705912E-03, 0.41726095E-03, 0.12684194E-03,
     +  0.47007774E-03,
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
      x33 = x32*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
c
c                  function
c
      y_s4_sext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x21            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)*x11                
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x23    *x41    
      y_s4_sext   =y_s4_sext   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21*x33        
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)    *x21*x31        
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)    *x21*x31*x41    
     8  +coeff( 17)*x11        *x41    
      y_s4_sext   =y_s4_sext   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x24            
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)    *x24    *x41    
     5  +coeff( 23)    *x23*x33        
     6  +coeff( 24)            *x42    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x21    *x43    
      y_s4_sext   =y_s4_sext   
     9  +coeff( 27)    *x22    *x42    
     1  +coeff( 28)    *x22*x33        
     2  +coeff( 29)    *x23*x31*x41    
c
      return
      end
      function p_s4_sext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5018020E-03/
      data xmin/
     1 -0.49937E-02,-0.40148E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.40930E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35297172E-03, 0.34210648E-01,-0.64850152E-02, 0.29084990E-02,
     +  0.45526507E-02, 0.94588124E-03,-0.30382045E-02,-0.50321780E-03,
     + -0.49843360E-03,-0.83587470E-03, 0.16720511E-02, 0.81405882E-03,
     +  0.23506994E-03,-0.79612795E-03,-0.40352668E-03, 0.28655084E-02,
     +  0.15876956E-03,-0.44894469E-03, 0.29571072E-03, 0.65807346E-03,
     + -0.13347723E-02,-0.65653789E-04, 0.13454786E-02,-0.97653648E-03,
     +  0.86676696E-03,-0.15967514E-02,-0.96104821E-04, 0.43759643E-03,
     + -0.20363573E-03, 0.18173406E-02,-0.83137857E-04,-0.14690048E-03,
     +  0.40877637E-03, 0.41267279E-03, 0.10713200E-02,-0.16821337E-03,
     + -0.48537360E-03,-0.67729875E-03, 0.26365030E-02, 0.16593296E-04,
     + -0.95448588E-04,-0.30697510E-03,-0.66691900E-04,-0.10168946E-03,
     + -0.81479899E-04, 0.63839584E-03, 0.10189680E-03, 0.55369579E-04,
     +  0.39284569E-04, 0.10864886E-03,-0.41132761E-03, 0.46788482E-05,
     +  0.97977382E-03,-0.11520822E-03,-0.40258886E-03, 0.54162494E-04,
     + -0.17093269E-03, 0.17908518E-03, 0.53455900E-04,-0.19943270E-03,
     + -0.74843463E-03,-0.37466318E-03, 0.38623410E-04, 0.18886756E-03,
     + -0.12206996E-03, 0.14505698E-03, 0.10542122E-03,-0.73725518E-04,
     +  0.85504515E-04,-0.54352457E-03,-0.20445722E-03,-0.23521531E-03,
     +  0.20295102E-03, 0.34604859E-03,-0.91462262E-03, 0.36049378E-03,
     +  0.12777167E-03,-0.35296491E-03, 0.14937717E-03, 0.14832256E-03,
     + -0.19918090E-03,-0.19985756E-03,-0.17281536E-03,-0.21399214E-03,
     +  0.36855767E-03, 0.48613313E-03,-0.48840203E-03, 0.62770961E-03,
     + -0.25711374E-03,-0.29331626E-03,
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
      x34 = x33*x3
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_s4_sext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x23    *x41    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x42    
     8  +coeff(  8)            *x41    
      p_s4_sext   =p_s4_sext   
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x24            
     7  +coeff( 16)    *x23    *x42    
     8  +coeff( 17)*x11*x21            
      p_s4_sext   =p_s4_sext   
     9  +coeff( 18)*x11        *x41    
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)    *x23*x31        
     3  +coeff( 21)    *x21    *x43    
     4  +coeff( 22)            *x44    
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x24    *x41    
     7  +coeff( 25)    *x23*x31*x41    
     8  +coeff( 26)    *x24    *x42    
      p_s4_sext   =p_s4_sext   
     9  +coeff( 27)        *x31        
     1  +coeff( 28)    *x22*x31        
     2  +coeff( 29)            *x43    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)    *x21*x31*x42    
     5  +coeff( 32)*x11*x23            
     6  +coeff( 33)*x11*x21*x32*x41    
     7  +coeff( 34)*x11*x22    *x42    
     8  +coeff( 35)    *x23    *x43    
      p_s4_sext   =p_s4_sext   
     9  +coeff( 36)*x11        *x44    
     1  +coeff( 37)*x11*x23    *x44    
     2  +coeff( 38)*x11*x23*x32*x43    
     3  +coeff( 39)*x12*x23*x34*x43*x52
     4  +coeff( 40)                *x51
     5  +coeff( 41)        *x31*x41    
     6  +coeff( 42)            *x42    
     7  +coeff( 43)*x11    *x31        
     8  +coeff( 44)        *x31*x41*x51
      p_s4_sext   =p_s4_sext   
     9  +coeff( 45)    *x22*x32        
     1  +coeff( 46)    *x22*x31*x41    
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)        *x31*x42*x51
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)*x11*x22*x31        
     6  +coeff( 51)    *x24*x31        
     7  +coeff( 52)*x11*x21    *x42    
     8  +coeff( 53)    *x22    *x43    
      p_s4_sext   =p_s4_sext   
     9  +coeff( 54)    *x21*x31*x43    
     1  +coeff( 55)    *x21    *x44    
     2  +coeff( 56)*x11*x21*x31    *x51
     3  +coeff( 57)    *x23    *x41*x51
     4  +coeff( 58)        *x31*x43*x51
     5  +coeff( 59)    *x21*x31*x41*x52
     6  +coeff( 60)*x11*x24            
     7  +coeff( 61)    *x24*x31*x41    
     8  +coeff( 62)    *x21*x31*x44    
      p_s4_sext   =p_s4_sext   
     9  +coeff( 63)*x11*x23        *x51
     1  +coeff( 64)*x11*x22*x31    *x51
     2  +coeff( 65)*x11*x22    *x41*x51
     3  +coeff( 66)*x11*x21*x31*x41*x51
     4  +coeff( 67)*x11        *x43*x51
     5  +coeff( 68)*x11*x21        *x53
     6  +coeff( 69)    *x21*x31*x41*x53
     7  +coeff( 70)*x11*x24    *x41    
     8  +coeff( 71)*x12*x21*x31*x41    
      p_s4_sext   =p_s4_sext   
     9  +coeff( 72)    *x24*x32*x41    
     1  +coeff( 73)*x11*x21*x32*x42    
     2  +coeff( 74)    *x23*x32*x42    
     3  +coeff( 75)    *x24    *x43    
     4  +coeff( 76)*x11*x23    *x41*x51
     5  +coeff( 77)*x12    *x31*x41*x51
     6  +coeff( 78)*x11*x22*x31*x41*x51
     7  +coeff( 79)*x11    *x31*x43*x51
     8  +coeff( 80)*x11*x22*x31    *x52
      p_s4_sext   =p_s4_sext   
     9  +coeff( 81)*x11*x21    *x41*x53
     1  +coeff( 82)*x11*x24*x32        
     2  +coeff( 83)*x12*x21*x31*x42    
     3  +coeff( 84)*x11*x23    *x43    
     4  +coeff( 85)    *x23*x32*x43    
     5  +coeff( 86)    *x23*x31*x44    
     6  +coeff( 87)*x11*x24*x31    *x51
     7  +coeff( 88)*x11*x23*x32    *x51
     8  +coeff( 89)*x11*x21*x34    *x51
      p_s4_sext   =p_s4_sext   
     9  +coeff( 90)*x11*x22*x32*x41*x51
c
      return
      end
      function l_s4_sext   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.1905026E-02/
      data xmin/
     1 -0.49937E-02,-0.40148E-01,-0.49973E-02,-0.26996E-01,-0.44954E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.49980E-02, 0.40930E-01, 0.49977E-02, 0.17284E-01, 0.44909E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10813222E-02,-0.12356169E-03, 0.16955681E-03, 0.38320699E-02,
     +  0.42848627E-03,-0.26600705E-02,-0.60609134E-03,-0.97635333E-04,
     + -0.28573326E-03, 0.23989081E-03, 0.27235330E-04,-0.70119670E-04,
     + -0.50151510E-04,-0.28380982E-04, 0.42556854E-04, 0.10437780E-03,
     +  0.14411876E-03,-0.28289758E-04, 0.21565717E-03,-0.38325124E-05,
     + -0.26250917E-04, 0.31731633E-04, 0.34571472E-04, 0.91446542E-04,
     +  0.36619411E-04, 0.10077246E-03, 0.82831481E-04,-0.55438391E-05,
     +  0.68853928E-05,
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
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_s4_sext   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)                *x51
     6  +coeff(  6)    *x22            
     7  +coeff(  7)            *x42    
     8  +coeff(  8)            *x41*x51
      l_s4_sext   =l_s4_sext   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)            *x43    
     2  +coeff( 11)        *x31*x41    
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x23            
     7  +coeff( 16)        *x31*x42    
     8  +coeff( 17)    *x24            
      l_s4_sext   =l_s4_sext   
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)            *x44    
     2  +coeff( 20)*x11                
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)    *x23    *x41    
     6  +coeff( 24)        *x31*x43    
     7  +coeff( 25)*x11*x23            
     8  +coeff( 26)    *x24    *x41    
      l_s4_sext   =l_s4_sext   
     9  +coeff( 27)    *x22    *x43    
     1  +coeff( 28)    *x21*x31        
     2  +coeff( 29)        *x32        
c
      return
      end

