C Forward transfer functions for right hrs with septum based on s5a_dir.dat
c HRS + PREX room temperature septum (right side)
c                     -JJL 11/17/08
c
c This is "tune X" = original tune used to design collimator

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
      function x_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.8618783E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.37417248E-01, 0.52946794E-03,-0.59743272E-03, 0.69543862E+00,
     + -0.11777457E-01, 0.29095370E-01, 0.23988811E-01, 0.18568523E+00,
     +  0.47014922E-01,-0.16980179E-01, 0.76017608E-02,-0.71987608E-02,
     + -0.90300795E-02, 0.21475352E-01,-0.17968938E-01,-0.17889466E-01,
     + -0.30139007E-02, 0.69464813E-02, 0.63501513E-02,-0.64095343E-02,
     +  0.16491048E-01, 0.26438578E-02,-0.29326640E-02, 0.12476296E-02,
     + -0.61341273E-02, 0.86364746E-02, 0.77172485E-02, 0.57302518E-02,
     + -0.15046716E-01,-0.72924918E-02,-0.54498115E-02,-0.91911536E-02,
     +  0.14785999E-01, 0.18972494E-02,-0.32236602E-02,-0.47872835E-02,
     + -0.32195977E-02, 0.24417001E-02, 0.95738526E-02, 0.86964620E-02,
     + -0.91280258E-03,-0.10838992E-02, 0.36019439E-03, 0.14546118E-03,
     +  0.27545157E-02,-0.37949148E-02,-0.16389930E-02,-0.17572411E-02,
     + -0.48408299E-02,-0.70347176E-02, 0.12159898E-01,-0.39486568E-02,
     +  0.37038847E-03,-0.23147953E-02, 0.87913490E-04, 0.68094458E-04,
     + -0.37657137E-02, 0.93632945E-04,-0.11516156E-02,-0.43008034E-02,
     +  0.15879923E-02, 0.24352197E-02, 0.12773297E-01,-0.75864064E-03,
     + -0.11401121E-02, 0.60194633E-02, 0.19226380E-02, 0.62626586E-02,
     + -0.17013619E-03, 0.92473095E-02,-0.58753295E-02, 0.22998229E-02,
     +  0.48764003E-02,-0.16336044E-02, 0.16475671E-02,-0.29945099E-02,
     + -0.10714880E-01, 0.39909491E-02,-0.62374459E-02, 0.35026136E-02,
     +  0.12367755E-02,-0.27203655E-01, 0.26362160E-01,-0.61001494E-04,
     +  0.33745798E-03, 0.95790027E-04,-0.51746349E-03, 0.26275995E-02,
     +  0.24772721E-03, 0.54882839E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      x_sp_fp     =x_sp_fp     
     9  +coeff(  9)                *x53
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)            *x41*x51
      x_sp_fp     =x_sp_fp     
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x24            
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)*x11            *x52
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)    *x22    *x42    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)    *x21    *x41*x52
     4  +coeff( 31)    *x21        *x53
     5  +coeff( 32)    *x24        *x51
     6  +coeff( 33)    *x23    *x41*x51
     7  +coeff( 34)*x11*x22            
     8  +coeff( 35)            *x41*x52
      x_sp_fp     =x_sp_fp     
     9  +coeff( 36)    *x24    *x41    
     1  +coeff( 37)            *x41*x53
     2  +coeff( 38)    *x23    *x42    
     3  +coeff( 39)    *x22    *x42*x51
     4  +coeff( 40)    *x22    *x41*x52
     5  +coeff( 41)*x11        *x41    
     6  +coeff( 42)            *x43    
     7  +coeff( 43)*x11*x21        *x51
     8  +coeff( 44)*x11        *x42    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 45)*x11*x22    *x41    
     1  +coeff( 46)            *x42*x52
     2  +coeff( 47)    *x23        *x52
     3  +coeff( 48)*x11*x24            
     4  +coeff( 49)    *x21    *x42*x52
     5  +coeff( 50)    *x24        *x52
     6  +coeff( 51)    *x23    *x43    
     7  +coeff( 52)    *x23    *x41*x53
     8  +coeff( 53)*x12*x21            
      x_sp_fp     =x_sp_fp     
     9  +coeff( 54)    *x21    *x43    
     1  +coeff( 55)*x11        *x42*x51
     2  +coeff( 56)*x11        *x41*x52
     3  +coeff( 57)    *x22        *x53
     4  +coeff( 58)*x12*x21        *x51
     5  +coeff( 59)    *x21    *x43*x51
     6  +coeff( 60)    *x24    *x41*x51
     7  +coeff( 61)            *x43*x52
     8  +coeff( 62)*x11*x21    *x43    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 63)    *x23    *x42*x51
     1  +coeff( 64)*x11*x21    *x41*x52
     2  +coeff( 65)*x11        *x43*x51
     3  +coeff( 66)    *x22    *x41*x53
     4  +coeff( 67)*x11*x23    *x42    
     5  +coeff( 68)*x11*x22    *x42*x51
     6  +coeff( 69)*x11*x22        *x53
     7  +coeff( 70)    *x23    *x42*x52
     8  +coeff( 71)*x11*x23    *x43    
      x_sp_fp     =x_sp_fp     
     9  +coeff( 72)*x12*x21    *x41*x52
     1  +coeff( 73)    *x21    *x43*x53
     2  +coeff( 74)*x12*x24    *x41    
     3  +coeff( 75)*x12        *x42*x52
     4  +coeff( 76)*x12*x23    *x42    
     5  +coeff( 77)    *x23    *x42*x53
     6  +coeff( 78)*x11*x22    *x43*x52
     7  +coeff( 79)*x11*x23    *x43*x52
     8  +coeff( 80)*x12*x21    *x42*x53
      x_sp_fp     =x_sp_fp     
     9  +coeff( 81)*x12        *x43*x53
     1  +coeff( 82)*x11*x24    *x42*x53
     2  +coeff( 83)*x11*x24    *x43*x53
     3  +coeff( 84)*x12                
     4  +coeff( 85)*x11*x21    *x41    
     5  +coeff( 86)*x12        *x41    
     6  +coeff( 87)*x11*x21    *x42    
     7  +coeff( 88)    *x22    *x43    
     8  +coeff( 89)*x11            *x53
      x_sp_fp     =x_sp_fp     
     9  +coeff( 90)*x11*x23        *x51
c
c      print *, x_sp_fp, "This is a sample to show we are in tune X"
      return
      end
      function t_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.3418943E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.53827497E-02,-0.17868737E-01, 0.23440915E-03, 0.11793327E+00,
     +  0.69333962E-02,-0.10600019E-01,-0.13667849E-02, 0.24625247E-02,
     + -0.14178583E-02,-0.18676152E-02,-0.97170536E-03,-0.57082187E-03,
     +  0.54593798E-03,-0.13544374E-02, 0.11259691E-02, 0.15495217E-02,
     +  0.46711744E-03, 0.10561070E-03, 0.18520495E-03,-0.63700526E-03,
     +  0.50658529E-03, 0.92911761E-03, 0.17783439E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)    *x22    *x41*x51
      t_sp_fp     =t_sp_fp     
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)    *x23        *x51
     4  +coeff( 22)    *x22    *x42*x51
     5  +coeff( 23)    *x23    *x43    
c
      return
      end
      function y_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.2932074E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.68480624E-02,-0.56617539E-01,-0.16990418E-01, 0.31690341E-02,
     + -0.56022597E-02,-0.64988539E-01,-0.21009030E-01, 0.22905756E-01,
     +  0.79020327E-02, 0.12036243E-01,-0.38173571E-01,-0.93722381E-02,
     +  0.32223329E-01, 0.12644882E-02, 0.12438965E-01, 0.15725441E-01,
     + -0.11947018E-01,-0.74405479E-02, 0.19345969E-02, 0.10380307E-01,
     +  0.23037637E-01, 0.10482660E-01,-0.93858568E-02,-0.25369012E-03,
     +  0.53583696E-02, 0.37323008E-02,-0.28354174E-02, 0.11979284E-01,
     +  0.10656082E-02,-0.40965988E-02, 0.10357367E-01, 0.20305406E-01,
     + -0.15237493E-02, 0.65014599E-03, 0.18730487E-02,-0.28591489E-02,
     + -0.45684422E-02, 0.19703570E-02, 0.30134874E-02,-0.14551069E-02,
     +  0.31668083E-02, 0.40095327E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_fp     =y_sp_fp     
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)    *x22            
     2  +coeff( 11)            *x41*x52
     3  +coeff( 12)                *x53
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)            *x41*x53
      y_sp_fp     =y_sp_fp     
     9  +coeff( 18)    *x23            
     1  +coeff( 19)            *x45    
     2  +coeff( 20)    *x21    *x41*x52
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)    *x22        *x52
     5  +coeff( 23)    *x23        *x51
     6  +coeff( 24)    *x21        *x54
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)            *x43*x51
      y_sp_fp     =y_sp_fp     
     9  +coeff( 27)                *x54
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)            *x42*x53
     4  +coeff( 31)    *x22    *x41*x52
     5  +coeff( 32)    *x22    *x42*x53
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)*x11            *x51
     8  +coeff( 35)            *x44    
      y_sp_fp     =y_sp_fp     
     9  +coeff( 36)    *x23    *x41    
     1  +coeff( 37)            *x41*x54
     2  +coeff( 38)*x11*x21    *x42    
     3  +coeff( 39)    *x21    *x42*x52
     4  +coeff( 40)*x11*x22        *x51
     5  +coeff( 41)*x11*x21    *x42*x51
     6  +coeff( 42)    *x23        *x53
c
      return
      end
      function p_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.4905075E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.27091493E-03, 0.11227905E-02,-0.25758455E-01,-0.11475112E-01,
     +  0.69090333E-02, 0.18649099E-01,-0.28616502E-02, 0.36001338E-02,
     + -0.16552670E-01,-0.19107126E-02, 0.87901630E-03,-0.46968702E-02,
     +  0.10949041E-01, 0.72224758E-03, 0.34443515E-02, 0.56009786E-02,
     + -0.56316698E-03, 0.19584596E-02,-0.14860746E-02, 0.83409750E-03,
     +  0.60910889E-03, 0.17309185E-02, 0.93105587E-03, 0.46916073E-02,
     +  0.16170974E-03, 0.49214219E-06,-0.44008941E-03,-0.14974739E-02,
     +  0.46211528E-03,-0.68295328E-03,-0.15818692E-02, 0.29010829E-03,
     +  0.53502532E-03, 0.51235437E-03,-0.62096195E-04, 0.57608431E-03,
     + -0.51731477E-03, 0.41901018E-03, 0.61245497E-04, 0.15018310E-03,
     +  0.77474362E-03, 0.20053520E-03,-0.99398848E-03,-0.62305556E-03,
     + -0.10256209E-03, 0.35367941E-03, 0.95220975E-03, 0.62720943E-03,
     + -0.85544871E-03,-0.98583079E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_fp     =p_sp_fp     
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)*x11        *x41    
      p_sp_fp     =p_sp_fp     
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)*x11*x21    *x42    
     6  +coeff( 24)    *x22    *x43    
     7  +coeff( 25)*x11                
     8  +coeff( 26)*x11            *x51
      p_sp_fp     =p_sp_fp     
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)            *x41*x53
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)                *x53
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)    *x21    *x42*x51
     8  +coeff( 35)*x11*x23            
      p_sp_fp     =p_sp_fp     
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)*x11*x22        *x51
     2  +coeff( 38)*x11*x21    *x41*x51
     3  +coeff( 39)*x11*x21        *x52
     4  +coeff( 40)*x12        *x42    
     5  +coeff( 41)    *x22    *x42*x51
     6  +coeff( 42)    *x23        *x52
     7  +coeff( 43)*x11*x23    *x41    
     8  +coeff( 44)*x11*x23        *x51
      p_sp_fp     =p_sp_fp     
     9  +coeff( 45)*x11*x22        *x52
     1  +coeff( 46)*x11*x21        *x53
     2  +coeff( 47)    *x23    *x41*x52
     3  +coeff( 48)    *x23        *x53
     4  +coeff( 49)*x11*x23        *x52
     5  +coeff( 50)*x11*x22    *x41*x52
c
      return
      end
      function l_sp_fp     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.2604719E+00/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10870708E+00,-0.26718977E+00, 0.49631214E-02,-0.19465520E+01,
     +  0.41215565E-01,-0.97497329E-01,-0.71382023E-01,-0.54067236E+00,
     +  0.83198756E-01,-0.13652910E+00, 0.19659646E-01,-0.31655043E-01,
     + -0.40956184E-01, 0.24763199E-01,-0.17007481E-01, 0.72442517E-01,
     +  0.51772542E-01,-0.72859630E-01,-0.11690352E-01, 0.81129642E-02,
     + -0.41335155E-02, 0.94488198E-02,-0.38131222E-01,-0.33236677E-02,
     +  0.14905673E-01,-0.30596966E-02, 0.93316454E-02,-0.17884996E-01,
     +  0.31609677E-01,-0.15244721E-01, 0.20416852E-01, 0.25673434E-01,
     +  0.54306551E-02,-0.41334815E-01,-0.33748250E-01,-0.24667159E-01,
     + -0.19347504E-01,-0.42343982E-01,-0.26015632E-01,-0.13088172E-01,
     + -0.21692906E-01,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      l_sp_fp     =l_sp_fp     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)                *x53
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x21    *x41*x51
      l_sp_fp     =l_sp_fp     
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x23        *x51
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11            *x52
     5  +coeff( 23)    *x22    *x42    
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)*x11*x22            
      l_sp_fp     =l_sp_fp     
     9  +coeff( 27)    *x21    *x43    
     1  +coeff( 28)    *x22    *x41*x51
     2  +coeff( 29)    *x21    *x42*x51
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)    *x21    *x41*x52
     5  +coeff( 32)    *x21        *x53
     6  +coeff( 33)*x11        *x43    
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x22    *x42*x51
      l_sp_fp     =l_sp_fp     
     9  +coeff( 36)    *x23        *x52
     1  +coeff( 37)*x11*x22    *x42    
     2  +coeff( 38)    *x23    *x43    
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)    *x22    *x43    
     5  +coeff( 41)    *x23        *x53
c
      return
      end
      function x_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1311138E-02/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.78544539E-03, 0.12058984E+00, 0.37662433E-02,-0.34358960E-02,
     +  0.17959137E-02, 0.11758781E-02, 0.21748135E-02,-0.55726717E-03,
     +  0.14353544E-02,-0.13432612E-02, 0.15823025E-03, 0.22543312E-03,
     + -0.58212271E-03,-0.78872294E-03, 0.12463300E-02, 0.14471050E-03,
     + -0.23153120E-03, 0.87196747E-03, 0.58907113E-03, 0.56269702E-04,
     +  0.36897265E-04,-0.23945920E-03, 0.14416910E-03,-0.16250474E-03,
     + -0.41805118E-03, 0.60775219E-05,-0.10167989E-03, 0.17298725E-03,
     +  0.15648766E-04, 0.21087499E-04,-0.31253169E-03, 0.98272511E-04,
     +  0.40120949E-03,-0.15482593E-03, 0.79641014E-03, 0.52254880E-03,
     + -0.25259121E-03,-0.12798267E-03, 0.13322629E-03,-0.23344489E-03,
     + -0.89187066E-04,-0.11839944E-03, 0.24108501E-03, 0.13904247E-03,
     +  0.59240341E-04, 0.86534492E-04, 0.67710411E-04, 0.12392132E-03,
     + -0.26264653E-03, 0.26326661E-04, 0.68062458E-04,-0.40470102E-03,
     + -0.15653810E-03,-0.55991342E-04,-0.10168086E-02,-0.49844930E-04,
     + -0.28813502E-03, 0.78081285E-06,-0.24237440E-03,-0.27235920E-03,
     + -0.23474934E-03, 0.39696076E-03,-0.84920664E-03, 0.12710229E-02,
     + -0.11433578E-02, 0.11376155E-03, 0.36124766E-03, 0.93964633E-03,
     + -0.11824804E-04,-0.35999756E-03, 0.19524177E-02, 0.87005366E-03,
     +  0.26716807E-03,-0.63231259E-04, 0.23589726E-04, 0.19213358E-04,
     +  0.24184736E-04, 0.13441971E-03, 0.16963796E-03, 0.32378703E-04,
     + -0.52969728E-04, 0.18207167E-03,-0.61545914E-04, 0.54336317E-04,
     +  0.32373573E-04,-0.10759479E-03,-0.15424807E-03,-0.86358916E-04,
     + -0.32022715E-03,-0.22510075E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23    *x41    
     8  +coeff(  8)            *x41    
      x_sp_col    =x_sp_col    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21        *x51
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x24    *x41    
     6  +coeff( 15)    *x23    *x42    
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)*x11        *x41    
      x_sp_col    =x_sp_col    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)    *x21    *x43*x51
     3  +coeff( 21)                *x51
     4  +coeff( 22)            *x42    
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)    *x24    *x42    
     8  +coeff( 26)            *x41*x51
      x_sp_col    =x_sp_col    
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x12*x21            
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)    *x23    *x41*x51
     6  +coeff( 33)    *x22    *x43    
     7  +coeff( 34)*x11*x23    *x41    
     8  +coeff( 35)*x11*x22    *x42    
      x_sp_col    =x_sp_col    
     9  +coeff( 36)    *x23    *x43    
     1  +coeff( 37)*x11*x24    *x41    
     2  +coeff( 38)*x12*x21    *x41*x51
     3  +coeff( 39)*x11*x21    *x42*x52
     4  +coeff( 40)*x11*x22    *x41*x53
     5  +coeff( 41)            *x43    
     6  +coeff( 42)*x11        *x42    
     7  +coeff( 43)    *x22    *x41*x51
     8  +coeff( 44)    *x21    *x42*x51
      x_sp_col    =x_sp_col    
     9  +coeff( 45)    *x24        *x51
     1  +coeff( 46)            *x43*x51
     2  +coeff( 47)*x11*x21    *x41*x51
     3  +coeff( 48)*x11*x24            
     4  +coeff( 49)    *x21    *x42*x52
     5  +coeff( 50)    *x21    *x41*x53
     6  +coeff( 51)*x11*x22    *x41*x51
     7  +coeff( 52)    *x24    *x41*x51
     8  +coeff( 53)    *x23    *x42*x51
      x_sp_col    =x_sp_col    
     9  +coeff( 54)*x11*x21    *x41*x52
     1  +coeff( 55)    *x22    *x43*x51
     2  +coeff( 56)*x12*x21    *x42    
     3  +coeff( 57)*x11*x23    *x41*x51
     4  +coeff( 58)*x12*x24            
     5  +coeff( 59)    *x24    *x43    
     6  +coeff( 60)*x11*x21    *x43*x51
     7  +coeff( 61)    *x23    *x43*x51
     8  +coeff( 62)    *x23    *x42*x52
      x_sp_col    =x_sp_col    
     9  +coeff( 63)*x11*x24    *x42    
     1  +coeff( 64)    *x24    *x43*x51
     2  +coeff( 65)*x11*x22    *x42*x52
     3  +coeff( 66)    *x24    *x41*x53
     4  +coeff( 67)*x11*x21    *x43*x52
     5  +coeff( 68)*x11*x23    *x43*x51
     6  +coeff( 69)*x12*x24    *x42    
     7  +coeff( 70)*x12*x24    *x41*x51
     8  +coeff( 71)*x11*x24    *x42*x52
      x_sp_col    =x_sp_col    
     9  +coeff( 72)*x12*x24    *x43*x51
     1  +coeff( 73)*x12*x24    *x42*x52
     2  +coeff( 74)    *x21    *x41*x52
     3  +coeff( 75)    *x21        *x53
     4  +coeff( 76)*x11        *x41*x52
     5  +coeff( 77)    *x22        *x53
     6  +coeff( 78)*x12*x21    *x41    
     7  +coeff( 79)*x11*x22        *x52
     8  +coeff( 80)*x12*x23            
      x_sp_col    =x_sp_col    
     9  +coeff( 81)*x11*x21    *x43    
     1  +coeff( 82)    *x23    *x41*x52
     2  +coeff( 83)    *x23        *x53
     3  +coeff( 84)*x11        *x42*x52
     4  +coeff( 85)*x11*x23        *x52
     5  +coeff( 86)    *x21    *x43*x52
     6  +coeff( 87)    *x23    *x41*x53
     7  +coeff( 88)*x12*x22    *x42    
     8  +coeff( 89)*x11*x24        *x52
      x_sp_col    =x_sp_col    
     9  +coeff( 90)*x12*x21    *x43    
c
      return
      end
      function t_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 29)
      data ncoeff/ 28/
      data avdat/ -0.7940747E-03/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11199360E-03, 0.45059815E-01,-0.43417178E-02, 0.66507841E-03,
     +  0.16468079E-02, 0.29239592E-02,-0.54701063E-03,-0.30943364E-03,
     +  0.94239437E-03, 0.28425953E-03,-0.13678051E-02, 0.30608851E-03,
     +  0.13620697E-02, 0.30640684E-04, 0.54033815E-04,-0.23536278E-03,
     +  0.16821267E-03, 0.85861335E-04, 0.34786024E-03,-0.19872622E-03,
     +  0.57588306E-04,-0.12328198E-03, 0.45584943E-03, 0.24899992E-03,
     + -0.13203612E-03, 0.87506505E-05, 0.97942429E-04, 0.27265667E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x23            
     6  +coeff(  6)    *x23    *x41    
     7  +coeff(  7)            *x41    
     8  +coeff(  8)*x11                
      t_sp_col    =t_sp_col    
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x23    *x42    
     5  +coeff( 14)                *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x21    *x41*x51
      t_sp_col    =t_sp_col    
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)            *x43*x51
     4  +coeff( 22)            *x42*x52
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x22    *x42*x52
     7  +coeff( 25)            *x42    
     8  +coeff( 26)*x11*x21    *x42    
      t_sp_col    =t_sp_col    
     9  +coeff( 27)*x11*x22    *x42    
     1  +coeff( 28)*x11*x21    *x42*x52
c
      return
      end
      function y_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.8805549E-02/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.16876624E-04, 0.55183716E-01, 0.46256855E-02,-0.13548728E-02,
     + -0.38068090E-02,-0.22121805E-02,-0.45134278E-03, 0.97014237E-03,
     + -0.36950698E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      y_sp_col    =y_sp_col    
     9  +coeff(  9)*x11*x21            
c
      return
      end
      function p_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4142276E-02/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.99295238E-03,-0.15301852E-02, 0.24275590E-01, 0.66562868E-02,
     + -0.44205156E-02,-0.20383941E-02,-0.81859151E-03,-0.24356825E-03,
     +  0.11996370E-02, 0.54814771E-03,-0.44894757E-03,-0.10363871E-02,
     + -0.37700319E-03,-0.22728855E-03, 0.65137638E-03, 0.14558263E-03,
     + -0.76470656E-04,-0.48325982E-03,-0.10149550E-04, 0.22390636E-03,
     + -0.76042699E-04, 0.57015580E-03,-0.74383820E-03, 0.26675794E-03,
     +  0.65552705E-03,-0.13083864E-02,-0.28240051E-04, 0.98277023E-03,
     + -0.18054897E-03,-0.23285307E-02, 0.86822227E-03,-0.64294873E-04,
     + -0.36339171E-03,-0.22236223E-03, 0.23895298E-04, 0.21135331E-03,
     +  0.22865123E-03,-0.18181969E-03, 0.35659800E-03, 0.93156312E-04,
     +  0.65877888E-03, 0.98608485E-04, 0.12341668E-02,-0.23767505E-03,
     + -0.39961172E-03,-0.47963450E-03,-0.13007275E-03,-0.18988832E-03,
     +  0.51805045E-03,-0.28558119E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      p_sp_col    =p_sp_col    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)            *x42    
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)    *x22        *x53
     8  +coeff( 17)            *x41*x51
      p_sp_col    =p_sp_col    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)                *x53
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)*x11            *x52
     4  +coeff( 22)            *x42*x52
     5  +coeff( 23)*x11*x21    *x42    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)            *x43*x52
     8  +coeff( 26)    *x22    *x42*x52
      p_sp_col    =p_sp_col    
     9  +coeff( 27)    *x21    *x42*x53
     1  +coeff( 28)*x11*x23    *x42    
     2  +coeff( 29)*x11        *x42*x53
     3  +coeff( 30)    *x23    *x43*x51
     4  +coeff( 31)*x12*x22    *x43    
     5  +coeff( 32)*x11                
     6  +coeff( 33)    *x21    *x41*x51
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)*x11*x21        *x51
      p_sp_col    =p_sp_col    
     9  +coeff( 36)    *x21    *x42*x51
     1  +coeff( 37)    *x22        *x52
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)*x11*x21    *x41*x51
     4  +coeff( 40)*x12        *x41*x51
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)    *x22    *x42*x51
     7  +coeff( 43)    *x21    *x43*x51
     8  +coeff( 44)            *x42*x53
      p_sp_col    =p_sp_col    
     9  +coeff( 45)*x12*x22    *x41    
     1  +coeff( 46)*x11*x21    *x43*x51
     2  +coeff( 47)*x11*x22    *x41*x52
     3  +coeff( 48)*x11*x21    *x42*x52
     4  +coeff( 49)    *x22    *x42*x53
     5  +coeff( 50)*x11*x21    *x43*x52
c
      return
      end
      function l_sp_col    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 37)
      data ncoeff/ 36/
      data avdat/ -0.4966813E-03/
      data xmin/
     1 -0.39997E-02,-0.51091E-01, 0.00000E+00,-0.25597E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18534E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13066409E-02, 0.50993895E-04,-0.48968908E-02,-0.12383491E-03,
     + -0.29814767E-02,-0.63069927E-03,-0.10579259E-03, 0.25094760E-03,
     +  0.84305721E-04,-0.54829132E-04,-0.82574363E-04, 0.83013838E-04,
     +  0.10381773E-04, 0.23721870E-04,-0.24100341E-05, 0.43254331E-04,
     + -0.84987123E-05,-0.49617600E-04, 0.23211323E-05,-0.30928529E-05,
     +  0.20870939E-05, 0.15275478E-04, 0.12295272E-04, 0.46022110E-05,
     + -0.37346312E-04, 0.37218186E-04, 0.83369687E-05,-0.83029909E-05,
     + -0.26677555E-04,-0.62671766E-05, 0.47956106E-04,-0.33427848E-04,
     + -0.31495685E-04, 0.49246661E-04,-0.19266969E-04, 0.33548720E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_col    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_col    =l_sp_col    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)    *x22    *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11*x21    *x41    
     6  +coeff( 15)                *x52
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)*x11*x22            
      l_sp_col    =l_sp_col    
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)*x11            *x52
     7  +coeff( 25)            *x42*x52
     8  +coeff( 26)*x11*x21    *x42    
      l_sp_col    =l_sp_col    
     9  +coeff( 27)    *x23    *x41*x51
     1  +coeff( 28)    *x21    *x42*x52
     2  +coeff( 29)            *x43*x52
     3  +coeff( 30)    *x22    *x43*x51
     4  +coeff( 31)    *x22    *x42*x52
     5  +coeff( 32)*x11*x23    *x42    
     6  +coeff( 33)*x11*x21    *x41*x53
     7  +coeff( 34)    *x22    *x43*x52
     8  +coeff( 35)*x12*x22    *x43    
      l_sp_col    =l_sp_col    
     9  +coeff( 36)*x11*x23    *x43*x51
c
      return
      end
      function x_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.1178958E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.60442043E-03, 0.13789745E+00,-0.77032307E-02, 0.17741575E-02,
     +  0.34695738E-02, 0.33930277E-02,-0.12408126E-02, 0.24525081E-02,
     +  0.28543314E-02, 0.42339256E-02,-0.33483685E-02, 0.60772354E-03,
     + -0.11905099E-02, 0.29182862E-02,-0.13995435E-03, 0.75045187E-04,
     + -0.46295318E-03, 0.31302913E-03, 0.20103648E-02,-0.98454487E-03,
     +  0.10452260E-02,-0.43171959E-03, 0.88438799E-04, 0.27185696E-03,
     + -0.16568507E-03, 0.39209539E-03,-0.20085127E-03,-0.33850607E-03,
     +  0.72748220E-03, 0.64276257E-02, 0.49996539E-04,-0.21315601E-02,
     + -0.37109881E-03,-0.16072350E-03,-0.11691867E-03, 0.93875373E-04,
     + -0.63899432E-04, 0.96915982E-03,-0.86456863E-03,-0.14161444E-02,
     +  0.88862848E-03, 0.56567957E-03,-0.34033973E-03,-0.34039628E-02,
     + -0.15225678E-02, 0.10463741E-02, 0.21457113E-02,-0.22451193E-02,
     + -0.16254548E-04, 0.99772689E-04,-0.23998928E-05, 0.64453241E-04,
     + -0.17409510E-03,-0.84188861E-04, 0.26898485E-03,-0.15654891E-03,
     +  0.15899873E-02,-0.22181367E-03,-0.38292617E-03, 0.92668397E-05,
     +  0.72716830E-04,-0.30629952E-04, 0.15418225E-02, 0.38689238E-03,
     + -0.94208539E-04,-0.38554909E-03,-0.11880986E-02,-0.80535730E-03,
     + -0.82647090E-03, 0.18269155E-03,-0.19597792E-03, 0.91770315E-03,
     + -0.44305835E-03, 0.33553868E-04,-0.18009649E-04, 0.85202577E-04,
     +  0.12341135E-02,-0.37658142E-03,-0.47908456E-03, 0.44652046E-04,
     +  0.96797559E-03, 0.42325849E-03,-0.14555685E-02,-0.51163521E-03,
     +  0.29650286E-02,-0.65756415E-03,-0.94679015E-03,-0.43113131E-03,
     +  0.22221196E-02, 0.36627731E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23            
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22            
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)    *x23    *x42    
     6  +coeff( 15)*x11        *x43    
     7  +coeff( 16)                *x51
     8  +coeff( 17)            *x42    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x24    *x41    
     3  +coeff( 21)*x11*x22    *x42    
     4  +coeff( 22)*x11        *x41    
     5  +coeff( 23)*x11            *x51
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)*x11*x21    *x41    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)            *x41*x51
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x24    *x41    
     7  +coeff( 34)            *x43    
     8  +coeff( 35)    *x22    *x41*x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 36)*x12*x21            
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)    *x22    *x43    
     3  +coeff( 39)*x11*x23    *x41    
     4  +coeff( 40)    *x24    *x42    
     5  +coeff( 41)    *x23    *x42*x51
     6  +coeff( 42)*x11*x23    *x42    
     7  +coeff( 43)*x11*x22    *x42*x51
     8  +coeff( 44)*x12*x24    *x41    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 45)*x12*x22    *x43    
     1  +coeff( 46)*x11*x24    *x42*x51
     2  +coeff( 47)*x12*x24    *x41*x51
     3  +coeff( 48)*x11*x22    *x43*x52
     4  +coeff( 49)                *x53
     5  +coeff( 50)*x11*x22        *x51
     6  +coeff( 51)*x11*x21    *x41*x51
     7  +coeff( 52)*x12*x22            
     8  +coeff( 53)    *x21    *x43*x51
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 54)    *x24        *x52
     1  +coeff( 55)*x11*x21    *x41*x52
     2  +coeff( 56)    *x23        *x53
     3  +coeff( 57)*x12*x22    *x41    
     4  +coeff( 58)*x11*x24        *x51
     5  +coeff( 59)    *x22    *x42*x52
     6  +coeff( 60)*x12*x21    *x42    
     7  +coeff( 61)*x11*x23    *x41*x51
     8  +coeff( 62)*x12        *x43    
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 63)*x11*x22    *x43    
     1  +coeff( 64)*x11*x22    *x41*x52
     2  +coeff( 65)    *x24    *x41*x52
     3  +coeff( 66)*x12*x23        *x51
     4  +coeff( 67)    *x23    *x43*x51
     5  +coeff( 68)*x11*x24    *x42    
     6  +coeff( 69)*x12*x22    *x41*x51
     7  +coeff( 70)*x11        *x43*x52
     8  +coeff( 71)    *x22    *x42*x53
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 72)*x11*x23    *x43    
     1  +coeff( 73)*x11*x23    *x42*x51
     2  +coeff( 74)*x12*x21        *x53
     3  +coeff( 75)*x11*x23        *x53
     4  +coeff( 76)*x12        *x43*x51
     5  +coeff( 77)    *x24    *x42*x52
     6  +coeff( 78)*x12*x23    *x42    
     7  +coeff( 79)*x11*x21    *x43*x52
     8  +coeff( 80)*x12        *x42*x53
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 81)*x12*x23    *x42*x51
     1  +coeff( 82)*x12*x23        *x53
     2  +coeff( 83)*x11*x23    *x43*x52
     3  +coeff( 84)*x12*x21    *x42*x53
     4  +coeff( 85)*x12*x24    *x43    
     5  +coeff( 86)*x12*x24    *x42*x51
     6  +coeff( 87)*x12*x23    *x41*x53
     7  +coeff( 88)*x12*x22    *x42*x53
     8  +coeff( 89)*x12*x23    *x43*x53
      x_sp_q1ex   =x_sp_q1ex   
     9  +coeff( 90)            *x42*x51
c
      return
      end
      function t_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.2296385E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.17086115E-03,-0.11445392E-01,-0.16230084E-02,-0.17007065E-02,
     +  0.27853509E-02,-0.20392342E-03, 0.23055117E-03, 0.33846172E-03,
     +  0.31141509E-03,-0.77067123E-03, 0.12253247E-02, 0.12208959E-04,
     +  0.70921997E-04,-0.16538042E-03,-0.14599014E-03, 0.79462967E-04,
     +  0.58499090E-05, 0.40890114E-04, 0.26085609E-03, 0.11892527E-03,
     +  0.17654149E-05, 0.20916294E-03, 0.46330187E-03, 0.20122397E-03,
     + -0.78510733E-04,-0.84984174E-04,-0.61102379E-04, 0.30399580E-04,
     +  0.33406297E-04,-0.36804377E-04, 0.17748460E-03, 0.99636638E-03,
     +  0.14431231E-04,-0.29219029E-03, 0.41039661E-04,-0.28588318E-04,
     +  0.53267904E-04, 0.60911389E-04,-0.39998005E-04,-0.92564595E-04,
     + -0.17278839E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)*x11        *x43    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x23        *x51
     3  +coeff( 21)            *x42*x52
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)            *x42    
     8  +coeff( 26)*x11        *x41    
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 27)            *x43    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)*x11        *x42    
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)*x12*x21            
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)    *x21    *x42*x51
      t_sp_q1ex   =t_sp_q1ex   
     9  +coeff( 36)*x11*x23            
     1  +coeff( 37)    *x22    *x42*x51
     2  +coeff( 38)*x11*x21    *x43    
     3  +coeff( 39)*x11*x22    *x41*x52
     4  +coeff( 40)*x11*x23    *x41*x52
     5  +coeff( 41)*x11*x22    *x42*x52
c
      return
      end
      function y_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/ -0.9726789E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.64733545E-02, 0.12335366E+00, 0.20201284E-01,-0.50231484E-02,
     +  0.31878429E-02,-0.13822552E-01,-0.39442605E-02,-0.18338546E-02,
     + -0.12011703E-02,-0.34992825E-02, 0.76466275E-03,-0.17884821E-02,
     + -0.11384519E-02, 0.37266628E-02,-0.98684700E-02, 0.90176164E-03,
     +  0.35941552E-02,-0.19222457E-02,-0.76694121E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)            *x41*x51
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21    *x41    
     2  +coeff( 11)            *x43    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x23    *x41    
      y_sp_q1ex   =y_sp_q1ex   
     9  +coeff( 18)*x11*x21    *x42    
     1  +coeff( 19)    *x22    *x43    
c
      return
      end
      function p_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.4125908E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.23677766E-02,-0.25220478E-02, 0.52058153E-01, 0.10646850E-01,
     + -0.69617056E-02,-0.16357367E-02,-0.61652693E-03, 0.18543025E-02,
     + -0.25028684E-02,-0.11401881E-02, 0.15348699E-02,-0.80809987E-03,
     +  0.11193370E-03,-0.49845180E-02,-0.10910979E-02, 0.24923493E-03,
     + -0.48634896E-03,-0.26541881E-02,-0.10869525E-02, 0.51932823E-03,
     +  0.37741080E-04,-0.84595702E-03, 0.12741069E-02, 0.14122792E-02,
     + -0.38763466E-02,-0.10161002E-03, 0.21862591E-03,-0.26612106E-03,
     +  0.25396177E-03, 0.98540600E-04,-0.28696188E-03,-0.26905007E-03,
     + -0.20565759E-03, 0.88942534E-03,-0.40864304E-03, 0.17419183E-02,
     + -0.19268335E-03,-0.49333815E-02,-0.19061794E-02,-0.20199455E-04,
     +  0.29419398E-05,-0.91826463E-04, 0.24583982E-03, 0.25173009E-03,
     + -0.13910454E-03, 0.78641780E-03, 0.16042786E-02, 0.86963858E-03,
     +  0.29880987E-03,-0.46526137E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q1ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x23            
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x43    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)                *x52
     4  +coeff( 13)            *x43    
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x21    *x41    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)    *x21    *x42    
     2  +coeff( 20)    *x23    *x41    
     3  +coeff( 21)*x11*x23            
     4  +coeff( 22)*x11*x21    *x42    
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)*x11*x22    *x41*x52
     7  +coeff( 25)*x11*x23    *x43*x51
     8  +coeff( 26)*x11                
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 27)    *x21        *x51
     1  +coeff( 28)            *x41*x52
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)    *x23        *x51
     5  +coeff( 32)*x11        *x41*x52
     6  +coeff( 33)            *x42*x53
     7  +coeff( 34)*x11*x23    *x41    
     8  +coeff( 35)*x12*x21    *x42    
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 36)    *x23    *x43    
     1  +coeff( 37)*x12*x22    *x41*x51
     2  +coeff( 38)    *x22    *x43*x52
     3  +coeff( 39)*x12*x22    *x42*x51
     4  +coeff( 40)*x11            *x51
     5  +coeff( 41)*x12                
     6  +coeff( 42)*x11            *x52
     7  +coeff( 43)            *x42*x52
     8  +coeff( 44)*x11*x21    *x41*x51
      p_sp_q1ex   =p_sp_q1ex   
     9  +coeff( 45)*x12*x22            
     1  +coeff( 46)    *x22    *x42*x51
     2  +coeff( 47)    *x22    *x41*x52
     3  +coeff( 48)            *x43*x52
     4  +coeff( 49)    *x22        *x53
     5  +coeff( 50)*x11*x21    *x43    
c
      return
      end
      function l_sp_q1ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1097887E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17678315E-02,-0.52619924E-06,-0.43437011E-02,-0.35442226E-02,
     +  0.18778631E-03,-0.21435309E-02,-0.77088195E-03, 0.59150683E-03,
     + -0.83512314E-04, 0.37717840E-03,-0.66443499E-04, 0.14782736E-03,
     + -0.10223038E-03,-0.53775770E-04, 0.40740149E-04, 0.36890101E-05,
     +  0.16418136E-03,-0.11351909E-03, 0.60676222E-04, 0.81221195E-04,
     +  0.39361792E-04, 0.63180318E-03,-0.88863535E-05, 0.46329178E-05,
     +  0.85588617E-04,-0.15622686E-03,-0.18005622E-05, 0.21200765E-03,
     +  0.26038216E-03, 0.24709610E-04, 0.16993930E-03,-0.24862431E-04,
     + -0.64744454E-04,-0.46957430E-04, 0.10297156E-03,-0.34506404E-03,
     + -0.22432543E-04,-0.12767279E-03,-0.20309875E-04, 0.48371377E-04,
     +  0.25191413E-04, 0.58661462E-04, 0.23687142E-04, 0.25893783E-03,
     + -0.87889683E-04,-0.19576932E-03, 0.10384072E-04,-0.43972900E-05,
     +  0.62245695E-05,-0.22878581E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)                *x51
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x21    *x42    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 18)            *x43    
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)            *x41*x52
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)*x11*x22            
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)*x11*x21    *x42    
     8  +coeff( 26)    *x23    *x42    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 27)*x11*x23    *x41    
     1  +coeff( 28)*x11*x23    *x42*x51
     2  +coeff( 29)*x12*x22    *x42*x51
     3  +coeff( 30)                *x53
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)    *x23        *x51
     6  +coeff( 33)            *x43*x52
     7  +coeff( 34)    *x22        *x53
     8  +coeff( 35)*x11*x21    *x43    
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 36)    *x23    *x43    
     1  +coeff( 37)*x12*x22        *x51
     2  +coeff( 38)*x11*x23    *x42    
     3  +coeff( 39)*x11*x22    *x41*x52
     4  +coeff( 40)*x11        *x43*x52
     5  +coeff( 41)*x11        *x42*x53
     6  +coeff( 42)*x12*x21    *x43    
     7  +coeff( 43)    *x23    *x43*x51
     8  +coeff( 44)    *x22    *x43*x52
      l_sp_q1ex   =l_sp_q1ex   
     9  +coeff( 45)*x12*x22    *x41*x52
     1  +coeff( 46)*x11*x22    *x43*x52
     2  +coeff( 47)    *x21    *x41*x51
     3  +coeff( 48)    *x21        *x52
     4  +coeff( 49)*x12        *x41    
     5  +coeff( 50)    *x22        *x52
c
      return
      end
      function x_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2799580E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.30031707E-03, 0.18376021E+00,-0.31237402E-02,-0.18908719E-01,
     +  0.11702695E-01, 0.55474397E-02, 0.68165972E-02,-0.50728652E-03,
     +  0.11721591E-01, 0.42695282E-02,-0.28177719E-02,-0.10338395E-02,
     +  0.64015216E-02,-0.86730253E-02, 0.13593561E-02,-0.27397338E-02,
     +  0.23337854E-03, 0.72300708E-03,-0.10392672E-02, 0.40014583E-03,
     + -0.82490867E-03, 0.18993375E-02,-0.20366558E-02, 0.53801420E-02,
     +  0.84977743E-03,-0.78655372E-03, 0.85344305E-04,-0.16966078E-03,
     +  0.28724503E-02,-0.27904178E-02, 0.13915183E-01,-0.46602741E-04,
     + -0.47375390E-03,-0.93074734E-04, 0.23619678E-03,-0.44840691E-02,
     +  0.30579575E-03,-0.28113637E-03,-0.18252867E-02, 0.18400744E-02,
     +  0.19478707E-02,-0.39672389E-03, 0.56797033E-03,-0.28128116E-03,
     + -0.49712523E-02,-0.82889171E-02,-0.42303945E-02,-0.37848481E-03,
     + -0.31099387E-02, 0.37228945E-03, 0.84525404E-04,-0.58021698E-04,
     +  0.14601945E-03,-0.14841123E-03, 0.72893490E-04, 0.67748821E-04,
     +  0.76424330E-03, 0.19355552E-03, 0.15125390E-03,-0.13124827E-04,
     +  0.68395439E-03, 0.45758896E-02,-0.79557428E-03, 0.13557334E-03,
     + -0.19005944E-03, 0.29273593E-03, 0.40899757E-02,-0.21650372E-02,
     +  0.28635929E-02,-0.15274342E-02,-0.28723417E-03, 0.11763999E-02,
     + -0.11958897E-02,-0.97347563E-03,-0.34443936E-02, 0.93599170E-04,
     + -0.13391091E-02, 0.76123991E-03, 0.35450472E-02,-0.52682091E-02,
     +  0.63799392E-02, 0.25015092E-02, 0.22000629E-02, 0.22178318E-02,
     +  0.27016702E-04, 0.83094550E-04, 0.22604031E-03,-0.76779281E-04,
     + -0.25342978E-03, 0.45327284E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23            
     8  +coeff(  8)            *x43    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff(  9)    *x23    *x41    
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)            *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x24            
     8  +coeff( 17)                *x51
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x23            
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 27)            *x43*x51
     1  +coeff( 28)*x11*x24            
     2  +coeff( 29)    *x22    *x43    
     3  +coeff( 30)    *x24    *x42    
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)    *x21    *x41*x51
     6  +coeff( 33)*x11        *x42    
     7  +coeff( 34)    *x22        *x52
     8  +coeff( 35)*x12*x21            
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)*x11        *x43    
     3  +coeff( 39)*x11*x23    *x41    
     4  +coeff( 40)*x11*x22    *x42    
     5  +coeff( 41)    *x23    *x42*x51
     6  +coeff( 42)*x11*x24    *x41    
     7  +coeff( 43)*x11*x23    *x42    
     8  +coeff( 44)*x11*x22    *x42*x51
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 45)    *x23    *x43*x51
     1  +coeff( 46)*x12*x24    *x41    
     2  +coeff( 47)*x12*x22    *x43    
     3  +coeff( 48)*x12*x21    *x42*x52
     4  +coeff( 49)*x11*x22    *x43*x52
     5  +coeff( 50)*x12*x24        *x53
     6  +coeff( 51)            *x41*x51
     7  +coeff( 52)                *x53
     8  +coeff( 53)*x11*x21        *x51
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 54)    *x22    *x41*x51
     1  +coeff( 55)    *x21    *x41*x52
     2  +coeff( 56)*x11*x21    *x41*x51
     3  +coeff( 57)    *x23    *x41*x51
     4  +coeff( 58)*x11        *x41*x52
     5  +coeff( 59)    *x22    *x41*x52
     6  +coeff( 60)            *x42*x53
     7  +coeff( 61)*x11*x21    *x41*x52
     8  +coeff( 62)*x12*x22    *x41    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 63)*x12*x22        *x51
     1  +coeff( 64)    *x22    *x41*x53
     2  +coeff( 65)*x11*x23    *x41*x51
     3  +coeff( 66)*x12        *x43    
     4  +coeff( 67)*x11*x22    *x43    
     5  +coeff( 68)    *x23    *x41*x53
     6  +coeff( 69)*x11*x23    *x43    
     7  +coeff( 70)*x11*x23    *x42*x51
     8  +coeff( 71)*x11*x23        *x53
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 72)*x12*x24        *x51
     1  +coeff( 73)*x12*x23    *x42    
     2  +coeff( 74)*x11*x21    *x43*x52
     3  +coeff( 75)*x11*x24    *x43    
     4  +coeff( 76)*x11        *x43*x53
     5  +coeff( 77)    *x22    *x43*x53
     6  +coeff( 78)*x11*x23    *x42*x52
     7  +coeff( 79)    *x23    *x43*x53
     8  +coeff( 80)*x11*x23    *x43*x52
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 81)*x12*x24    *x43    
     1  +coeff( 82)*x11*x24    *x42*x53
     2  +coeff( 83)*x12*x23    *x43*x52
     3  +coeff( 84)*x12*x23    *x43*x53
     4  +coeff( 85)*x12                
     5  +coeff( 86)            *x42*x51
     6  +coeff( 87)    *x23        *x51
     7  +coeff( 88)*x11        *x41*x51
     8  +coeff( 89)*x12        *x41    
      x_sp_q2ex   =x_sp_q2ex   
     9  +coeff( 90)*x12            *x51
c
      return
      end
      function t_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.9147560E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.41209649E-04, 0.43665167E-01,-0.18765669E-02,-0.62038233E-02,
     +  0.30942238E-02, 0.19976443E-02,-0.80104510E-03, 0.10477643E-02,
     +  0.10919506E-02,-0.32140969E-02, 0.40172329E-02, 0.77714743E-04,
     + -0.38007082E-03, 0.16812912E-03,-0.34956104E-03, 0.45544389E-03,
     +  0.18030654E-02, 0.19712944E-02,-0.36824640E-03, 0.12063439E-03,
     + -0.28395082E-03, 0.12697576E-02, 0.60203130E-03, 0.38113147E-04,
     +  0.29538266E-03,-0.47896501E-04, 0.45383000E-02,-0.13275849E-03,
     +  0.58118796E-04,-0.13293105E-02, 0.18159770E-03, 0.23280240E-03,
     +  0.22862603E-03,-0.27184759E-03, 0.17934522E-03,-0.15674681E-03,
     + -0.12182838E-03, 0.84228435E-03,-0.65241975E-03,-0.28237476E-03,
     +  0.36770670E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)    *x23            
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22            
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)                *x51
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)*x11            *x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x23    *x42    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 18)    *x22    *x43    
     1  +coeff( 19)            *x42    
     2  +coeff( 20)*x11*x21            
     3  +coeff( 21)            *x43    
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)            *x41*x51
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x23            
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 27)    *x23    *x43    
     1  +coeff( 28)    *x21    *x41*x51
     2  +coeff( 29)*x12*x21            
     3  +coeff( 30)    *x21    *x43    
     4  +coeff( 31)    *x21    *x42*x51
     5  +coeff( 32)*x11*x21    *x42    
     6  +coeff( 33)*x11        *x41*x52
     7  +coeff( 34)*x11*x23    *x41    
     8  +coeff( 35)*x11*x22    *x42    
      t_sp_q2ex   =t_sp_q2ex   
     9  +coeff( 36)*x11*x21    *x41*x52
     1  +coeff( 37)*x11        *x42*x52
     2  +coeff( 38)*x11*x22    *x43    
     3  +coeff( 39)*x11*x22    *x41*x52
     4  +coeff( 40)*x11        *x43*x52
     5  +coeff( 41)*x12*x22    *x42*x51
c
      return
      end
      function y_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 16)
      data ncoeff/ 15/
      data avdat/ -0.1351821E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.87683955E-02, 0.17318302E+00, 0.32577906E-01,-0.89050047E-02,
     +  0.78803683E-02,-0.21249520E-01,-0.21268055E-02,-0.99063050E-02,
     + -0.36304118E-02, 0.72962800E-02,-0.21331929E-01,-0.13281680E-02,
     + -0.23829478E-02,-0.31766489E-02,-0.15710382E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x44 = x43*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      y_sp_q2ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)    *x22            
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x22    *x41    
      y_sp_q2ex   =y_sp_q2ex   
     9  +coeff(  9)    *x21    *x43    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)            *x44    
     6  +coeff( 15)*x11*x21    *x41    
c
      return
      end
      function p_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1620335E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10316764E-02, 0.57135249E-03,-0.19113310E-01,-0.25831643E-02,
     +  0.16602469E-02, 0.30121480E-02, 0.67648443E-03,-0.43228475E-03,
     +  0.24437727E-03, 0.15443469E-03,-0.39867542E-03, 0.41431689E-04,
     + -0.32827919E-03,-0.17622176E-03,-0.20211123E-03, 0.11890146E-02,
     +  0.37953482E-03, 0.19097277E-03,-0.53636904E-04, 0.15130373E-03,
     +  0.26041432E-05,-0.26648177E-04, 0.10419448E-02, 0.11249373E-04,
     + -0.73952069E-04,-0.24908417E-04, 0.31889434E-04, 0.30972678E-03,
     + -0.23745105E-03,-0.10763758E-03,-0.24631237E-04,-0.76460412E-04,
     + -0.59291691E-04, 0.19542508E-03,-0.37421277E-03,-0.57552557E-03,
     + -0.26286216E-03,-0.69747859E-03, 0.39054707E-03, 0.74176473E-03,
     +  0.37444952E-04, 0.57807491E-04,-0.11511111E-03, 0.17714254E-04,
     +  0.13733302E-03,-0.12250640E-03, 0.87966597E-04,-0.88323068E-04,
     + -0.15741210E-03, 0.18207208E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x22        *x51
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)    *x23            
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x21        *x51
     6  +coeff( 15)            *x41*x52
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)    *x21    *x43    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x23        *x51
     3  +coeff( 21)    *x22    *x41*x51
     4  +coeff( 22)*x11*x23            
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)    *x21    *x42*x52
     7  +coeff( 25)            *x43*x52
     8  +coeff( 26)    *x22    *x42*x53
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 27)*x11                
     1  +coeff( 28)    *x21    *x42    
     2  +coeff( 29)            *x43    
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)                *x53
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)*x11*x21    *x42    
     8  +coeff( 35)    *x23    *x42    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 36)    *x22    *x42*x51
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)    *x23    *x43    
     3  +coeff( 39)*x11*x23    *x41*x51
     4  +coeff( 40)*x12*x22    *x42*x51
     5  +coeff( 41)    *x22        *x52
     6  +coeff( 42)            *x41*x53
     7  +coeff( 43)*x11*x21    *x41*x51
     8  +coeff( 44)*x12*x21    *x41    
      p_sp_q2ex   =p_sp_q2ex   
     9  +coeff( 45)    *x23    *x41*x51
     1  +coeff( 46)    *x22        *x53
     2  +coeff( 47)*x12*x21    *x42    
     3  +coeff( 48)*x12*x22        *x51
     4  +coeff( 49)    *x22    *x43*x51
     5  +coeff( 50)    *x22    *x42*x52
c
      return
      end
      function l_sp_q2ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/ -0.1724131E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23768211E-02, 0.10288813E-04,-0.39478862E-02,-0.46742936E-02,
     +  0.31390658E-03,-0.36260078E-02,-0.13866260E-02,-0.14494325E-03,
     +  0.10027507E-02, 0.14889907E-03, 0.16218548E-04, 0.68680674E-04,
     + -0.93907656E-04, 0.63523900E-03, 0.74316595E-04, 0.13125628E-03,
     +  0.84790423E-04,-0.18814454E-03, 0.11849107E-02,-0.43862168E-04,
     +  0.26184105E-03,-0.22132449E-03, 0.10289949E-03, 0.12716031E-04,
     + -0.25232701E-03,-0.18377139E-03, 0.32368605E-03,-0.95365958E-05,
     +  0.28863145E-03,-0.51855091E-04, 0.88645451E-04, 0.37636077E-04,
     + -0.18510127E-04,-0.70566188E-04,-0.40240811E-05,-0.59041451E-03,
     +  0.90512745E-04, 0.35788002E-03, 0.11787856E-04, 0.21402308E-04,
     + -0.11259644E-04, 0.41806587E-04, 0.76613296E-05,-0.47230067E-04,
     + -0.29509843E-04,-0.25299796E-04, 0.73840944E-04, 0.88246103E-04,
     + -0.77179138E-04, 0.11638137E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)                *x52
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff(  9)    *x22    *x41    
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)                *x51
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x42    
     6  +coeff( 15)    *x21        *x51
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)*x11*x21    *x41    
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)    *x22    *x43    
     2  +coeff( 20)            *x43*x52
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)            *x43    
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)*x11*x22    *x41*x52
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 27)*x12*x22    *x42*x51
     1  +coeff( 28)*x11                
     2  +coeff( 29)    *x21    *x43    
     3  +coeff( 30)    *x23        *x51
     4  +coeff( 31)*x11*x21    *x42    
     5  +coeff( 32)*x11        *x41*x52
     6  +coeff( 33)*x12        *x42    
     7  +coeff( 34)*x11*x23    *x41    
     8  +coeff( 35)*x11*x21    *x42*x51
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 36)    *x23    *x43    
     1  +coeff( 37)*x12*x21    *x43    
     2  +coeff( 38)*x11*x23    *x42*x51
     3  +coeff( 39)*x11            *x51
     4  +coeff( 40)    *x21    *x41*x51
     5  +coeff( 41)    *x21        *x52
     6  +coeff( 42)                *x53
     7  +coeff( 43)*x11            *x52
     8  +coeff( 44)    *x22        *x52
      l_sp_q2ex   =l_sp_q2ex   
     9  +coeff( 45)            *x42*x52
     1  +coeff( 46)*x11*x22        *x51
     2  +coeff( 47)    *x22    *x42*x51
     3  +coeff( 48)    *x22    *x41*x52
     4  +coeff( 49)    *x22        *x53
     5  +coeff( 50)*x11*x21    *x43    
c
      return
      end
      function x_sp_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5097530E+01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69499807E-03,-0.10571390E+00, 0.30871115E-02, 0.12849805E-01,
     + -0.69489330E-02, 0.19715701E-02,-0.43440107E-02,-0.52349479E-02,
     +  0.62422631E-02,-0.81493137E-02,-0.61762112E-03,-0.83780556E-03,
     +  0.20382439E-02, 0.23262178E-03, 0.61069598E-03,-0.31460606E-03,
     +  0.65217813E-03,-0.86024258E-04,-0.31945196E-02,-0.33362249E-02,
     + -0.11949282E-02, 0.11582772E-03,-0.11577542E-02, 0.67257154E-03,
     +  0.20542603E-02,-0.70692436E-03, 0.30071431E-03, 0.45684527E-03,
     + -0.16030094E-02,-0.93866978E-02,-0.21167405E-03,-0.15019755E-03,
     +  0.38140509E-03,-0.16947386E-03, 0.28824634E-02,-0.44597851E-03,
     +  0.12547377E-02, 0.26501340E-03,-0.39755451E-03,-0.15146365E-03,
     +  0.43763773E-03,-0.16628893E-02, 0.12027519E-02, 0.98891545E-03,
     +  0.19432332E-02, 0.96909789E-04,-0.73378196E-03,-0.66678499E-05,
     + -0.13172836E-03,-0.12587713E-02, 0.15546411E-04, 0.42259321E-02,
     + -0.96810708E-03, 0.99576231E-04, 0.25415921E-03, 0.35498049E-02,
     +  0.68471345E-04,-0.53804837E-04, 0.39652019E-03,-0.88758410E-04,
     +  0.26695645E-04,-0.15789429E-03,-0.42118867E-04, 0.34536817E-04,
     +  0.49191993E-04, 0.60263075E-04,-0.24422107E-03,-0.29671166E-03,
     +  0.13815094E-03,-0.19181902E-02, 0.56566880E-03, 0.14374753E-03,
     + -0.36223566E-04,-0.12950100E-03, 0.10550899E-03,-0.19671589E-02,
     + -0.68587641E-03,-0.40758701E-03,-0.10481178E-03, 0.13155093E-03,
     +  0.31033374E-03, 0.84908510E-03,-0.31882705E-03, 0.18592909E-03,
     +  0.72863436E-03,-0.19628328E-03, 0.12489119E-03,-0.60167362E-03,
     +  0.85483619E-03, 0.56971330E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      x_sp_den    =x_sp_den    
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24    *x41    
     5  +coeff( 14)*x11        *x43    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21        *x52
      x_sp_den    =x_sp_den    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)*x11*x22    *x42    
     4  +coeff( 22)    *x24        *x52
     5  +coeff( 23)    *x22            
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x24            
     8  +coeff( 26)*x11*x21    *x41    
      x_sp_den    =x_sp_den    
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)                *x51
     5  +coeff( 32)            *x41*x51
     6  +coeff( 33)    *x22        *x51
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)    *x21    *x43    
      x_sp_den    =x_sp_den    
     9  +coeff( 36)    *x21    *x42*x51
     1  +coeff( 37)*x11*x24    *x41    
     2  +coeff( 38)            *x43    
     3  +coeff( 39)    *x23        *x51
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x23    *x41    
     8  +coeff( 44)    *x21    *x43*x51
      x_sp_den    =x_sp_den    
     9  +coeff( 45)    *x24    *x42    
     1  +coeff( 46)*x11*x21    *x43    
     2  +coeff( 47)    *x23    *x42*x51
     3  +coeff( 48)*x11        *x42*x52
     4  +coeff( 49)*x11*x22    *x42*x51
     5  +coeff( 50)*x11*x23    *x43    
     6  +coeff( 51)*x11*x23    *x41*x52
     7  +coeff( 52)*x12*x24    *x41    
     8  +coeff( 53)*x12*x24        *x51
      x_sp_den    =x_sp_den    
     9  +coeff( 54)*x12*x22    *x43    
     1  +coeff( 55)*x12*x21    *x42*x52
     2  +coeff( 56)*x11*x22    *x43*x52
     3  +coeff( 57)                *x53
     4  +coeff( 58)*x11*x21        *x51
     5  +coeff( 59)    *x22    *x41*x51
     6  +coeff( 60)    *x21    *x41*x52
     7  +coeff( 61)*x12            *x51
     8  +coeff( 62)*x11*x21    *x42    
      x_sp_den    =x_sp_den    
     9  +coeff( 63)*x11*x21    *x41*x51
     1  +coeff( 64)*x11        *x42*x51
     2  +coeff( 65)*x11        *x41*x52
     3  +coeff( 66)    *x22    *x41*x52
     4  +coeff( 67)*x12*x21    *x41    
     5  +coeff( 68)    *x24    *x41*x51
     6  +coeff( 69)    *x23        *x53
     7  +coeff( 70)*x12*x22    *x41    
     8  +coeff( 71)*x12*x22        *x51
      x_sp_den    =x_sp_den    
     9  +coeff( 72)*x11*x23    *x41*x51
     1  +coeff( 73)*x12*x21        *x52
     2  +coeff( 74)*x12*x24            
     3  +coeff( 75)*x12        *x43    
     4  +coeff( 76)*x11*x22    *x43    
     5  +coeff( 77)*x11*x22    *x41*x52
     6  +coeff( 78)    *x24    *x41*x52
     7  +coeff( 79)*x12            *x53
     8  +coeff( 80)*x12*x23    *x41    
      x_sp_den    =x_sp_den    
     9  +coeff( 81)*x12*x23        *x51
     1  +coeff( 82)    *x23    *x41*x53
     2  +coeff( 83)*x11        *x43*x52
     3  +coeff( 84)*x12*x21    *x43    
     4  +coeff( 85)*x11*x23    *x42*x51
     5  +coeff( 86)*x12*x21        *x53
     6  +coeff( 87)*x11*x23        *x53
     7  +coeff( 88)    *x21    *x43*x53
     8  +coeff( 89)*x12*x23    *x42    
      x_sp_den    =x_sp_den    
     9  +coeff( 90)*x11*x21    *x43*x52
c
      return
      end
      function t_sp_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/  0.3603638E+01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.41653894E-01, 0.61587846E+00, 0.75405012E-02,-0.26042875E-01,
     +  0.17275757E+00,-0.94125852E-01, 0.43867867E-01, 0.69286987E-01,
     + -0.10068014E-02,-0.17872421E-01, 0.12134173E-01, 0.80567542E-02,
     + -0.69379732E-02,-0.44120770E-01, 0.18911004E-01, 0.58756132E-01,
     + -0.47258181E-02, 0.24379414E-01,-0.42992570E-02, 0.21755805E-02,
     +  0.26804563E-02, 0.86496538E-02, 0.99820830E-02,-0.12971380E-02,
     +  0.45698094E-02, 0.19297046E-02,-0.57454794E-02,-0.31754598E-02,
     +  0.27595987E-02, 0.65209411E-01,-0.21046434E-01, 0.38078502E-02,
     + -0.15565638E-02,-0.42118491E-02,-0.40964084E-02, 0.56080287E-03,
     +  0.61563626E-02,-0.23743336E-02, 0.63098604E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      t_sp_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      t_sp_den    =t_sp_den    
     9  +coeff(  9)            *x43    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21        *x52
      t_sp_den    =t_sp_den    
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)*x11*x21    *x41    
      t_sp_den    =t_sp_den    
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)            *x43*x52
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)    *x22        *x52
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x22    *x42*x51
      t_sp_den    =t_sp_den    
     9  +coeff( 36)*x12            *x52
     1  +coeff( 37)*x11*x22    *x42    
     2  +coeff( 38)*x11        *x42*x52
     3  +coeff( 39)*x11*x23    *x43    
c
      return
      end
      function y_sp_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/ -0.7348734E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.44901045E-02, 0.99475138E-01, 0.22400493E-01,-0.63730804E-02,
     +  0.56264075E-02, 0.11914965E-01,-0.26342151E-03,-0.15184742E-01,
     + -0.14983970E-02,-0.90736682E-02, 0.41951705E-02, 0.13742591E-02,
     + -0.13790208E-01,-0.81911578E-03,-0.19263171E-02,-0.21686191E-02,
     + -0.13823749E-02,-0.10150283E-02, 0.37228724E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22            
      y_sp_den    =y_sp_den    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)                *x54
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x44    
     8  +coeff( 17)    *x22        *x51
      y_sp_den    =y_sp_den    
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x23    *x41    
c
      return
      end
      function p_sp_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.5418017E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25268176E-02, 0.38742465E-02,-0.55461004E-01,-0.59286216E-02,
     +  0.61736177E-02,-0.18566305E-01,-0.29717816E-04,-0.37588002E-02,
     +  0.12132913E-01, 0.24796359E-02,-0.94683133E-02,-0.36140098E-02,
     +  0.32209881E-03, 0.73140085E-03, 0.32282320E-02,-0.10602347E-02,
     +  0.36975616E-02, 0.14668208E-02,-0.10758217E-02,-0.34381705E-03,
     +  0.96127979E-03,-0.17082912E-03,-0.74099359E-03, 0.16265608E-03,
     + -0.30013444E-03, 0.21756240E-02,-0.23243048E-02, 0.15802846E-02,
     + -0.13693585E-02, 0.10605219E-03,-0.86288404E-04,-0.10180486E-03,
     +  0.37663389E-03, 0.14936439E-03,-0.10384311E-03,-0.25409847E-03,
     + -0.94329217E-03,-0.25378354E-02, 0.78017532E-04,-0.70071736E-04,
     +  0.19229551E-03,-0.16741667E-04, 0.39735064E-03, 0.17907593E-03,
     + -0.25196481E-03,-0.11677594E-02, 0.30445692E-03,-0.47576541E-03,
     +  0.53170905E-03, 0.71945647E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_den    =p_sp_den    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)    *x22    *x42    
      p_sp_den    =p_sp_den    
     9  +coeff( 18)    *x23            
     1  +coeff( 19)            *x43    
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)                *x53
     8  +coeff( 26)    *x21    *x43    
      p_sp_den    =p_sp_den    
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)    *x22    *x42*x51
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)*x11*x22            
     5  +coeff( 32)*x11        *x41*x51
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)    *x22        *x52
     8  +coeff( 35)*x11*x23            
      p_sp_den    =p_sp_den    
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)    *x23    *x43    
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)*x11*x21        *x51
     5  +coeff( 41)    *x21    *x41*x52
     6  +coeff( 42)            *x41*x53
     7  +coeff( 43)*x11*x21    *x42    
     8  +coeff( 44)            *x42*x53
      p_sp_den    =p_sp_den    
     9  +coeff( 45)*x12*x22    *x41    
     1  +coeff( 46)    *x22    *x43*x51
     2  +coeff( 47)            *x43*x53
     3  +coeff( 48)*x11*x22    *x41*x52
     4  +coeff( 49)*x12*x22    *x41*x51
     5  +coeff( 50)    *x22    *x43*x52
c
      return
      end
      function l_sp_den    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 50)
      data ncoeff/ 49/
      data avdat/ -0.2278157E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51585832E-02, 0.20447290E+00,-0.68634884E-02,-0.59955809E-02,
     + -0.24418786E-01, 0.13927435E-01,-0.95742205E-02, 0.82017137E-02,
     + -0.56372047E-02, 0.66158632E-02,-0.12146716E-01, 0.15258301E-01,
     + -0.14507415E-02, 0.11037112E-02, 0.14679035E-02, 0.64611072E-02,
     + -0.30588533E-03, 0.77636275E-02, 0.65789890E-03, 0.58396929E-03,
     + -0.12920818E-02, 0.89207236E-02,-0.32589567E-04, 0.33565145E-02,
     + -0.21015985E-03,-0.12947259E-02,-0.14255648E-02,-0.60919626E-03,
     +  0.21124879E-03, 0.14598644E-02,-0.56535110E-03, 0.19573784E-02,
     +  0.11362380E-02,-0.88227965E-03, 0.16668793E-01, 0.24813847E-03,
     +  0.22724843E-03,-0.50082086E-02,-0.70450806E-04,-0.14866653E-02,
     +  0.34253139E-03,-0.78132330E-03, 0.95232004E-04, 0.26346687E-02,
     + -0.11885086E-02, 0.16762946E-02, 0.41148295E-02,-0.36757353E-02,
     + -0.26984522E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_den    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      l_sp_den    =l_sp_den    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)*x11        *x43    
      l_sp_den    =l_sp_den    
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)    *x22    *x42*x51
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)                *x52
     8  +coeff( 26)*x11        *x41    
      l_sp_den    =l_sp_den    
     9  +coeff( 27)            *x43    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)                *x53
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)    *x22        *x53
     8  +coeff( 35)    *x23    *x43    
      l_sp_den    =l_sp_den    
     9  +coeff( 36)                *x51
     1  +coeff( 37)*x12*x21            
     2  +coeff( 38)    *x21    *x43    
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x22        *x52
     6  +coeff( 42)*x11*x21    *x41*x52
     7  +coeff( 43)*x12*x22    *x41    
     8  +coeff( 44)*x11*x22    *x43    
      l_sp_den    =l_sp_den    
     9  +coeff( 45)*x11*x22    *x41*x52
     1  +coeff( 46)*x12*x22    *x41*x51
     2  +coeff( 47)    *x22    *x43*x52
     3  +coeff( 48)*x11*x22    *x42*x52
     4  +coeff( 49)*x12*x22    *x41*x52
c
      return
      end
      function x_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3661939E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.18286741E-02, 0.27598858E+00,-0.53301449E-02, 0.13341466E+00,
     + -0.11502085E-01, 0.20360515E-01,-0.38431406E-01, 0.40398363E-01,
     +  0.15320870E-01,-0.35385662E-02, 0.13067214E-01,-0.21232285E-01,
     +  0.25515195E-01,-0.45640557E-02,-0.28484932E-02, 0.27492533E-02,
     + -0.59276195E-02,-0.10933421E-02,-0.20149248E-02,-0.11499448E-02,
     +  0.82819397E-02, 0.44278507E-02,-0.39077913E-02, 0.10380959E-01,
     +  0.96271746E-03, 0.29436310E-03, 0.16031506E-02, 0.17608901E-02,
     +  0.27361480E-02,-0.12861944E-02, 0.28891638E-01,-0.88876169E-02,
     + -0.53914811E-03,-0.97667403E-03,-0.43417615E-03,-0.21944773E-03,
     +  0.31388152E-03,-0.83279790E-03, 0.54292579E-03, 0.49194053E-03,
     + -0.46185418E-04,-0.66581141E-03, 0.55723656E-02,-0.30726814E-02,
     + -0.30404686E-02, 0.27817057E-02,-0.41866223E-02,-0.49455941E-03,
     +  0.36621098E-02, 0.15993206E-02,-0.11289186E-02,-0.29025106E-02,
     +  0.11327727E-02,-0.15111127E-02, 0.15348545E-02, 0.24207521E-02,
     + -0.11648608E-01,-0.36343506E-02, 0.32943310E-03, 0.52649793E-02,
     + -0.13373429E-01, 0.29366596E-02, 0.63284686E-04,-0.20753537E-03,
     + -0.12176429E-03, 0.40673962E-03,-0.15819807E-02,-0.17268940E-03,
     + -0.98028158E-05,-0.10833485E-03,-0.96048490E-03, 0.52895624E-03,
     + -0.31672738E-04, 0.48223492E-02, 0.12526581E-02,-0.24393153E-03,
     +  0.70641264E-02, 0.18620364E-02,-0.45732353E-03,-0.21747437E-02,
     +  0.10027175E-02, 0.49012606E-02, 0.93247980E-03,-0.23673302E-02,
     + -0.72557555E-03, 0.17169903E-02, 0.40088792E-03, 0.18327709E-02,
     +  0.50752779E-03,-0.31236033E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
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
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)            *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x24            
      x_sp_dex    =x_sp_dex    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)*x11*x22    *x41    
     5  +coeff( 23)    *x24    *x41    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)*x11            *x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 27)    *x22        *x51
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x24            
     7  +coeff( 34)            *x43    
     8  +coeff( 35)            *x42*x51
      x_sp_dex    =x_sp_dex    
     9  +coeff( 36)            *x41*x52
     1  +coeff( 37)*x11*x21        *x51
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)    *x22    *x41*x51
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)    *x21    *x42*x51
     6  +coeff( 42)*x11        *x43    
     7  +coeff( 43)    *x22    *x43    
     8  +coeff( 44)*x11*x23    *x41    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 45)    *x21    *x43*x51
     1  +coeff( 46)*x11*x22    *x42    
     2  +coeff( 47)    *x24    *x42    
     3  +coeff( 48)    *x24        *x52
     4  +coeff( 49)    *x23    *x42*x51
     5  +coeff( 50)*x11*x21    *x41*x52
     6  +coeff( 51)    *x23        *x53
     7  +coeff( 52)*x11*x24    *x41    
     8  +coeff( 53)*x11*x23    *x42    
      x_sp_dex    =x_sp_dex    
     9  +coeff( 54)*x11*x22    *x42*x51
     1  +coeff( 55)    *x24    *x42*x51
     2  +coeff( 56)*x11*x22    *x41*x52
     3  +coeff( 57)*x12*x24    *x41    
     4  +coeff( 58)*x12*x23    *x42    
     5  +coeff( 59)*x12*x22    *x43    
     6  +coeff( 60)*x11*x24    *x42*x51
     7  +coeff( 61)*x11*x22    *x43*x52
     8  +coeff( 62)*x12*x24    *x41*x52
      x_sp_dex    =x_sp_dex    
     9  +coeff( 63)*x12                
     1  +coeff( 64)*x11        *x41*x51
     2  +coeff( 65)    *x21    *x41*x52
     3  +coeff( 66)    *x21        *x53
     4  +coeff( 67)    *x24        *x51
     5  +coeff( 68)            *x41*x53
     6  +coeff( 69)    *x23    *x41*x51
     7  +coeff( 70)*x11        *x41*x52
     8  +coeff( 71)    *x22        *x53
      x_sp_dex    =x_sp_dex    
     9  +coeff( 72)*x12*x21    *x41    
     1  +coeff( 73)*x11*x23        *x51
     2  +coeff( 74)*x12*x22    *x41    
     3  +coeff( 75)    *x21    *x43*x52
     4  +coeff( 76)*x12        *x43    
     5  +coeff( 77)*x11*x22    *x43    
     6  +coeff( 78)    *x24        *x53
     7  +coeff( 79)*x12*x23        *x51
     8  +coeff( 80)    *x23    *x41*x53
      x_sp_dex    =x_sp_dex    
     9  +coeff( 81)*x11        *x43*x52
     1  +coeff( 82)*x11*x23    *x43    
     2  +coeff( 83)*x12*x21    *x42*x51
     3  +coeff( 84)*x11*x23    *x41*x52
     4  +coeff( 85)*x11*x23        *x53
     5  +coeff( 86)    *x21    *x43*x53
     6  +coeff( 87)*x12*x24        *x51
     7  +coeff( 88)*x11*x22    *x42*x52
     8  +coeff( 89)*x11*x22    *x41*x53
      x_sp_dex    =x_sp_dex    
     9  +coeff( 90)*x11*x21    *x43*x52
c
      return
      end
      function t_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/  0.5934969E+00/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.83046296E-04,-0.89400105E-01, 0.15304978E-02, 0.30802112E-01,
     +  0.25619741E-02, 0.11017520E-01,-0.76823733E-02,-0.43482999E-02,
     + -0.16176520E-02,-0.13723398E-02,-0.17352600E-02, 0.60400078E-02,
     + -0.73421714E-02,-0.58029062E-03, 0.82417147E-03, 0.64056477E-03,
     + -0.59665740E-03,-0.18056580E-02, 0.14252292E-06,-0.20506204E-03,
     +  0.30754163E-03,-0.76382066E-03,-0.47041560E-03,-0.66760386E-03,
     + -0.31948891E-02,-0.13028310E-02, 0.61336020E-03, 0.24787404E-03,
     + -0.56138629E-03, 0.53505629E-03, 0.35784108E-03,-0.14417594E-02,
     + -0.29617406E-02,-0.80622071E-02, 0.34862183E-03,-0.52998128E-03,
     +  0.22116120E-02,-0.22249608E-03, 0.62917330E-03, 0.79544890E-03,
     +  0.76900976E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      t_sp_dex    =t_sp_dex    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x21        *x52
     8  +coeff( 17)*x11*x22            
      t_sp_dex    =t_sp_dex    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)*x11        *x43    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)*x11*x21    *x42    
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)*x11*x22    *x42    
      t_sp_dex    =t_sp_dex    
     9  +coeff( 27)*x11        *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)    *x22    *x41*x51
     4  +coeff( 31)    *x22        *x52
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)    *x22    *x43    
     7  +coeff( 34)    *x23    *x43    
     8  +coeff( 35)            *x43    
      t_sp_dex    =t_sp_dex    
     9  +coeff( 36)*x11*x21    *x41    
     1  +coeff( 37)    *x21    *x43    
     2  +coeff( 38)            *x42*x52
     3  +coeff( 39)*x11*x23    *x41    
     4  +coeff( 40)*x11*x22    *x41*x52
     5  +coeff( 41)*x11*x21    *x43*x52
c
      return
      end
      function y_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 38)
      data ncoeff/ 37/
      data avdat/ -0.1641307E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24390474E-02, 0.55390947E-01, 0.25746845E-01,-0.22952133E-02,
     +  0.70375181E-02, 0.40280823E-01, 0.63739880E-02,-0.38221683E-01,
     + -0.82015703E-02,-0.18274758E-02,-0.74991869E-03,-0.13959253E-01,
     +  0.10918544E-02,-0.93960465E-03,-0.38660013E-02,-0.16299116E-02,
     + -0.24665991E-01,-0.77589871E-02, 0.10498346E-01,-0.12719003E-01,
     +  0.62851788E-04,-0.48522637E-02,-0.30422306E-02, 0.64256520E-03,
     + -0.11960989E-02,-0.43000592E-03,-0.11289956E-02,-0.22988362E-03,
     +  0.10120182E-02, 0.23827227E-02,-0.20596914E-02,-0.33541143E-03,
     + -0.72637121E-02,-0.10293805E-02,-0.33420851E-02, 0.24291317E-02,
     +  0.56503001E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_dex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_dex    =y_sp_dex    
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x43    
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)    *x22            
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x22    *x41    
      y_sp_dex    =y_sp_dex    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)            *x43*x52
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)            *x44    
     8  +coeff( 26)*x11                
      y_sp_dex    =y_sp_dex    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11*x21        *x51
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)*x11*x21    *x42    
     5  +coeff( 32)            *x45*x51
     6  +coeff( 33)    *x22    *x43    
     7  +coeff( 34)*x11*x21    *x41*x51
     8  +coeff( 35)    *x22    *x42*x51
      y_sp_dex    =y_sp_dex    
     9  +coeff( 36)*x11*x23    *x41    
     1  +coeff( 37)    *x23    *x44    
c
      return
      end
      function p_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.8351380E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.98304248E-04, 0.23286650E-03,-0.21766075E-02, 0.28718282E-02,
     + -0.10936948E-02,-0.85634757E-02, 0.87465759E-03,-0.20824361E-02,
     +  0.87524690E-02, 0.16306038E-02,-0.19628665E-03, 0.29222906E-03,
     +  0.18417140E-02,-0.49886787E-02, 0.55735966E-03,-0.39557891E-03,
     + -0.21270043E-02,-0.18646098E-02,-0.10404852E-02,-0.30924124E-03,
     + -0.39529556E-03,-0.12387605E-02,-0.23302258E-04,-0.23872851E-03,
     +  0.88449731E-03, 0.48305839E-03,-0.90991281E-03,-0.67707129E-04,
     +  0.12539592E-04, 0.16709980E-03,-0.27184535E-04, 0.27059208E-03,
     +  0.10101016E-03,-0.25523725E-03,-0.12729679E-03,-0.99742634E-03,
     +  0.46922316E-03,-0.36608217E-04, 0.51892614E-04, 0.15315181E-03,
     + -0.10973218E-03,-0.13759447E-03,-0.73872099E-04, 0.70054404E-04,
     +  0.18831000E-03, 0.21636339E-03, 0.21031447E-03,-0.16590949E-03,
     + -0.12986833E-03,-0.45427226E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_dex    =p_sp_dex    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x22        *x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x22    *x41*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)                *x53
     7  +coeff( 25)    *x23    *x41    
     8  +coeff( 26)    *x23        *x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 27)    *x22    *x42*x51
     1  +coeff( 28)*x11                
     2  +coeff( 29)*x11            *x51
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)            *x42*x52
     7  +coeff( 34)*x11*x21    *x42    
     8  +coeff( 35)*x11*x21    *x41*x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)*x11*x21    *x41    
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)            *x43*x51
     5  +coeff( 41)    *x21    *x41*x52
     6  +coeff( 42)*x11*x22    *x41    
     7  +coeff( 43)*x11        *x43    
     8  +coeff( 44)*x11*x22        *x51
      p_sp_dex    =p_sp_dex    
     9  +coeff( 45)    *x22    *x41*x52
     1  +coeff( 46)    *x22        *x53
     2  +coeff( 47)*x11*x23    *x41    
     3  +coeff( 48)*x11*x21    *x42*x51
     4  +coeff( 49)*x12*x21    *x42    
     5  +coeff( 50)    *x22    *x43*x51
c
      return
      end
      function l_sp_dex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 44)
      data ncoeff/ 43/
      data avdat/ -0.1068454E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.69820900E-02,-0.40051723E+00, 0.32255938E-02,-0.98872170E-01,
     +  0.14727497E-01,-0.34658611E-01, 0.54101016E-01,-0.41393116E-01,
     + -0.21849278E-01,-0.75505236E-02, 0.28750094E-01,-0.35715811E-01,
     +  0.17665175E-02,-0.30824104E-02,-0.95893843E-02,-0.10478101E-03,
     + -0.15902253E-01,-0.90350577E-03,-0.88557479E-03,-0.20894941E-02,
     +  0.32933219E-02, 0.20457201E-02,-0.12749022E-01,-0.56452597E-02,
     +  0.29289960E-02, 0.18731805E-02,-0.21065250E-02, 0.11163820E-02,
     + -0.26220921E-02,-0.65775374E-02,-0.38102105E-01, 0.47602956E-03,
     +  0.73184422E-03, 0.17351280E-03,-0.40987600E-03, 0.11230943E-01,
     + -0.13870271E-02, 0.55141235E-03,-0.16596274E-02, 0.22613206E-02,
     +  0.33151587E-02,-0.23238345E-02, 0.27013232E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)                *x52
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)*x11        *x43    
     8  +coeff( 17)    *x23    *x42    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x21    *x41*x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)            *x43    
      l_sp_dex    =l_sp_dex    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)            *x42    
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)*x12*x21            
      l_sp_dex    =l_sp_dex    
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x21    *x42*x51
     2  +coeff( 38)*x11*x23            
     3  +coeff( 39)*x11*x21    *x42    
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x22    *x41*x52
     6  +coeff( 42)*x12*x22    *x41*x51
     7  +coeff( 43)*x12*x22    *x41*x52
c
      return
      end
      function x_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.2241619E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.77176250E-04, 0.17129853E+00,-0.33838567E-02, 0.13863553E+00,
     + -0.80186594E-02, 0.21411659E-01,-0.25195856E-01, 0.31680658E-01,
     +  0.11036162E-01,-0.60518007E-02,-0.14448079E-01, 0.17132448E-01,
     + -0.36344258E-02, 0.84247533E-02, 0.38140228E-02,-0.27352751E-02,
     +  0.19096520E-02,-0.34575462E-02,-0.10607891E-02,-0.14896691E-02,
     +  0.37582505E-02, 0.63586826E-02, 0.22981255E-02, 0.13646749E-02,
     +  0.19039235E-02, 0.51847496E-02,-0.76581887E-03,-0.33435426E-02,
     +  0.18116102E-01, 0.17248567E-03, 0.76449265E-04,-0.45711838E-03,
     +  0.41468797E-03,-0.54605869E-02,-0.58976305E-03,-0.67852897E-03,
     + -0.48813535E-03,-0.88053297E-04,-0.50366652E-03,-0.11006631E-03,
     +  0.22416802E-03, 0.28650547E-03,-0.23948750E-03, 0.42333610E-04,
     +  0.10641884E-02,-0.19678327E-02,-0.11048076E-02, 0.22728259E-02,
     + -0.27351764E-02, 0.10932691E-03, 0.30906955E-02, 0.13399225E-03,
     + -0.24769802E-02, 0.12526014E-02, 0.10880196E-02, 0.33122355E-02,
     + -0.83546073E-03,-0.94174634E-03,-0.11234049E-01,-0.14166973E-01,
     + -0.19769827E-02,-0.17733108E-03,-0.20258807E-03,-0.22802327E-02,
     +  0.93395654E-04,-0.18772018E-02,-0.18140886E-03, 0.65455053E-04,
     + -0.11678237E-03, 0.23512109E-02, 0.78477297E-03, 0.19162660E-02,
     + -0.45794793E-03,-0.26074825E-02, 0.64547517E-03, 0.71012510E-04,
     +  0.14356417E-02,-0.40632654E-02, 0.33709427E-03, 0.49921230E-03,
     + -0.48979017E-03,-0.23666107E-03, 0.25551254E-02, 0.37356333E-02,
     +  0.70675270E-03,-0.17360987E-02,-0.16747863E-02, 0.78677554E-02,
     +  0.17566262E-02, 0.55934652E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)*x11*x22            
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 18)    *x24            
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11        *x41    
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x21    *x41    
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x22    *x42    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 27)*x11*x23            
     1  +coeff( 28)    *x24    *x41    
     2  +coeff( 29)    *x23    *x43    
     3  +coeff( 30)*x11*x21            
     4  +coeff( 31)*x11            *x51
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)                *x53
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)*x11*x24            
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 36)            *x43    
     1  +coeff( 37)            *x42*x51
     2  +coeff( 38)            *x41*x52
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)*x11        *x41*x51
     5  +coeff( 41)    *x22    *x41*x51
     6  +coeff( 42)*x12*x21            
     7  +coeff( 43)    *x21    *x42*x51
     8  +coeff( 44)    *x23    *x41*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 45)    *x22    *x42*x51
     1  +coeff( 46)*x11*x23    *x41    
     2  +coeff( 47)    *x21    *x43*x51
     3  +coeff( 48)*x11*x22    *x42    
     4  +coeff( 49)    *x24    *x42    
     5  +coeff( 50)*x11*x22    *x41*x51
     6  +coeff( 51)    *x23    *x42*x51
     7  +coeff( 52)    *x23    *x41*x52
     8  +coeff( 53)*x11*x24    *x41    
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 54)*x11*x23    *x42    
     1  +coeff( 55)*x11*x22    *x41*x52
     2  +coeff( 56)    *x24        *x53
     3  +coeff( 57)*x11*x23    *x41*x52
     4  +coeff( 58)*x12*x21    *x42*x52
     5  +coeff( 59)*x11*x24    *x43*x52
     6  +coeff( 60)*x12*x24    *x43*x52
     7  +coeff( 61)    *x24        *x51
     8  +coeff( 62)            *x41*x53
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 63)*x11        *x42*x51
     1  +coeff( 64)    *x22    *x41*x52
     2  +coeff( 65)*x11            *x53
     3  +coeff( 66)    *x22        *x53
     4  +coeff( 67)*x11*x23    *x41*x51
     5  +coeff( 68)*x12*x24            
     6  +coeff( 69)*x12        *x43    
     7  +coeff( 70)*x11*x22    *x43    
     8  +coeff( 71)*x11*x22    *x42*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 72)    *x24    *x41*x52
     1  +coeff( 73)*x12*x23        *x51
     2  +coeff( 74)    *x23    *x43*x51
     3  +coeff( 75)*x11*x21    *x41*x53
     4  +coeff( 76)    *x22    *x43*x52
     5  +coeff( 77)*x11*x23    *x43    
     6  +coeff( 78)*x12*x24    *x41    
     7  +coeff( 79)*x12*x24        *x51
     8  +coeff( 80)    *x24    *x41*x53
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 81)*x12*x23    *x42    
     1  +coeff( 82)*x11*x21    *x43*x52
     2  +coeff( 83)    *x23    *x43*x52
     3  +coeff( 84)*x12*x22    *x43    
     4  +coeff( 85)*x12*x22    *x41*x52
     5  +coeff( 86)*x11*x23    *x42*x52
     6  +coeff( 87)*x11*x23    *x41*x53
     7  +coeff( 88)    *x24    *x43*x52
     8  +coeff( 89)*x12*x23    *x42*x51
      x_sp_q3en   =x_sp_q3en   
     9  +coeff( 90)*x12*x24    *x41*x52
c
      return
      end
      function t_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 43)
      data ncoeff/ 42/
      data avdat/ -0.1090703E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.58552553E-03,-0.66510260E-01, 0.11316757E-02, 0.22735398E-01,
     +  0.19440851E-02,-0.48979660E-02, 0.82121557E-02,-0.48080888E-02,
     + -0.33330738E-02,-0.14109076E-02,-0.12713238E-02,-0.52802456E-02,
     +  0.37979316E-02,-0.56126103E-03, 0.50930952E-03, 0.45018535E-03,
     + -0.15875178E-03,-0.40477407E-03, 0.56543091E-03,-0.10412288E-02,
     + -0.87240350E-03,-0.21874006E-02, 0.15771613E-03, 0.24777433E-03,
     + -0.19512497E-03,-0.20776249E-02,-0.56201443E-02, 0.32144101E-03,
     +  0.11105150E-03, 0.10874112E-03,-0.38770400E-03,-0.62405015E-04,
     +  0.16790463E-02,-0.49245090E-03,-0.42155071E-03, 0.10457728E-03,
     +  0.40131519E-03,-0.31815475E-03, 0.43492537E-03, 0.15229256E-03,
     + -0.81212237E-03, 0.45369141E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)            *x42    
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)*x11            *x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)            *x42*x51
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)*x11*x21    *x42    
     8  +coeff( 26)    *x22    *x43    
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 27)    *x23    *x43    
     1  +coeff( 28)            *x43    
     2  +coeff( 29)            *x41*x52
     3  +coeff( 30)                *x53
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)*x12*x21            
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)    *x23        *x51
     8  +coeff( 35)    *x21    *x42*x51
      t_sp_q3en   =t_sp_q3en   
     9  +coeff( 36)            *x43*x51
     1  +coeff( 37)    *x23    *x41*x51
     2  +coeff( 38)    *x23        *x52
     3  +coeff( 39)*x11*x23    *x41    
     4  +coeff( 40)*x11        *x42*x52
     5  +coeff( 41)*x11*x22    *x43    
     6  +coeff( 42)*x11*x22    *x41*x52
c
      return
      end
      function y_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 39)
      data ncoeff/ 38/
      data avdat/ -0.1132265E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22210993E-02, 0.55784378E-01, 0.29288331E-01,-0.17964842E-02,
     +  0.79352437E-02, 0.49855839E-01, 0.77981488E-02,-0.47073279E-01,
     + -0.10513500E-01,-0.24170768E-02,-0.15721044E-01, 0.13436724E-02,
     + -0.41417675E-02,-0.14084657E-02,-0.57913857E-02,-0.18980344E-02,
     + -0.29430533E-01,-0.99547794E-02, 0.12049542E-01,-0.12076062E-02,
     + -0.14052457E-01,-0.56522503E-02, 0.15942573E-02,-0.13520994E-02,
     + -0.47334481E-03,-0.78545493E-03, 0.58818306E-03,-0.13882595E-02,
     + -0.23850155E-03, 0.12240253E-02, 0.34645444E-02,-0.23154889E-02,
     + -0.77363797E-02,-0.10795661E-02,-0.45442530E-02, 0.32626418E-02,
     +  0.75231860E-02, 0.25838881E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_q3en   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_q3en   =y_sp_q3en   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x43    
     2  +coeff( 11)    *x22            
     3  +coeff( 12)*x11        *x41    
     4  +coeff( 13)            *x41*x52
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)*x11*x21            
     8  +coeff( 17)    *x22    *x41    
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)    *x23            
     2  +coeff( 20)            *x44*x51
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)    *x23    *x42    
     6  +coeff( 24)            *x44    
     7  +coeff( 25)*x11                
     8  +coeff( 26)    *x21        *x52
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 27)            *x42*x52
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11*x21        *x51
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)*x11*x21    *x42    
     6  +coeff( 33)    *x22    *x43    
     7  +coeff( 34)*x11*x21    *x41*x51
     8  +coeff( 35)    *x22    *x42*x51
      y_sp_q3en   =y_sp_q3en   
     9  +coeff( 36)*x11*x23    *x41    
     1  +coeff( 37)    *x23    *x44    
     2  +coeff( 38)*x11*x22    *x42*x51
c
      return
      end
      function p_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.6651885E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.52559688E-04, 0.86911675E-03, 0.37569855E-02,-0.15671937E-02,
     + -0.92360424E-02, 0.10397243E-02,-0.20316101E-02, 0.99436194E-02,
     +  0.18405362E-02,-0.25018596E-03, 0.28205037E-03, 0.19915884E-02,
     + -0.59794812E-02, 0.59009070E-03,-0.47950778E-03,-0.24373804E-02,
     + -0.14376020E-02,-0.13135556E-02,-0.40298319E-03,-0.15524648E-02,
     + -0.10828734E-03,-0.22239565E-03,-0.24762092E-03,-0.11092883E-04,
     +  0.34734825E-03,-0.12105353E-02,-0.66910154E-03, 0.12368022E-03,
     + -0.60979837E-04, 0.14384354E-03, 0.69424015E-03, 0.26165240E-03,
     + -0.36046366E-03, 0.86056200E-04,-0.15603349E-03, 0.19942256E-03,
     + -0.64437278E-03, 0.59030725E-04,-0.83654304E-04, 0.73078743E-04,
     + -0.35696561E-04,-0.62205618E-04,-0.55636065E-04, 0.32526924E-03,
     +  0.24042596E-03, 0.17398442E-03,-0.22010850E-03, 0.18670452E-03,
     + -0.11648704E-03,-0.26439782E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x22            
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff(  9)                *x52
     1  +coeff( 10)*x11*x21            
     2  +coeff( 11)*x11        *x41    
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)    *x21    *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 18)    *x22    *x41*x51
     1  +coeff( 19)            *x41*x52
     2  +coeff( 20)    *x22    *x42    
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)                *x53
     6  +coeff( 24)*x11*x21        *x51
     7  +coeff( 25)    *x23        *x51
     8  +coeff( 26)    *x22    *x43    
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 27)    *x22    *x42*x51
     1  +coeff( 28)    *x21            
     2  +coeff( 29)*x11                
     3  +coeff( 30)*x11*x22            
     4  +coeff( 31)    *x23    *x41    
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)*x11*x22        *x51
     8  +coeff( 35)*x11*x21    *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 36)*x11*x23    *x41    
     1  +coeff( 37)*x12*x22    *x42*x51
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)*x11            *x52
     4  +coeff( 40)            *x42*x52
     5  +coeff( 41)*x11        *x43    
     6  +coeff( 42)*x12*x22            
     7  +coeff( 43)*x12*x21    *x41    
     8  +coeff( 44)    *x23    *x41*x51
      p_sp_q3en   =p_sp_q3en   
     9  +coeff( 45)    *x22    *x41*x52
     1  +coeff( 46)    *x22        *x53
     2  +coeff( 47)*x11*x21    *x42*x51
     3  +coeff( 48)*x11*x22        *x52
     4  +coeff( 49)*x12*x21    *x42    
     5  +coeff( 50)*x11*x22    *x43    
c
      return
      end
      function l_sp_q3en   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.6745730E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.46294918E-02,-0.26224601E+00,-0.32164164E-01, 0.89591555E-02,
     + -0.29312681E-01, 0.34823027E-01,-0.19452387E-01,-0.14373324E-01,
     + -0.19650739E-02, 0.17944790E-01,-0.22749176E-01,-0.15533245E-02,
     + -0.37074981E-02, 0.97046158E-03,-0.21436177E-02,-0.20337049E-02,
     + -0.49963370E-02, 0.44920153E-03,-0.10156746E-01, 0.98348816E-03,
     + -0.74151234E-03, 0.12952310E-02,-0.24613726E-02,-0.59158570E-04,
     +  0.15908332E-02, 0.13121965E-02, 0.50396839E-03,-0.13100025E-02,
     +  0.71347208E-03, 0.95371774E-03,-0.29744238E-02,-0.69293529E-02,
     + -0.25666339E-01,-0.29184078E-03,-0.26607074E-03, 0.75962855E-02,
     + -0.16322270E-02,-0.10206326E-02,-0.35614222E-02, 0.22092920E-02,
     +  0.26229087E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     8  +coeff(  8)    *x23            
      l_sp_q3en   =l_sp_q3en   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)*x11*x22            
     8  +coeff( 17)    *x22    *x42    
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 18)*x11        *x43    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)            *x41    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)*x11*x22    *x42    
     6  +coeff( 24)*x11*x21            
     7  +coeff( 25)*x11        *x41    
     8  +coeff( 26)    *x21    *x41*x51
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 27)            *x42*x51
     1  +coeff( 28)*x11*x21    *x41    
     2  +coeff( 29)*x11        *x42    
     3  +coeff( 30)*x11*x23            
     4  +coeff( 31)*x11*x22    *x41    
     5  +coeff( 32)    *x22    *x43    
     6  +coeff( 33)    *x23    *x43    
     7  +coeff( 34)                *x52
     8  +coeff( 35)*x12*x21            
      l_sp_q3en   =l_sp_q3en   
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x11*x22    *x43    
     4  +coeff( 40)*x11*x22    *x41*x52
     5  +coeff( 41)*x11*x23    *x41*x52
c
      return
      end
      function x_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.3048534E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12270153E-02, 0.86975597E-01,-0.18126414E-02, 0.19257694E+00,
     + -0.59663206E-02, 0.15642000E-01,-0.15411070E-01, 0.27440697E-01,
     + -0.99735837E-02, 0.69274884E-02,-0.32949233E-02, 0.51163388E-02,
     +  0.33560295E-02,-0.10358118E-01,-0.32468573E-02,-0.28873954E-02,
     +  0.11081462E-01,-0.66389548E-04,-0.11279774E-02, 0.87412895E-03,
     +  0.45693141E-05, 0.11173950E-01, 0.74860518E-03, 0.11811148E-02,
     +  0.38356765E-02,-0.32407180E-02,-0.21545708E-02, 0.19997947E-02,
     + -0.95778622E-03,-0.22602521E-03, 0.85984543E-03, 0.46706237E-03,
     +  0.34798728E-02,-0.20680975E-02, 0.42172619E-02,-0.12093229E-02,
     + -0.19434292E-03,-0.30459528E-03,-0.56026457E-03, 0.66212303E-03,
     +  0.48223184E-03, 0.75806968E-03,-0.12340718E-02,-0.28388731E-02,
     + -0.87936700E-03,-0.48376195E-03,-0.27298348E-03, 0.82044375E-04,
     +  0.16356101E-03,-0.23539433E-03, 0.44843880E-03,-0.44931224E-03,
     +  0.19891378E-02, 0.14138785E-02,-0.88347500E-03,-0.18041342E-02,
     +  0.17968570E-02,-0.69575111E-03,-0.10565884E-02, 0.11951183E-03,
     + -0.63685940E-04, 0.26837341E-02, 0.22552130E-02,-0.62633277E-03,
     + -0.22079829E-03, 0.30224290E-03,-0.19731375E-02,-0.42629879E-03,
     + -0.95675880E-03, 0.11147660E-02, 0.11765623E-02, 0.93282078E-03,
     +  0.93527266E-03, 0.14245776E-02,-0.27522072E-03, 0.24222925E-03,
     + -0.38655165E-02,-0.53609598E-04, 0.22378697E-04, 0.11693266E-03,
     + -0.59559883E-04,-0.82606194E-03,-0.15585637E-03,-0.14675580E-03,
     +  0.23651193E-03, 0.10908224E-03, 0.48564161E-04,-0.13678727E-02,
     +  0.17003968E-03,-0.13659128E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      x_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x23            
     2  +coeff( 11)            *x42    
     3  +coeff( 12)    *x22    *x41    
     4  +coeff( 13)    *x22        *x51
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x21    *x41*x51
     7  +coeff( 16)    *x24            
     8  +coeff( 17)    *x23    *x41    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 18)            *x43*x51
     1  +coeff( 19)            *x41*x51
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)*x11        *x43    
     4  +coeff( 22)    *x23    *x43    
     5  +coeff( 23)                *x53
     6  +coeff( 24)    *x23        *x51
     7  +coeff( 25)    *x22    *x42    
     8  +coeff( 26)    *x21    *x43    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 27)    *x24    *x41    
     1  +coeff( 28)*x11*x22    *x42    
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)    *x21        *x52
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)*x11*x22    *x41    
     7  +coeff( 34)    *x24        *x51
     8  +coeff( 35)    *x23    *x42    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 36)*x11*x23    *x41    
     1  +coeff( 37)*x11            *x51
     2  +coeff( 38)*x11        *x42    
     3  +coeff( 39)    *x22        *x52
     4  +coeff( 40)*x11*x22        *x51
     5  +coeff( 41)*x11*x21    *x42    
     6  +coeff( 42)    *x23    *x41*x51
     7  +coeff( 43)    *x21    *x43*x51
     8  +coeff( 44)*x11*x24    *x41    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 45)*x11*x23        *x52
     1  +coeff( 46)            *x43    
     2  +coeff( 47)            *x42*x51
     3  +coeff( 48)*x11            *x52
     4  +coeff( 49)*x12*x21            
     5  +coeff( 50)    *x21        *x53
     6  +coeff( 51)*x11*x21        *x52
     7  +coeff( 52)    *x23        *x52
     8  +coeff( 53)    *x22    *x43    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 54)    *x22    *x42*x51
     1  +coeff( 55)    *x22    *x41*x52
     2  +coeff( 56)    *x24    *x42    
     3  +coeff( 57)    *x23    *x42*x51
     4  +coeff( 58)*x11*x24        *x51
     5  +coeff( 59)*x11*x23    *x41*x51
     6  +coeff( 60)    *x21    *x42*x53
     7  +coeff( 61)*x11*x22    *x41*x52
     8  +coeff( 62)    *x24    *x41*x52
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 63)    *x24        *x53
     1  +coeff( 64)*x11*x21    *x42*x52
     2  +coeff( 65)    *x22    *x42*x53
     3  +coeff( 66)    *x21    *x43*x53
     4  +coeff( 67)*x12*x24    *x41    
     5  +coeff( 68)*x11*x22    *x42*x52
     6  +coeff( 69)*x12*x23    *x42    
     7  +coeff( 70)    *x23    *x43*x52
     8  +coeff( 71)*x12*x22    *x43    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 72)*x12*x22    *x41*x52
     1  +coeff( 73)*x12*x24    *x41*x51
     2  +coeff( 74)*x12*x23    *x42*x51
     3  +coeff( 75)*x12*x23        *x53
     4  +coeff( 76)*x11*x21    *x43*x53
     5  +coeff( 77)*x11*x24    *x43*x52
     6  +coeff( 78)*x11*x21            
     7  +coeff( 79)*x12                
     8  +coeff( 80)*x11*x21        *x51
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 81)*x11        *x41*x51
     1  +coeff( 82)    *x21    *x42*x51
     2  +coeff( 83)            *x42*x52
     3  +coeff( 84)            *x41*x53
     4  +coeff( 85)*x11*x21    *x41*x51
     5  +coeff( 86)*x11        *x41*x52
     6  +coeff( 87)*x11            *x53
     7  +coeff( 88)    *x22        *x53
     8  +coeff( 89)*x12*x21    *x41    
      x_sp_q3m    =x_sp_q3m    
     9  +coeff( 90)*x11*x23        *x51
c
      return
      end
      function t_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.1536579E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.15644664E-02,-0.37594311E-01, 0.66807045E-03, 0.67421712E-01,
     + -0.58168294E-02, 0.37434231E-02,-0.10310899E-02, 0.10503181E-02,
     + -0.12773447E-02, 0.27320752E-03,-0.38231257E-03, 0.68091200E-03,
     + -0.24703548E-02,-0.22289461E-03,-0.64064102E-03, 0.96619484E-03,
     + -0.56987064E-03,-0.77688746E-04, 0.21590156E-03,-0.22739499E-03,
     + -0.26633751E-03, 0.18489490E-03, 0.22948418E-03,-0.16987039E-02,
     +  0.12992181E-03,-0.25190890E-03,-0.20960651E-03,-0.50384505E-03,
     + -0.46196819E-03,-0.21880081E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)                *x52
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x21        *x51
      t_sp_q3m    =t_sp_q3m    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)*x11                
     2  +coeff( 11)            *x42    
     3  +coeff( 12)                *x53
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)    *x21    *x41*x51
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 18)*x11            *x51
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)    *x21        *x52
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11*x23            
     5  +coeff( 23)*x11        *x43    
     6  +coeff( 24)    *x23    *x42    
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)    *x21    *x42*x51
      t_sp_q3m    =t_sp_q3m    
     9  +coeff( 27)    *x22        *x52
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)*x11*x22    *x42    
     3  +coeff( 30)*x11*x21    *x43    
c
      return
      end
      function y_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.1697584E-04/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.21410794E-02, 0.50222524E-01, 0.30922478E-01,-0.15491729E-02,
     +  0.88455575E-02, 0.58031656E-01,-0.52621635E-03, 0.97232079E-02,
     + -0.54106336E-01,-0.12423373E-01,-0.29119451E-02,-0.89432637E-03,
     + -0.15985914E-01, 0.15245290E-02,-0.13868940E-02,-0.72115506E-02,
     + -0.19698939E-02,-0.34159537E-01,-0.12136832E-01, 0.13561646E-01,
     + -0.14734507E-01, 0.72768435E-03,-0.71854698E-02,-0.71241177E-03,
     + -0.40222001E-02, 0.15609673E-02,-0.14930953E-02,-0.10045166E-02,
     +  0.13292016E-02, 0.40505058E-02,-0.30376256E-03, 0.84044681E-04,
     + -0.31619036E-03,-0.22441226E-02, 0.13432353E-02,-0.82291458E-02,
     + -0.11715306E-02,-0.50718221E-02, 0.83701797E-02, 0.26297432E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)*x11                
     8  +coeff(  8)                *x52
      y_sp_q3m    =y_sp_q3m    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)            *x43    
     3  +coeff( 12)            *x42*x51
     4  +coeff( 13)    *x22            
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)*x11*x21            
      y_sp_q3m    =y_sp_q3m    
     9  +coeff( 18)    *x22    *x41    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)    *x23            
     3  +coeff( 21)    *x22    *x42    
     4  +coeff( 22)            *x43*x52
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x21    *x44    
     7  +coeff( 25)            *x41*x52
     8  +coeff( 26)    *x21    *x42    
      y_sp_q3m    =y_sp_q3m    
     9  +coeff( 27)            *x44    
     1  +coeff( 28)    *x21        *x52
     2  +coeff( 29)*x11*x22            
     3  +coeff( 30)    *x23    *x41    
     4  +coeff( 31)            *x43*x51
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)*x11*x21        *x51
     7  +coeff( 34)*x11*x21    *x42    
     8  +coeff( 35)    *x23        *x51
      y_sp_q3m    =y_sp_q3m    
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)*x11*x21    *x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x23    *x44    
     4  +coeff( 40)*x11*x22    *x42*x51
c
      return
      end
      function p_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.8995809E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.37753026E-03, 0.62870217E-03,-0.13525983E-01,-0.47308281E-02,
     +  0.28819013E-02, 0.52643539E-02,-0.99184341E-03, 0.12218861E-02,
     + -0.45552985E-02,-0.17244612E-02, 0.29005543E-02, 0.14704511E-02,
     +  0.21674572E-02,-0.22647969E-03, 0.29762718E-03, 0.54110761E-03,
     + -0.20640412E-03, 0.21980562E-03, 0.25037571E-03, 0.29406298E-03,
     + -0.32986628E-03, 0.16962297E-02, 0.69103502E-04,-0.17138661E-03,
     + -0.18555226E-03, 0.33514283E-03, 0.30094496E-03,-0.29142320E-04,
     + -0.13631663E-02, 0.21664693E-02,-0.32310149E-04,-0.18312714E-03,
     + -0.20299065E-03,-0.28092597E-03,-0.12612289E-03, 0.40576578E-03,
     + -0.55868452E-03, 0.22334716E-03, 0.70740766E-03, 0.96604822E-03,
     +  0.13812751E-04,-0.96952470E-04,-0.23659572E-04, 0.34997007E-04,
     +  0.34077126E-04, 0.23234606E-03,-0.16598499E-03,-0.12171536E-03,
     + -0.11607955E-03, 0.10024323E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_q3m    =p_sp_q3m    
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)            *x41*x52
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)                *x52
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11        *x41    
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)            *x42*x51
     2  +coeff( 20)                *x53
     3  +coeff( 21)    *x23    *x41    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)*x11                
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)*x11*x22            
     8  +coeff( 26)*x11*x21    *x41    
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)    *x23    *x43    
     3  +coeff( 30)*x11*x23    *x43*x53
     4  +coeff( 31)*x11            *x51
     5  +coeff( 32)    *x21    *x41*x51
     6  +coeff( 33)    *x22        *x52
     7  +coeff( 34)    *x21    *x41*x52
     8  +coeff( 35)            *x41*x53
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)*x12*x21    *x42    
     3  +coeff( 39)*x12*x22    *x42*x51
     4  +coeff( 40)    *x23    *x43*x52
     5  +coeff( 41)*x12                
     6  +coeff( 42)    *x21    *x42    
     7  +coeff( 43)*x11        *x42    
     8  +coeff( 44)*x11        *x41*x51
      p_sp_q3m    =p_sp_q3m    
     9  +coeff( 45)*x11            *x52
     1  +coeff( 46)    *x21    *x43    
     2  +coeff( 47)    *x22    *x41*x51
     3  +coeff( 48)            *x43*x51
     4  +coeff( 49)            *x42*x52
     5  +coeff( 50)*x11        *x43    
c
      return
      end
      function l_sp_q3m    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.7771165E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56771142E-02,-0.26221728E+00,-0.32086298E-01, 0.89663928E-02,
     + -0.31968836E-01, 0.34885924E-01,-0.16627824E-01,-0.14725816E-01,
     + -0.20120936E-02, 0.17944967E-01,-0.22747541E-01,-0.16212424E-02,
     + -0.14508030E-02,-0.34589607E-02, 0.96190028E-03,-0.22290177E-02,
     + -0.20191823E-02,-0.48606941E-02, 0.43729346E-03,-0.10146136E-01,
     +  0.10084007E-02,-0.80208009E-03, 0.12199904E-02,-0.24478503E-02,
     +  0.42402546E-04, 0.15924966E-02,-0.12916395E-02, 0.69833372E-03,
     +  0.96690207E-03,-0.29571562E-02,-0.66571212E-02,-0.25974169E-01,
     +  0.10746599E-02, 0.52687008E-03,-0.25876763E-03, 0.77153994E-02,
     + -0.14920792E-02,-0.10108766E-02,-0.36355646E-02, 0.21347404E-02,
     +  0.24853407E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_q3m    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      l_sp_q3m    =l_sp_q3m    
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x23    *x41    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)                *x52
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)            *x43    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x22            
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)*x11        *x43    
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)            *x41    
     4  +coeff( 22)*x11            *x51
     5  +coeff( 23)    *x21        *x52
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)*x11        *x41    
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)*x11        *x42    
     2  +coeff( 29)*x11*x23            
     3  +coeff( 30)*x11*x22    *x41    
     4  +coeff( 31)    *x22    *x43    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)    *x21    *x41*x51
     7  +coeff( 34)            *x42*x51
     8  +coeff( 35)*x12*x21            
      l_sp_q3m    =l_sp_q3m    
     9  +coeff( 36)    *x21    *x43    
     1  +coeff( 37)    *x23        *x51
     2  +coeff( 38)    *x21    *x42*x51
     3  +coeff( 39)*x11*x22    *x43    
     4  +coeff( 40)*x11*x22    *x41*x52
     5  +coeff( 41)*x11*x23    *x41*x52
c
      return
      end
      function x_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1004723E-01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.13131803E-01,-0.14425130E-02, 0.32378420E+00,-0.72478461E-02,
     + -0.14702863E-01, 0.33028215E-01,-0.21014007E-01, 0.70569613E-02,
     + -0.12112223E-01,-0.36374775E-02, 0.80806596E-04, 0.40835589E-02,
     +  0.48901617E-04, 0.53960856E-01, 0.18668631E-01,-0.49187080E-02,
     +  0.53073317E-02,-0.50160550E-02, 0.11664884E-01,-0.20613745E-02,
     + -0.17139859E-02, 0.26776236E-02,-0.75362041E-03, 0.40120301E-02,
     + -0.67984866E-03, 0.95331913E-03, 0.11022643E-02, 0.47911606E-02,
     +  0.20397843E-02,-0.24483720E-03, 0.39283070E-03, 0.11531514E-02,
     + -0.24750249E-02,-0.26392790E-02,-0.11556695E-02, 0.95677674E-02,
     + -0.52545383E-03,-0.23429096E-02, 0.22768797E-02,-0.98157884E-03,
     + -0.94472831E-04, 0.16595483E-02, 0.17700443E-02,-0.35005314E-02,
     +  0.11339571E-03,-0.16442458E-02,-0.12500635E-02, 0.73069596E-03,
     +  0.22535570E-03, 0.24350693E-02,-0.21901648E-02,-0.21144433E-02,
     +  0.23060307E-03, 0.97118056E-03, 0.20462547E-02, 0.87668793E-03,
     + -0.12994049E-02, 0.54609217E-02, 0.24328998E-02,-0.11862969E-02,
     + -0.13450255E-02,-0.33237837E-04,-0.54982364E-04, 0.21985757E-03,
     +  0.21013341E-03,-0.20979762E-03,-0.12945784E-03, 0.32126985E-03,
     +  0.68196625E-03,-0.17651544E-03, 0.11739741E-03, 0.10620505E-02,
     +  0.94483834E-03, 0.38422417E-03,-0.12132803E-02,-0.11190464E-02,
     + -0.72599761E-03,-0.11707329E-02,-0.48607343E-03,-0.13803059E-02,
     + -0.10632703E-03, 0.15086857E-03, 0.12776052E-02, 0.28699794E-03,
     + -0.12623727E-02, 0.33760027E-03,-0.68603043E-03,-0.81066077E-03,
     + -0.10033435E-02, 0.37855585E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x23            
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x24            
     2  +coeff( 11)*x12*x21            
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)*x12*x23            
     5  +coeff( 14)    *x21            
     6  +coeff( 15)    *x22            
     7  +coeff( 16)            *x42    
     8  +coeff( 17)    *x22    *x41    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 18)    *x21    *x41*x51
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x22        *x53
     3  +coeff( 21)            *x41*x51
     4  +coeff( 22)                *x53
     5  +coeff( 23)*x11        *x41    
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)    *x21        *x52
     8  +coeff( 26)*x11*x22            
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)    *x22    *x42    
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)*x11*x21    *x41    
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)    *x24    *x41    
     7  +coeff( 34)    *x24        *x51
     8  +coeff( 35)    *x21    *x43*x51
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 36)    *x23    *x43    
     1  +coeff( 37)            *x43    
     2  +coeff( 38)    *x21    *x43    
     3  +coeff( 39)    *x23    *x41*x51
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)*x11*x24            
     6  +coeff( 42)    *x22    *x42*x51
     7  +coeff( 43)*x11*x22    *x43    
     8  +coeff( 44)*x11*x23    *x41*x52
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 45)*x11*x21            
     1  +coeff( 46)    *x22        *x52
     2  +coeff( 47)    *x21    *x42*x51
     3  +coeff( 48)*x11*x22        *x51
     4  +coeff( 49)*x11*x21    *x42    
     5  +coeff( 50)    *x22    *x43    
     6  +coeff( 51)    *x22    *x41*x52
     7  +coeff( 52)    *x24    *x42    
     8  +coeff( 53)*x11*x22        *x52
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 54)    *x24        *x52
     1  +coeff( 55)    *x23    *x42*x51
     2  +coeff( 56)*x11*x21    *x41*x52
     3  +coeff( 57)*x12*x24            
     4  +coeff( 58)    *x24    *x41*x52
     5  +coeff( 59)    *x24        *x53
     6  +coeff( 60)*x11*x23        *x53
     7  +coeff( 61)*x11*x24    *x41*x52
     8  +coeff( 62)*x12                
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 63)            *x42*x51
     1  +coeff( 64)            *x41*x52
     2  +coeff( 65)*x11*x21        *x51
     3  +coeff( 66)*x11        *x42    
     4  +coeff( 67)            *x43*x51
     5  +coeff( 68)*x11*x21    *x41*x51
     6  +coeff( 69)*x12*x22            
     7  +coeff( 70)*x11        *x43    
     8  +coeff( 71)*x12*x21    *x41    
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 72)*x11*x22    *x42    
     1  +coeff( 73)    *x24    *x41*x51
     2  +coeff( 74)*x11*x21    *x43    
     3  +coeff( 75)*x11*x24    *x41    
     4  +coeff( 76)*x11*x24        *x51
     5  +coeff( 77)    *x22    *x41*x53
     6  +coeff( 78)*x11*x23    *x41*x51
     7  +coeff( 79)*x11*x23        *x52
     8  +coeff( 80)    *x21    *x43*x52
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 81)    *x23    *x42*x52
     1  +coeff( 82)*x12*x22    *x42    
     2  +coeff( 83)*x11*x24    *x41*x51
     3  +coeff( 84)*x12*x22        *x52
     4  +coeff( 85)*x11*x24        *x52
     5  +coeff( 86)*x12*x21    *x42*x51
     6  +coeff( 87)*x12*x24    *x41    
     7  +coeff( 88)    *x24    *x43*x51
     8  +coeff( 89)*x11*x22    *x41*x53
      x_sp_q3ex   =x_sp_q3ex   
     9  +coeff( 90)    *x23    *x43*x52
c
      return
      end
      function t_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.4031088E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.47715148E-02,-0.17860426E-01, 0.23415944E-03, 0.11793287E+00,
     +  0.69330540E-02,-0.10603162E-01,-0.13678235E-02, 0.24538687E-02,
     + -0.14194965E-02,-0.18902665E-02,-0.97893435E-03,-0.57344855E-03,
     +  0.52677817E-03,-0.13727887E-02, 0.11256016E-02, 0.15889012E-02,
     +  0.50864165E-03, 0.10040021E-03, 0.20547293E-03,-0.69061952E-03,
     +  0.49777876E-03, 0.83788193E-03, 0.19467856E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
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
     1  +coeff( 10)    *x21    *x42    
     2  +coeff( 11)    *x21        *x52
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x23            
     5  +coeff( 14)    *x21    *x41*x51
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)    *x22    *x41*x51
      t_sp_q3ex   =t_sp_q3ex   
     9  +coeff( 18)*x11*x21            
     1  +coeff( 19)    *x23    *x41    
     2  +coeff( 20)    *x21    *x43    
     3  +coeff( 21)    *x23        *x51
     4  +coeff( 22)    *x22    *x42*x51
     5  +coeff( 23)    *x23    *x43    
c
      return
      end
      function y_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/  0.2089276E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.62500429E-03, 0.20918041E-01, 0.16048195E-01,-0.99204434E-03,
     +  0.48197466E-02, 0.34534376E-01, 0.59046312E-02,-0.33444162E-01,
     + -0.63860966E-02,-0.16809599E-02,-0.39286696E-03,-0.85243760E-02,
     +  0.99929061E-03,-0.74069074E-03,-0.57430926E-03,-0.58617294E-02,
     + -0.98549819E-03, 0.16114787E-03,-0.21247607E-01,-0.82234628E-02,
     +  0.81512704E-02,-0.83808498E-02,-0.61510918E-02, 0.16737833E-02,
     + -0.77282905E-03,-0.60653850E-03, 0.20831155E-02,-0.15742725E-02,
     + -0.25720522E-02,-0.22329474E-03, 0.12885532E-03, 0.42373539E-03,
     + -0.91356918E-03,-0.10056215E-03, 0.66417939E-03,-0.46171425E-02,
     + -0.72419061E-03,-0.20902101E-02, 0.20422356E-02,-0.13616104E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_q3ex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)    *x21    *x41    
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff(  9)    *x21        *x51
     1  +coeff( 10)            *x43    
     2  +coeff( 11)            *x42*x51
     3  +coeff( 12)    *x22            
     4  +coeff( 13)*x11        *x41    
     5  +coeff( 14)            *x41*x52
     6  +coeff( 15)                *x53
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)*x11*x21            
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 18)    *x21        *x52
     1  +coeff( 19)    *x22    *x41    
     2  +coeff( 20)    *x22        *x51
     3  +coeff( 21)    *x23            
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x21    *x44    
     7  +coeff( 25)            *x44    
     8  +coeff( 26)            *x43*x51
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)*x11*x21    *x42    
     2  +coeff( 29)    *x23        *x52
     3  +coeff( 30)*x11                
     4  +coeff( 31)*x11            *x51
     5  +coeff( 32)    *x21    *x43    
     6  +coeff( 33)            *x41*x53
     7  +coeff( 34)*x11*x21        *x51
     8  +coeff( 35)*x11*x22            
      y_sp_q3ex   =y_sp_q3ex   
     9  +coeff( 36)    *x22    *x43    
     1  +coeff( 37)*x11*x21    *x41*x51
     2  +coeff( 38)    *x22    *x42*x51
     3  +coeff( 39)    *x23    *x42    
     4  +coeff( 40)*x11*x22    *x43    
c
      return
      end
      function p_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1815546E-03/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.58655121E-03, 0.13232497E-02,-0.25753461E-01,-0.11368356E-01,
     +  0.68953177E-02, 0.18566763E-01,-0.28084770E-02, 0.35024153E-02,
     + -0.16627993E-01,-0.19760965E-02, 0.82565413E-03,-0.50850888E-02,
     +  0.11050425E-01, 0.75986428E-03, 0.29493258E-02, 0.57380428E-02,
     + -0.51283761E-03, 0.18920074E-02,-0.14620476E-02, 0.69550081E-03,
     +  0.61863574E-03, 0.28859412E-02, 0.10502596E-02, 0.42015882E-02,
     +  0.17422737E-03, 0.25286686E-05,-0.33505831E-03,-0.74349908E-03,
     +  0.96413726E-03,-0.55710919E-03,-0.19282463E-02,-0.76703541E-03,
     +  0.55345398E-03, 0.96945798E-04, 0.42833362E-03, 0.11195466E-03,
     + -0.18661185E-03,-0.16259060E-03,-0.72084695E-04, 0.16259945E-02,
     +  0.10957822E-02,-0.10633556E-02,-0.83032955E-03, 0.42765035E-03,
     + -0.73752843E-03,-0.18887380E-02, 0.15965124E-02,-0.13592870E-02,
     + -0.93956332E-03,-0.69166720E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x41    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)*x11        *x41    
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 18)            *x41*x52
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)    *x21    *x41*x51
     3  +coeff( 21)            *x42*x51
     4  +coeff( 22)    *x22    *x41*x51
     5  +coeff( 23)*x11*x21    *x42    
     6  +coeff( 24)    *x22    *x43    
     7  +coeff( 25)*x11                
     8  +coeff( 26)*x11            *x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)    *x23    *x41    
     2  +coeff( 29)    *x23        *x51
     3  +coeff( 30)            *x41*x53
     4  +coeff( 31)    *x23    *x43    
     5  +coeff( 32)    *x21        *x52
     6  +coeff( 33)*x11*x21    *x41    
     7  +coeff( 34)*x12            *x51
     8  +coeff( 35)    *x21    *x42*x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 36)*x11*x23            
     1  +coeff( 37)*x11*x22        *x51
     2  +coeff( 38)*x11*x21    *x41*x51
     3  +coeff( 39)*x11*x21        *x52
     4  +coeff( 40)    *x23        *x52
     5  +coeff( 41)    *x22        *x53
     6  +coeff( 42)*x11*x23    *x41    
     7  +coeff( 43)*x11*x23        *x51
     8  +coeff( 44)*x11*x21    *x42*x51
      p_sp_q3ex   =p_sp_q3ex   
     9  +coeff( 45)*x11*x22        *x52
     1  +coeff( 46)    *x22    *x41*x53
     2  +coeff( 47)*x11*x23    *x41*x51
     3  +coeff( 48)*x11*x23        *x52
     4  +coeff( 49)*x11*x22        *x53
     5  +coeff( 50)*x12*x23    *x41    
c
      return
      end
      function l_sp_q3ex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/ -0.5387782E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51419E-01, 0.00000E+00, 0.18179E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.34236271E-03,-0.26401848E+00,-0.31262416E-01, 0.89247422E-02,
     + -0.32926541E-01, 0.35062764E-01,-0.11912325E-01,-0.14792903E-01,
     + -0.89173466E-02,-0.23006622E-02,-0.17873172E-02, 0.18434949E-01,
     + -0.23012523E-01, 0.19057473E-02,-0.36721751E-02, 0.85712626E-03,
     + -0.20105646E-02,-0.50447825E-02,-0.11425902E-01, 0.95343235E-03,
     + -0.59386430E-03,-0.18322294E-02, 0.15404999E-02,-0.35709378E-02,
     +  0.62422856E-03, 0.95358590E-03,-0.12217276E-02,-0.15443751E-02,
     + -0.60972366E-02,-0.26132226E-01,-0.36766457E-02, 0.46789986E-02,
     +  0.34907236E-03, 0.66739821E-03, 0.79031158E-02,-0.11126164E-02,
     + -0.27234962E-02, 0.13069379E-03, 0.27844480E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)            *x42    
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x21    *x42    
     4  +coeff( 13)    *x23    *x41    
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x22    *x41    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)*x11*x22            
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)            *x41    
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)    *x22        *x51
     5  +coeff( 23)                *x53
     6  +coeff( 24)*x11*x22    *x41    
     7  +coeff( 25)*x11*x21            
     8  +coeff( 26)    *x21    *x41*x51
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 27)*x11*x21    *x41    
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)    *x22    *x43    
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)*x11*x22    *x43    
     5  +coeff( 32)*x11*x23    *x41*x52
     6  +coeff( 33)            *x42*x51
     7  +coeff( 34)*x11        *x42    
     8  +coeff( 35)    *x21    *x43    
      l_sp_q3ex   =l_sp_q3ex   
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)*x11*x22    *x42    
     2  +coeff( 38)*x11*x22        *x52
     3  +coeff( 39)*x11*x22    *x41*x52
c
      return
      end
      function x_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/  0.1179949E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12986526E-03,-0.10567220E-04,-0.29915418E-01,-0.19968271E-04,
     + -0.39090617E-04, 0.52949516E-04, 0.28759205E-04, 0.82994311E-05,
     + -0.61949931E-05,-0.11855483E-05, 0.32094135E-05,-0.73929964E-05,
     + -0.10290519E-04, 0.38495041E-05,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      x_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)    *x22    *x41    
     8  +coeff(  8)    *x21    *x41    
      x_sp_sen    =x_sp_sen    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)*x11                
     2  +coeff( 11)            *x41*x51
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x24            
     5  +coeff( 14)*x11*x21    *x41    
c
      return
      end
      function t_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.9716649E-01/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.11254171E-03,-0.80555474E-04,-0.24252066E-01,-0.24961354E-03,
     + -0.54003886E-03, 0.42976276E-03, 0.20554484E-03,-0.14445740E-03,
     + -0.24330846E-04,-0.12553362E-03, 0.18348423E-03, 0.34224442E-04,
     +  0.91741676E-04,-0.32975262E-04, 0.99505254E-04,-0.15012773E-03,
     +  0.23738110E-04,-0.55873574E-04,-0.86484019E-04,-0.83188745E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      t_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)            *x42    
      t_sp_sen    =t_sp_sen    
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x23            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)    *x21    *x42    
     5  +coeff( 14)*x11*x22            
     6  +coeff( 15)*x11*x21    *x41    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x22        *x51
      t_sp_sen    =t_sp_sen    
     9  +coeff( 18)*x11*x23            
     1  +coeff( 19)    *x23    *x42    
     2  +coeff( 20)*x11*x23    *x41    
c
      return
      end
      function y_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  6)
      data ncoeff/  5/
      data avdat/ -0.1015245E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.94468438E-03,-0.17227860E-07, 0.15608674E-05, 0.66365875E-01,
     +  0.39879410E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
c
      return
      end
      function p_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 15)
      data ncoeff/ 14/
      data avdat/ -0.8532951E-03/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.56719035E-03, 0.50793618E-01,-0.10544643E-02,-0.15845473E-03,
     + -0.98367287E-04, 0.18315678E-03, 0.36813351E-03, 0.65436921E-04,
     + -0.43348176E-04, 0.15503194E-03,-0.16648673E-03, 0.28580654E-03,
     +  0.40326311E-04, 0.68247435E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x21        *x51
      p_sp_sen    =p_sp_sen    
     9  +coeff(  9)*x11        *x41    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)    *x21    *x41*x51
     5  +coeff( 14)*x11*x22            
c
      return
      end
      function l_sp_sen    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1274153E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73294551E-03, 0.49782002E-05, 0.27576557E-02,-0.17160119E-02,
     + -0.35103349E-03, 0.17931454E-05,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x21 = x2
      x22 = x21*x2
      x41 = x4
      x42 = x41*x4
      x51 = x5
c
c                  function
c
      l_sp_sen    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)    *x22            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)                *x51
c
      return
      end
      function x_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.1790629E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.41872545E-03, 0.24351390E-03,-0.41963819E-01,-0.84961049E-03,
     +  0.83617860E-03, 0.99541433E-03, 0.27145250E-03,-0.19102635E-03,
     + -0.36448083E-03, 0.91489125E-03,-0.44140685E-03,-0.15067335E-03,
     +  0.81644233E-04, 0.20704635E-03,-0.21710174E-03, 0.57749738E-04,
     +  0.44142373E-04,-0.63115323E-04,-0.39976181E-04, 0.14774222E-03,
     + -0.63655258E-04,-0.21875856E-03, 0.48749059E-03,-0.12648918E-03,
     + -0.57508011E-03,-0.68995547E-06,-0.90231326E-04,-0.37557405E-03,
     +  0.14258396E-04, 0.23789060E-04,-0.81685039E-05,-0.89066722E-04,
     +  0.46403752E-05, 0.10992764E-03,-0.80099941E-04, 0.13918013E-03,
     +  0.68650363E-04, 0.17421842E-04,-0.13010632E-03,-0.13742522E-03,
     + -0.72441195E-04,-0.65056025E-04, 0.60243216E-04, 0.42430303E-04,
     + -0.29937761E-04, 0.36353646E-04,-0.27270602E-04, 0.15549853E-04,
     +  0.35011262E-05, 0.70825146E-04, 0.37343871E-04,-0.28461072E-04,
     +  0.20622416E-04,-0.55551723E-05,-0.19808767E-04, 0.15952509E-03,
     + -0.12565264E-03, 0.92507864E-04, 0.10219603E-03, 0.12388194E-04,
     + -0.53156273E-05, 0.44022505E-06,-0.11882915E-04,-0.82512342E-05,
     +  0.75895787E-05,-0.84208214E-05,-0.16978971E-04, 0.51590978E-05,
     +  0.11362852E-04,-0.14006022E-04,-0.12590411E-04,-0.82198549E-05,
     + -0.14411634E-04,-0.12895235E-04,-0.28051147E-04,-0.11178549E-04,
     + -0.16649639E-04,-0.11981630E-04,-0.46411842E-04, 0.96180884E-05,
     +  0.13438520E-04, 0.73679832E-04,-0.17443532E-04, 0.24674118E-04,
     + -0.71967079E-05,-0.30635227E-04, 0.71567033E-05,-0.24767543E-04,
     +  0.35314926E-04,-0.11461236E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x23            
      x_sp_sm     =x_sp_sm     
     9  +coeff(  9)    *x24            
     1  +coeff( 10)    *x22    *x42    
     2  +coeff( 11)    *x24    *x41    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x23    *x41    
     7  +coeff( 16)            *x41*x51
     8  +coeff( 17)                *x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)            *x43    
     2  +coeff( 20)*x11*x21    *x41    
     3  +coeff( 21)*x11*x23            
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)*x11*x23    *x41    
     7  +coeff( 25)    *x24    *x42    
     8  +coeff( 26)    *x21        *x51
      x_sp_sm     =x_sp_sm     
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)    *x24    *x43    
     2  +coeff( 29)*x11*x24    *x42    
     3  +coeff( 30)*x11                
     4  +coeff( 31)    *x21    *x41*x51
     5  +coeff( 32)*x11*x22            
     6  +coeff( 33)            *x41*x52
     7  +coeff( 34)    *x21    *x43    
     8  +coeff( 35)*x11*x22    *x41    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 36)*x11*x21    *x42    
     1  +coeff( 37)*x11*x24            
     2  +coeff( 38)*x11        *x43    
     3  +coeff( 39)    *x23    *x43    
     4  +coeff( 40)*x11*x23    *x42    
     5  +coeff( 41)*x11*x22    *x43    
     6  +coeff( 42)*x11*x21    *x42*x52
     7  +coeff( 43)*x12*x23    *x41*x51
     8  +coeff( 44)*x11*x24    *x43    
      x_sp_sm     =x_sp_sm     
     9  +coeff( 45)    *x21    *x42*x51
     1  +coeff( 46)    *x24        *x51
     2  +coeff( 47)            *x43*x51
     3  +coeff( 48)            *x41*x53
     4  +coeff( 49)*x11*x21        *x52
     5  +coeff( 50)    *x24    *x41*x51
     6  +coeff( 51)    *x23    *x42*x51
     7  +coeff( 52)*x11*x21    *x41*x52
     8  +coeff( 53)*x11*x23        *x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 54)*x12*x23    *x41    
     1  +coeff( 55)*x12*x24        *x51
     2  +coeff( 56)*x12*x23    *x41*x52
     3  +coeff( 57)*x12*x21    *x43*x52
     4  +coeff( 58)*x11*x22    *x43*x53
     5  +coeff( 59)*x12*x24    *x42*x52
     6  +coeff( 60)*x11        *x41    
     7  +coeff( 61)    *x21        *x52
     8  +coeff( 62)*x12                
      x_sp_sm     =x_sp_sm     
     9  +coeff( 63)*x11*x21        *x51
     1  +coeff( 64)    *x23        *x51
     2  +coeff( 65)*x11        *x41*x51
     3  +coeff( 66)*x12        *x41    
     4  +coeff( 67)            *x42*x52
     5  +coeff( 68)*x11*x21    *x41*x51
     6  +coeff( 69)    *x23    *x41*x51
     7  +coeff( 70)*x12*x22            
     8  +coeff( 71)    *x22    *x41*x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 72)*x12*x21        *x51
     1  +coeff( 73)    *x21    *x43*x51
     2  +coeff( 74)    *x21    *x41*x53
     3  +coeff( 75)*x11*x22    *x42    
     4  +coeff( 76)*x12        *x41*x51
     5  +coeff( 77)*x11*x22    *x41*x51
     6  +coeff( 78)            *x43*x52
     7  +coeff( 79)    *x23    *x41*x52
     8  +coeff( 80)*x11*x21        *x53
      x_sp_sm     =x_sp_sm     
     9  +coeff( 81)*x12*x22    *x41    
     1  +coeff( 82)*x11*x24    *x41    
     2  +coeff( 83)*x11        *x43*x51
     3  +coeff( 84)    *x22    *x43*x51
     4  +coeff( 85)*x11        *x42*x52
     5  +coeff( 86)    *x22    *x41*x53
     6  +coeff( 87)*x12*x21    *x42    
     7  +coeff( 88)*x12*x21    *x41*x51
     8  +coeff( 89)    *x21    *x43*x52
      x_sp_sm     =x_sp_sm     
     9  +coeff( 90)    *x21    *x42*x53
c
      return
      end
      function t_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1624449E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.95637824E-03, 0.86431816E-03,-0.25759574E-01,-0.35579421E-02,
     +  0.22492434E-02, 0.17317634E-02, 0.60065137E-03, 0.20468248E-03,
     + -0.52858429E-03,-0.32615641E-03, 0.12979809E-02,-0.37072896E-03,
     +  0.12801499E-03, 0.21405191E-03, 0.56379515E-03, 0.52027032E-03,
     + -0.48865756E-03, 0.46420759E-04,-0.63614927E-04,-0.73189885E-05,
     +  0.69799391E-03,-0.63393789E-03, 0.77396404E-03,-0.28273065E-03,
     + -0.87020523E-03,-0.22113070E-03,-0.51802286E-03,-0.22049077E-04,
     + -0.57789221E-04, 0.11035889E-03,-0.35366767E-04, 0.75608499E-04,
     +  0.47079060E-03,-0.93440474E-04,-0.12107458E-03,-0.16692215E-03,
     + -0.84118056E-03,-0.21984296E-03,-0.25059245E-03, 0.94862346E-03,
     +  0.36627028E-03, 0.83052815E-03,-0.11640072E-02,-0.59693831E-03,
     +  0.14616903E-04, 0.80594793E-04,-0.26692784E-04, 0.34028690E-04,
     + -0.58410053E-04, 0.54865497E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      t_sp_sm     =t_sp_sm     
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)            *x42    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)*x11*x21    *x41    
     8  +coeff( 17)*x11*x23    *x41    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x22        *x51
     2  +coeff( 20)*x11*x23            
     3  +coeff( 21)*x11*x21    *x42    
     4  +coeff( 22)    *x23    *x42    
     5  +coeff( 23)    *x22    *x43    
     6  +coeff( 24)            *x43*x52
     7  +coeff( 25)*x11*x23    *x42    
     8  +coeff( 26)*x12*x22    *x42    
      t_sp_sm     =t_sp_sm     
     9  +coeff( 27)*x12*x22    *x43    
     1  +coeff( 28)    *x21        *x51
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)            *x41*x52
     4  +coeff( 31)*x11*x21        *x51
     5  +coeff( 32)*x12        *x41    
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)    *x22    *x41*x51
     8  +coeff( 35)            *x42*x52
      t_sp_sm     =t_sp_sm     
     9  +coeff( 36)*x11*x21    *x41*x51
     1  +coeff( 37)    *x23    *x43    
     2  +coeff( 38)*x12*x23    *x41    
     3  +coeff( 39)    *x22    *x42*x53
     4  +coeff( 40)*x11*x23    *x43*x51
     5  +coeff( 41)*x11*x22    *x43*x52
     6  +coeff( 42)*x12*x23    *x43    
     7  +coeff( 43)*x11*x23    *x43*x53
     8  +coeff( 44)*x12*x23    *x43*x52
      t_sp_sm     =t_sp_sm     
     9  +coeff( 45)*x11            *x51
     1  +coeff( 46)            *x42*x51
     2  +coeff( 47)    *x21        *x52
     3  +coeff( 48)                *x53
     4  +coeff( 49)*x11*x22            
     5  +coeff( 50)*x11            *x52
c
      return
      end
      function y_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  7)
      data ncoeff/  6/
      data avdat/ -0.1419671E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12375814E-02,-0.77891826E-04, 0.13288310E-04, 0.90704642E-01,
     +  0.39470121E-02,-0.10026846E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
c
      return
      end
      function p_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 30)
      data ncoeff/ 29/
      data avdat/ -0.8581974E-03/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.39444552E-03, 0.49655240E-01,-0.36798054E-02, 0.16759750E-02,
     +  0.61067479E-03, 0.23285388E-02,-0.45775610E-03, 0.71447727E-03,
     + -0.14209021E-02,-0.18950296E-03, 0.24895702E-03, 0.13766111E-02,
     +  0.83047911E-04, 0.12708911E-03,-0.19315984E-03, 0.15327320E-03,
     +  0.43144665E-03,-0.46314566E-04, 0.33039376E-03, 0.66502536E-04,
     +  0.15504547E-04,-0.16976280E-03, 0.15584986E-03,-0.12229964E-03,
     +  0.11543134E-03,-0.15929306E-03, 0.11665253E-03,-0.90040900E-04,
     + -0.21325753E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x23    *x41    
     7  +coeff(  7)            *x41    
     8  +coeff(  8)    *x22    *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)*x11                
     2  +coeff( 11)*x11*x22            
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)    *x21        *x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11        *x41    
     7  +coeff( 16)    *x21    *x41*x51
     8  +coeff( 17)    *x22    *x42    
      p_sp_sm     =p_sp_sm     
     9  +coeff( 18)            *x42*x52
     1  +coeff( 19)*x11*x22    *x41    
     2  +coeff( 20)    *x22    *x42*x52
     3  +coeff( 21)                *x51
     4  +coeff( 22)            *x42    
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)*x11*x23            
     7  +coeff( 25)    *x22    *x43    
     8  +coeff( 26)*x11*x23    *x41    
      p_sp_sm     =p_sp_sm     
     9  +coeff( 27)*x11*x22    *x42    
     1  +coeff( 28)*x11*x21    *x41*x53
     2  +coeff( 29)*x12*x21    *x41*x53
c
      return
      end
      function l_sp_sm     (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/ -0.1875851E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10341090E-02,-0.25060746E-04, 0.42738244E-02, 0.11528481E-03,
     + -0.24196862E-02,-0.47784348E-03,-0.25713744E-04,-0.38276707E-04,
     + -0.27809718E-04, 0.20855836E-04,-0.61618994E-05,-0.12735727E-04,
     + -0.16118365E-04, 0.21037529E-04,-0.30942763E-04,-0.69387997E-05,
     + -0.57234065E-05, 0.47039825E-05, 0.79077809E-05, 0.15574085E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_sm     =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_sm     =l_sp_sm     
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)    *x23            
     2  +coeff( 11)                *x52
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x11*x21    *x41    
     5  +coeff( 14)    *x23    *x41    
     6  +coeff( 15)    *x22    *x42    
     7  +coeff( 16)    *x21    *x42    
     8  +coeff( 17)            *x43    
      l_sp_sm     =l_sp_sm     
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x23            
     2  +coeff( 20)*x11*x23    *x41    
c
      return
      end
      function x_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.2741545E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12635069E-02, 0.86288073E-03,-0.54637704E-01,-0.34126665E-02,
     +  0.34458342E-02, 0.25412063E-02, 0.68535120E-03, 0.31935587E-03,
     + -0.54658140E-03,-0.11071911E-02, 0.23424951E-02,-0.11730590E-02,
     + -0.47031883E-03, 0.18081391E-03, 0.51118806E-03,-0.55446109E-03,
     +  0.15560605E-03,-0.21818513E-03, 0.25770042E-03,-0.19842171E-03,
     + -0.54135820E-03, 0.11702139E-02,-0.17358171E-03,-0.14226430E-02,
     +  0.64434265E-04, 0.12966428E-05,-0.19748167E-03, 0.30999476E-03,
     + -0.29809409E-03,-0.78409488E-04, 0.11773873E-04,-0.21207395E-03,
     +  0.23383144E-03,-0.28562467E-03, 0.11104241E-03,-0.75191138E-05,
     + -0.30944229E-03,-0.85321342E-03, 0.56897665E-04,-0.12315289E-03,
     + -0.42347543E-04, 0.97294331E-04,-0.72634983E-04,-0.49219354E-04,
     +  0.36259520E-04, 0.10418652E-03, 0.13946752E-03,-0.41030467E-04,
     + -0.60417173E-04, 0.23928849E-03, 0.65184468E-04,-0.10890573E-03,
     +  0.26344109E-03,-0.29614122E-03,-0.34803768E-04,-0.11610422E-03,
     + -0.71887232E-04,-0.37776557E-04,-0.37598363E-05,-0.19663337E-04,
     + -0.10103000E-04,-0.12017749E-04,-0.57704565E-04,-0.25549358E-04,
     +  0.17356694E-03, 0.23927834E-04, 0.14575476E-03, 0.24226560E-04,
     + -0.93858929E-04, 0.75149146E-04, 0.47729541E-05, 0.13507674E-04,
     + -0.48986847E-04, 0.39613355E-04,-0.26205656E-04, 0.48596958E-04,
     +  0.13680174E-03, 0.54041044E-04,-0.66825091E-05, 0.67426256E-04,
     +  0.39561182E-04, 0.57620102E-04,-0.39179577E-04, 0.80743928E-04,
     + -0.26919673E-04, 0.34365901E-04,-0.22398149E-03, 0.15038677E-03,
     + -0.37434096E-04,-0.17108044E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      x_sp_sex    =x_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)    *x24            
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)    *x24    *x41    
     4  +coeff( 13)            *x42    
     5  +coeff( 14)                *x52
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)            *x41*x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 18)    *x22        *x51
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)*x11*x23            
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)*x11*x23    *x41    
     6  +coeff( 24)    *x24    *x42    
     7  +coeff( 25)*x11                
     8  +coeff( 26)    *x21        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 27)*x11*x22            
     1  +coeff( 28)*x11*x21    *x42    
     2  +coeff( 29)*x11*x23    *x42    
     3  +coeff( 30)    *x21    *x41*x51
     4  +coeff( 31)            *x41*x52
     5  +coeff( 32)    *x22    *x41*x51
     6  +coeff( 33)    *x21    *x43    
     7  +coeff( 34)*x11*x22    *x41    
     8  +coeff( 35)*x11*x24            
      x_sp_sex    =x_sp_sex    
     9  +coeff( 36)*x11*x22    *x42    
     1  +coeff( 37)    *x23    *x43    
     2  +coeff( 38)    *x24    *x43    
     3  +coeff( 39)*x11        *x41    
     4  +coeff( 40)            *x43    
     5  +coeff( 41)    *x23        *x51
     6  +coeff( 42)    *x24        *x51
     7  +coeff( 43)            *x43*x51
     8  +coeff( 44)            *x42*x52
      x_sp_sex    =x_sp_sex    
     9  +coeff( 45)            *x41*x53
     1  +coeff( 46)    *x23    *x41*x51
     2  +coeff( 47)    *x24    *x41*x51
     3  +coeff( 48)            *x43*x52
     4  +coeff( 49)*x11*x21    *x41*x52
     5  +coeff( 50)*x11*x24    *x41    
     6  +coeff( 51)    *x22    *x42*x52
     7  +coeff( 52)*x11*x21    *x42*x52
     8  +coeff( 53)*x12*x23    *x41*x52
      x_sp_sex    =x_sp_sex    
     9  +coeff( 54)*x12*x21    *x43*x52
     1  +coeff( 55)*x11*x21        *x51
     2  +coeff( 56)    *x21    *x42*x51
     3  +coeff( 57)    *x21    *x41*x52
     4  +coeff( 58)*x12        *x41    
     5  +coeff( 59)*x12            *x51
     6  +coeff( 60)*x11*x21    *x41*x51
     7  +coeff( 61)*x11*x21        *x52
     8  +coeff( 62)    *x22        *x53
      x_sp_sex    =x_sp_sex    
     9  +coeff( 63)*x12*x21        *x51
     1  +coeff( 64)*x12        *x41*x51
     2  +coeff( 65)*x11*x21    *x43    
     3  +coeff( 66)*x11*x21    *x42*x51
     4  +coeff( 67)    *x23    *x42*x51
     5  +coeff( 68)*x11*x21        *x53
     6  +coeff( 69)*x12*x22    *x41    
     7  +coeff( 70)    *x22    *x43*x51
     8  +coeff( 71)*x11        *x42*x52
      x_sp_sex    =x_sp_sex    
     9  +coeff( 72)*x11        *x41*x53
     1  +coeff( 73)    *x22    *x41*x53
     2  +coeff( 74)*x12*x21    *x42    
     3  +coeff( 75)*x12*x21        *x52
     4  +coeff( 76)*x11*x23        *x52
     5  +coeff( 77)    *x21    *x43*x52
     6  +coeff( 78)*x12        *x43    
     7  +coeff( 79)*x11*x22    *x43    
     8  +coeff( 80)*x12*x23        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 81)    *x23    *x42*x52
     1  +coeff( 82)*x11*x21    *x41*x53
     2  +coeff( 83)    *x23    *x41*x53
     3  +coeff( 84)*x11*x24        *x52
     4  +coeff( 85)*x11        *x43*x52
     5  +coeff( 86)*x12*x21    *x43    
     6  +coeff( 87)*x11*x23    *x43    
     7  +coeff( 88)*x12*x24    *x41    
     8  +coeff( 89)*x12*x24        *x51
      x_sp_sex    =x_sp_sex    
     9  +coeff( 90)*x11*x22    *x42*x52
c
      return
      end
      function t_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.2259666E+00/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22448760E-02, 0.17288085E-02,-0.26710913E-01,-0.69060116E-02,
     +  0.58903890E-02, 0.17755352E-02, 0.83520665E-03, 0.48364021E-03,
     + -0.11576822E-02, 0.38524190E-03, 0.16625368E-02,-0.64541215E-04,
     + -0.59726374E-03, 0.11512393E-03, 0.62210025E-03,-0.85887252E-04,
     +  0.66032860E-03,-0.72584959E-03,-0.17099232E-03,-0.71582652E-03,
     + -0.15590014E-03, 0.10276820E-03,-0.16028169E-03, 0.95218234E-03,
     + -0.66001265E-03, 0.10288426E-02,-0.40127482E-03,-0.54447871E-03,
     + -0.19219480E-04,-0.11323089E-02,-0.51124726E-03, 0.11694485E-02,
     + -0.10128248E-02, 0.10701626E-03,-0.42816617E-04, 0.21416230E-03,
     +  0.15233314E-03, 0.50950181E-04,-0.18279652E-03,-0.65246801E-04,
     +  0.74283518E-04, 0.12978043E-03, 0.13167430E-03, 0.11959457E-03,
     + -0.84362124E-04,-0.93550807E-04, 0.80797123E-04, 0.40009145E-04,
     + -0.95615615E-04,-0.34878490E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x22    *x41    
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)*x11*x21            
      t_sp_sex    =t_sp_sex    
     9  +coeff(  9)    *x23            
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x42    
     3  +coeff( 12)            *x42*x52
     4  +coeff( 13)            *x42    
     5  +coeff( 14)            *x41*x51
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)*x11*x21    *x41    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 18)    *x23    *x41    
     1  +coeff( 19)*x11*x23            
     2  +coeff( 20)*x11*x23    *x41    
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)*x11            *x52
     5  +coeff( 23)    *x21    *x42*x51
     6  +coeff( 24)*x11*x21    *x42    
     7  +coeff( 25)    *x23    *x42    
     8  +coeff( 26)    *x22    *x43    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 27)            *x43*x52
     1  +coeff( 28)*x11*x22        *x52
     2  +coeff( 29)*x12*x22    *x41    
     3  +coeff( 30)*x11*x23    *x42    
     4  +coeff( 31)    *x22    *x42*x53
     5  +coeff( 32)*x11*x22    *x42*x52
     6  +coeff( 33)*x12*x22    *x43    
     7  +coeff( 34)*x11                
     8  +coeff( 35)*x11        *x41    
      t_sp_sex    =t_sp_sex    
     9  +coeff( 36)            *x42*x51
     1  +coeff( 37)            *x41*x52
     2  +coeff( 38)                *x53
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)*x11*x21        *x51
     5  +coeff( 41)*x12*x21            
     6  +coeff( 42)*x12        *x41    
     7  +coeff( 43)    *x21    *x43    
     8  +coeff( 44)    *x23        *x51
      t_sp_sex    =t_sp_sex    
     9  +coeff( 45)    *x21    *x41*x52
     1  +coeff( 46)    *x21        *x53
     2  +coeff( 47)            *x41*x53
     3  +coeff( 48)*x11            *x53
     4  +coeff( 49)*x12*x22            
     5  +coeff( 50)    *x22    *x42*x51
c
      return
      end
      function y_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 10)
      data ncoeff/  9/
      data avdat/ -0.1852641E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.12859759E-02,-0.18188675E-03, 0.29027566E-04, 0.11363956E+00,
     +  0.38901742E-02,-0.35011182E-02, 0.68296766E-03, 0.18543696E-02,
     +  0.20657571E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x51 = x5
c
c                  function
c
      y_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      y_sp_sex    =y_sp_sex    
     9  +coeff(  9)    *x23    *x41    
c
      return
      end
      function p_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 31)
      data ncoeff/ 30/
      data avdat/ -0.1063462E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.35744937E-03, 0.47568046E-01,-0.42276108E-02, 0.16149377E-02,
     +  0.25674829E-02,-0.53457444E-03, 0.61357929E-03, 0.88637817E-03,
     + -0.28721563E-03, 0.21459232E-03,-0.15085072E-02, 0.14201850E-02,
     +  0.28745527E-03, 0.33199602E-04, 0.99134988E-04,-0.22678013E-03,
     +  0.13268647E-03, 0.45717543E-03,-0.66486944E-03,-0.56059929E-04,
     +  0.41866335E-03, 0.81187766E-03, 0.53503038E-04,-0.17825339E-03,
     +  0.17351232E-03,-0.11077688E-03, 0.76486736E-04, 0.11270749E-03,
     + -0.16767849E-03,-0.12687256E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)    *x21    *x41    
     4  +coeff(  4)    *x23            
     5  +coeff(  5)    *x23    *x41    
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x22    *x41    
      p_sp_sex    =p_sp_sex    
     9  +coeff(  9)*x11                
     1  +coeff( 10)    *x21        *x51
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x42    
     4  +coeff( 13)*x11*x22            
     5  +coeff( 14)                *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11        *x41    
     8  +coeff( 17)    *x21    *x41*x51
      p_sp_sex    =p_sp_sex    
     9  +coeff( 18)    *x22    *x42    
     1  +coeff( 19)    *x21    *x43    
     2  +coeff( 20)            *x42*x52
     3  +coeff( 21)*x11*x22    *x41    
     4  +coeff( 22)    *x23    *x43    
     5  +coeff( 23)    *x22    *x42*x52
     6  +coeff( 24)            *x42    
     7  +coeff( 25)*x11*x21    *x41    
     8  +coeff( 26)*x11*x23            
      p_sp_sex    =p_sp_sex    
     9  +coeff( 27)*x11*x21    *x42    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)*x11*x23    *x41    
     3  +coeff( 30)*x11*x21    *x41*x53
c
      return
      end
      function l_sp_sex    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 37)
      data ncoeff/ 36/
      data avdat/ -0.2791778E-02/
      data xmin/
     1 -0.39997E-02,-0.52018E-01, 0.00000E+00,-0.27908E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51938E-01, 0.00000E+00, 0.18551E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.14539291E-02,-0.13385137E-03, 0.67096828E-02, 0.61908388E-03,
     + -0.34270466E-02,-0.58458175E-03,-0.10971865E-03,-0.10888673E-03,
     + -0.62348699E-04,-0.39499919E-04,-0.40136718E-04, 0.72229341E-04,
     + -0.14000913E-03,-0.56794266E-04, 0.36445806E-04, 0.38032995E-04,
     +  0.75118519E-05,-0.36612125E-05,-0.55042554E-04, 0.73827623E-05,
     +  0.51155395E-04,-0.80256264E-04, 0.18863284E-04, 0.53333068E-04,
     + -0.85455231E-05, 0.15348789E-04, 0.58718624E-05,-0.44900771E-04,
     +  0.12368737E-04,-0.73477313E-04,-0.22489183E-04, 0.51814623E-04,
     +  0.86266977E-04, 0.21623315E-04, 0.16649918E-04, 0.68361238E-04,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
c
c                  function
c
      l_sp_sex    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)            *x42    
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)    *x22    *x41    
      l_sp_sex    =l_sp_sex    
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)                *x52
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)    *x23            
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21        *x51
      l_sp_sex    =l_sp_sex    
     9  +coeff( 18)            *x43    
     1  +coeff( 19)*x11*x21    *x41    
     2  +coeff( 20)*x11*x23            
     3  +coeff( 21)    *x23    *x42    
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)            *x43*x52
     6  +coeff( 24)*x11*x23    *x41    
     7  +coeff( 25)*x11                
     8  +coeff( 26)*x11*x22            
      l_sp_sex    =l_sp_sex    
     9  +coeff( 27)*x11*x21        *x51
     1  +coeff( 28)    *x21    *x43    
     2  +coeff( 29)    *x22    *x41*x51
     3  +coeff( 30)*x11*x21    *x42    
     4  +coeff( 31)*x12*x22    *x41    
     5  +coeff( 32)    *x23    *x43    
     6  +coeff( 33)*x11*x23    *x42    
     7  +coeff( 34)*x12*x22    *x42    
     8  +coeff( 35)*x12        *x43*x51
      l_sp_sex    =l_sp_sex    
     9  +coeff( 36)*x12*x22    *x43    
c
      return
      end
      function x_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 11)
      data ncoeff/ 10/
      data avdat/ -0.1178958E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23279034E-02, 0.69028623E-01, 0.70177950E-01, 0.54796779E-03,
     +  0.23404500E-02,-0.30422109E-03,-0.31480953E-03,-0.14190061E-03,
     +  0.22051652E-04, 0.84731008E-04,
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
      x_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)*x11                
     4  +coeff(  4)    *x21        *x51
     5  +coeff(  5)*x11            *x51
     6  +coeff(  6)    *x23            
     7  +coeff(  7)    *x21*x31*x41    
     8  +coeff(  8)    *x21        *x52
      x_sp_cq1x   =x_sp_cq1x   
     9  +coeff(  9)                *x51
     1  +coeff( 10)    *x21*x31        
c
      return
      end
      function t_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 22)
      data ncoeff/ 21/
      data avdat/ -0.2296385E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.23178289E-03, 0.34202360E-01,-0.35257199E-04,-0.44595029E-01,
     +  0.44348562E-04,-0.33143046E-04, 0.55217155E-03, 0.12512024E-04,
     +  0.20883109E-02,-0.14903080E-05,-0.11238839E-03,-0.12738830E-03,
     +  0.43741173E-04,-0.17607704E-03, 0.11620803E-02, 0.19958668E-03,
     + -0.41353254E-04, 0.20028750E-03,-0.42379345E-03,-0.19945168E-03,
     + -0.12857181E-02,
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
      t_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21*x31        
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x51
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff(  9)*x11            *x51
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)            *x41    
     5  +coeff( 14)    *x23            
     6  +coeff( 15)    *x23*x32        
     7  +coeff( 16)*x12*x21    *x41    
     8  +coeff( 17)    *x21*x32*x42*x51
      t_sp_cq1x   =t_sp_cq1x   
     9  +coeff( 18)    *x21    *x43    
     1  +coeff( 19)    *x23*x31*x41    
     2  +coeff( 20)    *x21*x33*x41    
     3  +coeff( 21)*x11*x22*x32        
c
      return
      end
      function y_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff(  5)
      data ncoeff/  4/
      data avdat/ -0.9726789E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.97567430E-02, 0.71794532E-01, 0.52567977E-01,-0.14819610E-02,
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
c
      return
      end
      function p_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 13)
      data ncoeff/ 12/
      data avdat/ -0.4125908E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.36679730E-02, 0.36403010E-05, 0.21243203E-01, 0.32554373E-01,
     + -0.10782221E-02, 0.23118513E-03,-0.21361593E-03,-0.43441812E-03,
     +  0.75386488E-04, 0.15394075E-04, 0.72749513E-04,-0.27677312E-03,
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
     7  +coeff(  7)    *x22*x31        
     8  +coeff(  8)            *x41*x51
      p_sp_cq1x   =p_sp_cq1x   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)            *x42*x51
     2  +coeff( 11)        *x31    *x52
     3  +coeff( 12)*x11*x23    *x41    
c
      return
      end
      function l_sp_cq1x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 42)
      data ncoeff/ 41/
      data avdat/ -0.5100312E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.47756752E-03, 0.27182582E-03, 0.17973357E-03, 0.27314358E-04,
     +  0.33856782E-05,-0.15080121E-02,-0.18202079E-04,-0.72198833E-03,
     +  0.17119064E-04, 0.54792478E-03,-0.15098462E-02, 0.45797438E-04,
     +  0.80532445E-04,-0.34055291E-03, 0.20686258E-02,-0.11049407E-02,
     + -0.16211432E-04, 0.37829936E-03,-0.13746567E-04,-0.46627360E-04,
     +  0.13154068E-04,-0.25917001E-04, 0.43054777E-04,-0.13453490E-02,
     + -0.89057052E-04, 0.25867799E-04,-0.83080886E-06, 0.22340888E-03,
     + -0.28971193E-03, 0.25175787E-02, 0.57973234E-05,-0.10808886E-04,
     + -0.11027701E-02, 0.10658231E-04, 0.12417911E-04,-0.34658566E-04,
     + -0.10567222E-04, 0.14969837E-04,-0.90013682E-05,-0.12721206E-04,
     +  0.24600145E-04,
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
      l_sp_cq1x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)        *x32    *x51
     4  +coeff( 13)        *x31*x41*x51
     5  +coeff( 14)        *x31    *x51
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x12                
     8  +coeff( 17)    *x21            
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 18)            *x41*x51
     1  +coeff( 19)        *x31*x42*x51
     2  +coeff( 20)                *x52
     3  +coeff( 21)    *x22*x31        
     4  +coeff( 22)        *x33        
     5  +coeff( 23)            *x43    
     6  +coeff( 24)    *x22        *x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)            *x41*x52
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 27)                *x53
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)*x11*x21        *x51
     4  +coeff( 31)        *x33*x41    
     5  +coeff( 32)        *x32*x42    
     6  +coeff( 33)*x12            *x51
     7  +coeff( 34)    *x22*x31    *x51
     8  +coeff( 35)        *x33    *x51
      l_sp_cq1x   =l_sp_cq1x   
     9  +coeff( 36)    *x22    *x41*x51
     1  +coeff( 37)*x11*x23            
     2  +coeff( 38)*x11*x22*x31        
     3  +coeff( 39)*x11*x21*x32        
     4  +coeff( 40)*x11*x22    *x41    
     5  +coeff( 41)*x11*x21*x31*x41    
c
      return
      end
      function x_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/ -0.5097530E+01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.69499807E-03,-0.10571390E+00, 0.30871115E-02, 0.12849805E-01,
     + -0.69489330E-02, 0.19715701E-02,-0.43440107E-02,-0.52349479E-02,
     +  0.62422631E-02,-0.81493137E-02,-0.61762112E-03,-0.83780556E-03,
     +  0.20382439E-02, 0.23262178E-03, 0.61069598E-03,-0.31460606E-03,
     +  0.65217813E-03,-0.86024258E-04,-0.31945196E-02,-0.33362249E-02,
     + -0.11949282E-02, 0.11582772E-03,-0.11577542E-02, 0.67257154E-03,
     +  0.20542603E-02,-0.70692436E-03, 0.30071431E-03, 0.45684527E-03,
     + -0.16030094E-02,-0.93866978E-02,-0.21167405E-03,-0.15019755E-03,
     +  0.38140509E-03,-0.16947386E-03, 0.28824634E-02,-0.44597851E-03,
     +  0.12547377E-02, 0.26501340E-03,-0.39755451E-03,-0.15146365E-03,
     +  0.43763773E-03,-0.16628893E-02, 0.12027519E-02, 0.98891545E-03,
     +  0.19432332E-02, 0.96909789E-04,-0.73378196E-03,-0.66678499E-05,
     + -0.13172836E-03,-0.12587713E-02, 0.15546411E-04, 0.42259321E-02,
     + -0.96810708E-03, 0.99576231E-04, 0.25415921E-03, 0.35498049E-02,
     +  0.68471345E-04,-0.53804837E-04, 0.39652019E-03,-0.88758410E-04,
     +  0.26695645E-04,-0.15789429E-03,-0.42118867E-04, 0.34536817E-04,
     +  0.49191993E-04, 0.60263075E-04,-0.24422107E-03,-0.29671166E-03,
     +  0.13815094E-03,-0.19181902E-02, 0.56566880E-03, 0.14374753E-03,
     + -0.36223566E-04,-0.12950100E-03, 0.10550899E-03,-0.19671589E-02,
     + -0.68587641E-03,-0.40758701E-03,-0.10481178E-03, 0.13155093E-03,
     +  0.31033374E-03, 0.84908510E-03,-0.31882705E-03, 0.18592909E-03,
     +  0.72863436E-03,-0.19628328E-03, 0.12489119E-03,-0.60167362E-03,
     +  0.85483619E-03, 0.56971330E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x24 = x23*x2
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
     4  +coeff(  4)    *x21    *x41    
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)            *x41    
     7  +coeff(  7)    *x23            
     8  +coeff(  8)    *x22    *x41    
      x_sp_cden   =x_sp_cden   
     9  +coeff(  9)    *x21    *x42    
     1  +coeff( 10)    *x23    *x41    
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)*x11*x22            
     4  +coeff( 13)    *x24    *x41    
     5  +coeff( 14)*x11        *x43    
     6  +coeff( 15)            *x42    
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21        *x52
      x_sp_cden   =x_sp_cden   
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x22    *x42    
     2  +coeff( 20)    *x23    *x42    
     3  +coeff( 21)*x11*x22    *x42    
     4  +coeff( 22)    *x24        *x52
     5  +coeff( 23)    *x22            
     6  +coeff( 24)*x11        *x41    
     7  +coeff( 25)    *x24            
     8  +coeff( 26)*x11*x21    *x41    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 27)*x11        *x42    
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)*x11*x22    *x41    
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)                *x51
     5  +coeff( 32)            *x41*x51
     6  +coeff( 33)    *x22        *x51
     7  +coeff( 34)    *x21    *x41*x51
     8  +coeff( 35)    *x21    *x43    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 36)    *x21    *x42*x51
     1  +coeff( 37)*x11*x24    *x41    
     2  +coeff( 38)            *x43    
     3  +coeff( 39)    *x23        *x51
     4  +coeff( 40)*x12*x21            
     5  +coeff( 41)    *x23    *x41*x51
     6  +coeff( 42)    *x22    *x43    
     7  +coeff( 43)*x11*x23    *x41    
     8  +coeff( 44)    *x21    *x43*x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 45)    *x24    *x42    
     1  +coeff( 46)*x11*x21    *x43    
     2  +coeff( 47)    *x23    *x42*x51
     3  +coeff( 48)*x11        *x42*x52
     4  +coeff( 49)*x11*x22    *x42*x51
     5  +coeff( 50)*x11*x23    *x43    
     6  +coeff( 51)*x11*x23    *x41*x52
     7  +coeff( 52)*x12*x24    *x41    
     8  +coeff( 53)*x12*x24        *x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 54)*x12*x22    *x43    
     1  +coeff( 55)*x12*x21    *x42*x52
     2  +coeff( 56)*x11*x22    *x43*x52
     3  +coeff( 57)                *x53
     4  +coeff( 58)*x11*x21        *x51
     5  +coeff( 59)    *x22    *x41*x51
     6  +coeff( 60)    *x21    *x41*x52
     7  +coeff( 61)*x12            *x51
     8  +coeff( 62)*x11*x21    *x42    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 63)*x11*x21    *x41*x51
     1  +coeff( 64)*x11        *x42*x51
     2  +coeff( 65)*x11        *x41*x52
     3  +coeff( 66)    *x22    *x41*x52
     4  +coeff( 67)*x12*x21    *x41    
     5  +coeff( 68)    *x24    *x41*x51
     6  +coeff( 69)    *x23        *x53
     7  +coeff( 70)*x12*x22    *x41    
     8  +coeff( 71)*x12*x22        *x51
      x_sp_cden   =x_sp_cden   
     9  +coeff( 72)*x11*x23    *x41*x51
     1  +coeff( 73)*x12*x21        *x52
     2  +coeff( 74)*x12*x24            
     3  +coeff( 75)*x12        *x43    
     4  +coeff( 76)*x11*x22    *x43    
     5  +coeff( 77)*x11*x22    *x41*x52
     6  +coeff( 78)    *x24    *x41*x52
     7  +coeff( 79)*x12            *x53
     8  +coeff( 80)*x12*x23    *x41    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 81)*x12*x23        *x51
     1  +coeff( 82)    *x23    *x41*x53
     2  +coeff( 83)*x11        *x43*x52
     3  +coeff( 84)*x12*x21    *x43    
     4  +coeff( 85)*x11*x23    *x42*x51
     5  +coeff( 86)*x12*x21        *x53
     6  +coeff( 87)*x11*x23        *x53
     7  +coeff( 88)    *x21    *x43*x53
     8  +coeff( 89)*x12*x23    *x42    
      x_sp_cden   =x_sp_cden   
     9  +coeff( 90)*x11*x21    *x43*x52
c
      return
      end
      function t_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 40)
      data ncoeff/ 39/
      data avdat/  0.3603638E+01/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.41653894E-01, 0.61587846E+00, 0.75405012E-02,-0.26042875E-01,
     +  0.17275757E+00,-0.94125852E-01, 0.43867867E-01, 0.69286987E-01,
     + -0.10068014E-02,-0.17872421E-01, 0.12134173E-01, 0.80567542E-02,
     + -0.69379732E-02,-0.44120770E-01, 0.18911004E-01, 0.58756132E-01,
     + -0.47258181E-02, 0.24379414E-01,-0.42992570E-02, 0.21755805E-02,
     +  0.26804563E-02, 0.86496538E-02, 0.99820830E-02,-0.12971380E-02,
     +  0.45698094E-02, 0.19297046E-02,-0.57454794E-02,-0.31754598E-02,
     +  0.27595987E-02, 0.65209411E-01,-0.21046434E-01, 0.38078502E-02,
     + -0.15565638E-02,-0.42118491E-02,-0.40964084E-02, 0.56080287E-03,
     +  0.61563626E-02,-0.23743336E-02, 0.63098604E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)    *x23            
      t_sp_cden   =t_sp_cden   
     9  +coeff(  9)            *x43    
     1  +coeff( 10)            *x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)            *x41*x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)    *x22        *x51
     7  +coeff( 16)    *x23    *x41    
     8  +coeff( 17)    *x21        *x52
      t_sp_cden   =t_sp_cden   
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11        *x41    
     2  +coeff( 20)*x11            *x51
     3  +coeff( 21)*x11*x22            
     4  +coeff( 22)    *x23        *x51
     5  +coeff( 23)*x11*x22    *x41    
     6  +coeff( 24)    *x21    *x41*x51
     7  +coeff( 25)            *x42*x51
     8  +coeff( 26)*x11*x21    *x41    
      t_sp_cden   =t_sp_cden   
     9  +coeff( 27)    *x22    *x41*x51
     1  +coeff( 28)*x11*x23            
     2  +coeff( 29)            *x43*x52
     3  +coeff( 30)    *x23    *x43    
     4  +coeff( 31)    *x21    *x43    
     5  +coeff( 32)    *x21    *x42*x51
     6  +coeff( 33)    *x22        *x52
     7  +coeff( 34)    *x23    *x41*x51
     8  +coeff( 35)    *x22    *x42*x51
      t_sp_cden   =t_sp_cden   
     9  +coeff( 36)*x12            *x52
     1  +coeff( 37)*x11*x22    *x42    
     2  +coeff( 38)*x11        *x42*x52
     3  +coeff( 39)*x11*x23    *x43    
c
      return
      end
      function y_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 20)
      data ncoeff/ 19/
      data avdat/ -0.7348734E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.44901045E-02, 0.99475138E-01, 0.22400493E-01,-0.63730804E-02,
     +  0.56264075E-02, 0.11914965E-01,-0.26342151E-03,-0.15184742E-01,
     + -0.14983970E-02,-0.90736682E-02, 0.41951705E-02, 0.13742591E-02,
     + -0.13790208E-01,-0.81911578E-03,-0.19263171E-02,-0.21686191E-02,
     + -0.13823749E-02,-0.10150283E-02, 0.37228724E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
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
      y_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)    *x21            
     5  +coeff(  5)            *x42    
     6  +coeff(  6)            *x41*x51
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x22            
      y_sp_cden   =y_sp_cden   
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x23            
     3  +coeff( 12)                *x54
     4  +coeff( 13)    *x22    *x42    
     5  +coeff( 14)            *x43    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x44    
     8  +coeff( 17)    *x22        *x51
      y_sp_cden   =y_sp_cden   
     9  +coeff( 18)*x11*x21    *x41    
     1  +coeff( 19)    *x23    *x41    
c
      return
      end
      function p_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.5418017E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.25268176E-02, 0.38742465E-02,-0.55461004E-01,-0.59286216E-02,
     +  0.61736177E-02,-0.18566305E-01,-0.29717816E-04,-0.37588002E-02,
     +  0.12132913E-01, 0.24796359E-02,-0.94683133E-02,-0.36140098E-02,
     +  0.32209881E-03, 0.73140085E-03, 0.32282320E-02,-0.10602347E-02,
     +  0.36975616E-02, 0.14668208E-02,-0.10758217E-02,-0.34381705E-03,
     +  0.96127979E-03,-0.17082912E-03,-0.74099359E-03, 0.16265608E-03,
     + -0.30013444E-03, 0.21756240E-02,-0.23243048E-02, 0.15802846E-02,
     + -0.13693585E-02, 0.10605219E-03,-0.86288404E-04,-0.10180486E-03,
     +  0.37663389E-03, 0.14936439E-03,-0.10384311E-03,-0.25409847E-03,
     + -0.94329217E-03,-0.25378354E-02, 0.78017532E-04,-0.70071736E-04,
     +  0.19229551E-03,-0.16741667E-04, 0.39735064E-03, 0.17907593E-03,
     + -0.25196481E-03,-0.11677594E-02, 0.30445692E-03,-0.47576541E-03,
     +  0.53170905E-03, 0.71945647E-03,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      p_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)    *x22            
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)            *x42    
     8  +coeff(  8)    *x21        *x51
      p_sp_cden   =p_sp_cden   
     9  +coeff(  9)            *x41*x51
     1  +coeff( 10)                *x52
     2  +coeff( 11)    *x22    *x41    
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)*x11*x21            
     5  +coeff( 14)*x11        *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x41*x52
     8  +coeff( 17)    *x22    *x42    
      p_sp_cden   =p_sp_cden   
     9  +coeff( 18)    *x23            
     1  +coeff( 19)            *x43    
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)*x11*x21    *x41    
     4  +coeff( 22)    *x23    *x41    
     5  +coeff( 23)    *x22    *x41*x51
     6  +coeff( 24)    *x21        *x52
     7  +coeff( 25)                *x53
     8  +coeff( 26)    *x21    *x43    
      p_sp_cden   =p_sp_cden   
     9  +coeff( 27)    *x23    *x42    
     1  +coeff( 28)    *x22    *x43    
     2  +coeff( 29)    *x22    *x42*x51
     3  +coeff( 30)*x11            *x51
     4  +coeff( 31)*x11*x22            
     5  +coeff( 32)*x11        *x41*x51
     6  +coeff( 33)    *x23        *x51
     7  +coeff( 34)    *x22        *x52
     8  +coeff( 35)*x11*x23            
      p_sp_cden   =p_sp_cden   
     9  +coeff( 36)*x11*x22    *x41    
     1  +coeff( 37)*x11*x23    *x41    
     2  +coeff( 38)    *x23    *x43    
     3  +coeff( 39)*x11        *x42    
     4  +coeff( 40)*x11*x21        *x51
     5  +coeff( 41)    *x21    *x41*x52
     6  +coeff( 42)            *x41*x53
     7  +coeff( 43)*x11*x21    *x42    
     8  +coeff( 44)            *x42*x53
      p_sp_cden   =p_sp_cden   
     9  +coeff( 45)*x12*x22    *x41    
     1  +coeff( 46)    *x22    *x43*x51
     2  +coeff( 47)            *x43*x53
     3  +coeff( 48)*x11*x22    *x41*x52
     4  +coeff( 49)*x12*x22    *x41*x51
     5  +coeff( 50)    *x22    *x43*x52
c
      return
      end
      function l_sp_cden   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 50)
      data ncoeff/ 49/
      data avdat/ -0.2278157E-02/
      data xmin/
     1 -0.39997E-02,-0.50464E-01, 0.00000E+00,-0.24130E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.39960E-02, 0.51179E-01, 0.00000E+00, 0.18179E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.51585832E-02, 0.20447290E+00,-0.68634884E-02,-0.59955809E-02,
     + -0.24418786E-01, 0.13927435E-01,-0.95742205E-02, 0.82017137E-02,
     + -0.56372047E-02, 0.66158632E-02,-0.12146716E-01, 0.15258301E-01,
     + -0.14507415E-02, 0.11037112E-02, 0.14679035E-02, 0.64611072E-02,
     + -0.30588533E-03, 0.77636275E-02, 0.65789890E-03, 0.58396929E-03,
     + -0.12920818E-02, 0.89207236E-02,-0.32589567E-04, 0.33565145E-02,
     + -0.21015985E-03,-0.12947259E-02,-0.14255648E-02,-0.60919626E-03,
     +  0.21124879E-03, 0.14598644E-02,-0.56535110E-03, 0.19573784E-02,
     +  0.11362380E-02,-0.88227965E-03, 0.16668793E-01, 0.24813847E-03,
     +  0.22724843E-03,-0.50082086E-02,-0.70450806E-04,-0.14866653E-02,
     +  0.34253139E-03,-0.78132330E-03, 0.95232004E-04, 0.26346687E-02,
     + -0.11885086E-02, 0.16762946E-02, 0.41148295E-02,-0.36757353E-02,
     + -0.26984522E-02,
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
      x4 =1.+(x(  4)-xmax(  4))*scale(  4)
      x5 =1.+(x(  5)-xmax(  5))*scale(  5)
c          set up monomials   functions
      x11 = x1
      x12 = x11*x1
      x21 = x2
      x22 = x21*x2
      x23 = x22*x2
      x41 = x4
      x42 = x41*x4
      x43 = x42*x4
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      l_sp_cden   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21    *x41    
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)    *x22            
     8  +coeff(  8)    *x23            
      l_sp_cden   =l_sp_cden   
     9  +coeff(  9)            *x42    
     1  +coeff( 10)    *x22    *x41    
     2  +coeff( 11)    *x21    *x42    
     3  +coeff( 12)    *x23    *x41    
     4  +coeff( 13)            *x41*x51
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11*x22            
     7  +coeff( 16)    *x22    *x42    
     8  +coeff( 17)*x11        *x43    
      l_sp_cden   =l_sp_cden   
     9  +coeff( 18)    *x23    *x42    
     1  +coeff( 19)*x11            *x51
     2  +coeff( 20)            *x42*x51
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)    *x22    *x43    
     5  +coeff( 23)    *x22    *x42*x51
     6  +coeff( 24)*x11*x22    *x42    
     7  +coeff( 25)                *x52
     8  +coeff( 26)*x11        *x41    
      l_sp_cden   =l_sp_cden   
     9  +coeff( 27)            *x43    
     1  +coeff( 28)    *x22        *x51
     2  +coeff( 29)                *x53
     3  +coeff( 30)*x11*x21    *x41    
     4  +coeff( 31)*x11        *x42    
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)*x11*x21    *x42    
     7  +coeff( 34)    *x22        *x53
     8  +coeff( 35)    *x23    *x43    
      l_sp_cden   =l_sp_cden   
     9  +coeff( 36)                *x51
     1  +coeff( 37)*x12*x21            
     2  +coeff( 38)    *x21    *x43    
     3  +coeff( 39)*x11*x23            
     4  +coeff( 40)*x11*x23    *x41    
     5  +coeff( 41)*x11*x22        *x52
     6  +coeff( 42)*x11*x21    *x41*x52
     7  +coeff( 43)*x12*x22    *x41    
     8  +coeff( 44)*x11*x22    *x43    
      l_sp_cden   =l_sp_cden   
     9  +coeff( 45)*x11*x22    *x41*x52
     1  +coeff( 46)*x12*x22    *x41*x51
     2  +coeff( 47)    *x22    *x43*x52
     3  +coeff( 48)*x11*x22    *x42*x52
     4  +coeff( 49)*x12*x22    *x41*x52
c
      return
      end
      function x_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 71)
      data ncoeff/ 70/
      data avdat/ -0.2241619E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.59165666E-02, 0.36186656E+00, 0.51541597E-03, 0.13882004E+00,
     + -0.18025248E+00, 0.36847979E-01, 0.25620686E-02, 0.24792604E-01,
     + -0.59285038E-02,-0.11920809E-02, 0.76509028E-03, 0.14334091E-02,
     + -0.35492757E-02,-0.25669588E-01, 0.53009973E-02,-0.34840652E-02,
     + -0.14183666E-02, 0.32226151E-03, 0.80357576E-02,-0.21317350E-02,
     + -0.29964070E-03, 0.85403473E-03, 0.42045262E-03, 0.46779758E-04,
     +  0.17531583E-03,-0.17205998E-02, 0.10508833E-02, 0.10468675E-02,
     +  0.28542315E-02,-0.70028019E-03, 0.12413575E-02, 0.10616746E-02,
     + -0.31207174E-02,-0.51434673E-02, 0.46113270E-03, 0.16802028E-02,
     +  0.15635346E-03, 0.10739246E-02, 0.13902239E-02,-0.21976498E-02,
     + -0.14775348E-02, 0.44633038E-02,-0.71728085E-02,-0.18924864E-02,
     + -0.30718838E-04, 0.11620937E-03,-0.21301852E-03,-0.25481361E-03,
     + -0.91253396E-03, 0.19472007E-04,-0.58969985E-04,-0.36379899E-03,
     +  0.12916960E-02,-0.38503672E-02,-0.97965682E-03, 0.19658513E-03,
     + -0.26469672E-04, 0.77901676E-03, 0.21896597E-02, 0.11597572E-03,
     +  0.15468587E-02,-0.25371037E-03,-0.32796459E-04, 0.24183889E-03,
     + -0.19715170E-02, 0.49740169E-03,-0.31929120E-03, 0.16767005E-02,
     +  0.53677353E-03,-0.16956910E-02,
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
     1  +coeff( 10)*x11    *x31*x41    
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x22        *x51
     4  +coeff( 13)        *x31*x41    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)*x11            *x51
     7  +coeff( 16)    *x24            
     8  +coeff( 17)    *x21*x31        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 18)            *x42    
     1  +coeff( 19)    *x23            
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x21        *x52
     4  +coeff( 22)            *x43    
     5  +coeff( 23)                *x53
     6  +coeff( 24)        *x31*x43*x51
     7  +coeff( 25)        *x31        
     8  +coeff( 26)    *x23*x31        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 27)    *x23    *x41    
     1  +coeff( 28)    *x21*x32*x41    
     2  +coeff( 29)    *x23*x32        
     3  +coeff( 30)*x11*x24            
     4  +coeff( 31)    *x22*x31*x42    
     5  +coeff( 32)    *x21*x33*x41    
     6  +coeff( 33)    *x21*x32*x42    
     7  +coeff( 34)*x11*x22*x31*x41    
     8  +coeff( 35)*x12*x23            
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 36)*x11    *x31        
     1  +coeff( 37)    *x22*x31        
     2  +coeff( 38)        *x33        
     3  +coeff( 39)        *x31*x41*x51
     4  +coeff( 40)            *x42*x51
     5  +coeff( 41)*x11        *x41    
     6  +coeff( 42)*x12                
     7  +coeff( 43)*x11*x22            
     8  +coeff( 44)        *x32*x41    
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 45)        *x31*x43    
     1  +coeff( 46)        *x32    *x52
     2  +coeff( 47)*x11*x21        *x52
     3  +coeff( 48)    *x23        *x52
     4  +coeff( 49)    *x22*x33        
     5  +coeff( 50)*x11    *x31*x42    
     6  +coeff( 51)*x11            *x53
     7  +coeff( 52)*x11*x23        *x51
     8  +coeff( 53)*x11*x24*x31        
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 54)*x11*x24*x32        
     1  +coeff( 55)    *x21*x32        
     2  +coeff( 56)        *x31    *x52
     3  +coeff( 57)*x11*x21*x31        
     4  +coeff( 58)*x11*x21        *x51
     5  +coeff( 59)*x11*x23            
     6  +coeff( 60)    *x21*x33        
     7  +coeff( 61)*x11*x22    *x41    
     8  +coeff( 62)        *x31*x41*x52
      x_sp_cq3e   =x_sp_cq3e   
     9  +coeff( 63)*x11*x21*x32        
     1  +coeff( 64)    *x21*x33*x41*x51
     2  +coeff( 65)*x12*x21*x33*x41    
     3  +coeff( 66)*x12*x21*x33    *x51
     4  +coeff( 67)*x12*x24*x31*x41    
     5  +coeff( 68)*x11*x22*x33*x42    
     6  +coeff( 69)*x12*x23*x33        
     7  +coeff( 70)    *x23*x33*x43    
c
      return
      end
      function t_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 28)
      data ncoeff/ 27/
      data avdat/ -0.1090703E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.12818176E-02,-0.10903450E+00,-0.15420499E-02, 0.23033269E-01,
     +  0.39464697E-01,-0.95034670E-02,-0.19375517E-03,-0.14187286E-02,
     +  0.13442170E-02, 0.37007514E-03, 0.58256965E-02,-0.38293833E-02,
     +  0.35635308E-04,-0.22560731E-03, 0.38742027E-03, 0.12605260E-03,
     + -0.88532304E-03,-0.43064470E-05, 0.29325817E-03,-0.19262353E-03,
     +  0.89906232E-03,-0.16224127E-03, 0.15138628E-03, 0.10093275E-03,
     + -0.38181484E-03,-0.26158008E-03, 0.74840739E-03,
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
     9  +coeff(  9)        *x31        
     1  +coeff( 10)    *x21*x31        
     2  +coeff( 11)*x11*x21            
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x21*x32        
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x21        *x52
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)*x11        *x41    
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 18)            *x43    
     1  +coeff( 19)        *x32    *x51
     2  +coeff( 20)    *x22*x31        
     3  +coeff( 21)    *x21*x31*x41    
     4  +coeff( 22)            *x42*x51
     5  +coeff( 23)            *x41*x52
     6  +coeff( 24)                *x53
     7  +coeff( 25)    *x21*x31*x41*x51
     8  +coeff( 26)*x11*x22        *x51
      t_sp_cq3e   =t_sp_cq3e   
     9  +coeff( 27)    *x23*x32        
c
      return
      end
      function y_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 26)
      data ncoeff/ 25/
      data avdat/ -0.1132265E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.17355980E-02,-0.98027058E-01, 0.17361742E+00,-0.88169752E-02,
     +  0.11920517E-01, 0.70459969E-01,-0.14687970E+00, 0.77528164E-01,
     +  0.49456630E-01,-0.52189147E-02, 0.29237957E-02,-0.43061092E-01,
     + -0.30976532E-01,-0.34608108E-02, 0.25087180E-01,-0.37413093E-02,
     + -0.15292321E-02,-0.15093255E-01, 0.10244740E-03, 0.38130695E-02,
     + -0.48909837E-03,-0.43874360E-02,-0.12374750E-01,-0.14341141E-02,
     +  0.29510392E-02,
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
     1  +coeff( 10)*x11                
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)        *x32*x41    
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)        *x31    *x52
     8  +coeff( 17)    *x21*x31    *x51
      y_sp_cq3e   =y_sp_cq3e   
     9  +coeff( 18)    *x22*x31        
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)    *x21    *x41*x51
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)    *x22*x31    *x51
     7  +coeff( 25)    *x23*x31        
c
      return
      end
      function p_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.6651885E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.26953581E-03,-0.22440212E-01, 0.26648762E-01,-0.14808230E-02,
     + -0.88825921E-03, 0.16834768E-02,-0.99387262E-02, 0.85468916E-02,
     + -0.46559642E-02,-0.17258378E-01, 0.88811675E-02, 0.69884383E-02,
     +  0.49000103E-02,-0.18218661E-02, 0.47044893E-03,-0.33084376E-04,
     +  0.23120434E-02, 0.56796372E-02,-0.94928802E-03, 0.31837393E-02,
     +  0.46418165E-03, 0.38882932E-04, 0.15608731E-02,-0.34884419E-03,
     +  0.10056791E-04,-0.28413020E-02,-0.68067359E-02,-0.26601265E-03,
     +  0.29680051E-02, 0.32135120E-03,-0.56176406E-03, 0.10254564E-02,
     + -0.16388871E-02,-0.45142617E-03,-0.82533013E-04,-0.47383143E-03,
     +  0.29534218E-03,-0.38257835E-03, 0.38979203E-03,-0.90778602E-03,
     +  0.10731557E-01,-0.14182670E-01,-0.51726418E-03, 0.25626406E-03,
     + -0.12820818E-04,-0.35293611E-04,-0.10676390E-02, 0.12499826E-03,
     + -0.21089235E-03, 0.65998291E-04,
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
      p_sp_cq3e   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)        *x31        
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)        *x32        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff(  9)    *x21    *x41    
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)            *x42    
     3  +coeff( 12)        *x31    *x51
     4  +coeff( 13)*x11    *x31        
     5  +coeff( 14)    *x22*x31        
     6  +coeff( 15)    *x21*x31    *x51
     7  +coeff( 16)        *x33*x42    
     8  +coeff( 17)    *x21            
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 18)            *x43    
     1  +coeff( 19)    *x21    *x41*x51
     2  +coeff( 20)            *x41*x51
     3  +coeff( 21)*x11            *x51
     4  +coeff( 22)        *x31    *x52
     5  +coeff( 23)    *x23*x31        
     6  +coeff( 24)    *x22*x31    *x51
     7  +coeff( 25)*x11*x23        *x51
     8  +coeff( 26)        *x33        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 27)    *x22    *x41    
     1  +coeff( 28)*x11*x21*x31        
     2  +coeff( 29)*x11*x21    *x41    
     3  +coeff( 30)    *x23*x31    *x51
     4  +coeff( 31)    *x21        *x51
     5  +coeff( 32)    *x21*x32        
     6  +coeff( 33)    *x21    *x42    
     7  +coeff( 34)            *x41*x52
     8  +coeff( 35)    *x21*x33        
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)    *x22*x31    *x52
     2  +coeff( 38)*x11*x23    *x41*x51
     3  +coeff( 39)                *x52
     4  +coeff( 40)*x11*x21            
     5  +coeff( 41)        *x32*x41    
     6  +coeff( 42)        *x31*x42    
     7  +coeff( 43)*x11    *x31    *x51
     8  +coeff( 44)    *x21*x32*x41    
      p_sp_cq3e   =p_sp_cq3e   
     9  +coeff( 45)            *x42*x52
     1  +coeff( 46)        *x31    *x53
     2  +coeff( 47)*x11*x22*x31        
     3  +coeff( 48)    *x21*x31    *x53
     4  +coeff( 49)*x11*x23*x31        
     5  +coeff( 50)*x11*x22        *x52
c
      return
      end
      function l_sp_cq3e   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.6170426E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.44671833E-02,-0.46846446E+00, 0.36266856E-02,-0.32566607E-01,
     +  0.19294158E+00,-0.46390049E-01, 0.31190638E-02, 0.35115407E-03,
     +  0.26243638E-01, 0.94423806E-02, 0.66301011E-03,-0.37512914E-02,
     + -0.19458598E-01,-0.11064994E-03, 0.13578399E-01,-0.15240341E-02,
     + -0.21853028E-02,-0.29957248E-02,-0.16791675E-01, 0.65076660E-03,
     + -0.13324132E-03, 0.56434079E-03, 0.25855487E-02,
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
     9  +coeff(  9)*x11*x21            
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31        
     3  +coeff( 12)        *x31*x41    
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x22    *x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)        *x31        
      l_sp_cq3e   =l_sp_cq3e   
     9  +coeff( 18)    *x21    *x41    
     1  +coeff( 19)    *x21*x31*x41    
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)                *x52
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)    *x23*x31*x41    
c
      return
      end
      function x_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 69)
      data ncoeff/ 68/
      data avdat/  0.1004723E-01/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.17325476E-01, 0.12895938E-02, 0.32400450E+00,-0.18144159E+00,
     +  0.33491679E-01,-0.20915778E-01,-0.45777247E-02, 0.24164346E+00,
     +  0.66552553E-02,-0.70896940E-02, 0.73459772E-02, 0.25260563E-02,
     + -0.64136373E-03, 0.18666065E-02, 0.20305028E-01,-0.19031361E-01,
     + -0.27597544E-03,-0.27687874E-01, 0.28210136E-03, 0.90174488E-03,
     + -0.19243590E-02, 0.82352096E-02, 0.29075434E-02,-0.19101713E-03,
     +  0.73438715E-02,-0.31526482E-02,-0.10779889E-02,-0.26672124E-02,
     + -0.18751795E-02,-0.76238031E-03, 0.41819684E-01, 0.84307184E-03,
     + -0.87281357E-03,-0.40107057E-02,-0.14988729E-02,-0.95365470E-03,
     + -0.10798717E-03,-0.52996917E-03, 0.19633416E-03,-0.22439932E-03,
     + -0.34372141E-02,-0.39703636E-02, 0.25658414E-04, 0.23850610E-02,
     +  0.61739893E-02, 0.35610567E-02,-0.68852762E-02,-0.97159779E-03,
     + -0.43459781E-02, 0.55512384E-03, 0.34837960E-02, 0.77068963E-03,
     +  0.39089904E-02, 0.44967915E-03, 0.43747965E-02,-0.38933987E-02,
     +  0.32770657E-03, 0.43422947E-03, 0.15540083E-02,-0.10224914E-02,
     +  0.19815916E-02,-0.24504107E-02,-0.16981479E-01,-0.17018210E-01,
     +  0.11298588E-01,-0.98773045E-03, 0.14299456E-02,-0.10480640E-02,
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
     2  +coeff(  2)            *x41    
     3  +coeff(  3)                *x51
     4  +coeff(  4)*x11                
     5  +coeff(  5)    *x21        *x51
     6  +coeff(  6)                *x52
     7  +coeff(  7)*x11*x21            
     8  +coeff(  8)    *x21            
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)    *x21*x31*x41    
     3  +coeff( 12)                *x53
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x22        *x51
     6  +coeff( 15)    *x22            
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)    *x21        *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 18)    *x24            
     1  +coeff( 19)    *x22*x32        
     2  +coeff( 20)    *x21*x31*x41*x51
     3  +coeff( 21)    *x24        *x51
     4  +coeff( 22)            *x42    
     5  +coeff( 23)*x11    *x31*x41    
     6  +coeff( 24)    *x21        *x53
     7  +coeff( 25)    *x22*x32    *x51
     8  +coeff( 26)    *x22*x31        
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 27)*x11            *x51
     1  +coeff( 28)    *x23        *x51
     2  +coeff( 29)*x11    *x32        
     3  +coeff( 30)    *x22        *x52
     4  +coeff( 31)*x11*x23            
     5  +coeff( 32)            *x42*x52
     6  +coeff( 33)    *x23*x32        
     7  +coeff( 34)*x11*x21*x31*x41    
     8  +coeff( 35)    *x23        *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 36)*x11*x24            
     1  +coeff( 37)*x11        *x42*x51
     2  +coeff( 38)    *x21*x33*x41    
     3  +coeff( 39)*x11*x22        *x52
     4  +coeff( 40)    *x23*x33        
     5  +coeff( 41)    *x22*x32    *x52
     6  +coeff( 42)*x12*x21*x31*x41    
     7  +coeff( 43)    *x21*x33*x42    
     8  +coeff( 44)    *x21    *x41    
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 45)    *x23            
     1  +coeff( 46)    *x22    *x41    
     2  +coeff( 47)    *x21    *x42    
     3  +coeff( 48)    *x21*x31    *x51
     4  +coeff( 49)*x11*x22            
     5  +coeff( 50)    *x23*x31        
     6  +coeff( 51)    *x21*x33        
     7  +coeff( 52)*x11*x21*x32        
     8  +coeff( 53)*x11*x21*x32    *x52
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 54)*x11    *x31*x43*x53
     1  +coeff( 55)    *x22    *x42    
     2  +coeff( 56)    *x21    *x43    
     3  +coeff( 57)    *x21*x32    *x51
     4  +coeff( 58)*x12            *x51
     5  +coeff( 59)*x11*x22        *x51
     6  +coeff( 60)        *x31*x41*x52
     7  +coeff( 61)        *x31    *x53
     8  +coeff( 62)            *x41*x53
      x_sp_cq3x   =x_sp_cq3x   
     9  +coeff( 63)*x12*x22            
     1  +coeff( 64)    *x22*x31*x41*x51
     2  +coeff( 65)    *x22    *x42*x51
     3  +coeff( 66)    *x22        *x53
     4  +coeff( 67)*x12*x21*x31        
     5  +coeff( 68)    *x21*x32*x42    
c
      return
      end
      function t_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.4031088E-02/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50792629E-02, 0.22505622E-01,-0.54356433E-02, 0.11892874E+00,
     + -0.39894432E-01, 0.76388191E-02,-0.10637320E-01, 0.52072448E-02,
     +  0.21180375E-02, 0.21468622E-03, 0.71490282E-03,-0.95874054E-03,
     +  0.31165313E-02, 0.11342067E-02,-0.31627899E-02,-0.15046949E-02,
     + -0.28736377E-02,-0.36745775E-03, 0.99775509E-03,-0.51065831E-03,
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
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)        *x31        
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)    *x21*x31*x41    
      t_sp_cq3x   =t_sp_cq3x   
     9  +coeff( 18)*x12    *x31        
     1  +coeff( 19)    *x22*x32    *x51
     2  +coeff( 20)*x11*x23            
c
      return
      end
      function y_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 32)
      data ncoeff/ 31/
      data avdat/  0.2089276E-02/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.10050718E-04,-0.77753574E-01, 0.10877363E+00,-0.56203278E-02,
     +  0.80173966E-02, 0.16621327E-01,-0.20623675E-01, 0.31381799E-02,
     +  0.13834178E-01, 0.23651076E-01,-0.30886610E-02, 0.10513347E-03,
     + -0.31560570E-01,-0.19522395E-01,-0.11313975E-01, 0.24547914E-01,
     + -0.15470510E-01,-0.14580353E-03, 0.16562333E-01, 0.19923689E-01,
     + -0.23424413E-01,-0.16525558E-02,-0.36922350E-03,-0.15118278E-01,
     + -0.30271249E-03,-0.31029587E-02, 0.38013014E-02,-0.11153836E-01,
     + -0.83989440E-03, 0.52528610E-02,-0.33383719E-02,
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
     3  +coeff( 12)                *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)    *x21    *x41    
     6  +coeff( 15)        *x33        
     7  +coeff( 16)        *x32*x41    
     8  +coeff( 17)        *x31*x42    
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff( 18)        *x32    *x51
     1  +coeff( 19)*x11    *x31        
     2  +coeff( 20)    *x21*x31*x41    
     3  +coeff( 21)    *x21    *x42    
     4  +coeff( 22)    *x21*x31    *x51
     5  +coeff( 23)    *x21    *x41*x51
     6  +coeff( 24)    *x22*x31        
     7  +coeff( 25)*x12                
     8  +coeff( 26)    *x22*x31    *x51
      y_sp_cq3x   =y_sp_cq3x   
     9  +coeff( 27)    *x22            
     1  +coeff( 28)    *x22    *x41    
     2  +coeff( 29)        *x31    *x53
     3  +coeff( 30)*x11*x21*x31        
     4  +coeff( 31)    *x22*x33*x42    
c
      return
      end
      function p_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.1815546E-03/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.72440592E-03,-0.38713661E-02, 0.35281129E-01,-0.66653885E-01,
     +  0.13765773E-02,-0.20956800E-02, 0.75985445E-02, 0.22751724E-01,
     +  0.75907018E-02,-0.20833118E-02,-0.12816302E-01, 0.57024159E-02,
     + -0.22296410E-02, 0.15598638E-04, 0.17280400E-02,-0.55778073E-02,
     +  0.68113493E-03, 0.23727313E-03,-0.95726887E-03,-0.10918212E-02,
     +  0.28311133E-02,-0.36723868E-02, 0.62468732E-02,-0.11042327E-02,
     + -0.11793843E-02,-0.45280319E-03,-0.61263965E-03,-0.91316318E-02,
     + -0.50133828E-03, 0.10169160E-02,-0.21449269E-02,-0.29044873E-02,
     +  0.15249336E-02, 0.73776866E-03, 0.28893125E-03,-0.51857537E-03,
     + -0.23217402E-03, 0.79962402E-03, 0.22487265E-03,-0.25245568E-03,
     +  0.11373776E-02,-0.11564692E-03, 0.12493231E-02, 0.32227065E-02,
     + -0.81569917E-03,-0.71803411E-03,-0.21521737E-03,-0.20133627E-03,
     +  0.16585190E-02, 0.19594977E-03,
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
      p_sp_cq3x   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)*x11    *x31        
     5  +coeff( 14)        *x33        
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)        *x32        
     8  +coeff( 17)*x11*x21            
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x21*x31    *x51
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)                *x51
     4  +coeff( 22)            *x41*x51
     5  +coeff( 23)    *x22    *x41    
     6  +coeff( 24)    *x23*x31        
     7  +coeff( 25)            *x42*x52
     8  +coeff( 26)        *x31    *x53
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 27)                *x52
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)    *x21*x32        
     3  +coeff( 30)    *x21    *x42    
     4  +coeff( 31)            *x42*x51
     5  +coeff( 32)*x11*x21*x31        
     6  +coeff( 33)*x11    *x31    *x51
     7  +coeff( 34)*x11            *x52
     8  +coeff( 35)    *x21*x33        
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 36)    *x22    *x42    
     1  +coeff( 37)        *x33    *x51
     2  +coeff( 38)        *x32    *x52
     3  +coeff( 39)*x11*x21*x32        
     4  +coeff( 40)    *x22*x33        
     5  +coeff( 41)*x11*x23*x31        
     6  +coeff( 42)*x12*x23        *x51
     7  +coeff( 43)        *x31*x42    
     8  +coeff( 44)        *x31*x41*x51
      p_sp_cq3x   =p_sp_cq3x   
     9  +coeff( 45)    *x21        *x52
     1  +coeff( 46)    *x21*x32*x41    
     2  +coeff( 47)    *x21*x31    *x53
     3  +coeff( 48)            *x42*x53
     4  +coeff( 49)*x11*x21*x32*x41    
     5  +coeff( 50)*x11    *x31*x42*x51
c
      return
      end
      function l_sp_cq3x   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.4815066E-02/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.88770967E-02,-0.46839881E+00, 0.37983279E-02,-0.31425275E-01,
     +  0.19261482E+00,-0.50899643E-01, 0.59119305E-02,-0.88837929E-02,
     +  0.50423588E-02,-0.13604821E-01,-0.20431415E-02, 0.22769134E-05,
     +  0.48609413E-02, 0.27630433E-01, 0.15868621E-02,-0.25002989E-02,
     + -0.14811853E-01, 0.82494393E-02,-0.10167813E-02,-0.20099273E-02,
     + -0.13363789E-02, 0.35734361E-03, 0.23434770E-02,
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)                *x52
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff(  9)        *x32        
     1  +coeff( 10)        *x31*x41    
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)            *x42    
     5  +coeff( 14)*x11*x21            
     6  +coeff( 15)                *x53
     7  +coeff( 16)        *x31        
     8  +coeff( 17)*x11            *x51
      l_sp_cq3x   =l_sp_cq3x   
     9  +coeff( 18)    *x21*x31*x41    
     1  +coeff( 19)    *x23        *x51
     2  +coeff( 20)    *x21    *x41    
     3  +coeff( 21)    *x22        *x51
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)    *x23*x32        
c
      return
      end
      function x_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 91)
      data ncoeff/ 90/
      data avdat/  0.8992787E-01/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.35103913E-01, 0.31482238E+00,-0.15513309E-01, 0.70119715E+00,
     + -0.30714041E+00, 0.21207345E+00,-0.60078766E-01, 0.18903148E+00,
     +  0.16150771E-01, 0.43576296E-01, 0.54704789E-01,-0.18929555E+00,
     +  0.24050359E-01,-0.11376468E-02,-0.34824159E-01,-0.25288551E-02,
     +  0.90768166E-01,-0.77480726E-01, 0.10684884E+00, 0.59912126E-02,
     + -0.29735330E-02,-0.31356979E-01, 0.11192175E-01,-0.78434292E-02,
     +  0.55915225E-01,-0.37969463E-01, 0.67333598E-02, 0.28509375E+00,
     +  0.39217848E-01,-0.42648244E-04, 0.30003944E-02, 0.21208540E-03,
     + -0.65282714E-02, 0.16145296E-01,-0.25470888E-01,-0.14891231E-01,
     +  0.16824562E-01,-0.15529895E-01, 0.23223059E-01,-0.23450976E-01,
     + -0.50779231E-01, 0.14917074E-02,-0.12157117E-02, 0.91830352E-02,
     + -0.24801774E+00, 0.13094683E+00,-0.41420568E-01,-0.62533718E+00,
     +  0.34046730E+00, 0.15275083E-01, 0.12690094E-01,-0.46379301E-02,
     + -0.29222963E-01, 0.31248093E-01,-0.28137278E-02, 0.77322707E-01,
     + -0.17690538E+00,-0.88149361E-01,-0.59693335E-02, 0.49606871E-01,
     +  0.75471816E-02, 0.10127027E+00,-0.67528347E-02, 0.11499772E-01,
     + -0.70032380E-02, 0.23678573E-01,-0.59035039E-02, 0.55092532E-01,
     + -0.10228904E-01,-0.18661957E-01,-0.98327547E-02, 0.95315479E-01,
     + -0.79670614E-02, 0.27020862E-02,-0.19557219E-01, 0.31539134E-02,
     + -0.10051022E-02,-0.21078477E-02, 0.36266502E-01,-0.29597737E-01,
     + -0.48602582E-02,-0.76891389E-02, 0.53223048E-03, 0.46938166E-01,
     + -0.11072947E-01, 0.24383161E-01, 0.13330251E-01,-0.15758839E-01,
     +  0.95054144E-02,-0.48588201E-01,
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
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)            *x41*x51
     8  +coeff(  8)                *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff(  9)        *x31        
     1  +coeff( 10)                *x53
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)*x11            *x51
     4  +coeff( 13)    *x22            
     5  +coeff( 14)        *x32        
     6  +coeff( 15)    *x21*x32        
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)    *x21        *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 18)*x11            *x52
     1  +coeff( 19)    *x21*x32    *x51
     2  +coeff( 20)    *x21    *x41    
     3  +coeff( 21)        *x32    *x51
     4  +coeff( 22)    *x22        *x52
     5  +coeff( 23)    *x21        *x53
     6  +coeff( 24)    *x24        *x52
     7  +coeff( 25)    *x21*x31*x41    
     8  +coeff( 26)    *x21    *x42    
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 27)*x12        *x41    
     1  +coeff( 28)        *x32    *x52
     2  +coeff( 29)*x11*x21        *x52
     3  +coeff( 30)    *x22*x32    *x51
     4  +coeff( 31)*x12    *x32        
     5  +coeff( 32)*x11*x21*x32    *x51
     6  +coeff( 33)    *x23*x31*x43*x53
     7  +coeff( 34)    *x22    *x41    
     8  +coeff( 35)    *x24            
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 36)        *x31    *x52
     1  +coeff( 37)            *x41*x52
     2  +coeff( 38)*x11*x21*x31        
     3  +coeff( 39)    *x23*x31        
     4  +coeff( 40)    *x23    *x41    
     5  +coeff( 41)*x11*x21        *x51
     6  +coeff( 42)*x11    *x32        
     7  +coeff( 43)    *x22*x31    *x51
     8  +coeff( 44)*x11*x23            
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 45)    *x21*x31*x41*x51
     1  +coeff( 46)    *x21    *x42*x51
     2  +coeff( 47)    *x24        *x51
     3  +coeff( 48)        *x31*x41*x52
     4  +coeff( 49)            *x42*x52
     5  +coeff( 50)            *x41*x53
     6  +coeff( 51)*x12*x22            
     7  +coeff( 52)*x11*x24            
     8  +coeff( 53)    *x22*x31    *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 54)    *x22    *x41*x52
     1  +coeff( 55)    *x21*x32    *x52
     2  +coeff( 56)        *x32    *x53
     3  +coeff( 57)        *x31*x41*x53
     4  +coeff( 58)*x11*x21        *x53
     5  +coeff( 59)    *x23        *x53
     6  +coeff( 60)    *x22    *x42*x52
     7  +coeff( 61)    *x22    *x41*x53
     8  +coeff( 62)*x12            *x53
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 63)    *x24        *x53
     1  +coeff( 64)    *x22    *x42*x53
     2  +coeff( 65)        *x31*x41    
     3  +coeff( 66)    *x23            
     4  +coeff( 67)    *x22*x31        
     5  +coeff( 68)    *x22        *x51
     6  +coeff( 69)        *x31    *x53
     7  +coeff( 70)    *x23        *x52
     8  +coeff( 71)*x11            *x53
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 72)            *x42*x53
     1  +coeff( 73)    *x22*x31    *x53
     2  +coeff( 74)    *x21*x32    *x53
     3  +coeff( 75)*x11*x22            
     4  +coeff( 76)*x11    *x31*x41    
     5  +coeff( 77)    *x23*x32        
     6  +coeff( 78)    *x21*x33*x41    
     7  +coeff( 79)*x11*x23        *x51
     8  +coeff( 80)*x11*x22        *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 81)    *x23*x33        
     1  +coeff( 82)    *x23*x32    *x51
     2  +coeff( 83)    *x21*x33*x42    
     3  +coeff( 84)*x12*x21        *x52
     4  +coeff( 85)*x11*x22*x31*x42    
     5  +coeff( 86)*x11*x22    *x43    
     6  +coeff( 87)        *x31*x42*x53
     7  +coeff( 88)            *x43*x53
     8  +coeff( 89)*x11*x21*x32    *x52
      x_sp_cfp    =x_sp_cfp    
     9  +coeff( 90)*x11*x21*x31*x41*x52
c
      return
      end
      function t_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 21)
      data ncoeff/ 20/
      data avdat/  0.4031016E-02/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.50795036E-02, 0.22504391E-01,-0.54329629E-02, 0.11892929E+00,
     + -0.39893061E-01, 0.76381871E-02,-0.10637998E-01, 0.52046305E-02,
     +  0.21185761E-02, 0.21616720E-03, 0.71355357E-03,-0.95982692E-03,
     +  0.31155243E-02, 0.11309916E-02,-0.31614560E-02,-0.15067544E-02,
     + -0.28716899E-02,-0.36709348E-03, 0.99568965E-03,-0.51103521E-03,
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
      x51 = x5
      x52 = x51*x5
      x53 = x52*x5
c
c                  function
c
      t_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21        *x51
     7  +coeff(  7)                *x52
     8  +coeff(  8)        *x31        
      t_sp_cfp    =t_sp_cfp    
     9  +coeff(  9)    *x22            
     1  +coeff( 10)        *x32        
     2  +coeff( 11)    *x21*x32        
     3  +coeff( 12)    *x21        *x52
     4  +coeff( 13)    *x21*x31        
     5  +coeff( 14)                *x53
     6  +coeff( 15)    *x21    *x41    
     7  +coeff( 16)        *x31*x41    
     8  +coeff( 17)    *x21*x31*x41    
      t_sp_cfp    =t_sp_cfp    
     9  +coeff( 18)*x12    *x31        
     1  +coeff( 19)    *x22*x32    *x51
     2  +coeff( 20)*x11*x23            
c
      return
      end
      function y_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.3275587E-02/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73451074E-02, 0.12852822E+00,-0.20273080E+00, 0.24079900E-01,
     + -0.58775726E-02,-0.70780760E+00, 0.14394141E+01,-0.72541785E+00,
     + -0.19730256E+00, 0.12591812E+00,-0.50010980E-03, 0.55452090E-01,
     + -0.49973298E-01,-0.13283254E+02, 0.61112614E+01,-0.13357754E+01,
     +  0.39830174E+01,-0.28075316E+01,-0.13390298E-01,-0.25692382E+00,
     +  0.33987597E+00,-0.99448577E-01, 0.11187131E+00, 0.16377339E-01,
     +  0.16804362E-02, 0.91584958E-02, 0.79170624E-02,-0.19849132E+01,
     +  0.92118797E+01,-0.14318369E-01, 0.15075479E-01,-0.89380750E-02,
     +  0.17247731E-01,-0.16768189E-01, 0.96540954E-02,-0.14286755E-01,
     + -0.34739566E-02,-0.18270519E-03, 0.37682645E-02,-0.15772936E-02,
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
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)            *x43    
     7  +coeff( 16)        *x32    *x51
     8  +coeff( 17)        *x31*x41*x51
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 18)            *x42*x51
     1  +coeff( 19)    *x22            
     2  +coeff( 20)        *x31    *x52
     3  +coeff( 21)            *x41*x52
     4  +coeff( 22)    *x21*x31*x41    
     5  +coeff( 23)    *x21    *x42    
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)        *x31    *x53
     8  +coeff( 26)    *x21*x31    *x52
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 27)    *x21        *x51
     1  +coeff( 28)        *x33        
     2  +coeff( 29)        *x32*x41    
     3  +coeff( 30)            *x41*x53
     4  +coeff( 31)    *x22*x31    *x51
     5  +coeff( 32)*x12            *x51
     6  +coeff( 33)*x11*x23            
     7  +coeff( 34)                *x52
     8  +coeff( 35)*x11    *x31        
      y_sp_cfp    =y_sp_cfp    
     9  +coeff( 36)                *x53
     1  +coeff( 37)    *x21    *x41*x51
     2  +coeff( 38)        *x34        
     3  +coeff( 39)        *x33*x41    
     4  +coeff( 40)    *x21        *x52
c
      return
      end
      function p_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 49)
      data ncoeff/ 48/
      data avdat/  0.1814489E-03/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.71535131E-03,-0.39349985E-02, 0.35330325E-01,-0.66698298E-01,
     +  0.14464370E-02,-0.17392813E-02, 0.79319850E-02, 0.22318365E-01,
     +  0.13102639E-01,-0.48238761E-02,-0.13662448E-01, 0.22262069E-02,
     + -0.21560057E-02, 0.98313915E-03, 0.19823068E-02,-0.83001982E-02,
     +  0.33420403E-03, 0.10635600E-03,-0.62012591E-03, 0.45262679E-03,
     + -0.33055339E-03, 0.28338043E-02,-0.67211041E-03,-0.23040783E-02,
     + -0.28996230E-02, 0.97063435E-02,-0.36397309E-03, 0.82130995E-04,
     + -0.91272807E-02, 0.19576028E-03, 0.66682213E-03, 0.15207418E-02,
     +  0.44226105E-03,-0.14457773E-03, 0.19678276E-02, 0.11274287E-02,
     + -0.79841568E-03, 0.73903898E-03, 0.12338631E-02,-0.61776012E-03,
     + -0.35574229E-02,-0.95388445E-03,-0.14829639E-02,-0.23296096E-03,
     + -0.18229196E-03, 0.23424483E-03, 0.12791717E-03,-0.32833585E-03,
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
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21*x31        
     8  +coeff(  8)    *x21    *x41    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)    *x22*x31        
     4  +coeff( 13)*x11    *x31        
     5  +coeff( 14)        *x33        
     6  +coeff( 15)        *x31    *x52
     7  +coeff( 16)        *x32        
     8  +coeff( 17)*x11*x21            
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 18)    *x23            
     1  +coeff( 19)    *x21*x31    *x51
     2  +coeff( 20)        *x32    *x51
     3  +coeff( 21)            *x41*x53
     4  +coeff( 22)                *x51
     5  +coeff( 23)    *x23*x31        
     6  +coeff( 24)        *x31*x41*x52
     7  +coeff( 25)            *x41*x51
     8  +coeff( 26)    *x22    *x41    
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 27)            *x41*x52
     1  +coeff( 28)    *x21        *x51
     2  +coeff( 29)*x11        *x41    
     3  +coeff( 30)    *x21*x32        
     4  +coeff( 31)*x11*x21*x31        
     5  +coeff( 32)*x11    *x31    *x51
     6  +coeff( 33)    *x21*x33        
     7  +coeff( 34)        *x33    *x51
     8  +coeff( 35)        *x32    *x52
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 36)    *x22*x33        
     1  +coeff( 37)*x12*x21    *x41    
     2  +coeff( 38)*x12    *x31*x41    
     3  +coeff( 39)*x11*x23*x31        
     4  +coeff( 40)                *x52
     5  +coeff( 41)*x11*x21    *x41    
     6  +coeff( 42)    *x21*x32*x41    
     7  +coeff( 43)*x11*x21*x31*x41    
     8  +coeff( 44)    *x21*x31    *x53
      p_sp_cfp    =p_sp_cfp    
     9  +coeff( 45)            *x42*x53
     1  +coeff( 46)*x11*x21*x32    *x51
     2  +coeff( 47)*x11    *x33    *x51
     3  +coeff( 48)            *x43*x53
c
      return
      end
      function l_sp_cfp    (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 41)
      data ncoeff/ 40/
      data avdat/ -0.2701343E+00/
      data xmin/
     1 -0.11919E+00,-0.45368E-01,-0.59439E-01,-0.27342E-01,-0.42633E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.21947E-01, 0.49837E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.96587010E-01,-0.13505145E+01, 0.41452520E-01,-0.19606249E+01,
     +  0.10528346E+01,-0.93669802E-01,-0.59395683E+00, 0.28014690E+00,
     + -0.59458995E+00,-0.13443172E+00, 0.52955550E+00, 0.16942929E-01,
     + -0.24789049E+00, 0.27797652E-01,-0.42088367E-01,-0.98114563E-02,
     + -0.23643617E+00,-0.48908904E-01, 0.20155773E+00, 0.29845797E-01,
     +  0.12310249E-01, 0.89808358E-02, 0.18015454E-01, 0.66377781E-01,
     + -0.80602011E-02,-0.10768568E-01,-0.54349229E-02, 0.16708244E+00,
     + -0.38140547E-02,-0.12972830E+00,-0.46198700E-01,-0.98637104E-01,
     +  0.16400295E+00, 0.79014981E-02, 0.21100126E-02, 0.29761350E+00,
     + -0.32226339E+00, 0.10036074E-01,-0.42327628E-01, 0.32512676E-01,
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
      l_sp_cfp    =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x41*x51
      l_sp_cfp    =l_sp_cfp    
     9  +coeff(  9)                *x52
     1  +coeff( 10)                *x53
     2  +coeff( 11)*x11            *x51
     3  +coeff( 12)        *x32        
     4  +coeff( 13)        *x31    *x51
     5  +coeff( 14)    *x21*x32        
     6  +coeff( 15)        *x31        
     7  +coeff( 16)    *x21*x31        
     8  +coeff( 17)    *x21        *x52
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 18)*x11        *x42    
     1  +coeff( 19)*x11            *x52
     2  +coeff( 20)    *x21*x32    *x51
     3  +coeff( 21)    *x21        *x53
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)*x11*x23            
     6  +coeff( 24)    *x21    *x42    
     7  +coeff( 25)    *x23*x31        
     8  +coeff( 26)    *x22*x32        
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 27)    *x23        *x51
     1  +coeff( 28)    *x22        *x52
     2  +coeff( 29)            *x41*x53
     3  +coeff( 30)*x11*x21        *x52
     4  +coeff( 31)    *x21    *x41*x53
     5  +coeff( 32)*x11    *x31    *x53
     6  +coeff( 33)*x11        *x41*x53
     7  +coeff( 34)        *x33*x41*x52
     8  +coeff( 35)*x12            *x53
      l_sp_cfp    =l_sp_cfp    
     9  +coeff( 36)    *x22*x31    *x53
     1  +coeff( 37)    *x22    *x41*x53
     2  +coeff( 38)*x11*x22*x31    *x52
     3  +coeff( 39)*x11*x23        *x53
     4  +coeff( 40)*x12*x21*x32    *x52
c
      return
      end
      function x_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 61)
      data ncoeff/ 60/
      data avdat/ -0.3661939E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.73624090E-02, 0.54408592E+00, 0.12991597E-02, 0.13331600E+00,
     + -0.25333229E+00,-0.19299382E-02, 0.28006369E-01,-0.43372321E-03,
     +  0.19683193E-01,-0.30577199E-02, 0.13395125E-02,-0.89188088E-02,
     + -0.19269234E-03,-0.35047766E-02,-0.94378386E-02, 0.10281038E-01,
     +  0.13868400E-01,-0.73659522E-02, 0.24886462E-02,-0.85966202E-03,
     +  0.21292618E-02,-0.55529363E-02,-0.23617770E-02,-0.10597468E-02,
     + -0.36377823E-02, 0.12307853E-01,-0.12758012E-02, 0.61965511E-05,
     +  0.71204454E-02,-0.91621783E-02, 0.13322660E-02, 0.34553681E-02,
     + -0.28105325E-03, 0.60220109E-03, 0.39563395E-03,-0.62230336E-02,
     + -0.21036442E-02, 0.78278519E-02,-0.82336711E-02,-0.74323928E-02,
     + -0.40811789E-02, 0.81864878E-03,-0.70517824E-03,-0.17356533E-02,
     +  0.20598262E-02,-0.21846138E-02, 0.20539546E-02,-0.30908993E-03,
     + -0.68707210E-04, 0.35265679E-03,-0.41904617E-02, 0.18841500E-02,
     + -0.18986130E-03, 0.23434365E-02, 0.87841594E-03, 0.74678240E-03,
     +  0.10977738E-02,-0.19191147E-02, 0.43783667E-02,-0.28426992E-02,
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
      x_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)        *x31        
      x_sp_cdex   =x_sp_cdex   
     9  +coeff(  9)    *x22            
     1  +coeff( 10)                *x52
     2  +coeff( 11)        *x32        
     3  +coeff( 12)    *x21*x32        
     4  +coeff( 13)    *x23    *x42*x52
     5  +coeff( 14)        *x31*x41    
     6  +coeff( 15)*x11*x21            
     7  +coeff( 16)*x11            *x51
     8  +coeff( 17)    *x21*x31*x41    
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 18)    *x21    *x42    
     1  +coeff( 19)    *x21*x31        
     2  +coeff( 20)            *x42    
     3  +coeff( 21)        *x31    *x51
     4  +coeff( 22)        *x32    *x51
     5  +coeff( 23)    *x22        *x51
     6  +coeff( 24)    *x21*x31    *x51
     7  +coeff( 25)    *x24            
     8  +coeff( 26)        *x31*x41*x51
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 27)    *x23*x31        
     1  +coeff( 28)*x11            *x52
     2  +coeff( 29)    *x21*x33        
     3  +coeff( 30)    *x21*x32*x41    
     4  +coeff( 31)    *x21*x32    *x51
     5  +coeff( 32)*x11*x22    *x41    
     6  +coeff( 33)*x11*x21*x32        
     7  +coeff( 34)    *x21*x33*x41    
     8  +coeff( 35)*x11*x22*x32        
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 36)*x12*x21*x31*x41    
     1  +coeff( 37)            *x41*x51
     2  +coeff( 38)    *x23            
     3  +coeff( 39)*x11*x22            
     4  +coeff( 40)            *x42*x51
     5  +coeff( 41)    *x21*x32*x42    
     6  +coeff( 42)*x11    *x31        
     7  +coeff( 43)    *x21        *x52
     8  +coeff( 44)        *x31    *x52
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 45)*x11*x21        *x51
     1  +coeff( 46)*x11    *x32        
     2  +coeff( 47)    *x21*x31*x42    
     3  +coeff( 48)*x11*x22        *x51
     4  +coeff( 49)            *x43*x51
     5  +coeff( 50)    *x22*x32    *x51
     6  +coeff( 51)*x11*x24    *x42    
     7  +coeff( 52)            *x41*x52
     8  +coeff( 53)    *x22        *x52
      x_sp_cdex   =x_sp_cdex   
     9  +coeff( 54)*x11*x23            
     1  +coeff( 55)    *x23    *x41*x51
     2  +coeff( 56)*x11*x24    *x41    
     3  +coeff( 57)    *x21*x32*x43    
     4  +coeff( 58)*x11*x22*x33        
     5  +coeff( 59)*x11*x22*x32*x41    
     6  +coeff( 60)    *x23*x33*x41    
c
      return
      end
      function t_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/  0.5934969E+00/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.24040234E-02,-0.14557219E+00,-0.30459273E-02, 0.31407263E-01,
     +  0.52291416E-01,-0.75922109E-03,-0.13214811E-02,-0.61241048E-03,
     + -0.12179046E-02, 0.22951143E-02, 0.91115251E-03, 0.27204363E-02,
     + -0.54518715E-02,-0.38851725E-02, 0.39490527E-02, 0.32644827E-03,
     + -0.18206682E-02,-0.79409027E-03, 0.37386510E-03, 0.92563780E-04,
     +  0.17415590E-03, 0.35935882E-03, 0.75522315E-03,
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
      t_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21    *x41    
     7  +coeff(  7)    *x21        *x51
     8  +coeff(  8)            *x42    
      t_sp_cdex   =t_sp_cdex   
     9  +coeff(  9)                *x52
     1  +coeff( 10)    *x21*x32        
     2  +coeff( 11)        *x32    *x51
     3  +coeff( 12)        *x31        
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)    *x22        *x51
     8  +coeff( 17)        *x31*x41*x51
      t_sp_cdex   =t_sp_cdex   
     9  +coeff( 18)            *x43*x51
     1  +coeff( 19)    *x21        *x52
     2  +coeff( 20)*x11*x22            
     3  +coeff( 21)    *x21*x32    *x51
     4  +coeff( 22)        *x33    *x51
     5  +coeff( 23)    *x23*x31*x41    
c
      return
      end
      function y_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.1641307E-02/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.22031993E-02,-0.77386558E-01, 0.14951643E+00,-0.73454655E-02,
     +  0.10073172E-01, 0.70882425E-01,-0.16234212E+00, 0.94544090E-01,
     +  0.37688684E-01,-0.44621546E-02, 0.13002012E-02,-0.35638258E-01,
     + -0.23015786E-01,-0.32701287E-02, 0.20566547E-01,-0.32914667E-02,
     + -0.17936911E-01,-0.33407140E-03, 0.54100511E-03, 0.37431810E-02,
     + -0.24195353E-02,-0.47009822E-03,-0.38314755E-02,
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
     1  +coeff( 10)*x11                
     2  +coeff( 11)                *x52
     3  +coeff( 12)    *x21*x31        
     4  +coeff( 13)    *x21    *x41    
     5  +coeff( 14)        *x31*x42    
     6  +coeff( 15)*x11    *x31        
     7  +coeff( 16)    *x21*x31    *x51
     8  +coeff( 17)    *x22*x31        
      y_sp_cdex   =y_sp_cdex   
     9  +coeff( 18)        *x33    *x52
     1  +coeff( 19)    *x21        *x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)        *x31    *x52
     4  +coeff( 22)    *x21        *x52
     5  +coeff( 23)*x11*x23    *x41    
c
      return
      end
      function p_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 51)
      data ncoeff/ 50/
      data avdat/  0.8351380E-03/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     +  0.46990227E-03, 0.20078470E-02,-0.22617333E-01, 0.23522906E-01,
     + -0.69950084E-03,-0.40764301E-02, 0.31498375E-02,-0.10559782E-01,
     + -0.52957064E-02, 0.21972433E-02, 0.36572840E-02, 0.16220970E-02,
     + -0.11679428E-03, 0.57225344E-02,-0.57736854E-03, 0.96684624E-03,
     + -0.20320185E-03, 0.32827454E-02,-0.14811368E-02, 0.81553351E-03,
     +  0.20016274E-03,-0.33808967E-04, 0.47759102E-02, 0.79011028E-04,
     +  0.28091847E-03, 0.13495512E-04, 0.56760572E-02, 0.36725653E-02,
     +  0.19349396E-02,-0.20956795E-02,-0.13239670E-03,-0.45924040E-03,
     +  0.30306240E-02,-0.75770142E-02,-0.89327684E-02,-0.30018936E-03,
     +  0.18547758E-03, 0.84126682E-03,-0.27068792E-03,-0.10219258E-03,
     +  0.18521794E-03,-0.33981658E-02,-0.17660451E-02, 0.11710961E-02,
     + -0.19163834E-02, 0.11995891E-02,-0.25985423E-02, 0.15450551E-03,
     +  0.53237472E-04,-0.53072633E-03,
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
      p_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)        *x31        
     4  +coeff(  4)            *x41    
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x21*x31        
     7  +coeff(  7)        *x32        
     8  +coeff(  8)    *x21    *x41    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff(  9)        *x31*x41    
     1  +coeff( 10)            *x42    
     2  +coeff( 11)        *x31    *x51
     3  +coeff( 12)*x11    *x31        
     4  +coeff( 13)    *x22*x31        
     5  +coeff( 14)    *x21    *x42    
     6  +coeff( 15)        *x31*x42    
     7  +coeff( 16)    *x21*x31    *x51
     8  +coeff( 17)    *x23            
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 18)    *x21*x32        
     1  +coeff( 19)                *x51
     2  +coeff( 20)    *x22            
     3  +coeff( 21)    *x21        *x51
     4  +coeff( 22)        *x31    *x52
     5  +coeff( 23)*x11*x21    *x41    
     6  +coeff( 24)*x11            *x52
     7  +coeff( 25)    *x23*x31        
     8  +coeff( 26)*x11*x23        *x51
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 27)            *x41*x51
     1  +coeff( 28)*x11        *x41    
     2  +coeff( 29)    *x22*x31*x41    
     3  +coeff( 30)    *x22    *x42    
     4  +coeff( 31)    *x21*x31    *x52
     5  +coeff( 32)    *x22*x33        
     6  +coeff( 33)    *x23    *x42    
     7  +coeff( 34)    *x22    *x41    
     8  +coeff( 35)    *x21*x31*x41    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 36)            *x41*x52
     1  +coeff( 37)    *x21*x31*x42    
     2  +coeff( 38)        *x33    *x51
     3  +coeff( 39)    *x22    *x41*x51
     4  +coeff( 40)        *x31    *x53
     5  +coeff( 41)    *x23*x32        
     6  +coeff( 42)    *x21    *x41*x51
     7  +coeff( 43)*x11*x21*x31        
     8  +coeff( 44)    *x23    *x41    
      p_sp_cdex   =p_sp_cdex   
     9  +coeff( 45)        *x32*x41*x51
     1  +coeff( 46)        *x31*x42*x51
     2  +coeff( 47)    *x23*x31*x41    
     3  +coeff( 48)    *x22*x31    *x52
     4  +coeff( 49)    *x22        *x53
     5  +coeff( 50)*x12*x22    *x41    
c
      return
      end
      function l_sp_cdex   (x,m)
      dimension x(m)
      dimension xmin(10),xmax(10),scale(10),xmean(10)
      dimension coeff( 24)
      data ncoeff/ 23/
      data avdat/ -0.1009623E-01/
      data xmin/
     1 -0.11827E+00,-0.45368E-01,-0.59439E-01,-0.28391E-01,-0.49971E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data xmax/
     1  0.11997E+00, 0.46230E-01, 0.40983E-01, 0.22340E-01, 0.49982E-01,
     2  0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00, 0.00000E+00/
      data scale /10*0./
      data coeff/
     + -0.64331507E-02,-0.74186116E+00, 0.87732906E-02,-0.10018444E+00,
     +  0.32092661E+00,-0.51410973E-01,-0.36355702E-02,-0.12361572E-01,
     +  0.26264332E-01, 0.11578904E-02, 0.16254442E-01, 0.29233132E-01,
     + -0.23848247E-01,-0.52483547E-01, 0.37819222E-01, 0.18863921E-02,
     + -0.12081645E-03,-0.69632567E-02,-0.40301576E-01, 0.22423444E-01,
     + -0.11934834E-02, 0.19607572E-02, 0.35778824E-02,
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
      l_sp_cdex   =avdat
     1  +coeff(  1)                    
     2  +coeff(  2)    *x21            
     3  +coeff(  3)            *x41    
     4  +coeff(  4)                *x51
     5  +coeff(  5)*x11                
     6  +coeff(  6)    *x22            
     7  +coeff(  7)    *x21    *x41    
     8  +coeff(  8)    *x21        *x51
      l_sp_cdex   =l_sp_cdex   
     9  +coeff(  9)    *x21*x32        
     1  +coeff( 10)                *x52
     2  +coeff( 11)        *x32        
     3  +coeff( 12)*x11*x21            
     4  +coeff( 13)*x11            *x51
     5  +coeff( 14)    *x21*x31*x41    
     6  +coeff( 15)    *x21    *x42    
     7  +coeff( 16)            *x43    
     8  +coeff( 17)    *x22        *x51
      l_sp_cdex   =l_sp_cdex   
     9  +coeff( 18)        *x31        
     1  +coeff( 19)        *x31*x41    
     2  +coeff( 20)            *x42    
     3  +coeff( 21)        *x33        
     4  +coeff( 22)    *x22    *x42    
     5  +coeff( 23)    *x23*x31*x41    
c
      return
      end
