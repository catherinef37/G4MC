      function dfringe(x,y,z,fact,a,fx,fy,fz)
      real x,y,z,fx,fy,fz,fact,a,delta
      real c0(0:5),s,b(-2:2,-2:2)
      data c0/0.2383,1.7395,-0.4768,0.5288,-0.1299,0.0222/
      data delta/80./
c      print *, x,y,z,fact,a,fx,fy,fz
c
c   calculates  entrance fringe field of dipole magnet
c   in the same manner as raytrace, for use in snake. jjl 2/17/87
c   n.b. can only be used for bending in the x-y plane.
c
c          fact = central field
c             a = gap size (radius)
c         x,y,z = coordinates relative to magnet entrance
c                  (set entrance=0 in the relative ref frame)
c            fz = returned field value
c         c0(i) = fringe field coefficients
c
c            fz = fact/(1+exp(s))
c             s = c0(0) + c0(1)*(y/a) + c0(2)*((y/a)**2) ....etc.
c****************************************************************
c
c  set up grid of b's for expansion (see raytrace manual p. 11-12)
c
      do 20 j=-2,2
      do 20 k=-2,2
      if(abs(j)+abs(k).ge.3) go to 20
      s=c0(0)
      do 10 i=1,5
c
c  when curved efb's are introduced y will depend on x in the following
c
 10          s=s+(c0(i)*(((y+(j*delta))/a)**i))
 20             b(j,k)=fact/(1+exp(s))
c
c  calculate fields
c        1         2         3         4         5         6         7
      fz=b(0,0)
     &   -((z**2/delta**2)*
     &      (((2./3.)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(4.*b(0,0))))
     &     -((1./24.)*(b(2,0)+b(-2,0)+b(0,2)+b(0,-2)-(4.*b(0,0)))))
     &   +((z**4/delta**4)*
     &      (((-1./6.)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(4.*b(0,0))))
     &      +((1./24.)*(b(2,0)+b(-2,0)+b(0,2)+b(0,-2)-(4.*b(0,0))))
     &      +((1./12.)*(b(1,1)+b(-1,1)+b(1,-1)+b(-1,-1)))
     &      -((1./6.0)*(b(1,0)+b(-1,0)+b(0,1)+b(0,-1)-(2.*b(0,0)))))))
c        1         2         3         4         5         6         7
      fy=((z/delta)*
     &      (((2./3.)*(b(1,0)-b(-1,0)))
     &       -((1./12.)*(b(2,0)-b(-2,0)))))
     &   +(((z**3)/(delta**3))*
     &      (((1./6.)*(b(1,0)-b(-1,0)))
     &       -((1./12.)*(b(2,0)-b(-2,0)))
     &       -((1./12.)*(b(1,1)+b(1,-1)-b(-1,1)-b(-1,-1)
     &           -(2.*b(1,0))+(2.*b(-1,0))))))
c        1         2         3         4         5         6         7
      fx=((z/delta)*
     &      (((2./3.)*(b(0,1)-b(0,-1)))
     &       -((1./12.)*(b(0,2)-b(0,-2)))))
     &   +(((z**3)/(delta**3))*
     &      (((1./6.)*(b(0,1)-b(0,-1)))
     &       -((1./12.)*(b(0,2)-b(0,-2)))
     &       -((1./12.)*(b(1,1)+b(-1,1)-b(1,-1)-b(-1,-1)
     &           -(2.*b(0,1))+(2.*b(0,-1))))))
c      print *, fx, fy, fz
      return
      end
