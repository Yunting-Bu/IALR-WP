      subroutine inverse_real(nc,A)
      implicit none

      integer nc,LDA,IPIV(nc),LWORK,INFO
      real*8 work(nc)
      real*8 A(nc,nc)

      lwork=nc
      lda=nc

      call dgetrf(nc,nc,a,lda,ipiv,info)
      call dgetri(nc,a,lda,ipiv,work,lwork,info)

      return
      endsubroutine



      subroutine inverse_cmplx(nc,A)
      implicit none

      integer nc,LDA,IPIV(nc),LWORK,INFO
      complex*16 work(nc)
      complex*16 A(nc,nc)

      lwork=nc
      lda=nc

      call zgetrf(nc,nc,a,lda,ipiv,info)
      call zgetri(nc,a,lda,ipiv,work,lwork,info)

      return
      endsubroutine


        FUNCTION spgndr(l,mm,x) !normalized associated Legendre function
        implicit real*8 (a-h, o-z)

        INTEGER l,m,mm

        m=abs(mm)
        fact=(2*l+1)/2.d0
        do i=l-m+1, l+m
                fact=fact/i
        enddo
        fact=sqrt(fact)
        if(mm.lt.0.and.m/2*2.ne.m) fact=-fact
        spgndr=fact*plgndr(l,m,x)
        ENDFUNCTION



      FUNCTION plgndr(l,m,x)
      implicit real*8(a-h,o-z)
      INTEGER l,m
c      REAL plgndr,x
c      INTEGER i,ll
c      REAL fact,pll,pmm,pmmp1,somx2
      if(m.lt.0.or.m.gt.l.or.dabs(x).gt.1.d0) then
        write(*,"(2I4,E15.5)") l,m,x
        stop 'bad arguments in plgndr'
      endif
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
          plgndr=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
          plgndr=pll
        endif
      endif
      return
      ENDFUNCTION

      function spharmonic(theta,phi,l,m)
      implicit none
      complex*16 :: spharmonic
      real*8,parameter :: pi = dacos(-1.d0)
      complex*16,parameter :: ai = (0.d0,1.d0)
      integer :: l,m
      real*8 :: theta_use,phi_use
      real*8 :: theta,phi
      real*8 :: PARSGN,PLM
      real*8 :: spgndr

      theta_use = dcos(theta)
!      spharmonic = PARSGN(m)*PLM(l,m,theta_use)*
      spharmonic  = spgndr(l,m,theta_use)*
     &cdexp(ai*m*phi)/(2.d0*pi)**0.5d0

      return
      end function




      SUBROUTINE GAULEG(X1,X2,X,W,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(N),W(N)
      DATA EPS/3.0d-14/
      PI=DACOS(-1.0D0)
      M=(N+1)/2
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
      Z=DCOS(PI*(I-0.25D0)/(N+0.5D0))
1     CONTINUE
      P1=1.0D0
      P2=0.0D0
      DO 11 J=1,N
      P3=P2
      P2=P1
      P1=((2.0D0*J-1.0D0)*Z*P2-(J-1.0D0)*P3)/J
11    CONTINUE
      PP=N*(Z*P1-P2)/(Z*Z-1.0D0)
      Z1=Z
      Z=Z1-P1/PP
      IF(ABS(Z-Z1).GT.EPS) GOTO 1
      X(I)=XM-XL*Z
      X(N+1-I)=XM+XL*Z
      W(I)=2.0*XL/((1.0D0-Z*Z)*PP*PP)
      W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      ENDSUBROUTINE



c***********************************************************************
        SUBROUTINE DJMM(J,MP,M,BETA,R,ID)
        IMPLICIT REAL*8(A-H,O-Z)
        DIMENSION AF(100)

c when id.eq.1 , d matrix reduced to lengendre ploynomials(m'=0 and m=0)
c when id.eq.2 , d matrix reduced to associate ploynomials(m=0)
c when id.eq.3 , d matrix, all three cases are normalized for theta

        PI=DACOS(-1.D0)
        CALL NFACTO(50,AF,0)
        CALL DLEB(J,MP,M,BETA,R,AF)
        
        IF(ID.EQ.1) THEN
          IF(MP.NE.0.OR.M.NE.0) STOP 'BAD INPUT OF M'
          R=DSQRT((2.D0*J+1)/2.D0)*R
        END IF

        IF(ID.EQ.2) THEN
          IF(M.NE.0) STOP 'BAD INPUT OF M'
        
c there is a factor between the normal lengendre and reduced d ,
c it is  dsqrt((j+m)!/(j-m)!)

          RR=R*DSQRT(AF(J+MP+1)/AF(J-MP+1))
          FACTOR=DSQRT((2.D0*J+1)*AF(J-MP+1)/(4.D0*PI*AF(J+MP+1)))
          R=(-1.D0)**MP*FACTOR*RR
        
          R=DSQRT(2.D0*PI)*R
        END IF

        IF(ID.EQ.3) R=(2.D0*PI)*DSQRT((2.D0*J+1)/(8.D0*PI**2))*R

        RETURN
        END     
c=================================================================
         SUBROUTINE DLEB(J,MP,M,BETA,R,AF)
         IMPLICIT REAL*8(A-H,O-Z)
CHE05720
C***********************************************************************
C        THIS SUBROUTINE CALCULATES THE SMALL D COEF
C        CHE05740
C        EQ. 4.13 ON PAGE 52 OF M.E.ROSE
C        CHE05750
C***********************************************************************
      DIMENSION AF(1)
         C=DCOS(BETA*0.5D0)
         S=-DSIN(BETA*0.5D0)
      N1=J-MP
      N2=J+M
      N3=MP-M
      NC=2*J+M-MP
      NS=MP-M
      A1=DSQRT(AF(N2+1)*AF(J-M+1)*AF(J+MP+1)*AF(N1+1))
c
c          determine the summation index bound
c
       KST=MAX0(0,M-MP)+1
       KEND=MIN0(N1,N2)+1                       
   10  R=0.D0
       DO 1 KK=KST,KEND
        K=KK-1
        NC1=NC-2*K
         NS1=NS+2*K
          K=KK-1
           X=(-1)**K
          IF(NC1.NE.0) X=X*C**NC1
         IF(NS1.NE.0) X=X*S**NS1
        Y=AF(N1-K+1)*AF(N2-K+1)*AF(N3+KK)*AF(KK)
   1   R=R+X/Y
       R=A1*R
       RETURN
      END
c=====================================================================
      SUBROUTINE NFACTO(N,A,IOP)
      IMPLICIT REAL*8(A-H,O-Z)
C*************************************************************************
C  IF IOP=0: CALCULATES THE NFACTORIAL AN(I)=(N-1)!
C  IF IOP GREATER THAN 0: CALCULATES THE LOG OF FACTORIAL
C  AN(I)=DLOG((N-1)!)
C*************************************************************************
      DIMENSION A(N)
      IF(N.LT.1) STOP 'N IS LESS THAN 1 IN NFACTO'
      IF(IOP.GT.0) GO TO 3
      A(1)=1.D0
      DO 1 I=2,N
      A(I)=DFLOAT(I-1)*A(I-1)
   1  CONTINUE
      RETURN
   3  CONTINUE
      A(1)=0.D0
      DO 2 I=2,N
      A(I)=A(I-1)+DLOG(DFLOAT(I-1))
   2  CONTINUE
      RETURN
      END



        DOUBLE PRECISION FUNCTION CG(J1,M1,J2,M2,J3)
        IMPLICIT REAL*8(A-H,O-Z)
        M3=M1+M2
        M33=-M3
        CG=(-1)**(J1+J2+M33)*F3J(J1,J2,J3,M1,M2,M33,2)
     $  *DSQRT(2.D0*J3+1.D0)
        RETURN
        END

      DOUBLE PRECISION FUNCTION F3J(JD1,JD2,JD3,MD1,MD2,MD3,ITWO)
C FROM NBS TECHNICAL NOTE 409
C
C CALCULATE 3J SYMBOL. THIS FUNCTION WORKS WITH BOTH INTEGRAL AND HALF
C INTEGRAL ARGUMENTS. IF ALL INTEGRAL ARGUMENTS, USE ITWO=2 AND JD1,ETC
C EQUAL TO THE ARGUMENTS. IF SOME HALF INTEGRAL ARGUMENTS, USE ITWO=1
C AND JD1,ETC EQUAL TO TWICE THE ARGUMENTS.
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DOUBLE PRECISION FL(3220)
      INTEGER NCALL
      SAVE FL, NCALL
      DIMENSION MTRI(9)
      DATA ZERO,HALF,ONE/0.0D+0,0.5D+0,1.0D+0/
      DATA EPS,EPS2,EPS3,EPS4,EPS5,EPS6/1.0D-10,1.0D+30,1.0D-30,8.0D+1,
     $1.0D+10,23.02585092994046D+0/
      J1=JD1*ITWO
      J2=JD2*ITWO
      J3=JD3*ITWO
      M1=MD1*ITWO
      M2=MD2*ITWO
      M3=MD3*ITWO
      IF(NCALL+1867)5,15,5
5     NCALL=-1867
      FL(1)=ZERO
      FL(2)=ZERO
      DO 50 N=3,3220
      FN=N-1
50    FL(N)=FL(N-1)+DLOG(FN)
15    I=J1+J2-J3
      I1=I/2
      IF(I-2*I1)1000,1010,1000
1010  MTRI(1)=I1
      I=J1-J2+J3
      I1=I/2
      IF(I-2*I1)1000,1020,1000
1020  MTRI(2)=I1
      I=-J1+J2+J3
      I1=I/2
      IF(I-2*I1)1000,1030,1000
1030  MTRI(3)=I1
      IF(M1+M2+M3)1000,1040,1000
1040  I=J1+M1
      I1=I/2
      IF(I-2*I1)1000,1050,1000
1050  MTRI(4)=I1
      MTRI(5)=(J1-M1)/2
      I=J2+M2
      I1=I/2
      IF(I-2*I1)1000,1060,1000
1060  MTRI(6)=I1
      MTRI(7)=(J2-M2)/2
      I=J3+M3
      I1=I/2
      IF(I-2*I1)1000,1070,1000
1070  MTRI(8)=I1
      MTRI(9)=(J3-M3)/2
      DO 30 N=1,9
      IF(MTRI(N))1000,30,30
30    CONTINUE
      IF(J3-J2+M1)40,45,45
40    KMIN=-J3+J2-M1
      GOTO 60
45    KMIN=0
60    IF(-J3+J1+M2-KMIN)80,80,70
70    KMIN=-J3+J1+M2
80    KMIN=KMIN/2
      IF(J2-J3+M1)90,100,100
90    KMAX=J1+J2-J3
      GOTO 110
100   KMAX=J1-M1
110   IF(J2+M2-KMAX)120,130,130
120   KMAX=J2+M2
130   KMAX=KMAX/2
      MIN1=MTRI(1)-KMIN+1
      MIN2=MTRI(5)-KMIN+1
      MIN3=MTRI(6)-KMIN+1
      MIN4=(J3-J2+M1)/2+KMIN
      MIN5=(J3-J1-M2)/2+KMIN
      UK=EPS
      S=UK
      NCUT=0
      KMAX=KMAX-KMIN
      IF(KMAX)165,165,155
155   DO 160 K=1,KMAX
      UK=-UK*DFLOAT(MIN1-K)*DFLOAT(MIN2-K)*DFLOAT(MIN3-K)/(DFLOAT(KMIN+
     1K)*DFLOAT(MIN4+K)*DFLOAT(MIN5+K))
      IF(ABS(UK)-EPS2)158,157,157
157   UK=EPS*UK
      S=EPS*S
      NCUT=NCUT+1
158   IF(ABS(UK)-EPS3)165,160,160
160   S=S+UK
165   DELOG=ZERO
      DO 170 N=1,9
      NUM=MTRI(N)
170   DELOG=DELOG+FL(NUM+1)
      NUM=(J1+J2+J3)/2+2
      DELOG=HALF*(DELOG-FL(NUM))
      ULOG=-FL(KMIN+1)-FL(MIN1)-FL(MIN2)-FL(MIN3)-FL(MIN4+1)-FL(MIN5+1)
      PLOG=DELOG+ULOG
      IF(PLOG+EPS4)172,171,171
171   IF(NCUT)175,175,172
172   SIG=SIGN(ONE,S)
      S=ABS(S)
      SLOG=DLOG(S)+DFLOAT(NCUT+1)*EPS6
      F3J=SIG*EXP(SLOG+PLOG)
      GOTO 178
175   S=S*EPS5
      P=EXP(PLOG)
      F3J=P*S
178   NUM=KMIN+(J1-J2-M3)/2
      IF(MOD(NUM,2))180,190,180
180   F3J=-F3J
190   CONTINUE
      GOTO 2000
1000  F3J=ZERO
2000  RETURN
      END



      SUBROUTINE bessik(x,xnu,ri,rk,rip,rkp)
      INTEGER MAXIT
      REAL*8 ri,rip,rk,rkp,x,xnu,XMIN
      DOUBLE PRECISION EPS,FPMIN,PI
      PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=10000,XMIN=2.,
     *PI=3.141592653589793d0)
CU    USES beschb
      INTEGER i,l,nl
      DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,
     *gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,
     *ripl,ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      if(x.le.0.d0.or.xnu.lt.0.d0) pause 'bad arguments in bessik'
      nl=int(xnu+.5d0)
      xmu=xnu-nl
      xmu2=xmu*xmu
      xi=1.d0/x
      xi2=2.d0*xi
      h=xnu*xi
      if(h.lt.FPMIN)h=FPMIN
      b=xi2*xnu
      d=0.d0
      c=h
      do 11 i=1,MAXIT
        b=b+xi2
        d=1.d0/(b+d)
        c=b+1.d0/c
        del=c*d
        h=del*h
        if(abs(del-1.d0).lt.EPS)goto 1
11    continue
      pause 'x too large in bessik; try asymptotic expansion'
1     continue
      ril=FPMIN
      ripl=h*ril
      ril1=ril
      rip1=ripl
      fact=xnu*xi
      do 12 l=nl,1,-1
        ritemp=fact*ril+ripl
        fact=fact-xi
        ripl=fact*ritemp+ril
        ril=ritemp
12    continue
      f=ripl/ril
      if(x.lt.XMIN) then
        x2=.5d0*x
        pimu=PI*xmu
        if(abs(pimu).lt.EPS)then
          fact=1.d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS)then
          fact2=1.d0
        else
          fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        do 13 i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS)goto 2
13      continue
        pause 'bessk series failed to converge'
2       continue
        rkmu=sum
        rk1=sum1*xi2
      else
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
        do 14 i=2,MAXIT
          a=a-2*(i-1)
          c=-a*c/i
          qnew=(q1-b*q2)/a
          q1=q2
          q2=qnew
          q=q+c*qnew
          b=b+2.d0
          d=1.d0/(b+a*d)
          delh=(b*d-1.d0)*delh
          h=h+delh
          dels=q*delh
          s=s+dels
          if(abs(dels/s).lt.EPS)goto 3
14      continue
        pause 'bessik: failure to converge in cf2'
3       continue
        h=a1*h
        rkmu=sqrt(PI/(2.d0*x))*exp(-x)/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
      endif
      rkmup=xmu*xi*rkmu-rk1
      rimu=xi/(f*rkmu-rkmup)
      ri=(rimu*ril1)/ril
      rip=(rimu*rip1)/ril
      do 15 i=1,nl
        rktemp=(xmu+i)*xi2*rk1+rkmu
        rkmu=rk1
        rk1=rktemp
15    continue
      rk=rkmu
      rkp=xnu*xi*rkmu-rk1
      return
      END



      SUBROUTINE beschb(x,gam1,gam2,gampl,gammi)
      INTEGER NUSE1,NUSE2
      DOUBLE PRECISION gam1,gam2,gammi,gampl,x
      PARAMETER (NUSE1=5,NUSE2=5)
CU    USES chebev
      REAL*8 xx,c1(7),c2(8),chebev
      SAVE c1,c2
      DATA c1/-1.142022680371172d0,6.516511267076d-3,3.08709017308d-4,
     *-3.470626964d-6,6.943764d-9,3.6780d-11,-1.36d-13/
      DATA c2/1.843740587300906d0,-.076852840844786d0,1.271927136655d-3,
     *-4.971736704d-6,-3.3126120d-8,2.42310d-10,-1.70d-13,-1.d-15/
      xx=8.d0*x*x-1.d0
      gam1=chebev(-1.d0,1.d0,c1,NUSE1,xx)
      gam2=chebev(-1.d0,1.d0,c2,NUSE2,xx)
      gampl=gam2-x*gam1
      gammi=gam2+x*gam1
      return
      END



      FUNCTION chebev(a,b,c,m,x)
      INTEGER m
      REAL*8 chebev,a,b,x,c(m)
      INTEGER j
      REAL*8 d,dd,sv,y,y2
      if ((x-a)*(x-b).gt.0.) pause 'x not in range in chebev'
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
      return
      END




      subroutine rbesjy (ell,x,z,zp,zu,zup,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Riccati-Bessel functions of fractional order  
c     and their first derivatives with respect to x:
c
c     j(ell,x) = cj * exp(ej)
c     y(ell,x) = cy * exp(ey)
c     d/dx j(ell,x) = dj * exp(ej)
c     d/dx y(ell,x) = dy * exp(ey)
c     ----------------------------------------------------------------- 
c
      if (x.le.0.0d0 .or. ell.lt.-0.5d0) stop 'rbesjy 0'       
      v = ell+0.5d0
      call bessjy (v,x,cj,dj,ej,cy,dy,ey)
      pi = acos(-1.d0)
      ex = 0.5d0*dlog(pi*x/2.d0)
      dj = dj+cj/(2.d0*x)
      sj = sqrt(cj*cj+dj*dj)
      cj = cj/sj
      dj = dj/sj
      ej = ej+dlog(sj)+ex
      dy = dy+cy/(2.d0*x)
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = ey+dlog(sy)+ex
     
c add by ZDH
      z = cj * exp(ej)
      zu = cy * exp(ey)
      zp = dj * exp(ej)
      zup = dy * exp(ey)

      return
      end



      subroutine bessjy (v,x,cj,dj,ej,cy,dy,ey)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses a combination of methods (mostly due
c     to Temme) to calculate the Ordinary Bessel functions
c
c     J(v,x) = cj * exp(ej)
c     Y(v,x) = cy * exp(ey)
c
c     and their first derivatives with respect to x
c
c     d/dx J(v,x) = dj * exp(ej)
c     d/dx Y(v,x) = dy * exp(ey)
c
c     for a given real order v >= 0 and real argument x > 0.  
c     Note the exponential scaling, which is used to avoid
c     overflow of Y(v,x) and underflow of J(v,x) for v >> x.
c     ----------------------------------------------------------------- 
c
      parameter (eps = 1.d-15)! consistent with rgamma
      parameter (maxit = 1000)
c
      if (v.lt.0.d0 .or. x.le.0.d0) stop 'bessjy 0'
      pi = acos(-1.d0)
      xmin = 3.d0
      xmax = 5.d0-dlog10(eps)
c
c     begin by calculating Y(a,x) and Y(a+1,x) for |a| <= 1/2
c
      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then
c
c        using Temme's series (bessya) for small x 
c        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
c
         b = x/2.d0
         d = -dlog(b)
         e = a*d
         if (abs(a) .lt. eps) then
            c = 1.d0/pi
         else
            c = a/sin(a*pi)
         endif
         if (abs(e) .lt. eps) then
            s = 1.d0
         else
            s = sinh(e)/e
         endif
         e = exp(e)
         g = e*rgamma(a,p,q)
         e = (e+1.d0/e)/2.d0
         f = 2*c*(p*e+q*s*d)
         e = a*a
         p = g*c
         q = 1.d0/g/pi
         c = a*pi/2.d0
         if (abs(c) .lt. eps) then
            r = 1.d0
         else
            r = sin(c)/c
         endif
         r = pi*c*r*r
         c = 1.d0
         d = -b*b
         ya = f+r*q
         ya1 = p
         do n = 1,maxit
            f = (f*n+p+q)/(n*n-e)
            c = c*d/n
            p = p/(n-a)
            q = q/(n+a)
            g = c*(f+r*q)
            h = c*p-n*g
            ya = ya+g
            ya1 = ya1+h
            del = abs(g)/(1.d0+abs(ya))
            del1 = abs(h)/(1.d0+abs(ya1))
            if (del+del1 .lt. eps) go to 1
         enddo
         stop 'bessjy 1'
   1     f = -ya
         g = -ya1/b
      else if (x.ge.xmin .and. x.lt.xmax) then
c
c        Temme's PQ method (besspqa) for intermediate x  
c        [ N.M.Temme, J Comput Phys 21 (1976) 343-350 ]
c
         c = 0.25d0-a*a
         b = x+x
         p = pi
         e = (x*cos(a*pi)/pi/eps)**2
         p = 1.d0
         q = -x
         r = 1.d0+x*x
         s = r
         do n = 2,maxit
            d = (n-1+c/n)/s
            p = (2*n-p*d)/(n+1)
            q = (-b+q*d)/(n+1)
            s = p*p+q*q
            r = r*s
            if (r*n*n .gt. e) go to 2
         enddo
         stop 'bessjy 2'
   2     p = p/s
         f = p
         q = -q/s
         g = q
         do m = n,1,-1
            r = (m+1)*(2.d0-p)-2.d0
            s = b+(m+1)*q
            d = (m-1+c/m)/(r*r+s*s)
            p = d*r
            q = d*s
            e = f+1.d0
            f = p*e-g*q
            g = q*e+p*g
         enddo
         f = 1.d0+f
         d = f*f+g*g
         pa = f/d
         qa = -g/d
         d = a+0.5d0-p
         q = q+x
         pa1 = (pa*q-qa*d)/x
         qa1 = (qa*q+pa*d)/x
         b = x-pi*(a+0.5d0)/2.d0
         c = cos(b)
         s = sin(b)
         d = sqrt(2.d0/x/pi)
         f = d*(pa*s+qa*c)
         g = d*(qa1*s-pa1*c)
      else if (x .ge. xmax) then
c 
c        and Hankel's asymptotic expansions for large x           
c        [ Abramowitz and Stegun, Section 9.2 ]
c
         p = 0.d0
         q = 0.d0
         do ia = 0,1
            pa = p
            qa = q
            y = 4.d0*(a+ia)**2
            z = 8.d0*x 
            d = 0.d0
            w = -1.d0
            p = 1.d0
            q = 0.d0
            tp = 1.d0
            do k = 1,maxit
               d = d+z
               w = w+2.d0
               tq = +tp*(y-w*w)/d
               q = q+tq
               d = d+z
               w = w+2.d0
               tp = -tq*(y-w*w)/d
               p = p+tp   
               if (abs(tp)+abs(tq) .lt. eps) go to 3
            enddo
            stop 'bessjy 3'
   3        p = p-0.5d0*tp
            q = q-0.5d0*tq
         enddo
         pa1 = p
         qa1 = q
         b = x-pi*(a+0.5d0)/2.d0
         c = cos(b)
         s = sin(b)
         d = sqrt(2.d0/x/pi)
         f = d*(pa*s+qa*c)
         g = d*(qa1*s-pa1*c)
      endif
c
c     now recur upwards from Y(a,x) to Y(v,x),
c     scaling to avoid overflow along the way
c
      p = 0.d0
      if (na .gt. 0) then
         y = 2.d0/x
         do n = 1,na
            h = y*(a+n)*g-f
            f = g
            g = h
   4        if (abs(f) .gt. 4.d0) then 
               p = p+1.d0
               f = 0.0625d0*f
               g = 0.0625d0*g
               go to 4
            endif 
         enddo
      endif 
      cy = f
      dy = (v/x)*f-g
      sy = sqrt(cy*cy+dy*dy)
      cy = cy/sy
      dy = dy/sy
      ey = dlog(sy)+p*dlog(16.d0)
c
c     finally, calculate J(v,x) and dJ(v,x)/dx
c
      vv = max(xmin,v)
      if (x .ge. vv) then
c
c        using upward recursion in the classically allowed region
c
         f = d*(pa*c-qa*s)
         g = d*(qa1*c+pa1*s)
         if (na .gt. 0) then
            y = 2.d0/x
            do n = 1,na
               h = y*(a+n)*g-f
               f = g
               g = h
            enddo
         endif
         cj = f
         dj = (v/x)*f-g
         sj = sqrt(cj*cj+dj*dj)
         cj = cj/sj
         dj = dj/sj
         ej = dlog(sj)
      else
c
c        and CF1 in the classically forbidden region
c        [ Numerical Recipes, 2nd Edition, Section 6.7 ]
c
         ap = 1.d0
         a = v/x
         bp = 0.d0
         b = 1.d0
         f = 0.d0
         g = 0.d0
         y = 2.d0/x
         w = y/pi
         do n = 1,maxit
            an = y*(v+n)*a-ap
            ap = a
            a = an
            bn = y*(v+n)*b-bp
            bp = b
            b = bn
            if (abs(b) .gt. abs(a)) then
               ap = ap/b
               a = a/b
               bp = bp/b
               b = 1.d0
               if (abs(a-f) .lt. eps*abs(f)) then
                  cj = w/(dy-cy*a)
                  dj = a*cj
                  go to 5
               endif
               f = a
            else
               bp = bp/a
               b = b/a
               ap = ap/a
               a = 1.d0
               if (abs(b-g) .lt. eps*abs(g)) then
                  dj = w/(dy*b-cy)
                  cj = b*dj
                  go to 5
               endif
               g = b
            endif
         enddo
         stop 'bessjy 4'
   5     sj = sqrt(cj*cj+dj*dj)
         cj = cj/sj
         dj = dj/sj
         ej = dlog(sj)-ey   
      endif
      return
      end 



      subroutine mbessk (v,x,ck,dk,ek)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This subroutine uses Temme's method [ N.M.Temme, J Comput Phys
c     19 (1975) 324-337 ] to calculate the Modified Bessel function
c
c     K(v,x) = ck * exp(ek)
c
c     and its first derivative with respect to x
c
c     d/dx K(v,x) = dk * exp(ek)
c
c     for a given real order v >= 0 and real argument x > 0.  
c     Note the exponential scaling, which is used to avoid
c     overflow of K(v,x) for v >> x and underflow for v << x.
c     ----------------------------------------------------------------- 
c
      parameter (eps = 1.d-15)! consistent with rgamma
      parameter (maxit = 1000)
c
      if (v.lt.0.d0 .or. x.le.0.d0) stop 'mbessk 0'
      pi = acos(-1.d0)
      xmin = 1.d0
c
c     begin by calculating K(a,x) and K(a+1,x) for |a| <= 1/2
c
      na = int(v+0.5d0)
      a = v-na
      if (x .lt. xmin) then
c
c        using Temme's series for small x 
c
         b = x/2.d0
         d = -dlog(b)
         e = a*d
         c = a*pi
         if (abs(c) .lt. eps) then
            c = 1.d0
         else
            c = c/sin(c)
         endif
         if (abs(e) .lt. eps) then
            s = 1.d0
         else
            s = sinh(e)/e
         endif
         e = exp(e)
         g = e*rgamma(a,p,q)
         e = (e+1.d0/e)/2.d0
         f = c*(p*e+q*s*d)
         e = a*a
         p = 0.5d0*g*c
         q = 0.5d0/g
         c = 1.d0
         d = b*b
         ak = f
         ak1 = p
         do n = 1,maxit
            f = (f*n+p+q)/(n*n-e)
            c = c*d/n
            p = p/(n-a)
            q = q/(n+a)
            g = c*(p-n*f)
            h = c*f
            ak = ak+h
            ak1 = ak1+g
            if (h/ak+abs(g)/ak1 .lt. eps) go to 1
         enddo
         stop 'mbessk 1'
   1     f = ak
         g = ak1/b
         ex = 0.d0
      else if (x .ge. xmin) then
c
c        and Temme's PQ method for large x  
c
         c = 0.25d0-a*a
         g = 1.d0
         f = 0.d0
         e = x*cos(a*pi)/pi/eps
         do n = 1,maxit
            h = (2*(n+x)*g-(n-1+c/n)*f)/(n+1)
            f = g
            g = h
            if (h*n .gt. e) go to 2
         enddo
         stop 'mbessk 2'
   2     p = f/g
         q = p
         b = x+x
         e = b-2.d0
         do m = n,1,-1
            p = (m-1+c/m)/(e+(m+1)*(2.d0-p))
            q = p*(q+1.d0)
         enddo
         f = sqrt(pi/b)/(1.d0+q)
         g = f*(a+x+0.5d0-p)/x
         ex = x
      endif
c
c     now recur upwards from K(a,x) to K(v,x),
c     scaling to avoid overflow along the way
c
      p = 0.d0
      if (na .gt. 0) then
         y = 2.d0/x
         do n = 1,na
            h = y*(a+n)*g+f
            f = g
            g = h
   3        if (abs(f) .gt. 4.d0) then 
               p = p+1.d0
               f = 0.0625d0*f
               g = 0.0625d0*g
               go to 3
            endif 
         enddo
      endif 
      ck = f
      dk = (v/x)*f-g
      sk = sqrt(ck*ck+dk*dk)
      ck = ck/sk
      dk = dk/sk
      ek = dlog(sk)+p*dlog(16.d0)-ex
      return
      end




      function rgamma(x,odd,even)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Direct fortran translation of Temme's algol routine for computing
c     rgamma = 1/Gamma(1-x), along with its odd and even parts, for
c     abs(x) .le. 0.5. [ N.M.Temme, J Comput Phys 19 (1975) 324-337 ]
c     ----------------------------------------------------------------- 
c
      dimension b(12)
      data b / -0.283876542276024d0, -0.076852840844786d0,
     *         +0.001706305071096d0, +0.001271927136655d0,
     *         +0.000076309597586d0, -0.000004971736704d0,
     *         -0.000000865920800d0, -0.000000033126120d0,
     *         +0.000000001745136d0, +0.000000000242310d0,
     *         +0.000000000009161d0, -0.000000000000170d0 / 
      save b
c
      x2 = x*x*8.d0
      alfa = -0.000000000000001d0 
      beta = 0.d0
      do i = 12,2,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      even = (beta/2.d0+alfa)*x2-alfa+0.921870293650453d0
      alfa = -0.000000000000034d0
      beta = 0.d0
      do i = 11,1,-2
         beta = -(2*alfa+beta)
         alfa = -beta*x2-alfa+b(i)
      enddo
      odd = 2*(alfa+beta)
      rgamma = odd*x+even
      return
      end     


      subroutine swap(a1,a2)
      implicit none
      real*8 a1,a2,temp
      temp=a1
      a1=a2
      a2=temp
      endsubroutine

      subroutine swap_int(i,j)
      implicit none
      integer i,j,temp

      temp=i
      i=j
      j=temp

      return
      endsubroutine


      subroutine swap_real(a,b)
      implicit none
      real*8 a,b,temp

      temp=a
      a=b
      b=temp

      return
      endsubroutine


c=================================================================
      function dm(ij,ik,im,beta)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     This function uses eq. (4.1.23) of Edmonds 
c     to calculate the reduced rotation matrix element 
c     d(j,k,m;beta) = <jk|exp(+i*beta*Jy/hbar)|jm>.
c     ----------------------------------------------------------------- 
c     
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)

c
c     half integer angular momenta
c
      sj = half*nint(two*ij)
      sk = half*nint(two*ik)
      sm = half*nint(two*im)
c
c     projection ranges
c
      dm = zero
      if (sk.gt.sj .or. sk.lt.-sj)  return
      if (sm.gt.sj .or. sm.lt.-sj)  return
      if (mod(sj-sk,one) .ne. zero) return
      if (mod(sj-sm,one) .ne. zero) return      
c
c     reflection symmetries
c      
      if (sk+sm .ge. zero) then
        if (sk-sm .ge. zero) then
          tk = sk
          tm = sm
          isign = 0
        else
          tk = sm
          tm = sk
          isign = sk-sm
        endif
      else
        if (sk-sm .ge. zero) then
          tk = -sm
          tm = -sk
          isign = 0
        else
          tk = -sk
          tm = -sm
          isign = sk-sm
        endif
      endif
c
c     evaluation
c
      n = sj-tk
      ia = tk-tm
      ib = tk+tm
      a = ia
      b = ib
      beta2 = half*beta
      cosb2 = cos(beta2)
      sinb2 = sin(beta2)
      cosb = (cosb2-sinb2)*(cosb2+sinb2)
      d1 = pjacob(n,a,b,cosb)
      d2 = cosb2**ib*sinb2**ia
      d3 = d1*d2
      d4 = d3*d3
      ti = tm
      do i = 1,ia
         ti = ti+one
         d4 = d4*(sj+ti)/(sj-ti+one)
      enddo
      d4 = sqrt(d4)
      dm = sign(d4,d3)
      if (mod(isign,2) .ne. 0) dm = -dm
      return
      end

      function pjacob (n,a,b,x)
      implicit double precision (a-h,o-z)
c
c     ----------------------------------------------------------------- 
c     Jacobi polynomial p(n,a,b;x)
c     Abramowitz and Stegun eq. (22.7.1)
c     ----------------------------------------------------------------- 
c
      parameter (zero = 0.0d0)
      parameter (half = 0.5d0)
      parameter (one  = 1.0d0)
      parameter (two  = 2.0d0)
c
      if (n .eq. 0) then
        fp = one
      else
        f = one
        apa = a+a
        apb = a+b
        amb = a-b
        apbamb = apb*amb
        apbp1 = apb+one
        apbp2 = apb+two
        onek = zero
        twok = zero
        fp = half*(amb+apbp2*x)
        do k = 1,n-1
          onek = onek+one
          twok = twok+two 
          a1 = (twok+two)*(onek+apbp1)*(twok+apb)
          a2 = (twok+apbp1)*apbamb
          a3 = (twok+apb)*(twok+apbp1)*(twok+apbp2)
          a4 = (twok+apa)*(onek+b)*(twok+apbp2)
          fm = f
          f = fp
          fp = ((a2+a3*x)*f-a4*fm)/a1
        enddo
      endif
      pjacob = fp
      return
      end


