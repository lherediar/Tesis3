C234567
C Usa la subrutina AMOEBA para calcular el minimo de una funcion (version2)
C      Julio Tello
C      Set 2021
C      Calcula para NP=5 parametros de entrada: Razon de Masas,Radio1,Radio2,Luminosidad relativa,incl
C      Calcula el sigma 
       program be_amoeba
C     driver for routine amoeba 
       INTEGER NP,MP
       REAL FTOL
       PARAMETER(NP=5,MP=6,FTOL=1.0E-6)
C       PARAMETER(NP=6,MP=7,FTOL=1.0E-6)
       INTEGER i,iter,j,ndim
       REAL famoeb,p(MP,NP),x(NP),y(MP)
       real fhase(2000),mhag(2000),lft,luz(2000)
       integer NPD,NRO
       EXTERNAL famoeb
C ------------------------------------------------------
C      DATA (Matriz P de MP filas e NP colunas,  MP = NP + 1 y NP es el numero de variables a ajustar)
C      lee los elementos de la matriz que contiene los MP puntos de NP dimensiones
       open(12, FILE='matriz-be.txt', STATUS='OLD') 
       do i=1,MP
          read(12,10) p(i,1),p(i,2),p(i,3),p(i,4),p(i,5)
       end do 
 10    format (F7.3,F7.3,F7.3,F7.3,F7.3)
       close (12)
C      lee en orden: Q,R1,R2,L1R,INCL   (R1 > R2  y Q < = 1.0  L >= 1.0)
C       0.02 <  Q  < 10.0 , 0.01 < R/RSOL < 1,   0.1 < Lrel/LSOL < 1,  85 < incl < 90
C      Asigna valores de X1,X2,X3,X4,X5 a las componentes de P
       ndim=NP
       do 12 i=1,MP
          do 11 j=1,NP
             x(j)=p(i,j)
 11       continue
          y(i)=famoeb(x)
 12    continue
       print*,'  '       
       write(*,'(/1x,a)') 'Vertices of initial 3-D simplex and'
       write(*,'(1x,a)') 'function values at the vertices:'
C      Q = M2/M1  LR1 = L1/L2   Sr2= suma de cuadrados de residuales       
       write(*,16)' Q     R1     R2     L1R     INCL    SR2 '        
       do 14 i=1,MP 
          write(*,'(1x,i3,6f8.3)') i,(p(i,j),j=1,NP),y(i)   
 14    continue
       print*,'  '
 16    format (7X,a)

       call amoeba(p,y,MP,NP,ndim,FTOL,famoeb,iter)
       write(*,'(/1x,a,i3)') 'Number of iterations: ',iter
       write(*,'(/1x,a)') 'Vertices of final 3-D simplex and'
       write(*,'(1x,a)') 'function values at the vertices:'
       write(*,16)'    Q     R1     R2     L1R     INCL    SR2 '      
C     $Lum1     Lum2    Incl    Sr²'   
       do 18 i=1,MP
       write(*,'(1x,i3,6f8.4)') i,(p(i,j),j=1,NP),y(i)
 18    continue
C       write(*,'(/1x,a)') 'True minimum is at (0.5,0.6,0.7)' 



C   Lee el numero de filas del archivo datosentrada
       open(15, FILE='UOph-magU_sort1.dat', STATUS='OLD') 
       NRO=0 
       i=1
       NPD=1
       DO WHILE (NRO.eq.0)
         READ(15,*,IOSTAT=NRO)
         IF (NRO.eq.0) THEN
             i=i+1  
         END IF
       END DO
       NPD=i-1
C       print*,'  Total de filas del archivo de datos: ',NPD
C      Total de dados: ND
       REWIND(15)
       open(15, FILE='UOph-magU_sort1.dat', STATUS='OLD') 
       DO i=1,NPD
          READ(15,*)fhase(i),mhag(i)
          lft=-0.4*mhag(i)
          luz(i)=10**lft 
       END DO
       close(15)
 1200   format (F5.3,1x,F6.3)
       open(80, FILE='UOph-magnitudU.dat')        
       DO i=1,NPD
          write(UNIT=80,FMT=1200)fhase(i),luz(i)
       END DO
       close(80)

       END

       REAL FUNCTION famoeb(x)
       REAL x(5)       
       real*8 ang,qmas
       real*8 rad1,rad2,lumr,resp
C       print *,'Llama a la subrutina'
        qmas=x(1)
        rad1=x(2)
        rad2=x(3)
        lumr=x(4)          
        ang=x(5)
C        print*,qmas,' ',rad1,' ',rad2,' ',lumr,' ',ang 
        call lcbinaria(qmas,rad1,rad2,lumr,ang,resp)
C        close(15)
C        close(10)        
        print*,'      Sr²',resp     
       famoeb=resp
       end

      
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,iter)
      INTEGER iter,mp,ndim,np,NMAX,ITMAX
      REAL ftol,p(mp,np),y(mp),funk,TINY
      PARAMETER (NMAX=20,ITMAX=5000,TINY=1.e-10) 
C     Maximum allowed dimensions and function evaluations, and a small number.
      EXTERNAL funk

C     USES amotry,funk
C     Multidimensional minimization of the function funk(x) where x(1:ndim) is a vector
C     in ndim dimensions, by the downhill simplex method of Nelder and Mead. The matrix
C     p(1:ndim+1,1:ndim) is input. Its ndim+1 rows are ndim-dimensional vectors which are
C     the vertices of the starting simplex. Also input is the vector y(1:ndim+1), whose 
C     components must be pre-initialized to the values of funk evaluated at the ndim+1 vertices (rows)
C     of p; and ftol the fractional convergence tolerance to be achieved in the function value
C     (n.b.!). On output, p and y will have been reset to ndim+1 new points all within ftol of
C     a minimum function value, and iter gives the number of function evaluations taken.


      INTEGER i,ihi,ilo,inhi,j,m,n
      REAL rtol,sum,swap,ysave,ytry,psum(NMAX),amotry
      iter=0
1      do 12 n=1,ndim 
C     Enter here when starting or have just overall contracted.
       sum=0. 
C       Recompute psum.
        do 11 m=1,ndim+1
        sum=sum+p(m,n)
11       enddo
       psum(n)=sum
12      enddo
2      ilo=1 
C     Enter here when have just changed a single point.
      if (y(1).gt.y(2)) then 
C     Determine which point is the highest (worst), next-highest,
       ihi=1 
C       and lowest (best),
       inhi=2
      else
	ihi=2
	inhi=1
      endif
      do 13 i=1,ndim+1 
C     by looping over the points in the simplex.
       if(y(i).le.y(ilo)) ilo=i
       if(y(i).gt.y(ihi)) then
        inhi=ihi
        ihi=i
       else if(y(i).gt.y(inhi)) then
        if(i.ne.ihi) inhi=i
       endif
13      enddo
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
C     Compute the fractional range from highest to lowest and return if satisfactory.
      if (rtol.lt.ftol) then  
C     If returning, put best point and value in slot 1.
       swap=y(1)
       y(1)=y(ilo)
       y(ilo)=swap
       do 14 n=1,ndim
        swap=p(1,n)
        p(1,n)=p(ilo,n)
        p(ilo,n)=swap
14       end do
       return
      endif
      if (iter.ge.ITMAX) then
        pause 
      end if   
C      ’ITMAX exceeded in amoeba’
        iter=iter+2
C      Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex across
C     from the high point, i.e., reflect the simplex from the high point.
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0)
      if (ytry.le.y(ilo)) then
C     Gives a result better than the best point, so try an additional extrapolation by a factor 2.
       ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0)
      else if (ytry.ge.y(inhi)) then
C      The reflected point is worse than the second-highest, so look for an intermediate lower point,
C      i.e., do a one-dimensional contraction.
       ysave=y(ihi)
       ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5)
       if (ytry.ge.ysave) then 
C       Can’t seem to get rid of that high point. Better contract
        do i=1,ndim+1 
C       around the lowest (best) point.
         if(i.ne.ilo)then
          do j=1,ndim
           psum(j)=0.5*(p(i,j)+p(ilo,j))
           p(i,j)=psum(j)
          enddo
          y(i)=funk(psum)
         endif
        enddo
       iter=iter+ndim 
C      Keep track of function evaluations.
       goto 1 
C      Go back for the test of doneness and the next iteration.
       endif
      else
       iter=iter-1
C      Correct the evaluation count.
      endif
      goto 2
      END
       
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac)
      INTEGER ihi,mp,ndim,np,NMAX
      REAL amotry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)
      EXTERNAL funk

C     USES funk
C     Extrapolates by a factor fac through the face of the simplex across from the high point,
C     tries it, and replaces the high point if the new point is better.

      INTEGER j
      REAL fac1,fac2,ytry,ptry(NMAX)
      fac1=(1.-fac)/ndim
      fac2=fac1-fac
      do j=1,ndim
       ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
       if (ptry(j).lt.0) then
         ptry(j) = abs(ptry(j))
         print*,' cambia de signo j =',j  
       end if  
       if (ptry(1).lt.0.1 ) then
          ptry(1)=0.1
       end if
      if (ptry(5).gt.90.0 ) then
          ptry(5)=89.5
       end if
      end do

      ytry=funk(ptry) 
C     Evaluate the function at the trial point.
      if (ytry.lt.y(ihi)) then 
C     If it’s better than the highest, then replace the highest.
       y(ihi)=ytry
       do 12 j=1,ndim
        psum(j)=psum(j)-p(ihi,j)+ptry(j)
        p(ihi,j)=ptry(j)
12       end do
      end if
      amotry=ytry
      return
      END

C Método de calculo de curva de luz de una binaria eclipsante asumiendo estrellas esfericas y orbitas circulares
C Método de Dan Bruton  
C ingresa masas,radios,luminosidades y angulo de inclinacion
C Obtiene curva de luz y posiciones de las componentes
C1234567 
       subroutine lcbinaria(QQ,RB,RR,LREL,INC,SQR)
       real*8 RR,RB,LREL,INC
C       real*8 M1,M2
       real*8 R1,R2,L1,L2,QQ,IR
       real*8 LITE(2000),FASE(2000),LUZM,LUZ(2000)
       real*8 XR(2000),YR(2000),ZR(2000)
       real*8 XB(2000),YB(2000),ZB(2000)
       real*8 PHI,SMA,SQR
       real*8 X1,Y1,Z1,X2,Y2,Z2,X,Y,Z,RHO 
       real*8 A1,A2,CA1,CA2,T1,T2 
       integer PT,ND,NRO
       integer i,J,K
       real*8 phase(2000),LUX(2000),mag(2000)
       real*8 dif,suma,luxm,fasm,factor,desf
       real*8 maxluz,sumfase
       parameter (Pi=3.14159265) 
C       parameter (PT=1000) 
C      character*80 archivo,saida,errmsg

C Sub SimpleModel (MR, MB, RR, RB, LR, LB, INC, PT, LITE(), XR(), YR(), ZR(), XB(), YB(), ZB())

C ********************************************************
C This subroutine generates a light curve for an eclipsing
C binary star system assuming spherical stars and circular
C orbits.  - Dan Bruton, October 1995, astro@tamu.edu
C ********************************************************
C MR,MB = Star Masses in arbitrary units
C RR,RB = Star Radii in fractions of the semimajor axis length
C LR,LB = Star Luminosities in arbitrary units
C INC = Orbital Inclination
C PT = Number of points to compute for one orbit
C LITE = Array containing the light curve
C XR,YR,ZR,XB,YB,ZB = Star Positions
C
C Note: The second character,"R" or "B", means red or blue star
C which are the colors used for the animations of the stars in orbit.
C ********************************************************

C   Lee el numero de filas del archivo datosentrada
       open(15, FILE='UOph-magU_sort1.dat', STATUS='OLD') 
       NRO=0 
       i=1
       ND=1
       DO WHILE (NRO.eq.0)
         READ(15,*,IOSTAT=NRO)
         IF (NRO.eq.0) THEN
             i=i+1  
         END IF
       END DO
       ND=i-1
       if (ND.ge.1000) then
           print*,' el numero de datos ND es mayor o igual a 1000'
           stop 
       end if 
C       print*,'  Total de filas del archivo de datos: ',ND
C      Total de dados: ND
       REWIND(15)
C Calcula la luz maxima luxm y fase maxima maxphase de los datos observados      
       luxm=-1000.0
       sumfase=0.0
       k = 0
       open(15, FILE='UOph-magU_sort1.dat', STATUS='OLD') 
       DO i=1,ND
          READ(15,*)phase(i),mag(i)
          xft=-0.4*mag(i)
          lux(i)=10**xft 
           if (lux(i).gt.luxm) then
              luxm = lux(i)
              fasm = phase(i)
           end if
           if (phase(i).gt.0.2) then
              if(phase(i).lt.0.3) then
                 sumfase = sumfase + lux(i) 
                 k = k+ 1                
              end if 
           end if
            maxluz = sumfase/k
       END DO
       close(15)
 1100   format (F5.3,1x,F6.3)
       open(50, FILE='UOph-magnitudU.dat')        
       DO i=1,ND
          write(UNIT=50,FMT=1100)phase(i),lux(i)
       END DO

       If (RR .gt. RB) Then
CC         MR > MB
CC         M1 = MR  ;   M2 = MB
          Q = QQ  
          if (QQ .gt. 1.0) then
               Q =  1/QQ
          end if       
          R1 = RR
          R2 = RB
          L1 = LREL
          L2 = 1.0 
          PHI = Pi
        Else
CC         MB > MR        
C          M1 = MB   ;   M2 = MR 
          Q = QQ  
          if (QQ .gt. 1.0) then
               Q =  1/QQ
          end if
          R1 = RB
          R2 = RR
          L1 = LREL
          L2 = 1.0
          PHI = 0
       End If

C       print*,' Razon de Masas:  Q ',QQ
C       print*,' Radios:  R1 ',R1,' R2 ',R2
C       print*,' Luminosidad :  Lrel ',LREL
C       Print*,' PHI : ',PHI
      PRINT*, '---------------------------------------------------',
     $'------------------'
       print*,'Q ',Q,' R1 ',R1,' R2 ',R2,' Lrel ',LREL,' Inc ',INC
     $,' PHI ',PHI
      
C Semimajor Axis
       SMA = 1
       LUZM=-1000
C Orbital Inclination
       IR = INC * Pi / 180
C       open(75, FILE='salidabe.txt', ACCESS='sequential')
       open(75, FILE='salidabe.txt', STATUS='OLD')

C Compute Star Positions and Light Curve
       PT = 2*ND
       DO J = 1, PT +1

C T:Phase Angle
          T = PHI + (J - 1)*(2*Pi/PT) 
          X = SMA * Sin(T)
          Y = SMA * Cos(IR)*Cos(T)
          Z = SMA * Sin(IR)*cos(T)
  
          X1 = -1 * X/((1/Q) + 1)
          Y1 = -1 * Y/((1/Q) + 1)
          Z1 = -1 * Z/((1/Q) + 1)

          X2 = X/(1 + Q)
          Y2 = Y/(1 + Q)
          Z2 = Z/(1 + Q)

C Apparent Distance Between Stars
          RHO = sqrt((X2 - X1)**2 + (Y2 - Y1)**2)

C          if (J .eq. ND) then
C           print*,' R ',RHO,' Q ',Q,' ',Z1,' ',Z2
C           print*,X1,Y1,Z1,' ',X2,Y2,Z2
C          end if

    
C Area of Star Disks (not during an eclipse)

          A1 = Pi*(R1**2)
          A2 = Pi*(R2**2)
    
C Total and Annular Eclipses

          If (RHO .lt. (R1 + R2)) Then
             If (RHO .lt. (R1 - R2)) Then
                If (Z1 .gt. Z2) Then
C Total Eclipse
                   A2 = 0   
                Else
C Annular Eclipse
                   A1 = (A1 - A2) 
                End If
             Else

C Partial Eclipses
C   The area covered during partial eclipses is determined
C by considering and area of a circle cut by a line segment:
C Area = r**2 (theta - Sin(theta)) / 2 where r is the radius
C of the circle and sin(theta)=h/r.  See math handbooks.

                H = ((R1**2) * (R2**2))
                H = H - (.25 * (((RHO**2) - (R1**2) - (R2**2))**2))
                If ((H .NE. R2).and.(H .NE. R1).and.(H .ge. 0)) Then
                   H = sqrt(H / (RHO**2))
                   T1 = 2 * Atan(H / (sqrt((R1**2) - (H**2))))
                   T2 = 2 * Atan(H / (sqrt((R2**2) - (H**2))))
                End If
                CA1 = (R1**2) * (T1 - Sin(T1)) / 2
                CA2 = (R2**2) * (T2 - Sin(T2)) / 2
C Shallow Partial Eclipse
                If (RHO .gt. (sqrt((R1**2)-(R2**2)))) Then
                   If (Z1 .gt. Z2) Then
                      A2 = (A2 - CA1 - CA2)
                   Else
                      A1 = (A1 - CA1 - CA2)
                   End If
C' Deep Partial Eclipse
                Else
                   If (Z1 .gt. Z2) Then
                      A2 = (CA2 - CA1)
                   Else
                      A1 = (A1 - A2 + CA2 - CA1)
                   End If
                End If
             End If
          End If
C Approximate Light Intensity (no limb darkening)
C Calcula la luminosidad maxima del modelo LUZM para calular alli la fase
        LITE(J) = ((L1 * A1 / R1**2) + (L2 * A2 / R2**2)) / (4 * Pi)
        FASE(J)=T/(2*Pi)
        if (LITE(J).gt.LUZM) then
           LUZM = LITE(J)
        end if
        If (RR .gt. RB) Then
           XR(J) = X1
           YR(J) = Y1
           ZR(J) = Z1
           XB(J) = X2
           YB(J) = Y2
           ZB(J) = Z2
         Else
          XR(J) = X2
          YR(J) = Y2
          ZR(J) = Z2
          XB(J) = X1
          YB(J) = Y1
          ZB(J) = Z1
        End If
      END DO
C      PRINT *,'LUZ MAX ',LUZM
C Calcula la luz maxima luxm y fase maxima maxphase de los datos observados
      suma=0
      factor = maxluz/LUZM
C      factor=luxm/LUZM
      K=0
      DO J=1, PT + 1
         DO i = 1, ND
C           READ(15,*)phase(i),mag(i)
            LUZ(J)=LITE(J)*factor 
            desf=abs(phase(i)-FASE(J))
            if (LITE(J).ne.LITE(J)) then
                LUZ(J)=1000.
                 desf =0.1 
            end if
           if (desf.lt.0.0005) then
C             print*, i,' desf ',desf
              dif=(lux(i)-LUZ(J)) 
C             print*,i,'d ',phase(i),'m',J,' ',FASE(J),' ',desf,' ',dif
C             print*,i,'d',desf,'m',J,'',lux(i),'',LUZ(J),' ',dif,'',suma
              suma = suma + (dif*dif)
              K=K+1
           end if
         END DO
C       print*,J,' ',FASE(J),' ',LUZ(J)
C        LUZ(J)=LITE(J)/LUZM 
C        PRINT 2200,FASE(J),LUZ(J) 
 2200   format (F5.3,1x,F6.3)
        write(UNIT=75,FMT=2200)FASE(J),LUZ(J)
      END DO
      SQR=suma/K
      print*,' SQR =',SQR,' puntos ',K 
      close(75)
      close(50, status='delete') 
      END


