C Método de calculo de curva de luz de una binaria eclipsante asumiendoestrellasesfericas y orbitas circulares
C Método de Dan Bruton  
C ingresa
C Obtiene
C12345 
      program lcbinaria1

      real*8 MR,MB,RR,RB,LR,LB,INC
      real*8 M1,M2,R1,R2,L1,L2,Q,IR
      real*8 LITE(2000),FASE(2000)
      real*8 phase(2000),LUZ(2000)
      real*8 lux(2000),mag(2000)
      real*8 XR(2000),YR(2000),ZR(2000)
      real*8 XB(2000),YB(2000),ZB(2000)
      real*8  PHI,SMA,factor
      real*8 X1,Y1,Z1,X2,Y2,Z2,X,Y,Z,RHO 
      real*8 A1,A2,CA1,CA2,T1,T2 
      real*8 xft,luxm,fasm
      integer NP,i,J,K,NRO,ND
      real*8 suma,desf,dif,SQR,LUZM
      parameter (Pi=3.14159265) 
C      parameter (NP=400) 
      
     
C      character*80 archivo,saida,errmsg

C Sub SimpleModel (MR, MB, RR, RB, LR, LB, INC, NP, LITE(), XR(), YR(), ZR(), XB(), YB(), ZB())

C ********************************************************
C This subroutine generates a light curve for an eclipsing
C binary star system assuming spherical stars and circular
C orbits.  - Dan Bruton, October 1995, astro@tamu.edu
C ********************************************************
C MR,MB = Star Masses in arbitrary units
C RR,RB = Star Radii in fractions of the semimajor axis length
C LR,LB = Star Luminosities in arbitrary units
C INC = Orbital Inclination
C NP = Number of points to compute for one orbit
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

C       print*,'  Total de filas del archivo de datos: ',ND
C      Total de dados: ND
       REWIND(15)
C Calcula la luz maxima luxm y fase maxima maxphase de los datos observados      
       luxm=-1000.0
       maxf=-1000.0
       open(15, FILE='UOph-magU_sort1.dat', STATUS='OLD') 
       DO i=1,ND
          READ(15,*)phase(i),mag(i)
          xft=-0.4*mag(i)
          lux(i)=10**xft 
           if (lux(i).gt.luxm) then
              luxm = lux(i)
              fasm = phase(i)
           end if
       END DO
       close(15)
  900   format (F5.3,1x,F6.3)
       open(50, FILE='UOph-magnitudU1.dat')        
       DO i=1,ND
          write(UNIT=50,FMT=900)phase(i),lux(i)
       END DO


C Require that R1>R2 to simplify calculations
C PHI = Starting angle for orbit
      print*,' Ingresar los datos con tantas cifras como se indica '
      print *,' Ingrese las masas de las componentes (Masas solares) '
      print *, '  M AZUL, M ROJA ( XX.XXX_XX.XXX)'
      read(*,1000)MB,MR
C      MR=0.6  
C      MB=0.1
 1000 format (F6.3,1x,F6.3)

      print *,' Ingrese los radios de las componentes (Radios solares)'
      print *, 'R AZUL, R ROJA (XX.XXX_XX.XXX)'
      read(*,1100)RB,RR
 1100 format (F6.3,1x,F6.3)
C      RR=0.4  
C      RB=0.3

      print *,' Ingrese las luminosidades de las componentes (en Luminosidad
     $  solar)' 
      print *,'  L AZUL, L ROJA (X.XXX_X.XXX) '
      read(*,1200)LB,LR
 1200 format (F5.3,1x,F5.3)
C      LR=0.6  
C      LB=0.2

      print *,' Ingrese el angulo de inclinacion (grados sexagesimales)' 
      print *,'  INC (XX.XXX)'
      read(*,1300)INC
 1300 format (F6.3)
C      INC=85.0

      print*,''
       If (RR .gt. RB) Then
          M1 = MR
          M2 = MB
          R1 = RR
          R2 = RB
          L1 = LR
          L2 = LB
          PHI = Pi
        Else
          M1 = MB
          M2 = MR
          R1 = RB
          R2 = RR
          L1 = LB
          L2 = LR
          PHI = 0
        End If

      print*,' ------------------------------------------'
      print*,'  '
      print*,' Masas:  M1 ',M1,' M2 ',M2
      print*,' Radios:  R1 ',R1,' R2 ',R2
      print*,' Masas:  L1 ',L1,' L2 ',L2
      print*,' Inclinacion : I ',INC  
      Print*,' PHI : ',PHI
      print*,'  '
      print*,' Calculo de la razon de masas Q=M2/M1 '
      print*,'  '
      Q=M2/M1 
      print*,'  Q = ',Q 
      print*,'  '   
      print*,' ------------------------------------------'
      print*,'  '   

C Semimajor Axis
      SMA = 1
      LUZM=-1000
C Orbital Inclination
      IR = INC * Pi / 180
      open(75, FILE='salidalc1.txt', ACCESS='sequential')
        print*,'  '
C        print*,'FASE X  Y Z  X1  Y1 Z1 X2 Y2 Z2'
C Compute Star Positions and Light Curve
      NP=2*ND
      DO J=1, NP+1

C T:Phase Angle
         T = PHI + (J - 1)*(2*Pi/NP) 
         X = SMA * sin(T)
         Y = SMA * cos(IR)*cos(T)
         Z = SMA * sin(IR)*cos(T)
  
         X1 = -1 * X/((1/Q) + 1)
         Y1 = -1 * Y/((1/Q) + 1)
         Z1 = -1 * Z/((1/Q) + 1)

         X2 = X/(1 + Q)
         Y2 = Y/(1 + Q)
         Z2 = Z/(1 + Q)


C         print*,T,' ',X,Y,Z,' ',X1,Y1,Z1,' ',X2,Y2,Z2

C Apparent Distance Between Stars
         RHO = sqrt((X2 - X1)**2 + (Y2 - Y1)**2)

C         print*,' fase ',T,' RHO ',RHO
    
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
C      factor=luxm/LUZM  
C      DO J=1, NP+1
C        LUZ(J)=LITE(J)*factor 
C        PRINT 2200,FASE(J),LUZ(J) 
C 2200   format (F5.3,1x,F6.3)
C        write(UNIT=75,FMT=2200)FASE(J),LUZ(J)
C      END DO

         
C Calcula la luz maxima luxm y fase maxima maxphase de los datos observados
      suma=0
      factor=luxm/LUZM
      K=0
      DO J=1, NP + 1
         DO i = 1, ND
C           READ(15,*)phase(i),mag(i)
            LUZ(J)=LITE(J)*factor 
            desf=abs(phase(i)-FASE(J))
            if (LITE(J).ne.LITE(J)) then
                LUZ(J)=1000.
                 desf =0.001 
            end if
           if (desf.lt.0.0005) then
C             print*, i,' desf ',desf
              dif=(lux(i)-LUZ(J)) 
              suma = suma + (dif*dif)
              K=K+1
             print*,i,'dato ',phase(i),'model',J,' ',FASE(J),' ',desf
             print*,i,'d',lux(i),' ',J,' ',LUZ(J),' ',dif,'',suma      
             print*,' '          
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
      close(50) 
      
      END


