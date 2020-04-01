C       ESTE PROGRAMA CONSTROI MODELOS COM NP PARTICULAS A PARTIR DE
C       UM MODELO DE UM ESFEROIDE SEGUINDO UM PERFIL DE NFW.

         PROGRAM NFW
         IMPLICIT NONE        
         INTEGER I,J,K,N,NMAX,NMODELS,ISEMA,ISEMB
         INTEGER JMAX,N200,NAUX,NRAD
         PARAMETER(NMAX=100000)
         PARAMETER(JMAX=400)
         DOUBLE PRECISION GASDEV
         EXTERNAL GASDEV
         DOUBLE PRECISION RAN1
         EXTERNAL RAN1
         DOUBLE PRECISION RATIO
         EXTERNAL RATIO
         DOUBLE PRECISION VC2
         EXTERNAL VC2
         DOUBLE PRECISION C,RC,R200,RMAX,DM,RS
         DOUBLE PRECISION CC,Ze,X1,X2,R(NMAX),PI
         DOUBLE PRECISION ASCRT,DECLI,DIST(NMAX),PHI
         DOUBLE PRECISION COSTHETA,SENTHETA,VROT
         DOUBLE PRECISION X(NMAX),Y(NMAX),Z(NMAX)
         DOUBLE PRECISION VX(NMAX),VY(NMAX),VZ(NMAX)
         DOUBLE PRECISION VIX,VIY,VIZ,VVX,VVY,VVZ
         DOUBLE PRECISION H0,ZRED,DC,M200,G,HZ
         DOUBLE PRECISION H,VDISP,Q1,V200,QZERO
         DOUBLE PRECISION RHOCRIT,RHOC,RTBIS,Q2
         DOUBLE PRECISION MASS,Mmin,Mmax,LOGN200
         DOUBLE PRECISION Zmin,Zmax,DZ,XACC,OMEGAL,OMEGAM
         DOUBLE PRECISION YASCRT(NMAX),ZDECL(NMAX),ZGAL(NMAX)
         DOUBLE PRECISION ANGAR(NMAX),ANGDECLI(NMAX),MPART1
         DOUBLE PRECISION KIN(NMAX),POT(NMAX),KINT,POTT
         DOUBLE PRECISION V1,V2,V3,V4,V5,MPART2
         DOUBLE PRECISION RR,RF,DR,ALPHA,BETA,GAMA
         DOUBLE PRECISION PARTE1,PARTE2,PARTE3,VLOS(NMAX)
CC        log(N200) = 0.47 * (log M200 - 14.5) + 1.58
C
CC        LOGN200 = 0.47*(LOG10(M200) - 14.5) + 1.58

         PRINT*, "DIGITE A SEMENTE ALEATORIA"
         READ*, ISEMA
         ISEMA = -ISEMA
         ISEMB = ISEMA - 20
         ASCRT = 24.0
         DECLI = 90.0
         OMEGAL = 0.7
         OMEGAM = 1. - OMEGAL
C OMEGARADIATION  IS VERY SMALL TO BE ACCOUNTED HERE.        
         H0 = 0.069 !km/s/kpc 
         C = 299792.458
         G = 43007.1
         PI = ACOS(-1.)
         ALPHA = 2.*PI*RAN1(ISEMA)
         BETA  = 2.*PI*RAN1(ISEMB)
         ISEMA = ISEMA+ISEMB
         GAMA  = 2.*PI*RAN1(ISEMA)
         Mmin = 4.0
         Mmax = 5.5
         DM = Mmax - Mmin
         Zmin = 0.03
         Zmax = 0.13
         DZ = (Zmax - Zmin)
C         Ze = Zmin + RAN1(ISEMA)*DZ
C         HZ = H0*SQRT(OMEGAM*(1.+Ze)**3 + OMEGAL)
         XACC = 1.0E-03
         PRINT*, "DIGITE O NUMERO DESEJADO DE MODELOS"
         READ*, NMODELS
        
         DM = DM/FLOAT(NMODELS)

         DO J = 1,NMODELS
         Ze = Zmin + RAN1(ISEMA)*DZ
         HZ = H0*SQRT(OMEGAM*(1.+Ze)**3 + OMEGAL) 
         QZERO = OMEGAM/2. - OMEGAL        
         M200 = Mmin + (J-1)*DM                
         M200 = 10.**(M200)   
C SEGUNDO A PRESCRICAO SEMIEMPIRICA DE ANDREON, 2008         
         LOGN200 = 0.47*(LOG10(M200*10.**(10)) - 14.5) + 1.58
         LOGN200 = 10.**(LOGN200) 
         N200 = NINT(LOGN200)                       
         V200 = 10*G*HZ*M200
         V200 = LOG(V200)/3.
         V200 = EXP(V200)
         R200 = V200/(10.*HZ)
cc COMERFORD AND NATARAYAN, 2008, MNRAS,379, Issue 1, July 2007 
         CC = (9./(1.+Ze))*(M200/1300.)**(-0.13)     
         VDISP = SQRT(G*M200/R200)
         RMAX = 2.5*R200
         RS = R200/CC
C         PRINT*, "N200,CC,RS,R200,M200,V200,VDISP"
C         PRINT*, N200,CC,RS,R200,M200,V200,VDISP
C         WRITE(9,*) J, N200, R200, M200, VDISP
         ASCRT = ASCRT*RAN1(ISEMA)
         DECLI = DECLI*(1.0 - 2.*RAN1(ISEMA))
         
         DO I=1,N200
         X1 = 0.0
         X2 = RMAX
 11      Q1 = RAN1(ISEMA)
         CALL RTBISS(X1,X2,XACC,I,Q1,RTBIS)
         R(I) = RS*RTBIS
         IF(R(I).GT.R200)THEN
         GO TO 11
         END IF        
         COSTHETA = 2.*(RAN1(ISEMA)-0.5)
         SENTHETA = SQRT(1.0 - COSTHETA**2)
         PHI = 2.*PI*RAN1(ISEMA)
         X(I) = R(I)*SENTHETA*COS(PHI)
         Y(I) = R(I)*SENTHETA*SIN(PHI)
         Z(I) = R(I)*COSTHETA
         RR = R(I)
         VROT = VC2(CC,R200,V200,G,RR)
         VX(I) = -VROT*SIN(PHI) 
         VY(I) =  VROT*COS(PHI)
         VZ(I) =  0.0                     
         WRITE(12,*) RR,VROT 
C-----------------------------------------------------
C - GERANDO ESFERAS COM ORIENTACOES DO EIXO DE ROTACAO ALEATORIOS
C - PARA PAPERS FUTUROS         

C CRIANDO UM MOVIMENTO DE ROTACAO DAS GALAXIAS-MEMBROS

C         VIX = -VROT*SIN(PHI)
C         VIY =  VROT*COS(PHI)
C         VIZ = 0.0                     
C MODELO EM ROTACAO COMPLETA A PARTIR DOS
C SEGUINTES PROCEDIMENTOS, COM EXCECAO DOS OBJETOS DO CAMPO
C CHAMAMOS A SUBROTINA COM OS ANGULOS DE EULER ATRIBUIDOS
C PREVIAMENTE.         
                    
C          CALL ROT(ALPHA,BETA,GAMA,VIX,VIY,VIZ,VVX,VVY,VVZ)
C          VX(I) = VVX
C          VY(I) = VVY
C          VZ(I) = VVZ
          

C         VX(I) = (1./SQRT(3.))*GASDEV(ISEMA,ISEMB)*VDISP
C         VY(I) = (1./SQRT(3.))*GASDEV(ISEMA,ISEMB)*VDISP
C         VZ(I) = (1./SQRT(3.))*GASDEV(ISEMA,ISEMB)*VDISP
C------------------------------------------------------------------------------------
          WRITE(11,*) J,X(I),Y(I),Z(I),VX(I),VY(I),VZ(I),0
C-----CALCULANDO A VELOCIDADE AO LONGO DA LINHA DE VISADA

C-----REVISAR VLOSS-PARA FUTUROS PAPERS-----------------------------
          
          
           
cC         WRITE(7,*) J,Y(I),Z(I),VX(I)

         DIST(I) = C*(Ze/H0)*(1.0 - ((1.0 + QZERO)/2.)*Ze) 
         PARTE1 = VX(I)*(X(I)-DIST(I))+VY(I)*Y(I)+VZ(I)*Z(I)
C         PARTE2 = SQRT(VX(I)**2 + VY(I)**2 + VZ(I)**2)
         PARTE3 = SQRT((X(I)-DIST(I))**2+Y(I)**2+Z(I)**2)
         VLOS(I) = PARTE1/(PARTE3)

         WRITE(7,*) J,Y(I),Z(I),VLOS(I),DIST(I),0


         ANGAR(I) = ATAN(Y(I)/DIST(I))*180./PI 
         ANGDECLI(I) =  ATAN(Z(I)/DIST(I))*180./PI
         YASCRT(I) = ASCRT + ANGAR(I)
         ZDECL(I) = DECLI + ANGDECLI(I)
         ZGAL(I) = Ze + (VLOS(I)/C)*(1.+Ze)          
         WRITE(8,*) J,YASCRT(I),ZDECL(I),ZGAL(I),0              
         ENDDO

C CRIANDO A CURVA DE ROTACAO PARA O MODELO         

         RR = 0.0
         RF = R200     
         WRITE(10,*) J,RR,RR
C PARA MONTAR A CURVA DE ROTACAO PARTIMOS DO PRINCIPIO QUE NA ORIGEM
C DO SISTEMA A VELOCIDADE CIRCULAR SERAH NULA.         
         NRAD = 100
         DR = (RF - RR)/FLOAT(NRAD)         
         DO I=1,NRAD
         RR = RR + DR
         VROT = VC2(CC,R200,V200,G,RR)
C ESCREVENDO EM FORT.10 O NUMERO DO MODELO E SUA CURVA DE ROTACAO
C CURVA DE ROTACAO RAIO EM kpc E VROT EM km/s. ADICIONAMOS O SOFTENING 
C PARAMETER, USADO EM N-CORPOS        
         WRITE(10,*) J,RR,VROT
         
         ENDDO

C REGIAO ALEM DE R200 ... =( SAO AS GALAXIAS DO CAMPO.
C SUGERIMOS QUE AS GALAXIAS DO CAMPO NAO ESTAO GRAVITACIONALMENTE
C LIGADAS AO CLUSTER.         

         
         Q1 = RATIO(CC)
         NAUX = NINT(N200*Q1)
         PRINT*, Q1, NAUX
         DO I = (N200+1),NAUX
         X1 = 0.0
         X2 = RMAX
         Q2 = RATIO(CC)
 12      Q1 = 1. + (Q2 - 1.0)*RAN1(ISEMA)
         CALL RTBISS(X1,X2,XACC,I,Q1,RTBIS)
         R(I) = RS*RTBIS
         IF((R(I).LT.R200).OR.(R(I).GT.RMAX))THEN
         GO TO 12
         END IF
         PRINT*, Q1,RAN1(ISEMA),R(I)
         COSTHETA = 2.*(RAN1(ISEMA)-0.5)
         SENTHETA = SQRT(1.0 - COSTHETA**2)
         PHI = 2.*PI*RAN1(ISEMA)
         X(I) = R(I)*SENTHETA*COS(PHI)
         Y(I) = R(I)*SENTHETA*SIN(PHI)
         Z(I) = R(I)*COSTHETA
         VX(I) = (1./SQRT(3.))*GASDEV(ISEMA,ISEMB)*RAN1(ISEMA)*VDISP
         VY(I) = (1./SQRT(3.))*GASDEV(ISEMA,ISEMB)*RAN1(ISEMA)*VDISP
         VZ(I) = (1./SQRT(3.))*GASDEV(ISEMA,ISEMB)*RAN1(ISEMA)*VDISP
         WRITE(11,*) J,X(I),Y(I),Z(I),VX(I),VY(I),VZ(I),1
                

         DIST(I) = C*(Ze/H0)*(1.0 - ((1.0 + QZERO)/2.)*Ze)
         PARTE1 = VX(I)*(X(I)-DIST(I))+VY(I)*Y(I)+VZ(I)*Z(I)
C         PARTE2 = SQRT(VX(I)**2 + VY(I)**2 + VZ(I)**2)
         PARTE3 = SQRT((X(I)-DIST(I))**2+Y(I)**2+Z(I)**2)
         VLOS(I) = PARTE1/(PARTE3)

         WRITE(7,*) J,Y(I),Z(I),VLOS(I),DIST(I),1


         ANGAR(I) = ATAN(Y(I)/DIST(I))*180./PI
         ANGDECLI(I) =  ATAN(Z(I)/DIST(I))*180./PI
         YASCRT(I) = ASCRT + ANGAR(I)
         ZDECL(I) = DECLI + ANGDECLI(I)
         ZGAL(I) = Ze + (VLOS(I)/C)*(1.+Ze)
         WRITE(8,*) J,YASCRT(I),ZDECL(I),ZGAL(I),1      

         ENDDO
         
C         MPART1 = M200/FLOAT(N200) 
C         DO I =1,NAUX
C         KIN(I) = VX(I)**2
C         POT(I) = G*MPART1/SQRT(Y(I)**2+Z(I)**2)         
C         ENDDO
         
          
C         KINT = 0.0
C         DO I=1,N200
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(0.5*R200))THEN
C         KINT = KIN(I) + KINT
C         ENDIF
C         ENDDO
                  
C         POTT = 0.0
C         DO I=1,N200
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(0.5*R200))THEN
C         POTT = POT(I) + POTT
C         ENDIF
C         ENDDO         
C         V1 = KINT/POTT

C         KINT = 0.0
C         DO I=1,N200
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(R200))THEN
C         KINT = KIN(I) + KINT
C         ENDIF
C         ENDDO

C         POTT = 0.0
C         DO I=1,N200
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(R200))THEN
C         POTT = POT(I) + POTT
C         ENDIF
C         ENDDO
C         V2 = KINT/POTT
         
C         KINT = 0.0
C         DO I=1,NAUX
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(1.5*R200))THEN
C         KINT = KIN(I) + KINT
C         ENDIF
C         ENDDO

C         POTT = 0.0
C         DO I=1,NAUX
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(1.5*R200))THEN
C         POTT = POT(I) + POTT
C         ENDIF
C         ENDDO
C         V3 = KINT/POTT

C         KINT = 0.0
C         DO I=1,NAUX
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(2.0*R200))THEN
C         KINT = KIN(I) + KINT
C         ENDIF
C         ENDDO

C         POTT = 0.0
C         DO I=1,N200
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(2.0*R200))THEN
C         POTT = POT(I) + POTT
C         ENDIF
C         ENDDO
C         V4 = KINT/POTT

C         KINT = 0.0
C         DO I=1,NAUX
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(2.5*R200))THEN
C         KINT = KIN(I) + KINT
C         ENDIF
C         ENDDO

C         POTT = 0.0
C         DO I=1,N200
C         IF((SQRT(Y(I)**2+Z(I)**2)).LE.(2.5*R200))THEN
C         POTT = POT(I) + POTT
C         ENDIF
c         ENDDO
C         V5 = KINT/POTT
         
C         WRITE(9,*) J,N200,R200,M200,VDISP,V1,V2,V3,V4,V5

          PRINT*,"DADOS EM FORT.9"
          PRINT*, "NUMERO DO MODELO, NUMERO TOTAL DE PARTICULAS,"
          PRINT*, "N200,R200,M200"
          WRITE(9,*) J,NAUX,N200,R200,M200         

         ENDDO


         END PROGRAM

C INICIO DAS SUBROTINAS E FUNCOES ANEXAS AO PROGRAMA
C=================================================================================
        
      FUNCTION gasdev(isema,isemb)
      INTEGER isema,isemb
      double precision gasdev
C     USES ran3
      INTEGER iset
      double precision fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran1(isema)-1.
        v2=2.*ran1(isemb)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 'Numerical Recipes Software 41$$'!)L.
        
C========================================================================

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 41$$'!)L.

C========================================================================      

        
      SUBROUTINE rtbiss(x1,x2,xacc,i,q1,rtbis)
      INTEGER JMAX
      DOUBLE PRECISION rtbis,x1,x2,xacc
      real ran3
      external ran3
      DOUBLE PRECISION  r
      DOUBLE PRECISION  func
      external func
      DOUBLE PRECISION q1,dens
      INTEGER i
      PARAMETER (JMAX=40)
      INTEGER j
      DOUBLE PRECISION dx,f,fmid,xmid


      fmid=func(x2,q1)
      f=func(x1,q1)


      if(f*fmid.ge.0.)then
      print*, i,f,fmid,q1
       print*, "root must be bracketed in rtbis"
	stop
        endif
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        raiz = rtbis
        fmid=func(xmid,q1)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      print*, "too many bisections in rtbis"
      stop
      END

C  (C) Copr. 1986-92 'Numerical Recipes Software 41$$'!)L.

        
C=======================================================================        
        
        FUNCTION func(xx,q1)
        double precision func,xx,q1,rs
        func = (1. + xx)*(log(1.+xx) - q1) - xx
        return
        end

C========================================================================

        FUNCTION RATIO(YY)
        DOUBLE PRECISION RATIO,YY,A1,A2
        A1 = (LOG(1. + 2.5*YY)-2.5*YY/(1.+2.5*YY))        
        A2 = (LOG(1. + YY)-YY/(1.+YY))
        RATIO = A1/A2
        RETURN
        END
        
C=============================================================================

        FUNCTION VC2(CC,R200,V200,G,RR)
        DOUBLE PRECISION CC,R200,V200,G,RR,VC2
        DOUBLE PRECISION A1,A2
        

        X1 = RR/R200
                
        A1 = (LOG(1. + CC*X1)-CC*X1/(1.+CC*X1))
        A2 = (LOG(1. + CC)-CC/(1.+CC))        
        VC2 = V200*SQRT((1./X1)*(A1/A2))
        RETURN
        END

C===============================================================================        
                
         SUBROUTINE ROT(a,b,c,vx,vy,vz,vvx,vvy,vvz)
         IMPLICIT NONE
         double precision mr(3,3)
         double precision a,b,c
	 double precision vx,vy,vz,vvx,vvy,vvz


         mr(1,1)= cos(a)*cos(b)*cos(c) - sin(a)*sin(c)
         mr(1,2)= sin(a)*cos(b)*cos(c)+cos(a)*sin(c)
         mr(1,3)= -sin(b)*cos(c)
         mr(2,1)= -cos(a)*cos(b)*sin(c) - cos(c)*sin(a)
         mr(2,2)= -sin(a)*cos(b)*sin(c) + cos(a)*cos(c)
         mr(2,3)= sin(b)*sin(c) 
         mr(3,1)= sin(b)*cos(a)
         mr(3,2)= sin(a)*sin(b)
         mr(3,3)= cos(b)

         vvx = mr(1,1)*vx + mr(1,2)*vy + mr(1,3)*vz         
         vvy = mr(2,1)*vx + mr(2,2)*vy + mr(2,3)*vz 
         vvz = mr(3,1)*vx + mr(3,2)*vy + mr(3,3)*vz 
        
         return
         end

