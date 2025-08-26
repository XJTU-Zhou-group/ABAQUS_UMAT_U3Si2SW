C********************************TOP***************************
C     SDVINI & UMAT FOR ABAQUS/STANDARD
C     MODIFIED FG RATE THEORY-BASED MODEL FOR AMORPHOUS U3SI2 SWELLING
C     C3D8
C     THIS CODE WAS DEVELOPED BY XIAOWEI YUE
C*******************************************************************
      
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC' 
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
C
      STATEV(1)=0.D0
      STATEV(2)=0.D0
      STATEV(3)=0.D0       
      STATEV(4)=1.D22
      STATEV(5)=1.D19
      STATEV(6)=1.D-1
      STATEV(7)=1.D-10
      STATEV(8)=1.D-17
      
!      WRITE(6,*)'SDVINI-STATEV',(STATEV(IXX),IXX=1,6)
      RETURN
      END

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
      
      INCLUDE 'ABA_PARAM.INC'
      
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      
      PARAMETER (M=3,N=3,ID=3,ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,
     1 FOUR=4.D0, FIVE=5.D0, SIX=6.D0, EIGHT=8.D0, NINE=9.D0,
     2 TOLER=1.D-5,TOLER2=100, PAI=3.141592654D0,BOCONS=1.380649D-23)

      DIMENSION DEPSW(6),EPSWOLD(6),DESTRAN(6),DEPTH(6),DSTRESS(6)
      INTEGER NProny,counter
      REAL*8 Radt,SW,THEXPL,EG,ELAM,E,XNUE,XTHEXPL,XTEMP,SURT,RP,RK,XN1,
     1 XN2,SWGBF,SWGAF,XXA,XXAKNEE,XMZE,BETA,FIR,PHI,VDWCONS,BZE,DG,FN,
     2 DZE,RG,ZETA,ZETAZE,ALPHA,XMRP,XMRT,RT,RI,XIRI,BETAE,BETAS,SWSD,
     3 FIRZE,BRS,XPHI1,XPHI2,XPHI3,RES1,RES2,RES3,CGI,CBI,XMBI,RBI,DBI,
     4 CGIT,CBIT,XMBIT,RBIT,DBIT,XPHI1PCG,XPHI2PCB,XPHI3PMB,DXDCG,DXDCB,
     5 DXDMB,XDCG,XDCB,XDMB,PHIR,SIGH,PHIRPRB,DRB,DTL,STEPNUM,XLSW

      
      E = PROPS(1)
      XNUE = PROPS(2)
      XTHEXPL = PROPS(3)
      FIR = PROPS(4)
      FIRZE = 1.25D20
      XTEMP = 373
      VDWCONS = (5.16D-5)/(6.023D23)
      BZE = 2.D-23
      BETA = 0.25
      RG = 0.216D-9
      SURT = 0.7
      BETAS = 1.75D-29
      
      ZETAZE = 2.D5*(FIRZE)
      ALPHA = 5.D-9*(FIRZE)
      ZETA=ZETAZE/FIR
      BRS=BZE*FIR
      DZE=(BOCONS*XTEMP)/(SIX*PAI*RG*ZETAZE)
      DG=(BOCONS*XTEMP)/(SIX*PAI*RG*ZETA)
      FN=FOUR*ALPHA*ZETAZE*RG
      SIGH=-(ONE/THREE)*(STRESS(1)+STRESS(2)+STRESS(3))*1.D6
!      SIGH=1000.D6
!      SIGH=0.D0
      
      
      IF(SIGH.LE.0.)THEN
          SIGH=0.
      END IF
      
      
      EG = E/(TWO*(ONE+XNUE))
      ELAM =(E*XNUE)/((ONE+XNUE)*(ONE-TWO*XNUE))

      
      DO I = 1,6
          DO J = 1,6
              DDSDDE(I,J) = 0.D0
          END DO
      END DO
      
      
      DO K1 = 1,3
          DO K2 = 1,3
              DDSDDE(K2,K1) = ELAM
          END DO
          DDSDDE(K1,K1) = TWO*EG+ELAM
      END DO
      
      DDSDDE(4,4)=EG
      DDSDDE(5,5)=EG
      DDSDDE(6,6)=EG
      
      DO K=1,3
          EPSWOLD(K)=STATEV(K)
      END DO
      
      DO K=4,6
          EPSWOLD(K)=ZERO
      END DO
      
      
      Radt=TIME(2)+DTIME
      
      CGI=STATEV(4)
      CBI=STATEV(5)
      XMBI=STATEV(6)
      RBI=STATEV(7)
      DBI=STATEV(8)
            
      
      STEPNUM=1.D3
      DTL=DTIME/STEPNUM
   
      
      DO KL=1,1000
      
      CGIT=CGI
      CBIT=CBI
      XMBIT=XMBI
      RBIT=RBI
      DBIT=DBI
      XDCG=ZERO
      XDCB=ZERO
      XDMB=ZERO
      
      
      DO KNEW=1,100
          XPHI1=BETA*FIR-(FOUR*PAI*FN*DG/((FIR)**TWO))*CGI*CGI
     1    -FOUR*PAI*DG*RBI*CBI*CGI+TWO*BRS*XMBI*CBI
          
          XPHI2=(FOUR*PAI*FN*DG/((FIR)**TWO))*CGI*CGI/XMBI-BRS*CBI
     1    -(EIGHT*TWO)*PAI*RBI*CBI*CBI*DBI
          
          XPHI3=FOUR*PAI*DG*RBI*CGI-BRS*XMBI
     1    +(EIGHT*TWO)*PAI*RBI*XMBI*CBI*DBI
          
          RES1=(XDCG/DTL)-XPHI1
          RES2=(XDCB/DTL)-XPHI2
          RES3=(XDMB/DTL)-XPHI3
          
          IF ((ABS(RES1/XPHI1).LT.TOLER).AND.
     1    (ABS(RES2/XPHI2).LT.TOLER).AND.
     2    (ABS(RES3/XPHI3).LT.TOLER)) GOTO 10
          
          XPHI1PCG=-(FOUR*PAI*FN*DG/((FIR)**TWO))*TWO*CGI
     1    -FOUR*PAI*DG*RBI*CBI
          
          XPHI2PCB=-BRS-(EIGHT*TWO)*PAI*RBI*DBI*TWO*CBI
          
          XPHI3PMB=-BRS+(EIGHT*TWO)*PAI*RBI*CBI*DBI
          
          DXDCG=-RES1/((ONE/DTL)-XPHI1PCG)
          DXDCB=-RES2/((ONE/DTL)-XPHI2PCB)
          DXDMB=-RES3/((ONE/DTL)-XPHI3PMB)
          
          XDCG=XDCG+DXDCG
          XDCB=XDCB+DXDCB
          XDMB=XDMB+DXDMB
          
          CGI=CGIT+XDCG
          CBI=CBIT+XDCB
          XMBI=XMBIT+XDMB
          
          DO KNEWNEW=1,10
              PHIR=(EIGHT/THREE)*PAI*SURT*(RBI**TWO)
     1        -TWO*SURT*VDWCONS*XMBI/RBI
     2        +(FOUR/THREE)*PAI*(RBI**THREE)*SIGH-VDWCONS*XMBI*SIGH
     3        -XMBI*BOCONS*XTEMP
              IF (ABS(PHIR).LT.1.D-30) GOTO 20
              PHIRPRB=(EIGHT*TWO/THREE)*PAI*SURT*RBI
     1        +TWO*SURT*VDWCONS*XMBI*(RBI**(-TWO))
     2        +FOUR*PAI*(RBI**TWO)*SIGH
              DRB=-PHIR/PHIRPRB
              RBI=RBI+DRB
          END DO
20        CONTINUE

          DBI=((RG/RBI)**THREE)*DG
          
          
      END DO
10    CONTINUE      
      
      
      END DO
      
      SWG=((FOUR/THREE)*PAI*(RBI**THREE))*CBI
      
      SWS=BETAS*FIR*Radt
      
      STATEV(10)=SWG
      STATEV(11)=SWS

      
      STATEV(12)=RES1
      STATEV(13)=RES2
      STATEV(14)=RES3
      STATEV(15)=KNEW
      STATEV(16)=KNEWNEW
      STATEV(17)=PHIR
      STATEV(18)=SIGH

      
      SW = SWG+SWS
      STATEV(9)=SW
      
      DO K = 1,6
          DEPSW(K) = 0.D0
          DEPTH(K) = 0.D0
          DSTRESS(K) = 0.D0
      END DO
  
      XLSW=(ONE/THREE)*LOG(ONE+SW)
      
!      DEPSW=DEPSW*ZERO
      DO K=1,3
          DEPSW(K)=XLSW-EPSWOLD(K)
      END DO
      
!      DEPTH=DEPTH*ZERO
      DO K=1,3
          DEPTH(K)=LOG(ONE+XTHEXPL*(TEMP+DTEMP))-LOG(ONE+XTHEXPL*(TEMP))
      END DO
     
      STATEV(19)=DEPSW(1)
      STATEV(20)=DEPSW(2)
      STATEV(21)=DEPSW(3)
      STATEV(22)=DEPSW(4)
      STATEV(23)=DEPSW(5)
      STATEV(24)=DEPSW(6)
      STATEV(25)=DEPTH(1)
      STATEV(26)=DEPTH(2)
      STATEV(27)=DEPTH(3)
      STATEV(28)=DEPTH(4)
      STATEV(29)=DEPTH(5)
      STATEV(30)=DEPTH(6)
      
      DO K=1,6
          DESTRAN(K)=DSTRAN(K)-DEPSW(K)-DEPTH(K)
      END DO
      
!      DO K=1,6
!      DESTRAN(K)=DSTRAN(K)
!      END DO
      
!      DSTRESS=DSTRESS*ZERO
      DO K1=1,NTENS
          DO K2=1,NTENS
              DSTRESS(K1)=DSTRESS(K1)+DDSDDE(K1,K2)*DESTRAN(K2)
          END DO
      END DO
      
      DO K=1,NTENS
          STRESS(K)=STRESS(K)+DSTRESS(K)
      END DO
      
      DO K=1,3
          STATEV(K)=EPSWOLD(K)+DEPSW(K)
      END DO      
      
      STATEV(4)=CGI
      STATEV(5)=CBI
      STATEV(6)=XMBI
      STATEV(7)=RBI
      STATEV(8)=DBI

     
      
      RETURN
      END
