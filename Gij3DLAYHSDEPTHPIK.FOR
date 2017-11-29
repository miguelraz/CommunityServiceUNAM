                                                                                                                                                                                                                     C
C     Programa Gij3DLAYHSDEPTHMIK 
C
C     Programa para calcular sismogramas sintéticos de Gij(XYZ,XfYfZf,t); 
C     Fuente y receptores en un semiespacio (HS).
C
C     La fuentes son una fuerza vertical y una horizontal en la interfaz
C     superior del estrato j (JFUEN = 1, o 2, ... o N, o N+1).
C
C     Por ahora los receptores están en superficie.  
C
C     A partir de la sol en coordenadas polares, Se hace integración en el 
C     número de onda radial con Jo y J1,Calculo con frecuencia imaginaria
C     Se calculan sismogramas sinteticos en hasta 25 estaciones.
C
C     Componentes del tensor de Green en coordenadas polares...
C 
C
      PARAMETER(NES=16)
      COMMON/ONE/B
      COMMON/MAT/A
      COMMON/ONEH/BH
      COMMON/MATH/AH

      COMMON/SINT/DT,TS,TP
	COMMON/ENU/NU, ENUD, ENUU
	COMMON/EGAM/GAM, EGAMD, EGAMU

      REAL K, L, AUX
      COMPLEX Gxx(25),Gyx(25),Gzx(25)
      COMPLEX Gxy(25),Gyy(25),Gzy(25)
      COMPLEX Gxz(25),Gyz(25),Gzz(25)
      COMPLEX INTEG,INTEGR,INTEGT
      COMPLEX COME, KPC, KSC, KPCHS, KSCHS
      COMPLEX NU, GAM, EGAMH, ENUH, NUHS, GAMHS, ETA
      COMPLEX UI, UR, ZER, C1, C2,  PIK
      COMPLEX A(4*NES+2,4*NES+2), B(4*NES+2,2)
      COMPLEX AH(2*NES+1,2*NES+1), BH(2*NES+2,1)

      COMPLEX SD(2,2),SU(2,2),UD(2,2),UU(2,2),  SHD,SHU,VD,VU
	REAL BETA(NES+1),ALFA(NES+1),AMU(NES+1),Z(NES+1), BBAL,ABAL
	REAL QDAMP(NES+1),QDAMS(NES+1)
	COMPLEX Vxx(512,25), Vyx(512,25), Vzx(512,25)
	COMPLEX Vxy(512,25), Vyy(512,25), Vzy(512,25)
      COMPLEX Vxz(512,25), Vyz(512,25), Vzz(512,25)
	COMPLEX CY(512),CR(512)                   
      COMPLEX ARG,H0,H1,CJ(101)

      DIMENSION TV(512),VT(512)                                                    
C
      OPEN(1,FILE="Gij3DLAYHS.DAT",FORM="FORMATTED")
C
      OPEN(11,FILE="SINTSxx.RES",FORM="FORMATTED") 
      OPEN(21,FILE="SINTSyx.RES",FORM="FORMATTED")
      OPEN(31,FILE="SINTSzx.RES",FORM="FORMATTED")
      OPEN(12,FILE="SINTSxy.RES",FORM="FORMATTED") 
      OPEN(22,FILE="SINTSyy.RES",FORM="FORMATTED")
      OPEN(32,FILE="SINTSzy.RES",FORM="FORMATTED")
      OPEN(13,FILE="SINTSxz.RES",FORM="FORMATTED")
      OPEN(23,FILE="SINTSyz.RES",FORM="FORMATTED")
      OPEN(33,FILE="SINTSzz.RES",FORM="FORMATTED")
C
      OPEN(30,FILE="MATRIZ.RES",FORM="FORMATTED")
C
C     COMIENZA DEFINICION DE CONSTANTES Y LECTURA DE PARAMETROS 
C

      UI =(0.0,1.0)
      UR =(1.0,0.0)
	ZER=(0.0,0.0)
      PI=4.0*ATAN(1.0)
C
      READ(1,*)
	READ(1,*)NL   !NUMBER OF LAYERS (<17) HALF-SPACE AT THE END
      READ(1,*)
C	PRINT*,"NL = ",NL
	Z(1)=0.0
C      PRINT*,"  H,  ALF,   BET,    RHO,   QP,   QS "
	DO J=1,NL
	   READ(1,*) H, ALF, BET, RHO, QP, QS
C        PRINT*,   H, ALF, BET, RHO, QP, QS
         Z(J+1)=Z(J)+H
         AMU(J)=RHO*BET**2
         BETA(J)=BET
         ALFA(J)=ALF
         QDAMP(J)=QP
         QDAMS(J)=QS
C        BEALF=SQRT((0.5-ANU)/(1.0-ANU)) !IF POISSON RATIO IS GIVEN
C        ALFA(J)=BET/BEALF
      END DO
      READ(1,*) H, ALF, BET, RHO, QP, QS
C      PRINT*, H, ALF, BET, RHO, QP, QS
      AMU(NL+1)=RHO*BET**2
      BETA(NL+1)=BET
      ALFA(NL+1)=ALF
	QDAMP(NL+1)=QP
	QDAMS(NL+1)=QS
C      BEALF=SQRT((0.5-ANU)/(1.0-ANU))   !IF POISSON RATIO IS GIVEN
C      ALFA(NL+1)=BET/BEALF
      READ(1,*)
	READ(1,*)L, FACTW, FACNMAX, JFUEN
C	PRINT*,"   L, FACTW, FACNMAX, JFUEN "
C     PRINT*,    L, FACTW, FACNMAX, JFUEN 
	READ(1,*)
	READ(1,*)XMIN, XMAX, YMIN, YMAX, ZRCPT, NEST
C      PRINT*,"  XMIN,XMAX,YMIN,YMAX,ZRCPT,NEST"
C      PRINT*,   XMIN,XMAX,YMIN,YMAX,ZRCPT,NEST
	READ(1,*)
	READ(1,*)TP, TS, DT, NT
C     PRINT*, "  TP, TS, DT, NT "
C     PRINT*,    TP, TS, DT, NT 
C
	NFREC=NT/2 
      DX=(XMAX-XMIN)/(NEST-1)
      DY=(YMAX-YMIN)/(NEST-1)	
	DUR=NT*DT
	TW=DUR * FACTW
	DFREC=1.0/DUR
C
C 	PRINT*,"  DX,DY,DUR,TW,DFREC"
C     PRINT*,   DX,DY,DUR,TW,DFREC              
C
C     COMIENZA CALCULO
C
      DO J=1,NFREC 
         FREC=DFREC*J
         OME=2.0*PI*FREC
	   OMEI=0.7*PI/TW
         COME=CMPLX(OME, -OMEI)    !AMORTIGUAMIENTO FUENTES PERIODICAS
C         COME=COME*(UR-UI/2.0/Q)  !AMORTIGUAMIENTO HISTERETICO
C	    KPC=COME/ALFA  ;  KPCHS=COME/ALFAHS 
C         KSC=COME/BETA  ;  KSCHS=COME/BETAHS
         DK=2.0*PI/L
         AUX=OME/BETA(1)
         NMAX=AUX/DK
	   NMAX=NMAX*FACNMAX+50    
C
	   PRINT*,' J = ',J,'  NMAX =',NMAX, ' SONRIE'
C
         DO IX=1,NEST
            Gxx(IX)=ZER; Gyx(IX)=ZER; Gzx(IX)=ZER
            Gxy(IX)=ZER; Gyy(IX)=ZER; Gzy(IX)=ZER
            Gxz(IX)=ZER; Gyz(IX)=ZER; Gzz(IX)=ZER
         END DO
C
C     WAVENUMBER LOOP
C
         DO IK= 1, NMAX
            K=0.5*DK*(2*IK-1)
            PIK=+UI*K
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
C
C           COEFFICIENTS MATRIX ---CASE P-SV
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      
	      DO ILIN=1,4*NL+2
               B(ILIN,1)=ZER
               B(ILIN,2)=ZER
               DO ICOL=1,4*NL+2
	            A(ILIN,ICOL)=ZER
               END DO
            END DO
C
c      PRINT*,'  EMPIEZO A ARMAR MATRIZ DE COEFS  P-SV '
            ZCAL=Z(1)
            BBAL=BETA(1);ABAL=ALFA(1);QBALP=QDAMP(1);QBALS=QDAMS(1) 
      CALL DESF(COME,K,BBAL,ABAL,QBALP,QBALS,ZCAL,Z,1,NL,UD,UU,SD,SU)
            A(1,1)=SD(1,1);  A(1,2)=SD(1,2)
		  A(1,3)=SU(1,1);  A(1,4)=SU(1,2)
C
C            PRINT*,'A11  A12  A13  A14', A(1,1),A(1,2),A(1,3),A(1,4) 
C            
            A(2,1)=SD(2,1);  A(2,2)=SD(2,2)
		  A(2,3)=SU(2,1);  A(2,4)=SU(2,2)
C
C            PRINT*,'A21  A22  A23  A24', A(2,1),A(2,2),A(2,3),A(2,4) 
C
	      DO IES=1,NL-1
               ZCAL=Z(IES+1)   
               BBAL=BETA(IES); ABAL=ALFA(IES)
               QBALP=QDAMP(IES); QBALS=QDAMS(IES)
      CALL DESF(COME,K,BBAL,ABAL,QBALP,QBALS,ZCAL,Z,IES,NL,UD,UU,SD,SU)
C               
               IR=4*(IES-1)+3; JR=4*(IES-1)+1  
C               
               A(IR+0,JR+0)=UD(1,1); A(IR+0,JR+1)=UD(1,2)
               A(IR+1,JR+0)=UD(2,1); A(IR+1,JR+1)=UD(2,2)
C                
               A(IR+0,JR+2)=UU(1,1); A(IR+0,JR+3)=UU(1,2)
               A(IR+1,JR+2)=UU(2,1); A(IR+1,JR+3)=UU(2,2)
C               
               A(IR+2,JR+0)=SD(1,1); A(IR+2,JR+1)=SD(1,2)
               A(IR+3,JR+0)=SD(2,1); A(IR+3,JR+1)=SD(2,2)
C               
               A(IR+2,JR+2)=SU(1,1); A(IR+2,JR+3)=SU(1,2)
               A(IR+3,JR+2)=SU(2,1); A(IR+3,JR+3)=SU(2,2)
C               
               BBAL=BETA(IES+1); ABAL=ALFA(IES+1)
               QBALP=QDAMP(IES+1); QBALS=QDAMS(IES+1)
C
      CALL DESF(COME,K,BBAL,ABAL,QBALP,QBALS,
     .                                     ZCAL,Z,IES+1,NL,UD,UU,SD,SU)
C                
               IR=4*(IES-1)+3; JR=4*(IES-1)+5 ; RAT=AMU(IES+1)/AMU(IES)
C               
               A(IR+0,JR+0)=-UD(1,1); A(IR+0,JR+1)=-UD(1,2)
               A(IR+1,JR+0)=-UD(2,1); A(IR+1,JR+1)=-UD(2,2)
C               
               A(IR+0,JR+2)=-UU(1,1); A(IR+0,JR+3)=-UU(1,2)
               A(IR+1,JR+2)=-UU(2,1); A(IR+1,JR+3)=-UU(2,2)
C               
               A(IR+2,JR+0)=-RAT*SD(1,1); A(IR+2,JR+1)=-RAT*SD(1,2)
               A(IR+3,JR+0)=-RAT*SD(2,1); A(IR+3,JR+1)=-RAT*SD(2,2)
C              
               A(IR+2,JR+2)=-RAT*SU(1,1); A(IR+2,JR+3)=-RAT*SU(1,2)
               A(IR+3,JR+2)=-RAT*SU(2,1); A(IR+3,JR+3)=-RAT*SU(2,2)
C
            END DO
            IES=NL
            ZCAL=Z(IES+1)
            BBAL=BETA(IES); ABAL=ALFA(IES)
            QBALP=QDAMP(IES); QBALS=QDAMS(IES)
C
      CALL DESF(COME,K,BBAL,ABAL,QBALP,QBALS,ZCAL,Z,IES,NL,UD,UU,SD,SU)
C               
            IR=4*(IES-1)+3; JR=4*(IES-1)+1
C
            A(IR+0,JR+0)=UD(1,1); A(IR+0,JR+1)=UD(1,2)
            A(IR+1,JR+0)=UD(2,1); A(IR+1,JR+1)=UD(2,2)
C                
            A(IR+0,JR+2)=UU(1,1); A(IR+0,JR+3)=UU(1,2)
            A(IR+1,JR+2)=UU(2,1); A(IR+1,JR+3)=UU(2,2)
C               
            A(IR+2,JR+0)=SD(1,1); A(IR+2,JR+1)=SD(1,2)
            A(IR+3,JR+0)=SD(2,1); A(IR+3,JR+1)=SD(2,2)
C               
            A(IR+2,JR+2)=SU(1,1); A(IR+2,JR+3)=SU(1,2)
            A(IR+3,JR+2)=SU(2,1); A(IR+3,JR+3)=SU(2,2)
C               
            BBAL=BETA(IES+1); ABAL=ALFA(IES+1)
            QBALP=QDAMP(IES+1); QBALS=QDAMS(IES+1)
C
      CALL DESF(COME,K,BBAL,ABAL,QBALP,QBALS,ZCAL,Z,
     .                                         IES+1,NL,UD,UU,SD,SU)
C                
            IR=4*(IES-1)+3; JR=4*(IES-1)+5 ; RAT=AMU(IES+1)/AMU(IES)
C               
            A(IR+0,JR+0)=-UD(1,1); A(IR+0,JR+1)=-UD(1,2)
            A(IR+1,JR+0)=-UD(2,1); A(IR+1,JR+1)=-UD(2,2)
C                           
            A(IR+2,JR+0)=-RAT*SD(1,1); A(IR+2,JR+1)=-RAT*SD(1,2)
            A(IR+3,JR+0)=-RAT*SD(2,1); A(IR+3,JR+1)=-RAT*SD(2,2)
C
c      PRINT*,'  YA ARME LA MATRIZ P-SV'            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C             
C           INDEPENDIENT TERMS --- CASE P-SV
C            
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
            IF(JFUEN.EQ.1)THEN
              SIGN=1.0
	        B(1,1)= -SIGN*UR/2./PI/AMU(1)*K *1.0 ! FZA VERT EN LA CIMA ESTRATO 1
	        B(2,2)= +SIGN*UI/2./PI/AMU(1)*K *1.0 ! FZA HOR  EN LA CIMA ESTRATO 1
              ELSE 
              SIGN=-1.0 
	        B(4*(JFUEN-1)+1,1)= -SIGN*UR/2./PI/AMU(JFUEN-1)*K*1.0   ! FZA VERT EN LA CIMA ESTRATO JFUEN
	        B(4*(JFUEN-1)+2,2)= +SIGN*UI/2./PI/AMU(JFUEN-1)*K*1.0   ! FZA HOR  EN LA CIMA ESTRATO JFUEN
            END IF
C           
            CALL IMP(A, 4*NL+2, 4*NL+2)
            CALL IMP(B, 4*NL+2, 2) 

C            CALL CATA(A,B,22,2)
C
C      PRINT*,'  YA VOY A SOLUCO P-SV'    
C
             CALL SOLUCO(4*NL+2,2)
C
C            CALL IMP(A, 22,22)
C            CALL IMP(B, 22,2) 
c      PRINT*,'  SALGO DE SOLUCO '
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC                  
C
C           COEFICIENTS MATRIX --- CASE SH
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
c      PRINT*,'  EMPIEZO A ARMAR MATRIZ DE COEFS SH '     
	      DO ILIN=1,2*NL+1
               BH(ILIN,1)=ZER
               DO ICOL=1,2*NL+1
	            AH(ILIN,ICOL)=ZER
               END DO
            END DO
C
            ZCAL=Z(1)
            BBAL=BETA(1); QBAL=QDAMS(1)
           CALL DESFH(COME,K,BBAL,QBAL,ZCAL,Z,1,NL,VD,VU,SHD,SHU)
            AH(1,1)=SHD;  AH(1,2)=SHU
C
	      DO IES=1,NL-1
               ZCAL=Z(IES+1)   
               BBAL=BETA(IES); QBAL=QDAMS(IES+1)
           CALL DESFH(COME,K,BBAL,QBAL,ZCAL,Z,IES,NL,VD,VU,SHD,SHU)
C               
               IR=2*(IES-1)+2; JR=2*(IES-1)+1 
C			 PRINT*, IR, JR 
C               
               AH(IR+0,JR+0)=VD
               AH(IR+0,JR+1)=VU
               AH(IR+1,JR+0)=SHD
               AH(IR+1,JR+1)=SHU
C               
               BBAL=BETA(IES+1); QBAL=QDAMS(IES+1)
           CALL DESFH(COME,K,BBAL,QBAL,ZCAL,Z,IES+1,NL,VD,VU,SHD,SHU)
C                
               IR=2*(IES-1)+2; JR=2*(IES-1)+3 ; RAT=AMU(IES+1)/AMU(IES)
C               
               AH(IR+0,JR+0)=-VD
               AH(IR+0,JR+1)=-VU
               AH(IR+1,JR+0)=-RAT*SHD
               AH(IR+1,JR+1)=-RAT*SHU
C
            END DO
            IES=NL
            ZCAL=Z(IES+1)
            BBAL=BETA(IES); QBAL=QDAMS(IES)
            CALL DESFH(COME,K,BBAL,QBAL,ZCAL,Z,IES,NL,VD,VU,SHD,SHU)
C               
            IR=2*(IES-1)+2; JR=2*(IES-1)+1
C            PRINT*, IR, JR
C
            AH(IR+0,JR+0)=VD
            AH(IR+0,JR+1)=VU
            AH(IR+1,JR+0)=SHD
            AH(IR+1,JR+1)=SHU
C               
            BBAL=BETA(IES+1); QBAL=QDAMS(IES+1)
            CALL DESFH(COME,K,BBAL,QBAL,ZCAL,Z,IES+1,NL,VD,VU,SHD,SHU)
C                
            IR=2*(IES-1)+2; JR=2*(IES-1)+3 ; RAT=AMU(IES+1)/AMU(IES)
C            PRINT*, IR, JR
C               
            AH(IR+0,JR+0)=-VD
            AH(IR+1,JR+0)=-RAT*SHD
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C             
C           INDEPENDIENT TERMS - CASE SH
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C   
            IF(JFUEN.EQ.1)THEN
              SIGN=1.0
	        BH(1,1)= +SIGN*UR/2./PI/AMU(1)*1.0  ! FZA HOR (SH) EN LA CIMA ESTRATO 1
              ELSE 
              SIGN=-1.0 
	        BH(2*JFUEN-1,1)=+SIGN*UR/2./PI/AMU(JFUEN-1)*1.0  ! FZA HOR EN LA CIMA ESTRATO JFUEN
            END IF
C
C            CALL IMPH(AH,2*NL+1,2*NL+1)  
c      PRINT*,'  YA VOY A SOLUCO SH'         
             CALL SOLUCOH(2*NL+1,1)
C
c      PRINT*,'  SALGO DE SOLUCOH '
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C           COEFICIENTES PARTE SH
C 
c            C2=-UI/ENUH/2.0/PI/AMU/NU/(UR/ENUH-ETA*ENUH)      !CdCOS(NU*H)
c            C1=-UI *ETA/2.0/PI/AMU/NU/(UR/ENUH-ETA*ENUH)      !CdCOS(NU*H)
C
C            C2=-UI/2.0/PI/AMU(1)/NU
C            C1=ZER  
c            PRINT*, 'NU = ', NU, ' C2 = ',C2
C
C           POSICIONES DE LOS RECEPTORES

            DO IX=1,NEST
               X=XMIN+DX*(IX-1)
               Y=YMIN+DY*(IX-1)
  		     R=SQRT(X*X+Y*Y)
	         CTETA=X/R
               STETA=Y/R
               ARG=CMPLX(K*R,0.0)    !ARGUMENTO DE LAS FUNCS DE BESSEL
C              CALL BESJC(ARG,2,CJ)  !CALCULO DE LAS FUNCS DE BESSEL   
c      PRINT*,'  VOY A CALCULAR H0 y H1' 
c	PRINT*,'  ARG = ',ARG, '  K = ', K, ' R = ', R     
	         CALL HANKEL(ARG,H0,H1) !CALCULO DE LAS FUNCS DE HANKEL --> BESSEL
c	PRINT*,'  ARG = ',ARG, ' H0 = ', H0,'  H1 = ', H1   
C
            ZCAL=Z(1)
            BBAL=BETA(1);ABAL=ALFA(1);QBALP=QDAMP(1);QBALS=QDAMS(1)
c      PRINT*,'  VOY A DESF P-SV'
      CALL  DESF(COME,K,BBAL,ABAL,QBALP,QBALS,ZCAL,Z,1,NL,UD,UU,SD,SU)
c      PRINT*,'  VOY A DESF SH'
      CALL DESFH(COME,K,BBAL,QBALS,ZCAL,Z,1,NL,VD,VU,SHD,SHU)
C           
C              CARGA VERTICAL; RECEPTORES EN SUPERFICIE  
C            
C              Ur:
C
               INTEG=UR*( -K *B(1,1)- K *B(3,1)*EGAMU                        ! P
     .           +PIK* NU*B(2,1)-PIK* NU*B(4,1)*ENUU )*REAL(H1) !CJ(2)=J1(kr)!SV
               Gxz(IX)=Gxz(IX)+INTEG*DK*CTETA
               Gyz(IX)=Gyz(IX)+INTEG*DK*STETA

C              Uz:

               INTEG=UI*(-GAM*B(1,1)+GAM*B(3,1)*EGAMU                     !P
     .          -PIK*K *B(2,1)-PIK*K *B(4,1)*ENUU)*REAL(H0)               !SV                                                                                      !CJ(1)=J0(kr)  !SV
               Gzz(IX)=Gzz(IX)+INTEG*DK
		  
C              CARGA HORIZONTAL RECEPTORES EN SUPERFICIE

C              Ur:
C
           INTEGR=+UR*(+K*B(1,2)*1+ K*B(3,2)*1*EGAMU-PIK*NU          !P
     .     *B(2,2)*1+PIK*NU*B(4,2)*1*ENUU)*(REAL(H0)-REAL(H1)/ARG)   !SV
     .     +UR*(VD*BH(1,1)*1+VU*BH(2,1)*1)*K*(REAL(H1)/ARG)          !SH
C    .         +UR*( C2+ C1*ENUU )        *K*(REAL(H1)/ARG)          !SH (SINGLE)
C   
C              Ut:              ^ OJO
C
          INTEGT=-UR*(+K*B(1,2)*1+ K*B(3,2)*1*EGAMU-PIK*NU          !P
     .    *B(2,2)*1+PIK*NU*B(4,2)*1*ENUU)*(REAL(H1)/ARG)            !SV
     .    -UR*(VD*BH(1,1)*1+VU*BH(2,1)*1)*K*(REAL(H0)-REAL(H1)/ARG)  !SH
C    .    -UR*( C2+ C1*ENUU)        *K*(REAL(H0)-REAL(H1)/ARG)  !SH (SINGLE)
C            
C            PRINT*, 'BH(1,1) = ', BH(1,1)
C
               Gxx(IX)=Gxx(IX)+(INTEGR*CTETA**2-INTEGT*STETA**2)*DK                       
               Gyx(IX)=Gyx(IX)+(INTEGR+INTEGT)*STETA*CTETA*DK 

               Gxy(IX)=Gxy(IX)+(INTEGR+INTEGT)*STETA*CTETA*DK                       
               Gyy(IX)=Gyy(IX)+(INTEGR*STETA**2-INTEGT*CTETA**2)*DK  

C              Uz:

               INTEG=+UI*(-GAM*B(1,2)+GAM*B(3,2)*EGAMU              !P
     .         -PIK* K*B(2,2) -PIK* K*B(4,2)*ENUU )*REAL(H1)        !SV      OK
               Gzx(IX)=Gzx(IX)+INTEG*DK * CTETA
               Gzy(IX)=Gzy(IX)+INTEG*DK * STETA

            END DO
	   END DO	
C
C         DO IX=1,NEST
C            X=XMIN+DX*(IX-1)
C            Y=YMIN+DY*(IX-1)
C            R=SQRT(X*X+Y*Y)
C            CTETA=X/R
C            STETA=Y/R
C            ARG=KSC*R
C            Gtx(IX)=CEXP(-UI*ARG)
C         END DO
C
C
C        WRITE(8,100)FREC, Grz, Gzz, Grx, Gtx, Gzx    !OJO
C
C         WRITE(13,100)FREC, (Grz(IX), IX=1,NEST)
C         WRITE(33,100)FREC, (Gzz(IX), IX=1,NEST)  
C         WRITE(11,100)FREC, (Grx(IX), IX=1,NEST)  
C         WRITE(21,100)FREC, (Gtx(IX), IX=1,NEST)  
C         WRITE(31,100)FREC, (Gzx(IX), IX=1,NEST)  
C

         DO IX=1,NEST
            Vxx(J,IX)=Gxx(IX); Vyx(J,IX)=Gyx(IX); Vzx(J,IX)=Gzx(IX) 
            Vxy(J,IX)=Gxy(IX); Vyy(J,IX)=Gyy(IX); Vzy(J,IX)=Gzy(IX) 
            Vxz(J,IX)=Gxz(IX); Vyz(J,IX)=Gyz(IX); Vzz(J,IX)=Gzz(IX) 
         END DO
	END DO
C
C     CALCULO DE SISMOGRAMAS
C     
      DO J=1,NT            
         FUN=RIC(J)*EXP(-OMEI*DT*(J-1))      
	   CR(J)=CMPLX(FUN,0.0) 
      END DO
C
      CALL FORK(NT,CR,-1.0)
C                    
C     CALCULO DE Grz, Gzz, Grx, Gtx, Gzx
C
c

      PRINT*,'  YA VOY A CALCULAR SINTETICOS'
      CALL SINTE(11,NEST,NT,DT,Vxx,OMEI,CR)
      CALL SINTE(21,NEST,NT,DT,Vyx,OMEI,CR)
      CALL SINTE(31,NEST,NT,DT,Vzx,OMEI,CR)
      CALL SINTE(12,NEST,NT,DT,Vxy,OMEI,CR)
      CALL SINTE(22,NEST,NT,DT,Vyy,OMEI,CR)
      CALL SINTE(32,NEST,NT,DT,Vzy,OMEI,CR)
      CALL SINTE(13,NEST,NT,DT,Vxz,OMEI,CR)
      CALL SINTE(23,NEST,NT,DT,Vyz,OMEI,CR)
      CALL SINTE(33,NEST,NT,DT,Vzz,OMEI,CR)
C
  200 FORMAT(2E15.6) 
  100 FORMAT(11e15.6)
      END
C
      SUBROUTINE DESF(COME,K,BETA,ALFA,QP,QS,ZCAL,Z,IES,N,UD,UU,SD,SU)
      PARAMETER(NES=16)
      COMMON/COM/UI,UR,ZER
	COMMON/ENU/NU, ENUD, ENUU
	COMMON/EGAM/GAM, EGAMD, EGAMU
      COMPLEX COME, KPC, KSC, NU, GAM 
      COMPLEX EGAMD, ENUD, EGAMU, ENUU
      COMPLEX SD(2,2),SU(2,2),UD(2,2),UU(2,2)
      COMPLEX UI, UR, ZER, PIK
	REAL K, ALFA, BETA, Z(NES+1)                              
C
      UI =CMPLX(0.0,1.0)
      UR =CMPLX(1.0,0.0)
      ZER=CMPLX(0.0,0.0)
	PIK=+UI*K
C
      PI=4.0*ATAN(1.0)
C
      KPC=COME/ALFA*(UR-UI/2.0/QP)
      KSC=COME/BETA*(UR-UI/2.0/QS )
C
C     COMPUTATION OF GAM AND NU USING COMPLEX FREQUENCY
C
      GAM=CSQRT(KPC*KPC-UR*K*K)
      NU= CSQRT(KSC*KSC-UR*K*K)
      IF(AIMAG(GAM).GT.0.0)GAM=-GAM
      IF(AIMAG(NU) .GT.0.0)NU =-NU
C
      EGAMD=CEXP( -UI*GAM*(ZCAL-Z(IES) ) )
      ENUD= CEXP( -UI* NU*(ZCAL-Z(IES) ) )
C
C     SUB-MATRICES UD,SD
C
      SD(1,1)=-(NU*NU-K*K*UR)*EGAMD
	SD(2,1)=-2.0*K*GAM*EGAMD
	SD(1,2)=-2.0*K*NU*ENUD       *PIK
	SD(2,2)=+(NU*NU-K*K*UR)*ENUD *PIK
C
      UD(1,1)=-GAM*EGAMD
	UD(2,1)=-K*EGAMD
	UD(1,2)=-K*ENUD              *PIK
	UD(2,2)=+NU*ENUD             *PIK   
C
      IF(IES.LE.N)THEN
	  EGAMU=CEXP(+UI*GAM*(ZCAL-Z(IES+1)))
	  ENUU= CEXP(+UI* NU*(ZCAL-Z(IES+1)))
C
C       SUB-MATRICES UU,SU
C
        SU(1,1)=-(NU*NU-K*K*UR)*EGAMU
	  SU(2,1)=+2.0*K*GAM*EGAMU
	  SU(1,2)=+2.0*K*NU*ENUU        *PIK
	  SU(2,2)=+(NU*NU-K*K*UR)*ENUU  *PIK
C         
        UU(1,1)=+GAM*EGAMU
	  UU(2,1)=-K*EGAMU
	  UU(1,2)=-K*ENUU            *PIK
        UU(2,2)=-NU*ENUU           *PIK
      END IF
C
      RETURN
	END
C
      SUBROUTINE DESFH(COME,K,BETA,QS,ZCAL,Z,IES,N,VD,VU,SHD,SHU)
      PARAMETER(NES=16)
      COMMON/COM/UI,UR,ZER
	COMMON/ENU/NU, ENUD, ENUU
      COMPLEX COME, KSC, NU 
      COMPLEX ENUD, ENUU
      COMPLEX SHD,SHU,VD,VU
      COMPLEX UI, UR, ZER
	REAL K, BETA, Z(NES+1)                              
      UI =CMPLX(0.0,1.0)
      UR =CMPLX(1.0,0.0)
      ZER=CMPLX(0.0,0.0)
C
      PI=4.0*ATAN(1.0)
      KSC=COME/BETA*(UR-UI/2.0/QS)
C     COMPUTAION OF NU USING COMPLEX FREQUENCY 
      NU= CSQRT(KSC*KSC-UR*K*K)
      IF(AIMAG(NU) .GT.0.0)NU =-NU
      ENUD= CEXP( -UI*NU*(ZCAL-Z(IES)) )
C     COMPUTATION OF VD,SHD
      SHD=+UI*NU*ENUD
      VD=+ENUD
      IF(IES.LE.N)THEN
	  ENUU= CEXP( +UI*NU*(ZCAL-Z(IES+1)))
C       SE CALCULAN VU,ShU
        SHU=-UI*NU*ENUU 
        VU=ENUU
      END IF
      RETURN
	END
C
	SUBROUTINE IMP (C,NR,NC )
	PARAMETER(NES=16)         
      COMPLEX C(4*NES+2,1)  
      IK =(NC+4)/5                                                     
      DO 10 K=1,IK                                                     
      KK=K-1                                                           
      IIN=KK*5+1                                                       
      IF(K.EQ.IK ) GO TO 20                                            
      IFN = IIN+4                                                      
      GO TO 30                                                         
   20 IFN=NC                                                           
   30 WRITE (30,100)(J,J=IIN,IFN)                                       
  100 FORMAT (/,11X,I3,5(13X,I3))                                      
      DO 10 I=1,NR                                                     
   10 WRITE (30,200) I,(C(I,J),J=IIN,IFN)                               
  200 FORMAT (I3,2X,10F8.4)                                   
      RETURN                                                           
      END 
C
	SUBROUTINE IMPH (C,NR,NC )                                       
	PARAMETER(NES=16)         
      COMPLEX C(2*NES+1,1)                                                 
      IK =(NC+4)/5                                                     
      DO 10 K=1,IK                                                     
      KK=K-1                                                           
      IIN=KK*5+1                                                       
      IF(K.EQ.IK ) GO TO 20                                            
      IFN = IIN+4                                                      
      GO TO 30                                                         
   20 IFN=NC                                                           
   30 WRITE (30,100)(J,J=IIN,IFN)                                       
  100 FORMAT (/,11X,I3,5(13X,I3))                                      
      DO 10 I=1,NR                                                     
   10 WRITE (30,200) I,(C(I,J),J=IIN,IFN)                               
  200 FORMAT (I3,2X,10F8.4)                                   
      RETURN                                                           
      END 
C
      SUBROUTINE CATA(A,B,N,M)               
C                                         
C----- A*TA  --> A   Y    A*TB  --> B         
C                                                     
      COMPLEX A(22,22),B(22,2)
      COMPLEX ASUM(22,22),BSUM(22,2)
      COMPLEX ZER,SUM,CON
      ZER=(0.,0.) 
      DO 20 I=1,N
      DO 20 J=1,N
         SUM=ZER 
         DO 10 K=1,N 
            CON=CONJG(A(K,I))
   10       SUM=SUM+CON*A(K,J)
   20 ASUM(I,J)=SUM 
      ZER=(0.,0.)  
      DO 40 I=1,N  
      DO 40 J=1,M 
         SUM=ZER 
         DO 30 K=1,N 
            CON=CONJG(A(K,I))     
   30       SUM=SUM+CON*B(K,J)   
   40 BSUM(I,J)=SUM 
      DO 50 I=1,N   
      DO 50 J=1,N   
   50 A(I,J)=ASUM(I,J)    
      DO 60 I=1,N    
      DO 60 J=1,M    
   60 B(I,J)=BSUM(I,J)  
      RETURN      
      END    
C        
      SUBROUTINE SOLUCO(N,M)   
C                                                                      
C     SOLUCION DEL SISTEMA DE ECUACIONES AX=B                          
C     A=MATRIZ COMPLEJA DE (N,N)                                       
C     B=MATRIZ COMPLEJA DE (N,M)                                       
C                               
      PARAMETER(NES=16)
      COMMON/ONE/B
      COMMON/MAT/A
      COMPLEX A(4*NES+2,4*NES+2), B(4*NES+2,2)
      COMPLEX SUM        
      EPS=1.0E-25 
	DO 10 I=1,N                                                      
      DO 9 J=I,N                                                       
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 5                                                
      DO 4 K=1,I-1                                                     
    4 SUM=SUM+A(I,K)*A(K,J)                                            
    5 A(I,J)=A(I,J)-SUM                                                
      IF(J.EQ.N)GO TO 9                                                
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 7                                                
      DO 6 K=1,I-1                                                     
    6 SUM=SUM+A(J+1,K)*A(K,I)                                          
    7 IF(CABS(A(I,I)).LE.EPS)GO TO 3                                   
      A(J+1,I)=(A(J+1,I)-SUM)/A(I,I)                                   
    9 CONTINUE                                                         
   10 CONTINUE    
      DO 100 ICAS=1,M                                                      
      DO 17 I=1,N                                                      
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 16                                               
      DO 15 K=1,I-1                                                    
   15 SUM=SUM+A(I,K)*B(K,ICAS)                                              
   16 B(I,ICAS)=B(I,ICAS)-SUM                              
   17 CONTINUE                                                         
      DO 22 I=1,N                                                      
      II=N-I+1                                                         
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 19                                               
      DO 20 K=1,I-1                                                    
      KK=N-K+1                                                         
   20 SUM=SUM+A(II,KK)*B(KK,ICAS)                                           
   19 IF(CABS(A(II,II)).LE.EPS)GO TO 3                                 
      B(II,ICAS)=(B(II,ICAS)-SUM)/A(II,II)                           
   22 CONTINUE      
  100 CONTINUE                                                   
      RETURN                                                           
    3 CONTINUE
C      WRITE(6,*)"***SISTEMA DE ECUACIONES SINGULAR***"                 
      RETURN                                                           
      END                                                              
C                                                                       
C        
      SUBROUTINE SOLUCOH(N,M)   
C                                                                      
C     SOLUCION DEL SISTEMA DE ECUACIONES AX=B                          
C     A=MATRIZ COMPLEJA DE (N,N)                                       
C     B=MATRIZ COMPLEJA DE (N,M)                                       
C                               
      PARAMETER(NES=16)
	COMMON/COM/UI,UR,ZER
      COMMON/ONEH/B
      COMMON/MATH/A
      COMPLEX A(2*NES+1,2*NES+1),B(2*NES+1,1)         
      COMPLEX SUM        
      EPS=1.0E-25 
	DO 10 I=1,N                                                      
      DO 9 J=I,N                                                       
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 5                                                
      DO 4 K=1,I-1                                                     
    4 SUM=SUM+A(I,K)*A(K,J)                                            
    5 A(I,J)=A(I,J)-SUM                                                
      IF(J.EQ.N)GO TO 9                                                
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 7                                                
      DO 6 K=1,I-1                                                     
    6 SUM=SUM+A(J+1,K)*A(K,I)                                          
    7 IF(CABS(A(I,I)).LE.EPS)GO TO 3                                   
      A(J+1,I)=(A(J+1,I)-SUM)/A(I,I)                                   
    9 CONTINUE                                                         
   10 CONTINUE    
      DO 100 ICAS=1,M                                                      
      DO 17 I=1,N                                                      
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 16                                               
      DO 15 K=1,I-1                                                    
   15 SUM=SUM+A(I,K)*B(K,ICAS)                                              
   16 B(I,ICAS)=B(I,ICAS)-SUM                              
   17 CONTINUE                                                         
      DO 22 I=1,N                                                      
      II=N-I+1                                                         
      SUM=(0.0,0.0)                                                    
      IF(I.EQ.1)GO TO 19                                               
      DO 20 K=1,I-1                                                    
      KK=N-K+1                                                         
   20 SUM=SUM+A(II,KK)*B(KK,ICAS)                                           
   19 IF(CABS(A(II,II)).LE.EPS)GO TO 3                                 
      B(II,ICAS)=(B(II,ICAS)-SUM)/A(II,II)                           
   22 CONTINUE      
  100 CONTINUE                                                   
      RETURN                                                           
    3 CONTINUE
C      WRITE(6,*)"***SISTEMA DE ECUACIONES SINGULAR***"                 
      RETURN                                                           
      END                                                              
C                                                                       C
C     CALCULO DE LAS FUNCIONES DE HANKEL          
C                                                                     
      SUBROUTINE HANKEL(Z,H0,H1)
C     Z = COMPLEX ARGUMENT
C
C     COMPUTE SECOND KIND HANKEL FUNCTIONS H0 AND H1
C
      COMPLEX Z,H0,H1,C,A,E,E2,ZH,P
      X=REAL(Z)
      Y=AIMAG(Z)
      R=SQRT(X*X+Y*Y)
      PHI=ATAN2(Y,X)
      IF(R.LE.10.0)GO TO 20
      J=2.0*R
      C=(0.0,0.1250)/Z
      K=2
      P=C*C
      A=4.5*P
      P=7.5*P
      H0=1.0+C+A
      H1=1.0-3.0*C-P
10    I=4*K
      K=K+1
      DI=I
      DK=K
      A=A*C*(DI+1.0/DK)
      P=P*C*(DI-3.0/DK)
      H0=H0+A
      H1=H1-P
      AR=ABS(REAL(P))+ABS(AIMAG(P))
      IF(AR.GT.1.E-16.AND.K.LT.J)GO TO 10
      AR=0.785398163397448-X-PHI/2.0
      E=0.0
      IF(Y.GT.-160.0)
     ^E=0.7978845608028650/SQRT(R)*EXP(Y)*CMPLX(COS(AR),SIN(AR))
      IF(X.EQ.0.0)E=CMPLX(0.0,AIMAG(E))
      H0=H0*E
      H1=H1*E*(0.0,1.0)
      GO TO 23
20    ZH=Z/2.0
      C=-ZH*ZH
      E=CMPLX(0.0,0.3183098861837910)
      E2=E*2.0
      A=1.0-E2*(0.5772156649015330+LOG(R/2.0))+PHI*0.636619772367582
      P=1.0
      K=1
      H0=A
      H1=A+E*(1.0-1.0/C)
25    A=A+E2/K
      P=P*C
      H0=H0+A*P
      K=K+1
      P=P/(K*K)
      H1=H1+(A*K+E)*P
      IF(ABS(REAL(P))+ABS(AIMAG(P)).GT.1.E-16)GO TO 25
      H1=H1*ZH
      IF(X.NE.0.0)GO TO 23
      H0=CMPLX(0.0,AIMAG(H0))
      H1=CMPLX(REAL(H1),0.0)
23    RETURN
      END
C                                                                          
C     CALCULO DE LAS FUNCIONES DE BESSEL          
C                                                
      SUBROUTINE BESJC(V,NMAX,U)              
      COMPLEX Z,R,S,SUM,U(101),UI(101),RR(101),LAMDA,I,C,V  
      INTEGER D                              
      REAL L                                 
      LOGICAL LL                             
C                                            
C     JN(V), N=0,1,...,NMAX                  
C                                            
C     V=ARGUMENTO COMPLEJO                   
C     NMAX=ORDEN MAXIMO                      
C     D=CIFRAS SIGNIFICATIVAS REQUERIDAS     
C     U(N)=FUNCIONES DE BESSEL DE PRIMERA ESPECIE  
C                                       
      D=10                              
C                                       
      EPS=0.5*10.0**(-D)                
      NMA=NMAX+1                        
      X=REAL(V)                         
      Y=AIMAG(V)                        
      I=(0.0,1.0)                       
      IF(((X.LE.0.0).AND.(Y.EQ.0.0)).OR.(NMAX.LT.0).OR.(NMAX.GT.100)) 
     .GO TO 200                         
      DO 10 K=1,NMA                     
   10 UI(K)=(0.0,0.0)                   
      LL=(Y.GE.0.0)                     
      Y=ABS(Y)                          
      Z=CMPLX(X,Y)                      
      R0=CABS(Z)                        
      R02=R0*R0                         
      SUM=CEXP(-I*Z)                    
      D1=2.3026*D+1.3863                
      KK=1                              
      S1=0.0                            
      IF(NMAX.EQ.0)GO TO 50             
      X=0.5*D1/NMAX                     
   20 IF(X.GT.10.0)GO TO 30             
      P=X*5.7941D-5-1.76148D-3          
      P=X*P+2.08645D-2                  
      P=X*P-1.29013D-1                  
      P=X*P+8.57770D-1                  
      S1=X*P+1.01250                    
      GO TO 40                          
   30 Q=ALOG(X)-0.775                   
      P=(0.775-ALOG(Q))/(1.0+Q)         
      P=1.0/(1.0+P)                     
      S1=X*P/Q                          
   40 IF(KK-2)50,60,60                  
   50 R1=S1*NMAX                        
      IF(Y-D1)52,51,51                  
   51 S1=1.3591*R0                      
      GO TO 61                          
   52 X=0.73576*(D1-Y)/R0               
      KK=KK+1                           
      GO TO 20                          
   60 S1=1.3591*R0*S1                   
   61 IF(R1-S1)62,62,63                 
   62 NU=1.0+S1                         
      GO TO 70                          
   63 NU=1.0+R1                         
   70 N=0                               
      L=1.0                             
      C=(1.0,0.0)                       
   80 N=N+1                             
      L=N*L/(N+1)                       
      C=-C*I                            
      IF(N.LT.NU)GO TO 80               
      R=(0.0,0.0)                       
      S=(0.0,0.0)                       
   81 R=1.0/(2.0*N/Z-R)                 
      L=(N+1)*L/N                       
      LAMDA=2.0*N*L*C                   
      C=C*I                             
      S=R*(LAMDA+S)                     
      IF(N.LE.NMAX)RR(N)=R              
      N=N-1                             
      IF(N.GE.1)GO TO 81                
      U(1)=SUM/(1.0+S)                  
      IF(NMAX.EQ.0)GO TO 90             
      DO 82 K=1,NMAX                    
   82 U(K+1)=U(K)*RR(K)                 
   90 IF(LL)GO TO 92                    
      DO 91 K=1,NMA                     
   91 U(K)=CONJG(U(K))                  
   92 CONTINUE                          
      DO 100 K=1,NMA                    
      IF(CABS(U(K)-UI(K))/CABS(U(K)).GT.EPS)GO TO 101  
  100 CONTINUE                          
      RETURN                            
  101 NU=NU+5                          
      DO 102 K=1,NMA                   
  102 UI(K)=U(K)                       
      GO TO 70                         
  200 WRITE(6,*)"ORDEN O ARGUMENTO INVALIDO"     
      RETURN                       
      END                          
C                        
      SUBROUTINE SINTE(IJSAL,NEST,N,DT,VIJ,OMEI,CR)
      COMPLEX VIJ(512,25)
      COMPLEX CR(512), CY(512)
      DIMENSION TV(512),VT(512)
      AN=N
	RAIZ=SQRT(AN)                      
      DO IX=1,NEST   
         CY(1)=CR(1)*RAIZ  
         DO J=2,N/2+1             
            CY(J)=CR(J)*VIJ(J-1,IX)*RAIZ 
            JN=N-J+2                     
            CY(JN)=CONJG(CY(J))           
         END DO                             
         CALL FORK(N,CY,1.0)              
         DO J=1,N                     
            CY(J)=CY(J)/RAIZ              
            TV(J)=DT*(J-1)            
            VT(J)=REAL(CY(J))
            WRITE(IJSAL,200)TV(J),VT(J)*EXP(OMEI*TV(J)) !CORRECCION POR TW
         END DO 
      END DO 
  200 FORMAT(2E15.6)
      RETURN
      END
C
      FUNCTION RIC(J)         
      COMMON/SINT/DT,TS,TP
      PI=4.0*ATAN(1.0)    
      A=PI*(DT*J-TS)/TP  
      A=A*A                  
      RIC=0.0                
      IF(A.GT.30.0)RETURN    
      RIC=(A-0.5)*EXP(-A)    
      RETURN              
      END          
C                                     
      SUBROUTINE FORK(LX,CX,SIGNI)      
      COMPLEX CX(LX),CARG,CW,CTEMP      
      J=1                               
      SC=SQRT(1.0/LX)                   
      DO 30 I=1,LX                      
      IF(I.GT.J)GO TO 10                
      CTEMP=CX(J)*SC                    
      CX(J)=CX(I)*SC                    
      CX(I)=CTEMP                       
   10 M=LX/2                            
   20 IF(J.LE.M)GO TO 30                
      J=J-M                             
      M=M/2                             
      IF(M.GE.1)GO TO 20                
   30 J=J+M                             
      L=1                               
   40 ISTEP=2*L                         
      DO 50 M=1,L                       
      CARG=(0.0,1.0)*(3.141592654*SIGNI*(M-1))/L  
	CW=CEXP(CARG)                
      DO 50 I=M,LX,ISTEP           
      CTEMP=CW*CX(I+L)             
      CX(I+L)=CX(I)-CTEMP          
   50 CX(I)=CX(I)+CTEMP            
      L=ISTEP                      
      IF(L.LT.LX)GO TO 40          
      RETURN                       
      END                        
C                                                         