	PROGRAM CARBYNE
	IMPLICIT NONE
	EXTERNAL BESS_NUL, BESS, F02HAF
	INTEGER, PARAMETER :: MAXSUM=50
C	?enei aaeaiee a naoea ii R
	INTEGER, Parameter:: NGRID = 300
C	?enei oeiia aoiiia
	INTEGER, Parameter:: NT = 1

C	Iaeneiaeuiia cia?aiea eaaioiaiai ?enea l

C	INTEGER, Parameter:: NLMAX = 12
      DOUBLE PRECISION, Parameter:: evs = 13.6058D0


	CHARACTER, PARAMETER :: JOB*1 = 'V'
C	Constraint: JOB = 'N' or 'V'.
	CHARACTER, PARAMETER :: UPLO*1 = 'L'

	DOUBLE PRECISION, Parameter:: 
	+LightSpeed = 274.0746D0
	DOUBLE PRECISION, Parameter:: 
	+PI = 3.14159265358979323846264338328D0

	DOUBLE COMPLEX, ALLOCATABLE :: AHA(:,:), ASMALL(:,:,:,:), 
	+BSMALL(:,:,:,:), HamAHA(:,:),WORK(:),
	+H(:,:), BigH(:,:,:)

	DOUBLE COMPLEX :: ResPoM, SIMP_I12, SUM_VAR,ResPoL, ResPoX
	DOUBLE PRECISION :: RES(MAXSUM), A, B, C, RMT, EPSI, BESSLJ,
	+BESSLJ1, K, Kstart, Kstop
	DOUBLE PRECISION, ALLOCATABLE :: KAPPA(:), RAD(:,:), BigD(:,:), 
	+RO(:), PHI(:), RSZC(:), PSIALL(:,:,:),PSIEALL(:,:,:),D(:),
	+PSIRALL(:,:,:),PSIERALL(:,:,:), V(:,:), RWORK(:), DE(:),EBig(:,:)
	INTEGER, ALLOCATABLE :: M(:), N(:), P(:)
	INTEGER :: COUNTX, I, I_, JRIS, NATOMS, X, L, MSMALL, J,
	+LENS, LENH, LEND, IT, NINTERVALS, A1, A2, A3, NPNTS, NLMAX,
	+ II, JJ, INFO, LWORK, I1, I1_

	CHARACTER(1) CHAR

C************************************************************************************************************
C	 N?eouaaiea COUNTX, M, N, P
C************************************************************************************************************
	OPEN(1, FILE = 'indefect9.dat', FORM = 'FORMATTED')
	READ(1,*) COUNTX
	ALLOCATE (KAPPA(COUNTX), M(COUNTX), N(COUNTX), P(COUNTX), 
	+AHA(COUNTX,COUNTX),HamAHA(CountX, CountX),D(COUNTX))
	READ(1,*) KAPPA
	READ(1,*) A, B, C
	DO I = 1, COUNTX
		READ(1,*) M(I), N(I), P(I)
	END DO
	CLOSE(1)
C	Ia?aiaiiua aey aeaaiiaeecaoee
	INFO=0
	LWORK = 64 * CountX * 2
	ALLOCATE(RWORK(7*CountX*2), WORK(LWORK))

C************************************************************************************************************
C	 ?AN?AO KAPPA
C************************************************************************************************************

      DO I=1,COUNTX	   
		CALL BESS_NUL(N(I),M(I),RES,MAXSUM)
		KAPPA(I)=RES(N(I))/A
c		WRITE(*,*) I, KAPPA(I)
	END DO

	OPEN(1, FILE = 'kStart_kStop.txt', FORM = 'FORMATTED')
	READ(1, *) kStart, kStop
	CLOSE(1)
C************************************************************************************************************
C	 N?eouaaiea JRIS-iiia? rmt a rad, RMT, RAD
C************************************************************************************************************
	ALLOCATE(RAD(NGRID,NT), V(NGRID,NT)) 
	OPEN(1, FILE = 'indefect2.dat', FORM = 'FORMATTED')
	READ(1,*) JRIS, RMT
	READ(1,*) RAD
c	WRITE(*,*) RAD
	CLOSE(1)

C************************************************************************************************************
C	 N?eouaaiea EPSI, NATOMS
C************************************************************************************************************

	OPEN(13, FILE = 'EPSI_NATOMSG.TXT', FORM = 'FORMATTED')
	READ(13,*) EPSI, NATOMS
	CLOSE(13)

C	IINEA N?EOUAAIE? ?ENEA AOIIIA A YEAIAIOA?IIE ??AEEA NATOMS II?II II?AAAEEOU IANNEAU OEEEIA?E?ANEEO EII?AEIAO

	ALLOCATE(RO(NATOMS), PHI(NATOMS), RSZC(NATOMS))

C************************************************************************************************************
C	 N?EOUAAIEA OEEEIA?E?ANEEO EII?AEIAO, RO - RHO, PHI, RSZC - Z; 
C************************************************************************************************************
	
   	OPEN(24, FILE = 'CYLINDRICALCOORDINATESG.TXT', FORM = 'FORMATTED')
	DO I=1,NATOMS
		READ(24, '(A1,3F20.10)') CHAR, RO(I),PHI(I),RSZC(I)
	END DO

	OPEN(15, FILE = 'NPNTS_LMAX_LENH_LENDG.TXT', FORM = 'FORMATTED')
	READ(15, *) NPNTS, NLMAX, LENS, LENH, LEND
C	NLMAX - YOI IAENEIAEUIIA CIA?AIEA L IAEAIUEIA(EAAIOIAIA ?ENEI)
	CLOSE(15)
	NLMAX=NLMAX-1
	NLMAX=1
	ALLOCATE(ASMALL(NLMAX,2*NLMAX+1,COUNTX,NPNTS),
	+BSMALL(NLMAX,2*NLMAX+1,COUNTX,NPNTS))
	ALLOCATE(BIGH(NPNTS,COUNTX,COUNTX), BIGD(NPNTS,COUNTX),
	+H(COUNTX,COUNTX),DE(CountX),
     +EBig(NPNTS,COUNTX))
	NINTERVALS=NPNTS-1
C************************************************************************************************************
C	N?eouaaiea PSIALL(NT, NLMAX, EPNTS), PSIEALL(NT, NLMAX, EPNTS), PSIRALL, PSIERALL
C************************************************************************************************************
	ALLOCATE(PSIALL(NT,NLMAX+1,JRIS),PSIEALL(NT,NLMAX+1,JRIS),
	+PSIRALL(NT,NLMAX+1,JRIS),PSIERALL(NT,NLMAX+1,JRIS))
  	OPEN(1, FILE = 'indefect4.dat', FORM = 'FORMATTED')

	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS
				READ(1, *) PSIALL(IT, L, X), A1, A2, A3

				PSIALL(IT,L,X) = PSIALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)

 	OPEN(1, FILE = 'indefect5.dat', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS
				READ(1, *) PSIEALL(IT, L, X), A1, A2, A3

				PSIEALL(IT, L, X)=PSIEALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)
	OPEN(1, FILE = 'PSIRALL.TXT', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS
				READ(1, *) PSIRALL(IT, L, X), A1, A2, A3

				PSIRALL(IT, L, X)=PSIRALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)
	OPEN(12, FILE='SO_EnergiesAll.txt')
	write(12,*)''
	close(12)
	OPEN(1, FILE = 'PSIERALL.TXT', FORM = 'FORMATTED')
	DO IT = 1, NT
		DO L = 1, NLMAX+1
			DO X = 1, JRIS
				READ(1, *) PSIERALL(IT, L, X), A1, A2, A3

				PSIERALL(IT, L, X)=PSIERALL(IT, L, X)/RAD(X,1)
			END DO
		END DO
	END DO
	CLOSE(1)

	OPEN(1, FILE = 'indefect3.dat', FORM = 'FORMATTED')
	DO IT = 1, NT
		READ(1,*) V(:,IT)
	END DO
	CLOSE(1)
	OPEN(22,FILE='H_AMNPG.TXT', ACCESS = 'DIRECT',
	+RECL = LENH+1000)
	open(101,file='coeff.txt')
	DO J = 1, NPNTS
	 READ(22, REC=J)((H(II,JJ),JJ=1,COUNTX),II=1,COUNTX)
	 DO JJ = 1, COUNTX
	 DO II = 1, COUNTX
	  BIGH(J,II,JJ) = H(II,JJ)
		write(101,*)II,JJ,H(II,JJ)
	  END DO
	 END DO
	END	DO
	CLOSE(22)
	close(101)
	open(32,FILE='d_en.txt')
	OPEN(23,FILE='D_ENERGIESG.TXT',
	+ACCESS = 'DIRECT', RECL = LEND+1000)
	DO J=1, NPNTS
	 READ(23,REC=J) D 
	 DO II = 1, COUNTX
	  BIGD(J,II) = D(II)
	write(32,*)II,D(II)
	 END DO
	END DO
	CLOSE(23)
	CLOSE(32)
C    ***************************************************?an?ao A e B iaeaiueeo************************************
	OPEN(213, FILE='new_temp.txt')
	DO I=1, COUNTX
		DO L=1,NLMAX
			DO MSMALL=-L, L
				DO J=1, NPNTS

	K=Kstart + (kStop - kStart) * (DBLE(J-1)/DBLE(NINTERVALS))
	ASMALL(L,MSMALL+L+1,I,J)=PSIEALL(1,L+1,JRIS)*
	+SIMP_I12(EPSI, RMT,K+(2*PI/C)*DBLE(P(I)),KAPPA(I), MSMALL,L,2)
     +-PSIERALL(1,L+1,JRIS)*
     +SIMP_I12(EPSI, RMT,K+(2*PI/C)*DBLE(P(I)), KAPPA(I), MSMALL,L,1)
	BSMALL(L,MSMALL+L+1,I,J)=PSIRALL(1,L+1,JRIS)*SIMP_I12(EPSI, RMT,K+ 
	+(2*PI/C)*DBLE(P(I)),KAPPA(I), MSMALL,L,1)-PSIALL(1,L+1,JRIS)*
     +SIMP_I12(EPSI, RMT,K+(2*PI/C)*DBLE(P(I)), KAPPA(I), MSMALL,L,2)
				END DO
			END DO
		END DO
		write(213,*)ASMALL(1,1,I,1),BSMALL(1,1,I,1)
	END DO
	close(213)

C    ****************************************ieii?aiea ?an?aoa A e B iaeaiueeo************************************
      DO J=NPNTS,1,-1
		DO I=1,COUNTX
			DO I_=1,COUNTX
				HamAHA(I,I_)=DCMPLX(0.0d0,0.0d0)
				IF (M(I).EQ.M(I_)) THEN
						SUM_VAR=RMT**4/(LightSpeed**2*C*A**2*2)
						CALL BESS(M(I),KAPPA(I)*A,BESSLJ,BESSLJ1)
						SUM_VAR=SUM_VAR/BESSLJ1
						CALL BESS(M(I_),KAPPA(I_)*A,BESSLJ,BESSLJ1)
						SUM_VAR=SUM_VAR/BESSLJ1
						MSMALL=M(I)
							ResPoL=0
						DO L=ABS(M(I)),NLMAX
							ResPoM=0
							DO X=1, NATOMS
								ResPoX=2*
	+int1_full1(V(:,1), PSIALL(1,L+1,:),PSIRALL(1,L+1,:), RAD(:,1),
     +JRIS)*DCONJG(ASMALL(L,MSMALL+L+1,I,J))*ASMALL(L,MSMALL+L+1,I_,J)+
	+2*int1_full1(V(:,1),PSIEALL(1,L+1,:),PSIERALL(1,L+1,:),
     +RAD(:,1),JRIS)*
     +DCONJG(BSMALL(L,MSMALL+L+1,I,J))*BSMALL(L,MSMALL+L+1,I_,J)
     ++(int1_full1(V(:,1), PSIALL(1,L+1,:), PSIERALL(1,L+1,:)
     +,RAD(:,1),JRIS)+int1_full1(V(:,1), PSIEALL(1,
     +L+1,:),PSIRALL(1,L+1,:),RAD(:,1),JRIS))
	+*(DCONJG(ASMALL(L,MSMALL+L+1,I,J))*BSMALL(L,MSMALL+L+1,I_,J)+
	+DCONJG(BSMALL(L,MSMALL+L+1,I,J))*ASMALL(L,MSMALL+L+1,I_,J))
		ResPoM=ResPoM+ResPoX*DCMPLX(DCOS((2*PI/C)*DBLE(P(I)-P(I_))*
	+	RSZC(X)),DSIN((2*PI/C)*DBLE(P(I)-P(I_))*RSZC(X)))
							END DO
							ResPoL=ResPoL+ResPoM*(2*L+1)
	+						*FCT(L-ABS(MSMALL))/FCT(L+ABS(MSMALL))
						END DO
			HamAHA(I,I_)=SUM_VAR*ResPoL
					END IF					  
			END DO
		END DO

		DO I1=1,COUNTX
			DO I1_=1,COUNTX
				AHA(I1,I1_)=DCMPLX(0.0D0,0.0D0)
				if (I1.GE.I1_)then
					DO I=1,COUNTX
						DO I_=1,COUNTX
							AHA(I1,I1_)=AHA(I1,I1_)-HamAHA(I,I_)*
	+						DCONJG(BIGH(J,I,I1))*BIGH(J,I_,I1_)
						END DO
					END DO
				END IF
			END DO
		END DO
	write(*,*)J,NPNTS
		DO I=1, CountX
				AHA(I,I) = AHA(I,I) + BIGD(J,I)
		END DO 
      info=0
		CALL F02HAF(JOB,UPLO,COUNTX,AHA,COUNTX,DE,
	+RWORK,WORK,LWORK,INFO)
		OPEN(100,FILE = 'SO_EnergiesAll.txt', STATUS = 'OLD',
	+FORM = 'FORMATTED', POSITION = 'APPEND')
		EBig(J,:) = DE
		DO I = 1, CountX
			EBig(J,I) = EBig(J,I)*evs
			K= kStart + (kStop - kStart) *(DBLE(J-1)/DBLE(NINTERVALS))
			WRITE(100,*) K, EBig(J,I)
	 
		END DO

		close(100)
	END DO
	CLOSE(12)
	CONTAINS
	
      INTEGER FUNCTION FCT(X)
	INTEGER :: X, I, RES
	RES=1
	DO I=2,X
		RES=RES*I
	END DO
	FCT=RES
	END FUNCTION

c     ************************************Eioaa?ae acaoa2 *********************************************************
	DOUBLE PRECISION FUNCTION dVdr(V1,V2,rad, rad1)
	DOUBLE PRECISION :: V1,V2,rad, rad1
		dVdr=(V1-V2)/(rad-rad1)
	END FUNCTION dVdr

	DOUBLE PRECISION FUNCTION pod_int1_full1(V1,V2,U1,U2,rad,rad1)
	DOUBLE PRECISION :: V1,V2,U1,U2,rad,rad1
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod_int1_full1 = U1*U2*rad**2*dVdr(V1,V2,rad,rad1)

	END FUNCTION pod_int1_full1

	DOUBLE PRECISION FUNCTION int1_full1(V, U1, U2, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), U1(:), U2(:), rad(:), integ_1_1

	integ_1_1 = (pod_int1_full1(V(1),V(2),U1(1),U2(1),rad(1),rad(2))+
	+			pod_int1_full1(V(2),V(1),U1(2),U2(2),rad(2),rad(1)))
	+			/2.D0 * (rad(2)-rad(1))
	DO i=2, JRIS-1
		integ_1_1 = integ_1_1 + (
	+	pod_int1_full1(V(i),V(i-1),U1(i),U2(i),rad(i),rad(i-1))+ 
	+	pod_int1_full1(V(i+1),V(i),U1(i+1),U2(i+1),rad(i+1),rad(i)))
     +	/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int1_full1=integ_1_1

	END FUNCTION int1_full1


	DOUBLE PRECISION FUNCTION pod1_1(V,U, UR, rad)
	DOUBLE PRECISION :: V,U,UR, rad
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod1_1 = V*U*UR*rad

	END FUNCTION pod1_1

	DOUBLE PRECISION FUNCTION pod1_2(V,U,UR,UR1,rad,rad1)
	DOUBLE PRECISION :: V,U,UR,UR1, rad,rad1
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod1_2 = V*U*(UR-UR1)/(rad-rad1)*rad**2

	END FUNCTION pod1_2

	DOUBLE PRECISION FUNCTION pod1_3(V, UR, rad)
	DOUBLE PRECISION :: V,UR, rad
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod1_3 = V*UR**2*rad**2

	END FUNCTION pod1_3

	DOUBLE PRECISION FUNCTION int1_1(V, u, ur, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), u(:), ur(:), rad(:), integ_1_1

	integ_1_1 = 0.d0

	DO i=1, JRIS-1
		integ_1_1 = integ_1_1 + (pod1_1(V(i),u(i),ur(i),rad(i))+ 
	+pod1_1(V(i+1),u(i+1),ur(i+1),rad(i+1)))/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int1_1=integ_1_1

	END FUNCTION int1_1

	DOUBLE PRECISION FUNCTION int1_2(V, u, ur, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), u(:), ur(:), rad(:), integ_1_2

	integ_1_2 = (-pod1_2(V(1),u(1),ur(2),ur(1),rad(1),rad(2))+
	+			pod1_2(V(2),u(2),ur(2),ur(1),rad(2),rad(1)))/
	+			2.D0 * (rad(2)-rad(1))


	DO i=2, JRIS-1
		integ_1_2 = integ_1_2 + 
	+			(pod1_2(V(i),u(i),ur(i),ur(i-1),rad(i),rad(i-1))+ 
	+			pod1_2(V(i+1),u(i+1),ur(i+1),ur(i),rad(i+1),rad(i)))
     +			/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int1_2=integ_1_2

	END FUNCTION int1_2

	DOUBLE PRECISION FUNCTION int1_3(V, ur, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), ur(:), rad(:), integ_1_3

	integ_1_3 = 0.d0

	DO i=1, JRIS-1
		integ_1_3 = integ_1_3 + (pod1_3(V(i),ur(i),rad(i))+ 
	+pod1_3(V(i+1),ur(i+1),rad(i+1)))/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int1_3=integ_1_3

	END FUNCTION int1_3

	DOUBLE PRECISION FUNCTION int1_full(V, U, UR, rad,JRIS)
	DOUBLE PRECISION :: V(:), UR(:), U(:), rad(:)
	INTEGER :: JRIS


	int1_full=V(JRIS)*U(JRIS)*UR(JRIS)*rad(JRIS)**2-
	+	V(1)*U(1)*UR(1)*rad(1)**2
	+	-(2*int1_1(V, U, UR, rad, JRIS)+
	+	int1_2(V, U, UR, rad, JRIS)+
	+	int1_3(V, UR, rad, JRIS))

	END FUNCTION int1_full


C************************************************************************************************************

	DOUBLE PRECISION FUNCTION int2_full(V, UE, UER, rad,JRIS)
	DOUBLE PRECISION :: V(:), UER(:), UE(:), rad(:)
	INTEGER :: JRIS


	int2_full=V(JRIS)*UE(JRIS)*UER(JRIS)*rad(JRIS)**2-
	+	V(1)*UE(1)*UER(1)*rad(1)**2
	+	-(2*int1_1(V, UE, UER, rad, JRIS)+
	+	int1_2(V, UE, UER, rad, JRIS)+
	+	int1_3(V, UER, rad, JRIS))

	END FUNCTION int2_full

c     ************************************Eioaa?ae acaoa2 *********************************************************

	DOUBLE PRECISION FUNCTION pod3_1(V,U, UER, rad)
	DOUBLE PRECISION :: V,U,UER, rad
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod3_1 = V*U*UER*rad

	END FUNCTION pod3_1

	DOUBLE PRECISION FUNCTION pod3_2(V,U,UER,UER1,rad,rad1)
	DOUBLE PRECISION :: V,U,UER,UER1, rad,rad1
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod3_2 = V*U*(UER-UER1)/(rad-rad1)*rad**2

	END FUNCTION pod3_2

	DOUBLE PRECISION FUNCTION pod3_3(V, UR, UER, rad)
	DOUBLE PRECISION :: V,UR, UER, rad
C	Ooieoey au?eneyao cia?aiea iiaeioaa?aeuiiai au?a?aiey Vu2dr(L,IT,i)

	pod3_3 = V*UR*UER*rad**2

	END FUNCTION pod3_3

	DOUBLE PRECISION FUNCTION int3_1(V, u, uer, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), u(:), uer(:), rad(:), integ_3_1

	integ_3_1 = 0.d0

	DO i=1, JRIS-1
		integ_3_1 = integ_3_1 + (pod3_1(V(i),u(i),uer(i),rad(i))+ 
	+pod3_1(V(i+1),u(i+1),uer(i+1),rad(i+1)))/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int3_1=integ_3_1

	END FUNCTION int3_1

	DOUBLE PRECISION FUNCTION int3_2(V, u, uer, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), u(:), uer(:), rad(:), integ_3_2

	integ_3_2 = (-pod3_2(V(1),u(1),uer(2),uer(1),rad(1),rad(2))+
	+			pod3_2(V(2),u(2),uer(2),uer(1),rad(2),rad(1)))/
	+			2.D0 * (rad(2)-rad(1))


	DO i=2, JRIS-1
		integ_3_2 = integ_3_2 + 
	+			(pod3_2(V(i),u(i),uer(i),uer(i-1),rad(i),rad(i-1))+ 
	+			pod3_2(V(i+1),u(i+1),uer(i+1),uer(i),rad(i+1),rad(i)))
     +			/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int3_2=integ_3_2

	END FUNCTION int3_2

	DOUBLE PRECISION FUNCTION int3_3(V, ur, uer, rad, JRIS)
	INTEGER :: i, JRIS
	DOUBLE PRECISION :: V(:), ur(:), uer(:), rad(:), integ_3_3

	integ_3_3 = 0.d0

	DO i=1, JRIS-1
		integ_3_3 = integ_3_3 + (pod3_3(V(i),ur(i),uer(i),rad(i))+ 
	+pod3_3(V(i+1),ur(i+1),uer(i+1),rad(i+1)))/2.D0 * (rad(i+1)-rad(i))
      ENDDO

	int3_3=integ_3_3

	END FUNCTION int3_3

	DOUBLE PRECISION FUNCTION int3_full(V, U, UR, UER, rad,JRIS)
	DOUBLE PRECISION :: V(:), UR(:), U(:), UER(:), rad(:)
	INTEGER :: JRIS


	int3_full=V(JRIS)*U(JRIS)*UER(JRIS)*rad(JRIS)**2-
	+	V(1)*U(1)*UER(1)*rad(1)**2
	+	-(2*int3_1(V, U, UER, rad, JRIS)+
	+	int3_2(V, U, UER, rad, JRIS)+
	+	int3_3(V, UR, UER, rad, JRIS))

	END FUNCTION int3_full


C************************************************************************************************************

	DOUBLE PRECISION FUNCTION int4_full(V, UE, UER, UR, rad,JRIS)
	DOUBLE PRECISION :: V(:), UER(:), UE(:), UR(:), rad(:)
	INTEGER :: JRIS


	int4_full=V(JRIS)*UE(JRIS)*UR(JRIS)*rad(JRIS)**2-
	+	V(1)*UE(1)*UE(1)*rad(1)**2
	+	-(2*int3_1(V, UE, UR, rad, JRIS)+
	+	int3_2(V, UE, UR, rad, JRIS)+
	+	int3_3(V, UER, UR, rad, JRIS))

	END FUNCTION int4_full

C************************************************************************************************************

	END PROGRAM