C----------------------------------------------------------------------
C     UMAT for compressible Neo-Hookean hyperelastity
C----------------------------------------------------------------------
      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
C     
C     abaqus job=Job-1.inp user=neo_hookean_my_compressible_umat.for int
C     LOCAL ARRAYS
C---------------------------------------------------------------------
C     EELAS  - LOGARITHMIC ELASTIC STRAINS
C     EELASP - PRINCIPAL ELASTIC STRAINS
C     bbar   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     bbarP  - PRINCIPAL VALUES OF bbar
C     bbarN  - PRINCIPAL DIRECTION OF bbar (AND EELAS)
C     DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C----------------------------------------------------------------------
C     STRSEE 11 22 33 12 13 23
      DIMENSION EELAS(6), EELASP(3), bbar(6), bbarP(3), bbarN(3,3),
     1DISTGR(3,3),bmat(6)
C
      PARAMETER(ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0,
     1          FOUR = 4.D0)
      
      real mu_0
      real lambda_0
      real mu
	real lambda
      real pressure
      
C     elastic preproties
      mu_0 = 2
      lambda_0 = 1000000.D0
      
C     calculate jacobian det(F)
      DET = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1    - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
      IF (NSHR .EQ. 3) THEN
          DET = DET + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     1              + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     2              - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     3              - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)
      END IF
      
c     calculate left cauchy-green tensor
      bmat(1) = DFGRD1(1,1)**2 + DFGRD1(1,2)**2 + DFGRD1(1,3)**2
      bmat(2) = DFGRD1(2,1)**2 + DFGRD1(2,2)**2 + DFGRD1(2,3)**2
      bmat(3) = DFGRD1(3,3)**2 + DFGRD1(3,1)**2 + DFGRD1(3,2)**2
      bmat(4) = DFGRD1(1,1)*DFGRD1(2,1) + DFGRD1(1,2)*DFGRD1(2,2) 
     1        + DFGRD1(1,3)*DFGRD1(2,3)
      IF (NSHR .EQ. 3) THEN
        bmat(5) = DFGRD1(1,1)*DFGRD1(3,1) + DFGRD1(1,2)*DFGRD1(3,2) 
     1          + DFGRD1(1,3)*DFGRD1(3,3)
        bmat(6) = DFGRD1(2,1)*DFGRD1(3,1) + DFGRD1(2,2)*DFGRD1(3,2) 
     1          + DFGRD1(2,3)*DFGRD1(3,3)
      END IF
      
C     calculate modified material perproties
      lambda = lambda_0
      mu = mu_0 - lambda*log(DET)

C     calculate cauchy stress
      EGXH   = mu_0 / DET
      PRXH   = -mu / DET 
      
      DO K1 = 1, NDI
        STRESS(K1) = EGXH * bmat(K1) + PRXH
      END DO 
      DO K1 = NDI+1, NDI+NSHR
        STRESS(K1) = EGXH * bmat(K1)
      END DO
      
C     calculate DDSDDE      
      DDSDDE(1,1) = (lambda + TWO * mu) / DET  + TWO * STRESS(1)
      DDSDDE(2,2) = (lambda + TWO * mu) / DET + TWO * STRESS(2)
      DDSDDE(3,3) = (lambda + TWO * mu) / DET + TWO * STRESS(3)
      DDSDDE(1,2) = lambda / DET
      DDSDDE(1,3) = lambda / DET
      DDSDDE(2,3) = lambda / DET
      DDSDDE(1,4) = ZERO + STRESS(4)
      DDSDDE(2,4) = ZERO + STRESS(4)
      DDSDDE(3,4) = ZERO
      DDSDDE(4,4) =  mu / DET + (STRESS(1) + STRESS(2)) / TWO 
      IF (NSHR .EQ. 3) THEN
        DDSDDE(1,5) = ZERO + STRESS(5)
        DDSDDE(2,5) = ZERO 
        DDSDDE(3,5) = ZERO + STRESS(5)
        DDSDDE(1,6) = ZERO
        DDSDDE(2,6) = ZERO + STRESS(6)
        DDSDDE(3,6) = ZERO + STRESS(6)
        DDSDDE(5,5) =  mu / DET + (STRESS(1) + STRESS(3)) / TWO 
        DDSDDE(6,6) =  mu / DET + (STRESS(2) + STRESS(3)) / TWO 
        DDSDDE(4,5) = ZERO + STRESS(6) / TWO
        DDSDDE(4,6) = ZERO + STRESS(5) / TWO
        DDSDDE(5,6) = ZERO + STRESS(4) / TWO
      END IF    
      DO K1 = 1, NTENS
        DO K2 = 1, K1 - 1
          DDSDDE(K1, K2) = DDSDDE(K2, K1)
        END DO
      END DO
      
      RETURN
      END