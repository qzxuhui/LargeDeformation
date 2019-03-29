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
C     LOCAL ARRAYS
C---------------------------------------------------------------------
C     EELAS  - LOGARITHMIC ELASTIC STRAINS
C     EELASP - PRINCIPAL ELASTIC STRAINS
C     bbar   - DEVIATORIC RIGHT CAUCHY-GREEN TENSOR
C     bbarP  - PRINCIPAL VALUES OF bbar
C     bbarN  - PRINCIPAL DIRECTION OF bbar (AND EELAS)
C     DISTGR - DEVIATORIC DEFORMATION GRADIENT (DISTORTION TENSOR)
C----------------------------------------------------------------------
C
      DIMENSION EELAS(6), EELASP(3), bbar(6), bbarP(3), bbarN(3,3),
     1DISTGR(3,3),bmat(6)
C
      PARAMETER(ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0,
     1          FOUR = 4.D0)
C
C----------------------------------------------------------------------
C     UMAT for compressible Neo-Hookean hyperelastity
C----------------------------------------------------------------------
      
C     elastic preproties
      C10  = PROPS(1)
      D1   = PROPS(2)
      mu_0 = PROPS(3)
      lambda_0 = PROPS(4)
      
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
      
c     calculate distortion tensor (will be no use)
      SCALE = DET**(-ONE /THREE)
      DO K1 = 1, 3
        DO K2 = 1, 3
          DISTGR(K2,K1) = SCALE * DFGRD1(K2,K1)
        END DO
      END DO
      
C     calculate left cauchy-green tensor's bar (will be no use)
      bbar(1) = DISTGR(1,1)**2 + DISTGR(1,2)**2 + DISTGR(1,3)**2
      bbar(2) = DISTGR(2,1)**2 + DISTGR(2,2)**2 + DISTGR(2,3)**2
      bbar(3) = DISTGR(3,3)**2 + DISTGR(3,1)**2 + DISTGR(3,2)**2
      bbar(4) = DISTGR(1,1)*DISTGR(2,1) + DISTGR(1,2)*DISTGR(2,2) 
     1        + DISTGR(1,3)*DISTGR(2,3)
      IF (NSHR .EQ. 3) THEN
        bbar(5) = DISTGR(1,1)*DISTGR(3,1) + DISTGR(1,2)*DISTGR(3,2) 
     1          + DISTGR(1,3)*DISTGR(3,3)
        bbar(6) = DISTGR(2,1)*DISTGR(3,1) + DISTGR(2,2)*DISTGR(3,2) 
     1          + DISTGR(2,3)*DISTGR(3,3)
      END IF

C     calculate cauchy stress
      TRbbar = (bbar(1) + bbar(2) + bbar(3)) / THREE
      EG     = TWO * C10 / DET
      EK     = TWO / D1 * (TWO * DET - ONE)
      PR     = TWO / D1 * (DET - ONE)
      DO K1 = 1, NDI
        STRESS(K1) = EG * (bbar(K1) - TRbbar) + PR
      END DO
      DO K1 = NDI+1, NDI+NSHR
        STRESS(K1) = EG * bbar(K1)
      END DO
     
C     calculate material jacobian
      EG23 = EG * TWO / THREE
      DDSDDE(1,1) = EG23 * (bbar(1) + TRbbar) + EK
      DDSDDE(2,2) = EG23 * (bbar(2) + TRbbar) + EK
      DDSDDE(3,3) = EG23 * (bbar(3) + TRbbar) + EK
      DDSDDE(1,2) =-EG23 * (bbar(1) + bbar(2) - TRbbar) + EK
      DDSDDE(1,3) =-EG23 * (bbar(1) + bbar(3) - TRbbar) + EK
      DDSDDE(2,3) =-EG23 * (bbar(2) + bbar(3) - TRbbar) + EK
      DDSDDE(1,4) = EG23 * bbar(4) / TWO
      DDSDDE(2,4) = EG23 * bbar(4) / TWO
      DDSDDE(3,4) =-EG23 * bbar(4)
      DDSDDE(4,4) = EG * (bbar(1) + bbar(2)) / TWO
      IF (NSHR .EQ. 3) THEN
        DDSDDE(1,5) = EG23 * bbar(5) / TWO
        DDSDDE(2,5) =-EG23 * bbar(5)
        DDSDDE(3,5) = EG23 * bbar(5) / TWO
        DDSDDE(1,6) =-EG23 * bbar(6)
        DDSDDE(2,6) = EG23 * bbar(6) / TWO
        DDSDDE(3,6) = EG23 * bbar(6) / TWO
        DDSDDE(5,5) = EG * (bbar(1) + bbar(3)) / TWO
        DDSDDE(6,6) = EG * (bbar(2) + bbar(3)) / TWO
        DDSDDE(4,5) = EG * bbar(6) / TWO
        DDSDDE(4,6) = EG * bbar(5) / TWO
        DDSDDE(5,6) = EG * bbar(4) / TWO
      END IF
      DO K1 = 1, NTENS
        DO K2 = 1, K1 - 1
          DDSDDE(K1, K2) = DDSDDE(K2, K1)
        END DO
      END DO

      RETURN
      END

    
