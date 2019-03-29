      program main
      
      
      CHARACTER*8 CMNAME
      DIMENSION STRESS(6),STATEV(10),
     1 DDSDDE(6,6),DDSDDT(6),DRPLDE(6),
     2 STRAN(6),DSTRAN(6),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(4),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      DIMENSION EELAS(6), EELASP(3), bbar(6), bbarP(3), bbarN(3,3),
     1DISTGR(3,3),bmat(6)
C
      PARAMETER(ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, THREE = 3.D0,
     1          FOUR = 4.D0)

      real(8) mu_0
      real(8) lambda_0
C     interface:
      NTENS=6
      NSHR=3
      NDI=3
      PROPS(1)=1.D0;
      PROPS(2)=1.D0;
c     F
      DFGRD1(1,1)=1.D0
      DFGRD1(1,2)=0.D0
      DFGRD1(1,3)=0.D0
      DFGRD1(2,1)=1.D0
      DFGRD1(2,2)=1.D0
      DFGRD1(2,3)=0.D0
      DFGRD1(3,1)=0.D0
      DFGRD1(3,2)=0.D0
      DFGRD1(3,3)=1.D0

C----------------------------------------------------------------------
C     UMAT for compressible Neo-Hookean hyperelastity
C----------------------------------------------------------------------
      
C     elastic preproties
      C10  = PROPS(1)
      D1   = PROPS(2)
      mu_0 = C10*2.D0
      lambda_0 = 0.D0

C     calculate jacobian det(F)
      DET = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1    - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
      IF (NSHR .EQ. 3) THEN
          DET = DET + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     1              + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     2              - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     3              - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)
      END IF
      print *,"det"
      print *,DET
      
c     calculate distortion tensor (will be no use)
      SCALE = DET**(-ONE /THREE)
      print *,"SCALE"
      print *,SCALE
      
      print *,"DISTGR"
      DO K1 = 1, 3
        DO K2 = 1, 3
          DISTGR(K2,K1) = SCALE * DFGRD1(K2,K1)
          print *,DISTGR(K2,K1)
        END DO
      END DO
      
      print *,"DFGRD1"
      DO K1 = 1, 3
        DO K2 = 1, 3
          print *,DFGRD1(K2,K1)
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
      print *,"bbar"
      call Print_Array_1D(bbar,NTENS)

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
      print *,"bmat"
      call Print_Array_1D(bmat,NTENS)

C     calculate cauchy stress
      TRbbar = (bbar(1) + bbar(2) + bbar(3)) / THREE
      print *,"trbbar"
      print *,TRbbar

      EG     = TWO * C10 / DET
      EK     = TWO / D1 * (TWO * DET - ONE)
      PR     = TWO / D1 * (DET - ONE)
      
      EGXH   = mu_0 / DET
      PRXH   = lambda_0 * LOG(DET) / DET
      print *,"EGXH",EGXH
      print *,"EG",EG
c---------------------------------------------------------------------
      DO K1 = 1, NDI
        STRESS(K1) = EGXH * (bmat(K1) - TRBBAR) + PRXH
      END DO
      DO K1 = NDI+1, NDI+NSHR
        STRESS(K1) = EGXH * bmat(K1)
      END DO

      print *,"stress my"
      call Print_Array_1D(STRESS,NTENS)
c---------------------------------------------------------------------
      DO K1 = 1, NDI
        STRESS(K1) = EG * (BBAR(K1) - TRBBAR) + PR
      END DO
      DO K1 = NDI+1, NDI+NSHR
        STRESS(K1) = EG * BBAR(K1)
      END DO
            
      print *,"stress abaqus"
      call Print_Array_1D(STRESS,NTENS)
c---------------------------------------------------------------------

      

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

      end program main




      subroutine Print_Array_1D(Array_1D,Dim)
          integer Dim
          real(8) Array_1D(Dim)
          integer i
          Print *,(Array_1D(i),i=1,Dim)
      end

      subroutine Print_Array_2D(Array_2D,Row,Col)
          integer Row,Col
          real(8) Array_2D(Row,Col)
          integer i
          do i=1,Row
             Print *,(Array_2D(i,j),J=1,Col)  
          end do
      end