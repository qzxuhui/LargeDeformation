c     程序说明：
c     各项同性等温弹性 UMAT
c     不能被用于平面应力问题
c     徐辉 20180517 于 清华大学
c     材料说明：
c     props(1) - E    弹性模量
c     props(2) - NU   泊松比
c     -----------------------------
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)

      INCLUDE 'ABA_PARAM.INC'

      CHARACTER*8 CMNAME

      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3)
     
      DIMENSION EELAS(6), EPLAS(6), FLOW(6)
      

      real E,nu
      real G,lambda
      integer i,j

c-----弹性常数计算----------------------------------
c     弹性模量
      E=props(1)
c     泊松比
      nu=props(2)
c     拉梅常数1：G剪切模量
      G=E/(2.D0*(1+nu))
c     拉梅常数2：lambda
      lambda=(E*nu)/(1.D0+nu)/(1-2.D0*nu)

c-----计算弹性刚度----------------------------------
      do i=1, NDI
          do j=1, NDI
              DDSDDE(j, i)=lambda
          end do
          DDSDDE(i, i)=2.D0*G+lambda
      end do
      do i=NDI+1, NTENS
          DDSDDE(i ,i)=G
      end do
      
c-----计算应力----------------------------------
c     sigma_new = sigma + dsigma
c     dsigma = ddsdde * dstrain
      do i=1, NTENS
          do j=1, NTENS
              STRESS(i)=STRESS(i)+DDSDDE(i, j)*DSTRAN(j)
          end do
      end do
      
      return
      end