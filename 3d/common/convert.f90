module convert

  implicit none
  private

  public :: convert__ptoc, convert__ctop


contains


  subroutine convert__ptoc(ix,jx,kx,gm,ro,pr,vx,vy,vz,bxc,byc,bzc &
                          ,dro,rx,ry,rz,ee)
!======================================================================
! Name :: convert_ptoc
!         convert to conserved
! Input ::
!          ix,jx,kx :: array size
!          gm :: specific heat retio
!          ro,pr,vx,vy,vz,bxc,byc,bzc
!             :: primitice variables
! Output ::
!          dro,rx,ry,rz,ee :: conserved varialbles
!
!======================================================================
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm
  real(8),dimension(ix,jx,kx),intent(in) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(in) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(in) :: bxc,byc,bzc
  real(8),dimension(ix,jx,kx),intent(out) :: dro,rx,ry,rz,ee

  integer :: i,j,k
  real(8) :: vsq,vb,lg,hh,pb,ilg ! lg : Lorentz factor hh : specific enthalpy

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j,vsq,lg,ilg,hh,pb,vb)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           vsq = vx(i,j,k)*vx(i,j,k) + vy(i,j,k)*vy(i,j,k) + vz(i,j,k)*vz(i,j,k)
           lg  = dsqrt(1d0 + vsq) ! Lorentz Factor
           ilg = 1d0/lg

           hh = 1d0 + gm*pr(i,j,k)/(ro(i,j,k)*(gm-1d0))

           dro(i,j,k) = lg*ro(i,j,k)
           pb = 0.5d0*(bxc(i,j,k)**2 + byc(i,j,k)**2 + bzc(i,j,k)**2) ! Magnetic Pressure

           vb = ( vx(i,j,k)*bxc(i,j,k) + vy(i,j,k)*byc(i,j,k) + vz(i,j,k)*bzc(i,j,k) )*ilg
           rx(i,j,k) = (dro(i,j,k)*lg*hh + 2d0*pb)*vx(i,j,k)*ilg - vb*bxc(i,j,k)
           ry(i,j,k) = (dro(i,j,k)*lg*hh + 2d0*pb)*vy(i,j,k)*ilg - vb*byc(i,j,k)
           rz(i,j,k) = (dro(i,j,k)*lg*hh + 2d0*pb)*vz(i,j,k)*ilg - vb*bzc(i,j,k)

           ee(i,j,k) =dro(i,j,k)*hh*lg - pr(i,j,k) + pb + 0.5d0*(vsq*ilg*ilg*2d0*pb -vb*vb )

           !! *** HD case
           ! vsq = vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2
           ! lg = 1d0/dsqrt(1.0d0-vsq)
           !
           ! hh =  1d0 + gm*pr(i,j,k)/(ro(i,j,k)*(gm-1d0))
           !
           ! dro(i,j,k) = lg*ro(i,j,k)
           !
           ! rx(i,j,k) = dro(i,j,k)*hh*lg*vx(i,j,k)
           ! ry(i,j,k) = dro(i,j,k)*hh*lg*vy(i,j,k)
           ! rz(i,j,k) = dro(i,j,k)*hh*lg*vz(i,j,k)
           !
           ! pb = 0.5d0*(bxc(i,j,k)**2 + byc(i,j,k)**2 + bzc(i,j,k)**2)
           ! vsq = vx(i,j,k)**2 + vy(i,j,k)**2 + vz(i,j,k)**2
           !
           ! ee(i,j,k) = dro(i,j,k)*hh*lg - pr(i,j,k)
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

  end subroutine convert__ptoc


  subroutine convert__ctop(ix,jx,kx,gm,ro,dro,ee,rx,ry,rz,bxc,byc,bzc &
                          ,vx,vy,vz,pr)
!======================================================================
! Name :: convert_ctop
!         convert to primitive
! Input ::
!          ix,jx,kx :: array size
!          gm :: specific heat retio
!          ro,ee,rx,ry,rz,bxc,byc,bzc
!             :: conserved variables
! Output ::
!           vx,vy,vz,pr :: primitive variable
!======================================================================
  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: gm
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(inout) :: ee
  real(8),dimension(ix,jx,kx),intent(in) :: dro,rx,ry,rz
  real(8),dimension(ix,jx,kx),intent(in) :: bxc,byc,bzc
  real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz

  integer :: i,j,k,count
  real(8), parameter :: pbeta_min=1d-3
  real(8) :: vsq,pb,roinverse,igm,temppr,signpr,temp1,temp2

  ! aloy et al.(1999) for HD
  !real(8) :: ps,vxs,vys,vzs,vvs,gls,eps,ros,css2,fp,dfp

  ! Mignone and Bodo et al. (2006) for MHD
  real(8) :: lg,ilg,ww,ss,bsq,rsq,pg,dlg,dpg,fw,dfw
  real(8) :: bb,cc,wwm,wwp

  real(8), parameter :: emin = 1d-10
  integer, parameter :: cmax = 100

  igm = 1d0/(gm-1d0)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j,count,ss,bsq,rsq,bb,cc,wwp,ww,lg,pg,dlg,dpg,fw,dfw)

  !***
  ! Primitve Ricovery for MHD case : Mignone and Bodo (2006)
  !***
  do k=1,kx
    do j=1,jx
      do i=1,ix
        ss  = rx(i,j,k)*bxc(i,j,k) + ry(i,j,k)*byc(i,j,k) + rz(i,j,k)*bzc(i,j,k)
        bsq = bxc(i,j,k)**2 + byc(i,j,k)**2 + bzc(i,j,k)**2
        rsq = rx(i,j,k)*rx(i,j,k)+ry(i,j,k)*ry(i,j,k)+rz(i,j,k)*rz(i,j,k)

        !*** first guess using previous stepes val. ***!
        ! vsq = vx(i,j,k)*vx(i,j,k) + vy(i,j,k)*vy(i,j,k) + vz(i,j,k)*vz(i,j,k)
        ! lg  = dsqrt(1d0 + vsq) ! Lorentz Factor
        ! ww = ( ro(i,j,k) + igm*gm*pr(i,j,k) )*lg*lg
        !
        ! ww = dsqrt(rsq+1d-8)
        ! ww = ww + ww*1d-4
        !**********************************************!

        !*** first guess using Mignone and McKinney ***!
        bb = 4.d0 *( bsq - ee(i,j,k) )
        cc = bsq*bsq + rsq - 2d0*bsq*ee(i,j,k)
        wwp = (-bb + dsqrt(bb*bb - 4d0*3d0*cc))/6d0
        ww  = wwp
        !***********************************************!
        do count=1,cmax
            lg  = 1d0 - (ss*ss*(bsq+ 2d0*ww)+ww*ww*rsq)/(ww*ww*(ww+bsq)*(ww+bsq))
            lg  = 1d0/dsqrt(lg)

            pg  = (ww-dro(i,j,k)*lg)/(igm*gm*lg*lg)

            dlg = lg*lg*lg*(rsq*ww*ww*ww+3d0*ss*ss*ww*ww+3d0*ss*ss*bsq*ww+ss*ss*bsq*bsq)
            dlg = -1d0*dlg/(ww*ww*ww*(ww+bsq)*(ww+bsq)*(ww+bsq))

            dpg = (lg*(1d0+dro(i,j,k)*dlg)-2d0*ww*dlg)/(gm*igm*lg*lg*lg)

            fw  = ww - pg + (1d0-0.5d0/(lg*lg))*bsq - 0.5d0*ss*ss/(ww*ww) - ee(i,j,k)

!            write(*,*) count, dabs((pg-pr(i,j,k))/pg),pr(i,j,k),pg
           if( dabs((pg-pr(i,j,k))/pg) .lt. emin)then
!            if( dabs(fw) .lt. emin)then
              exit
            endif
            dfw = 1d0 - dpg + dlg*bsq/(lg*lg*lg) + ss*ss/(ww*ww*ww)
            ww  = ww - fw/dfw
            pr(i,j,k) = pg
        end do

        if(count .eq. cmax+1)then
          write(*,*)"Stop Due to faile for Primitive recovery"
          write(*,*) i,j,k,pg,lg
          stop
        endif

        vx(i,j,k) = lg*(rx(i,j,k) + bxc(i,j,k)*ss/ww)/(ww+bsq)
        vy(i,j,k) = lg*(ry(i,j,k) + byc(i,j,k)*ss/ww)/(ww+bsq)
        vz(i,j,k) = lg*(rz(i,j,k) + bzc(i,j,k)*ss/ww)/(ww+bsq)
        pr(i,j,k) = pg
        ro(i,j,k) = dro(i,j,k)/lg
      enddo
    enddo
  enddo


  !***
  ! Primitive Ricovery for HD case : Aloy et al. (1999)
  !***
  ! do k=1,kx
  !    do j=1,jx
  !       do i=1,ix
  !           ! first guess for pressure
  !           ps = pr(i,j,k)
  !           do count=1,cmax
  !             vxs = rx(i,j,k)/(ee(i,j,k)+ps)
  !             vys = ry(i,j,k)/(ee(i,j,k)+ps)
  !             vzs = rz(i,j,k)/(ee(i,j,k)+ps)
  !             vvs = vxs*vxs+vys*vys+vzs*vzs
  !             gls = 1d0/dsqrt(1d0-vvs)
  !             eps = (ee(i,j,k)-dro(i,j,k)+dro(i,j,k)*(1d0-gls)+ps*(1d0-gls*gls))/(dro(i,j,k)*gls)
  !             ros = dro(i,j,k)/gls
  !
  !             fp = (gm-1d0)*ros*eps-ps
  !
  !             if( dabs(fp) .lt. emin)then
  !               pr(i,j,k) = ps
  !               exit
  !             endif
  !
  !             css2 = gm*(gm-1d0)*eps/(1+gm*eps)
  !             dfp = vvs*css2 - 1d0
  !
  !             ps = ps - fp/dfp
  !           enddo
  !
  !           if(count .eq. cmax)then
  !             write(*,*)"Stop Due to faile for Primitive recovery"
  !           endif
  !
  !           vx(i,j,k) = vxs*gls
  !           vy(i,j,k) = vys*gls
  !           vz(i,j,k) = vzs*gls
  !           ro(i,j,k) = ros
  !
  !       enddo
  !    enddo
  ! enddo

  !$OMP END PARALLEL DO
  end subroutine convert__ctop


end module convert
