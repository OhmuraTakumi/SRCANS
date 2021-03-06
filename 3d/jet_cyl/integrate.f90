module integrate

  implicit none
  private

  public :: integrate__TVDRK3


contains


  subroutine integrate__TVDRK3(margin,ix,jx,kx,gm,x,xm,y,ym,dx,dy,dz,dt          &
                              ,ro,pr,vx,vy,vz,bx,by,bz,phi,ch,cp &
                              ,eta,ccx,ccy,ccz)

  use convert
  use lr_state, only : lr_state__MP5, lr_state__MSCL2, lr_state__1st
  use flux_calc
  use bnd

  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: ch,cp
  real(8),intent(in) :: dt,gm
  real(8),dimension(ix),intent(in) :: x,dx
  real(8),dimension(jx),intent(in) :: y,dy
  real(8),dimension(kx),intent(in) :: dz
  real(8),dimension(0:ix),intent(in) :: xm
  real(8),dimension(0:jx),intent(in) :: ym
  real(8),dimension(5,2,ix),intent(in) :: ccx
  real(8),dimension(5,2,jx),intent(in) :: ccy
  real(8),dimension(5,2,kx),intent(in) :: ccz
  real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr,vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(inout) :: phi
  real(8),dimension(ix,jx,kx),intent(inout) :: eta

! !-- using flux
!   real(8),dimension(ix,jx,kx) :: ro1,pr1,vx1,vy1,vz1
!   real(8),dimension(ix,jx,kx) :: bx1,by1,bz1
!   real(8),dimension(ix,jx,kx) :: phi1
! !-conserved variable
!   real(8),dimension(ix,jx,kx) :: dro,rx,ry,rz,ee
!   real(8),dimension(ix,jx,kx) :: dro1,rx1,ry1,rz1,ee1
! !-surface variables
!   real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
!   real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw
!   real(8),dimension(ix,jx,kx) :: bx_m,by_m,bz_m,phi_m
! !-Numerical flux
! !x-component
!   real(8),dimension(ix,jx,kx) :: frox,feex,frxx,fryx,frzx
!   real(8),dimension(ix,jx,kx) :: fbyx,fbzx,fbxx,fphix
!
! !-other temporary variables
!   integer :: mdir
!   integer :: i,j,k,n
!   real(8), parameter :: fac=1.D0/12.D0
!   real(8) :: dtodx,k1,k2
!-- using flux
  real(8),dimension(ix,jx,kx) :: ro1,pr1,vx1,vy1,vz1
  real(8),dimension(ix,jx,kx) :: bx1,by1,bz1
  real(8),dimension(ix,jx,kx) :: phi1
!-conserved variable
  real(8),dimension(ix,jx,kx) :: dro,rx,ry,rz,ee
  real(8),dimension(ix,jx,kx) :: dro1,rx1,ry1,rz1,ee1
!-surface variables
  real(8),dimension(ix,jx,kx,2) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2) :: bxw,byw,bzw,phiw
  real(8),dimension(ix,jx,kx) :: bx_m,by_m,bz_m,phi_m
!-Numerical flux
!x-component
  real(8),dimension(ix,jx,kx) :: frox,feex,frxx,fryx,frzx
  real(8),dimension(ix,jx,kx) :: fbyx,fbzx,fbxx,fphix
  real(8),dimension(ix,jx,kx) :: fbyxr,fbzxr,fbxxr,feexr
!y-component
  real(8),dimension(ix,jx,kx) :: froy,feey,frxy,fryy,frzy
  real(8),dimension(ix,jx,kx) :: fbyy,fbzy,fbxy,fphiy
  real(8),dimension(ix,jx,kx) :: fbyyr,fbzyr,fbxyr,feeyr
!z-component
  real(8),dimension(ix,jx,kx) :: froz,feez,frxz,fryz,frzz
  real(8),dimension(ix,jx,kx) :: fbxz,fbyz,fbzz,fphiz
  real(8),dimension(ix,jx,kx) :: fbxzr,fbyzr,fbzzr,feezr

!-other temporary variables
  integer :: mdir
  integer :: i,j,k,n
  real(8), parameter :: fac=1.D0/12.D0
  real(8) :: dtodx,dtody,dtodz,k1,k2
  real(8),dimension(ix,jx,kx) :: curx,cury,curz

!-Source Term
  real(8) :: srmx,srmy
  real(8) :: lgl, bby,bbz
  
  
!-----Step 0.----------------------------------------------------------|
! primitive to conserve
  call convert__ptoc(ix,jx,kx,gm,ro,pr,vx,vy,vz,bx,by,bz &
                    ,dro,rx,ry,rz,ee)
  dro1=dro
  rx1=rx
  ry1=ry
  rz1=rz
  bx1=bx
  by1=by
  bz1=bz
  ee1=ee
  phi1=phi

  do n=1,3

!-----Step 1a.---------------------------------------------------------|
! Compute flux in x-direction
! set L/R state at x-direction
!
  mdir = 1

  ! call lr_state__MP5(mdir,ix,jx,kx,ro,pr &
  !       ,vx,vy,vz,bx,by,bz,phi &
  !       ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,ccx,ccy,ccz,dx,dy,dz)
  call lr_state__MSCL2(mdir,ix,jx,kx,ro,pr &
      ,vx,vy,vz,bx,by,bz,phi &
      ,ch,gm,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw,dx,dy,dz)
  ! call lr_state__1st(mdir,ix,jx,kx,ro,pr &
  !     ,vx,vy,vz,bx,by,bz,phi &
  !     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)

  call flux_calc__bp(ix,jx,kx,bxw,phiw &
       ,bx_m,phi_m,ch)
  call flux_calc__glm(bx_m,phi_m,ch,fbxx,fphix,ix,jx,kx)

  ! call flux_calc__hll(row,prw,vxw,vyw,vzw,bx_m,byw,bzw,gm,margin,ix,jx,kx &
  !                    ,frox,feex,frxx,fryx,frzx,fbyx,fbzx)
  call flux_calc__hllc(row,prw,vxw,vyw,vzw,bx_m,byw,bzw,gm,margin,ix,jx,kx &
                   ,frox,feex,frxx,fryx,frzx,fbyx,fbzx)

  ! mdir = 2
  ! ! call lr_state__MP5(mdir,ix,jx,kx,ro,pr &
  ! !      ,vy,vz,vx,by,bz,bx,phi &
  ! !      ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,ccx,ccy,ccz)
  ! call lr_state__MSCL2(mdir,ix,jx,kx,ro,pr &
  !     ,vy,vz,vx,by,bz,bx,phi &
  !     ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,dx,dy,dz)
  ! ! call lr_state__1st(mdir,ix,jx,kx &
  ! !     ,ro,pr,vx,vy,vz,bx,by,bz,phi &
  ! !     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)
  
  ! call flux_calc__bp(ix,jx,kx,byw,phiw &
  !      ,by_m,phi_m,ch)
  
  ! call flux_calc__glm(by_m,phi_m,ch,fbyy,fphiy,ix,jx,kx)
  
  ! call flux_calc__hll(row,prw,vyw,vzw,vxw,by_m,bzw,bxw,gm,margin,ix,jx,kx &
  !      ,froy,feey,fryy,frzy,frxy,fbzy,fbxy)
  ! call flux_calc__hllc(row,prw,vyw,vzw,vxw,by_m,bzw,bxw,gm,margin,ix,jx,kx &
  !      ,froy,feey,fryy,frzy,frxy,fbzy,fbxy)

  mdir = 3
  ! call lr_state__MP5(mdir,ix,jx,kx,ro,pr &
  !      ,vy,vz,vx,by,bz,bx,phi &
  !      ,ch,gm,row,prw,vyw,vzw,vxw,byw,bzw,bxw,phiw,ccx,ccy,ccz)
  call lr_state__MSCL2(mdir,ix,jx,kx,ro,pr &
      ,vz,vx,vy,bz,bx,by,phi &
      ,ch,gm,row,prw,vzw,vxw,vyw,bzw,bxw,byw,phiw,dx,dy,dz)
  ! call lr_state__1st(mdir,ix,jx,kx &
  !     ,ro,pr,vx,vy,vz,bx,by,bz,phi &
  !     ,row,prw,vxw,vyw,vzw,bxw,byw,bzw,phiw)
  
  call flux_calc__bp(ix,jx,kx,bzw,phiw &
       ,bz_m,phi_m,ch)
  
  call flux_calc__glm(bz_m,phi_m,ch,fbzz,fphiz,ix,jx,kx)
  
  ! call flux_calc__hll(row,prw,vzw,vxw,vyw,bz_m,bxw,byw,gm,margin,ix,jx,kx &
  !      ,froz,feez,frzz,frxz,fryz,fbxz,fbyz)
  call flux_calc__hllc(row,prw,vzw,vxw,vyw,bz_m,bxw,byw,gm,margin,ix,jx,kx &
                   ,froz,feez,frzz,frxz,fryz,fbxz,fbyz)    

!-----Step 2.---------------------------------------------------------|
! TVDRK substep
  k1 = fac*(-7.D0*n*n+30.D0*n-23.D0)
  k2 = fac*(+7.D0*n*n-30.D0*n+35.D0)
  do k=margin+1,kx-margin
     do j=margin+1,jx-margin
        do i=margin+1,ix-margin

           dtodx = dt/(x(i)*dx(i))
           dtodz = dt/dz(k)

           lgl = dsqrt(1d0+vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k))
           bby = by(i,j,k)/(lgl*lgl)+ &
                (vx(i,j,k)*bx(i,j,k)+vy(i,j,k)*by(i,j,k)+vz(i,j,k)*bz(i,j,k))*vy(i,j,k)/(lgl*lgl)
           ! bbz = bz(i,j,k)/(lgl*lgl)+ &
           !      (vx(i,j,k)*bx(i,j,k)+vy(i,j,k)*by(i,j,k)+vz(i,j,k)*bz(i,j,k))*vz(i,j,k)/(lgl*lgl)
           
           srmx = pr(i,j,k)/x(i) + (ry(i,j,k)*vy(i,j,k))/(lgl*x(i)) &
                -bby*by(i,j,k)/x(i) 
           srmy = 0d0
           
           dro(i,j,k) = k1*dro1(i,j,k)+k2*(+dro(i,j,k) &
                +dtodx*(xm(i-1)*frox(i-1,j,k)-xm(i)*frox(i,j,k))  &
                +dtodz*(froz(i,j,k-1)-froz(i,j,k))  &
                )
           ee(i,j,k) = k1*ee1(i,j,k)+k2*(+ee(i,j,k)  &
                +dtodx*(xm(i-1)*feex(i-1,j,k)-xm(i)*feex(i,j,k)) &
                +dtodz*(feez(i,j,k-1)-feez(i,j,k)) &
                )
           rx(i,j,k) = k1*rx1(i,j,k)+k2*(+rx(i,j,k) &
                +dtodx*(xm(i-1)*frxx(i-1,j,k)-xm(i)*frxx(i,j,k))  &
                +dtodz*(frxz(i,j,k-1)-frxz(i,j,k))  &
                +dt*srmx )
           ! ry(i,j,k) = k1*ry1(i,j,k)+k2*(+ry(i,j,k) &
           !      +dtodx*(xm(i-1)*fryx(i-1,j,k)-xm(i)*fryx(i,j,k))  &
           !      +dtodz*(fryz(i,j,k-1)-fryz(i,j,k))  &
           !      +dt*srmy )
           rz(i,j,k) = k1*rz1(i,j,k)+k2*(+rz(i,j,k) &
                +dtodx*(xm(i-1)*frzx(i-1,j,k)-xm(i)*frzx(i,j,k))  &
                +dtodz*(frzz(i,j,k-1)-frzz(i,j,k))  &
                 )
           bx(i,j,k) = k1*bx1(i,j,k)+k2*(+bx(i,j,k)  &
                +x(i)*dtodx*(fbxx(i-1,j,k)-fbxx(i,j,k))   &
                +dtodz*(fbxz(i,j,k-1)-fbxz(i,j,k)) &
                 )
           by(i,j,k) = k1*by1(i,j,k)+k2*(+by(i,j,k)  &
                +x(i)*dtodx*(fbyx(i-1,j,k)-fbyx(i,j,k)) &
                +dtodz*(fbyz(i,j,k-1)-fbyz(i,j,k))   &
                )
           bz(i,j,k) = k1*bz1(i,j,k)+k2*(+bz(i,j,k)  &
                +dtodx*(xm(i-1)*fbzx(i-1,j,k)-xm(i)*fbzx(i,j,k)) &
                +dtodz*(fbzz(i,j,k-1)-fbzz(i,j,k)) &
                )

           phi(i,j,k) = k1*phi1(i,j,k)+k2*(+phi(i,j,k) &
                +dtodx*(xm(i-1)*fphix(i-1,j,k)-xm(i)*fphix(i,j,k))   &
                +dtodz*(fphiz(i,j,k-1)-fphiz(i,j,k))   &
                )*exp(-dt*ch**2/cp**2)
        enddo
     enddo
  enddo

!-----Step 3.----------------------------------------------------------|
! conserved to primitive
!
  call convert__ctop(ix,jx,kx,gm,ro,dro,ee,rx,ry,rz,bx,by,bz &
                    ,vx,vy,vz,pr)
  call bnd__exec(margin,ix,jx,kx,x,y,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

  enddo

  end subroutine integrate__TVDRK3


end module integrate
