module model

  use const
  use mpi_setup, only : mpid

  implicit none
  private

  public :: model_setup


contains


  subroutine model_setup(ro,pr,vx,vy,vz,bx,by,bz,phi &
                        ,x,dx,xm,y,dy,ym,z,dz,zm,eta,min_dx)

!---Input & Output
  real(8),dimension(ix)      ,intent(out) :: x,dx
  real(8),dimension(0:ix)    ,intent(out) :: xm
  real(8),dimension(jx)      ,intent(out) :: y,dy
  real(8),dimension(0:jx)    ,intent(out) :: ym
  real(8),dimension(kx)      ,intent(out) :: z,dz
  real(8),dimension(0:kx)    ,intent(out) :: zm
  real(8),dimension(ix,jx,kx),intent(out) :: ro,pr
  real(8),dimension(ix,jx,kx),intent(out) :: vx,vy,vz
  real(8),dimension(ix,jx,kx),intent(out) :: bx,by,bz
  real(8),dimension(ix,jx,kx),intent(out) :: eta,phi
  real(8),                    intent(out) :: min_dx

  integer :: i,j,k
  integer :: ig,jg,kg
  integer :: izero,jzero,kzero
  real(8),dimension(igx) :: xg,dxg
  real(8),dimension(jgx) :: yg,dyg
  real(8),dimension(kgx) :: zg,dzg
  real(8),dimension(0:igx) :: xmg
  real(8),dimension(0:jgx) :: ymg
  real(8),dimension(0:kgx) :: zmg

  real(8) :: rr
  real(8) :: lg,sigp,sigt ! Lorentz Factor

! set global grid
  do i=1,igx
     dxg(i) = dxg0
  enddo
  do j=1,jgx
     dyg(j) = dyg0
  enddo
  do k=1,kgx
     dzg(k) = dzg0
  enddo

! XYZ position
  izero = margin+1
  xmg(izero) = xmin+dxg(izero)
  do i=izero,igx-1
     xmg(i+1) = xmg(i)+dxg(i+1)
  enddo
  do i=izero-1,0,-1
     xmg(i) = xmg(i+1)-dxg(i+1)
  enddo
  do i=1,igx
     xg(i) = 0.5d0*(xmg(i)+xmg(i-1))
  enddo

  jzero = margin+1
  ymg(jzero) = ymin+dyg(jzero)
  do j=jzero,jgx-1
     ymg(j+1) = ymg(j)+dyg(j+1)
  enddo
  do j=jzero-1,0,-1
     ymg(j) = ymg(j+1)-dyg(j+1)
  enddo
  do j=1,jgx
     yg(j) = 0.5d0*(ymg(j)+ymg(j-1))
  enddo

  kzero = margin+1
  zmg(kzero) = zmin+dzg(kzero)
  do k=kzero,kgx-1
     zmg(k+1) = zmg(k)+dzg(k+1)
  enddo
  do k=kzero-1,0,-1
     zmg(k) = zmg(k+1)-dzg(k+1)
  enddo
  do k=1,kgx
     zg(k) = 0.5d0*(zmg(k)+zmg(k-1))
  enddo

! set individual grid
  do i=1,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     x(i)=xg(ig)
     dx(i) = dxg(ig)
  enddo
  do j=1,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     y(j) = yg(jg)
     dy(j) = dyg(jg)
  enddo
  do k=1,kx
     kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
     z(k) = zg(kg)
     dz(k) = dzg(kg)
  enddo
  do i=0,ix
     ig = mpid%mpirank_3d(1)*(ix-2*margin)+i
     xm(i)=xmg(ig)
  enddo
  do j=0,jx
     jg = mpid%mpirank_3d(2)*(jx-2*margin)+j
     ym(j) = ymg(jg)
  enddo
  do k=0,kx
     kg = mpid%mpirank_3d(3)*(kx-2*margin)+k
     zm(k) = zmg(kg)
  enddo

! calculate min_dx
  min_dx = min(minval(dxg),minval(dyg),minval(dzg))

!----------------------------------------------------------------------|
  do k=1,kx
     do j=1,jx
        do i=1,ix
           rr = sqrt(x(i)*x(i)+y(j)*y(j)+z(k)*z(k))
           if( rr  < 0.1d0)then              
              ro(i,j,k) = 30d0
              pr(i,j,k) = 30d0
              vx(i,j,k) = 0.d0
              vy(i,j,k) = 0.d0
              vz(i,j,k) = 0.d0
              bx(i,j,k) = 5d0
              by(i,j,k) = 5d0
              bz(i,j,k) = 5d0
              phi(i,j,k) = 0.d0
              eta(i,j,k) = 0.d0

              lg = 1d0/dsqrt(1d0-(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))
              vx(i,j,k) = lg*vx(i,j,k)
              vy(i,j,k) = lg*vy(i,j,k)
              vz(i,j,k) = lg*vz(i,j,k)
           else
              ro(i,j,k) = 1d0
              pr(i,j,k) = 0.1d0
              vx(i,j,k) = 0.d0
              vy(i,j,k) = 0.d0
              vz(i,j,k) = 0.d0
              bx(i,j,k) = 5d0
              by(i,j,k) = 5d0
              bz(i,j,k) = 5d0
              phi(i,j,k) = 0.d0
              eta(i,j,k) = 0.d0

              lg = 1d0/dsqrt(1d0-(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))
              vx(i,j,k) = lg*vx(i,j,k)
              vy(i,j,k) = lg*vy(i,j,k)
              vz(i,j,k) = lg*vz(i,j,k)
           endif

        enddo
     enddo
  enddo

end subroutine model_setup


end module model
