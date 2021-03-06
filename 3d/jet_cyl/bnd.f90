module bnd

  implicit none
  private

  public :: bnd__exec
  
  
contains


  subroutine bnd__exec(margin,ix,jx,kx,x,y,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)
    
    use mpi_setup, only : mpid, mnull
    use boundary
    use const, only : r_jet,zin
    
    integer,intent(in) :: margin,ix,jx,kx
    real(8),dimension(ix),intent(in) :: x
    real(8),dimension(jx),intent(in) :: y    
    real(8),dimension(ix,jx,kx),intent(inout) :: ro,pr
    real(8),dimension(ix,jx,kx),intent(inout) :: vx,vy,vz
    real(8),dimension(ix,jx,kx),intent(inout) :: bx,by,bz
    real(8),dimension(ix,jx,kx),intent(inout) :: phi,eta
    integer :: i,j,k
    real(8) :: lg
    
!======================================================================
! inter-process communication by MPI
    call boundary__mpi(margin,ix,jx,kx,ro,pr,vx,vy,vz,bx,by,bz,phi,eta)

!----------------------------------------------------------------------|
    if(mpid%l == mnull) then
       ! call bd_frex(0,margin,ro,ix,jx,kx)
       ! call bd_frex(0,margin,pr,ix,jx,kx)
       ! call bd_frex(0,margin,vx,ix,jx,kx)
       ! call bd_frex(0,margin,vy,ix,jx,kx)
       ! call bd_frex(0,margin,vz,ix,jx,kx)
       ! call bd_frex(0,margin,bx,ix,jx,kx)
       ! call bd_frex(0,margin,by,ix,jx,kx)
       ! call bd_frex(0,margin,bz,ix,jx,kx)
       ! call bd_frex(0,margin,phi,ix,jx,kx)
       ! call bd_frex(0,margin,eta,ix,jx,kx)
       call bd_synpx_car(0,margin,ro,ix,jx,kx)
       call bd_synpx_car(0,margin,pr,ix,jx,kx)
       call bd_synnx_car(0,margin,vx,ix,jx,kx)
       call bd_synnx_car(0,margin,vy,ix,jx,kx)
       call bd_synpx_car(0,margin,vz,ix,jx,kx)
       call bd_synnx_car(0,margin,bx,ix,jx,kx)
       call bd_synnx_car(0,margin,by,ix,jx,kx)
       call bd_synpx_car(0,margin,bz,ix,jx,kx)
       call bd_synpx_car(0,margin,phi,ix,jx,kx)
       call bd_synpx_car(0,margin,eta,ix,jx,kx)       
    end if
    if(mpid%r == mnull) then
       call bd_frex(1,margin,ro,ix,jx,kx)
       call bd_frex(1,margin,pr,ix,jx,kx)
       call bd_frex(1,margin,vx,ix,jx,kx)
       call bd_frex(1,margin,vy,ix,jx,kx)
       call bd_frex(1,margin,vz,ix,jx,kx)
       call bd_frex(1,margin,bx,ix,jx,kx)
       call bd_frex(1,margin,by,ix,jx,kx)
       call bd_frex(1,margin,bz,ix,jx,kx)
       call bd_frex(1,margin,phi,ix,jx,kx)
       call bd_frex(1,margin,eta,ix,jx,kx)
    end if

!----------------------------------------------------------------------|
! free boundary in y
    if(mpid%b == mnull) then
       call bd_frey(0,margin,ro,ix,jx,kx)
       call bd_frey(0,margin,pr,ix,jx,kx)
       call bd_frey(0,margin,vx,ix,jx,kx)
       call bd_frey(0,margin,vy,ix,jx,kx)
       call bd_frey(0,margin,vz,ix,jx,kx)
       call bd_frey(0,margin,bx,ix,jx,kx)
       call bd_frey(0,margin,by,ix,jx,kx)
       call bd_frey(0,margin,bz,ix,jx,kx)
       call bd_frey(0,margin,phi,ix,jx,kx)
       call bd_frey(0,margin,eta,ix,jx,kx)
    end if
    if(mpid%f == mnull) then
       call bd_frey(1,margin,ro,ix,jx,kx)
       call bd_frey(1,margin,pr,ix,jx,kx)
       call bd_frey(1,margin,vx,ix,jx,kx)
       call bd_frey(1,margin,vy,ix,jx,kx)
       call bd_frey(1,margin,vz,ix,jx,kx)
       call bd_frey(1,margin,bx,ix,jx,kx)
       call bd_frey(1,margin,by,ix,jx,kx)
       call bd_frey(1,margin,bz,ix,jx,kx)
       call bd_frey(1,margin,phi,ix,jx,kx)
       call bd_frey(1,margin,eta,ix,jx,kx)
    end if

!----------------------------------------------------------------------|
! free boundary in z
    if(mpid%d == mnull) then
       ! call bd_synpz_car(0,margin,ro,ix,jx,kx)
       ! call bd_synpz_car(0,margin,pr,ix,jx,kx)
       ! call bd_synpz_car(0,margin,vx,ix,jx,kx)
       ! call bd_synpz_car(0,margin,vy,ix,jx,kx)
       ! call bd_synnz_car(0,margin,vz,ix,jx,kx)
       ! call bd_synnz_car(0,margin,bx,ix,jx,kx)
       ! call bd_synnz_car(0,margin,by,ix,jx,kx)
       ! call bd_synpz_car(0,margin,bz,ix,jx,kx)
       ! call bd_synpz_car(0,margin,phi,ix,jx,kx)
       ! call bd_synpz_car(0,margin,eta,ix,jx,kx)
       call bd_frez(0,margin,ro,ix,jx,kx)
       call bd_frez(0,margin,pr,ix,jx,kx)
       call bd_frez(0,margin,vx,ix,jx,kx)
       call bd_frez(0,margin,vy,ix,jx,kx)
       call bd_frez(0,margin,vz,ix,jx,kx)
       call bd_frez(0,margin,bx,ix,jx,kx)
       call bd_frez(0,margin,by,ix,jx,kx)
       call bd_frez(0,margin,bz,ix,jx,kx)
       call bd_frez(0,margin,phi,ix,jx,kx)
       call bd_frez(0,margin,eta,ix,jx,kx)
    end if
    if(mpid%t == mnull) then
       call bd_frez(1,margin,ro,ix,jx,kx)
       call bd_frez(1,margin,pr,ix,jx,kx)
       call bd_frez(1,margin,vx,ix,jx,kx)
       call bd_frez(1,margin,vy,ix,jx,kx)
       call bd_frez(1,margin,vz,ix,jx,kx)
       call bd_frez(1,margin,bx,ix,jx,kx)
       call bd_frez(1,margin,by,ix,jx,kx)
       call bd_frez(1,margin,bz,ix,jx,kx)
       call bd_frez(1,margin,phi,ix,jx,kx)
       call bd_frez(1,margin,eta,ix,jx,kx)
    end if
    
    do k=1,margin
       !$OMP PARALLEL DO &
       !$OMP PRIVATE(i)
       do j=1,jx
          do i=1,ix
             if (x(i).le.r_jet)then
                ro(i,j,k) = 1d0
                pr(i,j,k) = 0.015d0
                vx(i,j,k) = 0d0
                vy(i,j,k) = 0d0
                vz(i,j,k) = 0.995d0
                bx(i,j,k) = 0.0d0
                by(i,j,k) = 0.001d0
                bz(i,j,k) = 0.0d0
                phi(i,j,k) = 0.d0
                eta(i,j,k) = 0.d0

                lg = 1d0/dsqrt(1d0-(vx(i,j,k)*vx(i,j,k)+vy(i,j,k)*vy(i,j,k)+vz(i,j,k)*vz(i,j,k)))
                vx(i,j,k) = lg*vx(i,j,k)
                vy(i,j,k) = lg*vy(i,j,k)
                vz(i,j,k) = lg*vz(i,j,k)
             endif
          enddo
       enddo
          !$OMP END PARALLEL DO
    enddo
    
  end subroutine bnd__exec


end module bnd
