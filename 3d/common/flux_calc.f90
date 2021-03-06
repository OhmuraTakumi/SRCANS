module flux_calc

  implicit none
  private

  public :: flux_calc__hll, flux_calc__hllc, flux_calc__glm, flux_calc__bp, &
            flux_calc__fbres, flux_calc__feres

  real(8), parameter :: eps=1d-40


contains

  subroutine flux_calc__hll(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,margin,ix,jx,kx &
                           ,fro,fee,frx,fry,frz,fby,fbz)

  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: gm
! primitive variables :: 1 = left state , 2 = right state
  real(8),dimension(ix,jx,kx,2),intent(in) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(in) :: byw,bzw
  real(8),dimension(ix,jx,kx),intent(in) :: bx ! magnetic field
  real(8),dimension(ix,jx,kx),intent(out) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx),intent(out) :: fby,fbz

!----- U -----
! qql :: left state
! qqr :: right state
  real(8) :: rol,vxl,vyl,vzl,byl,bzl,ptl,eel
  real(8) :: ror,vxr,vyr,vzr,byr,bzr,ptr,eer
  real(8) :: rxl,ryl,rzl
  real(8) :: rxr,ryr,rzr
  real(8) :: bxs,bxsq
  real(8) :: pbl,pbr,prl,prr
  real(8) :: gmpl,gmpr,gpbl,gpbr
  real(8) :: vvl,lgl,hhl,drol
  real(8) :: vvr,lgr,hhr,dror
  real(8) :: vb

!----- convariant mag field
  real(8) :: bbl,b0xl,b0yl,b0zl,bbr,b0xr,b0yr,b0zr

!----- Wave Estimate
  real(8) :: lg,ilg,roh,cs2,bb,sig,hs,ca2,om2,rr,vsq,l1,l2,l3,lfl,lfr,lsl,lsr
!----- U* ----
! qqlst :: left state
! qqrst :: right state
  real(8) :: sl,sr
!----- flux ---
! fqql :: left physical flux
! fqqr :: right physical flux
! fluxqq :: intermediate HLLD flux (OUTPUT)
  real(8) :: frol,frxl,fryl,frzl
  real(8) :: fbyl,fbzl,feel
  real(8) :: fror,frxr,fryr,frzr
  real(8) :: fbyr,fbzr,feer
  real(8) :: bsql,bsqr
  real(8) :: cfl,cfr
  integer :: i,j,k

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j) &
  !$OMP PRIVATE(bxs,bxsq) &
  ! --- l state ---
  !$OMP PRIVATE(rol,vxl,vyl,vzl,byl,bzl,bsql,pbl,prl,ptl,eel,rxl,ryl,rzl) &
  !$OMP PRIVATE(vvl,lgl,hhl,drol,ilg,vb) &  
  ! --- r state ---
  !$OMP PRIVATE(ror,vxr,vyr,vzr,byr,bzr,bsqr,pbr,prr,ptr,eer,rxr,ryr,rzr) &
  !$OMP PRIVATE(vvr,lgr,hhr,dror) &    
  !--- step 1----
  !$OMP PRIVATE(gmpl,gmpr,gpbl,gpbr,cfl,cfr,sl,sr) &
  !$OMP PRIVATE(roh,cs2,bb,sig,hs,ca2,om2,rr,vsq,l1,l2,l3,lfl,lsl,lfr,lsr) &      
  !--- step 2 ---
  ! - left
  !$OMP PRIVATE(frol,frxl,fryl,frzl,feel,fbyl,fbzl) &
  !$OMP PRIVATE(bbl,b0xl,b0yl,b0zl) &  
  ! - right
  !$OMP PRIVATE(fror,frxr,fryr,frzr,feer,fbyr,fbzr)&
  !$OMP PRIVATE(bbr,b0xr,b0yr,b0zr)   
  do k=margin,kx-margin
     do j=margin,jx-margin
        do i=margin,ix-margin
!----- Step 0. ----------------------------------------------------------|
! set L/R-state
!
           bxs = bx(i,j,k)
           bxsq = bxs**2
!---- Left state

           rol = row(i,j,k,1)
           prl = prw(i,j,k,1)
           vxl = vxw(i,j,k,1)
           vyl = vyw(i,j,k,1)
           vzl = vzw(i,j,k,1)
           byl = byw(i,j,k,1)
           bzl = bzw(i,j,k,1)

           bsql = bxs**2+byl**2+bzl**2
           pbl = 0.5d0*(bxs**2 + byl**2 + bzl**2)
           ptl = prl + pbl

           vvl  = vxl*vxl+vyl*vyl+vzl*vzl
           lgl  = dsqrt(1d0+vvl)
           hhl  = 1d0 + gm*prl/(rol*(gm-1d0))
           drol = rol*lgl
           ilg = 1d0/lgl
           vb = ( vxl*bxs + vyl*byl + vzl*bzl )*ilg
           rxl = (drol*lgl*hhl + 2d0*pbl)*vxl*ilg - vb*bxs
           ryl = (drol*lgl*hhl + 2d0*pbl)*vyl*ilg - vb*byl
           rzl = (drol*lgl*hhl + 2d0*pbl)*vzl*ilg - vb*bzl
           eel =drol*hhl*lgl - prl + pbl + 0.5d0*(vvl*ilg*ilg*2d0*pbl -vb*vb )


!---- Right state

           ror = row(i,j,k,2)
           prr = prw(i,j,k,2)
           vxr = vxw(i,j,k,2)
           vyr = vyw(i,j,k,2)
           vzr = vzw(i,j,k,2)
           byr = byw(i,j,k,2)
           bzr = bzw(i,j,k,2)

           bsqr = bxs**2+byr**2+bzr**2
           pbr = 0.5d0*(bxs**2 + byr**2 + bzr**2)
           ptr = prr + pbr

           vvr  = vxr*vxr+vyr*vyr+vzr*vzr
           lgr  = dsqrt(1d0+vvr)
           hhr  = 1d0 + gm*prr/(ror*(gm-1d0))
           dror = ror*lgr
           ilg = 1d0/lgr
           vb = ( vxr*bxs + vyr*byr + vzr*bzr )*ilg
           rxr = (dror*lgr*hhr + 2d0*pbr)*vxr*ilg - vb*bxs
           ryr = (dror*lgr*hhr + 2d0*pbr)*vyr*ilg - vb*byr
           rzr = (dror*lgr*hhr + 2d0*pbr)*vzr*ilg - vb*bzr
           eer =dror*hhr*lgr - prr + pbr + 0.5d0*(vvr*ilg*ilg*2d0*pbr -vb*vb )

!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!

           ! ** Wave estimate Leismann et al. (2005) ** !
           ! Left hand Side
           roh = rol + prl*gm/(gm-1d0)
           cs2 = gm*prl/roh

           bb  = (bxs*bxs+byl*byl+bzl*bzl)/(lgl*lgl) &
              + ((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2
           sig = bb/rol
           hs  = roh/rol + sig

           ca2  = sig/hs
           om2 = cs2 + ca2 - cs2*ca2

           rr  = cs2*(((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2)/(rol*hs*lgl*lgl)

           vsq = (vxl*vxl+vyl*vyl+vzl*vzl)/(lgl*lgl)
           l1 = vxl*(1d0-om2)/lgl
           l2 = (1d0-om2*vsq-rr)
           l3 = ((vsq-1d0)*om2+rr)*( (vsq-vxl*vxl/(lgl*lgl))*om2 + (vxl*vxl)/(lgl*lgl) -1d0 + rr )
           lfl = l1/l2 + dsqrt(l3)/l2
           lsl = l1/l2 - dsqrt(l3)/l2

           ! right hand side
            roh = ror + prr*gm/(gm-1d0)
            cs2 = gm*prr/roh

            bb  = (bxs*bxs+byr*byr+bzr*bzr)/(lgr*lgr) &
                + ((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2
            sig = bb/ror
            hs  = roh/ror + sig

            ca2  = sig/hs
            om2 = cs2 + ca2 - cs2*ca2

            rr  = cs2*(((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2)/(ror*hs*lgr*lgr)

            vsq = (vxr*vxr+vyr*vyr+vzr*vzr)/(lgr*lgr)
            l1 = vxr*(1d0-om2)/lgr
            l2 = (1d0-om2*vsq-rr)
            l3 = ((vsq-1d0)*om2+rr)*( (vsq-vxr*vxr/(lgr*lgr))*om2 + (vxr*vxr)/(lgr*lgr) -1d0 + rr )
            lfr = l1/l2 + dsqrt(l3)/l2
            lsr = l1/l2 - dsqrt(l3)/l2

            sl = min(lsl,lsr)
            sr = max(lfl,lfr)

!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxes
!
! Left value
           ! Calc mangetic field four vector
           bbl  = bsql/(lgl*lgl) + ((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2
           b0xl = lgl*( bxs/(lgl*lgl) + vxl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))
           b0yl = lgl*( byl/(lgl*lgl) + vyl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))
           b0zl = lgl*( bzl/(lgl*lgl) + vzl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))

           frol = drol*vxl/lgl
           frxl = rxl*vxl/lgl + prl + 0.5d0*bbl - bxs*b0xl/(lgl)
           fryl = ryl*vxl/lgl - bxs*b0yl/(lgl)
           frzl = rzl*vxl/lgl - bxs*b0zl/(lgl)
           feel = rxl
           fbyl = byl*vxl/lgl - bxs*vyl/lgl
           fbzl = bzl*vxl/lgl - bxs*vzl/lgl

! Right value
           bbr  = bsqr/(lgr*lgr) + ((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2
           b0xr = lgr*( bxs/(lgr*lgr) + vxr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))
           b0yr = lgr*( byr/(lgr*lgr) + vyr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))
           b0zr = lgr*( bzr/(lgr*lgr) + vzr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))

           fror = dror*vxr/lgr
           frxr = rxr*vxr/lgr + prr +0.5d0*bbr - bxs*b0xr/(lgr)
           fryr = ryr*vxr/lgr - bxs*b0yr/(lgr)
           frzr = rzr*vxr/lgr - bxs*b0zr/(lgr)
           feer = rxr
           fbyr = byr*vxr/lgr - bxs*vyr/lgr
           fbzr = bzr*vxr/lgr - bxs*vzr/lgr

!----- Step 3. ----------------------------------------------------------|
! return upwind flux
!

           if (sl >= 0.0d0) then
              fro(i,j,k) = frol
              frx(i,j,k) = frxl
              fry(i,j,k) = fryl
              frz(i,j,k) = frzl
              fee(i,j,k) = feel
              fby(i,j,k) = fbyl
              fbz(i,j,k) = fbzl
              cycle
           endif

           if (sr <= 0.0d0) then
              fro(i,j,k) = fror
              frx(i,j,k) = frxr
              fry(i,j,k) = fryr
              frz(i,j,k) = frzr
              fee(i,j,k) = feer
              fby(i,j,k) = fbyr
              fbz(i,j,k) = fbzr
              cycle
           endif
!----- Step 4. ----------------------------------------------------------|
! compute HLL flux
!

           fro(i,j,k) = (sr*frol-sl*fror+sr*sl*(dror-drol))/(sr-sl)
           fee(i,j,k) = (sr*feel-sl*feer+sr*sl*(eer-eel))/(sr-sl)

           frx(i,j,k) = (sr*frxl-sl*frxr+sr*sl*(rxr-rxl))/(sr-sl)
           fry(i,j,k) = (sr*fryl-sl*fryr+sr*sl*(ryr-ryl))/(sr-sl)
           frz(i,j,k) = (sr*frzl-sl*frzr+sr*sl*(rzr-rzl))/(sr-sl)

           fby(i,j,k) = (sr*fbyl-sl*fbyr+sr*sl*(byr-byl))/(sr-sl)
           fbz(i,j,k) = (sr*fbzl-sl*fbzr+sr*sl*(bzr-bzl))/(sr-sl)
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  end subroutine flux_calc__hll

  !*****************************************************************!
  !      Calculation of Numerical flux by HLLC Method               !
  !             Honkkila and Janhunen 2007                          !
  !*****************************************************************!
  subroutine flux_calc__hllc(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,margin,ix,jx,kx &
                           ,fro,fee,frx,fry,frz,fby,fbz)

  integer,intent(in) :: ix,jx,kx,margin
  real(8),intent(in) :: gm
! primitive variables :: 1 = left state , 2 = right state
  real(8),dimension(ix,jx,kx,2),intent(in) :: row,prw,vxw,vyw,vzw
  real(8),dimension(ix,jx,kx,2),intent(in) :: byw,bzw
  real(8),dimension(ix,jx,kx),intent(in) :: bx ! magnetic field
  real(8),dimension(ix,jx,kx),intent(out) :: fro,fee,frx,fry,frz ! numerical flux
  real(8),dimension(ix,jx,kx),intent(out) :: fby,fbz

!----- U -----
! qql :: left state
! qqr :: right state
  real(8) :: rol,vxl,vyl,vzl,byl,bzl,ptl,eel
  real(8) :: ror,vxr,vyr,vzr,byr,bzr,ptr,eer
  real(8) :: rxl,ryl,rzl
  real(8) :: rxr,ryr,rzr
  real(8) :: bxs,bxsq
  real(8) :: pbl,pbr,prl,prr
  real(8) :: gmpl,gmpr,gpbl,gpbr
  real(8) :: vvl,lgl,hhl,drol
  real(8) :: vvr,lgr,hhr,dror
  real(8) :: vb

!----- convariant mag field
  real(8) :: bbl,b0xl,b0yl,b0zl,bbr,b0xr,b0yr,b0zr

!----- Wave Estimate
  real(8) :: lg,ilg,roh,cs2,bb,sig,hs,ca2,om2,rr,vsq,l1,l2,l3,lfl,lfr,lsl,lsr
!----- U* ----
! qqlst :: left state
! qqrst :: right state
  real(8) :: sl,sr
!----- flux ---
! fqql :: left physical flux
! fqqr :: right physical flux
! fluxqq :: intermediate HLLD flux (OUTPUT)
  real(8) :: frol,frxl,fryl,frzl
  real(8) :: fbyl,fbzl,feel
  real(8) :: fror,frxr,fryr,frzr
  real(8) :: fbyr,fbzr,feer
  real(8) :: bsql,bsqr
  real(8) :: cfl,cfr
  integer :: i,j,k

!---- HLLC ----
  real(8) :: u1,u2,u3,u4,u5,u6,u7,u8
  real(8) :: ros,prs,vxs,vys,vzs,lgs,bbs
  real(8) :: ees,bys,bzs,c1,g1,cy,cz,dros,rxs,rys,rzs,b0xs,b0ys,b0zs,hhs,bsq

  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j) &
  !$OMP PRIVATE(bxs,bxsq) &
  ! --- l state ---
  !$OMP PRIVATE(rol,vxl,vyl,vzl,byl,bzl,bsql,pbl,prl,ptl,eel,rxl,ryl,rzl) &
  !$OMP PRIVATE(vvl,lgl,hhl,drol,ilg,vb) &    
  ! --- r state ---
  !$OMP PRIVATE(ror,vxr,vyr,vzr,byr,bzr,bsqr,pbr,prr,ptr,eer,rxr,ryr,rzr) &
  !$OMP PRIVATE(vvr,lgr,hhr,dror) &      
  !--- step 1----
  !$OMP PRIVATE(gmpl,gmpr,gpbl,gpbr,cfl,cfr,sl,sr) &
  !$OMP PRIVATE(roh,cs2,bb,sig,hs,ca2,om2,rr,vsq,l1,l2,l3,lfl,lsl,lfr,lsr) &        
  !--- step 2 ---
  ! - left
  !$OMP PRIVATE(frol,frxl,fryl,frzl,feel,fbyl,fbzl) &
  !$OMP PRIVATE(bbl,b0xl,b0yl,b0zl) &    
  ! - right
  !$OMP PRIVATE(fror,frxr,fryr,frzr,feer,fbyr,fbzr)&
  !$OMP PRIVATE(bbr,b0xr,b0yr,b0zr) &
  !--- step 4---  
  !$OMP PRIVATE(u1,u2,u3,u4,u5,u6,u7,u8)&
  !$OMP PRIVATE(ros,prs,vxs,vys,vzs,lgs,bbs,ees,bys,bzs)&
  !$OMP PRIVATE(c1,g1,cy,cz)
  do k=margin,kx-margin
     do j=margin,jx-margin
        do i=margin,ix-margin
!----- Step 0. ----------------------------------------------------------|
! set L/R-state
!
           bxs = bx(i,j,k)
           bxsq = bxs**2
!---- Left state

           rol = row(i,j,k,1)
           prl = prw(i,j,k,1)
           vxl = vxw(i,j,k,1)
           vyl = vyw(i,j,k,1)
           vzl = vzw(i,j,k,1)
           byl = byw(i,j,k,1)
           bzl = bzw(i,j,k,1)

           bsql = bxs**2+byl**2+bzl**2
           pbl = 0.5d0*(bxs**2 + byl**2 + bzl**2)
           ptl = prl + pbl

           vvl  = vxl*vxl+vyl*vyl+vzl*vzl
           lgl  = dsqrt(1d0+vvl)
           hhl  = 1d0 + gm*prl/(rol*(gm-1d0))
           drol = rol*lgl
           ilg = 1d0/lgl
           vb = ( vxl*bxs + vyl*byl + vzl*bzl )*ilg
           rxl = (drol*lgl*hhl + 2d0*pbl)*vxl*ilg - vb*bxs
           ryl = (drol*lgl*hhl + 2d0*pbl)*vyl*ilg - vb*byl
           rzl = (drol*lgl*hhl + 2d0*pbl)*vzl*ilg - vb*bzl
           eel =drol*hhl*lgl - prl + pbl + 0.5d0*(vsq*ilg*ilg*2d0*pbl -vb*vb )


!---- Right state

           ror = row(i,j,k,2)
           prr = prw(i,j,k,2)
           vxr = vxw(i,j,k,2)
           vyr = vyw(i,j,k,2)
           vzr = vzw(i,j,k,2)
           byr = byw(i,j,k,2)
           bzr = bzw(i,j,k,2)

           bsqr = bxs**2+byr**2+bzr**2
           pbr = 0.5d0*(bxs**2 + byr**2 + bzr**2)
           ptr = prr + pbr

           vvr  = vxr*vxr+vyr*vyr+vzr*vzr
           lgr  = dsqrt(1d0+vvr)
           hhr  = 1d0 + gm*prr/(ror*(gm-1d0))
           dror = ror*lgr
           ilg = 1d0/lgr
           vb = ( vxr*bxs + vyr*byr + vzr*bzr )*ilg
           rxr = (dror*lgr*hhr + 2d0*pbr)*vxr*ilg - vb*bxs
           ryr = (dror*lgr*hhr + 2d0*pbr)*vyr*ilg - vb*byr
           rzr = (dror*lgr*hhr + 2d0*pbr)*vzr*ilg - vb*bzr
           eer =dror*hhr*lgr - prr + pbr + 0.5d0*(vsq*ilg*ilg*2d0*pbr -vb*vb )

!----- Step 1. ----------------------------------------------------------|
! Compute wave left & right wave speed
!

           ! ** Wave estimate Leismann et al. (2005) ** !
           ! Left hand Side
           roh = rol + prl*gm/(gm-1d0)
           cs2 = gm*prl/roh

           bb  = (bxs*bxs+byl*byl+bzl*bzl)/(lgl*lgl) &
              + ((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2
           sig = bb/rol
           hs  = roh/rol + sig

           ca2  = sig/hs
           om2 = cs2 + ca2 - cs2*ca2

           rr  = cs2*(((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2)/(rol*hs*lgl*lgl)

           vsq = (vxl*vxl+vyl*vyl+vzl*vzl)/(lgl*lgl)
           l1 = vxl*(1d0-om2)/lgl
           l2 = (1d0-om2*vsq-rr)
           l3 = ((vsq-1d0)*om2+rr)*( (vsq-vxl*vxl/(lgl*lgl))*om2 + (vxl*vxl)/(lgl*lgl) -1d0 + rr )
           lfl = l1/l2 + dsqrt(l3)/l2
           lsl = l1/l2 - dsqrt(l3)/l2

           ! right hand side
            roh = ror + prr*gm/(gm-1d0)
            cs2 = gm*prr/roh

            bb  = (bxs*bxs+byr*byr+bzr*bzr)/(lgr*lgr) &
                + ((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2
            sig = bb/ror
            hs  = roh/ror + sig

            ca2  = sig/hs
            om2 = cs2 + ca2 - cs2*ca2

            rr  = cs2*(((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2)/(ror*hs*lgr*lgr)

            vsq = (vxr*vxr+vyr*vyr+vzr*vzr)/(lgr*lgr)
            l1 = vxr*(1d0-om2)/lgr
            l2 = (1d0-om2*vsq-rr)
            l3 = ((vsq-1d0)*om2+rr)*( (vsq-vxr*vxr/(lgr*lgr))*om2 + (vxr*vxr)/(lgr*lgr) -1d0 + rr )
            lfr = l1/l2 + dsqrt(l3)/l2
            lsr = l1/l2 - dsqrt(l3)/l2

            sl = min(lsl,lsr)
            sr = max(lfl,lfr)

!----- Step 2. ----------------------------------------------------------|
! compute L/R fluxes
!
! Left value
           ! Calc mangetic field four vector
           bbl  = bsql/(lgl*lgl) + ((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2
           b0xl = lgl*( bxs/(lgl*lgl) + vxl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))
           b0yl = lgl*( byl/(lgl*lgl) + vyl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))
           b0zl = lgl*( bzl/(lgl*lgl) + vzl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))

           frol = drol*vxl/lgl
           frxl = rxl*vxl/lgl + prl + 0.5d0*bbl - bxs*b0xl/(lgl)
           fryl = ryl*vxl/lgl - bxs*b0yl/(lgl)
           frzl = rzl*vxl/lgl - bxs*b0zl/(lgl)
           feel = rxl
           fbyl = byl*vxl/lgl - bxs*vyl/lgl
           fbzl = bzl*vxl/lgl - bxs*vzl/lgl

! Right value
           bbr  = bsqr/(lgr*lgr) + ((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2
           b0xr = lgr*( bxs/(lgr*lgr) + vxr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))
           b0yr = lgr*( byr/(lgr*lgr) + vyr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))
           b0zr = lgr*( bzr/(lgr*lgr) + vzr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))

           fror = dror*vxr/lgr
           frxr = rxr*vxr/lgr + prr +0.5d0*bbr - bxs*b0xr/(lgr)
           fryr = ryr*vxr/lgr - bxs*b0yr/(lgr)
           frzr = rzr*vxr/lgr - bxs*b0zr/(lgr)
           feer = rxr
           fbyr = byr*vxr/lgr - bxs*vyr/lgr
           fbzr = bzr*vxr/lgr - bxs*vzr/lgr


!----- Step 3. ----------------------------------------------------------|
! return upwind flux
!

           if (sl >= 0.0d0) then
              fro(i,j,k) = frol
              frx(i,j,k) = frxl
              fry(i,j,k) = fryl
              frz(i,j,k) = frzl
              fee(i,j,k) = feel
              fby(i,j,k) = fbyl
              fbz(i,j,k) = fbzl
              cycle
           endif

           if (sr <= 0.0d0) then
              fro(i,j,k) = fror
              frx(i,j,k) = frxr
              fry(i,j,k) = fryr
              frz(i,j,k) = frzr
              fee(i,j,k) = feer
              fby(i,j,k) = fbyr
              fbz(i,j,k) = fbzr
              cycle
           endif

!----- Step.4 -----------------------------------------------------------!
! Calc. U_hll
!
          u1 = (sr*dror-sl*drol-fror+frol)/(sr-sl)
          u2 = (sr*rxr -sl*rxl -frxr+frxl)/(sr-sl)
          u3 = (sr*ryr -sl*ryl -fryr+fryl)/(sr-sl)
          u4 = (sr*rzr -sl*rzl -frzr+frzl)/(sr-sl)
          u5 = bxs
          u6 = (sr*byr -sl*byl -fbyr+fbyl)/(sr-sl)
          u7 = (sr*bzr -sl*bzl -fbzr+fbzl)/(sr-sl)
          u8 = (sr*eer -sl*eel -feer+feel)/(sr-sl)

          !! for initial guess
          ros = max(ror,rol)
          prs = max(prr,prl)
          vxs = max(vxr,vxl)
          vzs = max(vyr,vyl)
          vys = max(vzr,vzl)
          call convert_ctop(gm,ros,u1,u8,u2,u3,u4,u5,u6,u7,vxs,vys,vzs,prs)
          lgs = dsqrt(1d0+vxs*vxs+vys*vys+vzs*vzs)
          bbs = (u5*u5+u6*u6+u7*u7)/(lgs*lgs) + ((vxs*u5+vys*u6+vzs*u7)/lgs)**2
          bxs = u5
          vxs = vxs/lgs
          vys = vys/lgs
          vzs = vzs/lgs

!----- Step 4. ----------------------------------------------------------|
! compute HLL flux
!
           if( abs(bxs) .lt. 1d-8 )then
           if (vxs >= 0.0d0) then
             !! Left intermediate state
             ees = (sl*eel-rxl+(prs+0.5d0*bbs)*vxs)/(sl-vxs)
             bys = byl*(sl-vxl/lgl)/(sl-vxs)
             bzs = bzl*(sl-vxl/lgl)/(sl-vxs)

             c1  = (vxs-sl)*bys*bzs
             g1  = (sl*eel-feel+sl*(prs+0.5d0*bbs))/(sl-vxs) - (bxs*bxs+bys*bys+bzs*bzs)
             cy  = (sl-vxs)*(g1+bys*bys)
             cz  = (sl-vxs)*(g1+bzs*bzs)
             vys = (cy*(sl*ryl - fryl) - c1*(sl*rzl-frzl) )/(cy*cz-c1*c1)
             vzs = (cz*(sl*rzl - frzl) - c1*(sl*ryl-fryl) )/(cy*cz-c1*c1)

!             lgs = 1d0/dsqrt(1d0-(vxs*vxs+vys*vys+vzs*vzs))
!             ros = lgl*rol*(sl-vxl/lgl)/(lgs*(sl-vxs))
!             dros = ros*lgs

             dros = drol*(sl-vxl/lgl)/(sl-vxs)
             rxs = (ees+prs+0.5d0*bbs)*vxs - bxs*(vxs*bxs+vys*bys+vzs*bzs)
             rys = (ees+prs+0.5d0*bbs)*vys - bys*(vxs*bxs+vys*bys+vzs*bzs)
             rzs = (ees+prs+0.5d0*bbs)*vzs - bzs*(vxs*bxs+vys*bys+vzs*bzs)

             ! b0xs = lgs*( bxs/(lgs*lgs) + vxs*(vxs*bxs+bys*vys+bzs*vzs))
             ! b0ys = lgs*( bys/(lgs*lgs) + vys*(vxs*bxs+bys*vys+bzs*vzs))
             ! b0zs = lgs*( bzs/(lgs*lgs) + vzs*(vxs*bxs+bys*vys+bzs*vzs))

             fro(i,j,k) = frol + sl*(dros - drol )
             frx(i,j,k) = frxl + sl*( rxs - rxl  )
             fry(i,j,k) = fryl + sl*( rys - ryl  )
             frz(i,j,k) = frzl + sl*( rzs - rzl  )
             fee(i,j,k) = feel + sl*( ees - eel  )
             fby(i,j,k) = fbyl + sl*( bys - byl  )
             fbz(i,j,k) = fbzl + sl*( bzs - bzl  )

             ! fro(i,j,k) = dros*vxs
             ! frx(i,j,k) = rxs*vxs + prs +0.5d0*bbs - bxs*b0xs/(lgs)
             ! fry(i,j,k) = rys*vxs - bxs*b0ys/(lgs)
             ! frz(i,j,k) = rzs*vxs - bxs*b0zs/(lgs)
             ! fee(i,j,k) = rxs
             ! fby(i,j,k) = bys*vxs - bxs*vys
             ! fbz(i,j,k) = bzs*vxs - bxs*vzs
           else
             !! Right intermediate state
             ees = (sr*eer-rxr+(prs+0.5d0*bbs)*vxs)/(sr-vxs)
             bys = byr*(sr-vxr/lgr)/(sr-vxs)
             bzs = bzr*(sr-vxr/lgr)/(sr-vxs)

             c1  = (vxs-sr)*bys*bzs
             g1  = (sr*eer-feer+sr*(prs+0.5d0*bbs))/(sr-vxs) - (bxs*bxs+bys*bys+bzs*bzs)
             cy  = (sr-vxs)*(g1+bys*bys)
             cz  = (sr-vxs)*(g1+bzs*bzs)
             vys = (cy*(sr*ryr - fryr) - c1*(sr*rzr-frzr) )/(cy*cz-c1*c1)
             vzs = (cz*(sr*rzr - frzr) - c1*(sr*ryr-fryr) )/(cy*cz-c1*c1)

             dros = dror*(sr-vxr/lgr)/(sr-vxs)
             rxs = (ees+prs+0.5d0*bbs)*vxs - bxs*(vxs*bxs+vys*bys+vzs*bzs)
             rys = (ees+prs+0.5d0*bbs)*vys - bys*(vxs*bxs+vys*bys+vzs*bzs)
             rzs = (ees+prs+0.5d0*bbs)*vzs - bzs*(vxs*bxs+vys*bys+vzs*bzs)


             fro(i,j,k) = fror + sr*(dros - dror )
             frx(i,j,k) = frxr + sr*( rxs - rxr  )
             fry(i,j,k) = fryr + sr*( rys - ryr  )
             frz(i,j,k) = frzr + sr*( rzs - rzr  )
             fee(i,j,k) = feer + sr*( ees - eer  )
             fby(i,j,k) = fbyr + sr*( bys - byr  )
             fbz(i,j,k) = fbzr + sr*( bzs - bzr  )

             ! fro(i,j,k) = dros*vxs
             ! frx(i,j,k) = rxs*vxs + prs +0.5d0*bbs - bxs*b0xs/(lgs)
             ! fry(i,j,k) = rys*vxs - bxs*b0ys/(lgs)
             ! frz(i,j,k) = rzs*vxs - bxs*b0zs/(lgs)
             ! fee(i,j,k) = rxs
             ! fby(i,j,k) = bys*vxs - bxs*vys
             ! fbz(i,j,k) = bzs*vxs - bxs*vzs
           endif

         else
           ! case bx .neq. 0
           if (vxs >= 0.0d0) then
             !! Left intermediate state
             bys = u6
             bzs = u7
             ros = drol*(sl-vxl/lgl)/(lgs*(sl-vxs))

             dros = ros*lgs

             hhs = 1d0 + gm*prs/(ros*(gm-1d0))
             bsq = bxs*bxs+bys*bys+bzs*bzs
             vb  = vxs*bxs+vys*bys+vzs*bzs
             rxs = (dros*hhs*lgs+bsq)*vxs - bxs*vb
             rys = (dros*hhs*lgs+bsq)*vys - bys*vb
             rzs = (dros*hhs*lgs+bsq)*vzs - bzs*vb
             ees = dros*hhs*lgs - prs + 0.5d0*bsq+0.5*(bsq*(vxs*vxs+vys*vys+vzs*vzs) - vb*vb)

             fro(i,j,k) = frol + sl*(dros - drol )
             frx(i,j,k) = frxl + sl*( rxs - rxl  )
             fry(i,j,k) = fryl + sl*( rys - ryl  )
             frz(i,j,k) = frzl + sl*( rzs - rzl  )
             fee(i,j,k) = feel + sl*( ees - eel  )
             fby(i,j,k) = fbyl + sl*( bys - byl  )
             fbz(i,j,k) = fbzl + sl*( bzs - bzl  )
           else
             !! Right intermediate state
             bys = u6
             bzs = u7
             ros = dror*(sr-vxr/lgr)/(lgs*(sr-vxs))

             dros = ros*lgs

             hhs = 1d0 + gm*prs/(ros*(gm-1d0))
             bsq = bxs*bxs+bys*bys+bzs*bzs
             vb  = vxs*bxs+vys*bys+vzs*bzs
             rxs = (dros*hhs*lgs+bsq)*vxs - bxs*vb
             rys = (dros*hhs*lgs+bsq)*vys - bys*vb
             rzs = (dros*hhs*lgs+bsq)*vzs - bzs*vb
             ees = dros*hhs*lgs - prs + 0.5d0*bsq+0.5*(bsq*(vxs*vxs+vys*vys+vzs*vzs) - vb*vb)

             fro(i,j,k) = fror + sr*(dros - dror )
             frx(i,j,k) = frxr + sr*( rxs - rxr  )
             fry(i,j,k) = fryr + sr*( rys - ryr  )
             frz(i,j,k) = frzr + sr*( rzs - rzr  )
             fee(i,j,k) = feer + sr*( ees - eer  )
             fby(i,j,k) = fbyr + sr*( bys - byr  )
             fbz(i,j,k) = fbzr + sr*( bzs - bzr  )
           endif
         endif

        end do
     end do
  end do
  !$OMP END PARALLEL DO
end subroutine flux_calc__hllc

  !**************************************************************!
  !      Calculation of Numerical flux by HLLC Method            !
  !             J. Kim and D.S. Balsara  2014                    !
  !**************************************************************!
!   subroutine flux_calc__hllc(row,prw,vxw,vyw,vzw,bx,byw,bzw,gm,margin,ix,jx,kx &
!                            ,fro,fee,frx,fry,frz,fby,fbz)

!   integer,intent(in) :: ix,jx,kx,margin
!   real(8),intent(in) :: gm
! ! primitive variables :: 1 = left state , 2 = right state
!   real(8),dimension(ix,jx,kx,2),intent(in) :: row,prw,vxw,vyw,vzw
!   real(8),dimension(ix,jx,kx,2),intent(in) :: byw,bzw
!   real(8),dimension(ix,jx,kx),intent(in) :: bx ! magnetic field
!   real(8),dimension(ix,jx,kx),intent(out) :: fro,fee,frx,fry,frz ! numerical flux
!   real(8),dimension(ix,jx,kx),intent(out) :: fby,fbz

! !----- U -----
! ! qql :: left state
! ! qqr :: right state
!   real(8) :: rol,vxl,vyl,vzl,byl,bzl,ptl,eel
!   real(8) :: ror,vxr,vyr,vzr,byr,bzr,ptr,eer
!   real(8) :: rxl,ryl,rzl
!   real(8) :: rxr,ryr,rzr
!   real(8) :: bxs,bxsq
!   real(8) :: pbl,pbr,prl,prr
!   real(8) :: gmpl,gmpr,gpbl,gpbr
!   real(8) :: vvl,lgl,hhl,drol
!   real(8) :: vvr,lgr,hhr,dror
!   real(8) :: vb

! !----- convariant mag field
!   real(8) :: bbl,b0xl,b0yl,b0zl,bbr,b0xr,b0yr,b0zr

! !----- Wave Estimate
!   real(8) :: lg,ilg,roh,cs2,bb,sig,hs,ca2,om2,rr,vsq,l1,l2,l3,lfl,lfr,lsl,lsr
! !----- U* ----
! ! qqlst :: left state
! ! qqrst :: right state
!   real(8) :: sl,sr
! !----- flux ---
! ! fqql :: left physical flux
! ! fqqr :: right physical flux
! ! fluxqq :: intermediate HLLD flux (OUTPUT)
!   real(8) :: frol,frxl,fryl,frzl
!   real(8) :: fbyl,fbzl,feel
!   real(8) :: fror,frxr,fryr,frzr
!   real(8) :: fbyr,fbzr,feer
!   real(8) :: bsql,bsqr
!   real(8) :: cfl,cfr
!   integer :: i,j,k

! !---- HLLC ----
!   real(8) :: u1,u2,u3,u4,u5,u6,u7,u8
!   real(8) :: ros,prs,vxs,vys,vzs,lgs,bbs
!   real(8) :: ees,bys,bzs,c1,g1,cy,cz,dros,rxs,rys,rzs,b0xs,b0ys,b0zs,hhs,bsq
!   real(8) :: vbs,ff1,ff2,a11,a12,a21,a22,aaa,dvys,dvzs,ptot
!   real(8) ,parameter :: eps = 1d-10
!   integer,parameter :: cmax=100
!   integer :: ll

!   !$OMP PARALLEL DO &
!   !$OMP PRIVATE(i,j) &
!   !$OMP PRIVATE(bxs,bxsq) &
!   ! --- l state ---
!   !$OMP PRIVATE(rol,vxl,vyl,vzl,byl,bzl,bsql,pbl,prl,ptl,eel,rxl,ryl,rzl) &
!   ! --- r state ---
!   !$OMP PRIVATE(ror,vxr,vyr,vzr,byr,bzr,bsqr,pbr,prr,ptr,eer,rxr,ryr,rzr) &
!   !--- step 1----
!   !$OMP PRIVATE(gmpl,gmpr,gpbl,gpbr,cfl,cfr,sl,sr) &
!   !--- step 2 ---
!   ! - left
!   !$OMP PRIVATE(frol,frxl,fryl,frzl,feel,fbyl,fbzl) &
!   ! - right
!   !$OMP PRIVATE(fror,frxr,fryr,frzr,feer,fbyr,fbzr)
!   do k=margin,kx-margin
!      do j=margin,jx-margin
!         do i=margin,ix-margin
! !----- Step 0. ----------------------------------------------------------|
! ! set L/R-state
! !
!            bxs = bx(i,j,k)
!            bxsq = bxs**2
! !---- Left state

!            rol = row(i,j,k,1)
!            prl = prw(i,j,k,1)
!            vxl = vxw(i,j,k,1)
!            vyl = vyw(i,j,k,1)
!            vzl = vzw(i,j,k,1)
!            byl = byw(i,j,k,1)
!            bzl = bzw(i,j,k,1)

!            bsql = bxs**2+byl**2+bzl**2
!            pbl = 0.5d0*(bxs**2 + byl**2 + bzl**2)
!            ptl = prl + pbl

!            vvl  = vxl*vxl+vyl*vyl+vzl*vzl
!            lgl  = dsqrt(1d0+vvl)
!            hhl  = 1d0 + gm*prl/(rol*(gm-1d0))
!            drol = rol*lgl
!            ilg = 1d0/lgl
!            vb = ( vxl*bxs + vyl*byl + vzl*bzl )*ilg
!            rxl = (drol*lgl*hhl + 2d0*pbl)*vxl*ilg - vb*bxs
!            ryl = (drol*lgl*hhl + 2d0*pbl)*vyl*ilg - vb*byl
!            rzl = (drol*lgl*hhl + 2d0*pbl)*vzl*ilg - vb*bzl
!            eel =drol*hhl*lgl - prl + pbl + 0.5d0*(vvl*ilg*ilg*2d0*pbl -vb*vb )


! !---- Right state

!            ror = row(i,j,k,2)
!            prr = prw(i,j,k,2)
!            vxr = vxw(i,j,k,2)
!            vyr = vyw(i,j,k,2)
!            vzr = vzw(i,j,k,2)
!            byr = byw(i,j,k,2)
!            bzr = bzw(i,j,k,2)

!            bsqr = bxs**2+byr**2+bzr**2
!            pbr = 0.5d0*(bxs**2 + byr**2 + bzr**2)
!            ptr = prr + pbr

!            vvr  = vxr*vxr+vyr*vyr+vzr*vzr
!            lgr  = dsqrt(1d0+vvr)
!            hhr  = 1d0 + gm*prr/(ror*(gm-1d0))
!            dror = ror*lgr
!            ilg = 1d0/lgr
!            vb = ( vxr*bxs + vyr*byr + vzr*bzr )*ilg
!            rxr = (dror*lgr*hhr + 2d0*pbr)*vxr*ilg - vb*bxs
!            ryr = (dror*lgr*hhr + 2d0*pbr)*vyr*ilg - vb*byr
!            rzr = (dror*lgr*hhr + 2d0*pbr)*vzr*ilg - vb*bzr
!            eer =dror*hhr*lgr - prr + pbr + 0.5d0*(vvr*ilg*ilg*2d0*pbr -vb*vb )

! !----- Step 1. ----------------------------------------------------------|
! ! Compute wave left & right wave speed
! !

!            ! ** Wave estimate Leismann et al. (2005) ** !
!            ! Left hand Side
!            roh = rol + prl*gm/(gm-1d0)
!            cs2 = gm*prl/roh

!            bb  = (bxs*bxs+byl*byl+bzl*bzl)/(lgl*lgl) &
!               + ((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2
!            sig = bb/rol
!            hs  = roh/rol + sig

!            ca2  = sig/hs
!            om2 = cs2 + ca2 - cs2*ca2

!            rr  = cs2*(((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2)/(rol*hs*lgl*lgl)

!            vsq = (vxl*vxl+vyl*vyl+vzl*vzl)/(lgl*lgl)
!            l1 = vxl*(1d0-om2)/lgl
!            l2 = (1d0-om2*vsq-rr)
!            l3 = ((vsq-1d0)*om2+rr)*( (vsq-vxl*vxl/(lgl*lgl))*om2 + (vxl*vxl)/(lgl*lgl) -1d0 + rr )
!            lfl = l1/l2 + dsqrt(l3)/l2
!            lsl = l1/l2 - dsqrt(l3)/l2

!            ! right hand side
!             roh = ror + prr*gm/(gm-1d0)
!             cs2 = gm*prr/roh

!             bb  = (bxs*bxs+byr*byr+bzr*bzr)/(lgr*lgr) &
!                 + ((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2
!             sig = bb/ror
!             hs  = roh/ror + sig

!             ca2  = sig/hs
!             om2 = cs2 + ca2 - cs2*ca2

!             rr  = cs2*(((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2)/(ror*hs*lgr*lgr)

!             vsq = (vxr*vxr+vyr*vyr+vzr*vzr)/(lgr*lgr)
!             l1 = vxr*(1d0-om2)/lgr
!             l2 = (1d0-om2*vsq-rr)
!             l3 = ((vsq-1d0)*om2+rr)*( (vsq-vxr*vxr/(lgr*lgr))*om2 + (vxr*vxr)/(lgr*lgr) -1d0 + rr )
!             lfr = l1/l2 + dsqrt(l3)/l2
!             lsr = l1/l2 - dsqrt(l3)/l2

!             sl = min(lsl,lsr)
!             sr = max(lfl,lfr)

! !----- Step 2. ----------------------------------------------------------|
! ! compute L/R fluxes
! !
! ! Left value
!            ! Calc mangetic field four vector
!            bbl  = bsql/(lgl*lgl) + ((vxl*bxs+byl*vyl+bzl*vzl)/lgl)**2
!            b0xl = lgl*( bxs/(lgl*lgl) + vxl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))
!            b0yl = lgl*( byl/(lgl*lgl) + vyl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))
!            b0zl = lgl*( bzl/(lgl*lgl) + vzl/lgl*((vxl*bxs+byl*vyl+bzl*vzl)/lgl))

!            frol = drol*vxl/lgl
!            frxl = rxl*vxl/lgl + prl + 0.5d0*bbl - bxs*b0xl/(lgl)
!            fryl = ryl*vxl/lgl - bxs*b0yl/(lgl)
!            frzl = rzl*vxl/lgl - bxs*b0zl/(lgl)
!            feel = rxl
!            fbyl = byl*vxl/lgl - bxs*vyl/lgl
!            fbzl = bzl*vxl/lgl - bxs*vzl/lgl

! ! Right value
!            bbr  = bsqr/(lgr*lgr) + ((vxr*bxs+byr*vyr+bzr*vzr)/lgr)**2
!            b0xr = lgr*( bxs/(lgr*lgr) + vxr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))
!            b0yr = lgr*( byr/(lgr*lgr) + vyr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))
!            b0zr = lgr*( bzr/(lgr*lgr) + vzr/lgr*((vxr*bxs+byr*vyr+bzr*vzr)/lgr))

!            fror = dror*vxr/lgr
!            frxr = rxr*vxr/lgr + prr +0.5d0*bbr - bxs*b0xr/(lgr)
!            fryr = ryr*vxr/lgr - bxs*b0yr/(lgr)
!            frzr = rzr*vxr/lgr - bxs*b0zr/(lgr)
!            feer = rxr
!            fbyr = byr*vxr/lgr - bxs*vyr/lgr
!            fbzr = bzr*vxr/lgr - bxs*vzr/lgr


! !----- Step 3. ----------------------------------------------------------|
! ! return upwind flux
! !
!            if (sl >= 0.0d0) then
!               fro(i,j,k) = frol
!               frx(i,j,k) = frxl
!               fry(i,j,k) = fryl
!               frz(i,j,k) = frzl
!               fee(i,j,k) = feel
!               fby(i,j,k) = fbyl
!               fbz(i,j,k) = fbzl
!               cycle
!            endif

!            if (sr <= 0.0d0) then
!               fro(i,j,k) = fror
!               frx(i,j,k) = frxr
!               fry(i,j,k) = fryr
!               frz(i,j,k) = frzr
!               fee(i,j,k) = feer
!               fby(i,j,k) = fbyr
!               fbz(i,j,k) = fbzr
!               cycle
!            endif

! !----- Step.4 -----------------------------------------------------------!
! ! Calc. U_hll

!           u1 = (sr*dror-sl*drol-fror+frol)/(sr-sl)
!           u8 = (sr*eer -sl*eel -feer+feel)/(sr-sl)
!           u2 = (sr*rxr -sl*rxl -frxr+frxl)/(sr-sl)
!           u3 = (sr*ryr -sl*ryl -fryr+fryl)/(sr-sl)
!           u4 = (sr*rzr -sl*rzl -frzr+frzl)/(sr-sl)
!           u5 = bxs
!           u6 = (sr*byr -sl*byl -fbyr+fbyl)/(sr-sl)
!           u7 = (sr*bzr -sl*bzl -fbzr+fbzl)/(sr-sl)

!           ros = rol!0.1d0!min(ror,rol)
!           prs = prl!0.1d0!min(prr,prl)
!           vxs = vxl!0.1d0!min(abs(vxr),abs(vxl))!max(vxr,vxl)
!           vys = vyl!0.1d0!min(abs(vyr),abs(vyl))!max(vyr,vyl)
!           vzs = vzl!0.1d0!min(abs(vzr),abs(vzl))!max(vzr,vzl)
!           call convert_ctop(gm,ros,u1,u8,u2,u3,u4,u5,u6,u7,vxs,vys,vzs,prs)
!           lgs = dsqrt(1d0+vxs*vxs+vys*vys+vzs*vzs)
!           bbs = (u5*u5+u6*u6+u7*u7)/(lgs*lgs) + ((vxs*u5+vys*u6+vzs*u7)/lgs)**2
!           ptot= prs+0.5d0*bbs
!           bxs = u5
!           bys = u6
!           bzs = u7

! !_______________________ Calc. Three velocity ___________________________!
!           vxs = vxs/lgs
!           vys = vys/lgs
!           vzs = vzs/lgs
! !----- Step 4. ----------------------------------------------------------!
! ! compute HLLC flux                                                      !
! !------------------------------------------------------------------------!
!           if (vxs >= 0.0d0) then
!              !! Left intermediate state
!              do ll = 1, cmax
!                vbs = vxs*bxs + vys*bys + vzs*bzs

!                ff1 = (sl*eel-feel)*vys+ptot*vys*sl-bys*vbs*(sl-vxs) &
!                    + bxs*bys*(1d0-(vxs*vxs+vys*vys+vzs*vzs))-bxs*b0yl/lgl - ryl*(sl-vxl/lgl)
!                ff2 = (sl*eel-feel)*vzs+ptot*vzs*sl-bzs*vbs*(sl-vxs) &
!                    + bxs*bzs*(1d0-(vxs*vxs+vys*vys+vzs*vzs))-bxs*b0zl/lgl - rzl*(sl-vxl/lgl)

!                a11 = feel - sl*eel-ptot*sl+(bys*bys)*(sl-vxs)+2d0*bxs*bys*vys
!                a12 = bys*bzs*(sl-vxs) + 2d0*bxs*bys*vzs
!                a21 = bys*bzs*(sl-vxs) + 2d0*bxs*bzs*vys
!                a22 = feel - sl*eel-ptot*sl+(bzs*bzs)*(sl-vxs)+2d0*bxs*bzs*vzs
!                aaa = 1d0/(a11*a22-a12*a21)

!                dvys = aaa*( a22*ff1 - a12*ff2)
!                dvzs = aaa*(-a21*ff1 + a11*ff2)

!               if( (dabs(ff1) .lt. eps) .and. (dabs(ff2) .lt. eps) )then
!                   exit
!                endif
!                vys = vys+dvys
!                vzs = vzs+dvzs
!              enddo

!              if( ll .eq. cmax+1)then
!                write(*,*)'Stop Due to faile for tranvese velocity'
!                write(*,*)'Left hand',i,vys,vzs,dvys,dvzs
!                stop
!              endif

!              dros = drol*(sl-vxl/lgl)/(sl-vxs)

!              lgs = 1d0/dsqrt(1d0-(vxs*vxs+vys*vys+vzs*vzs))
!              bbs = (bxs*bxs+bys*bys+bzs*bzs)/(lgs*lgs) + (vxs*bxs+vys*bys+vzs*bzs)**2
!              prs = ptot - 0.5d0*bbs
             
!              hhs = 1d0 + gm*prs/(dros/lgs*(gm-1d0))
!              bsq = bxs*bxs+bys*bys+bzs*bzs
!              vb  = vxs*bxs+vys*bys+vzs*bzs
!              rxs = (dros*hhs*lgs+bsq)*vxs - bxs*vb
!              rys = (dros*hhs*lgs+bsq)*vys - bys*vb
!              rzs = (dros*hhs*lgs+bsq)*vzs - bzs*vb
!              ees = dros*hhs*lgs - prs + 0.5d0*bsq+0.5*(bsq*(vxs*vxs+vys*vys+vzs*vzs) - vb*vb)

!              fro(i,j,k) = frol + sl*(dros - drol)
!              frx(i,j,k) = frxl + sl*( rxs - rxl )
!              fry(i,j,k) = fryl + sl*( rys - ryl )
!              frz(i,j,k) = frzl + sl*( rzs - rzl )
!              fee(i,j,k) = feel + sl*( ees - eel )
!              fby(i,j,k) = fbyl + sl*( bys - byl )
!              fbz(i,j,k) = fbzl + sl*( bzs - bzl )
!            else
!              !! Right intermediate state
!              do ll = 1, cmax
!                vbs = vxs*bxs + vys*bys + vzs*bzs

!                ff1 = (sr*eer-feer)*vys+ptot*vys*sr-bys*vbs*(sr-vxs) &
!                    + bxs*bys*(1d0-(vxs*vxs+vys*vys+vzs*vzs))-bxs*b0yr/lgr - ryr*(sr-vxr/lgr)
!                ff2 = (sr*eer-feer)*vzs+ptot*vzs*sr-bzs*vbs*(sr-vxs) &
!                    + bxs*bzs*(1d0-(vxs*vxs+vys*vys+vzs*vzs))-bxs*b0zr/lgr - rzr*(sr-vxr/lgr)
!                ! if( (dabs(ff1) .lt. eps) .and. (dabs(ff2) .lt. eps) )then
!                !    exit
!                ! endif
!                a11 = rxr - sr*eer-ptot*sr+(bys*bys)*(sr-vxs)+2d0*bxs*bys*vys
!                a12 = bys*bzs*(sr-vxs) + 2d0*bxs*bys*vzs
!                a21 = bys*bzs*(sr-vxs) + 2d0*bxs*bzs*vys
!                a22 = rxr - sr*eer-ptot*sr+(bzs*bzs)*(sr-vxs)+2d0*bxs*bzs*vzs
!                aaa = 1d0/(a11*a22-a12*a21)

!                dvys = aaa*(a22*ff1  - a12*ff2)
!                dvzs = aaa*(-a21*ff1 + a11*ff2)

!                ! if( (dabs((dvys)/(vys+1d-15) ) .lt. eps) .and. (dabs((dvzs)/(vzs+1d-15)) .lt. eps) )then
!                if( (dabs(ff1) .lt. eps) .and. (dabs(ff2) .lt. eps) )then
!                  exit
!                endif
!                vys = vys+dvys
!                vzs = vzs+dvzs
!              enddo

!              if( ll .eq. cmax+1)then
!                write(*,*)'Stop Due to faile for tranvese velocity'
!                write(*,*)'Right hand',i,vys,vzs
!                stop
!              endif

!              dros = dror*(sr-vxr/lgr)/(sr-vxs)

!              lgs = 1d0/dsqrt(1d0-(vxs*vxs+vys*vys+vzs*vzs))
!              bbs = (bxs*bxs+bys*bys+bzs*bzs)/(lgs*lgs) + (vxs*bxs+vys*bys+vzs*bzs)**2
!              prs = ptot - 0.5d0*bbs

!              hhs = 1d0 + gm*prs/(dros/lgs*(gm-1d0))
!              bsq = bxs*bxs+bys*bys+bzs*bzs
!              vb  = vxs*bxs+vys*bys+vzs*bzs
!              rxs = (dros*hhs*lgs+bsq)*vxs - bxs*vb
!              rys = (dros*hhs*lgs+bsq)*vys - bys*vb
!              rzs = (dros*hhs*lgs+bsq)*vzs - bzs*vb
!              ees = dros*hhs*lgs - prs + 0.5d0*bsq+0.5*(bsq*(vxs*vxs+vys*vys+vzs*vzs) - vb*vb)

!              fro(i,j,k) = fror + sr*(dros - dror)
!              frx(i,j,k) = frxr + sr*( rxs - rxr )
!              fry(i,j,k) = fryr + sr*( rys - ryr )
!              frz(i,j,k) = frzr + sr*( rzs - rzr )
!              fee(i,j,k) = feer + sr*( ees - eer )
!              fby(i,j,k) = fbyr + sr*( bys - byr )
!              fbz(i,j,k) = fbzr + sr*( bzs - bzr )
!            endif
!         end do
!      end do
!   end do
!   !$OMP END PARALLEL DO

! end subroutine flux_calc__hllc


  subroutine flux_calc__glm(bx_m,phi_m,ch,fbx,fphi,ix,jx,kx)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch ! divergence B wave speed
! cell surface Magnetic field and divergence B
  real(8),dimension(ix,jx,kx), intent(in) :: bx_m,phi_m
! magnetic field flux @ cell surface of normal component
! divergence B flux
  real(8),dimension(ix,jx,kx), intent(out) :: fbx,fphi

  integer :: i,j,k


  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           fbx(i,j,k) = phi_m(i,j,k)
           fphi(i,j,k) = bx_m(i,j,k)*ch**2
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  end subroutine flux_calc__glm


  subroutine flux_calc__bp(ix,jx,kx,bxw,phiw &
                          ,bx_m,phi_m,ch)

  integer,intent(in) :: ix,jx,kx
  real(8),intent(in) :: ch
  real(8),dimension(ix,jx,kx,2),intent(in) :: bxw,phiw
  real(8),dimension(ix,jx,kx),intent(out) :: bx_m,phi_m

  integer :: i,j,k

  ! calcurate magnetic field & divergence B @ cell surface
  !$OMP PARALLEL DO &
  !$OMP PRIVATE(i,j)
  do k=1,kx
     do j=1,jx
        do i=1,ix
           bx_m(i,j,k) = bxw(i,j,k,1) &
                +(0.5d0*(bxw(i,j,k,2)-bxw(i,j,k,1)) &
                -0.5d0*(phiw(i,j,k,2)-phiw(i,j,k,1))/ch)

           phi_m(i,j,k) = phiw(i,j,k,1) &
                +(0.5d0*(phiw(i,j,k,2)-phiw(i,j,k,1)) &
                -0.5d0*ch*(bxw(i,j,k,2)-bxw(i,j,k,1)))
        end do
     end do
  end do
  !$OMP END PARALLEL DO

  end subroutine flux_calc__bp


  subroutine flux_calc__fbres(mdir,margin,ix,jx,kx,fbx,curx,eta,pm &
                             ,fbx_res)

  integer,intent(in) :: ix,jx,kx,mdir,margin
! +1 or -1, consistency between Electric fields and numerical flux
  real(8),intent(in) :: pm
  real(8),dimension(ix,jx,kx),intent(in) :: fbx,curx,eta
  real(8),dimension(ix,jx,kx),intent(out) :: fbx_res

  integer :: i,j,k
  real(8) :: fres

! Average eta*current at cell curface
! Flux at i+1/2
  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i+1,j,k)*curx(i+1,j,k))

              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at j+1/2
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j+1,k)*curx(i,j+1,k))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at k+1/2
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres = 0.5d0*pm*(eta(i,j,k)*curx(i,j,k) &
                   +eta(i,j,k+1)*curx(i,j,k+1))
              fbx_res(i,j,k) = fbx(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  end subroutine flux_calc__fbres


  subroutine flux_calc__feres(mdir,margin,ix,jx,kx,fee,curx,cury,curz,bx,by,bz,eta &
                            ,fee_res)

  integer,intent(in) :: ix,jx,kx,mdir,margin
  real(8),dimension(ix,jx,kx),intent(in) :: fee,curx,cury,curz,bx,by,bz,eta
  real(8),dimension(ix,jx,kx),intent(out) :: fee_res

  integer :: i,j,k
  real(8) :: fres

! Average eta*current at cell curface
! Flux at i+1/2
  if(mdir == 1)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres =0.5d0*(eta(i,j,k)*(cury(i,j,k)*bz(i,j,k)-curz(i,j,k)*by(i,j,k)) &
                   +eta(i+1,j,k)*(cury(i+1,j,k)*bz(i+1,j,k)-curz(i+1,j,k)*by(i+1,j,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at j+1/2
  else if(mdir == 2)then
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres =0.5d0*(eta(i,j,k)*(curz(i,j,k)*bx(i,j,k)-curx(i,j,k)*bz(i,j,k)) &
                   +eta(i,j+1,k)*(curz(i,j+1,k)*bx(i,j+1,k)-curx(i,j+1,k)*bz(i,j+1,k)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
! Flux at k+1/2
  else
     !$OMP PARALLEL DO &
     !$OMP PRIVATE(i,j,fres)
     do k=margin,kx-margin
        do j=margin,jx-margin
           do i=margin,ix-margin
              fres =0.5d0*(eta(i,j,k)*(curx(i,j,k)*by(i,j,k)-cury(i,j,k)*bx(i,j,k)) &
                   +eta(i,j,k+1)*(curx(i,j,k+1)*by(i,j,k+1)-cury(i,j,k+1)*bx(i,j,k+1)))
              fee_res(i,j,k) = fee(i,j,k)+fres
           end do
        end do
     end do
     !$OMP END PARALLEL DO
  end if

  end subroutine flux_calc__feres

  subroutine convert_ctop(gm,ro,dro,ee,rx,ry,rz,bxc,byc,bzc &
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
  real(8),intent(in) :: gm
  real(8),intent(inout) :: ro,pr
  real(8),intent(inout) :: ee
  real(8),intent(in) :: dro,rx,ry,rz
  real(8),intent(in) :: bxc,byc,bzc
  real(8),intent(inout) :: vx,vy,vz

  integer :: count
  real(8), parameter :: pbeta_min=1d-10
  real(8) :: vsq,pb,roinverse,igm,temppr,signpr,temp1,temp2

  ! aloy et al.(1999) for HD
  !real(8) :: ps,vxs,vys,vzs,vvs,gls,eps,ros,css2,fp,dfp

  ! Mignone and Bodo et al. (2006) for MHD
  real(8) :: lg,ilg,ww,ss,bsq,rsq,pg,dlg,dpg,fw,dfw
  real(8) :: bb,cc,wwm,wwp

  real(8), parameter :: emin = 1d-10
  integer, parameter :: cmax = 50

  igm = 1d0/(gm-1d0)
  !***
  ! Primitve Ricovery for MHD case : Mignone and Bodo (2006)
  !***
  ss  = rx*bxc + ry*byc + rz*bzc
  bsq = bxc**2 + byc**2 + bzc**2
  rsq = rx*rx+ry*ry+rz*rz

  !*** first guess using previous stepes val. ***!
  ! vsq = vx*vx + vy*vy + vz*vz
  ! lg  = dsqrt(1d0 + vsq) ! Lorentz Factor
  ! ww = ( ro + igm*gm*pr )*lg*lg
  !**********************************************!

  !*********** Hydrodynamics Case ***************!
  ! ww = dsqrt(rsq+1d-8)
  ! ww = ww + ww*1d-4
  !**********************************************!

  !*** first guess using Mignone and McKinney ***!
  bb = 4.d0 *( bsq - ee )
  cc = bsq*bsq + rsq - 2d0*bsq*ee
  wwp = (-bb + dsqrt(bb*bb - 4d0*3d0*cc))/6d0
  ww  = wwp
  !***********************************************!

  ! lg  = 1d0 - (ss*ss*(bsq+ 2d0*ww)+ww*ww*rsq)/(ww*ww*(ww+bsq)*(ww+bsq))
  ! lg  = 1d0/dsqrt(lg)
  ! if( isnan(lg) )then
  !   ww = 100d0
  !   lg  = 1d0 - (ss*ss*(bsq+ 2d0*ww)+ww*ww*rsq)/(ww*ww*(ww+bsq)*(ww+bsq))
  !   lg  = 1d0/dsqrt(lg)
  !   write(*,*) lg,ww,wwp
  ! endif
  
  do count=1,cmax
      lg  = 1d0 - (ss*ss*(bsq+ 2d0*ww)+ww*ww*rsq)/(ww*ww*(ww+bsq)*(ww+bsq))
      lg  = 1d0/dsqrt(lg)
      pg  = (ww-dro*lg)/(igm*gm*lg*lg)
      
      dlg = lg*lg*lg*(rsq*ww*ww*ww+3d0*ss*ss*ww*ww+3d0*ss*ss*bsq*ww+ss*ss*bsq*bsq)
      dlg = -1d0*dlg/(ww*ww*ww*(ww+bsq)*(ww+bsq)*(ww+bsq))

      dpg = (lg*(1d0+dro*dlg)-2d0*ww*dlg)/(gm*igm*lg*lg*lg)

      fw  = ww - pg + (1d0-0.5d0/(lg*lg))*bsq - 0.5d0*ss*ss/(ww*ww) - ee

!     if( dabs(fw) .lt. emin)then
      if( dabs((pg-pr)/pg) .lt. emin)then
        exit
      endif
      dfw = 1d0 - dpg + dlg*bsq/(lg*lg*lg) + ss*ss/(ww*ww*ww)
      ww  = ww - fw/dfw
      pr = pg
  end do
  if(count .eq. cmax+1)then
    write(*,*)"Stop Due to faile for Primitive recovery at HLLC"
    write(*,*) wwp,pg,lg,rx,ry,rz
    stop
  endif

  vx = lg*(rx + bxc*ss/ww)/(ww+bsq)
  vy = lg*(ry + byc*ss/ww)/(ww+bsq)
  vz = lg*(rz + bzc*ss/ww)/(ww+bsq)
  pr = pg
  ro = dro/lg

  end subroutine convert_ctop

end module flux_calc
