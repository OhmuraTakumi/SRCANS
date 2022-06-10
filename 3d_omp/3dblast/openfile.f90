module openfile
  use mpi_setup
  use dac_header
  use const, only : input_dir,output_dir

  implicit none
  private

  public :: file_input, file_output, file_output_param

  character :: cno*4
  character :: cnond*4
contains

  subroutine file_input(nd,mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,&
       ix,jx,kx,igx,jgx,kgx,margin)

    integer,intent(in) :: nd,mpirank,ix,jx,kx,igx,jgx,kgx,margin
    real(8),intent(inout),dimension(ix,jx,kx) :: ro,pr,vx,vy,vz,bx,by,bz,phi,eta
    integer :: is,ie,js,je,ks,ke
    integer :: fh,ftype
    integer, dimension(3) :: gridsize,subsize,start
    integer(kind=mkind) :: disp    
    disp=0
    
    write(cnond,'(i4.4)') nd

    gridsize = (/igx,jgx,kgx/); subsize= (/ix,jx,kx/)
    start    = (/mpid%mpirank_3d(1)*(ix-2*margin), &
                 mpid%mpirank_3d(2)*(jx-2*margin), &
                 mpid%mpirank_3d(3)*(kx-2*margin)/)

    ! ro
    call mpi_file_open(mcomw, input_dir//cnond//'_ro.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,ro(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! pr
    call mpi_file_open(mcomw, input_dir//cnond//'_pr.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,pr(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! vx
    call mpi_file_open(mcomw, input_dir//cnond//'_vx.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,vx(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! vy
    call mpi_file_open(mcomw, input_dir//cnond//'_vy.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,vy(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! vz
    call mpi_file_open(mcomw, input_dir//cnond//'_vz.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,vz(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! bx
    call mpi_file_open(mcomw, input_dir//cnond//'_bx.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,bx(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! by
    call mpi_file_open(mcomw, input_dir//cnond//'_by.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,by(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! bz
    call mpi_file_open(mcomw, input_dir//cnond//'_bz.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,bz(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! phi
    call mpi_file_open(mcomw, input_dir//cnond//'_phi.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,phi(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)

    ! eta
    call mpi_file_open(mcomw, input_dir//cnond//'_eta.dac',amoder,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize,subsize,start,order, &
         mdp,ftype,merr)
    call mpi_type_commit(ftype,merr)        
    call mpi_file_set_view(fh,disp,mdp,ftype,"native",minfonull,merr)
    call mpi_file_read_all(fh,eta(1,1,1),ix*jx*kx,mdp,mstat,merr)
    call mpi_type_free(ftype,merr)
    call mpi_file_close(fh,merr)
     
  end subroutine file_input

  subroutine file_output(nd,mpirank,ro,pr,vx,vy,vz,bx,by,bz,phi,eta,&
       ix,jx,kx,igx,jgx,kgx,margin,mpisize_x,mpisize_y,mpisize_z)

    integer,intent(in) :: nd,mpirank,ix,jx,kx,igx,jgx,kgx,margin
    real(8),intent(in),dimension(ix,jx,kx) :: ro,pr,vx,vy,vz,bx,by,bz,phi,eta
    integer,intent(in) :: mpisize_x,mpisize_y,mpisize_z
    integer :: is,ie,js,je,ks,ke
    integer :: fh,ftype1,ftype2
    integer, dimension(3) :: gridsize1,subsize1,start1
    integer, dimension(3) :: gridsize2,subsize2,start2    
    integer(kind=mkind) :: disp

    is = margin+1
    ie = ix-margin
    if( mpid%mpirank_3d(1) .eq. 0)  is= 1
    if( mpid%mpirank_3d(1) .eq. (mpisize_x-1)) ie = ix
    js = margin+1
    je = jx-margin
    if( mpid%mpirank_3d(2) .eq. 0)  js=1
    if( mpid%mpirank_3d(2) .eq. (mpisize_y-1)) je = jx
    ks = margin+1
    ke = kx-margin
    if( mpid%mpirank_3d(3) .eq. 0)  ks=1
    if( mpid%mpirank_3d(3) .eq. (mpisize_z-1)) ke = kx

    gridsize1 = (/igx,jgx,kgx/)
    subsize1  = (/ie-is+1,je-js+1,ke-ks+1/)    
    start1    = (/mpid%mpirank_3d(1)*(ix-2*margin)+(is-1), &
                 mpid%mpirank_3d(2)*(jx-2*margin)+(js-1), &
                 mpid%mpirank_3d(3)*(kx-2*margin)+(ks-1)/)
    gridsize2 = (/ix,jx,kx/)
    subsize2  = (/ie-is+1,je-js+1,ke-ks+1/)        
    start2    = (/is-1,js-1,ks-1/)    
    
    disp = 0
    
    write(cnond,'(i4.4)') nd
    ! ro
    call mpi_file_open(mcomw, output_dir//cnond//'_ro.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,ro(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

    ! pr
    call mpi_file_open(mcomw, output_dir//cnond//'_pr.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,pr(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

    ! vx
    call mpi_file_open(mcomw, output_dir//cnond//'_vx.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,vx(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

    ! vy
    call mpi_file_open(mcomw, output_dir//cnond//'_vy.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,vy(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

    ! vz
    call mpi_file_open(mcomw, output_dir//cnond//'_vz.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,vz(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)    

    ! bx
    call mpi_file_open(mcomw, output_dir//cnond//'_bx.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,bx(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

    ! by
    call mpi_file_open(mcomw, output_dir//cnond//'_by.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,by(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)    

    ! bz
    call mpi_file_open(mcomw, output_dir//cnond//'_bz.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,bz(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)
    
    ! phi
    call mpi_file_open(mcomw, output_dir//cnond//'_phi.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,phi(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

    ! eta
    call mpi_file_open(mcomw, output_dir//cnond//'_eta.dac',amode,minfonull, fh,merr)
    call mpi_type_create_subarray(3,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(3,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,eta(1,1,1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)
  end subroutine file_output

  subroutine file_output_param(nd,dtout,tend,ix,jx,kx,igx,jgx,kgx,margin,mpisize &
       ,mpirank,mpisize_x,mpisize_y,mpisize_z,gm,x,y,z,dx,dy,dz)
    
    integer,intent(in) :: nd,ix,jx,kx,igx,jgx,kgx,margin,mpisize,mpirank
    integer,intent(in) :: mpisize_x,mpisize_y,mpisize_z
    real(8),intent(in),dimension(ix) :: x,dx
    real(8),intent(in),dimension(jx) :: y,dy
    real(8),intent(in),dimension(kx) :: z,dz
    real(8),intent(in) :: dtout,tend
    real(8),intent(in) :: gm
    integer :: mf_params
    integer :: mf_x,mf_y,mf_z
    integer :: is,ie,js,je,ks,ke
    integer :: fh,ftype1,ftype2
    integer, dimension(3) :: gridsize1,subsize1,start1
    integer, dimension(3) :: gridsize2,subsize2,start2    
    integer(kind=mkind) :: disp
    
    write(cnond,'(i4.4)') nd
    write(cno,'(i4.4)') mpirank

    if(mpirank .eq. 0)then
       mf_params=9
       open (mf_params,file=output_dir//'params_rank='//cno//'.txt',form='formatted')
    
       call dacputparamc(mf_params,'comment','model_shktb')
       call dacputparami(mf_params,'ix',ix)
       call dacputparami(mf_params,'jx',jx)
       call dacputparami(mf_params,'kx',kx)
       call dacputparami(mf_params,'igx',igx)
       call dacputparami(mf_params,'jgx',jgx)
       call dacputparami(mf_params,'kgx',kgx)
       call dacputparami(mf_params,'margin',margin)
       call dacputparamd(mf_params,'tend',tend)
       call dacputparami(mf_params,'mpi',1)
       call dacputparamd(mf_params,'dtout',dtout)
       call dacputparami(mf_params,'mpisize',mpisize)
       call dacputparami(mf_params,'mpirank',mpirank)
       call dacputparami(mf_params,'mpix',mpisize_x)
       call dacputparami(mf_params,'mpiy',mpisize_y)
       call dacputparami(mf_params,'mpiz',mpisize_z)
       call dacputparamd(mf_params,'x(1)',x(1))
       call dacputparamd(mf_params,'y(1)',y(1))
       call dacputparamd(mf_params,'z(1)',z(1))
       call dacputparamd(mf_params,'dx(1)',dx(1))
       call dacputparamd(mf_params,'dy(1)',dy(1))
       call dacputparamd(mf_params,'dz(1)',dz(1))
       call dacputparamd(mf_params,'gm',gm)
       close(mf_params)
    endif
 
    is = margin+1
    ie = ix-margin
    if( mpid%mpirank_3d(1) .eq. 0)  is= 1
    if( mpid%mpirank_3d(1) .eq. (mpisize_x-1)) ie = ix
    js = margin+1
    je = jx-margin
    if( mpid%mpirank_3d(2) .eq. 0)  js=1
    if( mpid%mpirank_3d(2) .eq. (mpisize_y-1)) je = jx
    ks = margin+1
    ke = kx-margin
    if( mpid%mpirank_3d(3) .eq. 0)  ks=1
    if( mpid%mpirank_3d(3) .eq. (mpisize_z-1)) ke = kx

! --- x-coordinate ----
    gridsize1 = igx
    subsize1  = ie-is+1
    start1 = mpid%mpirank_3d(1)*(ix-2*margin)+(is-1)

    gridsize2 = ix
    subsize2  = ie-is+1
    start2    = is-1

    disp=0
    call mpi_file_open(mcomw, output_dir//'x.dac',amode & 
            ,minfonull, fh,merr)
    call mpi_type_create_subarray(1,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)    
    call mpi_type_create_subarray(1,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,x(1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)    
    call mpi_file_close(fh,merr)

! --- y-coordinate ----
    gridsize1 = jgx
    subsize1  = je-js+1
    start1 = mpid%mpirank_3d(2)*(jx-2*margin)+(js-1)

    gridsize2 = jx
    subsize2  = je-js+1
    start2    = js-1

    call mpi_file_open(mcomw, output_dir//'y.dac' & 
            ,amode,minfonull, fh,merr)
    call mpi_type_create_subarray(1,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)
    call mpi_type_create_subarray(1,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,y(1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)
    call mpi_file_close(fh,merr)

! --- z-coordinate ----
    gridsize1 = kgx
    subsize1  = ke-ks+1
    start1 = mpid%mpirank_3d(3)*(kx-2*margin)+(ks-1)

    gridsize2 = kx
    subsize2  = ke-ks+1
    start2    = ks-1

    call mpi_file_open(mcomw, output_dir//'z.dac' & 
            ,amode,minfonull, fh,merr)
    call mpi_type_create_subarray(1,gridsize1,subsize1,start1,order, &
         mdp,ftype1,merr)
    call mpi_type_commit(ftype1,merr)
    call mpi_type_create_subarray(1,gridsize2,subsize2,start2,order, &
         mdp,ftype2,merr)
    call mpi_type_commit(ftype2,merr)
    call mpi_file_set_view(fh,disp,mdp,ftype1,"native",minfonull,merr)
    call mpi_file_write_all(fh,z(1),1,ftype2,mstat,merr)
    call mpi_type_free(ftype1,merr)
    call mpi_type_free(ftype2,merr)
    call mpi_file_close(fh,merr)

  end subroutine file_output_param

end module openfile
