module mod_fft
  !
  use iso_c_binding , only: C_INT
  use mod_common_mpi, only: ierr
  use mod_fftw_param
  use mod_param     , only: dims
  !
  !$ use omp_lib
  private
  public fftini,fftend,fftd,ffti
  !
  contains
  !
  subroutine fftini(nx,ny,nz,bcxy,c_or_f,arrplan,normfft)
    !
    implicit none
    !
    integer         , intent(in )                   :: nx,ny,nz
    character(len=1), intent(in ), dimension(0:1,2) :: bcxy
    character(len=1), intent(in ), dimension(2)     :: c_or_f
    type(C_PTR)     , intent(out), dimension(2,2)   :: arrplan
    real(8)         , intent(out)                   :: normfft
    !
    real(8), dimension(nx,ny/dims(1),nz/dims(2)) :: arrx
    real(8), dimension(nx/dims(1),ny,nz/dims(2)) :: arry
    type(C_PTR) :: plan_fwd_x,plan_bwd_x, &
                   plan_fwd_y,plan_bwd_y
    type(fftw_iodim), dimension(1) :: iodim
    type(fftw_iodim), dimension(2) :: iodim_howmany
    integer :: kind_fwd,kind_bwd
    real(8), dimension(2) :: norm
    integer(C_INT) :: nx_x,ny_x,nz_x, &
                      nx_y,ny_y,nz_y
    integer :: ix,iy
    !$ call dfftw_init_threads(ierr)
    !$ call dfftw_plan_with_nthreads(omp_get_max_threads())
    !
    ! fft in x
    !
    ! prepare plans with guru interface
    !
    nx_x = nx
    ny_x = ny/dims(1)
    nz_x = nz/dims(2)
    nx_y = nx/dims(1)
    ny_y = ny
    nz_y = nz/dims(2)
    !
    normfft = 1.d0
    ix = 0 
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if(bcxy(0,1)//bcxy(1,1).eq.'DD'.and.c_or_f(1).eq.'f') ix = 1
    iodim(1)%n  = nx_x-ix
    iodim(1)%is = 1
    iodim(1)%os = 1
    iodim_howmany(1)%n  = ny_x
    iodim_howmany(1)%is = nx_x
    iodim_howmany(1)%os = nx_x
    iodim_howmany(2)%n  = nz_x
    iodim_howmany(2)%is = nx_x*ny_x
    iodim_howmany(2)%os = nx_x*ny_x
    call find_fft(bcxy(:,1),c_or_f(1),kind_fwd,kind_bwd,norm)
    plan_fwd_x=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrx,arrx,kind_fwd,FFTW_ESTIMATE)
    plan_bwd_x=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arrx,arrx,kind_bwd,FFTW_ESTIMATE)
    normfft = normfft*norm(1)*(nx_x+norm(2)-ix)
    !
    ! fft in y
    !
    ! prepare plans with guru interface
    !
    iy = 0
    ! size of transform reduced by 1 point with Dirichlet BC in face
    if(bcxy(0,2)//bcxy(1,2).eq.'DD'.and.c_or_f(2).eq.'f') iy = 1
    iodim(1)%n  = ny_y-iy
    iodim(1)%is = nx_y
    iodim(1)%os = nx_y
    iodim_howmany(1)%n  = nx_y
    iodim_howmany(1)%is = 1
    iodim_howmany(1)%os = 1
    iodim_howmany(2)%n  = nz_y
    iodim_howmany(2)%is = nx_y*ny_y
    iodim_howmany(2)%os = nx_y*ny_y
    call find_fft(bcxy(:,2),c_or_f(2),kind_fwd,kind_bwd,norm)
    plan_fwd_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arry,arry,kind_fwd,FFTW_ESTIMATE)
    plan_bwd_y=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arry,arry,kind_bwd,FFTW_ESTIMATE)
    normfft = normfft*norm(1)*(ny_y+norm(2)-iy)
    !
    normfft = normfft**(-1)
    arrplan(1,1) = plan_fwd_x
    arrplan(2,1) = plan_bwd_x
    arrplan(1,2) = plan_fwd_y
    arrplan(2,2) = plan_bwd_y
    !
    return
  end subroutine fftini
  !
  subroutine fftend(arrplan)
    !
    implicit none
    !
    type(C_PTR), intent(in), dimension(2,2) :: arrplan
    call dfftw_destroy_plan(arrplan(1,1))
    call dfftw_destroy_plan(arrplan(1,2))
    call dfftw_destroy_plan(arrplan(2,1))
    call dfftw_destroy_plan(arrplan(2,2))
    !$call dfftw_cleanup_threads(ierr)
    !
    return
  end subroutine fftend
  !
  subroutine fftd(plan,arr)
    !
    implicit none
    !
    type(C_PTR), intent(in   )                   :: plan 
    real(8)    , intent(inout), dimension(:,:,:) :: arr
    call dfftw_execute_r2r(plan,arr,arr)
    !
    return
  end subroutine fftd
  !
  subroutine ffti(plan,arr)
    !
    implicit none
    !
    type(C_PTR), intent(in   )                   :: plan 
    real(8)    , intent(inout), dimension(:,:,:) :: arr
    call dfftw_execute_r2r(plan,arr,arr)
    !
    return
  end subroutine ffti
  !
  subroutine find_fft(bc,c_or_f,kind_fwd,kind_bwd,norm)
    !
    implicit none
    !
    character(len=1), intent(in ), dimension(0:1) :: bc
    character(len=1), intent(in )                 :: c_or_f
    integer         , intent(out)                 :: kind_fwd,kind_bwd
    real(8)         , intent(out), dimension(2)   :: norm
    !
    if(c_or_f.eq.'c') then
      select case(bc(0)//bc(1))
      case('PP')
        kind_fwd = FFTW_R2HC
        kind_bwd = FFTW_HC2R
        norm = (/1.d0,0.d0/)
      case('NN')
        kind_fwd = FFTW_REDFT10
        kind_bwd = FFTW_REDFT01
        norm = (/2.d0,0.d0/)
      case('DD')
        kind_fwd = FFTW_RODFT10
        kind_bwd = FFTW_RODFT01
        norm = (/2.d0,0.d0/)
      case('ND')
        kind_fwd = FFTW_REDFT11
        kind_bwd = FFTW_REDFT11
        norm = (/2.d0,0.d0/)
      case('DN')
        kind_fwd = FFTW_RODFT11
        kind_bwd = FFTW_RODFT11
        norm = (/2.d0,0.d0/)
      end select
    elseif(c_or_f.eq.'f') then
      select case(bc(0)//bc(1))
      case('PP')
        kind_fwd = FFTW_R2HC
        kind_bwd = FFTW_HC2R
        norm = (/1.d0,0.d0/)
      case('NN')
        kind_fwd = FFTW_REDFT00
        kind_bwd = FFTW_REDFT00
        norm = (/2.d0,-1.d0/)
      case('DD')
        kind_fwd = FFTW_RODFT00
        kind_bwd = FFTW_RODFT00
        norm = (/2.d0,1.d0/)
      case('ND')
        kind_fwd = FFTW_REDFT10
        kind_bwd = FFTW_REDFT01
        norm = (/2.d0,0.d0/)
      case('DN')
        kind_fwd = FFTW_RODFT01
        kind_bwd = FFTW_RODFT10
        norm = (/2.d0,0.d0/)
      end select
    endif
    !
    return
  end subroutine find_fft
  !
end module mod_fft
