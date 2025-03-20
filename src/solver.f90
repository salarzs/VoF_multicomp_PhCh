module mod_solver
  !
  use iso_c_binding, only: C_PTR
  use decomp_2d
  use mod_fft      , only: fftd,ffti
  use mod_param    , only: dims
  !
  implicit none
  ! 
  private
  public solver
  !
  contains
  !
  subroutine solver(n,arrplan,normfft,lambdaxy,a,b,c,bcz,c_or_f,p)
    !
    implicit none
    !
    integer         , intent(in   ), dimension(3)         :: n
    type(C_PTR)     , intent(in   ), dimension(2,2)       :: arrplan
    real(8)         , intent(in   )                       :: normfft
    real(8)         , intent(in   ), dimension(n(1),n(2)) :: lambdaxy
    real(8)         , intent(in   ), dimension(n(3))      :: a,b,c
    character(len=1), intent(in   ), dimension(0:1)       :: bcz
    character(len=1), intent(in   ), dimension(3)         :: c_or_f
    real(8)         , intent(inout), dimension(0:,0:,0:)  :: p
    !
    real(8), dimension(n(1)*dims(1),n(2)*dims(2)/dims(1),n(3)/dims(2)) :: px
    real(8), dimension(n(1)*dims(1)/dims(1),n(2)*dims(2),n(3)/dims(2)) :: py
    real(8), dimension(n(1)        ,n(2)                ,n(3)        ) :: pz
    !real(8), allocatable, dimension(:,:,:) :: px,py
    integer, dimension(3) :: ng
    integer :: q
    ng(:) = n(:)
    ng(1:2) = ng(1:2)*dims(1:2)
    !allocate(px(ng(1),ng(2)/dims(1),ng(3)/dims(2)))
    !allocate(py(ng(1)/dims(1),ng(2),ng(3)/dims(2)))
    !
    !$OMP WORKSHARE
    pz(:,:,:) = p(1:n(1),1:n(2),1:n(3))
    !$OMP END WORKSHARE
    call transpose_z_to_x(pz,px)
    !call transpose_z_to_y(pz,py)
    !call transpose_y_to_x(py,px)
    call fftd(arrplan(1,1),px) ! fwd transform in x
    !
    call transpose_x_to_y(px,py)
    call fftd(arrplan(1,2),py) ! fwd transform in y
    !
    call transpose_y_to_z(py,pz)
    q = 0
    if(c_or_f(3).eq.'f'.and.bcz(1).eq.'D') q = 1
    if(bcz(0)//bcz(1).eq.'PP') then
      call gaussel_periodic(n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
    else
      call gaussel(         n(1),n(2),n(3)-q,a,b,c,lambdaxy,pz)
    endif
    !
    call transpose_z_to_y(pz,py)
    call ffti(arrplan(2,2),py) ! bwd transform in y
    !
    call transpose_y_to_x(py,px)
    call ffti(arrplan(2,1),px) ! bwd transform in x
    !
    call transpose_x_to_z(px,pz)
    !call transpose_x_to_y(px,py)
    !call transpose_y_to_z(py,pz)
    !$OMP WORKSHARE
    pz(:,:,:) = pz(:,:,:)*normfft
    p(1:n(1),1:n(2),1:n(3)) = pz(:,:,:)
    !$OMP END WORKSHARE
    !deallocate(px,py)
    !
    return
  end subroutine solver
  !
  subroutine gaussel(nx,ny,n,a,b,c,lambdaxy,p)
    !
    implicit none
    !
    integer, intent(in   )                   :: nx,ny,n
    real(8), intent(in   ), dimension(:)     :: a,b,c
    real(8), intent(in   ), dimension(nx,ny) :: lambdaxy
    real(8), intent(inout), dimension(:,:,:) :: p
    !
    real(8), dimension(n) :: bb
    integer :: i,j
    !
    !solve tridiagonal system
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP PRIVATE(i,j,bb) &
    !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        bb(:) = b(1:n) + lambdaxy(i,j)
        call dgtsv_homebrewed(n,a,bb,c,p(i,j,1:n))
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !
    return
  end subroutine gaussel
  !
  subroutine gaussel_periodic(nx,ny,n,a,b,c,lambdaxy,p)
    implicit none
    integer, intent(in) :: nx,ny,n
    real(8), intent(in), dimension(:) :: a,b,c
    real(8), intent(in), dimension(nx,ny) :: lambdaxy
    real(8), intent(inout), dimension(:,:,:) :: p
    real(8), dimension(n) :: bb,p1,p2
    integer :: i,j,info
    !
    !solve tridiagonal system
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP PRIVATE(i,j,bb,p1,p2) &
    !$OMP SHARED(nx,ny,n,a,b,c,lambdaxy,p)
    !$OMP DO COLLAPSE(2)
    do j=1,ny
      do i=1,nx
        bb(:)  = b(:) + lambdaxy(i,j)
        p1(1:n-1) = p(i,j,1:n-1)
        call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-1),p1(1:n-1))
        p2(:) = 0.d0
        p2(1  ) = -a(1  )
        p2(n-1) = -c(n-1)
        call dgtsv_homebrewed(n-1,a(1:n-1),bb(1:n-1),c(1:n-1),p2(1:n-1))
        p(i,j,n) = (p(i,j,n) - c(n)*p1(1) - a(n)*p1(n-1)) / &
                   (bb(   n) + c(n)*p2(1) + a(n)*p2(n-1))
        p(i,j,1:n-1) = p1(1:n-1) + p2(1:n-1)*p(i,j,n)
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    !
    return
  end subroutine gaussel_periodic
  !
  subroutine dgtsv_homebrewed(n,a,b,c,p)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in   ), dimension(:) :: a,b,c
    real(8), intent(inout), dimension(:) :: p
    real(8), dimension(n) :: d
    real(8) :: z
    integer :: l
    !
    ! Gauss elimination
    !
    z = 1.d0/b(1)
    d(1) = c(1)*z
    p(1) = p(1)*z
    do l=2,n-1
      z    = 1.d0/(b(l)-a(l)*d(l-1))
      d(l) = c(l)*z
      p(l) = (p(l)-a(l)*p(l-1))*z
    enddo
    z = b(n)-a(n)*d(n-1)
    if(z.ne.0.d0) then
      p(n) = (p(n)-a(n)*p(n-1))/z
    else
      p(n) = 0.d0
    endif
    !
    ! backward substitution
    !
    do l=n-1,1,-1
      p(l) = p(l) - d(l)*p(l+1)
    enddo
    return
  end subroutine dgtsv_homebrewed
  !
end module mod_solver
