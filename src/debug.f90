module mod_debug
  !
  use mpi
  use mod_common_mpi, only: myid,ierr,coord
  use mod_param     , only: dims
  !
  implicit none
  !
  private
  public chkmean,chk_helmholtz
  !
  contains
  !
  subroutine chkmean(n,qmin,dzlzi,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    !
    integer, intent(in ), dimension(3)                    :: n
    integer, intent(in )                                  :: qmin
    real(8), intent(in ), dimension(-2:)                  :: dzlzi
    real(8), intent(in ), dimension(-qmin:,-qmin:,-qmin:) :: p
    real(8), intent(out)                                  :: mean
    !
    integer :: i,j,k
    mean = 0.d0
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,p,dzlzi) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mean = mean + p(i,j,k)*dzlzi(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/(1.d0*n(1)*dims(1)*n(2)*dims(2))
    !
    return
  end subroutine chkmean
  !
  subroutine chk_helmholtz(n,dli,dzci,dzfi,alpha,fp,fpp,bc,c_or_f,diffmax)
    !
    ! this subroutine checks if the implementation of implicit diffusion is
    ! correct
    !
    implicit none
    !
    integer,          intent(in ), dimension(3)        :: n
    real(8),          intent(in ), dimension(2)        :: dli
    real(8),          intent(in )                      :: alpha
    real(8),          intent(in ), dimension(-2:)      :: dzfi,dzci
    real(8),          intent(in ), dimension(0:,0:,0:) :: fp,fpp
    character(len=1), intent(in ), dimension(0:1,3)    :: bc
    character(len=1), intent(in ), dimension(3)        :: c_or_f
    real(8),          intent(out)                      :: diffmax
    !
    real(8) :: val
    integer :: i,j,k,im,ip,jm,jp,km,kp
    integer :: idir
    integer, dimension(3) :: q
    !integer :: ii,jj
    !
    q(:) = 0
    do idir = 1,3
      if(bc(1,idir).ne.'P'.and.c_or_f(idir).eq.'f') q(idir) = 1
    enddo
    select case(c_or_f(3))
    !
    ! need to compute the maximum difference!
    !
    case('c')
      diffmax = 0.d0
      do k=1,n(3)-q(3)
        kp = k + 1
        km = k - 1
        do j=1,n(2)-q(2)
          jp = j + 1
          jm = j - 1
          do i=1,n(1)-q(1)
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(ip,j,k)-2.d0*fpp(i,j,k)+fpp(im,j,k))*(dli(1)**2) + &
                  (fpp(i,jp,k)-2.d0*fpp(i,j,k)+fpp(i,jm,k))*(dli(2)**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k ))*dzci(k ) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzci(km))*dzfi(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !ii = coord(1)*n(1)+i
            !jj = coord(2)*n(2)+j
            !if(abs(val-fp(i,j,k)).gt.1.e-8) print*, 'Large difference : ', val-fp(i,j,k),ii,jj,k
          enddo
        enddo
      enddo
    case('f')
      diffmax = 0.d0
      do k=1,n(3)-q(3)
        kp = k + 1
        km = k - 1
        do j=1,n(2)-q(2)
          jp = j + 1
          jm = j - 1
          do i=1,n(1)-q(1)
            ip = i + 1
            im = i - 1
            val =  fpp(i,j,k)+(1./alpha)*( &
                  (fpp(ip,j,k)-2.d0*fpp(i,j,k)+fpp(im,j,k))*(dli(1)**2) + &
                  (fpp(i,jp,k)-2.d0*fpp(i,j,k)+fpp(i,jm,k))*(dli(2)**2) + &
                 ((fpp(i,j,kp)-fpp(i,j,k ))*dzfi(kp) - &
                  (fpp(i,j,k )-fpp(i,j,km))*dzfi(k ))*dzci(k) )
            val = val*alpha
            diffmax = max(diffmax,abs(val-fp(i,j,k)))
            !ii = coord(1)*n(1)+i
            !jj = coord(2)*n(2)+j
            !if(abs(val-fp(i,j,k)).gt.1.e-8) print*, 'Large difference : ', val,fp(i,j,k),ii,jj,k
          enddo
        enddo
      enddo
    end select
    call mpi_allreduce(MPI_IN_PLACE,diffmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
    !
    return
  end subroutine chk_helmholtz
  !
end module mod_debug
