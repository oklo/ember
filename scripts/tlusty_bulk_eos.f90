module ember_atmosphere_bulk
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: dp=real64
  integer :: nt=0,np=0
  real(dp) :: xh
  real(dp),allocatable :: ts(:),ps(:),table(:,:,:,:)
  logical,allocatable :: cells(:,:)
  logical :: loaded=.false.,enabled=.false.
contains
  subroutine load_bulk
    character(len=1024) :: path
    character(len=80) :: header
    integer :: u,i,j,k,status,version
    integer,allocatable :: mask(:,:)
    if(loaded) return
    loaded=.true.
    call get_environment_variable('EMBER_ATMOSPHERE_BULK_EOS',path,status=status)
    if(status/=0.or.len_trim(path)==0) return
    open(newunit=u,file=trim(path),status='old',action='read',iostat=status)
    if(status/=0) error stop 'EMBER bulk EOS: cannot open table'
    read(u,'(a)')header
    select case(trim(header))
    case('EMBER_ATMOSPHERE_BULK_EOS 1'); version=1
    case('EMBER_ATMOSPHERE_BULK_EOS 2'); version=2
    case default; error stop 'EMBER bulk EOS: table format'
    end select
    read(u,*)nt,np,xh
    if(nt<2.or.np<2.or.xh<=0.or.xh>=1) error stop 'EMBER bulk EOS: invalid dimensions'
    allocate(ts(nt),ps(np),table(4,3,nt,np),cells(nt-1,np-1))
    cells=.true.
    read(u,*)ts
    read(u,*)ps
    do j=1,np
      do i=1,nt
        read(u,*)((table(k,status,i,j),k=1,4),status=1,3)
      end do
    end do
    if(version==2)then
      read(u,'(a)')header
      if(trim(header)/='cells')error stop 'EMBER bulk EOS: missing cell mask'
      allocate(mask(nt-1,np-1));read(u,*)mask
      if(any(mask<0).or.any(mask>1))error stop 'EMBER bulk EOS: invalid cell mask'
      cells=mask==1
    end if
    close(u)
    if(.not.all(ieee_is_finite(ts)).or..not.all(ieee_is_finite(ps)).or. &
       .not.all(ieee_is_finite(table)))error stop 'EMBER bulk EOS: nonfinite table'
    if(any(ts(2:)<=ts(:nt-1)).or.any(ps(2:)<=ps(:np-1))) error stop 'EMBER bulk EOS: unordered axes'
    enabled=.true.
    write(6,*)'EMBER NONIDEAL BULK EOS: ',trim(path),xh
  end subroutine
  subroutine cubic(u,b)
    real(dp),intent(in)::u
    real(dp),intent(out)::b(4)
    b=[1-3*u*u+2*u**3,3*u*u-2*u**3,u-2*u*u+u**3,-u*u+u**3]
  end subroutine
  subroutine bulk_value(t,p,rho,s,e)
    real(dp),intent(in)::t,p
    real(dp),intent(out)::rho,s,e
    real(dp)::lt,lp,dt,dpp,bt(4),bp(4),c(4,4),v(3)
    integer :: i,j,k,ii,jj
    call load_bulk
    if(.not.enabled) error stop 'EMBER bulk EOS: provider disabled'
    if(.not.ieee_is_finite(t).or..not.ieee_is_finite(p).or.t<=0.or.p<=0) &
      error stop 'EMBER bulk EOS: invalid state'
    lt=log(t);lp=log(p)
    if(lt<ts(1).or.lt>ts(nt).or.lp<ps(1).or.lp>ps(np))then
      write(6,*)'EMBER bulk EOS outside table: ',t,p
      error stop 'EMBER bulk EOS: unsupported state'
    end if
    i=1;j=1
    do while(i<nt-1.and.lt>ts(i+1));i=i+1;end do
    do while(j<np-1.and.lp>ps(j+1));j=j+1;end do
    ! At a shared edge or corner, use any containing supported cell.
    if(.not.cells(i,j))then
      search: do jj=j,min(j+1,np-1)
        if(lp<ps(jj).or.lp>ps(jj+1))cycle
        do ii=i,min(i+1,nt-1)
          if(lt<ts(ii).or.lt>ts(ii+1))cycle
          if(cells(ii,jj))then
            i=ii;j=jj
            exit search
          end if
        end do
      end do search
    end if
    if(.not.cells(i,j))then
      write(6,*)'EMBER bulk EOS masked cell: ',t,p
      error stop 'EMBER bulk EOS: unsupported cell'
    end if
    dt=ts(i+1)-ts(i);dpp=ps(j+1)-ps(j)
    call cubic((lt-ts(i))/dt,bt);call cubic((lp-ps(j))/dpp,bp)
    do k=1,3
      c(:2,:2)=table(1,k,i:i+1,j:j+1)
      c(3:,:2)=table(2,k,i:i+1,j:j+1)*dt
      c(:2,3:)=table(3,k,i:i+1,j:j+1)*dpp
      c(3:,3:)=table(4,k,i:i+1,j:j+1)*dt*dpp
      v(k)=exp(dot_product(bt,matmul(c,bp)))
    end do
    rho=v(1);s=v(2);e=v(3)
  end subroutine
end module
