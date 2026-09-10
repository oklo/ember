! Offline reference to the unmodified Ioffe EOS EIP (June 2022).
! Outputs retain the source's inferred density and classical phase choice.
! This is not an Ember runtime EOS or a treatment of partially ionized matter.
program eip_probe
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: n, i, ios, phase
  real(real64) :: x(32), z(32), a(32), rho, temperature, tau
  real(real64) :: prad, ne, zmean, amean, z2, gamma, chi, quantum, p, u, s, cv, chir, chit
  real(real64), parameter :: kb=1.380649e-16_real64, mu=1.66053906660e-24_real64
  real(real64), parameter :: bohr=5.29177210903e-9_real64, hartree_over_k=315775.02480407_real64
  real(real64) :: ni, actual_rho, gas_constant
  do
    read(*,*,iostat=ios) n, rho, temperature
    if (ios<0) exit
    if (ios/=0 .or. n<1 .or. n>32 .or. .not.ieee_is_finite(rho+temperature)) error stop 'invalid query'
    if (rho<1e-19_real64 .or. rho>1e15_real64 .or. temperature<=0) error stop 'outside source density range'
    do i=1,n
      read(*,*) z(i),a(i),x(i)
      if (.not.ieee_is_finite(z(i)+a(i)+x(i))) error stop 'nonfinite mixture'
      if (z(i)<=0 .or. a(i)<z(i) .or. x(i)<=0) error stop 'invalid ion'
    end do
    if (abs(sum(x(:n))-1)>1e-12_real64) error stop 'number fractions must sum to one'
    tau=temperature/hartree_over_k
    call MELANGE9(n,x,z,a,rho,tau,prad,ne,zmean,amean,z2,gamma,chi,quantum,phase,p,u,s,cv,chir,chit)
    ni=ne/(bohr**3*zmean)
    actual_rho=ni*amean*mu
    gas_constant=kb/(amean*mu)
    ! rho_requested, rho_inferred, T, Pmaterial, Ematerial, Smaterial, cv,
    ! chi_rho, chi_T, Gamma_effective, Tplasma/T, electron_eta, phase(0/1).
    write(*,'(12(es25.16e3,1x),i2)') rho,actual_rho,temperature,p*ni*kb*temperature, &
      u*gas_constant*temperature,s*gas_constant,cv*gas_constant,chir,chit,gamma,quantum,chi,phase
    flush(6)
  end do
end program eip_probe
