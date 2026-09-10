! Offline, explicitly forced-phase audit of the unmodified EOSFI22 source.
! This isolates the ionic + electron interaction terms from ELECT11's
! approximate density inversion. It excludes ideal electrons and radiation.
! Input per line: phase[0/1], Z, A, rho[g/cm3], T[K].
program eip_phase_probe
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: phase, ios
  real(real64) :: z,a,rho,t,ne,ni,rs,gamma,tpt,rsp
  real(real64) :: f1,u1,p1,s1,cv1,pt1,pr1,f2,u2,p2,s2,cv2,pt2,pr2
  real(real64), parameter :: kb=1.380649e-16_real64,mu=1.66053906660e-24_real64
  real(real64), parameter :: bohr=5.29177210903e-9_real64,hartree_k=315775.02480407_real64
  real(real64), parameter :: pi=acos(-1._real64),aum=1822.88848_real64
  do
    read(*,*,iostat=ios) phase,z,a,rho,t
    if(ios<0) exit
    if(ios/=0 .or. .not.ieee_is_finite(z+a+rho+t)) error stop 'invalid query'
    if(phase<0 .or. phase>1 .or. z<1 .or. a<z .or. rho<=0 .or. t<=0) error stop 'invalid state'
    ni=rho/(a*mu);ne=ni*z*bohr**3
    rs=(3/(4*pi*ne))**(1._real64/3)
    gamma=z**(5._real64/3)*hartree_k/(rs*t)
    tpt=gamma/sqrt(rs)*sqrt(3/(aum*a))/z**(7._real64/6)
    rsp=kb/(a*mu)
    call EOSFI22(phase,a,z,rs,gamma, f1,u1,p1,s1,cv1,pt1,pr1, f2,u2,p2,s2,cv2,pt2,pr2)
    ! F,U,S,P,cv,dP/dlnT,dP/dlnrho in cgs; r_s,Gamma,Tp/T.
    ! The chosen phase is an input diagnostic, not a stable-phase solution.
    write(*,'(10(es25.16e3,1x))') f2*rsp*t,u2*rsp*t,s2*rsp,p2*ni*kb*t, &
      cv2*rsp,pt2*ni*kb*t,pr2*ni*kb*t,rs,gamma,tpt
    flush(6)
  end do
end program eip_phase_probe
