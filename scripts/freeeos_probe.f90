program probe
  use mod_free_eos, only: free_eos
  use mod_free_eos_types, only: fp_kind
  implicit none
  real(fp_kind) :: eps(20), lr, lt, fl, t, rho, rl, p, pl, cf, cp, s, sf, st
  real(fp_kind) :: grada, rtp, qe, qv, rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta
  real(fp_kind) :: gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2
  real(fp_kind) :: pressure(3), energy(3), entropy(3)
  integer :: ios, iterations, info, opt, modif, ion
  read(*,*) eps
  read(*,*) opt, modif, ion
  do
    read(*,*,iostat=ios) lr,lt
    if (ios /= 0) exit
    call free_eos(0, opt, modif, ion, 2, eps, lr, lt, fl, &
      t, rho, rl, p, pl, cf, cp, s, sf, st, grada, rtp, &
      qe, qv, rmue, fh2, fhe2, fhe3, xmu1, xmu3, eta, &
      gamma1, gamma2, gamma3, h2rat, h2plusrat, lambda, gamma_e, sound2, &
      iterations, info, pressure=pressure,energy=energy,entropy=entropy)
    write(*,'(2(i10,1x),20(es25.16e3,1x))') info,iterations,rho,t,p,qe,s, &
      pressure(2:3),energy(2:3),entropy(2:3),cp,grada,rtp,xmu1,xmu3,h2rat,fh2,fhe2,fhe3
    flush(6)
  end do
end program probe
