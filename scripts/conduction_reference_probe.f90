! Optional external Ioffe conduct21.f comparison; no Fortran in normal builds.
program probe
implicit none
real(8) :: t,rho,z,a,k,s,q,st,kt,qt,sh,kh,qh,enh
integer :: ios,mode,iz
do
 read(*,*,iostat=ios) t,rho,z,a,mode
 if(ios/=0)exit
 call condconv(t*1d-6,rho,0d0,z,a,0d0,s,k,q,st,kt,qt,sh,kh,qh)
 enh=1d0;iz=nint(z)
 if(iz==1.or.iz==2)then
  if(mode==0)call blouin20(iz,rho,t*1d-6,enh)
  if(mode==1)call blouin1(iz,rho,t*1d-6,enh)
 endif
 write(*,'(es25.16e3)')k*enh
enddo
end program
