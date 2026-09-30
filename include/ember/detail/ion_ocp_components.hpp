#pragma once
// Liquid Coulomb free energies from Potekhin & Chabrier (2013),
// https://arxiv.org/abs/1212.3405, appendices A, B.1 and C.1.
// Derivatives come from the common free energy. Quantum ions are supplied by ion_quantum.cpp.
#include "ember/detail/taylor3.hpp"
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
namespace ember::detail::ioffe {
// cgs atomic constants.
inline constexpr double bohr=5.29177210903e-9,hartree_k=315775.02480407,aum=1822.88848;
// x_r = p_F/(m_e c) = (9 pi/4)^(1/3) alpha / r_s.
inline const double relativity_x_rs=std::cbrt(9*std::numbers::pi/4)/137.035999084;
template<class J>J pow(const J& x,double p){return exp(p*log(x));}

// One-component-plasma parameters of an ion species (A,Z) at electron density n_e and temperature T.
template<class J> struct Plasma { J rs,gami,tpt,rsi,x; double Z,A; };
template<class J> Plasma<J> plasma(const J& logT,const J& logne,double A,double Z) {
  constexpr double pi=std::numbers::pi;
  using std::exp,std::sqrt;   // also usable with J=double (scalar domain checks)
  const J rs=exp((std::log(3/(4*pi))-logne)/3)/bohr;           // a_e/a_0
  const J gami=std::pow(Z,5./3)*hartree_k/rs*exp(-logT);
  const J tpt=gami/sqrt(rs)*(std::sqrt(3/(aum*A))/std::pow(Z,7./6));
  const J ratio=gami/tpt;
  return {rs,gami,tpt,3*ratio*ratio,relativity_x_rs/rs,Z,A};
}
// FITION9: classical OCP liquid (ii part).
template<class J> J fition9(const J& g) {
  constexpr double A1=-.907347,A2=.62849,C1=.0045,G1=170.,C2=-8.4e-5,SQ32=.8660254038;
  constexpr double G2=.0037;   // PC2013 eq. B.1: A1,A2,A3; B1=C1, B2=G1, B3=C2, B4=G2
  const double A3=-SQ32-A1/std::sqrt(A2);
  const J sg=sqrt(g);
  const J F0=A1*(sqrt(g*(A2+g))-A2*log(sqrt(g/A2)+sqrt(1+g/A2)))+2*A3*(sg-atan(sg));
  return F0+C1*(g-G1*log1p(g/G1))+C2/2*log1p(g*g/G2);
}
// Ideal ion gas (without mixing entropy), as in EOSFI22.
template<class J> J ideal_ion(const Plasma<J>& p){return 1.5*log(p.tpt*p.tpt/p.gami)-1.323515;}
// FSCRliq8: electron-ion screening (polarization) in the liquid.
template<class J> J fscr_liquid(const Plasma<J>& p) {
  const double Z=p.Z,ZLN=std::log(Z),Z13=std::cbrt(Z);
  const J game=p.gami/std::pow(Z,5./3);
  const J& rs=p.rs;const J& x=p.x;
  // PC2013 eq. C.1: c_DH, c_TF, a (P01), b (P03), nu (PTX), g1 (cor0), g2 (cor1), h1.
  const double CDH=Z/std::sqrt(3.)*(std::pow(1+Z,1.5)-std::pow(Z,1.5)-1);
  const double CTF=Z*Z*(18./175)*std::cbrt(144/(std::numbers::pi*std::numbers::pi))*(Z13-1+.2/std::sqrt(Z13)),
      P01=1.11*std::exp(.475*ZLN),P03=.2+.078*ZLN*ZLN,PTX=1.16+.08*ZLN,P1=(Z-1)/9.,
      Q1=.18/std::sqrt(std::sqrt(Z)),Q2=.2+.37/std::sqrt(Z);
  const J sqg=sqrt(game),tx=pow(game,PTX);
  const J ty1=1/(1e-3*Z*Z+2*game),ty2=1+6*rs*rs,rs3=rs*rs*rs;
  const J cor1=1+P1*rs3/ty2*(1+ty1);
  const J u0=.78*sqrt(game/Z)*rs3,d0=game*(Z*Z*Z)+21*rs3,cor0=1+u0/d0;
  const J h1=(1+x*x/5.)/(1+Q1*x+Q2*x*x);
  const J up=CDH*sqg+P01*CTF*tx*cor0*h1;
  const J dn=1+(P03*sqg+P01/rs*tx*cor1)/sqrt(1+x*x);
  return -up/dn*game;
}
// EXCOR7: electron exchange-correlation per ELECTRON (units of kT), identical in both ion phases.
// Perrot & Dharma-wardana non-relativistic Hartree-Fock fit with Ichimaru-type correlation.
template<class J> J excor7(const J& rs,const J& game) {
  // PC2013 eq. A.1 (a, b, d, f, g; code names A, B, D, E and C = g f). theta = T/T_F.
  const J th=(2/std::pow(9*std::numbers::pi/4,2./3))*rs/game,sqth=sqrt(th),th2=th*th,th3=th2*th,th4=th3*th;
  auto tanh_=[](const J& v){const J e=exp(-2*v);return (1-e)/(1+e);};
  const J t1=th.value()>.005?tanh_(1/th):J(1),t2=th.value()>.005?tanh_(1/sqth):J(1);
  const J A=std::cbrt(9/(4*std::numbers::pi*std::numbers::pi))*(.75+3.04363*th2-.09227*th3+1.7035*th4)/(1+8.31051*th2+5.1105*th4)*t1;
  const J B=sqth*t2*(.341308+12.0708*th2+1.148889*th4)/(1+10.495346*th2+1.326623*th4);
  const J D=sqth*t2*(.614925+16.996055*th2+1.489056*th4)/(1+10.10935*th2+1.22184*th4);
  const J E=th*t1*(.539409+2.522206*th2+.178484*th4)/(1+2.555501*th2+.146319*th4);
  const J C=(.872496+.025248*exp(-1/th))*E;
  const J discr=sqrt(4*E-D*D);
  const J S1=-C/E*game,B2=B-C*D/E,sqge=sqrt(game);
  const J S2=-2/E*B2*sqge;
  const J R3=E*game+D*sqge+1,B3=A-C/E,C3=(D/E*B2-B3)/E;
  const J S3=C3*log(R3);
  const J B4=2-D*D/E,C4=2*E*sqge+D;
  const J S4=2/E/discr*(D*B3+B4*B2)*(atan(C4/discr)-atan(D/discr));
  return S1+S2+S3+S4;
}
template<class J> J classical_liquid(const Plasma<J>& p) {
  return fition9(p.gami)+ideal_ion(p)+fscr_liquid(p);
}
} // namespace ember::detail::ioffe
