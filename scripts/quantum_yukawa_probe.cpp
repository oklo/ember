// Independent radial Schrodinger integration for attractive Yukawa scattering.
// Radius is in screening lengths. a=k*lambda; c=2*m*Z*e^2*lambda/hbar^2.
// u'' + [a^2-l(l+1)/r^2+c*exp(-r)/r]u=0.
// Propagate ratios of Numerov amplitudes to avoid exponential growth/overflow.
// Output two adjacent normalized amplitudes for matching to free waves in Python.
// This is an offline numerical probe, not a selected stellar transport model.
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

double regular_factor(double radius, double a, double c, int ell) {
    // u=r^(ell+1) h. Frobenius series for h regular at the Coulomb singularity.
    std::array<long double, 40> coefficients{};
    coefficients[0]=1;
    long double sum=1, power=1;
    for(int n=1;n<40;++n){
        long double convolution=0, exponential=1;
        for(int j=0;j<n;++j){
            convolution+=exponential*coefficients[n-1-j];
            exponential*=-1.0L/(j+1);
        }
        coefficients[n]=-(c*convolution+(n>=2?a*a*coefficients[n-2]:0))
                        /(static_cast<long double>(n)*(n+2*ell+1));
        power*=radius;sum+=coefficients[n]*power;
    }
    if(!std::isfinite(sum))throw std::runtime_error("regular series overflow");
    return static_cast<double>(sum);
}

int main(){
    try{
        std::cout<<std::setprecision(17);
        double a,c,step,rmax;int lmax;
        while(std::cin>>a>>c>>step>>rmax>>lmax){
            if(!(std::isfinite(a)&&a>0&&std::isfinite(c)&&c>=0&&c<=200&&
                 std::isfinite(step)&&step>0&&std::isfinite(rmax)&&rmax>=12&&
                 lmax>=1&&lmax<=900))throw std::invalid_argument("unsupported scattering inputs");
            const int count=static_cast<int>(std::ceil(rmax/step));
            if(count<100||count>2000000||count<lmax+4)throw std::invalid_argument("unsupported radial grid");
            const double h=rmax/count,h2=h*h;
            std::vector<double> inverse_radius_squared(count+1),potential(count+1);
            for(int n=1;n<=count;++n){
                const double r=n*h;
                inverse_radius_squared[n]=1/(r*r);
                potential[n]=a*a+c*std::exp(-r)/r;
            }
            std::cout<<a<<' '<<c<<' '<<h<<' '<<rmax<<' '<<lmax;
            for(int ell=0;ell<=lmax;++ell){
                const double centrifugal=static_cast<double>(ell)*(ell+1);
                const int first=std::max(1,(ell+2)/2);
                auto factor=[&](int n){return 1+h2*(potential[n]-centrifugal*inverse_radius_squared[n])/12;};
                const double r1=first*h,r2=(first+1)*h;
                // Starting beyond the coarse centrifugal singularity keeps
                // the Numerov coefficient positive; the series supplies the
                // regular boundary condition at both starting points.
                double ratio=std::pow(static_cast<double>(first+1)/first,ell+1)
                             *regular_factor(r2,a,c,ell)/regular_factor(r1,a,c,ell)
                             *factor(first+1)/factor(first);
                for(int n=first+1;n<count;++n){
                    const double qh=h2*(potential[n]-centrifugal*inverse_radius_squared[n]);
                    ratio=2*(1-5*qh/12)/(1+qh/12)-1/ratio;
                }
                const double raw=ratio*factor(count-1)/factor(count);
                if(!std::isfinite(raw))throw std::runtime_error("nonfinite Numerov ratio");
                const double scale=std::max(1.,std::abs(raw));
                std::cout<<' '<<raw/scale<<' '<<1/scale;
            }
            std::cout<<'\n';
        }
        if(!std::cin.eof())throw std::invalid_argument("malformed scattering query");
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
