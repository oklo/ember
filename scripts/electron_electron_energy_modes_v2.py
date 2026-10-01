"""Extended energy intervals for the checked Born/Pauli moment operator.

Separate from the frozen three-mode comparison. The collision geometry is
evaluated directly on all four momenta, providing an independent check of its
reduced polynomial algebra. No stellar evolution prescription is selected.
"""
import math
from functools import lru_cache
import numpy as np
from scipy.special import expit,roots_jacobi
from electron_electron_collision import energy_basis, rule
from electron_ion_born import coulomb_bracket


class PolynomialBasis:
    def __init__(self,eta,degree,*,statistics='fermi',points=None,normalization_tail=80.):
        if not isinstance(degree,int) or not 2<=degree<=10:
            raise ValueError('polynomial degree 2 through 10 required')
        # Resolve the narrow Fermi surface at larger degeneracy.
        centered=statistics=='fermi' and eta>=128
        if points is None:
            points=256 if centered else 384 if eta<=64 else 768 if eta<=256 else 1536
            if not centered and eta>1024:points=math.ceil(1536*(eta+normalization_tail)/(1024+normalization_tail))
        if not isinstance(points,int) or points<32:raise ValueError("at least 32 normalization points required")
        self.base=energy_basis(eta,statistics=statistics)
        self.degree=degree;self.statistics=statistics
        if not np.isfinite(normalization_tail) or not 50<=normalization_tail<=120:
            raise ValueError('controlled polynomial energy tail required')
        self.normalization_tail=normalization_tail
        length=max(eta,0)+normalization_tail
        if centered:
            y,w=_three_panel_rule(points,[-normalization_tail,-8.,8.,normalization_tail])
            x=eta+y
            f=expit(-y)*expit(y)
            weight=w*x**1.5*f/self.base.normalization
        else:
            node,w=roots_jacobi(points,0,1.5)
            x=(node+1)*length/2
            f=(expit(eta-x)*expit(x-eta) if statistics=='fermi' else np.exp(eta-x))
            weight=w*(length/2)**2.5*f/self.base.normalization
        u=(x-eta-self.base.mean)/math.sqrt(self.base.variance)
        q=[np.ones_like(x),u]
        self.alpha=[0.];self.beta=[0.,1.]
        for n in range(1,degree):
            a=float(np.sum(weight*u*q[n]*q[n]));self.alpha.append(a)
            residual=(u-a)*q[n]-self.beta[n]*q[n-1]
            b=math.sqrt(float(np.sum(weight*residual*residual)))
            self.beta.append(b);q.append(residual/b)
        q=np.asarray(q)
        self.gram_error=float(np.max(abs((q*weight)@q.T-np.eye(degree+1))))
        self.points=points
        if self.gram_error>2e-9:raise ValueError('polynomial orthogonality check failed')

    def values(self,x):
        u=(np.asarray(x)-self.base.eta-self.base.mean)/math.sqrt(self.base.variance)
        q=[np.ones_like(u),u]
        for n in range(1,self.degree):
            q.append(((u-self.alpha[n])*q[n]-self.beta[n]*q[n-1])/self.beta[n+1])
        result=np.asarray(q)
        result[1]/=math.sqrt(self.base.variance)
        return result

    def ion_matrix(self,bthermal,*,points=None,tail=None):
        if not np.isfinite(bthermal) or bthermal<=0:raise ValueError('positive screening required')
        tail=self.normalization_tail if tail is None else tail
        if not np.isfinite(tail) or not 50<=tail<=120:raise ValueError('controlled ion energy tail required')
        centered=self.statistics=='fermi' and self.base.eta>=128
        if points is None:
            points=256 if centered else 512 if self.base.eta<=256 else 1536
            if not centered and self.base.eta>1024:
                points=math.ceil(1536*(self.base.eta+tail)/(1024+tail))
        if not isinstance(points,int) or points<32:raise ValueError('at least 32 ion quadrature points required')
        if centered:
            y,w=_three_panel_rule(points,[-tail,-8.,8.,tail])
            x=self.base.eta+y
            occupation=expit(-y)*expit(y)
        else:
            momentum,weight=rule(points,0,math.sqrt(max(self.base.eta,0)+tail))
            x,w=momentum**2,2*momentum*weight
            occupation=(expit(self.base.eta-x)*expit(x-self.base.eta)
                        if self.statistics=='fermi' else np.exp(self.base.eta-x))
        p=self.values(x)
        return (p*(w*occupation*coulomb_bracket(bthermal*x)))@p.T


@lru_cache(maxsize=None)
def _unit_rule(n):
    return np.polynomial.legendre.leggauss(n)

def _three_panel_rule(n, edges):
    counts = (n//4, n-2*(n//4), n//4)
    nodes, weights = [], []
    for count, low, high in zip(counts, edges[:-1], edges[1:]):
        x, w = _unit_rule(count)
        nodes.append(low+(x+1)*(high-low)/2)
        weights.append(w*(high-low)/2)
    return np.concatenate(nodes), np.concatenate(weights)


def energy_mode_matrix(eta,b_thermal,*,degree=7,orders=(64,32,32,40,40),
                       statistics='fermi',tail=50.,normalization_tail=80.,quadrature_layout='uniform',split_delta=8.,
                       normalization_points=None):
    if quadrature_layout not in ['uniform','energy','angle','both']:raise ValueError('unknown quadrature layout')
    if not np.isfinite(split_delta) or split_delta<=0 or split_delta>=tail:raise ValueError('invalid split width')
    if not np.isfinite(b_thermal) or b_thermal<=0:raise ValueError('positive screening required')
    if len(orders)!=5 or any(not isinstance(n,int) or not 8<=n<=160 for n in orders):
        raise ValueError('bounded quadrature orders required')
    if not np.isfinite(tail) or not 35<=tail<=60:raise ValueError('controlled tail required')
    basis=PolynomialBasis(eta,degree,statistics=statistics,normalization_tail=normalization_tail,points=normalization_points)
    na,npair,nu,nt,nphi=orders
    energies,wa=rule(na,0,max(eta,0)+tail)
    if statistics=='fermi' and eta>split_delta and quadrature_layout in ['energy','both']:
        energies,wa=_three_panel_rule(na,[0.,eta-split_delta,eta+split_delta,eta+tail])
    angles,wpair=rule(npair,0,math.pi/2)
    original_u,original_wu=rule(nu,-1,1)
    u=original_u[:,None,None];wu=original_wu;su=np.sqrt(1-u*u)
    phi=(np.arange(nphi)+.5)*(2*math.pi/nphi)
    cp=np.cos(phi)[None,None,:];sp=np.sin(phi)[None,None,:]
    unit,wt=rule(nt,0,1);unit=unit[None,:,None]
    angular=wu[:,None,None]*wt[None,:,None]*(2*math.pi/nphi)
    screening=4/b_thermal
    result=np.zeros((degree,degree))
    for a,w_a in zip(energies,wa):
        r=math.sqrt(a)
        for psi,w_psi in zip(angles,wpair):
            co,si=math.cos(psi),math.sin(psi);G,g=r*co,2*r*si;g2=g*g
            if statistics=='fermi' and eta>split_delta and quadrature_layout in ['angle','both']:
                width=split_delta/(G*g)
                un,uw=_three_panel_rule(nu,[-1.,-width,width,1.]) if width<1 else (original_u,original_wu)
                u=un[:,None,None];wu=uw;su=np.sqrt(1-u*u)
                angular=wu[:,None,None]*wt[None,:,None]*(2*math.pi/nphi)
            length=math.log1p(g2/(2*screening));t=unit*length
            y=(screening/g2)*np.expm1(t);ct=1-2*y;st=2*np.sqrt(y*(1-y))
            direct=np.exp(-t)/screening;other=1/(screening+g2*(1-y))
            jac=2*screening/g2*np.exp(t)*length
            b=G*g*u;bp=G*g*(u*ct+su*st*cp)
            pa,ma=basis.values(a+b)[1:],basis.values(a-b)[1:]
            pb,mb=basis.values(a+bp)[1:],basis.values(a-bp)[1:]
            even=pa+ma-pb-mb;pre=(pa-ma)/2;post=-(pb-mb)/2
            if statistics=='fermi':
                occupation=expit(eta-a-b)*expit(eta-a+b)*expit(a+bp-eta)*expit(a-bp-eta)
            else:occupation=math.exp(2*(eta-a))
            radial=8*a**2.5*co*co*si**3*w_a*w_psi
            weight=np.sqrt(radial*angular*jac*(direct*direct+other*other-direct*other)*occupation)
            components=[G*su*even+g*st*cp*post,
                        g*st*sp*post,
                        G*u*even+g*pre+g*ct*post]
            for component in components:
                v=(component*weight).reshape(degree,-1)
                result+=v@v.T
    matrix=np.zeros((degree+1,degree+1))
    matrix[1:,1:]=3/(8*math.pi*basis.base.normalization)*result
    np.linalg.cholesky(matrix[1:,1:])
    return dict(eta=eta,b_thermal=b_thermal,degree=degree,statistics=statistics,
                orders=list(orders),tail=tail,normalization_tail=normalization_tail,normalization_points=basis.points,quadrature_layout=quadrature_layout,split_delta=split_delta,collision_matrix=matrix.tolist(),
                electron_ion_matrix=basis.ion_matrix(b_thermal).tolist(),
                orthogonal_basis_gram_error=basis.gram_error,
                basis_alpha=basis.alpha,basis_beta=basis.beta)
