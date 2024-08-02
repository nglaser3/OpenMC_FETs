from math import *
import scipy
import matplotlib.pyplot as plt
import numpy as np
class zernike:
    '''
    n : order of Zernike Polynomial
    m : sub-order of Zernike Polynomial
    
    '''
    def __init__(self,n,m,r0=1,ri=0):
        self.n = n
        self.m = m
        self.r0=r0
        self.ri=ri
        self.func = None
        self.rad_func = None
        self.powers=[]
        self.coeffs=[]
        N = (2*(n+1)/(1+1*(m==0)))**.5
        
    def get_c_p(self):
        '''
        Obtains the coefficients and powers of zernike polynomial
        '''
        _m = self.m
        _n = self.n
        
        N = (2*(_n+1)/(1+1*(_m==0)))**.5
        upper_bound = int((_n-abs(_m))/2)
        for k in range(upper_bound+1):
            coeff = N*(-1)**k * factorial(_n-k) / factorial(k) / factorial(int((_n+abs(_m))/2) - k) / factorial(int((_n-abs(_m))/2)-k)
            self.coeffs.append(coeff)
            power = _n-2*k
            self.powers.append(power)
    def get_function(self):
        def zern_function(r,theta):
            _p,_c,_m = self.powers,self.coeffs,self.m
            _r0,_ri = self.r0,self.ri
            ans = 0
            for p,c in zip(_p,_c):
                ans+= c*((r-_ri)/(_r0-_ri))**p
            angular = np.cos(_m*theta)*(_m>=0) + np.sin(_m*theta)*(_m<0)
            return ans*angular*(r-_ri>=0)
        self.func = zern_function

    def get_radial(self):
        def zern_radial(r):
            _p,_c,_m = self.powers,self.coeffs,self.m
            _r0,_ri = self.r0,self.ri
            ans = 0
            for p,c in zip(_p,_c):
                ans+= c*((r-_ri)/(_r0-_ri))**p
            return ans*(r-_ri>=0)
        self.rad_func = zern_radial
    def plot(self):
        _f = self.func
        _r0,_ri=self.r0,self.ri
        azimuths = np.linspace(0,2*pi , 1000)
        zeniths = np.linspace(0, _r0, 500)
        _r, _theta = np.meshgrid(zeniths, azimuths)
        values = _f(_r,_theta)
        fig,ax = plt.subplots(subplot_kw=dict(projection='polar'))
        contour = ax.contourf(_theta,_r,values,cmap='jet')
        fig.colorbar(contour)
        plt.show()
    def plot_rad(self,vline=False,xloc = 0,show=True):
        _fr = self.rad_func
        _n,_m = self.n,self.m
        _r0,_ri=self.r0,self.ri
        _r = np.linspace(0,_r0,1000)
        plt.plot(_r,_fr(_r),label='R$^{}_{}$'.format(_m,_n))
        if vline:
            plt.axvline(xloc,color='k',linestyle=(0,(5,3)))
        plt.axhline(0,color='k')
        plt.legend()
        if show:
            plt.grid()
            plt.show()
        
def get_two(n1,n2,m1,m2,rmin=0,rmax=1):
    a = zernike(n1,m1,ri=rmin,r0=rmax)
    b = zernike(n2,m2,ri=rmin,r0=rmax)
    a.get_c_p()
    b.get_c_p()
    a.get_function()
    b.get_function()
    return a.func,b.func
    
def check_ortho(z1,z2,ri=0,r0=1,theta_bounds=(0,2*pi),r_bounds=(0,1)):
    t_l,t_u = theta_bounds
    r_l,r_u = r_bounds
    function = lambda r,theta : z1((r-ri)/r0,theta)*z2((r-ri)/r0,theta)*(r-ri)/pi/(r0-ri)
    integration = scipy.integrate.dblquad(function,t_l,t_u,r_l,r_u)
    return integration
    
def check_ortho_r(z1,z2,ri=0,r0=1,r_bounds=(0,1)):
    r_l,r_u = r_bounds
    function = lambda r : z1(r)*z2(r)*r/pi
    integration = scipy.integrate.quad(function,r_l,r_u)
    return integration

def randomizer(max_order):
    n_1,n_2 = randint(1,max_order),randint(1,max_order)
    def m(n):
        m = [-n]
        param = 1
        while param>0:
            m.append(m[-1]+2)
            param = n-m[-1]
            
        return m
    _m1,_m2 = m(n_1),m(n_2)
    m_1,m_2 = _m1[randint(0,len(_m1)-1)],_m2[randint(0,len(_m2)-1)]
    return n_1,n_2,m_1,m_2