from math import *
import scipy
import matplotlib.pyplot as plt
import numpy as np
class zernike:
    def __init__(self,n,m,r0=1,ri=0):
        self.n = n
        self.m = m
        self.r0=r0
        self.ri=ri
        self.func = None
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
            ans = 0
            for p,c in zip(_p,_c):
                ans+= c*r**p
            angular = np.cos(_m*theta)*(_m>=0) + np.sin(_m*theta)*(_m<0)
            return ans*angular
        self.func = zern_function
        
    def plot(self):
        _f = self.func
        azimuths = np.radians(np.linspace(0, 360, 1000))
        zeniths = np.linspace(0, 1, 500)
        _r, _theta = np.meshgrid(zeniths, azimuths)
        values = _f(_r,_theta)
        fig,ax = plt.subplots(subplot_kw=dict(projection='polar'))
        ax.contourf(_theta,_r,values,cmap='jet')
        plt.show()
        
def get_two(n1,n2,m1,m2):
    a = zernike(n1,m1)
    b = zernike(n2,m2)
    a.get_c_p()
    b.get_c_p()
    a.get_function()
    b.get_function()
    return a.func,b.func
    
def check_ortho(z1,z2,ri=0,r0=1,theta_bounds=(0,2*pi),r_bounds=(0,1)):
    t_l,t_u = theta_bounds
    r_l,r_u = r_bounds
    function = lambda r,theta : z1((r-ri)/r0,theta)*z2((r-ri)/r0,theta)*r/pi
    integration = scipy.integrate.dblquad(function,t_l,t_u,r_l,r_u)
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