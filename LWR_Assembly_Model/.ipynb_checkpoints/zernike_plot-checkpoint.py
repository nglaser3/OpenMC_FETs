import openmc, numpy as np, matplotlib.pyplot as plt, os, glob
from matplotlib.animation import FuncAnimation
from assembly_build import fuel_outter_radius


statepoint = glob.glob('state*')
if len(statepoint) != 1:
    raise Exception('Make sure only one statepoint file is present')
else:
    statepoint = statepoint[0]


with openmc.StatePoint(statepoint) as sp:
    df = sp.get_tally()
    df = df.get_pandas_dataframe()

zern_orders = np.array(df['zernike'])
zern_mean = df['mean']

def get_max_index(orders,desired):
    for i,order in enumerate(orders):
        order = order.split(',')[0][1:]
        if int(order) <= desired:
            pass
        else:
            break
    return i
max_order = int(zern_orders[-1].split(',')[0][1:])
max_indices = [get_max_index(zern_orders,n) for n in range(max_order+1)]

'''
Plotting
'''
fig,ax = plt.subplots(subplot_kw=dict(projection='polar'))

azimuths = np.radians(np.linspace(0, 360, 50))
zeniths = np.linspace(0, fuel_outter_radius, 100)
r, theta = np.meshgrid(zeniths, azimuths)

def funcanimate(order):
    max_index = max_indices[order-1]
    z_n = openmc.Zernike(zern_mean[:max_index], radius=fuel_outter_radius)
    values = z_n(zeniths,azimuths)
    ax.contourf(theta,r,values, cmap='jet')

animation = FuncAnimation(fig = fig,func = funcanimate, frames = [n for n in range(1,max_order+1)])
animation.save('Zernike_plot.mp4')
    
    






    
    