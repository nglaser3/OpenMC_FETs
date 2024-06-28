import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import glob
compact_radius = 566.7e-4*15
radius_fuel = 400e-4
os.system('echo "Make sure compact radius in this file matches triso_build.py"')
statepoint = glob.glob('statepoint*')
if len(statepoint) == 1:
    statepoint = statepoint[0]
else:
    raise Exception('Wrong number of statepoint files, should only have one')

#there is only one tally so its id is one, I need to look into a better way to do this though
with openmc.StatePoint(statepoint) as sp:
    df = sp.tallies[1].get_pandas_dataframe()


#heads up the plotting takes a very, very long time, order 100 takes about 10 minutes, order 200 takes significantly longer
#this is because I am passing a pandas dateframe which is super slow, and also FuncAnimation is very slow
'''
Plotting Functions
'''
def zern_order(n,df):
    '''
    n : max order to include in analysis
    df : dataframe of zernike filter
    '''
    zerns = df['zernike']
    for i in range(len(zerns)):
        order = zerns[i]
        order = int(order.split(',')[0][1:])
        if order <= n:
            pass
        else:
            return i
            break
def xy(r,theta):
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return x,y  

def zern_plot(n):
    index = zern_order(n,df)
    z_n = df['mean'][:index]
    zz = openmc.Zernike(z_n, compact_radius)

    values = zz(zeniths, azimuths)
    ax.set_title('Order {} Zernike Polynomial'.format(n))
    
    contour_f = ax.contourf(x, y, values, cmap='jet')
    #contour_f = ax.contourf(theta, r, values, cmap='jet')
    fig.colorbar(contour_f,cax=cax)
def find_centroids(geo_path,fuel_radius):
    geom = openmc.Geometry.from_xml(geo_path)
    cell_dict = geom.get_all_cells()
    centroids = []
    for cell_id in cell_dict.keys():
        cell = cell_dict[cell_id]
        if cell.fill.name == 'Fuel':
            surf_dic = cell.region.get_surfaces()
            surf = [surf_dic[x] for x in surf_dic.keys()][0]
            x,y,z = surf.x0, surf.y0, surf.z0
            if abs(z) < fuel_radius:
                centroids.append([x,y,z])
    cents = np.array(centroids)
    x,y,z = cents[:,0],cents[:,1],cents[:,2]
    return x,y,z


def circle(x,y,z,R,order):
    '''
    x = x0
    y = y0
    z = z0
    R = Sphere radius
    '''
    r = (R**2-z**2)**.5
    theta = np.linspace(0,2*np.pi,order)
    x_ = r*np.cos(theta)+x
    y_ = r*np.sin(theta)+y
    return x_,y_
  


'''
Plotting
'''



# this plotting scheme is a little bit wonky, but I couldn't get the cross geometry to plot on the same plot when in polar 



max_order = df['zernike'].iloc[-1]
max_order = int(max_order.split(',')[0][1:])
azimuths = np.radians(np.linspace(0, 360, 1000)) 
zeniths = np.linspace(0, compact_radius, 500)    
r, theta = np.meshgrid(zeniths, azimuths)

#x,y = xy(r,theta)

fig, ax = plt.subplots()

div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')
x,y,z = find_centroids('geometry.xml',radius_fuel)
circs = [circle(x0,y0,z0,radius_fuel,30) for x0,y0,z0 in zip(x,y,z)]
for circ in circs: 
    ax.plot(circ[0],circ[1],color = 'k',linewidth = .75)
ax.set_aspect('equal')
animation = FuncAnimation(fig, func=zern_plot, frames=np.arange(1,max_order+1), interval=250)
animation.save('TRISO_short_10000_Z_150.mp4')


os.system('echo "Done! Yay!"')
