import openmc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
import glob

for file in glob.glob('*.xml'):     #removes all current .xml files just to start from scratch
    os.remove(file)

save_path = r'/Users/nglaser3/OpenMC_FETs/TRISO_Model/'


os.system('echo "Starting Model Construction"')

TRISO = openmc.Model()




''' 
Materials 
'''

SiC = openmc.Material(1,name='SiC')
SiC.add_elements_from_formula('SiC')
SiC.set_density('g/cm3',3.2)

PyC = openmc.Material(2,name = 'PyC')
PyC.add_element('C',1)
PyC.set_density('g/cm3',1.9)

Carbon = openmc.Material(3,name = 'Carbon Buffer')
Carbon.add_element('C',1)
Carbon.set_density('g/cm3',1)

UCO = openmc.Material(4,name='Fuel')
UCO.add_element('U',3,enrichment=19.75,)
UCO.add_element('C',3)
UCO.add_element('O',2)
UCO.set_density('g/cm3',10.75)

TRISO.materials = openmc.Materials([SiC,PyC,Carbon,UCO])




'''
Geometry 
'''
#triso universe generation
outer_radius = 566.7e-4 #outer radius of OPyC
pf = .3
seed = 1
radii = [400e-4, 475e-4, 510e-4, 546.7e-4] #outer radii for fuel, buffer, ipyc, and SiC
spheres = [openmc.Sphere(r = radius) for radius in radii]
cells = [openmc.Cell(fill=UCO, region=-spheres[0]),
        openmc.Cell(fill = Carbon,region = +spheres[0] & -spheres[1]),
        openmc.Cell(fill = PyC, region = +spheres[1]&-spheres[2]),
        openmc.Cell(fill=SiC, region = +spheres[2]&-spheres[3]),
        openmc.Cell(fill = PyC, region = +spheres[3])]
triso_univ = openmc.Universe(cells = cells)

#compact generation and sphere packing
compact_radius,compact_height= outer_radius *15 , 1   #eyeball guess from pictures
compact_surface = openmc.ZCylinder(r = compact_radius)
top,bottom = openmc.ZPlane(z0 = compact_height/2,boundary_type ='reflective'), openmc.ZPlane(z0 = -compact_height/2,boundary_type = 'reflective')
compact_region =  -compact_surface & - top & + bottom
centroids = openmc.model.pack_spheres(radius = outer_radius,region=compact_region, pf = pf, seed = seed)

#bounding cells (hexagon / cylinder)
outterhex = openmc.model.HexagonalPrism(edge_length=2*compact_radius,boundary_type = 'reflective',corner_radius=0)
hexagonal_compact = openmc.Cell(region = -outterhex &+bottom &-top&+compact_surface,fill=SiC)
cylinder_compact = openmc.Cell(region = compact_region)

#lattice generation 
trisos = [openmc.model.TRISO(outer_radius=outer_radius,fill = triso_univ, center=centroid) for centroid in centroids]
shape = (4,4,4)
lowerleft,upperright = cylinder_compact.region.bounding_box
pitch = (upperright-lowerleft)/shape
lattice = openmc.model.create_triso_lattice(trisos = trisos, lower_left=lowerleft, pitch=pitch,shape = shape,background=Carbon)
cylinder_compact.fill = lattice

universe = openmc.Universe(cells=[cylinder_compact,hexagonal_compact])
TRISO.geometry = openmc.Geometry(universe)




'''
Settings 
'''
trigger = openmc.Trigger('std_dev',.0001)
flux_tally_zernike = openmc.Tally()
flux_tally_zernike.scores = ['fission']
zernike_filter = openmc.ZernikeFilter(order=5, x=0.0, y=0.0, r=compact_radius)
flux_tally_zernike.filters = [zernike_filter]
flux_tally_zernike.triggers = [trigger]
TRISO.tallies = openmc.Tallies([flux_tally_zernike])

TRISO.settings.source = openmc.IndependentSource(space=openmc.stats.Box(
    [-compact_radius, -compact_radius, -compact_height], [compact_radius, compact_radius, compact_height],only_fissionable=True))
TRISO.settings.keff_trigger = {'type':'std_dev', 'threshold': .0001}
TRISO.settings.trigger_active = True
TRISO.settings.trigger_batch_interval = 5
TRISO.settings.trigger_max_batches = 1000
TRISO.settings.inactive = 500
TRISO.settings.batches = 1000
TRISO.settings.particles = 500000



TRISO.export_to_xml()
universe.plot(color_by = 'cell',width = (compact_radius*4,compact_radius*4))
plt.savefig('geometry.png')
