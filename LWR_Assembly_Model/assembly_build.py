import openmc, numpy, os, matplotlib.pyplot as plt

model = openmc.Model()

fuel_outter_radius = .5
zirc_outter_radius = .025+fuel_outter_radius
topz, bottomz = 100,-100
pitch = .55
num_pins = 8 # 8x8
zern_ord = 100
pinx,piny = 7,0 #0th column, 8th row going from left -right, down-up 

'''
Materials
'''

uo2 = openmc.Material(1,'Fuel')
uo2.add_element('U',1,enrichment= 5.)
uo2.add_element('O',2)
uo2.set_density('g/cm3',10.)

h2o = openmc.Material(2, 'Water')
h2o.add_elements_from_formula('H2O')
h2o.set_density('g/cm3',1.)

zirc = openmc.Material(3,'Zirconium')
zirc.add_element('Zr',1)
zirc.set_density('g/cm3',5.)

model.materials = openmc.Materials([uo2,h2o,zirc])


'''
Geometry
'''
#top and bottom planes
top = openmc.ZPlane(z0 = topz, boundary_type = 'reflective')
bottom = openmc.ZPlane(z0 = bottomz, boundary_type = 'reflective')

inner = openmc.ZCylinder(r=fuel_outter_radius)
outter = openmc.ZCylinder(r=zirc_outter_radius)
fuel = openmc.Cell(region = -inner,fill=uo2)
zirconium = openmc.Cell(region = -outter&+inner,fill = zirc)
water = openmc.Cell(region = +outter,fill = h2o)

pin_universe = openmc.Universe(cells = [fuel,zirconium,water])

lattice = openmc.RectLattice()

lattice.lower_left = (-num_pins*pitch,-num_pins*pitch)
lattice.pitch = (pitch*2,pitch*2)
lattice.universes = [[pin_universe]*num_pins]*num_pins
left,right = openmc.XPlane(x0=-pitch*num_pins,boundary_type = 'reflective'),openmc.XPlane(x0=pitch*num_pins,boundary_type = 'reflective')
back,front = openmc.YPlane(y0=-pitch*num_pins,boundary_type = 'reflective'),openmc.YPlane(y0=pitch*num_pins,boundary_type = 'reflective')
boundingrectangle = openmc.Cell(region = +left &-right &+back &-front &-top &+bottom, fill = lattice)
universe = openmc.Universe(cells = [boundingrectangle])
model.geometry = openmc.Geometry(universe)
ll,ur = universe.bounding_box

_furthest = num_pins*pitch*2 - pitch

_x0,_y0 = .55,.55 
'''
Tallies
'''

flux_tally_zern = openmc.Tally()
flux_tally_zern.scores = ['kappa-fission']
zernike_filter = openmc.ZernikeFilter(order = zern_ord,r = fuel_outter_radius,x=_x0,y=_y0)
flux_tally_zern.filters = [zernike_filter]
model.tallies = openmc.Tallies([flux_tally_zern])
'''
Settings 
'''
model.settings.source = openmc.IndependentSource(space = openmc.stats.Box(lower_left = ll,upper_right = ur, only_fissionable = True))
model.settings.particles = 500000
model.settings.batches = 1250
model.settings.inactive = 500

model.export_to_xml()





