from solcore import material
import solcore.poisson_drift_diffusion as PDD
import numpy as np
from solcore.structure import Layer, Structure
from solcore.structure import Junction
import solcore.quantum_mechanics as QM
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
T = 300


wl = np.linspace(350, 1050, 1001) * 1e-9

# First, we create the materials of the QW
QWmat = material('InGaN')(T=T, strained=True)
Bmat =  material('GaN')(T=T)
Bulk = material('GaN')(T=T)

# The final device will have 4 of these QWs.
QW = PDD.QWunit([Layer(width=5e-9, material=Bmat,  role="barrier"),
                 Layer(width=3e-9, material=QWmat, role="well"),
                 Layer(width=5e-9, material=Bmat,  role="barrier")],
                 T=T, repeat=4, substrate=Bulk)

# We solve the quantum properties of the QW, leaving the default values of all parameters
QW_list = QW.GetEffectiveQW(wavelengths=wl)

# Materials for the junction
n_GaN = material('GaN')(T=T, Nd=1e24)

p_GaN = material('GaN')(T=T, Na=8e22)


MyJunction = Junction([Layer(width=10e-9, material=n_GaN, role= "Emitter")] +
                          QW_list +
                         [Layer(width=50e-9, material=p_GaN, role="Base")],
                          sn=1e6, sp=1e6, T=T, kind='PDD')

#=======================

MgF2 = material('MgF2')()
ZnS = material('ZnScub')()

#============================

my_solar_cell = SolarCell([Layer(width = 100e-9, material=MgF2, role="ARC1"),
                           Layer(width = 60e-9,  material=ZnS,  role="ARC2"),
                           MyJunction],T=T, substrate=Bulk)

#====================================

light_source = LightSource(source_type='standard', version='AM1.5g', x=wl,output_units='photon_flux_per_m', concentration=1)

solar_cell_solver(my_solar_cell, 'qe',user_options={'light_source': light_source, 'wavelength': wl, 'optics_method': 'TMM'})

# Calculating the IV characteristics

con = np.logspace(0, 3, 19)
vint = np.linspace(-3.5, 4, 600)
V = np.linspace(-3.5, 0, 300)

allI = []
isc = []
voc = []
FF = []
pmpp = []

for c in con:
    light_source.options['concentration'] = c

    solar_cell_solver(my_solar_cell, 'iv',user_options={'light_source': light_source, 'wavelength': wl,
                                    'optics_method': None,'light_iv': True, 'mpp': True, 'voltages': V,
                                    'internal_voltages': vint})
    isc.append(my_solar_cell.iv['Isc'])
    voc.append(my_solar_cell.iv['Voc'])
    FF.append(my_solar_cell.iv['FF'])
    pmpp.append(my_solar_cell.iv['Pmpp'])
    allI.append(my_solar_cell.iv['IV'][1])

