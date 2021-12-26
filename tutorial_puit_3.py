import matplotlib.pyplot as plt
from solcore import si, material
from solcore.constants import vacuum_permittivity, q
import solcore.poisson_drift_diffusion as PDD
import numpy as np
from solcore.structure import Layer, Structure
from solcore.structure import Junction
import solcore.quantum_mechanics as QM
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource

T = 300
num_energy = 300

wl = np.linspace(350, 1050, 1001) * 1e-9
E = np.linspace(2.5, 3.25, num_energy) * q

# First, we create the materials of the QW
p_GaN = material('GaN')(T=T, Na=8e22)
QWmat =  material('InGaN')(T=T, strained=True)
Bmat  =  material('GaN')(T=T)
n_GaN = material('GaN')(T=T, Nd=1e24)
Bulk  =  material('GaN')(T=T)

# As well as some of the layers
top_layer     = Layer(width=si("500nm"), material=p_GaN)
barrier_layer = Layer(width=5e-9,       material=Bmat)
well_layer    = Layer(width=3e-9,       material=QWmat)
bottom_layer  = Layer(width=si("100nm"), material=n_GaN)

# We define some parameters need to calculate the shape of the excitonic absorption
alpha_params = {
    "well_width": si("3nm"),
    "theta": 0,
    "eps": 12.9 * vacuum_permittivity,
    "espace": E,
    "hwhm": si("6meV"),
    "dimensionality": 0.16,
    "line_shape": "Gauss"}

#  A single QW with interlayers

structure_1 = Structure([barrier_layer, well_layer, barrier_layer], substrate=Bulk)

output_1 = QM.schrodinger(structure_1, quasiconfined=0, graphtype='potentials', num_eigenvalues=20, show=True)

#  Multiple QWs with interlayers
structure_2 = Structure([top_layer, barrier_layer] + 4 * [well_layer, barrier_layer] + [bottom_layer],substrate=Bulk)

output_2 = QM.schrodinger(structure_2, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200,show=True)

# ========= absorption in QW

output = QM.schrodinger(structure_2, quasiconfined=0, num_eigenvalues=20, alpha_params=alpha_params, calculate_absorption=True)

alfa = output[0]['alphaE'](E)
plt.plot(1240 / (E / q), alfa / 100, label='{}'.format('GaN/GaInN'))

plt.xlabel('Wavelength (nm)')
plt.ylabel('$\\alpha$ cm$^{-1}$')
plt.legend(loc='upper right', frameon=False)
plt.tight_layout()

plt.show()




