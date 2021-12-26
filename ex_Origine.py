import matplotlib.pyplot as plt
import numpy as np
from solcore.constants import vacuum_permittivity, q
from solcore import si, material
from solcore.structure import Layer, Structure
import solcore.quantum_mechanics as QM

T=300
# First we create the materials we need
bulk = material("GaN")(T=T, strained=False)
barrier = material("GaN")(T=T, strained=False)

# As well as some of the layers
top_layer = Layer(width=si("10nm"), material=barrier)
barrier_layer = Layer(width=si("4nm"), material=barrier)
bottom_layer = top_layer
#=========
# The absorption coefficients will be calculated at these energies and stored in alfas
num_energy = 300

E = np.linspace(2.5, 3.25, num_energy) * q

# We define some parameters need to calculate the shape of the excitonic absorption
alpha_params = {
    "well_width": si("3nm"),
    "theta": 0,
    "eps": 12.9 * vacuum_permittivity,
    "espace": E,
    "hwhm": si("6meV"),
    "dimensionality": 0.16,
    "line_shape": "Gauss"
}
#=========
# We create the QW material at the given composition

#QW =material("InGaN")(T=T, strained=False)
QW =material("GaInN")(T=T, In=0.6, strained=False)

# And the layer
well_layer = Layer(width=si("3nm"), material=QW)

# The following lines create the QW structure, with different number of QWs and interlayers. Indicating the substrate
# material with the keyword "substrate" is essential in order to calculate correctly the strain.

# A single QW with interlayers
test_structure_1 = Structure([top_layer, well_layer, bottom_layer], substrate=bulk)
output_1 = QM.schrodinger(test_structure_1, quasiconfined=0, graphtype='potentials', num_eigenvalues=20, show=True)

test_structure_2 = Structure([top_layer, barrier_layer] + 10 * [well_layer, barrier_layer] + [bottom_layer],
                             substrate=bulk)
output_2 = QM.schrodinger(test_structure_2, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200,
                          show=True)
# ========= absorption in QW

output = QM.schrodinger(test_structure_2, quasiconfined=0, num_eigenvalues=20, alpha_params=alpha_params, calculate_absorption=True)

alfa = output[0]['alphaE'](E)
plt.plot(1240 / (E / q), alfa / 100, label='{}'.format('GaN/GaInN'))
#int(0.6 * 100)
#plt.xlim(826, 1100)
#plt.ylim(0, 23000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('$\\alpha$ cm$^{-1}$')
plt.legend(loc='upper right', frameon=False)
plt.tight_layout()

plt.show()




