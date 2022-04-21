import matplotlib.pyplot as plt, matplotlib as mpl, numpy as np
from scipy.integrate import quad as integrate

from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpl_patches
from scipy.optimize import root_scalar

V0 = -100
hbar = 1
a = 1
m = 1
def energy_equation(E):
    try:
        return np.tan(np.sqrt(2*m*(E-V0)) * a / (2 * hbar)) - np.sqrt(-E / (E - V0))
    except ZeroDivisionError:
        return 3
E0_info = root_scalar(energy_equation, bracket = (-97, -96) ,x0 = -96)
E0 = E0_info.root
print(E0)