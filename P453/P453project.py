from turtle import clear
import matplotlib.pyplot as plt, matplotlib as mpl, numpy as np, os
from scipy.integrate import quad as integrate
 
from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpl_patches
from scipy.optimize import root_scalar
from scipy import integrate
 
V0 = -100.
hbar = 1.
a = 1.
m = 1.
def energy_equation(E):
   try:
       return np.tan(np.sqrt(2*m*(E-V0)) * a / (2 * hbar)) - np.sqrt(-E / (E - V0))
   except ZeroDivisionError:
       return 3
E0_info = root_scalar(energy_equation, bracket = (-97, -96) ,x0 = -96)
E0 = E0_info.root
 
q = np.sqrt(-2*m*E0) / hbar
k = np.sqrt(2*m*(E0-V0)) / hbar
# def test(x):
#     return np.exp(q*x/2)
# print(test(-5))
n1 = np.cos(k*a/2)/ (np.exp(-q*a/2)) * integrate.quad(lambda x: (np.exp(q*x))**2, -500*a, -a/2)[0]
n2 = integrate.quad(lambda x: (np.cos(k*x))**2, -a/2, a/2)[0]
norm = 2*n1 + n2
H = np.cos(k*a/2) / (np.exp(-q*a/2))
 
def psi_t0(x):
   if x <= -a/2:
       return H*np.exp(q*x) / norm**0.5
   elif x >= a/2:
       return H*np.exp(-q*x) / norm**0.5
   else:
       return np.cos(k*x) / norm**0.5
 
## Custom potentials:
N = 500
x = np.linspace(-5*a, 5*a, N+1)
dx = x[1] - x[0]
def V1(x):
   if x < -a/2 or x > a/2:
       return 0
   else:
       return V0
 
def V2(x):
   if x < -a or x > a:
       return 0
   else:
       return 2*V0
V_used = V1
used_potential = "V1"
V = [V_used(x) for x in x]
##
 
imag = complex(0,1)
 
dt = 0.1
T = 10
times = np.arange(dt, T, dt)
 
 
psi_initial = np.array([complex(psi_t0(value), 0) for value in x], dtype = "complex")
 
def dpsi_dt_func(t, psi):
   N = len(psi) - 1
   dpsi_dt = np.zeros(N+1, dtype = "complex")
 
   dpsi_dt[0] = 0
   for i in range(1, N):
       dpsi_dt[i] = imag*hbar/(2*m) * (psi[i+1] - 2*psi[i] + psi[i-1]) / dx**2 - imag * (V[i] * psi[i]) / hbar
   dpsi_dt[N] = 0
 
   return dpsi_dt
if f"{T}{used_potential}{N}{dt}{dx}_data.npz" in os.listdir(): # checks if data from this exact simulation is already saved
   data = np.load(f"{T}{used_potential}{N}{dt}{dx}_data.npz")
   time = data["time"]
   psi = data["psi"]
   norm = data["norm"]
else:
   sol = integrate.solve_ivp(dpsi_dt_func, (0, T), psi_initial)
   time = sol.t
   psi = sol.y
   psi_squared = psi * np.conjugate(psi)
   norm = np.array([integrate.cumtrapz(psi_squared[:, i], x) for i in range(len(time))])[:,-1]
   np.savez(f"{T}{used_potential}{N}{dt}{dx}_data", time = time, psi = psi, norm = norm) # saves the simulation so that lengthy simulations need only be done once
#print(psi[:, 0])
psi_real = np.real(psi)
psi_imag = np.imag(psi)
psi_squared = psi * np.conjugate(psi)
 
fig, ax = plt.subplots()
fig.set_size_inches(8.5, 8.5)
 
time_text = ax.text(0.025, np.min(psi_real)+0.01, "")
 
ln_real, = ax.plot([], [], label = "Real")
ln_imag, = ax.plot([], [], label = "Imaginary")
ln_squared, = ax.plot([], [], label = "Mod Squared")
 
def init():
   time_text.set_text('0')
   ax.set_xlim(min(x), max(x))
   ax.set_ylim(np.min([np.min(psi_real), np.min(psi_imag), np.min(psi_squared)]), np.max([np.max(psi_real), np.max(psi_imag), np.max(psi_squared)]))
   return time_text, ln_real, ln_imag, ln_squared
 
###
times = time
###
def update(frame):
   time_text.set_text(f"Time = {round(times[frame], 2)}")
   ln_real.set_data(x, psi_real[:,frame] / norm[frame]**0.5)
   ln_imag.set_data(x, psi_imag[:,frame]/ norm[frame]**0.5)
   ln_squared.set_data(x, psi_squared[:,frame] / norm[frame])
   return time_text, ln_real, ln_imag, ln_squared
 
interval = 1 / dt / dx
wait_time = interval / 250#; print(wait_time)
ani = FuncAnimation(fig, update, frames = len(times), init_func=init, blit=True, interval = wait_time, cache_frame_data = True)
 
ax.set_xlabel("x")
ax.set_ylabel("Psi")
ax.set_title("Particle in 1D Finite Square Well")
#plt.plot(x, V_used(x), legend = "Potential")
ax.legend()
plt.show()

