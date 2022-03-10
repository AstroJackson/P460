import matplotlib.pyplot as plt, matplotlib as mpl, numpy as np
from scipy.integrate import quad as integrate

from matplotlib.animation import FuncAnimation
import matplotlib.patches as mpl_patches

a = 1.
c = 1.
k = .1
pi = np.pi
s = .01
points = 500

def gaussian(x, sigma = 1):
    value_at_0 = 1 /(np.sqrt(2*pi)) * np.exp( -(0-a/2)**2/ (2*sigma**2)) # ensure boundary conditions are satisfied, though this is not required
    return 1 /(np.sqrt(2*pi)) * np.exp( -(x-a/2)**2/ (2*sigma**2)) - value_at_0

def flat_then_peak(x):
    if x <= 0.4 or x >= .6:
        return 0
    elif x <= 0.5:
        return x - .4
    else:
        return -x + .6

def quad(x):
    return x*x

def order3(x):
    return 15 * (x**3-1.42*x**2 + .5*x)

def phi_sin(x, n, *args):
    phi = order3(x)
    return phi * np.sin(n*pi*x/a)

#plt.plot(np.linspace(0, a, points), gaussian(np.linspace(0, a, points), s))

def psi_sin(x, n, *args):
    psi = 0
    return psi * np.sin(n*pi*x/a)

def An(n, *args):
    integration = integrate(phi_sin, 0, a, args = (n, *args))
    return 2. / a * integration[0]

def Bn(n, *args):
    integration = integrate(psi_sin, 0, a, args = (n, *args))
    return 2. / (n*pi*c) * integration[0]

# integration = integrate(phi_sin, 0, a, args = (1, 2))
# print(integration)

#points = 500
x = np.linspace(0, a, points)

max_time = 5
interval = 0.009 # timestep, in seconds
total_frames = int(max_time/interval)
times = np.arange(0, max_time, interval)
u = [0 for i in range(len(times))]

wave_equation = False
dif_equation = True
if wave_equation:
    for i, t in enumerate(times):
        value = np.zeros(points)
        for m in range(1, 50): 
            value += np.sin(m*pi*x/a) * (An(m, s)*np.cos(m*pi*c*t/a) + Bn(m, s)*np.sin(m*pi*c*t/a))

        u[i] = value

if dif_equation:
    for i, t in enumerate(times):
        value = np.zeros(points)
        for m in range(1, 50): 
            value += np.sin(m*pi*x/a) * (An(m, s)*np.exp(-k*m*m*t))

        u[i] = value
#print(u)

u_array = np.array(u)

#plt.plot(x, u[0])


# plt.figure(figsize=(15,8))
# plt.plot(x, u[0])
# plt.xlabel("x")
# plt.ylabel("Height")
# plt.show()

fig, ax = plt.subplots()
fig.set_size_inches(8.5, 8.5)

time_text = ax.text(0.025, np.min(u)+0.01, "")

ln, = ax.plot([], [])

def init():
    time_text.set_text('0')
    ax.set_xlim(0, a)
    ax.set_ylim(np.min(u), np.max(u))
    return time_text, ln

def update(frame):
    time_text.set_text(f"Time = {round(times[frame], 2)}")
    ln.set_data(x, u[frame])
    return time_text, ln

wait_time = interval*1000 *.5
ani = FuncAnimation(fig, update, frames = len(times), init_func=init, blit=True, interval = wait_time, cache_frame_data = True)

ax.set_xlabel("x")
ax.set_ylabel("Height")
ax.set_title("String Over Time")

plt.show()
