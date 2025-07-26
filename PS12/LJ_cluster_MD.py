# import packages for basic math, plotting, linear algebra, etc.
from numpy import *
from numpy.linalg import *
from numpy.random import *
from matplotlib.pyplot import *
from scipy.special import binom, erf, erfc

class histogram():
    def __init__(self,limits,binwidth):
        self.limits = limits
        self.binwidth = binwidth
        self.vals = arange(self.limits[0] + self.binwidth / 2, self.limits[1], self.binwidth)
        self.histo = 0 * self.vals
        self.N_samples = 0

    def add_sample(self,dat):
        self.N_samples += 1
        if dat > self.limits[0] and dat < self.limits[1]:
            bin_index = int((dat - self.limits[0]) / self.binwidth)
            self.histo[bin_index] += 1

    def normalize(self):
        self.histo = self.histo / (self.N_samples * self.binwidth)

    def barplot(self):
        bar(self.vals, self.histo, width=0.95 * self.binwidth, color='k')
        
    def lineplot(self):
        plot(self.vals, self.histo)

def plot_circle(center,radius):
    npoints = 100
    theta = arange(0,2*pi + 1e-7,2*pi/npoints)
    x = center[0] + radius*cos(theta)
    y = center[1] + radius*sin(theta)
    plot(x,y,'k',linewidth=2)

def draw_config():
    clf()
    for i in range(N):
        plot_circle(r[i, :], 0.5)
    axis('equal')
    gca().set_adjustable("box")
    view_scale = 4
    xlim(-view_scale, view_scale)
    ylim(-view_scale, view_scale)
    pause(0.01)

def init_config():
    r = zeros((N, 2))
    n_side = int(sqrt(N) + 0.99)
    count = 0
    for row in range(n_side):
        for column in range(n_side):
            if count < N:
                r[count, :] = [row, column]
                count += 1
    return r

N = 7
delta_t = 0.01
total_time = 100
N_steps = int(total_time/delta_t)
T = 0.05
k_coll = 1

r = init_config()
v = zeros((N,2))
draw_config()

from numba import jit
@jit(nopython=True)
def compute_forces_and_potential(r):
    forces = zeros((N,2))
    potential = 0
    for i in range(N):
        for j in range(i+1,N):
            dr = r[i,:] - r[j,:]
            dr2 = dr @ dr
            force_factor = 48 * ( dr2**(-7) - 0.5 * dr2**(-4) )
            forces[i,:] += force_factor * dr
            forces[j,:] += force_factor * (-dr)

            potential += 4 * ( dr2**(-6) - dr2**(-3) )

    return forces, potential

forces, potential_energy = compute_forces_and_potential(r)

kinetic_traj = zeros(N_steps)
potential_traj = zeros(N_steps)
for step in range(N_steps):
    v = v + 0.5 * delta_t * forces
    r = r + delta_t * v
    forces, potential_energy = compute_forces_and_potential(r)
    v = v + 0.5 * delta_t * forces

    kinetic_energy = 0.5 * sum( v**2 )
    kinetic_traj[step] = kinetic_energy
    potential_traj[step] = potential_energy

    for i in range(N):
        if rand() < k_coll * delta_t:
            speed = sqrt(-2 * T * log(rand()) )
            angle = 2 * pi * rand()
            v[i,:] = speed * array([cos(angle), sin(angle)])

    if step % 1000 == 0:
        draw_config()

clf()
time_traj = arange(N_steps)*delta_t
plot(time_traj,kinetic_traj)
plot(time_traj,potential_traj)
plot(time_traj,kinetic_traj + potential_traj)
