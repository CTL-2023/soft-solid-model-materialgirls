import numpy as np
from numpy.random import normal
from matplotlib import pyplot as plt
from tqdm import tqdm

N = 30
g = 3.5
k = 0.05
r_c = 3.4
T = 0.1 #up to 0.6

def create_initial_configuration(g, N):
    L = g * (N-1)
    X, Y = np.zeros((N, N)), np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            x, y = j*g - L/2, i*g - L/2,
            X[i, j], Y[i, j] = x, y

    return X, Y, L

def create_initial_velocities(N, T): 
    vx = normal(loc=0, scale=np.sqrt(T/2), size=(N, N))
    vy = normal(loc=0, scale=np.sqrt(T/2), size=(N, N))

    return vx, vy

def thermostat(N, T, vx, vy):
    E_k = 0

    for i in tqdm(range(N), desc='Calculating E_k, func:\'thermostat()\''):
        for j in range(N):
            E_k += vx[i][j]**2 + vy[i][j]**2

    return E_k, N**2 * T, np.mean(vx), np.mean(vy)

def visualize_configuration(X, Y, order_parameter=None):
    plt.scatter(X, Y)
    plt.axis('equal')
    plt.title(f'The order parameter is \'{order_parameter}\'.')
    plt.show()

def get_dist(L, x1, x2, y1, y2):
    dx, dy = x2-x1, y2-y1
    dx, dy = dx-L*round(dx/L), dy-L*round(dy/L)
    dist_vect, abs_dist = np.array((dx, dy)), np.sqrt(dx**2+dy**2)

    return dist_vect, abs_dist

def Lennard_Jones(r, r_c, k_LJ):
    if r >= r_c:
        value = 0
    else:
        value = k_LJ * ((r_c/r)**12 - (r_c/r)**6)
        
    

# Fuck dis! Epair(r) = 4( r**(-12) - r**(-6) - r_c**(-12) + r_c**(-6 ))


def get_forces(N, L, X, Y):
    k_S, k_LJ = 1, 1
    spring_forces_X, spring_forces_Y, LJ_forces_X, LJ_forces_Y = np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N)), np.zeros((N, N))

    for i in tqdm(range(N), desc='get_forces: spring'):
        for j in range(N):
            dist_vect, abs_dist = get_dist(L, X[i], X[(i+1)%N], Y[j], Y[j])
            spring_forces_X[i][j] += dist_vect[0]*k_S
            spring_forces_X[(i+1)%N][j] -= dist_vect[0]*k_S

            dist_vect, abs_dist = get_dist(L, X[i], X[i], Y[j], Y[(j+1)%N])
            spring_forces_Y[i][j] += dist_vect[1]*k_S
            spring_forces_Y[i][(j+1)%N] -= dist_vect[1]*k_S

    for i in tqdm(range(N), desc = 'get_forces: lennart jones'):
        for j in range(N):
            dist_vect_x, abs_dist = get_dist(L, X[i], X[(i+1)%N], Y[j], Y[j])
            dist_vect_y, abs_dist = get_dist(L, X[i], X[i], Y[j], Y[(j+1)%N])
            if abs_dist <= r_c:
                ...
              



            
            for p in range(N):
                for q in range(N):
                    dist_vect, abs_dist = get_dist(L, X[i], X[(i+1)%N], Y[j], Y[j])



            
def main():
    N, T = 10_000, 273.15
    vx, vy = create_initial_velocities(N, T)
    print(thermostat(N, T, vx, vy))

    

    