import numpy as np
from numpy.random import normal
from matplotlib import pyplot as plt
from tqdm import tqdm
import matplotlib.patches as patches

def create_initial_configuration(g, N):
    L = g * N
    X, Y = np.zeros((N, N)), np.zeros((N, N))
    
    for i in range(N):
        for j in range(N):
            x, y = j*g - L/2 + g/2, i*g - L/2 + g/2
            X[i, j], Y[i, j] = x, y

    return X, Y, L

def create_initial_velocities(N, T): 
    vx, vy = [normal(loc=0, scale=np.sqrt(T/2), size=(N, N)) for _ in range(2)]

    return vx, vy

def thermostat(N, T, vx, vy):
    E_k = 0

    for i in tqdm(range(N), desc='Calculating E_k, func:\'thermostat()\''):
        for j in range(N):
            E_k += vx[i][j]**2 + vy[i][j]**2

    return E_k, N**2 * T, np.mean(vx), np.mean(vy)

def visualize_configuration(X, Y, L, order_parameter=None):
    plt.scatter(X, Y)
    rectangle = patches.Rectangle((-L/2, -L/2), L, L, edgecolor='r', facecolor='none')
    plt.gca().add_patch(rectangle)
    plt.axis('equal')
    plt.title(f'The order parameter is \'{order_parameter}\'.')
    plt.show()

def get_dist(L, x1, x2, y1, y2):
    dx, dy = x2 - x1, y2 - y1
    dx -= L * np.round(dx / L)
    dy -= L * np.round(dy / L)

    dist_vect = np.array([dx, dy])
    r = np.linalg.norm(dist_vect)
    norm_vect = [entry/r for entry in dist_vect]

    return dist_vect, r, norm_vect


def Lennard_Jones(r, r_c, k_LJ, norm_vect): #r = float, norm_vect from get_dist to get direction of force
    delta = 1e-5*r_c
    if r >= r_c:
        return [0, 0]
    
    else:
        function_value = k_LJ * ((r_c/r)**12 - (r_c/r)**6)
        gradient_value = (k_LJ * ((r_c/(r+delta))**12 - (r_c/(r+delta))**6) - function_value) / delta
        gradient = [entry*gradient_value for entry in norm_vect]
    
    return gradient

def get_forces(N, L, X, Y):
    k_S, k_LJ, r_c = 1, 1, 3.4
    spring_forces_X, spring_forces_Y, LJ_forces_X, LJ_forces_Y = [np.zeros((N, N)) for _ in range(4)]

    for i in range(N):
        for j in range(N):
            
            dist_vect, r, norm_vect = get_dist(L, X[i, j], X[i, (j+1)%N], Y[i, j], Y[i, j])
            spring_force_X, spring_force_Y = k_S*dist_vect[0], k_S*dist_vect[1]
            spring_forces_X[i, j] += spring_force_X
            spring_forces_X[i, (j+1)%N] -= spring_force_X
            spring_forces_Y[i, j] += spring_force_Y
            spring_forces_Y[i, (j+1)%N] -= spring_force_Y
            
            dist_vect, r, norm_vect = get_dist(L, X[i, j], X[i, j], Y[i, j], Y[(i+1)%N, j])
            spring_force_X, spring_force_Y = k_S*dist_vect[0], k_S*dist_vect[1]
            spring_forces_X[i, j] += spring_force_X
            spring_forces_X[(i+1)%N, j] -= spring_force_X
            spring_forces_Y[i, j] += spring_force_Y
            spring_forces_Y[(i+1)%N, j] -= spring_force_Y

    for h in tqdm(range(N), desc='getting LJ-forces: '):
        for i in range(h+1, N): #ranges are wrong, think...
            x1, y1 = X[h, i], Y[h, i]
            for j in range(h, N):
                for k in range(i+1, N):
                    x2, y2 = X[j, k], Y[j, k]
                    dist_vect, r, norm_vect = get_dist(L, x1, x2, y1, y2) #r = abs distance
                    LJ_force = Lennard_Jones(r, r_c, k_LJ, norm_vect)
                    
                    LJ_forces_X[h, i] += LJ_force[0]
                    LJ_forces_Y[h, i] += LJ_force[1]
                    LJ_forces_X[j, k] -= LJ_force[0]
                    LJ_forces_Y[j, k] -= LJ_force[1]
    
    total_forces_X = [spring_forces_X[i] + LJ_forces_X[i] for i in range(len(spring_forces_X))] 
    total_forces_Y = [spring_forces_Y[i] + LJ_forces_Y[i] for i in range(len(spring_forces_Y))]   

    return LJ_forces_X, LJ_forces_Y           
                    
def main():
    g, N = 0.1, 5
    X, Y, L = create_initial_configuration(g, N)
    total_forces_X, total_forces_Y = get_forces(N, L, X, Y)
    print(total_forces_X)

main()
