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

#flo solution
"""
def get_forces(N,L,X,Y):
    N = N-1
    X = np.reshape(X,(N,-1))
    Y = np.reshape(Y,(N,-1))
    M = np.shape(X)[1]
    FX = np.zeros((N,M))
    FY = np.zeros((N,M))
    # storage_LJ = np.zeros((N,M,N,M))
    # storage_spring = np.zeros((N,M,N,M))
    storage_direction = np.zeros((N,M,N,M), dtype=object)
    storage_R = np.zeros((N,M,N,M))
    rc = 3.4
    neighbours = np.array([[1,0],[-1,0],[0,1],[0,-1]])
    if M == 1:
        neighbours = np.array([[1,0],[-1,0]])
    for n in range(N):
        for m in range(M):
            df = np.zeros((2))
            R = np.zeros((N,M))
            direction = np.zeros((N,N), dtype=object)
            for i in range(N):
                for u in range(M):
                    if storage_R[n,m,i,u] == 0:
                        direction[i,u] = get_dist(L,X[n,m],Y[n,m],X[i,u],Y[i,u])[0]
                        R[i,u] = get_dist(L,X[n,m],Y[n,m],X[i,u],Y[i,u])[1]
                        storage_R [n,m,i,u] = R[i,u]
                        storage_direction[n,m,i,u] = direction[i,u]
                        storage_R [i,u,n,m] = R[i,u]
                        storage_direction[i,u,n,m] = -direction[i,u]
                    else:
                        R[i,u] = storage_R[n,m,i,u]
                        direction[i,u] = storage_direction[n,m,i,u]
                #     if R[i,u] <= rc and R[i,u] != 0:
                #         flennard = (24*R[i,u]**6-48)/(R[i,u]**13)
                # # print(direction[i[0],i[1]])
                #         df += flennard * direction[i,u]/np.linalg.norm(direction[i,u])
            for neighbour in neighbours:
                df -= k * direction[(n+neighbour[0])%N,(m+neighbour[1])%N]
                if np.linalg.norm(k * direction[(n+neighbour[0])%N,(m+neighbour[1])%N]) > 1000:
                    print("spring")
                
            close_points = np.argwhere(R <= rc)
            for i in close_points:
                if R[i[0],i[1]] <= 0.0000000001:
                   continue 
                flennard = (24 * R[i[0],i[1]]**6 - 48)/(R[i[0],i[1]]**13)
                # if np.linalg.norm(flennard * direction[i[0],i[1]]/np.linalg.norm(direction[i[0],i[1]])) > 10**15:
                #     print("lennard")
                #     print(R[i[0],i[1]])
                #     plt.scatter([X[n,m],X[i[0],i[1]]],[Y[n,m],Y[i[0],i[1]]], color = ["b","r"],alpha=0.5)
                #     plt.ylim(-L/2-0.1,L/2+0.1)
                #     plt.xlim(-L/2-0.1,L/2+0.1)
                #     plt.show()
                # print(direction[i[0],i[1]])
                df += flennard * direction[i[0],i[1]]/np.linalg.norm(direction[i[0],i[1]])
            if np.linalg.norm(df) > 10**15:
                print(np.linalg.norm(df))
                # df[1] += flennard * direction[i[0],i[1]][1]/np.linalg.norm(direction[i[0],i[1]])
                # print(flennard * direction[i[0],i[1]][0]/np.linalg.norm(direction[i[0],i[1]]),flennard * direction[i[0],i[1]][1]/np.linalg.norm(direction[i[0],i[1]]))
            FX[n,m] = df[0]
            FY[n,m] = df[1]
    return FX, FY
"""



def main():
    N, T = 10_000, 273.15
    vx, vy = create_initial_velocities(N, T)
    print(thermostat(N, T, vx, vy))

    

    