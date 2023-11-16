import numpy as np
import matplotlib as plt
import random as rn
from numpy.linalg import norm
from numpy.random import normal

N = 30
g = 3.5
k = 0.05
rc = 3.4
T = 0.1 #up to 0.6


def create_initial_configuration(g, N):
    L = g * N
    X, Y = np.zeros((N, N)), np.zeros((N, N))  
    for i in range(N):
        for j in range(N):
            x, y = j * g - L/2, i * g - L/2,
            X[i, j], Y[i, j] = x, y

    return X, Y, L


def create_initial_velocities(N,T): 
    #normal uses gaussian distribution and takes rn values, but probabilities for each arent the same, problem?
    vx = normal(loc=0, scale=np.sqrt(T/2), size=(N,N))
    vy = normal(loc=0, scale=np.sqrt(T/2), size=(N,N))
    return vx, vy


def thermostat(N,T,vx,vy):
    vxflat = vx.flatten() #create 1D arrays of vx, vy
    vyflat = vy.flatten()
    ekin1 = (N**2) * T
    ekin2 = np.dot(vxflat,vxflat) + np.dot(vyflat,vyflat) #dot product = squared velocities
    placeholder = ekin1/ekin2 #calculate factor so ekin1 and ekin2 are equal
    vxnew = vx * np.sqrt(placeholder) #times placeholder, so that ekin1 =ekin2 are the same
    vynew = vy * np.sqrt(placeholder)
    vxnew -= np.mean(vxnew)
    vynew -= np.mean(vynew)  
    return vxnew, vynew


def fold(L,x): #not sure if right like this
    x -= L*round(x/L)
    return x


def get_dist(L,x1,y1,x2,y2):
    distance_vector = fold(L,[x1,y1] - [x2,y2])
    distance = norm(distance_vector)
    return distance_vector,distance


def get_forces(N,L,X,Y): #feather force und pot force nedd be equal => equilibrium
    n,m = create_initial_configuration(g, N)
    fx = []
    fy = []
    r = get_dist(L, x1, y1, x2, y2)
    F_lennart = (rc/r)**12 - (rc/r)**(6)
   

    F = -k * r
    return fx, fy

def visualize_configuration(X, Y, order_parameter=None):
    plt.scatter(X, Y)
    plt.axis('equal')
    plt.title(f'The order parameter is \'{order_parameter}\'.')
    plt.show()

print(thermostat(N,T,*create_initial_velocities(N,T)))
#print(create_initial_velocities(N,T))
#grid_nodes = create_initial_configuration(g, N)

# Printing the coordinates of the nodes
#for node in grid_nodes:
    #print(node)
