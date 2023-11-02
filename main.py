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

def create_initial_configuration(g,N):
    L = g * N #g = distance between nodes
    nodes = []
    for i in range(N):
        for j in range(N):
            x = -L/2 + j * g
            y = -L/2 + i *g
            nodes.append((x, y))

    return nodes, L

def create_initial_velocities(N,T): 
    #normal uses gaussian distribution and takes rn values, but probabilities for each arent the same, problem?
    vx = [normal(loc=0, scale=1, size=N) for element in range(N)]
    vy = [normal(loc=0, scale=1, size=N) for element in range(N)]
    vx -= np.mean(vx)
    vy -= np.mean(vy)
    """
    #to check if sum of all vx/vy elements is 0 (don't want drift)
    sum = 0
    sum2 = 0
    for i in vx:
        sum += np.sum(i)
    for j in vy:
        sum2 += np.sum(j)

    if round(sum, 10) != 0 or round(sum2,10) != 0:
        print("drift in x or y", sum2)
    else:
        print("no drift")"""
    return vx, vy

def thermostat(N,T,vx,vy):
    vxflat = vx.flatten() #create 1D arrays of vx, vy
    vyflat = vy.flatten()
    ekin1 = (N**2) * T
    ekin2 = np.dot(vxflat,vxflat) + np.dot(vyflat,vyflat) #dot product = squared velocities
    placeholder = ekin1/ekin2 #calculate factor so ekin1 and ekin2 are equal
    ekinnew = ekin2 * placeholder 
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

def get_forces(N,L,X,Y):


    return


#print(thermostat(N,T,*create_initial_velocities(N,T)))
#print(create_initial_velocities(N,T))
#grid_nodes = create_initial_configuration(g, N)

# Printing the coordinates of the nodes
#for node in grid_nodes:
    #print(node)
