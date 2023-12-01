import numpy as np
from numpy.random import normal
from matplotlib import pyplot as plt
from tqdm import tqdm
import imageio.v2 as imageio
import os
import shutil
from datetime import datetime

'''
Manual:
1.  Install all the necessary libraries
2.  Decide for a directory and make the png and mp4 folders, include these folders in the according constants
3.  In molecular_dynamics(), set up all the preferred parameters
4.  Run the code
'''

png_folder = '/Users/nicol/Documents/Python_Projects/CTL_II/Soft_Solid/png/' # Make sure that these folder exist, consider the differences between mac and PC
mp4_folder = '/Users/nicol/Documents/Python_Projects/CTL_II/Soft_Solid/mp4/'

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
    vx -= np.mean(vx)
    vy -= np.mean(vy)
    return vx, vy

def thermostat(N, T, vx, vy):
    E_k = 0

    for i in tqdm(range(N), desc='Calculating E_k, func:\'thermostat()\''):
        for j in range(N):
            E_k += vx[i][j]**2 + vy[i][j]**2

    return E_k, N**2 * T, np.mean(vx), np.mean(vy)

def get_dist(L, x1, x2, y1, y2):
    dx, dy = x2 - x1, y2 - y1
    dx -= L * np.round(dx / L)
    dy -= L * np.round(dy / L)
    dist_vect = [dy, dx]
    r = np.linalg.norm(dist_vect)
    
    if r == 0:
        dist_vect = norm_vect = [0, 0]
        return dist_vect, r, norm_vect

    norm_vect = [entry/r for entry in dist_vect]

    return dist_vect, r, norm_vect

def get_spring_forces(N, L, X, Y, k_S):
    spring_forces_X, spring_forces_Y = np.zeros((N, N)), np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            dist_vect = get_dist(L, X[i, j], X[i, (j+1)%N], Y[i, j], Y[i, (j+1)%N])[0]
            spring_force_X, spring_force_Y = k_S*dist_vect[1], k_S*dist_vect[0]
            spring_forces_X[i, j] += spring_force_X
            spring_forces_X[i, (j+1)%N] -= spring_force_X
            spring_forces_Y[i, j] += spring_force_Y
            spring_forces_Y[i, (j+1)%N] -= spring_force_Y
            
            dist_vect = get_dist(L, X[i, j], X[(i+1)%N, j], Y[i, j], Y[(i+1)%N, j])[0]
            spring_force_X, spring_force_Y = k_S*dist_vect[1], k_S*dist_vect[0]
            spring_forces_X[i, j] += spring_force_X
            spring_forces_X[(i+1)%N, j] -= spring_force_X
            spring_forces_Y[i, j] += spring_force_Y
            spring_forces_Y[(i+1)%N, j] -= spring_force_Y

    
    return spring_forces_X, spring_forces_Y

def Lennard_Jones(r, r_c, k_LJ, norm_vect, sigma):
    if r >= r_c or r == 0:
        return [0, 0]

    else:
        gradient_value = 24 * k_LJ * ((sigma*r)**-7 - 2*(sigma*r)**-13)
        gradient = [entry*gradient_value for entry in norm_vect]

    return gradient

def get_LJ_forces(N, L, X, Y, k_LJ, r_c, sigma):
    LJ_forces_X, LJ_forces_Y = np.zeros((N, N)), np.zeros((N, N))

    for h in range(N):
        for i in range(N):
            x1, y1 = X[h, i], Y[h, i]
            for j in range(N):
                for k in range(N):
                    x2, y2 = X[j, k], Y[j, k]
                    dist_vect, r, norm_vect = get_dist(L, x1, x2, y1, y2)
                    LJ_vect = Lennard_Jones(r, r_c, k_LJ, norm_vect, sigma)
                    LJ_forces_X[h, i] += LJ_vect[1]
                    LJ_forces_Y[h, i] += LJ_vect[0]

    return LJ_forces_X, LJ_forces_Y

def total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma):
    if k_LJ == 0 or r_c >= L/2:
        return get_spring_forces(N, L, X, Y, k_S)

    elif k_S == 0:
        return get_LJ_forces(N, L, X, Y, k_LJ, r_c, sigma)
    
    else:
        LJ_forces_X, LJ_forces_Y = get_LJ_forces(N, L, X, Y, k_LJ, r_c, sigma)
        spring_forces_X, spring_forces_Y = get_spring_forces(N, L, X, Y, k_S)
        
        forces_x = LJ_forces_X + spring_forces_X
        forces_y = LJ_forces_Y + spring_forces_Y

        return forces_x, forces_y

def single_MD_step(N, L, T, dt, X, Y, vx, vy, forces_x, forces_y, k_S, k_LJ, r_c, sigma):
    vx += forces_x * dt/2
    vy += forces_y * dt/2

    X += vx * dt
    Y += vy * dt
    X = (X+L/2)%L - L/2
    Y = (Y+L/2)%L - L/2

    forces_x, forces_y = total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma)

    vx += forces_x * dt/2
    vy += forces_y * dt/2   

    vx *= np.sqrt(T/2) / np.std(vx)
    vy *= np.sqrt(T/2) / np.std(vy)
    vx -= np.mean(vx)
    vy -= np.mean(vy)

    return X, Y, vx, vy, forces_x, forces_y

def order_parameter(N, L, X, Y, cutoff_distance):
    phi = 0
    order_array = np.zeros((N, N))
    
    for h in range(N):
        for i in range(N):
            x1, y1 = X[h, i], Y[h, i]
            j = 0
            while j < N:
                k = 0
                while k < N:
                    x2, y2 = X[j, k], Y[j, k]
                    r = abs(get_dist(L, x1, x2, y1, y2)[1])
                    if r <= cutoff_distance and (h, i) != (j, k):
                        phi += 1
                        order_array[h, i] = 1
                        j, k = N, N
                    else:
                        k += 1
                j += 1
    
    phi /= N**2

    return phi, order_array

def visualize_configuration(X, Y, g, N, T, dt, k_S, k_LJ, r_c, L, cutoff_distance, step, r_0, show_grid):
    plt.clf()
    plt.cla()
    order_param, order_array = order_parameter(N, L, X, Y, cutoff_distance)
    
    if show_grid: # Don't ask
        for i in range(N):
            for j in range(N):
                if abs(X[i, j] - X[i, (j+1)%N]) <= L/2 and abs(Y[i, j] - Y[i, (j+1)%N]) <= L/2:
                    plt.plot([X[i, j], X[i, (j+1)%N]], [Y[i, j], Y[i, (j+1)%N]], color='grey')
                elif abs(X[i, j] - X[i, (j+1)%N]) > L/2 and abs(Y[i, j] - Y[i, (j+1)%N]) <= L/2:
                    if X[i, j] < X[i, (j+1)%N]:
                        plt.plot([X[i, j] + L, X[i, (j+1)%N]], [Y[i, j], Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N] - L], [Y[i, j], Y[i, (j+1)%N]], color='grey')
                    else:
                        plt.plot([X[i, j], X[i, (j+1)%N] + L], [Y[i, j], Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j] - L, X[i, (j+1)%N]], [Y[i, j], Y[i, (j+1)%N]], color='grey')
                elif abs(X[i, j] - X[i, (j+1)%N]) <= L/2 and abs(Y[i, j] - Y[i, (j+1)%N]) > L/2:
                    if Y[i, j] < Y[i, (j+1)%N]:
                        plt.plot([X[i, j], X[i, (j+1)%N]], [Y[i, j] + L, Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N]], [Y[i, j], Y[i, (j+1)%N] - L], color='grey')
                    else:
                        plt.plot([X[i, j], X[i, (j+1)%N]], [Y[i, j], Y[i, (j+1)%N] + L], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N]], [Y[i, j] - L, Y[i, (j+1)%N]], color='grey')
                else:
                    if X[i, j] < X[i, (j+1)%N] and Y[i, j] < Y[i, (j+1)%N]:
                        plt.plot([X[i, j] + L, X[i, (j+1)%N]], [Y[i, j] + L, Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N] - L], [Y[i, j], Y[i, (j+1)%N] - L], color='grey')
                    elif X[i, j] < X[i, (j+1)%N] and Y[i, j] >= Y[i, (j+1)%N]:
                        plt.plot([X[i, j] + L, X[i, (j+1)%N]], [Y[i, j] - L, Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N] - L], [Y[i, j], Y[i, (j+1)%N] + L], color='grey')
                    elif X[i, j] >= X[i, (j+1)%N] and Y[i, j] < Y[i, (j+1)%N]:
                        plt.plot([X[i, j] - L, X[i, (j+1)%N]], [Y[i, j] + L, Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N] + L], [Y[i, j], Y[i, (j+1)%N] - L], color='grey')
                    elif X[i, j] >= X[i, (j+1)%N] and Y[i, j] >= Y[i, (j+1)%N]:
                        plt.plot([X[i, j] - L, X[i, (j+1)%N]], [Y[i, j] - L, Y[i, (j+1)%N]], color='grey')
                        plt.plot([X[i, j], X[i, (j+1)%N] + L], [Y[i, j], Y[i, (j+1)%N] + L], color='grey')
                
                if abs(X[i, j] - X[(i+1)%N, j]) <= L/2 and abs(Y[i, j] - Y[(i+1)%N, j]) <= L/2:
                    plt.plot([X[i, j], X[(i+1)%N, j]], [Y[i, j], Y[(i+1)%N, j]], color='grey')
                elif abs(X[i, j] - X[(i+1)%N, j]) > L/2 and abs(Y[i, j] - Y[(i+1)%N, j]) <= L/2:
                    if X[i, j] < X[(i+1)%N, j]:
                        plt.plot([X[i, j] + L, X[(i+1)%N, j]], [Y[i, j], Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j] - L], [Y[i, j], Y[(i+1)%N, j]], color='grey')
                    else:
                        plt.plot([X[i, j], X[(i+1)%N, j] + L], [Y[i, j], Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j] - L, X[(i+1)%N, j]], [Y[i, j], Y[(i+1)%N, j]], color='grey')
                elif abs(X[i, j] - X[(i+1)%N, j]) <= L/2 and abs(Y[i, j] - Y[(i+1)%N, j]) > L/2:
                    if Y[i, j] < Y[(i+1)%N, j]:
                        plt.plot([X[i, j], X[(i+1)%N, j]], [Y[i, j] + L, Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j]], [Y[i, j], Y[(i+1)%N, j] - L], color='grey')
                    else:
                        plt.plot([X[i, j], X[(i+1)%N, j]], [Y[i, j], Y[(i+1)%N, j] + L], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j]], [Y[i, j] - L, Y[(i+1)%N, j]], color='grey')
                else:
                    if X[i, j] < X[(i+1)%N, j] and Y[i, j] < Y[(i+1)%N, j]:
                        plt.plot([X[i, j] + L, X[(i+1)%N, j]], [Y[i, j] + L, Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j] - L], [Y[i, j], Y[(i+1)%N, j] - L], color='grey')
                    elif X[i, j] < X[(i+1)%N, j] and Y[i, j] >= Y[(i+1)%N, j]:
                        plt.plot([X[i, j] + L, X[(i+1)%N, j]], [Y[i, j] - L, Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j] - L], [Y[i, j], Y[(i+1)%N, j] + L], color='grey')
                    elif X[i, j] >= X[(i+1)%N, j] and Y[i, j] < Y[(i+1)%N, j]:
                        plt.plot([X[i, j] - L, X[(i+1)%N, j]], [Y[i, j] + L, Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j] + L], [Y[i, j], Y[(i+1)%N, j] - L], color='grey')
                    elif X[i, j] >= X[(i+1)%N, j] and Y[i, j] >= Y[(i+1)%N, j]:
                        plt.plot([X[i, j] - L, X[(i+1)%N, j]], [Y[i, j] - L, Y[(i+1)%N, j]], color='grey')
                        plt.plot([X[i, j], X[(i+1)%N, j] + L], [Y[i, j], Y[(i+1)%N, j] + L], color='grey')                
                
    plt.scatter(X[order_array == 1], Y[order_array == 1], color='red')
    plt.scatter(X[order_array == 0], Y[order_array == 0], color='green')
    
    plt.axis('equal')
    plt.xlim(-L/2, L/2)
    plt.ylim(-L/2, L/2)
    
    plt.gca().set_aspect('equal', adjustable='box')
    plt.title(f'$g={g},\;N={N},\;T={T},\;\Delta t={dt},\;t={step*dt:.2f},$ \n $k_S={k_S},\;k_{{LJ}}={k_LJ},\;r_c={r_c},\;r_0^{{\;LJ}}={r_0:.2f},\;\Phi={order_param:.2f}$')
    plt.xticks([])
    plt.yticks([])

def create_video(fps):
    images = []
    images.clear()

    for filename in tqdm(sorted(os.listdir(png_folder), key=lambda x: int(x.split('_')[1].split('.')[0])), desc='Generating video'):
        if filename.endswith('.png'):
            file_path = os.path.join(png_folder, filename)
            images.append(imageio.imread(file_path))
    
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime("%Y.%m.%d_%H.%M.%S")
    video_path = f'{mp4_folder}output_{formatted_datetime}.mp4'

    if fps > 30:
        frame_skip = round(fps / 30)
        effective_fps = round(fps / frame_skip)
        images = images[::frame_skip]  # Include every n-th frame

        imageio.mimsave(video_path, images, fps=effective_fps)
    else:
        imageio.mimsave(video_path, images, fps=fps)

    print(f"Video created successfully: {video_path}")

def delete_folder_contents():
    folder_path = png_folder
    if os.path.exists(folder_path):
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")
        print(f"Contents of folder '{folder_path}' deleted.")
    else:
        print(f"Folder '{folder_path}' does not exist.")

def molecular_dynamics():
    g = 3.5                 # Initial lattice constant
    N = 10                  # N^2 is the particle count
    T = 0.2                 # Temperature, usually <1
    dt = 0.02               # Time step, usually < 0.05
    k_S = 0.05              # Spring constant
    steps = 10_000          # Number of simulated steps
    k_LJ = 1                # Lennard-Jones-prefactor
    r_c = 3.4               # For all r>r_c, the LJ-potential is set to 0
    cutoff_distance = 1.5   # Particles, who have at least one neighbouring particle closer than cutoff_distance contribute to the order parameter
    sigma = 1               # Influences the equilibrium distance (r_0) of the LJ-potential, r_0 = 2^(1/6)/sigma
    r_0 = 2**(1/6)/sigma    # The gradient of the LJ-potential is 0 at r = r_0
    show_grid = True        # Sets if the grid is shown or not

    # Initialize system:
    X, Y, L = create_initial_configuration(g, N)
    vx, vy = create_initial_velocities(N, T)
    forces_x, forces_y = total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma)

    # Run the simulation:
    for step in tqdm(range(steps), desc='Generating frames'):
        X, Y, vx, vy, forces_x, forces_y = single_MD_step(N, L, T, dt, X, Y, vx, vy, forces_x, forces_y, k_S, k_LJ, r_c, sigma)
        visualize_configuration(X, Y, g, N, T, dt, k_S, k_LJ, r_c, L, cutoff_distance, step, r_0, show_grid)
        plt.savefig(f'{png_folder}plot_{step}.png', dpi=100)

if __name__ == "__main__":
    delete_folder_contents()    # Deletes all the frames from a previous simulation
    molecular_dynamics()        # Runs the simulation
    create_video(fps=600)       # Higher fps speeds up the video
