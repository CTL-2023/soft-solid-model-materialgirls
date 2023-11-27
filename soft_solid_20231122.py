import numpy as np
from numpy.random import normal
from matplotlib import pyplot as plt
from tqdm import tqdm
import imageio.v2 as imageio
import os
import shutil
from datetime import datetime

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

def visualize_configuration(X, Y, N, L, order_param):
    plt.clf()
    plt.cla()
    
    # Plot lines connecting points
    for i in range(N):
        for j in range(N-1):
            plt.plot([X[i, j], X[i, j+1]], [Y[i, j], Y[i, j+1]], color='grey')
            plt.plot([X[j, i], X[j+1, i]], [Y[j, i], Y[j+1, i]], color='grey')
    
    # Plot additional lines
    for i in range(N):
        if abs(get_dist(L, X[i, 0]+L, X[i, N-1], Y[i, 0], Y[i, N-1])[1]) < L/2:
            plt.plot([X[i, 0]+L, X[i, N-1]], [Y[i, 0], Y[i, N-1]], color='grey')
        if abs(get_dist(L, X[i, 0], X[i, N-1]-L, Y[i, 0], Y[i, N-1])[1]) < L/2:
            plt.plot([X[i, 0], X[i, N-1]-L], [Y[i, 0], Y[i, N-1]], color='grey')
        if abs(get_dist(L, X[0, i], X[N-1, i], Y[0, i], Y[N-1, i]-L)[1]) < L/2:
            plt.plot([X[0, i], X[N-1, i]], [Y[0, i], Y[N-1, i]-L], color='grey')
        if abs(get_dist(L, X[0, i], X[N-1, i], Y[0, i]+L, Y[N-1, i])[1]) < L/2:
            plt.plot([X[0, i], X[N-1, i]], [Y[0, i]+L, Y[N-1, i]], color='grey')

    # Scatter plot of points
    plt.scatter(X, Y, color='black')
    
    # Set axis properties
    plt.axis('equal')
    plt.xlim(-L/2, L/2)
    plt.ylim(-L/2, L/2)
    
    plt.gca().set_aspect('equal', adjustable='box')
    #plt.title(f'The order parameter is \'{order_param}\'.')
    plt.xticks([])
    plt.yticks([])



def get_dist(L, x1, x2, y1, y2):
    dx, dy = x2 - x1, y2 - y1
    dx -= L * np.round(dx / L)
    dy -= L * np.round(dy / L)
    dist_vect = np.array([dy, dx])
    r = np.linalg.norm(dist_vect)
    
    if r == 0:
        dist_vect = norm_vect = [0, 0]
        return dist_vect, r, norm_vect

    norm_vect = [entry/r for entry in dist_vect]

    return dist_vect, r, norm_vect

def Lennard_Jones(r, r_c, k_LJ, norm_vect): #r = float, norm_vect from get_dist to get direction of force
    r, r_c = float(r), float(r_c)
    delta = 1e-5*r_c
    
    if r >= r_c or r == 0:
        return (0, 0)
    
    else:
        function_value = k_LJ * ((r_c/r)**12 - (r_c/r)**6)
        gradient_value = (k_LJ * ((r_c/(r+delta))**12 - (r_c/(r+delta))**6) - function_value) / delta
        gradient = tuple(round(entry*gradient_value, 15) for entry in norm_vect)

    return gradient

def get_spring_forces(N, L, X, Y, k_S):
    
    spring_forces_X, spring_forces_Y = np.zeros((N, N)), np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            dist_vect, r, norm_vect = get_dist(L, X[i, j], X[i, (j+1)%N], Y[i, j], Y[i, (j+1)%N])
            spring_force_X, spring_force_Y = round(k_S*dist_vect[1], 15), round(k_S*dist_vect[0], 15)
            spring_forces_X[i, j] += spring_force_X
            spring_forces_X[i, (j+1)%N] -= spring_force_X
            spring_forces_Y[i, j] += spring_force_Y
            spring_forces_Y[i, (j+1)%N] -= spring_force_Y
            
            dist_vect, r, norm_vect = get_dist(L, X[i, j], X[(i+1)%N, j], Y[i, j], Y[(i+1)%N, j])
            spring_force_X, spring_force_Y = round(k_S*dist_vect[1], 15), round(k_S*dist_vect[0], 15)
            spring_forces_X[i, j] += spring_force_X
            spring_forces_X[(i+1)%N, j] -= spring_force_X
            spring_forces_Y[i, j] += spring_force_Y
            spring_forces_Y[(i+1)%N, j] -= spring_force_Y

    #print(f'spring forces X:\n{spring_forces_X}\nspringY:\n{spring_forces_Y}\n')
    
    return spring_forces_X, spring_forces_Y

def get_LJ_forces(N, L, X, Y):
    # MUESS NO GMACHT WERDE
    k_LJ, r_c = 1, 0.1
    LJ_forces_X, LJ_forces_Y = np.zeros((N, N)), np.zeros((N, N))

    for h in range(N):
        for i in range(N):
            x1, y1 = X[h, i], Y[h, i]
            for j in range(N):
                for k in range(N):
                    x2, y2 = X[j, k], Y[j, k]
                    dist_vect, r, norm_vect = get_dist(L, x1, x2, y1, y2)
                    LJ_vect = Lennard_Jones(r, r_c, k_LJ, norm_vect)
                    LJ_forces_X[h, i] += LJ_vect[1]
                    LJ_forces_Y[h, i] += LJ_vect[0]

    return LJ_forces_X, LJ_forces_Y

def total_forces(N, L, X, Y, k_S):
    return get_spring_forces(N, L, X, Y, k_S) #remove dis

    LJ_forces_X, LJ_forces_Y = get_LJ_forces(N, L, X, Y)
    spring_forces_X, spring_forces_Y = get_spring_forces(N, L, X, Y)
    forces_x = LJ_forces_X + spring_forces_X
    forces_y = LJ_forces_Y + spring_forces_Y

    return forces_x, forces_y

def single_MD_step(N, L, dt, X, Y, vx, vy, forces_x, forces_y, k_S):
    forces_x, forces_y = total_forces(N, L, X, Y, k_S)
    vx += forces_x * dt
    vy += forces_y * dt
    X += vx * dt
    Y += vy * dt
    X = (X+L/2)%L - L/2
    Y = (Y+L/2)%L - L/2

    return X, Y, vx, vy, forces_x, forces_y

def order_parameter(N, L, X, Y):
    phi = 0
    cutoff_distance = 1.5

    for h in range(N):
        for i in range(N):
            x1, y1 = X[h, i], Y[h, i]
            j, k = 0, 0
            while j < N:
                while k < N:
                    x2, y2 = X[j, k], Y[j, k]
                    if get_dist(L, x1, x2, y1, y2)[1] <= cutoff_distance:
                        j = k = N
                        phi += 1
                    else:
                        j += 1
                        k += 1
    phi /= N**2

    return phi

def molecular_dynamics(N, L, g, k, T, MDsteps, dt):

    return

def create_video():
    images = []
    fps = 30
    
    image_folder = '/Users/nicol/Documents/Python_Projects/CTL_II/Soft_Solid/png/'
    video_folder = '/Users/nicol/Documents/Python_Projects/CTL_II/Soft_Solid/mp4/'

    # Clear the images list
    images.clear()

    # Read all image files in the specified folder, sorted alphabetically
    for filename in sorted(os.listdir(image_folder), key=lambda x: int(x.split('_')[1].split('.')[0])):
        if filename.endswith('.png'):
            file_path = os.path.join(image_folder, filename)
            images.append(imageio.imread(file_path))
    current_datetime = datetime.now()
    formatted_datetime = current_datetime.strftime("%Y.%m.%d_%H:%M:%S")
    # Create the video
    video_path = f'{video_folder}output_{formatted_datetime}.mp4'
    imageio.mimsave(video_path, images, fps=fps)

    print(f"Video created successfully: {video_path}")

def delete_folder_contents():
    folder_path = '/Users/nicol/Documents/Python_Projects/CTL_II/Soft_Solid/png'
    # Check if the folder exists
    if os.path.exists(folder_path):
        # Iterate over the folder contents and remove them
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

def frequency_fft(t, dt, y):
    plt.clf(), plt.cla()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    ax1.plot(t, y)
    
    fft_result = np.fft.fft(y)
    
    frequencies = np.fft.fftfreq(len(t), dt)
    frequencies = [float(frequencies[i]) for i in range(len(frequencies))]
    fft_result = [float(fft_result[i]) for i in range(len(fft_result))]

    index = frequencies.index(next(x for x in frequencies if x < 0))
    frequencies = frequencies[:index]
    fft_result = fft_result[:index]
    mode = frequencies[fft_result.index(max(fft_result))]

    ax2.plot(frequencies, np.abs(fft_result))
    ax2.set_title(f'mode: {mode:.2f}')

def main():
    delete_folder_contents()
    g, N, T, dt, k_S, dampening, steps = 0.1, 4, 20, 0.003, 1e4, 0.00, 1_000
    X, Y, L = create_initial_configuration(g, N)
    LJ_forces_X, LJ_forces_Y = get_LJ_forces(N, L, X, Y)
    print(LJ_forces_X)
    print()
    print(LJ_forces_Y)

if __name__ == "__main__":
    main()