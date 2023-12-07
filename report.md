# Group Report: Soft-Solid Particle Simulation Project

In this project, a simulation codebase for a 2D elastic Lennard-Jones model has been implemented using python. Most importantly, the code is able to generate videos of the particle systems and also plot the progression of the ordering parameter $\Phi$, which is defined as the fraction of particles in the system that are part of a so-called high-density phase (HDP). The particle forces, velocities and positions are calculated using the velocity-Verlet algorithm, a classical Newtonian approach. The initial position of the paricle is a square lattice, where each particle is connected to their 4 neighbours via elastic springs. Furthermore, each particle interacts with every other particle through Lennard-Jones forces. The square-shaped system follows periodic boundary conditions. As a result, the particles behave as if they were moving on the surface of a thorus. The initial conditions are highly tweakable, some of the parameters include:

- The initial lattice constant $g$
- The square root of the particle count $N$
- The temperature T
- The time interval of each simulation step $\Delta t$, resp. $dt$
- The spring constant $k_S$
- The number of total steps, $steps$
- A Lennard-Jones constant $k_{LJ}$
- The critical radius $r_c$, after which the LJ-potential is considered $0$
- The $cutoff\;distance$, that determines which particles are part of a HDP
- $\sigma$, that determines the equilibrium distance $r_0$ of the LJ-potential

Here are some formulas:
    $$F_s(r) = k_S\cdot r$$
    $$F_{LJ}(r) = 24\cdot k_{LJ}\cdot [ (\sigma r)^{-7}-2(\sigma r)^{-13}]$$
    $$r_0=\frac{\sqrt[6]{2}}{\sigma}$$
    $$F_{total} = F_S + F_{LJ}$$

The final version of the project can be found [here](https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/final_draft.py). A more thorough description of the task is given in the [readme](https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/README.md).

## Layout of the Algorithm and Functions
### Create Initial Configuration
Generates the initial $N\times N$-arrays of the particle positions of side length $L=N\cdot g$. The coordinates $[0,\; 0]$ are defined as the center of the system.

    create_initial_configuration(g, N):
        Args:
            g: Initial lattice constant
            N: Square root of the particle count
        Returns:
            numpy.ndarray: NxN-array of X-positions of particles
            numpy.ndarray: Y-positions of particles
            float: System side length L
### Create Initial Velocities
Creates two arrays of randomly distributed particle X-, and Y-speeds. The speeds are Gauss-distributed with $\mu_{x,y}=0$ and $\sigma_{x,y}=\sqrt{T/2}$, meaning that the total kinetic energy will be set by the temperature: $E_{tot} = \sum_{N^2}(v_x^2+v_y^2)=N^2\cdot T$.

    create_initial_velocities(N, T):
        Args:
            T: System temperature
        Returns:
            numpy.ndarray: NxN-array of X-velocities of particles
            numpy.ndarray: Y-velocities of particles
### Thermostat
This is a testing function that was important during the development of the code but has no application in the final product.

    thermostat(N, T, vx, vy):
        Args:
            vx: X-velocities
            vy: Y-velocities
        Returns:
            float: Kinetic energy of the system, calculated using vx/vy
            float: Kinetic energy, calculated using T
            float: Mean X-velocity
            float: Mean Y-velocity
### Get Distance
Takes the coordinates of two particles as an input and returns the distance vector between particles, the distance scalar and the normalized vector. It considers periodic boundary conditions:

$$r_x = 
\begin{cases}
    x_2-x_1 & \text{if } x_2-x_1 \leq L/2, \\
    L-(x_2-x_1) & \text{if } x_2-x_1>L/2.
\end{cases}$$

for the x-component of the distance vector.

    get_dist(L, x1, x2, y1, y2):
        Args:
            L: System side length
            x1: X-coordinate of particle 1
            x2: X-coord of part. 2
            y1: ...
            y2: ...
        Returns:
            list: Distance vector between particles
            float: Distance
            list: Normalized vector
### Get Spring Forces
Iterates over each particle and calculates the spring force to the particle next to it. It then adds this force to particle $n$ and substracts this force from particle $n-1$ (actio = reactio).

    get_spring_forces(N, L, X, Y, k_S):
        Args:
            X: NxN array of X-positions
            Y: ...
            k_S: spring constant
        Returns:
            numpy.ndarray: NxN array of spring X-forces
            numpy.ndarray: Y-forces
        Calls:
            get_dist(L, x1, x2, y1, y2)
### Lennard-Jones
Includes the definition of the Lennard-Jones-force:
$$F_{LJ}(r)=\frac{\partial}{\partial r}E_{LJ}(r)=24\cdot k_{LJ} [(\sigma r)^{-7}-2\cdot(\sigma r)^{-13}]$$
Returns the vector of the corresponding gradient.

    Lennard_Jones(r, r_c, k_LJ, norm_vect, sigma):
        Args: 
            r: Distance between two particles
            r_c: Critical radius
            k_LJ: Lennard_Jones constant
            norm_vect: Normalized vector between the two particles
            sigma: Determines eq. distance of the LJ-potential
        Returns:
            list: Directed LJ-force

### Get Lennard-Jones-Forces
Iterates over every particle and calculates the LJ-force to every other particle, then adds the result to the first particle.

    get_LJ_forces(N, L, X, Y, k_LJ, r_c, sigma):
        Returns:
            numpy.ndarray: NxN array of LJ X-forces
            numpy.ndarray: Y-forces
        Calls:
            get_dist(L, x1, x2, y1, y2)
            Lennard_Jones(r, r_c, k_LJ, norm_vect, sigma)
### Get Lennard-Jones-Forces Fast
Iterates over every particle and calculates the Lennard-Jones forces to all the following particles. Then adds this force to the first particle $i$ and substracts it from the second particle $i+j>i$ (actio = reactio). In theory, this should speed up the function by a factor of two for large $N$.

    get_LJ_forces_fast(N, L, X, Y, k_LJ, r_c, sigma):
        ...
### Get Total Forces
Adds the LJ-forces and the spring forces together.

    total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma, fast):
        Args:
            fast: determines which LJ-function to use
        Returns:
            numpy.ndarray: NxN array of total X-forces
            numpy.ndarray: Y-forces
        Calls:
            get_spring_forces(N, L, X, Y, k_S)
            get_LJ_forces_fast(N, L, X, Y, k_LJ, r_c, sigma)
            get_LJ_forces(N, L, X, Y, k_LJ, r_c, sigma)
### Single MD Step
Does one single timestep using the Verlet-algorithm, then makes sure that $\sigma_v$ stays 0 and $\mu_v$ stays $\sqrt{T/2}$. Then returns the updated position, velocity and force arrays.

    single_MD_step(N, L, T, dt, X, Y, vx, vy, forces_x, forces_y, k_S, k_LJ, r_c, sigma, fast)
        Args:
            dt: Length of single time step
            forces_x: NxN array of total X-forces
            forces_y: Y-forces
        Returns:
            numpy.ndarray: NxN array of particle X-positions
            numpy.ndarray: Y-positions
            numpy.ndarray: NxN array of X-velocities
            numpy.ndarray: Y-velocities
            numpy.ndarray: NxN array of X-forces
            numpy.ndarray: Y-forces
        Calls:
            total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma, fast)
### Order Parameter
Calculates the order parameter using the position arrays and outputs this parameter and also an array with entries $1$ or $0$, which states which particles are part of a HDP.

    order_parameter(N, L, X, Y, cutoff_distance):
        Args:
            cutoff_distance: Distance that defines a HDP
        Returns:
            float: Order parameter
            numpy.ndarray: NxN array that states which particles belong to a HDP
        Calls:
            get_dist(L, x1, x2, y1, y2)
### Visualize Configuration
Plots the system and the used parameters with or without grid lines and saves the plot to the *frames folder*. Particles outside of a HDP are plotted green, inside they are plotted red.

    visualize_configuration(X, Y, g, N, T, dt, k_S, k_LJ, r_c, L, step, r_0, show_grid, phi, order_array):
        Args:
            step: Number of the step that is being plotted
            r_0: Equilibrium distance of the LJ-potential
            show_grid: Sets whether or not the lattice grid is shown
            phi: Order parameter
            order_array: Array that states which particles belong to a HDP
### Plot $\Phi$
Plots $\Phi$ against $t$ using two provided lists. Also includes the other provided parameters in the plot. Saves the plot to the *plot folder*.

    plot_phi(t_list, phi_list, g, N, T, dt, k_S, k_LJ, r_c, r_0):
        Args:
            t_list: Sets X-axis
            phi_list: Phi values, Y-axis
Takes all the frames from the frames folder and zips them into a video, which is saved as an mp4 to the *mp4 folder*.

    create_video(fps):
        Args:
            fps: Sets how fast the video is played back (not actual frame rate of the resulting file)
### Delete Folder Contents
Deletes all the frames from the *frame folder*.

    delete_folder_contents():
        ...
### Calculate $\Phi$ vs. $T$
Runs multiple simulations at several different temperatures. Optionally, per temperature step several simulations can be run and the mean value and the standard deviation is then plotted. Finally, one curve is plotted per temperature step and all curves are merged into one plot, which is saved in the *plots folder*.

    phi_vs_T():
        Calls:
            create_initial_configuration(g, N)
            create_initial_velocities(N, T)
            total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma, fast)
            single_MD_step(N, L, T, dt, X, Y, vx, vy, forces_x, forces_y, k_S, k_LJ, r_c, sigma, fast)
            order_parameter(N, L, X, Y, cutoff_distance)
### Molecular Dynamics
Initializes the system using the provided parameters and runs the simulation. After each time step, a frame is saved in the *frames folder*. Finally, the plot_phi() function is called and a plot of the time progression of the order parameter is saved to the *plots folder*.

    molecular_dynamics():
        Calls:
            create_initial_configuration(g, N)
            create_initial_velocities(N, T)
            total_forces(N, L, X, Y, k_S, k_LJ, r_c, sigma, fast)
            single_MD_step(N, L, T, dt, X, Y, vx, vy, forces_x, forces_y, k_S, k_LJ, r_c, sigma, fast)
            order_parameter(N, L, X, Y, cutoff_distance)
            visualize_configuration(X, Y, g, N, T, dt, k_S, k_LJ, r_c, L, step, r_0, show_grid, phi, order_array)
            plot_phi(t_list, phi_list, g, N, T, dt, k_S, k_LJ, r_c, r_0)






## How to Run the Code

1.  Make sure, all the neccessary libraries are installed and that you run python3.  The code was successfully run on an M1 MacBook, using MacOS 14.0 and Python 3.11.5. The used non-standard libraries include:
    - numpy         (simplifies array manipulation)
    - matplotlib    (to generate plots)
    - tqdm          (loading bars for time-intensive simulations)
    - imageio.v2    (uses ffmpeg to generate videos)
2.  Create the three folders and reference them in the according constants *'frames_folder'*, *'plots_folder'*, and *'mp4_folder'*.
3.  Decide for a function you want to run, either *molecular_dynamics()* or *phi_vs_T()*. Enter all the parameters of interest. For example, for the first function, the parameters could be chosen as follows:

        g = 3.5                 # Initial lattice constant
        N = 20                  # N^2 is the particle count
        T = 0.2                 # Temperature, usually < 1
        dt = 0.01               # Time step, usually < 0.05
        k_S = 0.05              # Spring constant
        steps = 25_000          # Number of simulated steps
        k_LJ = 1                # Lennard-Jones-prefactor
        r_c = 3.4               # For all r>r_c, the LJ-potential is set to 0
        cutoff_distance = 1.5   # Particles, who have at least one neighbouring particle closer than cutoff_distance contribute to the order parameter
        sigma = 1               # Influences the equilibrium distance (r_0) of the LJ-potential, r_0 = 2^(1/6)/sigma
        r_0 = 2**(1/6)/sigma    # The gradient of the LJ-potential is 0 at r = r_0
        show_grid = True        # Sets if the grid is shown or not
        fast = True             # Sets which get_LJ_forces() function is being used (only makes a small difference for small N). 
4. Call this function and/or any other functions at the very end of the code. Example:

        if __name__ == "__main__":
            delete_folder_contents()     # Deletes all the frames from a previous simulation
            molecular_dynamics()         # Runs the simulation
            create_video(fps=1_200)      # Higher fps speeds up the video
5.  Run the code.

In this example code, 25.000 time steps of a $20\times20$ particle system are simulated, which leads to the functions single_MD_step(), order_parameter(), and visualize_configuration() being called 25.000 times. Finally, the frames are wrapped into a 30fps-mp4, with a simulated 1.200 fps, resulting in a approximate video length of 20 seconds.

## Problems & Limitations

The execution time of the two functions *get_LJ_forces_fast()* and *order_parameter()* scale with the fourth power of *N*. Some different approaches to optimize this would be:
-   Combining the *order_parameter()*-function with the *get_LJ_forces_fast()*-function, so that the distance between each particle to any other particle has to be only calculated once. This could roughly half the computation time for larger *N*, while still maintaining a time complexity of $N^4$.
- Lowering the dpi of the plot would lower the computation time by a constant value.
- The system could be subdivided into smaller square bins, with side lengths $r_c$, and $cutoff\_ value$ respectively. The *get_LJ_forces_fast()*-function could then compare the particle of interest with only the particles that are in the surrounding 4 bins with side length $r_c$ if we consider actio = reactio. The *order_parameter()*-function could do the same with the second binning system. This would optimize code efficiency to a computation time that scales approximately linearly with the particle count, $N^2$ (Suggested by MK).

Some aspects of the simulation seem unphysical, for example in the *single_MD_step()* function, the velocity of the particles is held constant, no matter how much energy is stored as potential energy in the Lennard-Jones-potetntial. This may imply, that the total energy of the system changes over time. Also, the periodic boundary conditions lead to a thoroidal topology of the system, which probably has some implications on the behavior of the system, especially for smaller *N*.

## Results
The example code in chapter *How to Run the Code* provided the following results:

<table>
  <tr>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/demo_17.png" alt="Image 1"></td>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/demo_17-1.png" alt="Image 2"></td>
  </tr>
</table>
<table>
  <tr>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/demo_17-2.png" alt="Image 1"></td>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/demo_17-3.png" alt="Image 2"></td>
  </tr>
</table>
<table>
  <tr>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/demo_17-4.png" alt="Image 1"></td>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/demo_17-5.png" alt="Image 2"></td>
  </tr>
</table>

The whole simulation took approximately 7 hours using the previously mentioned setup. A video of the same simulation can be found [here](https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_videos/demo_17-1.mp4). Here you find more example [videos](https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_videos/) and [images/plots](https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/).

The order parameter has been tracked for systems of different temperatures, using *phi_vs_T()*. It shows, that a higher $T$ has a significant effect on how fast the first HDPs form, but only a negligible effect on $\Phi_{max}$. This could be due to some errors in the code, since smaller $\Phi_{max}$ for increasing $T$ were expected.
<table>
  <tr>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/phi_vs_T_1.png" alt="Image 1"></td>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/phi_vs_T_2.png" alt="Image 2"></td>
  </tr>
</table>

The execution time has been measured for three different cases: 
-   Repeatedly calling *single_MD_step()*
-   Repeatedly calling *single_MD_step()* & *order_parameter()*
-   Repeatedly calling *single_MD_step()*, *order_parameter()* & *visualize_configuration()*

<table>
  <tr>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/runtime_N.png" alt="Image 1"></td>
    <td><img src="https://github.com/CTL-2023/soft-solid-model-materialgirls/blob/main/Demo_images/runtime_N2.png" alt="Image 2"></td>
  </tr>
</table>

As the code is right now, some of the simulations suggested in the README ($30\times 30$ particles, 100.000 time steps) would take more than four days to compute, which was not feasible.
