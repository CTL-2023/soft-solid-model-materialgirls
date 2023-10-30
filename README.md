General information about the CTL course available at https://ctl.polyphys.mat.ethz.ch/

# :wave: PROJECT SoftSolid

<img src="https://www.complexfluids.ethz.ch/images/PROJECT-soft-solid.PNG" width="50%">

In this project you implement a known 2D model for a soft solid, the so-called elastic Lennard-Jones model [1], and explore its features. For movies see [2,3], snapshots are given below. The model consists of a number of particles, whose positions are initially coinciding with the node positions of a regular *N* $\times$ *N* grid with grid constant *g*. The size of the simulation box is therefore *L* $\times$ *L* with *L=Ng*, and periodic boundary conditions apply. The center of the simulation box is at the origin, i.e., at [0,0]. Each of the particles is *permanently* spring-connected with spring coefficient *k* to those four particles that are nearest neighbors in the *initial* configuration. If **r** is the relative vector between two spring-connected particles, the potential energy of the connecting harmonic spring is 

*E<sub>spring</sub>(r)=* $\frac{1}{2}$*kr<sup>2</sup>* 

where *r* is the extension of the spring, i.e., *r=|**r**|*. Moreover, each particle interacts with all those particles whose centers are located within a circle of radius *r<sub>c</sub>* around its own position (note that *L>2r<sub>c</sub>* must hold). This interaction is captured by the radially symmetric Lennard-Jones potential

*E<sub>pair</sub>(r) = 4( r<sup>-12</sup> - r<sup>-6</sup> - r<sub>c</sub><sup>-12</sup> + r<sub>c</sub><sup>-6</sup> )*

which reaches zero at *r=r<sub>c</sub>* by construction. The contribution to the total Lennard-Jones energy from a single pair of particles separated by distance vector **r** is *E<sub>pair</sub>(r)*. As obvious from these expressions for energies, all quantities are considered dimensionless, while dimensions may also be introduced. 

Qualitatively, particles experience a long-range repulsion caused by the elastic springs to the original neighbors, and they attract each other due to the Lennard-Jones potential at small separation distances. This competition gives rise to a phase transition, if parameters are suitably chosen. The four model parameters are *g* (grid constant), *k* (spring coefficient), *T* (temperature), and integer-valued *N=L/g*. Typical parameters to be used in this project are: *N*=30, *g=3.5*, *k=0.05*, cutoff distance *r<sub>c</sub>=3.4* and temperatues *T* in the range from 0.1 to 0.6.

The model just described is solved applying Newton's equations using a certain integration time step &Delta;*t*, supplemented by a temperature control that ensures that the mean squared velocity of particles equals *T/2* (the particle mass is *m=1*). Initially, coordinates **x** are placed on a regular square grid inside the simulation box, velocities **v** are chosen randomly while ensuring $\langle$*v<sup>2</sup>*$\rangle$*=T/2*, initial forces **F** on all particles are calculated, time is set to *t=0*, and the time step is set to &Delta;*t=0.01*. The new coordinates **x** of particles at time *t+&Delta;t* are given in terms of the existing coordinates **x** and existing velocities **v** via the so-called temperature-controlled velocity Verlet algorithm, which consists of six successive operations

1. **v** +=  **F** &Delta;t/2, i.e., change velocities (first part)
2. **x** +=  **v** &Delta;t, i.e., change coordinates
3. **x** -= *L* round(**x**/*L*), i.e., fold particle coordinates back to the central simulation box
4. calculate new forces **F** using the new coordinates **x**
5. **v** +=  **F** &Delta;t/2, i.e., change velocities (second part)
6. multiply all **v** with a common factor so that $\langle$*v<sup>2</sup>*$\rangle$*=T/2*; this is the temperature control

The above 6 commands constitute a single time step during which time *t* increased by &Delta;*t*.

Implementation of this algorithm requires knowledge of all forces **F** on all particles for given coordinates **x**, which is the main part of the algorithm. Using the above expressions for energies, one has to first write down the corresponding forces analytically via force = negative gradient of the potential energy 

     you should test your force routine using the following: If you have only two particles in the system, 
     the first particle at position [-0.1,0.0] and the second particle at position [0.8,0.2], then the Lennard-Jones force 
     on the first particle due to the presence of the 2nd particle within the cutoff distance should be [-93.3782,20.7507]. 
     Actio = reactio applies, so that the sum over all forces on all particles vanishes. If the two particle are moreover 
     spring-connected, the additional spring force on particle 1 is spring coeffient times [0.9,0.2].
     
Running the molecular dynamics for some time generates particle trajectories (coordinates + velocities) for all particles at equidistantly spaced times *t*, while all particles effectively stay within the central simulation box. 

The coordinates contain all structural information about the system. As an example, you will calculate the numbers of particles belonging to the so-called high-density phase (HDP). Each particle that has at least one other particle at distance *1.5* or smaller belongs to the HDP. The fraction of particles belonging to the HDP defines the order parameter &Phi;.
Other order parameters are for example the fraction of particles that belong to the largest cluster of HDP particles. For a completely random configuration of particles inside the simulation box, one has &Phi;*=1-*exp*(-&pi;r<sub>c</sub><sup>2</sup>/g<sup>2</sup>)*. 
     
     This formula can be used to test your routine that calculates the order parameter.

The following table provides some snapshots at different times, for later comparison with your own code. Particles belonging to the HDP phase are red, others green, and bonds (except those that cross boundaries) are represented by black lines. Click on image to zoom in. 

| N| *g* | *k* | *T* | *t=0* | *t=20* | *t=100* | *t=1000* | 
|---|---|---|---|---|---|---|---|
| 30 | 3.5 | 0 | 0.3 | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0-T = 0.3-dt=0.01-Phi=0-t=0.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0-T = 0.3-dt=0.01-Phi=0.62-t=20.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0-T = 0.3-dt=0.01-Phi=0.9-t=99.99.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0-T = 0.3-dt=0.01-Phi=0.97-t=999.99.png" width="100%"> |
| 30 | 3.5 | 0.01 | 0.2 | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0-T = 0.3-dt=0.01-Phi=0-t=0.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.01-T = 0.2-dt=0.01-Phi=0.7-t=20.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.01-T = 0.2-dt=0.01-Phi=0.94-t=99.99.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.01-T = 0.2-dt=0.01-Phi=0.98-t=999.99.png" width="100%"> | 
| 30 | 3.5 | 0.05 | 0.1 | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.1-dt=0.01-Phi=0-t=0.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.1-dt=0.01-Phi=0.36-t=20.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.1-dt=0.01-Phi=0.66-t=99.99.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.1-dt=0.01-Phi=0.82-t=999.99.png" width="100%">
| 30 | 3.5 | 0.05 | 0.2 | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.2-dt=0.01-Phi=0-t=0.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.2-dt=0.01-Phi=0.43-t=20.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.2-dt=0.01-Phi=0.71-t=99.99.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.2-dt=0.01-Phi=0.87-t=999.99.png" width="100%"> |
| 30 | 3.5 | 0.05 | 0.5 | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.5-dt=0.01-Phi=0-t=0.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.5-dt=0.01-Phi=0.48-t=20.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.5-dt=0.01-Phi=0.6-t=99.99.png" width="100%"> | <img src="https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/fig-MD-N=30-g=3.5-k=0.05-T = 0.5-dt=0.01-Phi=0.87-t=999.99.png" width="100%"> |

# Code

The following tasks are meant to become python-functions, with some input like the set of *g,k,T,N* and some output like *X,Y,L". The reason for creating such functions becomes most evident lateron in the molecular_dynamics function, where you may start to read. 

## def create_initial_configuration(*g,k,T,N*)

1. Place *N* $\times$ *N* nodes on a regular two-dimensional grid with lattice spacing *g* in a square box of size *L* $\times$ *L* where *L = gN*. The box is centered at the origin, i.e., the coordinates are all in [-*L*/2,*L*/2]. 
2. Collect the <em>x</em>- and <em>y</em>-components of the node positions in *N* $\times$ *N* arrays *X* and *Y*. The component *X*[*n*,*m*] then contains the *x*-component of the position of the node (*n,m*). In python it is most convenient to use *n,m* $\in$ {0,..,*N*-1}.
3. Return *X* and *Y* and *L*

## def create_initial_velocities(*N,T*)

This function creates and returns two *N* $\times$ *N* random velocity arrays (VX and VY). The element (*n,m*) of VX carries the *x*-component of the velocity of the node (*n,m*). Make sure that the mean(VX)=0 and mean(VY)=0, so that there is no drift of your systems center of mass. 

## def thermostat(*N,T*,VX,VY)

This function 
1. multiplies both VX and VY with a unique factor so that the kinetic energy is *N* $\times$ *T*. The kinetic energy is the sum of all squared velocity components.
2. Calculate and print the new kinetic energy and check that it now indeed equals *N* $\times$ *T*. This command can be removed once the test has been successfully passed. 
3. Calculates the mean of all velocity vectors meanVX, meanVY and afterwards sets VX -= meanVX and VY -= meanVY to make sure the mean velocity is zero.
2. returns the modified VX and VY.

## def fold(*L,x*)

Because coordinates may leave the central simulation box, one has to fold them back to the central simulation box 

     def fold(L,x):
        x -= - L*round(x/L)
        return x
        
## def get_dist(L,x1,y1,x2,y2)

This function calculates the true vector between two nodes at positions (x1,y1) and (x2,y2). Because the simulation box has periodic boundaries, the true distance vector and its length are 

    def get_dist(L,x1,y1,x2,y2):
        distance_vector = fold(L,[x1,y1] - [x2,y2])
        distance = norm(distance_vector)
        return distance_vector,distance
    
## def get_forces(N,L,X,Y)

This function calculates force components FX and FY for all nodes. The force on node (*n,m*) has two qualitatively different contributions
1. an elastic force from the nodes that have been direct neighbors of (*n,m*) directly after create_initial_configuration(). This means, that node (*n,m*) is spring-connected (permanently bonded) to nodes (*n*+1,*m*), (*n*-1,*m*), (*n*,*m*+1), and (*n*,*m*-1). Because of periodic boundary conditions, *n*+1 stands for 1, if *n=N*, and *n*-1 stands for *N*-1, if *n*=0. 
2. a Lennard-Jones-type repulsive force from all those nodes that are currently closer to (*n,m*) than distance *r<sub>c</sub>*. 

This is the difficult part. To calculate distances and distance vectors, always use get_dist. If the distance vector **r** between two bonded particles is known, the elastic force is **F** = -*k***r** on one of the two particles (and -**F** on the other), where *k* is the spring coefficient. The Lennard-Jones force between any pair of particles (two particles form a pair if their distance is &le;*r<sub>c</sub*) is **F** = -&nabla;*E<sub>pair</sub>*.

## def single_MD_step(....)

    t += dt
    etc. as described in text

## def order_parameter(*N,L,X,Y*)

For each of the nodes check if it has at least one other node in its neighborhood, i.e., at a distance smaller 1.5. Count the number of nodes belonging to this species and divide by N<sup>2</sup> (the total number of nodes). This defines the order parameter &Phi; which is returned by this function.

## def visulize_configuration(....)

    visualize the configuration, mention the order parameter in the figure.

## def molecular_dynamics(*N,L,g,k,T*,MDsteps,*dt*)

Now we are ready to implement a full molecular dynamics with MDsteps steps, each of time duration $dt$. The overall structure looks like this, using the above functions. The integration scheme implemented here is known as the velocity-Verlet algorithm. 

    create_initial_configuration
    create_initial_velocities
    thermostat
    get_forces
    t = 0 

    # Loop over the following for MDsteps steps: 
    single_MD_step

    # eventually, at some intervals, visualize the configuration and eventually save a picture. 
    
return the final t,coordinates, velocities, and forces. 

# Tasks

1. reproduce one or more of the above configurations qualitatively (you have different initial velocities)
2. plot &Phi; versus time to check if &Phi; reaches a stationary value
3. create a set of snapshots visualizing the dynamical evolution of one of the systems shown above.
4. check if and how the result depends on initial velocities
5. plot steady-state values of &Phi; versus temperature *T* using *N=30*, *g=3.5*, *k=0.05* over the range of temperatures *T* in [0.1,0.6].

## Optional tasks

1. Calculate the time your code needs to do a single molecular dynamics step, and divide it by the number of nodes, *N*<sup>2</sup>. Report this value. 
2. Create a routine that calculates all forces **F** with a computational effort that is proportional to the number of nodes, and not quadratic in the number of nodes (hint: neighbor lists). This way, your code could run at a larger *N*.
3. Create a movie
4. Calculate the size of the largest cluster of HDP particles, which defines another order parameter, and run the above application 5. for this new order parameter.

### References

[1] http://doi.org/10.1209/0295-5075/77/58007

[2] Movie: https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/Media1.mp4

[3] Movie: https://ctl.polyphys.mat.ethz.ch/CTL-I-PUBLIC/soft-solid/Media2.mp4


