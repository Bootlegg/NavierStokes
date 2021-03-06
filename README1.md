#Fluid Dynamics
This sections contains various different subjects.  

1. [Shallow Water Equations](https://github.com/mintDan/FluidDynamics#shallow-water-equations)
2. [NavierStokes](https://github.com/mintDan/FluidDynamics#navierstokes)
3. [Newton's Method for system of nonlinear equations](https://github.com/mintDan/FluidDynamics#newtons-method-for-system-of-nonlinear-equations)
4. [Delaunay Triangulation for mesh generation](https://github.com/mintDan/FluidDynamics#delaunay-triangulation-for-mesh-generation)

# Shallow Water Equations
Solutions based on different schemes.  
Here is Semi-Lagrangian bi-cubic interpolation, with filter for short wavelengths.
Periodic grid with a mountain on the bottom of the ocean. The governing equations are 

![SWE.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/SWE.png)

As seen these are the Lagrangian form of the equations, following trajectories of water particles.
## Water Waves
![SWE2D.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/SW2D.png)

## Mountain at the seabed
![Bottomtopo.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/bottomtopog.png)


# NavierStokes
Based on Lorena Barba's excellect Navier-Stokes steps.
Output is written to files and then animated.

The momentum equations are given by  

![Momentum.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/Momentum.png)

The pressure equation is given by  

![Peq.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/Peq.png)

The pressure equation can be solved by pseudotime methods, such as

![Pseudotime.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/Pseudotime.png)

The result of velocity(arrows) and pressure(contour) fields is given below 

![NV.png](https://github.com/mintDan/FluidDynamics/blob/master/figs/NV.png)

# Newton's Method for system of nonlinear equations
Solves n equations in n variables using Newton's method. At the moment uses constant stepping dx, an idea for improvement would be a dynamic stepping.
Calls on NumPy to invert a matrix, this could be changed to be done with Jacobi method instead.  

The iteration process to solve the system of equations is

![NewtonMethod](https://github.com/mintDan/FluidDynamics/blob/master/figs/NewtonMethod.png)

Where J is the Jacobian given below.

![Jacobian](https://github.com/mintDan/FluidDynamics/blob/master/figs/Jacobian.png)

# Delaunay Triangulation for mesh generation
Self made mesh generation from Delaunay Triangulation. A direct method and not efficient, but gets the job done. Should try more efficient algos like divide and conquer.

![Jacobian](https://github.com/mintDan/FluidDynamics/blob/master/figs/DT.png)