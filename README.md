#Fluid Dynamics
This sections contains various different subjects.  

1. [Shallow Water Equations](https://github.com/mintDan/Fluiddynamicss#shallowwaterequations)
2. [NavierStokes](https://github.com/mintDan/Fluiddynamicss#navierstokes)


# ShallowWaterEquations
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

![NV.png](https://github.com/Bootlegg/FluidDynamics/blob/master/figs/NV.png)