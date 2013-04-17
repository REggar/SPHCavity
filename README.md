# SPHCavity

Java code for simulating lid driven cavity flows using SPH. Note that the moving is moving positively in the x direction at the y=0 instead of y=H. The data is outputted as ParaView files. ParaView is an open-source data analysis and visualization program and is available at http://www.paraview.org/.

## Getting Started
Edit the `SPHPoiseuille.java` file with your required input setting currently you are able to set:
Smoothing length
Mass of particle
Speed of sound in the fluid
Time step
Acceleration due to gravity
Kinematic viscosity of fluid
Velocity of plate at y=0

Set the initial box size using `GenerateParticles.java`

## Improvements
As with any project some improvements could be made:
- Better method of initialising particles
- Nicer method of writing ParaView files
- Implement more weighting functions
- Saving state to allow resuming

## Notes
Used for my undergraduate disseration
	

