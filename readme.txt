--------------------------------------------------------------------------------------------
-- LibreGrowth
--------------------------------------------------------------------------------------------

	Here we introduce LibreGrowth: a libre tumor growth code for simulating spheroid and rim development at benign and malignant stages over a conditioned media. It models each medium as a three-dimensional domain with an spatially variable diffusion coefficient, through the resolution of reaction-diffusion equations. The code is implemented in C++ for GNU/Linux systems, and optimized through the shared memory technology OpenMP. LibreGrowth provides a flexible implementation in order to facilitate upcoming studies concerning to the impact of the environment over the infiltration patterns. We expect that this kind of novel research tools help in the promotion of standard cancer therapies optimization, particularly in the context of personalized medicine.


--------------------------------------------------------------------------------------------
-- Code file descriptions
--------------------------------------------------------------------------------------------

LibreGrowth code is organized in four main source files and two makefiles:

	src/main.cpp
		Main simulation process (finite difference iterative scheme) is implemented in this file.

	src/Params.h
		Model parameters (i.e. proliferation coefficient) are defined in this file.

	src/ScalarField.h
		ScalarField class is defined in this file.
		An instance of this class is used in main.h for depicting the tumor concentration.

	src/Mesh.h
		Mesh class keeps information about geometry and domain discretization.

	src/makefile
		Compilation process is implemented in this file.

	makefile
		File cleaning, compilation, execution and loggin are implemented in this file.


--------------------------------------------------------------------------------------------
-- How to modify the code
--------------------------------------------------------------------------------------------

	You can modify src/Params.h, where you will find a list of parameters to be adapted to your particular problem.

		- Initial core radius. Units: um. (r_min_core_dim)
		- Domain dim. Units: um. (x_min_dim, x_max_dim, etc.)
		- Max. concentration of tumoral cells. Units: cells/um^3. (u_max_dim)
		- Max. concentration of core cells. Units: cells/um^3. (u_core_max_dim)
		- Invasion diffusion coefficient. Units: um^2/h. (D_inv_dim)
		- Core velocity. Units: um/h. (vc_dim)
		- Proliferation coefficient. Units: 1/h. (p_coef_dim)
		- Source coefficient. Units: cell/(um^3 h). (s_dim)
		- Max. time after invasion begins: 5 days. Units: h. (t_max_inv_dim)
		- Max. time of benign stage. 10 days. Units: h. (t_max_benign_dim)
		- How often you save an output file. (save_step)


--------------------------------------------------------------------------------------------
-- Compilation and execution instructions
--------------------------------------------------------------------------------------------

	make
		File cleaning and standard compilation

	make clean
		Last compilation and execution file cleaning

	make run
		File cleaning, standard compilation and execution

	make run-omp
		File cleaning, OpenMP compilation and multi-thread execution

	make run-dbg
		File cleaning, debug compilation and execution


--------------------------------------------------------------------------------------------
-- Output
--------------------------------------------------------------------------------------------

	In the "data" folder you will find vtk and/or csv output files.

	Each file represents the tumor concentration of a particular time.

	Output formats were selected to be compatible with the powerful visualization tool: Paraview.


