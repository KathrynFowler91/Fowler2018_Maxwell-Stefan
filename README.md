# Maxwell–Stefan diffusion: a framework for predicting condensed phase diffusion and phase separation in atmospheric aerosol (supporting code)

Kathryn Fowler, Paul J. Connolly, David O. Topping, and Simon O'Meara

School of Earth and Environmental Sciences, The University of Manchester, Manchester, UK

Correspondence: kathryn.fowler@manchester.ac.uk / p.connolly@manchester.ac.uk

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1156928.svg)](https://doi.org/10.5281/zenodo.1156928)

This is the first time the Maxwell-Stefan framework has been applied to an atmospheric aerosol core-shell model and shows that there is a complex interplay between the viscous and solubility effects on aerosol composition. Understanding aerosol composition is essential to accurately model their interactions within atmospheric systems. We use simple binary systems to demonstrate how viscosity and solubility both play a role in affecting the rate of diffusion through aerosol particles.

The complete article can be found: <https://www.atmos-chem-phys-discuss.net/acp-2017-424/#discussion>

### Description

This repository contains a function called 'implicit_maxwell_stefan03', which solves the diffusion equation on a spherical core shell model with neumann boundary conditions. Fick's second law is solved using the backward Euler method of finite differences and the non-ideal effects of diffusion are included using the Maxwell-Stefan framework, relating diffusion flux to gradient in activity coefficients. The UNIFAC group contribution model is used to estimate the activity coefficients.

Functions included by the repository, which are called by 'implicit_maxwell_stefan03' include:
 - outer_shell_equilibration.m
 - outer_shell_redistribution.m
 - UFC_datamain.m
 - UFC_dataqi_v1.m
 - UFC_interact_params_v1.m
 - UNIFAC_gamma.m

### Getting Started

Call the model by typing [u_save, RN] = implicit_maxwell_stefan03()

Input parameters maybe changed in the 'setting up parameters' section in the main function or using the 'implicit_maxwell_stefan_runscript' function.

Input varibles are:
 - fickian =  whether the fickian or Maxwell-Stefan framework of diffusion is used.
 - diffusion = the diffusion coefficient equation ('darken' or 'vignes').
 - n_components = the number of components.
 - dt = time-step [seconds].
 - ip = number of grid points.
 - ntm = number of time steps.
 - DSelf = self-diffusion coefficent of non-volatile component.
 - aerosol = the Maxwell-Stefan variable switch (see MS variable switch).
 - CNum = number of carbon atoms in monocarboxylic acid chain (see MS variable switch).
 - Xw_init = initial water mole fraction of particle (0-1).
 - Xw_shell(:) = the equilibrated water mole fraction for the outer shell, defined for each time step (0-1).
 - r0 = lowest radius edge [m].
 - rN = upper radius edge [m].
 - R = initial radius of aerosol particle [m].
 - T = temperature [Kelvin]

Outputs variables are:
 - u_save = concentration (or molar density) for a given shell, timestep and component.
 - RN = inital shell boundaries as determined from the number of shells and size of grid.

### Contact

If you have any further issues using the model please contact either Kathryn Fowler (kathryn.fowler@manchester.ac.uk) or Dr. Paul Connolly (p.connolly@manchester.ac.uk).

### License

Copyright (C) 2018  Kathryn Fowler

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>
