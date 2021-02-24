# Sequential community assembly with dFBAwMC

## Introduction
This tool, conveniently packaged in a Docker image, simulates a community of microorganisms using dynamic flux balance analysis with molecular crowding (dFBAwMC). Microorganisms in this system are defined as a metabolic network with random parts knocked out. To form a community, microorganisms are added one-by-one and simulated until a steady state or extiction occurs, after which a new microorganism is added. The organisms live in a chemostat-like environment with inflow of only a single limiting carbon resource, glucose.

## Installation
To install this tool, install Docker on your system and pull from Docker Hub using the following command:
```
docker pull gmteunisse/fba-community-assembly:0.1.0
```

Alternatively, you can build a local version of this tool by cloning this repository and building a local docker image:

```
git clone https://github.com/gmteunisse/fba-community-assembly.git
cd fba-community-assembly
docker build <username>/fba-community-assembly:0.1.0
```

Lastly, the code can downloaded and run locally. This requires the following dependencies:
* Python 2.x (will not run on Python 3.x)
	* NumPy
	* SciPy
	* lxml
	* libsbml
	* swiglpk
	* cobrapy v0.3.2
* GLPK


## Simulation
To run a simulation with the default settings, run:

```
docker run -v <results_path>:/output gmteunisse/fba-invasion:0.0.1 -o /output
```

To get an overview of parameters, add `--help` to his command, which displays the help file below:

```
usage: assembly_model.py [-h] [-g growth threshold] [-i increase factor]
                         [-n N] [-o Output path] [-p Knockout probability]
                         [-s simulation name] [-S seed] [-r removal threshold]
                         [-d dt] [-I initial size] [-u] [-v protein volume]
                         [-V chemostat volume] [-f inflow rate]

Simulate a microbial community in a chemostat.

optional arguments:
  -h, --help            show this help message and exit
  -g growth threshold   If the growth of each species in the community is below -g, a steady
                        state has been reached and a new species is added.
  -i increase factor    Factor by which crowding coefficients are increased
                        when a process is knocked-out
  -n N                  Number of metabacteria species to be added to
                        chemostat
  -o Output path        subdirectory to store data in
  -p Knockout probability
                        Knockout probabilty for every reaction
  -s simulation name    outputpath>simulation name directory
  -S seed               seed for random number generator
  -r removal threshold  Threshold community density for removinga
                        metabacterium from the community: invasion failure
  -d dt                 length of timestep (dt) in hours
  -I initial size       Initial size of a metabacterium in gDW
  -u                    Only knock-out exchange reactions; default all
                        reactions
  -v protein volume     Fraction of cellular volume occupied by protein
  -V chemostat volume   Chemostat volume in liters
  -f inflow rate        Inflow rate of chemostat in liters/hour

This script simulates a chemostat environment containing a community of
microbial species that feed on and excrete up to 115 free metabolites. For
every timestep (dt), fluxes are calculated using FBAwMC, implemented in
cobrapy using glpk as solver. Iteratively, species are addedto the chemostat
when a steady state it reached.
```