#!/usr/bin/python

# Title: 				A computational model of microbial community assembly in a gut-like environment.
# Author: 				Guus Martijn Teunisse
# Date last changed: 	March 27 2016

#Import libraries
import cobra
from numpy import random 	#Numpy's random has the same random number generator algorithm as MATLAB.
import timeit
import os
from copy import deepcopy	#Use deepcopy to copy entire objects and their attributes.
import argparse

#Import classes
from class_chemostat import Chemostat
from class_community import Community
from class_metabacterium import Metabacterium
from class_data import Data


#Get input arguments from command line or use default
def ArgumentParse():
	parser = argparse.ArgumentParser(description = "Simulate a microbial community in a chemostat.",
										epilog = "This script simulates a chemostat environment"
										" containing a community of microbial species that feed on and excrete up to 115"
										" free metabolites. For every timestep (dt), fluxes are calculated using FBAwMC, "
										"implemented in cobrapy using glpk as solver. Iteratively, species are added"
										"to the chemostat when a steady state it reached.")

	parser.add_argument('-g', metavar = 'growth threshold', dest = 'growthThreshold', help = 'If the growth of each species is below -g, '
										'a steady state has been reached and a new species is added.', default = 0.0000001, type = float)
	#increasefactor of 10,000 has been chosen so that even when a process has the lowest crowding coeff, it is still higher than the highest
	#default crowding coeff.
	parser.add_argument('-i', metavar = 'increase factor', dest = 'increaseFactor', help = 'Factor by which crowding coefficients '
										'are increased when a process is knocked-out', default = 10000, type = float)
	parser.add_argument('-n', metavar = 'N', dest = 'ntot', help = 'Number of metabacteria species to be added to chemostat', default = 1,
										 type = int)
	parser.add_argument('-o', metavar = 'Output path', dest = 'outputPath', help = 'subdirectory to store data in', default = 'output',
										 type = str)
	parser.add_argument('-p', metavar = 'Knockout probability', dest = 'pKO', help = 'Knockout probabilty for every reaction', 
										 default = 0, type = float)
	parser.add_argument('-s', metavar = 'simulation name', dest = 'simName', help = 'outputpath>simulation name directory',
										default = 'data', type = str)
	parser.add_argument('-S', metavar = 'seed', dest = 'seed', help = 'seed for random number generator',
										default = 1, type = int)
	parser.add_argument('-r', metavar = 'removal threshold', dest = 'removalThreshold', help = 'Threshold community density for removing'
										'a metabacterium from the community: invasion failure', default = 0.0001, type = float)
	parser.add_argument('-d', metavar = 'dt', dest = 'timeStep', help = 'length of timestep (dt) in hours', default = 0.1, type = float)
	parser.add_argument('-I', metavar = 'initial size', dest = 'initSize', help = 'Initial size of a metabacterium in gDW',
										 default = 0.0001, type = float)
	parser.add_argument('-u',  dest = 'upEx', action = 'store_true', help = 'Only knock-out '
										'exchange reactions; default all reactions', default = False)
	parser.add_argument('-v', metavar = 'protein volume', dest = 'VProt', help = 'Fraction of cellular volume occupied by protein',
										default = 0.2, type = float)
	parser.add_argument('-V', metavar = 'chemostat volume', dest = 'V_chemostat', help = 'Chemostat volume in liters', default = 1,
										 type = float)
	parser.add_argument('-f', metavar = 'inflow rate', dest = 'F_in', help = 'Inflow rate of chemostat in liters/hour', default = 0.1,
										 type = float)

	args = parser.parse_args()

	del parser 		#Clear some space in working memory
	return(args)

#When steady state occurs: assess whether invasion was successful or not, whether extinction occurred or not, and write to file. 
def AssessInvasion(community, metabact, initialMembers, simName, pKOStr, simNStr, time, pKO, communitySize, last_change):
	
	#Get an overview of all members in community
	finalMembers = list()
	for member in community.members:
		finalMembers.append(member.n)

	#Assess whether the current invader was succesfull
	if metabact.n in finalMembers:	#Invasion success
		invasion = 1
		invasion_time = metabact.n - last_change

	else:							#Invasion failed
		invasion = 0
		invasion_time = 0
	
	#Count the number of species that went extinct
	extinct = 0
	for member in initialMembers:
		if member not in finalMembers:
			extinct += 1

	#Calculate how many invaders were required to cause extinction of one or more species
	if extinct > 0:
		extinct_time = metabact.n - last_change
	else:
		extinct_time = 0

	#Write information about this invasion attempt to file
	WriteInvasion(simName, pKOStr, simNStr, time, metabact.n, initialMembers, invasion, extinct, pKO, communitySize, invasion_time, extinct_time)

	#Store at which invader the last change in community size occurred
	if invasion > 0 or extinct > 0:
		last_change = metabact.n
		if invasion > 0:
			invasibilityFile = open(os.path.join(simName, pKOStr, 'invasibility.txt'), 'a')
			invasibilityFile.write("%.1f\t%d\t%d\n" % (metabact.pKO, len(initialMembers), invasion_time))
			invasibilityFile.close()
		if extinct > 0:
			stabilityFile = open(os.path.join(simName, pKOStr, 'stability.txt'), 'a')
			stabilityFile.write("%.1f\t%d\t%d\n" % (metabact.pKO, len(initialMembers), extinct_time))
			stabilityFile.close()

	return last_change

#Check whether simulation directories exist and if not, create new ones
def CreateDirectories(simName, simN, pKOStr):
	#Check whether folder for simulation exists
	if not os.path.exists(simName):
		os.makedirs(simName)								#Create simulation dir

	#Check whether directories for pKO and the current repeat exist, and if not, 
	#create new directories.
	if not os.path.exists(os.path.join(simName, pKOStr)):
		if simN == None:
			simN = 1
		simNStr = '%04d' % simN 							#Store repeat number as string i.e. 0001
		os.makedirs(os.path.join(simName, pKOStr))			#Create subdir for pKO
		os.makedirs(os.path.join(simName, pKOStr, simNStr))	#Create subdir for first repeat
		simN = 1 											#Store repeat number
		simNStr = '%04d' % simN 							#Store repeat number as string
	else:
		if simN == None:
			#List all directories
			simList = os.listdir(os.path.join(os.getcwd(), simName, pKOStr))
			#See whether directories are repeat numbers, and if so, get the last repeat number
			for i in range(len(simList)):
				try:
					simList[i] = int(simList[i])
				except:
					simList[i] = 0
			simN = max(simList)
		
		#Create subdirectory for this repeat
		simN += 1
		simNStr = '%04d' % simN
		try:
			os.makedirs(os.path.join(simName, pKOStr, simNStr))
		except OSError:
			print 'OSError with seed:', simN
	return(simN, simNStr)

#Creates a text file that stores basic information about the simulation
def CreateSupplement(simName, pKOStr, simNStr, increaseFactor, ntot, pKO, initSize, timeStep, upEx, VProt, F_in, V_chemostat, seed):
	pKO = 'Variable' 	#pKO can be anything in this simulation

	path = os.path.join(simName, pKOStr, simNStr)
	supplFile = open(os.path.join(path, 'supplementary_information.txt'), 'w')
	supplFile.write('increase factor\t%s\n' % increaseFactor)
	if upEx:
		supplFile.write('KO reactions\texchange\n')
	else:
		supplFile.write('KO reactions\tall\n')
	supplFile.write('KO probability\t%s\n' % pKO)
	supplFile.write('seed\t%d\n' % seed)
	supplFile.write('dt (hours)\t%s\n' % timeStep)
	supplFile.write('initial community size (g DW)\t%s\n' % initSize)
	supplFile.write('initial number of species\t%d\n' % ntot)
	supplFile.write('Protein fraction of cellular volume\t%s\n' % VProt)
	supplFile.write('Inflow rate (liter/hour)\t%s\n' % F_in)
	supplFile.write('chemostat volume (liter)\t%s\n' % V_chemostat)
	supplFile.close()

#perform FBA for timestep, update environment and sizes of metabacteria
def Simulate(time, time_added, community, converged, chemostat, growthThreshold, removalThreshold):
	
	#Assess whether the community has reached a steady state
	if time - time_added > 10:									#Give invader at least one hour before assessing its growth
		
		steadyState = 0 										#Number of metabacteria that are growing at a steady state rate
		
		for n in range(len(community.members)-1, -1, -1): 		#Count down to prevent indexing issues when removing metabacteria
			
			member = community.members[n]
			
			#Check whether change in community density is lower than threshold
			if abs(member.size - member.prevsize) < growthThreshold:
				steadyState += 1
			
			#Metabacteria with a size close to 0 are removed to save resources --> invasion failed!
			if member.size < removalThreshold:
				community.members.pop(n)

		if steadyState == len(community.members):	#All members are in steady state
			converged = True

	#Write to files when: a new species is added; when one hour of simtime had passed; or when converged.
	written = False
	if time == time_added:
		community.WriteEnvironment(time)	#Write environment to file
		community.WriteSizes(time)
		community.WriteCrossfeeding(time)
		written = True

	if time % 1 < 0.1 and not written:		#Only write once every 1 hours of simtime
		community.WriteEnvironment(time)
		community.WriteSizes(time)
		community.WriteCrossfeeding(time)						
		written = True

	if converged and not written:			#Write when steady state occurs
		community.WriteEnvironment(time)
		community.WriteSizes(time)
		community.WriteCrossfeeding(time)

	#If chemostat is empty, stop and add a new bacterium
	if len(community.members) == 0:
		converged = True
		return(converged)
	
	#Add metabolites to environment by inflow
	chemostat.AddMetabolites()						#Inflow of metabolites to environment

	#Perform FBA
	community.SolveAllLPs()							#Perform FBA for every metabacterium. FBA is just an application of linear programming (LP)
	
	#Update metbacteria and environment using FBA solution
	community.UpdateEnvironment()					#Update environment
	chemostat.calcConcentrations()					#Update environmental concentrations
	community.UpdateSizes()							#Update size of each metabacterium

	#Update chemostat concentrations and metabacterium sizes due to outflow
	chemostat.RemoveMetabolites()
	chemostat.RemoveBacteria()

	#Print fluxes every now and then, or when steady state is reached, or when a new species is added
	#Fluxes are only written at this point, because they require LPs to be solved first.
	if time % 1 < 0.1 or converged or time == time_added:
		community.WriteFluxes(time)

	return(converged)

#Write invasion success and extinctions causes to file
def WriteInvasion(simName, pKOStr, simNStr, time, n, initialMembers, invasion, extinct, pKO, communitySize, invasion_time, extinct_time):
	#Specify path of file
	invasionFilePath = os.path.join(simName, pKOStr, simNStr, 'invasion_full.txt')

	#If new file, first write header to file
	if not os.path.exists(invasionFilePath):
		invasionFile = open(invasionFilePath,  'w')
		invasionFile.write('time\tpKO\tinvader\tcommunity\tcommunity_size\tsuccess\textinct\tsize\tinvasion_time\textinction_time\n')
		invasionFile.write('\n')				#Add new line for values
		invasionFile.close()

	invasionFile = open(invasionFilePath,  'a')
	invasionFile.write('%.1f\t%.1f\t%d\t%s\t%d\t%d\t%d\t%f\t%d\t%d\n' % (time, pKO, n, initialMembers, len(initialMembers), invasion, extinct, communitySize, invasion_time, extinct_time))
	invasionFile.close()

#This is where the magic happens
def main():
	start_time = timeit.default_timer()

	#Get parameters, optionally from command line
	parameters = ArgumentParse()
	
	#Simulation parameters
	growthThreshold 	= parameters.growthThreshold	#If growth is lower than threshold, assume steady state (default = 0.0000001)
	increaseFactor 		= parameters.increaseFactor		#Factor by which crowding coefficients are increased by knockout (default = 10000)
	ntot 				= parameters.ntot 				#Pool size (default = 1). Required for initializing community(stems from batch model)
	outputPath			= parameters.outputPath			#output path (default = 'output')
	pKO					= parameters.pKO 			 	#Initial knockout probability (default = 0)
	pKOStr				= 'pKO_%03d' % (pKO*100)		#pKO in string format as percentage. Used for creating directories
	#removalThreshold 	= parameters.removalThreshold	#Size at which a metabacterium will be removed from the community: invasion failure 
	removalThreshold 	= parameters.initSize 			#In case of assembly experiments, this is the initial size.
	simN 				= None							#Number of the current simulation. I chose to use seeds for this
	simName 			= parameters.simName			#Simulation name (default = data)
	initSize 			= parameters.initSize			#Initial metabcaterium size (default = 0.0001 gDW)
	timeStep 			= parameters.timeStep			#length of dt (default = 0.1 hour)
	upEx 				= parameters.upEx				#Only uptake and excretion reactions eligible for knockout (default = False)
	VProt 				= parameters.VProt				#Fraction of cellular volume occupied by protein. Constraint for FBAwMC (default = 0.2)
	wallTime 			= 10							#hours. Required for terminating simulation before walltime on Lisa runs out.

	#chemostat parameters
	F_in				= parameters.F_in				#Inflow rate in liters/hour (default = 0.1 liter/hour)
	V_chemostat 		= parameters.V_chemostat		#Chemostat volume in liters (default = 1 liter)

	#Get static data
	simName 			= os.path.join(outputPath, simName)		#Create output path name
	data 				= Data()								#Read data from metabacterium folder using class Data	
	random.seed(parameters.seed)								#set seed for random number generator

	#Create directories to store data
	simN, simNStr = CreateDirectories(simName, simN, pKOStr)		

	#Write simulation parameters to file
	CreateSupplement(simName, pKOStr, simNStr, increaseFactor, ntot, pKO, initSize, timeStep, upEx, VProt, F_in, V_chemostat, parameters.seed)

	#Create a community. This is required for creating a template bacterium, and for setting the environment.
	community = Community(data, simName, simNStr, initSize, 1, increaseFactor, VProt, upEx, pKO, pKOStr, timeStep)	#Create community
	community.InitCommunity(data, increaseFactor, initSize, VProt, upEx, pKO)										#InitCommunity stems from the batch model. Required for setting the environment
	community.InitEnv()
	community.members = list() #Remove all members from the community
	
	#Create a chemostat and its components. A community resides within a chemostat.
	chemostat = Chemostat(V_chemostat, F_in, timeStep)				#Create chemostat
	chemostat.community = community 								#Point to community so that they're linked
	chemostat.inflowMedium = deepcopy(community.environment)		#Set the inflow medium concentrations
	chemostat.inflowMedium['EX_glc(e)_back'] = 1 					#1 mM glucose in inflow medium
	chemostat.InitEnvironment()										#Initialize the chemostat environment
	
	#Point community environment to chemostat environment so that they're linked
	community.environment = chemostat.environmentmmol

	#Create files to store invasibility and stability
	if not os.path.exists(os.path.join(simName, 'invasibility.txt')):
		try:
			invasibilityFile = open(os.path.join(simName, pKOStr, 'invasibility.txt'), 'w')
			stabilityFile = open(os.path.join(simName, pKOStr, 'stability.txt'), 'w')
			invasibilityFile.write('pKO\tcommunity size\tinvasions\n')
			stabilityFile.write('pKO\tcommunity size\tinvasions\n')

			invasibilityFile.close()
			stabilityFile.close()
		except OSError:
			pass


	print('Started simulation with seed: %d' % (parameters.seed))

	#Iteratively add a metabacterium to the community and simulate until steady state
	time = 0
	last_change = 0 	#Last metabacterium that caused a change in community size
	for n in range(ntot):

		print("Invasion %d" % n)

		#Store initial community
		initialMembers = list()
		for member in community.members:
			initialMembers.append(member.n)
		communitySize = community.sizeTot

		#Add new member to community
		community.AddMetabacterium(n, data, increaseFactor, initSize, VProt, upEx, pKO)
		metabact = community.members[len(community.members)-1]

		#Write crowding coefficients to file
		metabact.WriteCCs(simName, pKOStr, simNStr)

		#Simulate until all metabacteria reach a steady state growth rate
		time_added = time
		converged = False

		while not converged:
			converged = Simulate(time, time_added, community, converged, chemostat, growthThreshold, removalThreshold)	#Simulate the current timestep

			time += timeStep 	#Prepare for next timestep

			#Quit simulation if walltime is about to expire
			if (timeit.default_timer() - start_time) > (wallTime * 60 * 60):
				last_change = AssessInvasion(community, metabact, initialMembers, simName, pKOStr, simNStr, time, pKO, communitySize, last_change)
				print('Ended simulation with seed: %d after %f seconds' % (parameters.seed, (timeit.default_timer() - start_time)))
				exit()
			
		last_change = AssessInvasion(community, metabact, initialMembers, simName, pKOStr, simNStr, time, pKO, communitySize, last_change)
		
	#Free up working memory
	del(community) 		
	del(chemostat)

#Only call main function when script is called as program from command line
if __name__ == '__main__':
	main()




