#!/usr/bin/python

#Import libraries
import cobra
from numpy import random 	#Use random package frmo numpy, as this is the same random number generator as in MATLAB.
import os
import timeit

class Metabacterium:
	def __init__(self, n, data, increaseFactor, initSize, VProt, upEx, pKO, template, metabacterium = False):

		#Set parameters
		self.genome_size 	= 1164
		self.n 				= n
		self.increaseFactor = increaseFactor
		self.pKO			= pKO
		self.initSize		= initSize
		self.size 			= initSize
		self.prevsize		= 0
		self.upEx 			= upEx
		self.VProt			= VProt
		self.solver 		= cobra.solvers.cglpk
		self.KO 			= dict()

		#Get parameters that have been read from files. Bounds, S and Z depend on metabolites and rxns
		self.CCdist			= data.CCdist 		#100,000 crowding coefficients drawn from distribution of E. coli (list)
		self.metabolites 	= data.metabolites 	#metabolites FBA model (list)
		self.rxns 			= data.rxns 			#reactions in FBA model (list)
		self.rxnBounds 		= data.rxnBounds 	#Bounds of reactions in FBA model (dict: ub & lb per rxn)
		self.S 				= data.S 			#Stoichiometry matrix in FBA model (dict: rxns per metabolite)
		self.Z 				= data.Z 			#Coefficients of objective function in FBA model

		#Define model:
		if template: 		#The first metabacterium to be initiliazed is the template for other metabacteria. Preferably read from json file.
			if not os.path.exists(os.path.join('metabacterium_files','metabacterium.json')):	#If no file exists, initialize from text files and create json file.
				self.model 		= self.InitializeMetabacterium() 								#Initialize a template metabacterium from text files
				cobra.io.save_json_model(self.model, os.path.join('metabacterium_files','metabacterium.json'))	#Save as .json file. SMBL saving leads to loss
			else:																								#of 1 metabolite for unknown reasons
				self.model 		= cobra.io.load_json_model(os.path.join('metabacterium_files','metabacterium.json'))
		else:
			if not metabacterium:
				raise RuntimeError('Please provide a metabacterium template.')
			else:
				self.model = metabacterium.model.copy() 	#Copy the model from template. This is much faster than reading from file(s)

		#Create linear programming (LP) problem from model. LP problems are created and solved in C++.
		self.problem 		= self.solver.create_problem(self.model, objective_sense='maximize')
		self.solution 		= None 		#Attribute to store solution at time (t)
		self.growthRate 	= None		#Attribute to store growth rate at time (t)
		self.fluxes 		= dict()	#attribute to store fluxes are time (t)


	#Create a cobra model that holds all reactions and metabolites of a metabacterium.
	#Also sets the constraints of the reaction: steady state, lower and upper bounds of reactions
	def InitializeMetabacterium(self):
		#Define LP problem object
		metabact = cobra.Model('metabacterium')
		reactions = list()
		metabolites = list()

		#Define all metabolites
		for i in range(len(self.metabolites)):
			met = self.metabolites[i]
			metabolites.append(cobra.Metabolite(met))				#Add cobra metabolite object
			if met == 'crowding_coeff':								#the sum of crowding coefficients should be < VProt
				metabolites[i]._bound = self.VProt					#for all other metabolites, _bound == 0.
				metabolites[i]._constraint_sense = "L"				# "L" specifies "lower than"

		#Define all reactions, their stoichiometric coefficients per metabolite, upper bounds, lower bounds
		#and objective coefficient (z)
		for i in range(len(self.rxns)):
			rxn = self.rxns[i]
			reactions.append(cobra.Reaction(rxn))					#Add cobra reaction object
			for j in range(len(metabolites)):
				met = metabolites[j]
				s_coeff = self.S[self.metabolites[j]][rxn]			#Obtain the stoichiometric coefficients 
				reactions[i].add_metabolites({met:s_coeff})			#Store for every reaction	

			reactions[i].lower_bound = self.rxnBounds[rxn]['lb']	#Set reaction lower- (lb) and upper bounds (ub)
			reactions[i].upper_bound = self.rxnBounds[rxn]['ub']
			reactions[i].objective_coefficient = self.Z[rxn]		#Set coefficients of objective function (Z) for every reaction
																	#Currently every index is 0, except for 'biomass_LPL6.0', which is 1.
			metabact.add_reaction(reactions[i])						#Add all reactions incl. metabolites to metabact object

		return(metabact)

	#Set random crowding coefficients for every reaction, except for exchange reactions
	def SetCCs(self):
		#Update CCs in stoichiometry matrix (S)
		for rxn in self.rxns:
			KO = 0 														#Reaction knocked out: False (0); True (1).
			
			#Assign random crowding coefficients to every reactions
			if self.S['crowding_coeff'][rxn] > 0.000000001: 			#EX_fluxes (and a few other reactions) have CC = 0.000000001
				randNum = int(random.random()*99999+0.5)				#10000 random crowding coefficients, add 0.5
				self.S['crowding_coeff'][rxn] = self.CCdist[randNum]	#so that int(0.5) == 1 and not 0.
			
			#Knockout exchange reactions. Only if self.upEx = True
			if self.upEx:
				if 'EX_' in rxn and random.random() < self.pKO:					#Knock out in- & effluces with prob. pKO by
					randNum = int(random.random()*99999+0.5)					#assigning a random CC > 0.000000001 * increaseFactor
					self.S['crowding_coeff'][rxn] = self.CCdist[randNum] * self.increaseFactor
					KO = 1
					self.genome_size -= 1 										#Update genome size
			
			#Knockout non-exchange reactions
			else:
				if random.random() < self.pKO:									#Knock out all processes with prob. pKO
					self.S['crowding_coeff'][rxn] *= self.increaseFactor
					KO = 1
					self.genome_size -= 1 										#Update genome size
			
			self.KO[rxn] = KO 													#Store whether reaction was knocked out
		
		#Change CCs in cgplk object. 
		ccInd = self.model.metabolites.index("crowding_coeff") 					#Get index for crowding coefficient metabolite in cgplk object.
		for rxn in self.rxns:
			rxnInd = self.model.reactions.index(rxn)	
			self.solver.change_coefficient(self.problem, ccInd, rxnInd, self.S['crowding_coeff'][rxn])	#Update CCs for metabacterium problem

	#Set concentrations of essential compounds in the environment to infinity
	def SetEssentialEnv(self):
		#Dict containing essential metabolites and their values
		essentialMets = {'EX_h2(e)_back':1000000}
						 #These metabolites are not essential in the current environment, but may be essential in more realistic environments.
						 #{'EX_h2o(e)_back':1000000,
						 #'EX_na1(e)_back':1000000,
						 #'EX_so4(e)_back':1000000,
						 #'EX_pi(e)_back':1000000,
						 #'EX_h(e)_back':1000000,
						 #'EX_nh4(e)_back':1000000}
		
		for met in essentialMets:
			ind = self.model.reactions.index(met)											#Get  index for environmental metabolite
			self.model.reactions[ind].upper_bound = essentialMets[met]						#Update to dict value in model
			self.solver.change_variable_bounds(self.problem, ind, 0, essentialMets[met])	#Update to dict value in problem

	#Perform FBA
	def SolveLP(self):

		status = self.solver.solve_problem(self.problem)						#Solve LP
		flux = dict()

		#Optimal means a solution has been found
		if status == 'optimal':

			self.solution = self.solver.format_solution(self.problem, self.model)	#Store solution
			self.growthRate = self.solution.f 										#Store solution for objective function
			
			#Store exchange fluxes
			for x in self.solution.x_dict:
				if 'EX_' in x:
					flux[x] = self.solution.x_dict[x]
					if flux[x] < 0:		#sometimes a very small negative solution is returned
						flux[x] = 0 	#which must be set to 0, as the lower bound of each rxn = 0.
		
		#If no optimal solution is found, this usually means that the solution is the null space,
		#which means that no flux is possible for the metabacterium in the current environment.
		#However, the ATP sink in this model makes this an infeasible solution.
		#For now, I assume that such an infeasible solution leads to 0 flux for all reactions.
		else:
			self.growthRate = 0
			for x in self.rxns:
				if 'EX_' in x:
					flux[x] = 0
		
		#Store as attribute
		self.fluxes = flux


	#Update bounds of uptake reactions from environment
	def UpdateBounds(self, environment, envMets, sizeTot, timestep):
		for met in envMets:
			ind = self.model.reactions.index(met) 															#Obtain metabolite index
			self.model.reactions[ind].upper_bound = environment[met]/(timestep*sizeTot)						#Set upper bound:
			self.solver.change_variable_bounds(self.problem, ind, self.model.reactions[ind].lower_bound, environment[met]/(timestep*sizeTot))	#mmol/(population gDW*hr)
		
	#Write crowding coefficients and KO reactions to files for a single metabacterium
	def WriteCCs(self, simName, pKOStr, simNStr):
		#Define filepaths
		ccFilePath = os.path.join(simName, pKOStr, simNStr, 'crowding_coefficients.txt')
		KORxnsFilePath = os.path.join(simName, pKOStr, simNStr, 'KO_reactions.txt')

		#If file does not exist, create file and add header
		if not os.path.exists(ccFilePath) or not os.path.exists(KORxnsFilePath):
			ccFile = open(ccFilePath,  'w')				#File to store crowding coefficients
			KORxnsFile = open(KORxnsFilePath,  'w')		#File to store KO reactions
			KORxnsFile.write('metabacterium\t')
			ccFile.write('metabacterium\t')
			for rxn in self.rxns:					#Write reaction names
				KORxnsFile.write('%s\t' % rxn)
				ccFile.write('%s\t' % rxn)

			KORxnsFile.write('\n')
			ccFile.write('\n')
			ccFile.close()
			KORxnsFile.close()

		#Open files
		ccFile = open(ccFilePath,  'a')				#File to store crowding coefficients
		KORxnsFile = open(KORxnsFilePath,  'a')		#File to store KO reactions

		#Write metabacterium number
		KORxnsFile.write('%d\t' % self.n)
		ccFile.write('%d\t' % self.n)
		
		#Write CCs and KOs
		for rxn in self.rxns:
			KORxnsFile.write('%d\t' % self.KO[rxn])	#Write KO status of rxn
			ccFile.write('%f\t' % self.S['crowding_coeff'][rxn])	#Write crowding coeff.

		#Add a new line
		KORxnsFile.write('\n')
		ccFile.write('\n')

		#Close files
		ccFile.close()
		KORxnsFile.close()












