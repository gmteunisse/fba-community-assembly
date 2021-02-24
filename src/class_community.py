#!/usr/bin/python

#Import libraries
import timeit
import os
import copy

#Import classes
from class_metabacterium import Metabacterium

class Community:

	def __init__(self, data, simName, simNStr, initSize, ntot, increaseFactor, VProt, upEx, pKO, pKOStr, timestep):
		#Parameters for file writing
		self.simName 		= simName
		self.simNStr		= simNStr
		self.pKOStr			= pKOStr
		self.timestep 		= timestep	#hr

		#Community parameters
		self.sizeTot 		= initSize
		self.sizeTotPrev	= 0
		self.ntot 			= ntot
		#The template bacterium is used to copy other bacteria from. This speeds up the process of bacterium creation greatly.
		self.template 		= Metabacterium(0, data, increaseFactor, initSize, VProt, upEx, pKO, template = True)
		self.members 		= list()
		self.environment 	= dict()
		self.envMets 		= list()
		self.exchangeRxns 	= list()

		#carbon table for crossfeeding calculation
		self.carbon_table = data.C

	#Add a metabacterium to the community from template with random crowding coefficients and
	#set its environment. Note: these environmental values will be overwritten by the SolveAllLPs function
	#at the next cycle of FBA. These values are only used to set the initial environment.
	def AddMetabacterium(self, n, data, increaseFactor, initSize, VProt, upEx, pKO):
		metabact = Metabacterium(n, data, increaseFactor, initSize, VProt, upEx, pKO, template = False, metabacterium = self.template)
		glcInd = metabact.model.reactions.index('EX_glc(e)_back')					#Get problem index for environmental glucose
		metabact.model.reactions[glcInd].upper_bound = 0							#Env: 0 mmol glucose in model		
		metabact.SetEssentialEnv()													#Set essential metabolites to inf.
		metabact.SetCCs()															#Set random crowding coefficients
		
		metabact.solver.set_parameter(metabact.problem, 'lp_method', "auto")		#Set solver parameters equal to Matlab's COBRA toolbox
		metabact.solver.set_parameter(metabact.problem, 'tolerance_feasibility', 1e-9)
		
		self.members.append(metabact)												#Add to community


	#Create ntot metabacteria 
	def InitCommunity(self, data, increaseFactor, initSize, VProt, upEx, pKO):
		for n in range(self.ntot):
			self.AddMetabacterium(n, data, increaseFactor, initSize, VProt, upEx, pKO)

	#Initialize metabolite concentrations in the environment
	def InitEnv(self):
		metabact = self.members[0]								#Any metabacterium will do for initialization, as long as it's not the template:
																#the bacterium is used for its environmental values.
		for rxn in metabact.rxns:	
			if 'EX_' in rxn: 
				self.exchangeRxns.append(rxn)
				if '_back' in rxn:								#Only store external metabolites
					ind = metabact.model.reactions.index(rxn)	#Get index of reaction
					self.environment[rxn] = metabact.model.reactions[ind].upper_bound	#Store value
					self.envMets.append(rxn)					#Store reaction name

	#Update the bounds, then solve LP for every metabacterium
	def SolveAllLPs(self):
		for metabact in self.members:
			metabact.UpdateBounds(self.environment, self.envMets, self.sizeTot, self.timestep)
			metabact.SolveLP()

	#Calculate the in- and decreases of all environmental metabolites based on fluxes of all metabacteria
	def UpdateEnvironment(self):

		for metabact in self.members:

			for met in self.envMets:
				
				#Uptake by metabacteria
				self.environment[met] = self.environment[met] -  self.timestep*metabact.size*metabact.fluxes[met]		#Remove met from environment
				if self.environment[met] < 0:		#For some reason, some of the Inflow fluxes become negative, even though they are bounded by 0.
					if self.environment[met] < -0.1:	#If serious negative values occur, let user know.
						print 'influx metabact %s; metabolite %s yields negative environment! Setting to 0.' % (metabact.n, met)
					self.environment[met] = 0
				
				#Excretion by metabacteria
				metInflow = met[0:len(met)-5]			#Remove '_back' from met for flux into environment
				self.environment[met] = self.environment[met] + self.timestep*metabact.size*metabact.fluxes[metInflow]	#Add flux to environment
				if self.environment[met] < 0:		#For some reason, some of the Inflow fluxes become negative, even though they are bounded by 0.
					if self.environment[met] <-0.1:	#If serious negative values occur, let user know.
						print 'effflux metabact %s; metabolite %s yields negative environment! Setting to 0.' % (metabact.n, met)
					self.environment[met] = 0 			#Set back to 0

	#Update size of every community member and count how many bacteria are in the community
	def UpdateSizes(self):
		
		#Set community size and density to 0
		self.sizeTotPrev = self.sizeTot
		self.sizeTot = 0
		self.ntot = 0

		#Update individual sizes and sum
		for metabact in self.members:
			metabact.prevsize = metabact.size
			metabact.size += metabact.growthRate*self.timestep*metabact.size
			self.sizeTot += metabact.size
			self.ntot += 1


	#Write crowding coefficients and KO reactions for every metabacterium to files
	def WriteCCs(self):
		#Define filepaths
		ccFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'crowding_coefficients.txt')
		KORxnsFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'KO_reactions.txt')

		#If file does not exist, create file and add header
		if not os.path.exists(ccFilePath) or not os.path.exists(KORxnsFilePath):
			ccFile = open(ccFilePath,  'w')				#File to store crowding coefficients
			KORxnsFile = open(KORxnsFilePath,  'w')		#File to store KO reactions
			KORxnsFile.write('metabacterium\t')
			ccFile.write('metabacterium\t')
			for rxn in self.template.rxns:					#Write reaction names
				KORxnsFile.write('%s\t' % rxn)
				ccFile.write('%s\t' % rxn)

			KORxnsFile.write('\n')
			ccFile.write('\n')
			ccFile.close()
			KORxnsFile.close()

		#Open files
		ccFile = open(ccFilePath,  'a')				#File to store crowding coefficients
		KORxnsFile = open(KORxnsFilePath,  'a')		#File to store KO reactions


		for metabact in self.members:
			KORxnsFile.write('%d\t' % metabact.n)			#Write metabacterium number
			ccFile.write('%d\t' % metabact.n)
			for rxn in metabact.rxns:
				KORxnsFile.write('%d\t' % metabact.KO[rxn])	#Write KO status of rxn
				ccFile.write('%f\t' % metabact.S['crowding_coeff'][rxn])	#Write crowding coeff.

			KORxnsFile.write('\n')
			ccFile.write('\n')
		ccFile.close()
		KORxnsFile.close()

	#Calculate and write crossfeeding coefficients from Van Hoek & Merks (2016). Also used in Leonie's report.
	def WriteCrossfeeding(self, time):
		#Define empty dictionaries
		F_up = dict()
		F_ex = dict()
		C_crossfeeding = dict()
		C_net = dict()
		C_up = dict()
		CF_coef_ind = dict() 	#Crossfeeding coefficient per metabacterium

		C_crossfeeding_tot = 0 			#Crossfeeding coefficient for population
		C_up_tot = 0

		#Create crossfeeding file
		cfFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'crossfeeding.txt')
		if not os.path.exists(cfFilePath):
			cfFile = open(cfFilePath, 'w')
			cfFile.write('time\tmetabacterium\tCF coeff\n')
			cfFile.close()

		cfFile = open(cfFilePath, "a")

		for metabact in self.members:

			C_net[metabact.n] = 0
			C_up[metabact.n] = 0

			F_up[metabact.n] = dict()
			F_ex[metabact.n] = dict()

			for rxn in metabact.fluxes:

				if metabact.fluxes[rxn] > 0.00001: 										#only calculate crossfeeding for active fluxes

					if 'back' in rxn: 													#'back' indicates an uptake flux
						met = rxn[0:-5]													#Remove the 'back' part
						F_up[metabact.n][met] = metabact.fluxes[rxn]*metabact.size*0.1 	#Transform flux (mmol/(gDW*hr)) into mmol
						if met not in F_ex[metabact.n]:
							F_ex[metabact.n][met] = 0 									#Add 0 for excretion flux to avoid KeyErrors later on
					else:
						met = rxn
						F_ex[metabact.n][met] = metabact.fluxes[rxn]*metabact.size*0.1 	#Transform flux (mmol/(gDW*hr)) into mmol
						if rxn not in F_up[metabact.n]:
							F_up[metabact.n][met] = 0 									#Add 0 for uptake flux to avoid KeyErrors later on

					C_net[metabact.n] += float(self.carbon_table[met]) * max(float(0), (F_up[metabact.n][met]-F_ex[metabact.n][met])) #Calculate the net uptake of C in mmol by metabacterium
					C_up[metabact.n] += float(self.carbon_table[met]) * F_up[metabact.n][met]

			#Calculate the amount of carbon taken up by using glucose (NOT a crossfeeding metabolite)
			if 'EX_glc(e)' in F_up[metabact.n]:
				C_glucose = float(self.carbon_table['EX_glc(e)']) * F_up[metabact.n]['EX_glc(e)']
			else:
				C_glucose = 0

			#Calculate the amount of carbon obtained by crossfeeding
			C_crossfeeding[metabact.n] = C_net[metabact.n] - C_glucose

			#Normalize to obtain the crossfeeding coefficients for individual metabacteria
			if C_up[metabact.n] > 0:
				CF_coef_ind[metabact.n] = C_crossfeeding[metabact.n]/C_up[metabact.n]
			else:
				CF_coef_ind[metabact.n] = 0
			
			#Write to file
			cfFile.write('%f\t%d\t%f\n' % (time, metabact.n, CF_coef_ind[metabact.n]))

			#Sum the total amount of carbon obtained by crossfeeding
			C_crossfeeding_tot += C_crossfeeding[metabact.n]
			C_up_tot += C_up[metabact.n]

		#Calculate the community crossfeeding coefficient
		if C_up_tot > 0:
			CF_coef = C_crossfeeding_tot/C_up_tot
		else:
			CF_coef = 0

		#Write to file
		cfFile.write('%f\t%s\t%f\n' % (time, 'Total', CF_coef))
		cfFile.close()



	#Write environment to file
	def WriteEnvironment(self, time):
		envFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'metabolites_environment.txt')
		
		#If new file, first write header
		if not os.path.exists(envFilePath):
			envFile = open(envFilePath,  'w')
			envFile.write('time\tmetabolite\tmmol\n')
			envFile.close()

		#Write environmental values to file
		envFile = open(envFilePath,  'a')
		for met in self.envMets:
			if self.environment[met] > 0.000001:
				envFile.write('%f\t%s\t%f\n' % (time, met, self.environment[met]))
		envFile.close()

	#Write fluxes from and to the environment (exchange reactions) to file. 
	#The summed flux file sums all fluxes for a certain timepoint, but only for metabacteria that are larger than the initial size.
	#This is to prevent the fluxes from looking too "spiky"
	def WriteFluxes(self, time):
		fluxFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'fluxes.txt')
		sumfluxFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'fluxes_summed.txt')

		#If new file, first write header
		if not os.path.exists(fluxFilePath):
			fluxFile = open(fluxFilePath,  'w')
			fluxFile.write('time\tmetabacterium\treaction\tflux\n')
			fluxFile.close()
			sumfluxFile = open(sumfluxFilePath,  'w')
			sumfluxFile.write('time\treaction\tflux\n')
			sumfluxFile.close()

		fluxFile = open(fluxFilePath,  'a')
		sumfluxFile = open(sumfluxFilePath,  'a')
		
		for rxn in self.exchangeRxns:
			summed_value = 0
			for metabact in self.members:
				if metabact.fluxes[rxn] > 0.00001: 			#only write when flux is larger than 0 to save space.
					fluxFile.write('%f\t%d\t%s\t%f\n' % (time, metabact.n, rxn, metabact.fluxes[rxn]))
					if metabact.size > metabact.initSize:
						summed_value += metabact.fluxes[rxn]
			if summed_value > 0:
				sumfluxFile.write("%f\t%s\t%f\n" % (time, rxn, summed_value))
		fluxFile.close()
		sumfluxFile.close()


	#Write size for every metabacterium to file
	def WriteSizes(self, time):
		sizeFilePath = os.path.join(self.simName, self.pKOStr, self.simNStr, 'sizes_summed.txt')

		#If new file, first write metabacterium names to file
		if not os.path.exists(sizeFilePath):
			sizeFile = open(sizeFilePath,  'w')
			sizeFile.write('time\tmetabacterium\tsize\tpKO\n')
			sizeFile.close()

		#Write sizes
		sizeFile = open(sizeFilePath,  'a')
		for metabact in self.members:
			sizeFile.write('%s\t%d\t%f\t%.1f\n' % (time, metabact.n, metabact.size, metabact.pKO))

		sizeFile.write('%s\t%s\t%f\t%s\n' % (time, 'Total', self.sizeTot, 'NA'))
		#sizeFile.write('\n')			#Add new line for next timestep		
		sizeFile.close()




