#!/usr/bin/python

class Chemostat:
	def __init__(self, V, F_in, timestep):
		self.V 					= V 				#Chemostat volume (liter)
		self.F_in 				= F_in				#Inflow rate (liter/hour)
		self.D 					= self.F_in/self.V 	#Dilution rate (1/hour)
		self.timestep 			= timestep			#size of timestep (hour)
		self.inflowMedium		= dict()			#Concentrations of metabolites in inflow medium (mM)
		self.environmentmM  	= dict()			#Concentrations of metabolites in chemostat medium: environment for metabacteria (mM)
		self.environmentmmol 	= dict()			#Abundances of metabolites in mmol in chemostat medium: environment for metabacteria (mmol)
		self.community 	 		= None				#Community inhabiting chemostat

	#Inflow of metabolites from the inflow medium to the chemostat environment
	def AddMetabolites(self):
		for met in self.inflowMedium:
			self.environmentmM[met] += self.D*self.timestep*self.inflowMedium[met]	#update concentrations in chemostat
			self.environmentmmol[met] = self.environmentmM[met]*self.V 				#Update mmol in chemostat

	#Calculate mM from mmol
	def calcConcentrations(self):
		for met in self.environmentmmol:
			self.environmentmM[met] = self.environmentmmol[met]/self.V 				#Update mM in chemostat
	
	#Calculate mmol from mM
	def calcMoles(self):
		for met in self.environmentmM:
			self.environmentmmol[met] = self.environmentmM[met]*self.V 				#Update mmol in chemostat

	#Initialize an empty chemostat environment, this is NOT the inflow medium!
	def InitEnvironment(self):
		for met in self.inflowMedium:
			self.environmentmM[met] = 0
			self.environmentmmol[met] = self.environmentmM[met]*self.V

	#Print the metabolites in the chemostat environment
	def PrintEnvironment(self, time, metabolites = False):

		#determine whether to print all metabolites, or only metabolites specified by user
		if not metabolites:
			metList = self.environmentmM
		else:
			metList = metabolites

		#Print metabolites tab separated
		printline = '%s\t' % (time)
		for met in metList:
			printline += '%s\t%2.2f mM\t%2.2f mmol\t\t' % (met, self.environmentmM[met], self.environmentmmol[met])
		print printline

	#Decrease metabacterium size due to outflow with dilution rate
	def RemoveBacteria(self):
		for metabact in self.community.members:
			metabact.size -= metabact.size*self.D*self.timestep
			self.community.sizeTot -= metabact.size*self.D*self.timestep

	#Decrease metabolite concentration due to outflow with dilution rate
	def RemoveMetabolites(self):
		for met in self.inflowMedium:
			self.environmentmM[met] -= self.D*self.timestep*self.environmentmM[met]
			self.environmentmmol[met] = self.environmentmM[met]*self.V 				#Update mmol in chemostat
