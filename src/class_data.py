import os

class Data:
	
	#Get parameters from files. Bounds, S, Z and C depend on metabolites and rxns
	def __init__(self):
		self.CCdist 		= self.ReadCrowdingCoef()
		self.metabolites 	= self.ReadMetabolites()
		self.rxns 			= self.ReadRxns()

		self.S 				= self.ReadS(self.metabolites, self.rxns)
		self.rxnBounds 		= self.ReadBounds(self.rxns)
		self.Z 				= self.ReadZ(self.rxns)
		self.C				= self.ReadC(self.rxns)

	#Read file with 100000 randomly generated crowding coefficients
	def ReadCrowdingCoef(self):
		crowdFile = open(os.path.join('metabacterium_files', 'crowding_coli.txt'), 'r')
		CC = [float(line.rstrip('\n')) for line in crowdFile]
		crowdFile.close()
		return(CC)

	#Read number of carbons per exchange reaction
	def ReadC(self, rxns):
		CFile = open(os.path.join('metabacterium_files', 'carbon.txt'), 'r').read().split('\r')
		C = dict()
		
		for x in range(len(rxns)):
			if 'EX_' in rxns[x]:
				value = CFile[x]
				C[rxns[x]] = value

		return C

	#Read metabolites in model
	def ReadMetabolites(self):
		#Open file containing all metabolites (679)
		metFile = open(os.path.join('metabacterium_files', 'metabolites.txt'), 'r')
		metabolites = metFile.read().split('\r')
		metFile.close()
		metabolites.append('crowding_coeff')			#Add a row to contain crowding_coefficients
		return(metabolites)

	#Read reactions in model
	def ReadRxns(self):
		#Open file containing all reactions (1164)
		rxnFile = open(os.path.join('metabacterium_files', 'reactions.txt'), 'r')
		rxns = rxnFile.read().split('\r')
		rxnFile.close()
		return(rxns)

	#Read stoichiometry matrix
	def ReadS(self, metabolites, rxns):
		#Open file containing stoichiometry matrix (679 rows, 1164 cols)
		stoichFile = open(os.path.join('metabacterium_files', 'stoichiometry.txt'), 'r')
		stoich = stoichFile.readlines()
		stoichFile.close()

		#Store stoichiometry in a dictionary 'matrix' per metabolite and reaction
		S = dict()								#Store in dictionary per metabolite
		for i in range(0, len(stoich)):
			stoichLine = stoich[i].split()		#Separate all values
			S[metabolites[i]] = dict()			#Store in dictionary per reaction
			for j in range(0, len(stoichLine)):
				S[metabolites[i]][rxns[j]] = float(stoichLine[j])
		return(S)

	#Read lower (lb) and upper bounds (ub) for all reactions
	def ReadBounds(self, rxns):
		#Open files containing lower and upper bounds of all reactions
		lbFile = open(os.path.join('metabacterium_files', 'rxn_lower_bounds.txt'), 'r')
		lb = lbFile.readlines()
		lbFile.close()
		ubFile = open(os.path.join('metabacterium_files', 'rxn_upper_bounds.txt'), 'r')
		ub = ubFile.readlines()
		ubFile.close()

		rxnBounds = dict()
		for i in range(0, len(rxns)):
			rxnBounds[rxns[i]] = {'lb':float(lb[i].rstrip()), 'ub':float(ub[i].rstrip())}
		return(rxnBounds)

	#Read coefficients for every reaction for the objective function Z
	def ReadZ(self, rxns):
		#Open file containing coefficients for objective function Z for all reactions
		zFile = open(os.path.join('metabacterium_files', 'objective_function_coeff.txt'), 'r')
		z = zFile.readlines()
		zFile.close()
		Z = dict()
		for i in range(0, len(rxns)):
			Z[rxns[i]] = float(z[i].rstrip())
		return Z