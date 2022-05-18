##Python source code for ScrumPy3

from importlib import reload
from ScrumPy.Util import Set


import numpy as np
import BuildLP
reload(BuildLP)
import Media
reload(Media)
import SepiBiomass
reload(SepiBiomass)
media = Media.media

empform = {'ACETOIN':{'C': 4, 'H': 8, 'O': 2},
	   'BUTANEDIOL':{'C': 4, 'H': 10, 'O': 2}
	  }

class ModelResults():
	'''slecetion and lp constraint operations on reactions in the model'''
	def __init__(self,m):
		self.reactions = m.sm.cnames
		self.smx = m.smx

	def Reactions(self):
		'''Pre: True
		Post: return list of reactions excluding transporters'''
		return list(filter(lambda s:'_tx' not in s, self.reactions))
	
	def GetAllTx(self):
		'''Pre: True
		Post: return list of transporters excluding biomass transporters'''
		return Set.Complement(list(filter(lambda s:'_tx' in s, self.reactions)),list(filter(lambda s:'_bm_tx' in s, self.reactions)))
	
	def MediaTx(self):
		'''Pre: True
		Post: return list of media transporters'''
		return list(filter(lambda s:'_mm_tx' in s, self.reactions))

	def AAMediaTx(self):
		'''Pre: True
		Post: return list of amion acids media transporters'''
		return list(filter(lambda s:'_AA_mm_tx' in s, self.reactions))

	def GetAllSubstrates(self):
		'''Pre: True
		Post: return list of substrates in media'''
		rv = []
		for tx in self.MediaTx():
			subs = tx.split('_')[0]
			if subs not in rv:
				rv.append(subs)
		return rv

	def GetAllAASubstrates(self):
		'''Pre: True 
		Post: return list of amino acids in media'''
		rv = []
		for tx in self.AAMediaTx():
			subs = tx.split('_')[0]
			if subs not in rv:
				rv.append(subs)
		return rv

	def GetAllByProd(self):
		'''Pre: True
		Post: return list of by-products'''
		rv = []
		for tx in list(filter(lambda s:"_bp_tx" in s,self.reactions)):
			subs = tx.split('_')[0]
			if subs not in rv:
				rv.append(subs)
		return rv

	def CompToTx(self, Compound):
		'''Pre: True
		Post: return list of transporters for the compound'''
		if not Compound.startswith("x_"):
			Compound = "x_" + Compound
		txs = [tx for tx in self.smx.InvolvedWith(Compound)]
		if Compound == "x_GLC":
			txs.append('glycogen_bm_tx')
		return txs

	def CompToMediaTx(self,Compound):
		'''Pre: True
		Post: return list of import transporters for the compound'''
		return [tx for tx in self.CompToTx(Compound) if not "_bm_" in tx and not "_bp_" in tx]
	
	def LPBound(self,lp):
		'''Pre: lp
		Post: return lp with upper bound on substrate uptake''' 

		for Compound in media.keys():
			txs = self.CompToMediaTx(Compound)
			if len(txs)==1.0:
				lp.SetFluxBounds({txs[0]:(0.0,media[Compound])})
			if len(txs)>1.0:
				lp.SetSumFluxConstraint(txs,media[Compound],"Total"+Compound)
				lp.SetRowBounds("Total"+Compound,0,media[Compound])
		return lp

	def NetFlux(self,Compound,sol):
		"""Pre:True
		Post: Total substrate uptake flux"""
		transporters = self.CompToTx(Compound)
		influx = 0.0
		for tx in transporters:
			if tx in sol.keys():
				influx += sol[tx]
		return influx 

	def TotalUptakeFlux(self,Compound,sol):
		"""Pre:True
		Post: Total substrate uptake flux"""
		transporters = self.CompToMediaTx(Compound)
		influx = 0.0
		for tx in transporters:
			if tx in sol.keys():
				influx += sol[tx]
		return influx 


	def ElementalCont(self,db,sol):
		'''Pre: database, lp solution
		Post: return carbon and nitrogen contribution by each substrate in the lp solution'''
		rv = {}
		byprod = self.GetAllByProd()
		substrates = self.GetAllSubstrates()
		for subs in substrates+byprod:
			emp = {}
			if subs in db.dbs['Compound'].keys(): 
				if 'C' in db[subs].EmpForm.keys():
					emp['C'] = self.NetFlux(subs,sol)*db[subs].EmpForm['C']
				if 'N' in db[subs].EmpForm.keys():
					emp['N'] = self.NetFlux(subs,sol)*db[subs].EmpForm['N']

			else:
				if subs in empform.keys(): 
					if 'C' in empform[subs].keys():
						emp['C'] = self.NetFlux(subs,sol)*empform[subs]['C']
					if 'N' in empform[subs].keys():
						emp['N'] = self.NetFlux(subs,sol)*empform[subs]['N']
			rv[subs]=emp		
		return rv
	def TotalElementalCont(self,db,sol):
		'''Pre: True
		Post: return total carbon and nitrogen in the lp solution'''
		elecont = self.ElementalCont(db,sol)
		totalCin = 0.0
		totalbpC = 0.0
		totalNin = 0.0
		totalbpN = 0.0
		
		for subs in elecont.keys():
			if 'C' in elecont[subs].keys(): 
				if elecont[subs]['C']< 0.0:
					totalbpC+= abs(elecont[subs]['C'])
				else:
					totalCin += elecont[subs]['C']
			if 'N' in elecont[subs].keys():
				if elecont[subs]['N']< 0.0:
					totalbpN+= abs(elecont[subs]['N'])
				else:
					totalNin += elecont[subs]['N']

		totalbmN = totalNin-totalbpN
		totalbmC = totalCin-totalbpC
		print ('Total C and N input', totalCin, totalNin)		
		print ('Total C and N in biomass', totalbmC, totalbmN)
		print ('total C and N excretion',totalbpC,totalbpN)
		return totalCin,totalNin,totalbmC,totalbmN

	def PercentElementalCont(self,db,sol,Input='Total'):
		'''Pre: True
		Post: return percentage carbon and nitrogen contribution by each substrate in the lp solution'''
		rv = {}
		elecont = self.ElementalCont(db,sol)
		totalCin,totalNin,totalbmC,totalbmN = self.TotalElementalCont(db,sol)
		for subs in elecont.keys():
			pemp = {}
			if 'C' in elecont[subs]:
				if Input == 'Total':
					pemp['C'] = (elecont[subs]['C']/totalCin)*100
				if Input == 'Biomass':
					pemp['C'] = (elecont[subs]['C']/totalbmC)*100
			else:
				pemp['C'] = 0.0

			if 'N' in elecont[subs]:
				if Input == 'Total':
					pemp['N'] = (elecont[subs]['N']/totalNin)*100
				if Input == 'Biomass':
					pemp['N'] = (elecont[subs]['N']/totalbmN)*100
			else:
				pemp['N'] = 0.0

			if pemp['N'] != 0.0 or pemp['C'] != 0.0:
				rv[subs] = pemp
		return rv


def EssRxn(m,mu=SepiBiomass.DefaultMu,mediabound=True): # To generate essential reaction result in section Model simulation on defined rich media
	'''pre:True
	post: list of reactions essential for biomass production'''
	rv = []	
	mr = ModelResults(m)
	lp = BuildLP.BiomassLP(m,mu)
	if mediabound == True:
		lp = mr.LPBound(lp)
	for r in Set.Complement(mr.Reactions(),['ATPASE-RXN']):
		lp.SetFixedFlux({r:0})
		lp.Solve(False)
		if lp.GetStatusMsg() != 'optimal':
			rv.append(r)
		lp.ClearFluxConstraint(r)
	print ("No. of essential reactions",len(rv))	
	return rv

def ATPSol(m,mu=SepiBiomass.DefaultMu,mediabound=True): #To generate result in Energy metabolism section 
	"""Pre:True
	Post: lp sol for ATP synthesis"""
	print ("***********Generating results for Energy metabolism*************")
	mr = ModelResults(m)
	lp = BuildLP.BuildLP(m)
	AMaint = BuildLP.MaintCost(mu)
	lp.SetFixedFlux({'ATPASE-RXN':AMaint})
	if mediabound == True:
		lp = mr.LPBound(lp)
	lp.Solve()
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
		print ("Obj val",lp.GetObjVal(),"\n","no. of reacs",len(sol.keys()),"\n","Net import of aa",list(filter(lambda s:"_AA_mm_tx" in s,sol.keys())))
	return sol

def LPSol(m,mu=SepiBiomass.DefaultMu,mediabound=True): # To generate results in section Model simulation on defined rich media
	"""Pre:True
	Post: lp sol for biomass synthesis"""
	#print ("***********Generating results for Biomass production including maintenance cost*************")
	rv = []
	mr = ModelResults(m)
	lp = BuildLP.BiomassLP(m,mu)
	if mediabound == True:
		lp = mr.LPBound(lp)
	lp.Solve()
	sol = lp.GetPrimSol()
	#print (media,"Obj val",lp.GetObjVal(),"\n","no. of reactions",len(sol.keys()))
	return sol
	
def NetTxFlux(m,sol=None,mu=SepiBiomass.DefaultMu,mediabound=True):
	"""Pre:True
	Post: net flux in media transporters"""
	mr = ModelResults(m)
	rv = {}
	aa = []
	if sol==None:
		sol = LPSol(m,mu,mediabound)
	substrates = mr.GetAllSubstrates()
	for subs in substrates:
		netflux = mr.NetFlux(subs,sol)
		if netflux == 0:
			aa.append(subs)
		else:
			rv[mr.CompToMediaTx(subs)[0]] = netflux
			if len(Set.Intersect(sol.keys(),mr.CompToMediaTx(subs)))>1.0:
				print (subs,"imported from more than 1 transporters")
	print ("####No net import#####",aa)
	return rv

	
def PercentCont(m,db,sol,Input='Total'): #To generate results, C and N contribution by each substrate, in section Model simulation on defined rich media 
	'''Pre: model, database, lp solution
	Post: return percentage carbon and nitrogen contribution by each substrate in the lp solution'''
	print ("***********Generating results for C and N contribution*************")
	mr = ModelResults(m)
	pec = mr.PercentElementalCont(db,sol,Input)
	return pec


def DoAll(m,db,mu=SepiBiomass.DefaultMu,mediabound=True):
	'''Pre: model
	Post: Generates all results on MHHW media'''

	ess = EssRxn(m,mu) #uncomment to get essential reactions
	atpsol = ATPSol(m,mu,mediabound) #uncomment to get results on energy metabolism
	lpsol = LPSol(m,mu,mediabound)
	nettxflux = NetTxFlux(m,lpsol,mu,mediabound)
	percent_cont = PercentCont(m,db,lpsol)
	percent_netcont = PercentCont(m,db,{k: v for (k, v) in nettxflux.items() if v >= 0})
	return ess,atpsol,lpsol,nettxflux,percent_cont,percent_netcont


##To get all results for MHHW media in the manuscript 	
##ess,atpsol,lpsol,nettxflux,percent_cont,percent_netcont = ModelResults.DoAll(m,db)



