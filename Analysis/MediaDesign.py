##Python source code for ScrumPy3
from importlib import reload
from ScrumPy.Util import Set
import matplotlib.pyplot as plt
import BuildLP
reload(BuildLP)
import ModelResults
reload(ModelResults)
import SepiBiomass
reload(SepiBiomass)


def GetSol(m,lp,mediabound=True,substrate=None):
	'''Pre: m, lp and substrate to be blocked
	Post: return lp solution, obj value, no of reac with blocked substrate transporters '''
	
	if substrate != None:
		transporters = ModelResults.SubTx(m,substrate)
		for tx in transporters:
			lp.SetFixedFlux({tx:0.0})
		lp.Solve(False)
	
	else:
		lp.Solve(False)
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
		return (sol,lp.GetObjVal(),len(sol.keys()))
	else:
		return None

def SolMM(m,Biomass=True,mu=SepiBiomass.DefaultMu):
	'''Pre: True
	Post:  LP Sol for biomass production on inorganic N source'''
	rv = {}
	lp = BuildLP.BlockAAUptakeLP(m,Biomass,mu)
	lp.Solve()
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
	print ("Inputs on minimal media",list(filter(lambda s:'_mm_tx' in s, sol.keys())))
	print ("Obj val and no. of reacs",lp.GetObjVal(),len(sol.keys()))
	return rv


def SingleNSource(m,Biomass=True,mu=SepiBiomass.DefaultMu,Block='AMMONIUM_mm_tx'): #To generate results in section Single nitrogen sources
	'''Pre: True
	Post: returns dict with keys as single N source that can lead to Biomass production'''
	rv = {}
	nf = []
	nftx = []
	mr = ModelResults.ModelResults(m)
	substrates = mr.GetAllAASubstrates()
	lp = BuildLP.BlockAAUptakeLP(m,Biomass,mu)
	lp.SetFixedFlux({Block:0.0})
	
	for subs in substrates:
		transporters = mr.SubTx(subs)
		lp.ClearFluxConstraints(transporters)
		lp.Solve()
		if lp.IsStatusOptimal():
			sol = lp.GetPrimSol()
			rv[subs] = (sol,lp.GetObjVal(),len(sol.keys()))
		else:
			nf.append(subs)
			nftx.extend(transporters)
		for tx in transporters:
			lp.SetFixedFlux({tx:0.0})

	print ('Substrates not utilised as single N source', nf)
	lp.ClearFluxConstraints(nftx)
	lp.Solve()
	if lp.IsStatusOptimal():
		sol = lp.GetPrimSol()
		print ('feasible solution when substrates not utilised as single N source present in combination')
		rv['Suboptimal N sources'] = (sol,lp.GetObjVal(),len(sol.keys()))
	else:
		print ("no feasible solution when", nf, "present in combination either")
	return rv

def SingleSubsRemoval(m,mu=SepiBiomass.DefaultMu,mediabound=True):#To generate results in section Identification of substrate auxotrophies
	'''Pre: True
	Post: return dict with keys as substrate and values as tuple of LP sol, obj value and no. of reacs on individual substrate removal'''
	rv = {}
	aux = []
	mr = ModelResults.ModelResults(m)
	substrates = mr.GetAllSubstrates()
	
	lp = BuildLP.BiomassLP(m,mu)
	if mediabound==True:	
		lp = mr.MediaBound(lp)
	
	lp.Solve(False)
	sol = lp.GetPrimSol()
	rv['MHHW'] = (sol,lp.GetObjVal(),len(sol.keys()))
	
	for subs in substrates:
		lp = BuildLP.BiomassLP(m,mu)
		if mediabound==True:	
			lp = mr.MediaBound(lp)
		for tx in mr.SubTx(subs):
			lp.SetFixedFlux({tx:0.0})
		lp.Solve(False)
		if lp.IsStatusOptimal():
			sol = lp.GetPrimSol()
			rv[subs] = (sol,lp.GetObjVal(),len(sol.keys()))
		else:
			aux.append(subs)
				
		
	print ('Auxotrophy substrates', aux)
	return rv
		

