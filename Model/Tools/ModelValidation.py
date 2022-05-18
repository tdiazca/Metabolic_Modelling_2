
######
########
##########   This module is used to perform theoretical validation of the model
##########  
########        
######                
from importlib import reload
import BuildLP
reload(BuildLP)
import SepiBiomass
reload(SepiBiomass)


## E conservation

def ATPCons(m): # lp object checks m for energy conservation

    """pre: m = the model
       post: if energy is not conseved, returns a solution for
       producing ATP out of nothing"""

    lp = BuildLP.ClosedLP(m)   # constraints tx flux to zero
    lp.SetFixedFlux({"ATPASE-RXN":1}) # set flux through ATPase r to 1
    lp.Solve(False)
    
    if lp.IsStatusOptimal():
        sol=lp.GetPrimSol()
        return sol
    
    else:
        print ("ATP is not produced when medium components are not available")

## RedOx conservation

def NADHCons(m,redox='NADHOxid'): # lp object checks m for RedOx conservation

    """pre: NADHOxid reaction must be included in m (uncomment in top module)
            m = the model
       post: if RedOx balance is not conserved returns solution"""

    lp = BuildLP.ClosedLP(m)   # creates new LP, min flux as objective. Constraints tx flux to zero
    lp.SetFixedFlux({"NADHOxid":1}) # set flux through NADHOxid r to 1
    lp.Solve(False)
    
    if lp.IsStatusOptimal():
        sol=lp.GetPrimSol()
        return sol
    
    else:
        print  ("NAD(P)H is not produced when medium components are not available")


## Mass conservation

def MassCons(m): # lp object checks m for mass conservation

    """pre: m = the modelfilter(lambda s: "_tx" in s, m.sm.cnames)
       post: returns a dictionary with keys as transporters for biomass components
       and values as solution for their production out of 'nothing' if mass is not
       conserved, or value {} if conserved"""

    MassDic={}

    lp = BuildLP.ClosedLP(m)   # constraints tx flux to zero

    for bmtx in list(filter(lambda s: "_bm_tx" in s, m.sm.cnames)):
        lp.ClearFluxConstraint(bmtx)
        lp.SetFixedFlux({bmtx:-1}) # set flux through BM tx to -1
        lp.Solve(False)
        
        sol = lp.GetPrimSol()
        MassDic[bmtx] = sol
        
        lp.SetFixedFlux({bmtx:0.0})

    UnconComp = []

    for comp in MassDic.keys():
        if MassDic[comp] != {}:
            UnconComp.append(comp)

    print ("These compounds can be produced out of 'nothing': ", "\n", UnconComp)

    return MassDic

##Biomass Checks
def CheckUnitFlux(m,reac,neg=False,lp=None):
    """ pre: reac in m.sm.cnames
       post: the lp solution for a unit flux in reac or {} if none exists"""

    if lp ==None:
        lp = BuildLP.BuildLP(m)

    j = -1 if neg else 1
    lp.SetFixedFlux({reac:j})
    lp.Solve()
    lp.ClearFluxConstraint(reac)
    return lp.GetPrimSol()


			
def CheckBMProduction(m,neg=True,lp=None,fd=None):
	"""pre: True
	post: the lp solution for production of each biomass components"""
	rv = {}
	if fd == None:
		bm = list(filter(lambda s:'_bm_tx' in s, m.sm.cnames))
	else:
		fd = SepiBiomass.FluxDic(m)
		bm = list(fd.keys())
	if lp == None:
        	lp = BuildLP.BuildLP(m)
	for reac in bm:
		sol = CheckUnitFlux(m,reac,neg,lp)
		if sol == {}:
			print(reac)
			rv[reac] = 'NF'
		else:
			rv[reac] = sol

	return rv
	


