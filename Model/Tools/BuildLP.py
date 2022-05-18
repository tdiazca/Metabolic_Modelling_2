from importlib import reload
import SepiBiomass
reload(SepiBiomass)



#
##
### Growth rates are not a fundamental property of the model
### So need to be included in Analysis, not here.
##
#


def MaintCost(mu=1,GAM=SepiBiomass.GAM, NGAM=SepiBiomass.NGAM ):
    """ we don't care what the media was, just mu,GAM and NGAM.
        There ias quite a strong argument for putting this in SepiBiomass
        as calculating the fd depends on it (and nothing else) """

    """ Growth associated and non-growth associated maintenance cost """
    return mu*GAM+NGAM

  
def BuildLP(m):
    """ build lp with objective as minimisation of total flux """
    
    lp = m.GetLP() 
    lp.SetObjective(m.sm.cnames)   # minimise total flux

    return lp


def BiomassLP(m, mu=SepiBiomass.DefaultMu, GAM=SepiBiomass.GAM, NGAM=SepiBiomass.NGAM):
    """" Again, we don't care what the media was - what if we want to repeat the investigation
         under other conditions/media? Modules in Model/Tools should be *general*, code relating
         to  *specific* conditions goes in analysis (e.g. mapping media<>mu values)
         Removed fd from  the args - how should we process it if specified ?"""

    fd = SepiBiomass.FluxDic(m, mu)
    fd["ATPASE-RXN"] = MaintCost(mu,GAM, NGAM)

    lp = BuildLP(m)
    lp.SetFixedFlux(fd)

    return lp


def ClosedLP(m):

    """ lp object contraining flux through tx reacs to zero """

    lp = BuildLP(m)
    
    for tx in filter(lambda s: "_tx" in s, m.sm.cnames): # all tx in model
        lp.SetFixedFlux({tx:0.0}) # set flux to 0.0

    return lp

def BlockAAUptakeLP(m,Biomass=False, mu=SepiBiomass.DefaultMu):
	"""Pre:True
	Post: lp object, constrained  AA media transporter fluxes to 0  """
	if Biomass == False:
		lp = BuildLP(m)
	else:
		lp = BiomassLP(m,mu)
	for tx in filter(lambda s:"_AA_mm_tx" in s, m.sm.cnames):
		lp.SetFixedFlux({tx:0.0})
	return lp

