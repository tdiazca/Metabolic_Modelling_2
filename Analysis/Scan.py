import numpy
from ScrumPy.Data import DataSets
from importlib import reload
import SepiBiomass
import Media
media = Media.media
import BuildLP
reload(BuildLP)


CarbonSFs = {}

MultOpt = []
reacs = []


class ScanRes(DataSets.DataSet):

    def __init__(self, *args, **kwargs):

        DataSets.DataSet.__init__(self,*args, **kwargs)
    
    
    def Changers(self,lim=1e-7):
        rv = []
        for cname in self.cnames:
            col = self.GetCol(cname)
            if abs(max(col)-min(col))>lim:
                rv.append(cname)
        return rv
     
def CompToTx(m, Compound):

    if not Compound.startswith("x_"):
        Compound = "x_" + Compound
        
    txs = [tx for tx in m.smx.InvolvedWith(Compound) if not "_bm_" in tx and not "_bp_" in tx] # exclude biomass and by-product "tansporters"
    return txs


def BuildLPBound(m,mu):
    '''Pre: lp and substrate
    Post: return lp with upper bound on substrate uptake''' 
    lp = BuildLP.BiomassLP(m,mu)
    for Compound in media.keys():
        txs = CompToTx(m, Compound)
        if len(txs)==1.0:
            lp.SetFluxBounds({txs[0]:(0.0,media[Compound])})
        if len(txs)>1.0:
            lp.SetSumFluxConstraint(txs,media[Compound],"Total"+Compound)
            lp.SetRowBounds("Total"+Compound,0,media[Compound])
    return lp


def ATPScan(m,mu=SepiBiomass.DefaultMu,lo=20.0, hi=40, np=100, ForceATP=True, mediabound=False):
    """ pre: True
	post: returns ATP scan results"""

    rv = ScanRes()
    if mediabound == False:
        lp = BuildLP.BiomassLP(m,mu)
    else:
        lp = BuildLPBound(m,mu)

    if ForceATP:
        SetATP = lambda ATP: lp.SetFixedFlux({'ATPASE-RXN': ATP})
    else:
        SetATP = lambda ATP: lp.SetFluxBounds({'ATPASE-RXN':(ATP,None)})

    ranges = numpy.arange(lo,hi,(hi-lo)/(np-1))
    for ATPLim in ranges:
            SetATP(ATPLim)
            lp.Solve()
            if lp.IsStatusOptimal():
                    sol = lp.GetPrimSol()
                    sol['ObjVal'] = lp.GetObjVal()
                    sol["ATPLim"] = ATPLim
                    rv.UpdateFromDic(sol)
    rv.SetPlotX('ATPLim')
    return rv           

     
