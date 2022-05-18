from ScrumPy.Data import DataSets
from ScrumPy.Util import Sci
import SepiBiomass
import BuildLP
import Media
media = Media.media

AAs =[
    'x_L-ALPHA-ALANINE', 'x_ARG', 'x_ASN', 'x_L-ASPARTATE', 'x_CYS',
    'x_GLN', 'x_GLT', 'x_GLY', 'x_HIS', 'x_ILE', 'x_LEU',   'x_LYS',
    'x_MET', 'x_PHE', 'x_PRO', 'x_SER', 'x_THR', 'x_TRP', 'x_TYR', 'x_VAL'
]

def AABal(m, sol):
    """ pre: sol is a flux dict of m """

    rv = {}

    for aa in AAs:
        In = Out = 0
        txs = m.smx.InvolvedWith(aa)
        
        for tx in txs:
            if tx in sol:
                v = sol[tx]
        
                if v >0:
                    In += v
                else:
                    Out += v
        Rat = abs(In/Out)
                
        rv[aa] = In,Out,Rat

    return rv

def CompToTx(m, Compound):

    if not Compound.startswith("x_"):
        Compound = "x_" + Compound
        
    txs = [tx for tx in m.smx.InvolvedWith(Compound) if not "_bm_" in tx and not "_bp_" in tx] # exclude biomass and by-product "tansporters"
    return txs

def BlockDict(m, Compound):  
    txs = CompToTx(m, Compound)
    return dict(zip(txs, [0]*len(txs)))


def SetMediaConstraints(m, lp, MediaDic = Media.media):

    for Compound in MediaDic:
        txs = CompToTx(m, Compound)
        
        if len(txs) == 1.0:
            lp.SetFluxBounds({txs[0]:(0.0,MediaDic[Compound])})
            
        elif len(txs) > 1.0:
            lp.SetSumFluxConstraint(txs,media[Compound],"Total"+Compound)
            lp.SetRowBounds("Total"+Compound,0,media[Compound])
            
        else:
            print("no tx founds for", Compound) # defensive coding
            # we shouldn't encounter this but if we do we need to know


def BuildLPBound(m,mu, MediaDic = Media.media):
    '''Pre: lp and substrate
    Post: return lp with upper bound on substrate uptake'''
    
    lp = BuildLP.BiomassLP(m,mu)
    SetMediaConstraints(m, lp, MediaDic) 
  
    return lp


def ModelImpact(m):

    rv = {}

    ds = DataSets.DataSet()
    ds.ReadFile("../Analysis/GrowthImpacts.tsv")

    lp = m.GetLP()
    lp.SetObjective(m.sm.cnames)

    for aa in ds.rnames:
        fd = SepiBiomass.FluxDic(m, ds[aa,"mu"])
        fd["ATPASE-RXN"] = BuildLP.MaintCost(ds[aa,"mu"])
        fd.update(BlockDict(m, aa))
        lp.SetFixedFlux(fd)
        lp.Solve()
        rv[aa] = lp.GetPrimSol()
        lp.ClearFluxConstraints(fd)
    return rv


def ModelImpactBoundedScaled(m):

    rv = {}

    ds = DataSets.DataSet()
    ds.ReadFile("../Analysis/GrowthImpacts.tsv")

    for aa in ds.rnames:
        lp = BuildLPBound(m, ds[aa,"mu"])
        lp.SetFixedFlux(BlockDict(m, aa))
        lp.Solve()
        rv[aa] = lp.GetPrimSol()
    return rv


                   
def ModelImpactBoundedUnscaled(m):

    rv = {}
    
    ds = DataSets.DataSet()
    ds.ReadFile("../Analysis/GrowthImpacts.tsv")

    mu = 0.1 # BuildLP.MaintCost(BuildLP.SepiBiomass.DefaultMu)
   
    for aa in ds.rnames:
        lp = BuildLPBound(m, mu) # horrid building new lp at each iteration!
        Block = BlockDict(m, aa)
        lp.SetFixedFlux(Block)
        lp.Solve()
        rv[aa] = lp.GetPrimSol()
        
    return rv


def FluxDS(m,soldic):
    """ pre: soldic is a dictionary of of fluxdics from m. e.g from ModelImpact """

    keys = [k for k in soldic] # because d.keys() is not indexible in py3
    k0 = keys[0]
    rv = m.smx.RateVector(soldic[k0])
    rv.cnames[0] = k0
    
    for k in keys[1:]:
        v= m.smx.RateVector(soldic[k])
        rv.NewCol(v.GetCol(0), k)
    return rv

    
def ModelVExp(m, Remove=["VAL"]):

    Ref = "MHHW"

    rv = DataSets.DataSet()
    rv.ReadFile("../Analysis/GrowthImpacts.tsv")
    rv.NewCol(name="EDist")

    soldic = ModelImpactBoundedScaled(m)
    for r in Remove:
        del soldic[r]
        rv.DelRow(r)
    #soldic[Ref] = soldic["x_MHHW"]
    #del soldic["x_MHHW"]

    FluxMtx = FluxDS(m, soldic)
    FluxMtx.Transpose() # more convenient to iterate over rows than cols
    for aa in rv.rnames:
        rv[aa, "EDist"] = Sci.EucDist(FluxMtx[Ref], FluxMtx[aa])

    return rv
    
    

    

    

        
        
        
