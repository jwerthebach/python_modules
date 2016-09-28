#Import the icecube packages into python :
from icecube import icetray, dataclasses, dataio, common_variables, simclasses, phys_services
from icecube.common_variables import direct_hits, hit_multiplicity, hit_statistics, track_characteristics
from icecube.icetray import I3Units


def CheckEventType(frame, Bundle_Leading, Bundle_EinSum):
    """Applying cuts to set the event type.

    
    EventTypes:
    -1 -- Uncharacterized
     0 -- Background
     1 -- HEMuon
     2 -- Balloon
    """

    try:
        if Bundle_Leading/Bundle_EinSum > 0.5:
            if frame['HitStatisticsValues'].q_max_doms/frame['HitStatisticsValues'].q_tot_pulses > 0.2:
                return 2
            else:
                return 1
        else:
            return 0
    
    except:
        return -1



def BundleMultiplicity(frame):
    """Generates Tomasz bundle parameters."""

    MMCTrackList = frame['MMCTrackList']
    bundle = 0
    MaxMult = 200

    E_in = []
    E_lost = []
    E_prod = []
    zen = []
    sum_in= 0.0
    sum_prod = 0.0
    maxE = -1
    secE = -1
    trdE = -1
    prodE = -1
    for p in MMCTrackList:
        pdg = p.GetI3Particle().pdg_encoding
        if pdg in (-13,13):
            if(p.Ei>0):
                bundle+=1
                E_in.append(p.Ei)
                E_lost.append(p.Elost)
                E_prod.append(p.GetI3Particle().energy)
                zen.append(p.GetI3Particle().dir.zenith)
                sum_in += p.Ei
                sum_prod += p.GetI3Particle().energy
                if p.Ei>maxE:
                    trdE = secE
                    secE = maxE
                    maxE = p.Ei
                    prodE = p.GetI3Particle().energy
                elif p.Ei>secE:
                    trdE = secE
                    secE = p.Ei
                elif p.Ei>trdE:
                    trdE = p.Ei


    if(bundle>MaxMult):
        E_in = [0]
        E_lost = [0]
        E_prod = [0]
        zen = [-1]
    
    return bundle, E_in, E_lost, E_prod, zen, sum_in, sum_prod, maxE, secE, trdE, prodE




class TomaszReco(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('RecoName',                 # name
                          'RecoName',                 # doc
                          'TomaszAnalysis')           # default

    def Configure(self):
        self.RecoName = self.GetParameter('RecoName')
        return True
    
    def DAQ(self, frame):
        self.PushFrame(frame)
    
    def Physics(self, frame):

        bundle, E_in, E_lost, E_prod, zen, sum_in, sum_prod, maxE, secE, trdE, prodE = BundleMultiplicity(frame)

        I3_double_container = dataclasses.I3MapStringDouble()


        I3_double_container['MC_Bundle_Mult']          = float( bundle )
        I3_double_container['MC_Bundle_EinSum']        = float( sum_in )
        I3_double_container['MC_Bundle_EprodSum']      = float( sum_prod )
        I3_double_container['MC_Bundle_Leading']       = float( maxE )
        I3_double_container['MC_Bundle_LeadingSurf']   = float( prodE )
        I3_double_container['MC_Bundle_SecondLeading'] = float( secE )
        I3_double_container['MC_Bundle_ThirdLeading']  = float( trdE )
        I3_double_container['MC_Bundle_Type']          = float( CheckEventType(frame, maxE, sum_in) )

        frame['Tomasz_MC_Bundle_Ein']    = dataclasses.I3VectorFloat( E_in )
        frame['Tomasz_MC_Bundle_Elost']  = dataclasses.I3VectorFloat( E_lost )
        frame['Tomasz_MC_Bundle_Eprod']  = dataclasses.I3VectorFloat( E_prod )
        frame['Tomasz_MC_Bundle_Zenith'] = dataclasses.I3VectorFloat( zen )


        frame[self.RecoName] = I3_double_container
        
        self.PushFrame(frame)


