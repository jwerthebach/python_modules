#Import the icecube packages into python :
from icecube import icetray, dataclasses, dataio, common_variables, simclasses, phys_services
from icecube.common_variables import direct_hits, hit_multiplicity, hit_statistics, track_characteristics
from icecube.icetray import I3Units
import numpy as np
from icecube import ddddr

def PatrickEReco(ccasc_100 = 0, coszen_mpe=0, d_slant = 0):
    cascest_atmu = 10**np.exp(0.524+0.214*np.log10(ccasc_100))
    esurf_base = 3.44*cascest_atmu
    esurf_angcorr = 10**(2.92862
                         -22.0136*coszen_mpe
                         +71.4991*coszen_mpe**2
                         -118.7*coszen_mpe**3
                         +96.6*coszen_mpe**4
                         -30.55*coszen_mpe**5
                         +.241823/(1+np.exp(-9.63*(np.log10(d_slant)-3.88)))
                         )
    esurf_atmu = 10**(.553757+.884*np.log10(esurf_base*esurf_angcorr))
    
    return ccasc_100 , coszen_mpe , d_slant , cascest_atmu , esurf_base , esurf_angcorr , esurf_atmu

class GetPatrickReco(icetray.I3ConditionalModule):
    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter('RecoName',                 # name
                          'RecoName',                 # doc
                          'PatrickAnalysis')          # default
                          
        self.AddParameter('D4R_RecoName',             # name
                  'Name of D4R reconstruction',       # doc
                  'PatrickAnalysis_D4R')              # default

    def Configure(self):
        self.D4R_RecoName = self.GetParameter('D4R_RecoName')
        self.RecoName = self.GetParameter('RecoName')
        return True
    
    def DAQ(self, frame):
        self.PushFrame(frame)
    
    def Physics(self, frame):
        
        I3_double_container = dataclasses.I3MapStringDouble()

        try:
            coszen = np.cos( frame['SplineMPE'].dir.zenith )
            ccasc_100 = frame[self.D4R_RecoName +'100'+ 'SplineMPECascadeParams'].cascade_energy
            d_slant = frame[self.D4R_RecoName +'100' + 'SplineMPECascadeParams'].cascade_slant_depth
            
            ccasc_100 , coszen , d_slant , cascest_atmu , esurf_base , esurf_angcorr , esurf_atmu = PatrickEReco(ccasc_100 , coszen, d_slant)

            I3_double_container['IN_ccasc_100'] = ccasc_100
            I3_double_container['IN_coszen']    = coszen
            I3_double_container['IN_d_slant']   = d_slant
            
            I3_double_container['CALC_cascest_atmu']  = cascest_atmu
            I3_double_container['CALC_esurf_base']    = esurf_base
            I3_double_container['CALC_esurf_angcorr'] = esurf_angcorr
            
            I3_double_container['esurf_atmu'] = esurf_atmu
            
            I3_double_container['CUTS_Qtot1000']     = float( frame["HitStatisticsValues"].q_tot_pulses > 1000 )
            I3_double_container['CUTS_CosZen01']     = float( coszen > 0.1 )
            I3_double_container['CUTS_QmaxQtot05']   = float( frame["HitStatisticsValues"].q_max_doms/frame["HitStatisticsValues"].q_tot_pulses < 0.5 )
            I3_double_container['CUTS_Ndom40']       = float( frame[self.D4R_RecoName +'100'+ 'SplineMPEParams'].nDOMs > 40 )
            I3_double_container['CUTS_PeakMedian10'] = float( np.min( [ frame[self.D4R_RecoName +'70'+ 'SplineMPEParams'].peak_energy / frame[self.D4R_RecoName +'70'+ 'SplineMPEParams'].median,
                                                                        frame[self.D4R_RecoName +'100'+ 'SplineMPEParams'].peak_energy / frame[self.D4R_RecoName +'100'+ 'SplineMPEParams'].median,
                                                                        frame[self.D4R_RecoName +'150'+ 'SplineMPEParams'].peak_energy / frame[self.D4R_RecoName +'150'+ 'SplineMPEParams'].median]
                                                                     ) > 10 )
            I3_double_container['CUTS_Median02'] = float( frame[self.D4R_RecoName +'100'+ 'SplineMPEParams'].median > 0.2 )
            I3_double_container['CUTS_ECasc5']   = float( ccasc_100 > 5 )
            I3_double_container['CUTS_ALL']      = np.product( [I3_double_container['CUTS_Qtot1000'],
                                                                I3_double_container['CUTS_CosZen01'],
                                                                I3_double_container['CUTS_QmaxQtot05'],
                                                                I3_double_container['CUTS_Ndom40'],
                                                                I3_double_container['CUTS_PeakMedian10'],
                                                                I3_double_container['CUTS_Median02'],
                                                                I3_double_container['CUTS_ECasc5']]
                                                          )


        except:        
            I3_double_container['IN_ccasc_100'] = 0.0
            I3_double_container['IN_coszen']    = 0.0
            I3_double_container['IN_d_slant']   = 0.0
            
            I3_double_container['CALC_cascest_atmu']  = 0.0
            I3_double_container['CALC_esurf_base']    = 0.0
            I3_double_container['CALC_esurf_angcorr'] = 0.0
            
            I3_double_container['esurf_atmu'] = 0.0
            
            I3_double_container['CUTS_Qtot1000']     = 0.0
            I3_double_container['CUTS_CosZen01']     = 0.0
            I3_double_container['CUTS_QmaxQtot05']   = 0.0
            I3_double_container['CUTS_Ndom40']       = 0.0
            I3_double_container['CUTS_PeakMedian10'] = 0.0
            I3_double_container['CUTS_Median02']     = 0.0
            I3_double_container['CUTS_ECasc5']       = 0.0
            I3_double_container['CUTS_ALL']          = 0.0

        
        frame[self.RecoName] = I3_double_container
        
        self.PushFrame(frame)

# Using SRTHVInIcePulses instead of HVInIcePulses because it is a bit 'cleaner'
@icetray.traysegment
def PatrickReco(tray, Name='PatrickAnalysis',
                  pulse_map="SRTHVInIcePulses",
                  bin_width=50,
                  ref_spline='SplineMPE',
                  keep_D4R=False):
    D4R_Name = Name + '_D4R'
    
    tray.AddModule('I3MuonEnergy', 'Eme70',
       BinWidth      = bin_width,
       InputPulses   = pulse_map,
       MaxImpact     = 70,
       Seed          = ref_spline,
       SaveDomResults = False,
       Prefix = D4R_Name+str(70)
      )
    tray.AddModule('I3MuonEnergy', 'Eme100',
       BinWidth      = bin_width,
       InputPulses   = pulse_map,
       MaxImpact     = 100,
       Seed          = ref_spline,
       SaveDomResults = False,
       Prefix = D4R_Name+str(100)
		  )
    tray.AddModule('I3MuonEnergy', 'Eme150',
       BinWidth      = bin_width,
       InputPulses   = pulse_map,
       MaxImpact     = 150,
       Seed          = ref_spline,
       SaveDomResults = False,
       Prefix = D4R_Name+str(150)
		  )
        
    tray.AddModule(GetPatrickReco,
                    'PatrickReco_OutPut',
                    RecoName=Name,
                    D4R_RecoName=D4R_Name)