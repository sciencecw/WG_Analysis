import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import uuid
import time
import math
import os
from uncertainties import ufloat
import re
from collections import OrderedDict
import pickle

from SampleManager import SampleManager
from argparse import ArgumentParser
parser = ArgumentParser()

parser = ArgumentParser()
#parser.add_argument('--baseDirNoFilt' , default=None,           dest='baseDirNoFilt',  required=True, help='Path to base directory with no filter')
parser.add_argument('--baseDirLepFilt', default=None,           dest='baseDirLepFilt', required=True, help='Path to base directory with lepton filter')
parser.add_argument('--samplesConf'   , default=None,           dest='samplesConf',    required=True, help='Sample configuration' )
parser.add_argument('--outputFile'   , default=None,           dest='outputFile',    required=False, help='Path to output file' )

options = parser.parse_args()


_TREENAME = 'tupel/EventTree'
_FILENAME = 'tree.root'
_XSFILE   = 'cross_sections/photon15.py'
_LUMI     = 36200

rand = ROOT.TRandom3()
rand.SetSeed( int( time.time() ) )

def main() :

    #sampManNoFilt  = SampleManager( options.baseDir, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManLepFilt = SampleManager( options.baseDirLepFilt, _TREENAME, filename=_FILENAME )

    if options.samplesConf is not None :

        #sampManNoFilt.ReadSamples( options.samplesConf )
        sampManLepFilt.ReadSamples( options.samplesConf )
    else :
        print 'Must provide a sample configuration.  Exiting'
        return


    mode_cuts = { 'electron' : 'isWElDecay == 1', 'muon' : 'isWMuDecay == 1', 'tau' : 'isWTauDecay == 1' , 'tauel' : 'isWTauDecay == 1 && isWTauElDecay == 1'  , 'taumu' : 'isWTauDecay == 1 && isWTauMuDecay == 1' }
    cuts_nofilt = OrderedDict()

    #cuts_nofilt  = { 'electron' : OrderedDict(), 'muon' : OrderedDict() }
    cuts_lepfilt = { 'electron' : OrderedDict(), 'muon' : OrderedDict() }

    cuts_lepfilt['electron']['pass trigger'] = 'passTrig_HLT_Ele27_eta2p1_WPTight_Gsf'
    cuts_lepfilt['muon']    ['pass trigger'] = '(passTrig_HLT_IsoTkMu24 || passTrig_HLT_IsoMu24 )'

    cuts_lepfilt['electron']['At least one medium electron'] = 'el_n > 0'
    cuts_lepfilt['muon']    ['At least one tight muon']      = 'mu_n > 0'

    cuts_lepfilt['electron']['Electron trigger match'] = ' Sum$( el_hasTrigMatch ) > 0 && Sum$( el_pt > 30 ) > 0'
    cuts_lepfilt['muon']    ['Muon trigger match']      = ' Sum$( mu_hasTrigMatch ) > 0 && Sum$( mu_pt > 25 ) > 0'

    cuts_lepfilt['electron']['Second Lepton Veto'] = 'Sum$( el_pt > 10 ) == 1'
    cuts_lepfilt['muon']    ['Second Lepton Veto'] = 'Sum$( mu_pt > 10 ) == 1'

    cuts_lepfilt['electron']['At least one Medium Photon'] = 'Sum$( ph_passMedium ) > 0'
    cuts_lepfilt['muon']    ['At least one Medium Photon'] = 'Sum$( ph_passMedium ) > 0'

    cuts_lepfilt['electron']['At least one Photon with pt $>$ 30'] = 'Sum$( ph_pt > 30 ) > 0'
    cuts_lepfilt['muon']    ['At least one Photon with pt $>$ 30'] = 'Sum$( ph_pt > 30 ) > 0'

    cut_efficiencies = {}
    selected_samples = []

    for samp in sampManLepFilt.get_samples() :

        res = re.match( 'ResonanceMass(\d+)_width(\d+)', samp.name )

        print samp.name

        if res is not None :

            width = int(res.group(2))
            mass  = int(res.group(1))

            cut_efficiencies.setdefault( width, {})
            cut_efficiencies[width].setdefault( mass, {} )
            selected_samples.append( samp.name )


            eff_ele = get_efficiencies( sampManLepFilt, samp.name, mode_cuts, cuts_lepfilt['electron'] )
            eff_mu  = get_efficiencies( sampManLepFilt, samp.name, mode_cuts, cuts_lepfilt['muon'] )

            cut_efficiencies[width][mass] = {'muon' : eff_mu, 'electron' : eff_ele}

    if options.outputFile is not None :

        dirpath = os.path.dirname( options.outputFile )

        if not os.path.isdir( dirpath ) :
            os.makedirs( dirpath )

        ofile = open( options.outputFile, 'w' )

        pickle.dump( cut_efficiencies, ofile )

        ofile.close()


def get_efficiencies( sampManLepFilt, samp_name, mode_cuts, cuts ) :

    binning = ( 100, 0, 100 )
    cut_vals = {'counts' : {}, 'eff' : {}, 'ftau' : {}, 'effwithtau' : {}, 'effCorrftau' : {}}

    for modename, modec in mode_cuts.iteritems() :

        cut_vals['counts'].setdefault( modename, OrderedDict() )
        cut_vals['eff'].setdefault( modename, OrderedDict() )

        sampManLepFilt.create_hist( samp_name, 'ph_n', modec, binning )
        integral = sampManLepFilt.get_samples( name=samp_name )[0].hist.Integral()
        cut_vals['counts'][modename]['denominator'] = integral

        cut_str = ''
        for cut_name, cut in cuts.iteritems() :

            cut_str += ' && %s '%cut

            sampManLepFilt.create_hist( samp_name, 'ph_n',  modec + cut_str, binning )
            npass_v = sampManLepFilt.get_samples( name=samp_name )[0].hist.Integral()
            cut_vals['counts'][modename][cut_name] = npass_v

            tot = cut_vals['counts'][modename]['denominator']
            nfail_v = tot - npass_v

            npass = ufloat( npass_v, math.sqrt(npass_v ) )
            nfail = ufloat( nfail_v, math.sqrt(nfail_v ) )

            cut_vals['eff'][modename][cut_name] = npass / ( npass + nfail )

    # now make the tau fraction
    for modename, mode_dic in cut_vals['counts'].iteritems() :

        if modename.count('tau') > 0 :
            continue

        cut_vals['effwithtau'].setdefault( modename, OrderedDict() )
        cut_vals['ftau'].setdefault( modename, OrderedDict() )
        cut_vals['effCorrftau'].setdefault( modename, OrderedDict() )

        if modename == 'muon' :
            taumode = 'taumu'
        if modename == 'electron' :
            taumode = 'tauel'

        tau_total_v = cut_vals['counts'][taumode]['denominator']
        lep_total_v = mode_dic['denominator']

        for cut_name, cut_res in mode_dic.iteritems() :
            if cut_name == 'denominator' :
                continue

            pass_tau_v = cut_vals['counts'][taumode][cut_name]
            pass_lep_v = cut_res

            fail_tau_v = tau_total_v - pass_tau_v
            fail_lep_v = lep_total_v - pass_lep_v

            pass_tau = ufloat( pass_tau_v, math.sqrt( pass_tau_v ) ) 
            fail_tau = ufloat( fail_tau_v, math.sqrt( fail_tau_v ) ) 
            pass_lep = ufloat( pass_lep_v, math.sqrt( pass_lep_v ) ) 
            fail_lep = ufloat( fail_lep_v, math.sqrt( fail_lep_v ) ) 

            effwithtau = ( pass_tau + pass_lep ) / ( pass_tau + pass_lep + fail_tau + fail_lep )
            ftau = ( pass_tau + fail_tau ) / (  pass_tau + pass_lep + fail_tau + fail_lep  )

            effCorrftau = effwithtau / ( 1 - ftau )

            cut_vals['ftau'][modename][cut_name] = ftau
            cut_vals['effwithtau'][modename][cut_name] = effwithtau
            cut_vals['effCorrftau'][modename][cut_name] = effCorrftau


    # now make the totals
    cut_vals['counts']['total'] = {}
    cut_vals['eff']['total'] = {}

    cut_vals['counts']['total'].setdefault('denominator', 0 )

    for modename in mode_cuts.keys() :

        cut_vals['counts']['total']['denominator'] = cut_vals['counts']['total']['denominator'] + cut_vals['counts'][modename]['denominator']

        for cutname in cuts.keys() :

            cut_vals['counts']['total'].setdefault( cutname, 0 )
            cut_vals['counts']['total'][cutname] = cut_vals['counts']['total'][cutname] + cut_vals['counts'][modename][cutname]


    #now calculate efficiencies

    total = cut_vals['counts']['total']['denominator']

    for cutname, cutv in cut_vals['counts']['total'].iteritems() :
        npass_v = cutv
        nfail_v = total - npass_v

        npass = ufloat( npass_v, math.sqrt( npass_v ) )
        nfail = ufloat( nfail_v, math.sqrt( nfail_v ) )

        cut_vals['eff']['total'][cutname] = npass / ( npass + nfail )


    return cut_vals







main()


