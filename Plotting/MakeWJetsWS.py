import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import re
import os
import uuid
import math
import pickle
import selection_defs as defs
#from uncertainties import ufloat
from FitManager import FitManager
#ROOT.TVirtualFitter.SetMaxIterations( 100000 )
ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls( 100000)

from SampleManager import SampleManager
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('--baseDirMuG',      default=None,           dest='baseDirMuG',         required=False, help='Path to muon base directory')
parser.add_argument('--baseDirElG',      default=None,           dest='baseDirElG',         required=False, help='Path to electron base directory')
parser.add_argument('--outputDir',      default=None,           dest='outputDir',         required=False, help='Output directory to write histograms')
parser.add_argument('--useRooFit',       default=False,    action='store_true',      dest='useRooFit', required=False, help='Make fits using roostats' )
parser.add_argument('--doClosure',       default=False,   action='store_true',       dest='doClosure', required=False, help='make closure tests' )

options = parser.parse_args()

_TREENAME = 'UMDNTuple/EventTree'
_FILENAME = 'tree.root'
_XSFILE   = 'cross_sections/photon15.py'
_LUMI     = 36000
_BASEPATH = '/home/jkunkle/usercode/Plotting/LimitSetting/'
_SAMPCONF = 'Modules/Resonance.py'


#def get_cut_defaults( _var, ieta ) :
#
#    cut_defaults = {'sigmaIEIE' : { 'EB' : ( 0.012, 0.02 ), 'EE' : ( 0.0, 0.0 ) },
#                    'chIso'     : { 'EB' : ( 4, 10 ),       'EE' : (4, 10) },
#                   }
#
#    return cut_defaults[_var][ieta]


ROOT.gROOT.SetBatch(True)
if options.outputDir is not None :
#    ROOT.gROOT.SetBatch(True)
    if not os.path.isdir( options.outputDir ) :
        os.makedirs( options.outputDir )

def main() :

    sampManMuG= SampleManager( options.baseDirMuG, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManElG= SampleManager( options.baseDirElG, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )

    sampManMuG.ReadSamples( _SAMPCONF )
    sampManElG.ReadSamples( _SAMPCONF )

    sampManMuG.outputs = {}
    sampManElG.outputs = {}

    sel_base_mu = 'mu_pt30_n==1 && mu_n==1'
    #sel_base_el = 'el_pt30_n==1 && el_n==1'
    sel_base_el = 'ph_n==1 && el_n==1'

    #eta_cuts = ['EB', 'EE']
    eta_cuts = ['EB']

    workspaces_to_save = {}

    #xmin_m = 160
    #xmax_m = 2000
    #bin_width_m = 20

    #xmin_pt = xmin_m/2
    #if xmin_pt < 50 :
    #    xmin_pt = 50
    #xmax_pt = xmax_m/2
    #bin_width_pt = bin_width_m/2.

    #binning_m = ((xmax_m-xmin_m)/bin_width_m, xmin_m, xmax_m)

    #binning_pt = ( (xmax_pt - xmin_pt )/bin_width_pt, xmin_pt, xmax_pt )

    #xvar_m = ROOT.RooRealVar( 'x_m', 'x_m',xmin_m , xmax_m)

    #xvar_pt = ROOT.RooRealVar( 'x_pt', 'x_pt', xmin_pt, xmax_pt )

    #kine_vars = { #'mt_incl_lepph_z' : { 'var' : 'mt_lep_met_ph'   , 'xvar' : xvar_m  , 'binning' : binning_m},
    #              #'m_incl_lepph_z'  : { 'var' : 'm_lep_met_ph'    , 'xvar' : xvar_m  , 'binning' : binning_m},
    #              ##'mt_rotated'      : { 'var' : 'mt_rotated'      , 'xvar' : xvar_m  , 'binning' : binning_m},
    #              'mt_fulltrans'    : { 'var' : 'mt_res'          , 'xvar' : xvar_m  , 'binning' : binning_m},
    #              #'mt_constrwmass'  : { 'var' : 'recoM_lep_nu_ph' , 'xvar' : xvar_m  , 'binning' : binning_m},
    #              #'ph_pt'           : { 'var' : 'ph_pt[0]'        , 'xvar' : xvar_pt , 'binning' : binning_pt},
    #            }

    selections = { 'base'    : { 
                               # 'mu' : {'selection' : sel_base_mu }, 
                                'el' : { 'selection' : sel_base_el }, 
                               },
                   #'jetVeto' : { 'mu' : {'selection' : sel_jetveto_mu }, 
                   #              'el' : { 'selection' : sel_jetveto_el } ,
                   #            },
                 }

    elefake            = ROOT.RooWorkspace( 'elefake' )

    make_efake( sampManElG, 'Z+jets', sel_base_el,'EB', 'mt_res', suffix='efake', closure=False, workspace=elefake)
    #for seltag, chdic in selections.iteritems() : 

    #    for ch, seldic in chdic.iteritems() : 

    #        for name, vardata in kine_vars.iteritems() :
    # ###   definition: make_wjets_fit( sampMan, sample, sel_base, eta_cut, plot_var, _var, num_var, binning, xvar, suffix='', closure=False, workspace=None)                                 
    #            make_wjets_fit( sampManMuG, 'Data', seldic['selection'], 'EB', vardata['var'], 'chIso', 'sigmaIEIE', vardata['binning'], vardata['xvar'], suffix='wjets_%s_EB_%s_%s' %(ch,name,seltag), closure=False, workspace=wjets)

    #            if options.doClosure :

    #                closure_res_mu = make_wjets_fit( sampManMuG, 'Wjets', seldic['selection'], 'EB', 'mt_lep_met_ph', 'chIso', 'sigmaIEIE', binning_m, xvar_m, suffix='closure_%s_EB_%s' %( ch, seltag ), closure=True )


    if options.outputDir is not None :

        elefake.writeToFile( '%s/%s.root' %( options.outputDir,elefake.GetName() ) )

        for fileid, ws_list in workspaces_to_save.iteritems() :
            for idx, ws in enumerate(ws_list) :
                if idx == 0 :
                    recreate = True
                else  :
                    recreate = False

                ws.writeToFile( '%s/workspace_%s.root' %( options.outputDir, fileid ), recreate )

        for key, can in sampManMuG.outputs.iteritems() :
            can.SaveAs('%s/%s.pdf' %( options.outputDir, key ) )
        for key, can in sampManElG.outputs.iteritems() :
            can.SaveAs('%s/%s.pdf' %( options.outputDir, key ) )


def make_efake( sampMan, sample, sel_base, eta_cut, plot_var, suffix='', closure=False, workspace=None) :

	## el_n==1 && ph_n==1 && ph_passEleVeto[0] == 1 && met_pt<40 && abs(m_lep_ph-m_Z)<30 
	## el_n==1 && ph_n==1 && ph_passEleVeto[0] == 0 && met_pt<40 && abs(m_lep_ph-m_Z)<30 
    #---------------------------------------
    # Get the base selection for each region
    #---------------------------------------
    ph_selection_sr    = 'met_pt>40 && ph_passEleVeto[0]==1' #'%s==1' %defs.get_phid_selection('medium')
    ph_selection_B   = 'met_pt<40 && ph_passEleVeto[0]==0' #'%s==1'%defs.get_phid_selection( num_var, _var )
    #'el_n==1 && ph_n==1 && ph_passEleVeto[0] == 1 && met_pt<40 &&ph_pt[0]>160'
    ph_selection_A   = 'met_pt<40 && ph_passEleVeto[0]==1' #'%s==1' %defs.get_phid_selection( num_var )
    ph_selection_D = 'met_pt>40 && ph_passEleVeto[0]==0'#'%s==1' %defs.get_phid_selection( _var )

    full_sel_D = ' && '.join( [sel_base, ph_selection_D, ] )
    full_sel_A   = ' && '.join( [sel_base, ph_selection_A,] )
    full_sel_B   = ' && '.join( [sel_base, ph_selection_B,] )
    full_sel_sr    = ' && '.join( [sel_base, ph_selection_sr,] )

    label_D = 'D_%s' %suffix
    label_A   = 'A_%s' %suffix
    label_B   = 'B_%s' %suffix
    label_sr    = 'sr_%s' %suffix

    if workspace is None :
        ws = ROOT.RooWorkspace( 'ws') 
    else :
        ws = workspace

    #---------------------------------------
    # draw the histograms
    #---------------------------------------
    binning = (200,0,200)
    hist_D = clone_sample_and_draw( sampMan, sample, 'm_lep_ph', full_sel_D, binning )
    hist_A   = clone_sample_and_draw( sampMan, sample, 'm_lep_ph', full_sel_A  , binning )
    hist_B   = clone_sample_and_draw( sampMan, sample, 'm_lep_ph', full_sel_B  , binning )
    hist_sr     = clone_sample_and_draw( sampMan, sample, 'm_lep_ph', full_sel_sr  , binning )
    print hist_D
    c1  = ROOT.TCanvas('c1','c1')
    hist_A.Draw()
    c1.SaveAs("hist_a.pdf","pdf")
    hist_B.Draw()
    c1.SaveAs("hist_b.pdf","pdf")
    hist_sr.Draw()
    c1.SaveAs("hist_sr.pdf","pdf")
    hist_D.Draw()
    c1.SaveAs("hist_d.pdf","pdf")

    print  "Region A: ", hist_A.Integral(80,100)
    print  "Region B: ", hist_B.Integral(80,100)
    print  "Region Signal: ", hist_sr.Integral()
    print  "Region D: ", hist_D.Integral()



def clone_sample_and_draw( sampMan, samp, var, sel, binning ) :

    newSamp = sampMan.clone_sample( oldname=samp, newname=samp+str(uuid.uuid4()), temporary=True ) 
    sampMan.create_hist( newSamp, var, sel, binning )


main()

    




    
    




