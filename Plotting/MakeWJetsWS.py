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
    #sel_base_el = 'ph_n==1 && el_n==1'
    sel_base_el = 'ph_n==2 && el_n>=1'
    #sel_base_el = '1'

    #eta_cuts = ['EB', 'EE']
    eta_cuts = ['EB']

    workspaces_to_save = {}

    selections = { 'base'    : { 
                               # 'mu' : {'selection' : sel_base_mu }, 
                                'el' : { 'selection' : sel_base_el }, 
                               },
                   #'jetVeto' : { 'mu' : {'selection' : sel_jetveto_mu }, 
                   #              'el' : { 'selection' : sel_jetveto_el } ,
                   #            },
                 }

    elefake            = ROOT.RooWorkspace( 'elefake' )

    make_efake( sampManElG, 'Z+jets', sel_base_el,'EB', 'ph_phi', suffix='OvLrm', closure=False, workspace=elefake)

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
    ph_selection_sr    = 'met_pt>40 && ph_passEleVeto[1]==1' #'%s==1' %defs.get_phid_selection('medium')
    ph_selection_B   = 'met_pt<40 && ph_passEleVeto[1]==0' #'%s==1'%defs.get_phid_selection( num_var, _var )
    #'el_n==1 && ph_n==1 && ph_passEleVeto[0] == 1 && met_pt<40 &&ph_pt[0]>160'
    ph_selection_A   = 'met_pt<40 && ph_passEleVeto[1]==1' #'%s==1' %defs.get_phid_selection( num_var )
    ph_selection_D = 'met_pt>40 && ph_passEleVeto[1]==0'#'%s==1' %defs.get_phid_selection( _var )

    full_sel_D = ' && '.join( [sel_base, ph_selection_D, ] )
    full_sel_A   = ' && '.join( [sel_base, ph_selection_A,] )
    full_sel_B   = ' && '.join( [sel_base, ph_selection_B,] )
    full_sel_sr    = ' && '.join( [sel_base, ph_selection_sr,] )

    label_D = 'd_%s_' %suffix
    label_A   = 'a_%s_' %suffix
    label_B   = 'b_%s_' %suffix
    label_S    = 's_%s_' %suffix

    if workspace is None :
        ws = ROOT.RooWorkspace( 'ws') 
    else :
        ws = workspace

    #---------------------------------------
    # draw the histograms
    #---------------------------------------
    binning = (200,-5,5)
    hist_D   = clone_sample_and_draw( sampMan, sample, plot_var , full_sel_D, binning )
    print hist_D
    hist_A   = clone_sample_and_draw( sampMan, sample, plot_var , full_sel_A  , binning )
    hist_B   = clone_sample_and_draw( sampMan, sample, plot_var , full_sel_B  , binning )
    hist_sr  = clone_sample_and_draw( sampMan, sample, plot_var , full_sel_sr  , binning )
    c1  = ROOT.TCanvas('c1','c1')
    hist_A.Draw()
    c1.SaveAs("hist"+label_A+plot_var+".pdf","pdf")
    hist_B.Draw()
    c1.SaveAs("hist"+label_B+plot_var+".pdf","pdf")
    hist_sr.Draw()
    c1.SaveAs("hist"+label_S+plot_var+".pdf","pdf")
    hist_D.Draw()
    c1.SaveAs("hist"+label_D+plot_var+".pdf","pdf")
    zones(hist_A,hist_B,hist_sr,hist_D,suffix+"_"+plot_var)


    print  "Region A: ", hist_A.Integral(80,100)
    print  "Region B: ", hist_B.Integral(80,100)
    print  "full range:"
    print  "Region A: ", hist_A.Integral()
    print  "Region B: ", hist_B.Integral()
    print  "Region Signal: ", hist_sr.Integral()
    print  "Region D: ", hist_D.Integral()


def make_efake2( sampMan, sample, sel_base, eta_cut, plot_var, suffix='', closure=False, workspace=None) :

	#only single plot
    #---------------------------------------
    # Get the base selection for each region
    #---------------------------------------
    ph_selection_A   = 'met_pt<40 && ph_passEleVeto[0]==1' #'%s==1' %defs.get_phid_selection( num_var )
    full_sel_A   = ' && '.join( [sel_base, ph_selection_A,] )
    label_A   = 'a_%s_' %suffix

    #---------------------------------------
    # draw the histograms
    #---------------------------------------
    binning = (200,-5,5)
    hist_A   = clone_sample_and_draw( sampMan, sample, plot_var , full_sel_A  , binning )
    c1  = ROOT.TCanvas('c1','c1')
    hist_A.Draw()
    c1.SaveAs("hist"+label_A+plot_var+".pdf","pdf")

    print  "bin 80 -100: ", hist_A.Integral(80,100)
    print  "full integral: ", hist_A.Integral()



def clone_sample_and_draw( sampMan, samp, var, sel, binning ) :

    newSamp = sampMan.clone_sample( oldname=samp, newname=samp+str(uuid.uuid4()), temporary=True ) 
    sampMan.create_hist( newSamp, var, sel, binning )
    return newSamp.hist

def zones(h1,h3,h2,h4,suffix): 
   c1 = ROOT.TCanvas("c2","multipads",900,700)
   ROOT.gStyle.SetOptStat(0)
   c1.Divide(2,2,0,0)

   c1.cd(1)
   ROOT.gPad.SetTickx(2)
   h1.Draw()

   c1.cd(2)
   ROOT.gPad.SetTickx(2)
   ROOT.gPad.SetTicky(2)
   h2.GetYaxis().SetLabelOffset(0.01)
   h2.Draw()

   c1.cd(3)
   h3.Draw()

   c1.cd(4)
   ROOT.gPad.SetTicky(2)
   h4.Draw()
   c1.SaveAs("hist_multi_"+suffix+".pdf")

main()

    




    
    




