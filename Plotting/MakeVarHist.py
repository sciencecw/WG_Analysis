import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import re
import numpy as np
import os
import uuid
import math
import pickle
import selection_defs as defs
from itertools import product
from pprint import pprint
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
parser.add_argument('--data',           default=False,          dest='data',          required=False, help='Use data or MC')
parser.add_argument('--useRooFit',       default=False,    action='store_true',      dest='useRooFit', required=False, help='Make fits using roostats' )
parser.add_argument('--doClosure',       default=False,   action='store_true',       dest='doClosure', required=False, help='make closure tests' )

options = parser.parse_args()

_TREENAME = 'UMDNTuple/EventTree'
_FILENAME = 'tree.root'
_XSFILE   = 'cross_sections/photon15.py'
_LUMI     = 36000
_BASEPATH = '/home/jkunkle/usercode/Plotting/LimitSetting/'
_SAMPCONF = 'Modules/Resonance.py'

base = 'ph_n==1 && el_n==1'
iseb = ' ph_IsEB[0]'
baseeta = base+"&&"+iseb
passpix = 'ph_hasPixSeed[0]==0'  #Pixel seed
failpix = 'ph_hasPixSeed[0]==1'
ltmet = 'met_pt<25'
gtmet = 'met_pt>25'
unblind = "ph_hasPixSeed[0]==1 || met_pt<25"
weight = "PUWeight*NLOWeight"


outputdir = 'varhist/'

ROOT.gROOT.SetBatch(True)
if options.outputDir is not None :
#    ROOT.gROOT.SetBatch(True)
    if not os.path.isdir( options.outputDir ) :
        os.makedirs( options.outputDir )

def main() :
#    sampManMuG= SampleManager( options.baseDirMuG, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
#    sampManElG= SampleManager( options.baseDirElG, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )

#    sampManMuG.ReadSamples( _SAMPCONF )
#    sampManElG.ReadSamples( _SAMPCONF )

    #elefake            = ROOT.RooWorkspace( 'elefake' )
    #f1 = ROOT.TFile("%s/output.root"%(options.outputDir),"RECREATE")
    #f1.Write()
    #f1.Close()


    selectionset = []
    etalist = [-3,-2,-1.4,-1,-0.5,0.,0.5,1,1.4,2,3]
    etasellist = [("pheta%i%i" %r, "ph_eta[0]>%g&&ph_eta[0]<=%g" %r)\
                        for r in zip(etalist[:-1],etalist[1:])]
    ptlist = [0,30,50,80,100,200,1000]
    #ptlist = [0,30,50]
    ptsellist = [("phpt%i%i" %r, "ph_pt[0]>%g&&ph_pt[0]<=%g" %r)\
                        for r in zip(ptlist[:-1],ptlist[1:])]
    metsellist = [("metlt25",ltmet), ("metgt25",gtmet)]
    pixsellist = [("passpix",passpix),("failpix",failpix)]
    abcdlist = zip("ABCD",["&&".join((a,b)) for a in (ltmet, gtmet) for b in (passpix,failpix)])
    selectionproduct([abcdlist,abcdlist],True)

def selectionproduct(slist,talks=False):
    #slist2 = list(product(*slist))
    #slist3 = [zip(*x) for x in slist2]
    #slist2 = list(product(*slist))
    slist3 = [zip(*x) for x in product(*slist)]
    slist4 = [("_".join(i[0]), "&&".join(i[1])) for i in slist3]
    if talks:
        pprint(slist)
        print
        #print slist2
        #print 
        pprint(slist3)
        print 
        pprint( slist4)
    return  slist4

def makeplotsset(plotvar, selectionsets, xbins, hist_config,canvas=None):
    ## TODO: check no duplicate in selectionsets
    for filename, selections in selectionsets:
        outtext = makeplots(plotvar,selections, xbinx, hist_config, canvas, filename+".pdf")




def makeplots(plotvar, selection, xbins, hist_config,canvas=None,savename = "temp.pdf"):
    if canvas==None:
        canvas = ROOT.TCanvas("c1","c1")
    canvas.cd() 
    #samples.Draw("m_lep_ph","%s&&%s>%g&&%s<%g" %(basesel,var,ptrange[0],var,ptrange[1]), (200,0,200), { "weight": weight })
    samples.Draw(plotvar, selection, xbins, hist_config)
    canvas.SaveAs(savename)
    outtext = savename +"\n ===== \n"
    result = samples.get_stack_count()
    for line in result: outtext+="%15s %8.3g +/- %5.3g\n" %line
    return outtext



main()

    




    
    




