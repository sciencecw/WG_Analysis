"""
Interactive script to plot data-MC histograms out of a set of trees.
"""

# Parse command-line options
from argparse import ArgumentParser
p = ArgumentParser()
p.add_argument('--baseDir',      default=None,           dest='baseDir',         help='Path to base directory containing all ntuples')
p.add_argument('--baseDirModel',      default=None,           dest='baseDirModel', help='Path to base directory containing all ntuples for the model')
p.add_argument('--fileName',     default='ntuple.root',  dest='fileName',        help='( Default ntuple.root ) Name of files')
p.add_argument('--treeName',     default='events'     ,  dest='treeName',        help='( Default events ) Name tree in root file')
p.add_argument('--treeNameModel',     default='photons'     ,  dest='treeNameModel',help='( Default photons ) Name tree in root file')
p.add_argument('--samplesConf',  default=None,           dest='samplesConf',     help=('Use alternate sample configuration. '
                                                                                       'Must be a python file that implements the configuration '
                                                                                       'in the same manner as in the main() of this script.  If only '
                                                                                       'the file name is given it is assumed to be in the same directory '
                                                                                       'as this script, if a path is given, use that path' ) )

                                                                                       
p.add_argument('--xsFile',     default=None,  type=str ,        dest='xsFile',         help='path to cross section file.  When calling AddSample in the configuration module, set useXSFile=True to get weights from the provided file')
p.add_argument('--lumi',     default=None,  type=float ,        dest='lumi',         help='Integrated luminosity (to use with xsFile)')
p.add_argument('--mcweight',     default=None,  type=float ,        dest='mcweight',         help='Weight to apply to MC samples')
p.add_argument('--outputDir',     default=None,  type=str ,        dest='outputDir',         help='output directory for histograms')
p.add_argument('--readHists',     default=False,action='store_true',   dest='readHists',         help='read histograms from root files instead of trees')

p.add_argument('--quiet',     default=False,action='store_true',   dest='quiet',         help='disable information messages')
p.add_argument('--jupyt',     default=False,action='store_true',   dest='jupyt',         help='use setting for jupyter notebook')
p.add_argument('--reload',     default=False,action='store_true',   dest='reld',         help='reload sample manager')

options = p.parse_args()

import sys
import os
import re
import math
import uuid
import copy
import imp
import ROOT
from ROOT import RooFit
from array import array
if options.reld:
        import SampleManager;reload(SampleManager)
from SampleManager import SampleManager
if options.jupyt:
    try:
        from IPython.display import display, Math, Latex
        import rootnotes
    except ImportError: print "Fail to import jupyt modules"
else: ROOT.gROOT.SetBatch(False)


samples = None

def main() :

    global samples

    if not options.baseDir.count('/eos/') and not os.path.isdir( options.baseDir ) :
        print 'baseDir not found!'
        return

    samples = SampleManager(options.baseDir, options.treeName, mcweight=options.mcweight, treeNameModel=options.treeNameModel, filename=options.fileName, base_path_model=options.baseDirModel, xsFile=options.xsFile, lumi=options.lumi, readHists=options.readHists, quiet=options.quiet)


    if options.samplesConf is not None :

       samples.ReadSamples( options.samplesConf )

       # print 'Samples ready.\n'  

       # print 'The draw syntax follows that of TTree.Draw.  Examples : '
       # 
       # print 'samples.Draw(\'met_et\', \'EventWeight && passcut_ee==1\', \'(300, 0, 300)\'\n'

       # print 'The first argument is a branch in the tree to draw'
       # print 'The second argument is a set of cuts and/or weights to apply'
       # print 'The third argument are the bin limits to use \n'

       # print 'To see all available branches do ListBranches()'


#---------------------------------------
# User functions
#---------------------------------------

#---------------------------------------
def DrawFormatted(varexp, selection, histpars=None ) :
    """ Example wrapper function.   Make a stack and then add some extra goodies """

    global samples

    print 'DrawFormatted histpars ', histpars
    samples.MakeStack(varexp, selection, histpars)

    statuslabel = ROOT.TLatex(0.4, 0.8, 'Atlas Internal')
    statuslabel.SetNDC()
    luminosity = ROOT.TLatex(0.4, 0.7, ' 10.0 fb^{-1}')
    luminosity.SetNDC()
    samples.add_decoration(statuslabel)
    samples.add_decoration(luminosity)

    samples.DrawCanvas()

    statuslabel.Draw()
    luminosity.Draw()

#---------------------------------------
def WriteCurrentHists( filename='hist.root') :
    """ write all histograms in samples to a root file """

    file = ROOT.TFile.Open( filename, 'RECREATE')

    for hist, samp in samples.samples.iteritems() :
        newhist = samp.hist.Clone(hist)
        file.cd()
        newhist.Write()

    file.Close()
        
#---------------------------------------
def SaveStack( name, can=None ) :
    if options.outputDir is None :
        print 'No outputDir given!'
        return

    samples.SaveStack( name, outputDir=options.outputDir, canname=can )



#---------------------------------------


#graveyard.py
#def MakeTAndPHists( outputfile, tagprobe_min=0, tagprobe_max=1e9, normalize=1 ) :
#def MakeTAndPPlots( ) :
#def MakeTAndPCompPlots( ) :
#def MakeTAndPCompPlotsFull( ) :
#def FitTAndPComp( ) :

#
# The following is to get the history 
#
import atexit
historyPath = os.path.expanduser("~/.pyhistory")
try:
    import readline
except ImportError:
    print "Module readline not available."
else:
    import rlcompleter
    readline.parse_and_bind("tab: complete")
    if os.path.exists(historyPath):
        readline.read_history_file(historyPath)


# -----------------------------------------------------------------
def save_history(historyPath=historyPath):
    try:
        import readline
    except ImportError:
        print "Module readline not available."
    else:
        readline.write_history_file(historyPath)

atexit.register(save_history)

main()
