import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import os
import re
import math
import selection_defs as defs
from SampleManager import SampleManager
from argparse import ArgumentParser

parser = ArgumentParser()

parser.add_argument( '--baseDirMu'  ,  dest='baseDirMu'  , default=None, help='Path to single Muon ntuples')
parser.add_argument( '--baseDirEl'  ,  dest='baseDirEl'  , default=None, help='Path to single Electron ntuples')
parser.add_argument( '--baseDirMuG' ,  dest='baseDirMuG' , default=None, help='Path to Muon + Photon ntuples')
parser.add_argument( '--baseDirElG' ,  dest='baseDirElG' , default=None, help='Path to Electron + Photon ntuples')
parser.add_argument( '--baseDirMuGNoId' ,  dest='baseDirMuGNoId' , default=None, help='Path to Muon + Photon ntuples with no photon Id')
parser.add_argument( '--baseDirElGNoId' ,  dest='baseDirElGNoId' , default=None, help='Path to Electron + Photon ntuples with no photon Id')
parser.add_argument( '--baseDirMuMu',  dest='baseDirMuMu', default=None, help='Path to DiMuon ntuples')
parser.add_argument( '--baseDirElEl',  dest='baseDirElEl', default=None, help='Path to DiElectron ntuples')
parser.add_argument( '--baseDirMuEl',  dest='baseDirMuEl', default=None, help='Path to Muon + Electron ntuples')
parser.add_argument( '--baseDirNoFilt'  ,  dest='baseDirNoFilt'  , default=None, help='Path to single Muon ntuples')
parser.add_argument( '--baseDirReco'  ,  dest='baseDirReco'  , default=None, help='Path to single Muon ntuples')
parser.add_argument( '--outputDir'  ,  dest='outputDir'  , default=None, help='Path to output directory, save plots here')
parser.add_argument( '--makeZCR'  ,  dest='makeZCR'  , default=False, action='store_true', help='Make ZCR plots')
parser.add_argument( '--makeWCR'  ,  dest='makeWCR'  , default=False, action='store_true', help='Make WCR plots')
parser.add_argument( '--makeSR'  ,  dest='makeSR'  , default=False, action='store_true', help='Make signal region plots')
parser.add_argument( '--makeJetBkg'  ,  dest='makeJetBkg'  , default=False, action='store_true', help='Make jet background plots')
parser.add_argument( '--makeEleBkg'  ,  dest='makeEleBkg'  , default=False, action='store_true', help='Make electron background plots')
parser.add_argument( '--makeSigTruth'  ,  dest='makeSigTruth'  , default=False, action='store_true', help='Make signal truth plots')
parser.add_argument( '--makeBkgTruth'  ,  dest='makeBkgTruth'  , default=False, action='store_true', help='Make background truth plots')


options = parser.parse_args()

ROOT.gROOT.SetBatch(False)
if options.outputDir is not None :
    ROOT.gROOT.SetBatch(True)


_TREENAME = 'UMDNTuple/EventTree'
_FILENAME = 'tree.root'
_XSFILE   = 'cross_sections/photon15.py'
_LUMI     = 36000
_MODULE   = 'Modules/Resonance.py'

def main() :


    sampManReco     = SampleManager( options.baseDirReco  , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManNoFilt   = SampleManager( options.baseDirNoFilt  , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManMu       = SampleManager( options.baseDirMu  , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManEl       = SampleManager( options.baseDirEl  , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManMuG      = SampleManager( options.baseDirMuG , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManElG      = SampleManager( options.baseDirElG , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManMuGNoId  = SampleManager( options.baseDirMuGNoId , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManElGNoId  = SampleManager( options.baseDirElGNoId , _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManMuMu     = SampleManager( options.baseDirMuMu, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManElEl     = SampleManager( options.baseDirElEl, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManMuEl     = SampleManager( options.baseDirMuEl, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )

    sampManReco.ReadSamples( _MODULE )
    sampManNoFilt.ReadSamples( _MODULE )
    sampManMu.ReadSamples( _MODULE )
    sampManEl.ReadSamples( _MODULE )
    sampManMuG.ReadSamples( _MODULE )
    sampManElG.ReadSamples( _MODULE )
    sampManMuGNoId.ReadSamples( _MODULE )
    sampManElGNoId.ReadSamples( _MODULE )
    sampManMuMu.ReadSamples( _MODULE )
    sampManElEl.ReadSamples( _MODULE )
    sampManMuEl.ReadSamples( _MODULE )

    #---------------------------------------
    # Z and Z+photon CR plots
    #---------------------------------------
    if options.makeZCR :
        MakeZCRPlots( sampManMuMu, sampManElEl )

    #---------------------------------------
    # Single lepton (W) plots
    #---------------------------------------
    if options.makeWCR : 
        MakeWCRPlots( sampManMu, sampManEl )

    #---------------------------------------
    # lepton + photon plots
    #---------------------------------------
    if options.makeSR :
        MakeSignalRegionPlots( sampManMuG, sampManElG )

    #---------------------------------------
    # Plots for jet background
    #---------------------------------------
    if options.makeJetBkg :
        MakeWJetsPlots( sampManMuGNoId, sampManElGNoId )

    #---------------------------------------
    # Plots for electron background
    #---------------------------------------
    if options.makeEleBkg :
        MakeEFakePlots( sampManElG )

    #---------------------------------------
    # Make generator level distributions for signal
    #---------------------------------------
    if options.makeSigTruth :
        MakeSignalTruthPlots( sampManNoFilt )

    #---------------------------------------
    # for bkg studies at truth level
    #---------------------------------------
    if options.makeBkgTruth :
        MakeBkgTruthPlots( sampManReco )

    #---------------------------------------
    # for signal shape studies
    #---------------------------------------
    #Make2DSignalPlot( sampManMuG, sampManElG )


    print '^.^ Finished ^.^'

def MakeWCRPlots( sampManMu, sampManEl ) :

    name = 'WCRPlots'
    outdir = '%s/%s' %(options.outputDir, name )


    selection_mu = defs.get_base_selection( 'mu' )
    selection_el = defs.get_base_selection( 'el' )

    weight_str = defs.get_weight_str()

    plot_var = 'mt_lep_met'

    full_sel_str_mu = ' %s * ( %s ) ' %( weight_str, selection_mu )
    full_sel_str_el = ' %s * ( %s ) ' %( weight_str, selection_el )

    #sampManMu.Draw( 'mt_lep_met', full_sel_str_mu, (50, 0, 500), 
    #               hist_config={'logy' : True, 'xlabel' : 'Transverse Mass [GeV]', 'doratio' : 1}, 
    #               label_config={'labelStyle' : 'fancy13'},
    #              )

    #if options.outputDir is not None :
    #    sampManMu.SaveStack( 'mt_lep_met_mu.pdf' ,outdir, 'base' ) 
    #else :
    #    raw_input('cont')

    #sampManEl.Draw( 'mt_lep_met', full_sel_str_el, (50, 0, 500), 
    #               hist_config={'logy' : True, 'xlabel' : 'Transverse Mass [GeV]', 'doratio' : 1}, 
    #               label_config={'labelStyle' : 'fancy13'},
    #              )

    #if options.outputDir is not None :
    #    sampManEl.SaveStack( 'mt_lep_met_el.pdf' ,outdir, 'base' ) 
    #else :
    #    raw_input('cont')

    sampManMu.Draw( 'met_pt', full_sel_str_mu, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Missing Transverse Momentum [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManMu.SaveStack( 'met_pt_mu.pdf' ,outdir, 'base' ) 
    else :
        raw_input('cont')

    sampManEl.Draw( 'met_pt', full_sel_str_el, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Missing Transverse Momentum [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'met_pt_el.pdf' , outdir, 'base' ) 
    else :
        raw_input('cont')

    sampManMu.Draw( 'mu_pt[0]', full_sel_str_mu, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Lepton p_{T} [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManMu.SaveStack( 'lep_pt_mu.pdf' ,outdir, 'base' ) 
    else :
        raw_input('cont')

    sampManEl.Draw( 'el_pt[0]', full_sel_str_el, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Lepton p_{T} [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'lep_pt_el.pdf' , outdir, 'base' ) 
    else :
        raw_input('cont')


def MakeZCRPlots( sampManMu, sampManEl ) :

    name = 'ZCRPlots'
    outdir = '%s/%s' %(options.outputDir, name )

    selection_mu = defs.get_base_selection( 'mumu' )
    selection_el = defs.get_base_selection( 'elel' )

    weight_str = defs.get_weight_str()

    plot_var = 'm_ll'

    full_sel_str_mu = ' %s * ( %s ) ' %( weight_str, selection_mu )
    full_sel_str_el = ' %s * ( %s ) ' %( weight_str, selection_el )

    
    sampManMu.Draw( plot_var, full_sel_str_mu, (100, 0, 1000), 
                   hist_config={'logy' : True, 'xlabel' : 'Dilepton Mass [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManMu.SaveStack( 'm_ll_mumu.pdf' ,outdir, 'base' ) 
    else :
        raw_input('cont')

    sampManEl.Draw( plot_var, full_sel_str_el, (100, 0, 1000), 
                   hist_config={'logy' : True, 'xlabel' : 'Dilepton Mass [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'm_ll_elel.pdf' , outdir, 'base' ) 
    else :
        raw_input('cont')

    selection_ph = 'ph_n==1 && ph_pt[0] > 15 && ph_IsEB[0]'

    full_sel_str_mug = ' %s * ( %s ) ' %( weight_str, ' && '.join( [selection_mu, selection_ph] ) )
    full_sel_str_elg = ' %s * ( %s ) ' %( weight_str, ' && '.join( [selection_el, selection_ph] ) )

    #sampManMu.Draw( plot_var, full_sel_str_mug, (50, 0, 500), 
    #               hist_config={'logy' : True, 'xlabel' : 'Dilepton Mass [GeV]', 'doratio' : 1}, 
    #               label_config={'labelStyle' : 'fancy13'},
    #              )

    #if options.outputDir is not None :
    #    sampManMu.SaveStack( 'm_ll_mumug.pdf' ,outdir, 'base' ) 
    #else :
    #    raw_input('cont')

    #sampManEl.Draw( plot_var, full_sel_str_elg, (50, 0, 500), 
    #               hist_config={'logy' : True, 'xlabel' : 'Dilepton Mass [GeV]', 'doratio' : 1}, 
    #               label_config={'labelStyle' : 'fancy13'},
    #              )

    #if options.outputDir is not None :
    #    sampManEl.SaveStack( 'm_ll_elelg.pdf' , outdir, 'base' ) 
    #else :
    #    raw_input('cont')


    #sampManMu.Draw( 'dr_lep_ph', full_sel_str_mug, (50, 0, 5), 
    #               hist_config={'logy' : True, 'xlabel' : '#Delta R (#ell, #gamma)', 'doratio' : 1}, 
    #               label_config={'labelStyle' : 'fancy13'},
    #              )

    #if options.outputDir is not None :
    #    sampManMu.SaveStack( 'dr_lep_ph_mumug.pdf', outdir, 'base' ) 
    #else :
    #    raw_input('cont')

    #sampManEl.Draw( 'dr_lep_ph', full_sel_str_elg, (50, 0, 5), 
    #               hist_config={'logy' : True, 'xlabel' : '#Delta R (#ell, #gamma)', 'doratio' : 1}, 
    #               label_config={'labelStyle' : 'fancy13'},
    #              )

    #if options.outputDir is not None :
    #    sampManEl.SaveStack( 'dr_lep_ph_elelg.pdf' , outdir, 'base' )
    #else :
    #    raw_input('cont')

    sampManMu.Draw( 'mu_pt[0]', full_sel_str_mu, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Lead Lepton p_{T} [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManMu.SaveStack( 'leadlep_pt_mu.pdf' ,outdir, 'base' ) 
    else :
        raw_input('cont')

    sampManEl.Draw( 'el_pt[0]', full_sel_str_el, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Lead Lepton p_{T} [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'leadlep_pt_el.pdf' , outdir, 'base' ) 
    else :
        raw_input('cont')


    sampManMu.Draw( 'mu_pt[1]', full_sel_str_mu, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Sublead Lepton p_{T} [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManMu.SaveStack( 'lep_pt_mu.pdf' ,outdir, 'base' ) 
    else :
        raw_input('cont')

    sampManEl.Draw( 'el_pt[1]', full_sel_str_el, (50, 0, 500), 
                   hist_config={'logy' : True, 'xlabel' : 'Sublead Lepton p_{T} [GeV]', 'doratio' : 1}, 
                   label_config={'labelStyle' : 'fancy13'},
                  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'lep_pt_el.pdf' , outdir, 'base' ) 
    else :
        raw_input('cont')



def MakeEFakePlots( sampManEl ) :

    name = 'EFakePlots'
    outdir = '%s/%s' %(options.outputDir, name )

    binning = ( 100, 0, 1000 )
    plot_var = 'ph_pt[0]'
    selection_el = defs.get_base_selection( 'el' )

    weight_str = defs.get_weight_str()

    selection_ph = 'ph_n==1 && ph_passEleVeto[0] == 0 && !isBlinded'

    full_sel_str_el = ' %s * ( %s ) ' %( weight_str, ' && ' .join( [selection_el, selection_ph] ) )

    sampManEl.Draw( plot_var, full_sel_str_el, binning  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'ph_pt_elg_failVeto.pdf' , outdir, 'base' )
    else :
        raw_input('cont')

    selection_ph = 'ph_n==1 && ph_passEleVeto[0] == 1 && !isBlinded'

    full_sel_str_el = ' %s * ( %s ) ' %( weight_str, ' && ' .join( [selection_el, selection_ph] ) )

    sampManEl.Draw( plot_var, full_sel_str_el, binning  )

    if options.outputDir is not None :
        sampManEl.SaveStack( 'ph_pt_elg_passVeto.pdf' , outdir, 'base' )
    else :
        raw_input('cont')


def MakeWJetsPlots( sampManMu, sampManEl ) :

    name = 'WJetsPlots'
    outdir = '%s/%s' %(options.outputDir, name )

    #ph_id_cats = ['chIso', 'sigmaIEIE']
    ph_id_cats = [ 'chIso']
    ph_eta_cats = ['EB', 'EE']
    ph_pt_cats = [(15, 50), (50, 1e9 )]

    binning = {
               'chIso' : { 
                          'EB' : ( 50, 0, 25 ), 
                          'EE' : ( 50, 0, 25 ), 
                         },
               'sigmaIEIE' : { 
                          'EB' : ( 50, 0, 0.05 ), 
                          'EE' : ( 50, 0, 0.1 ), 
                         },
    }

    xlabels = {'chIso' : 'Charged hadron isolation', 'sigmaIEIE' : '#sigma i#eta i#eta' }
    yranges = {'chIso' : ( 10, 1e7), 'sigmaIEIE' : ( 1, 1e6 ) }

    extra_tags = {'EB' : 'Barrel Photons', 'EE' :  'Endcap Photons' }

    selection_mu = defs.get_base_selection( 'mu' )
    selection_el = defs.get_base_selection( 'el' )

    weight_str = defs.get_weight_str()

    for iid in ph_id_cats :

        phid_selection = defs.get_phid_selection( iid )

        phid_idx = defs.get_phid_idx( iid )

        plot_var = '%s[%s]' %( defs.get_phid_cut_var( iid ), phid_idx )

        for ipt in ph_pt_cats :

            ph_pt_sel  = 'ph_pt[%s] > %d && ph_pt[%s] < %d' %( phid_idx, ipt[0], phid_idx, ipt[1] )

            for ieta in ph_eta_cats :

                ph_eta_sel = 'ph_Is%s[%s]' %( ieta, phid_idx )

                all_selections_mu = [selection_mu, phid_selection, ph_pt_sel, ph_eta_sel]

                full_sel_str_mu = '%s * ( %s )' %( weight_str, ' && '.join( all_selections_mu ) )

                ymin = yranges[iid][0]
                ymax = yranges[iid][1]

                ylabel = 'Events / %.03f ' %( (binning[iid][ieta][2] - binning[iid][ieta][1]) / float( binning[iid][ieta][0] ) )

                sampManMu.Draw( plot_var, full_sel_str_mu, binning[iid][ieta], 
                               hist_config={'logy' : True, 'xlabel' : xlabels[iid], 'ylabel' : ylabel, 'ymin' : ymin, 'ymax' : ymax}, 
                               label_config={'labelStyle' : 'fancy13', 'extra_label': extra_tags[ieta], 'extra_label_loc':(0.2, 0.87)},
                              )

                if options.outputDir is not None :
                    sampManMu.SaveStack( 'ph_%s_mu_%s_phpt_%d_%d.pdf' %(iid, ieta, ipt[0], ipt[1]) ,outdir, 'base' )
                else :
                    raw_input('cont')


def MakeBkgTruthPlots( sampManReco ) :

    name = 'BkgTruthPlots'
    outdir = '%s/%s' %(options.outputDir, name )

    sampManReco.CompareSelections( 'trueph_pt[0]', ['trueph_n>0']*3, ['WGToLNuG-amcatnloFXFX', 'WGToLNuG_PtG-130-amcatnloFXFX', 'WGToLNuG_PtG-500-amcatnloFXFX'], ( 50, 0, 1000 ), hist_config={'colors' : [ROOT.kBlack, ROOT.kRed, ROOT.kBlue], 'logy' : True, 'xlabel' : 'Photon p_{T} [GeV]'  }, legend_config={'legend_entries' : ['W#gamma Inclusive', 'W#gamma, p_{T}^{#gamma} > 130 GeV', 'W#gamma, p_{T}^{#gamma} > 500 GeV'] } )

    if options.outputDir is not None :
        sampManReco.SaveStack( 'Wgamma_NLO_comp' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManReco.CompareSelections( 'trueph_pt[0]', ['trueph_n>0']*3, ['WGToLNuG-madgraphMLM', 'WGToLNuG_PtG-130-madgraphMLM', 'WGToLNuG_PtG-500-madgraphMLM'], ( 50, 0, 1000 ), hist_config={'colors' : [ROOT.kBlack, ROOT.kRed, ROOT.kBlue], 'logy' : True, 'xlabel' : 'Photon p_{T} [GeV]'  }, legend_config={'legend_entries' : ['W#gamma Inclusive', 'W#gamma, p_{T}^{#gamma} > 130 GeV', 'W#gamma, p_{T}^{#gamma} > 500 GeV'] } )

    if options.outputDir is not None :
        sampManReco.SaveStack( 'Wgamma_LO_comp' ,outdir, 'base' )
    else :
        raw_input('cont')

def MakeSignalTruthPlots( sampManNoFilt ) :

    name = 'TruthRecoCompPlots'
    outdir = '%s/%s' %(options.outputDir, name )

    sig_samples = []
    for samp in sampManNoFilt.get_samples() :

        print samp.name

        if samp.name.count('ResonanceMass' ) : 

            res = re.match('ResonanceMass(\d+)_width(\w+|\d+)', samp.name )
            if res is not None :
                mass = int( res.group(1) )
                width =  res.group(2) 

            sig_samples.append( (samp.name, mass, width) )

    for samp, mass, width in sig_samples :

        xmax = int(mass*1.5)
        xmin = 0

        while xmax%10 != 0 :
            xmax+=1

        bin_width = 10
        if mass < 500 : 
            bin_width = 5
        if mass >= 1000 :
            bin_width = 20
        if mass >= 2000 :
            bin_width = 50

        width_val = width 
        if width == '0p01' :
            width_val = '0.01'

        nbins = int((xmax - xmin)/bin_width)

        sampManNoFilt.CompareVars( ['truemt_res', 'mt_res','truemt_res', 'mt_res'], ['isWMuDecay==1 && trueph_n>0 ', 'mu_pt30_n==1 && ph_n==1','isWElDecay==1 && trueph_n>0 ', 'el_pt30_n==1 && ph_n==1',], [samp]*4, (nbins, xmin, xmax ), hist_config={'colors' : [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen], 'xlabel' : 'Transverse Mass [GeV]', 'ylabel' : 'Normalized Events', 'normalize' : 1, 'ymin' : 0.0001, 'ymax' : 1.0, 'logy' : 1 } , legend_config={'legendLoc' : 'TopLeft', 'legendWiden' : 1.2, 'legendCompress' : 1.2, 'legend_entries' : ['Truth, Muon channel', 'Reconstructed, Muon channel','Truth, Electron channel', 'Reconstructed, Electron channel',]}, label_config={'labelStyle' : 'fancy13', 'extra_label':'Mass = %d GeV, width = %s' %( mass, width_val ) + '%', 'extra_label_loc':(0.45, 0.87)})

        if options.outputDir is not None :
            sampManNoFilt.SaveStack( 'mt_res_trueRecoComp_M%d_width%s'%(mass, width) ,outdir, 'base' )
        else :
            raw_input('cont')

def MakeSignalRegionPlots( sampManMuG, sampManElG ) :

    sampManMuG.deactivate_sample( 'Data')
    sampManElG.deactivate_sample( 'Data')

    name = 'SignalRegionPlots'
    outdir = '%s/%s' %(options.outputDir, name )

    legend_conf = sampManMuG.config_legend( legendCompress = 1.0, legendWiden = 1.4 )
    legend_conf_s = sampManMuG.config_legend( legendCompress = 0.8, legendWiden = 1.0 )
    label_loc_conf = (0.2, 0.87)

    weight_str = defs.get_weight_str()
    sel_base_mu = defs.get_base_selection( 'mu' )
    sel_base_el = defs.get_base_selection( 'el' )

    ph_str = 'ph_n==1 && ph_IsEB[0] && ph_pt[0] > 50 && ph_passEleVeto[0]==1' 

    sel_mu_nominal = '%s * ( %s && %s )'%( weight_str, sel_base_mu, ph_str ) 
    sel_el_nominal = '%s * ( %s && %s )'%( weight_str, sel_base_el, ph_str ) 


    sampManMuG.Draw( 'mt_res', sel_mu_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mt_res_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'mt_res', sel_el_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'mt_res_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'mt_lep_met_ph', sel_mu_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mt_lep_met_ph_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'mt_lep_met_ph', sel_el_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'mt_lep_met_ph_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')



 # kak: Z_mass
    sampManElG.Draw( 'm_ll', sel_mu_nominal, (80, 0, 1000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'm_ll_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')
    sampManElG.Draw( 'm_lep_ph', sel_mu_nominal, (80, 0, 1000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )
    if options.outputDir is not None :
        sampManElG.SaveStack( 'm_lep_ph_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')



    sampManMuG.Draw( 'm_lep_met_ph', sel_mu_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'm_lep_met_ph_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'm_lep_met_ph', sel_el_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'm_lep_met_ph_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'recoM_lep_nu_ph', sel_mu_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Mass (W mass constraint) [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'constrained_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'recoM_lep_nu_ph', sel_el_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Mass (W mass constraint) [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'constrained_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'ph_pt[0]', sel_mu_nominal, (80, 0, 2000 ), hist_config={'xlabel' : 'Photon p_{T}[GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'ph_pt_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'ph_pt[0]', sel_el_nominal, (80, 0, 2000 ), hist_config={'xlabel' : 'Photon p_{T} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'ph_pt_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    theta = math.atan(1./2.0)

#    sampManMuG.Draw( 'ph_pt[0] *sin( %f ) + mt_res*cos(%f)' %(theta,theta), sel_mu_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Combined m_{T} and p_{T} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )
#
#    if options.outputDir is not None :
#        sampManMuG.SaveStack( 'rotated_mug_nodata' ,outdir, 'base' )
#    else :
#        raw_input('cont')
#
#    sampManElG.Draw( 'ph_pt[0] *sin( %f ) + mt_res*cos(%f)' %(theta,theta), sel_el_nominal, (80, 0, 4000 ), hist_config={'xlabel' : 'Combined m_{T} and p_{T} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )
#
#    if options.outputDir is not None :
#        sampManElG.SaveStack( 'rotated_elg_nodata' ,outdir, 'base' )
#    else :
#        raw_input('cont')

    sampManMuG.Draw( 'fabs(dphi_lep_ph)', sel_mu_nominal, (32, 0, 3.2 ), hist_config={'xlabel' : '#Delta #phi (l, #gamma)', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 200000000}, legend_config={'legendLoc' : 'TopLeft'}, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':(0.6, 0.87)} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'dphi_lep_ph_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'fabs(dphi_lep_ph)', sel_el_nominal, (32, 0, 3.2 ), hist_config={'xlabel' : '#Delta #phi (l, #gamma)', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 10000000000}, legend_config={'legendLoc' : 'TopLeft'}, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':(0.6, 0.87)} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'dphi_lep_ph_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'fabs(dphi_lep_met)', sel_mu_nominal, (32, 0, 3.2 ), hist_config={'xlabel' : '#Delta #phi (l, MET)', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 200000000}, legend_config={'legendLoc' : 'TopLeft'}, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':(0.6, 0.87)} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'dphi_lep_met_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'fabs(dphi_lep_met)', sel_el_nominal, (32, 0, 3.2 ), hist_config={'xlabel' : '#Delta #phi (l, MET)', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 10000000000}, legend_config={'legendLoc' : 'TopLeft'}, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':(0.6, 0.87)} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'dphi_lep_met_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'recoM_lep_nu_ph-mt_res', sel_mu_nominal, (100, -500, 500 ), hist_config={'xlabel' : 'Reconstructed Mass Difference [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mass_diff_ph_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'recoM_lep_nu_ph-mt_res', sel_el_nominal , (100, -500, 500 ), hist_config={'xlabel' : 'Reconstructed Mass Difference [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'mass_diff_ph_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'mu_pt[0]', sel_mu_nominal, (100, 0, 500 ), hist_config={'xlabel' : 'Muon p_{T} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mu_pt_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'el_pt[0]', sel_el_nominal, (100, 0, 500 ), hist_config={'xlabel' : 'Electron p_{T} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'el_pt_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'met_pt', sel_mu_nominal, (100, 0, 500 ), hist_config={'xlabel' : 'E_{T}^{miss} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'met_pt_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'met_pt', sel_el_nominal, (100, 0, 500 ), hist_config={'xlabel' : 'E_{T}^{miss} [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'met_pt_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'mt_lep_met', sel_mu_nominal, (100, 0, 500 ), hist_config={'xlabel' : 'W Transverse mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mt_lep_met_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'mt_lep_met', sel_el_nominal, (100, 0, 500 ), hist_config={'xlabel' : 'W Transverse mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf_s, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'mt_lep_met_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'jet_n', sel_mu_nominal, (10, 0, 10 ), hist_config={'xlabel' : 'Jet Multiplicity', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'jet_n_mug_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'jet_n', sel_el_nominal, (10, 0, 10 ), hist_config={'xlabel' : 'Jet Multiplicity', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 2000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'jet_n_elg_nodata' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.activate_sample( 'Data')
    sampManElG.activate_sample( 'Data')

    sampManMuG.Draw( 'mt_res', '%s * ( %s && mt_res < 500 && !isBlinded )'%(weight_str, sel_base_mu), (50, 0, 500 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 0, 'ymin' : 0.01, 'ymax' : 10000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mt_res_mug' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'mt_res', '%s * ( %s && mt_res < 500  && !isBlinded)'%(weight_str, sel_base_el), (50, 0, 500 ), hist_config={'xlabel' : 'Transverse Mass [GeV]', 'logy' : 0, 'ymin' : 0.01, 'ymax' : 10000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'mt_res_elg' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'mt_lep_met_ph', '%s * ( %s && %s )'%(weight_str, sel_base_mu, ph_str), (50, 0, 500 ), hist_config={'xlabel' : 'W Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'mt_lep_met_ph_mug' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'mt_lep_met_ph', '%s * ( %s && %s )'%(weight_str, sel_base_el, ph_str), (50, 0, 500 ), hist_config={'xlabel' : 'W Transverse Mass [GeV]', 'logy' : 1, 'ymin' : 0.01, 'ymax' : 1000000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'mt_lep_met_ph_elg' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManMuG.Draw( 'ph_pt[0]', '%s * ( %s && ph_pt[0] < 250 )'%(weight_str, sel_base_mu), (50, 0, 250 ), hist_config={'xlabel' : 'Photon p_{T}[GeV]', 'logy' : 0, 'ymin' : 0.01, 'ymax' : 100000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Muon Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManMuG.SaveStack( 'ph_pt_mug' ,outdir, 'base' )
    else :
        raw_input('cont')

    sampManElG.Draw( 'ph_pt[0]', '%s * ( %s && ph_pt[0] < 250 )'%(weight_str, sel_base_el), (50, 0, 250 ), hist_config={'xlabel' : 'Photon p_{T} [GeV]', 'logy' : 0, 'ymin' : 0.01, 'ymax' : 100000}, legend_config=legend_conf, label_config = {'labelStyle' : 'fancy13', 'extra_label':'Electron Channel', 'extra_label_loc':label_loc_conf} )

    if options.outputDir is not None :
        sampManElG.SaveStack( 'ph_pt_elg' ,outdir, 'base' )
    else :
        raw_input('cont')


def Make2DSignalPlot( sampManMuG, sampManElG ) :

    name = 'SignalBackground2DComp'
    out_dir = '%s/%s' %( options.outputDir, name )

    if options.outputDir is not None :
        if not os.path.isdir( out_dir ) :
            os.makedirs( out_dir )

    bkg_samp = 'MCBackground'
    var = 'mt_res:ph_pt[0]' #y:x
    selection_mu = 'mu_n==1 && ph_n==1'
    selection_el = 'el_n==1 && ph_n==1'

    sig_masses = {}
    for samp in sampManMuG.get_samples() :
        if samp.name.count('width5') :
            continue
        if samp.name.count('ResonanceMass') == 0 :
            continue

        res = re.match('ResonanceMass(\d+).*', samp.name )

        if res is not None :
            sig_masses[int( res.group(1))] = samp.name

    myf = ROOT.TF1( 'myf', '[0] + [1]*x', 0, 2000 )
    myf.SetParameter(0, 0. )
    myf.SetParameter(1, 2. )

    for mass, sig in sig_masses.iteritems() :

        mt_res_min = mass/2.
        mt_res_max = mass*1.2

        ph_pt_min = mt_res_min/2.
        ph_pt_max = mt_res_max/2.

        binning = ( 50, ph_pt_min, ph_pt_max, 50, mt_res_min, mt_res_max )

        selection_cut = 'mt_res > %f && mt_res < %f && ph_pt[0] > %f && ph_pt[0] < %f ' %( mt_res_min, mt_res_max, ph_pt_min, ph_pt_max ) 

        sampManMuG.create_hist(sig, var, selection_mu + ' && %s  ' %(selection_cut) ,binning  )
        sampManElG.create_hist(sig, var, selection_el + ' && %s  ' %(selection_cut) ,binning  )
        
        sampManMuG.create_hist(bkg_samp, var, selection_mu + ' && %s  ' %(selection_cut) , binning )
        sampManElG.create_hist(bkg_samp, var, selection_el + ' && %s  ' %(selection_cut) , binning )
        
        sig_hist_mu = sampManMuG.get_samples( name=sig )[0].hist
        sig_hist_el = sampManElG.get_samples( name=sig )[0].hist

        sig_hist_profx_mu = sig_hist_mu.ProfileX('%s_px' %sig_hist_mu.GetName() )
        sig_hist_profx_el = sig_hist_el.ProfileX('%s_px' %sig_hist_el.GetName() )

        can_sig_profx_mu = ROOT.TCanvas( 'can_sig_profx_mu_M%d' %mass, 'can_sig_profx_mu_M%d' %mass )
        sig_hist_profx_mu.Draw()

        sig_hist_profx_mu.Fit( myf, "+", "", ph_pt_min, (mass/2.)*0.9)

        sig_hist_profx_mu.GetXaxis().SetTitle( 'Photon pT [GeV]' )
        sig_hist_profx_mu.GetYaxis().SetTitle( 'Profiled Transverse Mass [GeV]' )
        sig_hist_profx_mu.SetStats( 0 )
        sig_hist_profx_mu.SetMarkerStyle(20)
        sig_hist_profx_mu.SetMarkerSize(1)
        sig_hist_profx_mu.SetMarkerColor(ROOT.kBlack)
        sig_hist_profx_mu.SetLineColor(ROOT.kBlack)


        latex_mu = ROOT.TLatex(0.15, 0.8, 'Slope = %.1f' %myf.GetParameter(1) )
        latex_mu.SetNDC()
        latex_mu.SetX(0.15)
        latex_mu.SetY(0.8)
        latex_mu.Draw()
        latex_mu_title = ROOT.TLatex(0.1, 0.91, 'Muon Channel, Mass = %d GeV' %mass )
        latex_mu_title.SetNDC()
        latex_mu_title.SetX(0.1)
        latex_mu_title.SetY(0.91)
        latex_mu_title.Draw()

        if options.outputDir is not None :
            can_sig_profx_mu.SaveAs( '%s/%s.pdf' %( out_dir, can_sig_profx_mu.GetName() ) )
        else :
            raw_input('cont')
            
        can_sig_profx_el = ROOT.TCanvas( 'can_sig_profx_el_M%d' %mass, 'can_sig_profx_el_M%d' %mass )
        sig_hist_profx_el.Draw()

        sig_hist_profx_el.Fit( myf, "+", "", ph_pt_min, (mass/2.)*0.9)

        sig_hist_profx_el.GetXaxis().SetTitle( 'Photon pT [GeV]' )
        sig_hist_profx_el.GetYaxis().SetTitle( 'Transverse Mass [GeV]' )
        sig_hist_profx_el.SetStats( 0 )
        sig_hist_profx_el.SetMarkerStyle(20)
        sig_hist_profx_el.SetMarkerSize(1)
        sig_hist_profx_el.SetMarkerColor(ROOT.kBlack)
        sig_hist_profx_el.SetLineColor(ROOT.kBlack)

        latex_el = ROOT.TLatex(0.15, 0.8, 'Slope = %.1f' %myf.GetParameter(1) )
        latex_el.SetNDC()
        latex_el.SetX(0.15)
        latex_el.SetY(0.8)
        latex_el.Draw()
        latex_el_title = ROOT.TLatex(0.1, 0.91, 'Electron Channel, Mass = %d GeV' %mass )
        latex_el_title.SetNDC()
        latex_el_title.SetX(0.1)
        latex_el_title.SetY(0.91)
        latex_el_title.Draw()

        if options.outputDir is not None :
            can_sig_profx_el.SaveAs( '%s/%s.pdf' %( out_dir, can_sig_profx_el.GetName() ) )
        else :
            raw_input('cont')

        bkg_hist_mu = sampManMuG.get_samples( name=bkg_samp )[0].hist
        bkg_hist_el = sampManElG.get_samples( name=bkg_samp )[0].hist

        sig_hist_mu.GetXaxis().SetTitle( 'Photon pT [GeV]' )
        sig_hist_mu.GetYaxis().SetTitle( 'Transverse mass [GeV]' )
        sig_hist_mu.GetYaxis().SetTitleOffset( 1.1 )

        sig_hist_el.GetXaxis().SetTitle( 'Photon pT [GeV]' )
        sig_hist_el.GetYaxis().SetTitle( 'Transverse mass [GeV]' )
        sig_hist_el.GetYaxis().SetTitleOffset( 1.1 )

        bkg_hist_mu.GetXaxis().SetTitle( 'Photon pT [GeV]' )
        bkg_hist_mu.GetYaxis().SetTitle( 'Transverse mass [GeV]' )
        bkg_hist_mu.GetYaxis().SetTitleOffset( 1.1 )

        bkg_hist_el.GetXaxis().SetTitle( 'Photon pT [GeV]' )
        bkg_hist_el.GetYaxis().SetTitle( 'Transverse mass [GeV]' )
        bkg_hist_el.GetYaxis().SetTitleOffset( 1.1 )

        can_sig_mu = ROOT.TCanvas( 'can_sig_mu_M%d' %mass, 'can_sig_mu_M%d' %mass )
        can_sig_mu.SetTitle()
        sig_hist_mu.SetStats(0)
        sig_hist_mu.Draw('colz')
        print sig_hist_mu.Integral()
        can_sig_mu.SetLogz()
        if options.outputDir is not None :
            can_sig_mu.SaveAs( '%s/%s.pdf' %( out_dir, can_sig_mu.GetName() ) )
        else :
            raw_input('cont')
            
        can_sig_el = ROOT.TCanvas( 'can_sig_el_M%d'%mass, 'can_sig_el_M%d'%mass )
        can_sig_el.SetTitle()
        sig_hist_el.SetStats(0)
        sig_hist_el.Draw('colz')
        can_sig_el.SetLogz()
        if options.outputDir is not None :
            can_sig_el.SaveAs( '%s/%s.pdf' %( out_dir, can_sig_el.GetName() ) )
        else :
            raw_input('cont')
            
        can_bkg_mu = ROOT.TCanvas( 'can_bkg_mu_M%d'%mass, 'can_bkg_mu_M%d'%mass )
        can_bkg_mu.SetTitle()
        bkg_hist_mu.SetStats(0)
        bkg_hist_mu.Draw('colz')
        can_bkg_mu.SetLogz()
        if options.outputDir is not None :
            can_bkg_mu.SaveAs( '%s/%s.pdf' %( out_dir, can_bkg_mu.GetName() ) )
        else :
            raw_input('cont')
            
        can_bkg_el = ROOT.TCanvas( 'can_bkg_el_M%d'%mass, 'can_bkg_el_M%d'%mass )
        can_bkg_el.SetTitle()
        bkg_hist_el.SetStats(0)
        bkg_hist_el.Draw('colz')
        can_bkg_el.SetLogz()
        if options.outputDir is not None :
            can_bkg_el.SaveAs( '%s/%s.pdf' %( out_dir, can_bkg_el.GetName() ) )
        else :
            raw_input('cont')
            

        
main()
