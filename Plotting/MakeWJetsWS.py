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

    #workspace_signal   = ROOT.RooWorkspace( 'workspace_signal' )
    #workspace_wgamma   = ROOT.RooWorkspace( 'workspace_wgamma' )
    #workspace_wgammalo = ROOT.RooWorkspace( 'workspace_wgammalo' )
    #workspace_top      = ROOT.RooWorkspace( 'workspace_top' )
    #workspace_zgamma   = ROOT.RooWorkspace( 'workspace_zgamma' )
    #wjets              = ROOT.RooWorkspace( 'workspace_wjets' )
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

   # ph_idx_sr =  defs.get_phid_idx( 'medium' )
   # ph_idx_den = defs.get_phid_idx( num_var, _var )
   # ph_idx_num = defs.get_phid_idx( num_var )
   # ph_idx_ = defs.get_phid_idx( shape_var )

   # xmin = xvar.getMin()
   # xmax = xvar.getMax()

   # #---------------------------------------
   # # Add eta cuts, (IsEB, IsEE)
   # #---------------------------------------
   # eta_str_ = 'ph_Is%s[%s]' %( eta_cut, ph_idx_shape )
   # eta_str_den = 'ph_Is%s[%s]' %( eta_cut, ph_idx_den)
   # eta_str_num = 'ph_Is%s[%s]' %( eta_cut, ph_idx_num )
   # eta_str_sr = 'ph_Is%s[%s]' %( eta_cut, ph_idx_sr )

   # #---------------------------------------
   # # Add additional cuts, mainly restricting
   # # the fitted variable to the plotting limits
   # #---------------------------------------
   # addtl_cuts_sr = 'ph_pt[%s] > 50 && %s > %d && %s < %d '     %( ph_idx_sr, plot_var, xmin, plot_var, xmax )
   # addtl_cuts_den = 'ph_pt[%s] > 50 && %s > %d && %s < %d  '   %( ph_idx_den, plot_var, xmin, plot_var, xmax )
   # addtl_cuts_num = 'ph_pt[%s] > 50 && %s > %d && %s < %d  '   %( ph_idx_num, plot_var, xmin, plot_var, xmax )
   # addtl_cuts_ = 'ph_pt[%s] > 50 && %s > %d && %s < %d  ' %( ph_idx_shape, plot_var, xmin, plot_var, xmax )


#    #---------------------------------------
#    # Get the cuts that define the photon
#    # sideband regions
#    #---------------------------------------
#    cut_str_base = ' {var}[{idx}] > {val_low} && {var}[{idx}] < {val_high}'
#
#    cut_var_ = defs.get_phid_cut_var( shape_var )
#    cut_var_num = defs.get_phid_cut_var( num_var )
#
#    cut_vals_ = get_cut_defaults( shape_var, eta_cut )
#    cut_vals_num = get_cut_defaults( num_var, eta_cut )
#
#    cut_str_ = cut_str_base.format(var=cut_var_shape, idx=ph_idx_shape, val_low=cut_vals_shape[0], val_high=cut_vals_shape[1] )
#    cut_str_num = cut_str_base.format(var=cut_var_num, idx=ph_idx_num, val_low=cut_vals_num[0], val_high=cut_vals_num[1] )
#
#    cut_str_den_1 = cut_str_base.format(var=cut_var_, idx=ph_idx_den, val_low=cut_vals_shape[0], val_high=cut_vals_shape[1] )
#    cut_str_den_2 = cut_str_base.format(var=cut_var_num, idx=ph_idx_den, val_low=cut_vals_num[0], val_high=cut_vals_num[1] )
#
#    cut_str_den = cut_str_den_1 + ' && ' + cut_str_den_2
#
    #---------------------------------------
    # put the cuts together
   # #---------------------------------------
   # full_sel_ = ' && '.join( [sel_base, ph_selection_shape, eta_str_shape, addtl_cuts_shape, cut_str_shape ] )
   # full_sel_num   = ' && '.join( [sel_base, ph_selection_num, eta_str_num, addtl_cuts_num, cut_str_num] )
   # full_sel_den   = ' && '.join( [sel_base, ph_selection_den, eta_str_den, addtl_cuts_den, cut_str_den] )
   # full_sel_sr    = ' && '.join( [sel_base, ph_selection_sr, eta_str_sr, addtl_cuts_sr] )
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

 #   #---------------------------------------
 #   # make fit managers
 #   #---------------------------------------
 #   fitMan_ = FitManager( 'dijet', 2, 'wjets_shape', hist_shape, plot_var, eta_cut, xvar, label_shape, useRooFit=False )
 #   fitMan_num   = FitManager( 'dijet', 2, 'wjets_num', hist_num, plot_var, eta_cut, xvar, label_num, useRooFit=False )
 #   fitMan_den   = FitManager( 'dijet', 2, 'wjets_den', hist_den, plot_var, eta_cut, xvar, label_den, useRooFit=False )

 #   #---------------------------------------
 #   # Do the fits
 #   #---------------------------------------
 #   result_= fitMan_shape.fit_histogram( ws )
 #   result_num = fitMan_num.fit_histogram( ws )
 #   result_den = fitMan_den.fit_histogram( ws )

 #   #---------------------------------------
 #   # save the results
 #   #---------------------------------------
 #   fitMan_den.save_fit( sampMan, ws, logy=True )
 #   fitMan_.save_fit( sampMan, ws, logy=True )
 #   fitMan_num.save_fit( sampMan, ws, logy=True )

 #   #---------------------------------------
 #   # calculate the function that describes
 #   # the ratio of num/den
 #   #---------------------------------------
 #   _integral_num   =  ws.pdf( 'dijet_%s' %label_num ).getNormIntegral(ROOT.RooArgSet( xvar ) )
 #   _integral_den   =  ws.pdf( 'dijet_%s' %label_den).getNormIntegral(ROOT.RooArgSet( xvar ) )
 #   _integral_ =  ws.pdf( 'dijet_%s' %label_shape).getNormIntegral(ROOT.RooArgSet( xvar ) )

 #   func_integral_num    = _integral_num.getValV()
 #   func_integral_den    = _integral_den.getValV()
 #   func_integral_  = _integral_shape.getValV()

 #   hist_integral_num = result_num['integral'].n
 #   hist_integral_den = result_den['integral'].n
 #   hist_integral_ = result_shape['integral'].n

 #   norm_num = hist_integral_num / func_integral_num
 #   norm_den = hist_integral_den / func_integral_den
 #   norm_ = hist_integral_shape / func_integral_shape

 #   print 'func integral Num = ', func_integral_num
 #   print 'hist integral Num = ', hist_integral_num
 #   print 'normalization Num = ', norm_num

 #   print 'func integral Den = ', func_integral_den
 #   print 'hist integral Den = ', hist_integral_den
 #   print 'normalization Den = ', norm_den

 #   print 'func integral Shape = ', func_integral_
 #   print 'hist integral Shape = ', hist_integral_
 #   print 'normalization Shape = ', norm_

 #   power_pred_name    = 'power_pred_%s' %suffix
 #   logcoef_pred_name  = 'logcoef_pred_%s' %suffix
 #   #power_ratio_name   = 'power_ratio_%s' %suffix
 #   #logcoef_ratio_name = 'logcoef_ratio_%s' %suffix

 #   #name_power_num     ='power' 
 #   #name_power_den     ='power' 
 #   #name_power_   ='power' 

 #   #name_logcoef_num   ='logcoef1' 
 #   #name_logcoef_den   ='logcoef1' 
 #   #name_logcoef_ ='logcoef1' 

 #   name_power_num     ='dijet_order1_%s' %label_num  
 #   name_power_den     ='dijet_order1_%s' %label_den  
 #   name_power_   ='dijet_order1_%s' %label_shape

 #   name_logcoef_num   ='dijet_order2_%s' %label_num  
 #   name_logcoef_den   ='dijet_order2_%s' %label_den  
 #   name_logcoef_ ='dijet_order2_%s' %label_shape

 #   val_power_num     = ws.var( name_power_num      )
 #   val_power_den     = ws.var( name_power_den      )
 #   val_power_   = ws.var( name_power_shape    )

 #   val_logcoef_num   = ws.var( name_logcoef_num    )
 #   val_logcoef_den   = ws.var( name_logcoef_den    )
 #   val_logcoef_ = ws.var( name_logcoef_shape  )

 #   power_pred    = ROOT.RooRealVar( power_pred_name   , 'power'  , (val_power_num.getValV() + val_power_.getValV() - val_power_den.getValV()) , -100, 100)
 #   logcoef_pred  = ROOT.RooRealVar( logcoef_pred_name , 'logcoef', (val_logcoef_num.getValV() + val_logcoef_.getValV() - val_logcoef_den.getValV()), -10, 10 )
 #   #power_ratio   = ROOT.RooRealVar( power_ratio_name  , 'power'  , val_power_num - val_power_den , -100, 100)
 #   #logcoef_ratio = ROOT.RooRealVar( logcoef_ratio_name, 'logcoef', val_logcoef_num - val_logcoef_den, -10, 10 )


 #   func = 'TMath::Power(@0/13000, @1+@2*TMath::Log10(@0/13000))'  
 #   prediction = ROOT.RooGenericPdf('dijet_prediction_%s' %suffix , 'prediction', func, ROOT.RooArgList(xvar,power_pred, logcoef_pred))
 #   norm_pred = ROOT.RooRealVar( 'dijet_prediction_%s_norm' %(suffix), 'prediction normalization', (hist_integral_ * hist_integral_num) / hist_integral_den )
 #   getattr( ws , 'import' ) ( norm_pred )
 #   getattr( ws , 'import' ) ( prediction )

 #   #---------------------------------------
    # Save the ratio hisotgram
    #---------------------------------------
    can_ratio = ROOT.TCanvas( str(uuid.uuid4()), '' )

   # ratiohist = ROOT.TH1F( 'ratio_%s' %(suffix), 'ratio', hist_num.GetNbinsX(), hist_num.GetXaxis().GetXmin(), hist_num.GetXaxis().GetXmax())
   # hist_num.Copy( ratiohist) 
   # #ratiohist = hist_num.Clone( 'ratio_%s' %(suffix) )
   # ROOT.SetOwnership( ratiohist, False )
   # ratiohist.Divide( hist_den )
   # ratiohist.SetMarkerStyle(20)
   # ratiohist.SetMarkerSize(1)
   # ratiohist.GetXaxis().SetTitle( 'Transverse Mass [GeV]' )
   # ratiohist.GetYaxis().SetTitle( 'ratio of passing to failing %s' %_var )
   # ratiohist.Draw()
   # ratiohist.SetStats(0)
   # ratiohist.SetMinimum( 0) 
   # ratiohist.SetMaximum( 5) 

   # ratio_func = ROOT.TF1( 'ratio_func', '( [2]*TMath::Power(x/13000, [0] + [1]*TMath::Log10(x/13000) ) ) ', xmin, xmax )

   # ROOT.SetOwnership( ratio_func, False )
   # ratio_func.SetParameter(0, result_num[name_power_num].n - result_den[name_power_den].n)
   # ratio_func.SetParameter(1, result_num[name_logcoef_num].n  - result_den[name_logcoef_den].n )
   # ratio_func.SetParameter(2, norm_num / norm_den )

   # ratio_func.Draw('same')

   # ratio_power_v = val_power_num.getValV() - val_power_den.getValV()
   # ratio_power_e = math.sqrt(val_power_num.getErrorHi()*val_power_num.getErrorHi() + val_power_den.getErrorHi()*val_power_den.getErrorHi() )
   # ratio_logcoef_v = val_logcoef_num.getValV() - val_logcoef_den.getValV()
   # ratio_logcoef_e = math.sqrt(val_logcoef_num.getErrorHi()*val_logcoef_num.getErrorHi() + val_logcoef_den.getErrorHi()*val_logcoef_den.getErrorHi() )

   # power_tex   = ROOT.TLatex(0, 0, 'power = %.01f #pm %.02f' %( ratio_power_v,  ratio_power_e ))
   # ROOT.SetOwnership( power_tex, False )
   # logcoef_tex = ROOT.TLatex(0, 0, 'logcoef = %.01f #pm %.02f' %( ratio_logcoef_v,  ratio_logcoef_e ))
   # ROOT.SetOwnership( logcoef_tex, False )

   # power_tex.SetNDC()
   # logcoef_tex.SetNDC()

   # power_tex  .SetX( 0.6 )
   # power_tex  .SetY( 0.84 )
   # logcoef_tex.SetX( 0.6 )
   # logcoef_tex.SetY( 0.78 )

   # power_tex.Draw()
   # logcoef_tex.Draw()

   # if closure :
   #     sampMan.outputs['wjetsclosure_ratio_%s' %suffix] = can_ratio 
   # else :
   #     sampMan.outputs['wjets_ratio_%s' %suffix] = can_ratio 
   # 
   # if closure :

   #     hist_sr    = clone_sample_and_draw( sampMan, sample, plot_var, full_sel_sr   , binning )
   #     fitMan_sr = FitManager( 'dijet', 2, samp.name, hist_sr, plot_var, eta_cut, xvar, full_suffix, useRooFit=False )
   #                                 
   #     result_sr  = fitMan_sr.fit_histogram( ws )

   #     can_sr = ROOT.TCanvas( str(uuid.uuid4()), '' )

   #     pred_val = (hist_integral_ * hist_integral_num) / hist_integral_den 

   #     tot_sr = hist_sr.Integral( hist_sr.FindBin( xmin ), hist_sr.FindBin( xmax ) )
   #     print 'SR integral = ', tot_sr
   #     sr_func = ROOT.TF1( 'sr_func', '( [2]*TMath::Power(x/13000, [0] + [1]*TMath::Log10(x/13000) ) ) ', xmin, xmax )
   #     sr_func.SetParameter(0, result_num['power'].n + result_['power'].n - result_den['power'].n )
   #     sr_func.SetParameter(1, result_num['logcoef'].n + result_['logcoef'].n - result_den['logcoef'].n )
   #     sr_func.SetParameter(2, 1 )
   #     sr_int = sr_func.Integral( xmin, xmax )
   #     print 'Normalization = ',(norm_*norm_num) / norm_den
   #     print 'func Integral SR Before = ', sr_func.Integral( xmin, xmax )
   #     sr_func.SetParameter(2, pred_val/sr_int )

   #     hist_sr.Draw()
   #     sr_func.Draw('same')

   #     sampMan.outputs['wjetsclosure_pred_%s' %suffix] = can_sr


   #     tot_ratio = result_num['integral']/result_den['integral']


def clone_sample_and_draw( sampMan, samp, var, sel, binning ) :

    newSamp = sampMan.clone_sample( oldname=samp, newname=samp+str(uuid.uuid4()), temporary=True ) 
    sampMan.create_hist( newSamp, var, sel, binning )
    return newSamp.hist

main()

    




    
    




