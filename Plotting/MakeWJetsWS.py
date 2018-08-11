import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import re
import numpy as np
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

parser.add_argument('--baseDirMuG',		 default=None,			 dest='baseDirMuG',			required=False, help='Path to muon base directory')
parser.add_argument('--baseDirElG',		 default=None,			 dest='baseDirElG',			required=False, help='Path to electron base directory')
parser.add_argument('--outputDir',		default=None,			dest='outputDir',		  required=False, help='Output directory to write histograms')
parser.add_argument('--useRooFit',		 default=False,    action='store_true',		 dest='useRooFit', required=False, help='Make fits using roostats' )
parser.add_argument('--doClosure',		 default=False,   action='store_true',		 dest='doClosure', required=False, help='make closure tests' )

options = parser.parse_args()

_TREENAME = 'UMDNTuple/EventTree'
_FILENAME = 'tree.root'
_XSFILE   = 'cross_sections/photon15.py'
_LUMI	  = 36000
_BASEPATH = '/home/jkunkle/usercode/Plotting/LimitSetting/'
_SAMPCONF = 'Modules/Resonance.py'


#def get_cut_defaults( _var, ieta ) :
#
#	 cut_defaults = {'sigmaIEIE' : { 'EB' : ( 0.012, 0.02 ), 'EE' : ( 0.0, 0.0 ) },
#					 'chIso'	 : { 'EB' : ( 4, 10 ),		 'EE' : (4, 10) },
#					}
#
#	 return cut_defaults[_var][ieta]


ROOT.gROOT.SetBatch(True)
if options.outputDir is not None :
#	 ROOT.gROOT.SetBatch(True)
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
	#sel_base_el = 'ph_n>=1 && el_n>=1'
	sel_base_el = '1'
	sel_base_el = 'ph_n>=1 && el_n==1 &&ph_pt[0]>50'
	sel_base_el = 'ph_n>=1 && el_n==1 '

	#eta_cuts = ['EB', 'EE']
	eta_cuts = ['EB']

	workspaces_to_save = {}

	selections = { 'base'	 : { 
							   # 'mu' : {'selection' : sel_base_mu }, 
								'el' : { 'selection' : sel_base_el }, 
							   },
				  }

	elefake			   = ROOT.RooWorkspace( 'elefake' )

	#make_efake( sampManElG, 'Z+jets', sel_base_el,'EB', 'ph_eta', suffix='noOvLrm', closure=False, workspace=elefake, overlaprm=0)


	f1 = ROOT.TFile("%s/output.root"%(options.outputDir),"RECREATE")

	#undo scaling by cross section
	#print "scale ", sampManElG.get_samples(name="DYJetsToLL_M-50")[0].scale
#	sampManElG.get_samples(name="DYJetsToLL_M-50")[0].scale=1.
	 #closure( sampManElG, 'Z+jets', sel_base_el,'EB', plot_var = 'm_lep_ph',var = "ph_eta",varbins=np.linspace(-2.5,2.5,51),xtitle="photon #eta")
#	closure( sampManElG, 'Z+jets', sel_base_el,'EB', plot_var = 'm_lep_ph',varbins=[-2.4,-1.8,-1.4,-1,-0.6,-0.3,0,0.3,0.6,1,1.4,1.8,2.4], var="ph_eta")
	closure( sampManElG, 'Z+jets', sel_base_el,'EB', plot_var = 'm_lep_ph',varbins=range(0,200,20))
#	closure( sampManElG, 'Wgamma', sel_base_el,'EB', plot_var = 'm_lep_ph',varbins=[0,50,250])

	#elefake.writeToFile( '%s/%s.root' %( options.outputDir,elefake.GetName() ) )
	#for key, can in sampManElG.outputs.iteritems() :
	f1.Write()
	f1.Close()



def make_efake( sampMan, sample, sel_base, eta_cut, plot_var, suffix='', workspace=None, overlaprm=0):
	#---------------------------------------
	# Get the base selection for each region
	#---------------------------------------
	passev = 'ph_mediumPassCSEV_n==1'
	failev = 'ph_mediumFailCSEV_n==1'
	ph_selection_sr  = 'met_pt>40 && ' + passev #'%s==1' %defs.get_phid_selection('medium')
	ph_selection_B	 = 'met_pt<40 && ' +failev #'%s==1'%defs.get_phid_selection( num_var, _var )
	ph_selection_A	 = 'met_pt<40 && '+passev #'%s==1' %defs.get_phid_selection( num_var )
	ph_selection_D	 = 'met_pt>40 && '+failev #'%s==1' %defs.get_phid_selection( _var )
	varfail = '[ptSorted_ph_mediumPassEleOlapFailCSEV_idx[0]]'
	varpass = '[ptSorted_ph_mediumPassEleOlapPassCSEV_idx[0]]'

	if overlaprm:
		ph_selection_sr  = 'met_pt>40 && ph_passEleVeto[0]==1'
		ph_selection_B	 = 'met_pt<40 && ph_passEleVeto[0]==0'
		ph_selection_A	 = 'met_pt<40 && ph_passEleVeto[0]==1'
		ph_selection_D	 = 'met_pt>40 && ph_passEleVeto[0]==0'
		varfail=varpass=''

	full_sel_D = ' && '.join( [sel_base, ph_selection_D, ] )
	full_sel_A	 = ' && '.join( [sel_base, ph_selection_A,] )
	full_sel_B	 = ' && '.join( [sel_base, ph_selection_B,] )
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
	binning = (160,-4,4)
	#binning = (200,0,200)
	hist_D	 = clone_sample_and_draw( sampMan, sample, plot_var+varfail , full_sel_D, binning )
	print hist_D
	hist_A	 = clone_sample_and_draw( sampMan, sample, plot_var+varpass , full_sel_A  , binning )
	hist_B	 = clone_sample_and_draw( sampMan, sample, plot_var+varfail , full_sel_B  , binning )
	hist_sr  = clone_sample_and_draw( sampMan, sample, plot_var+varpass , full_sel_sr  , binning )
	c1	= ROOT.TCanvas('c1','c1')
	c1.SetGridx()
	#c1.SetLogy()
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


def make_efake2( sampMan, sample, sel_base, eta_cut, plot_var, suffix='',  workspace=None) :

	#only single plot
	#---------------------------------------
	# Get the base selection for each region
	#---------------------------------------
	ph_selection_A	 = 'met_pt<40 && ph_passEleVeto[0]==1' #'%s==1' %defs.get_phid_selection( num_var )
	full_sel_A	 = ' && '.join( [sel_base, ph_selection_A,] )
	label_A   = 'a_%s_' %suffix

	#---------------------------------------
	# draw the histograms
	#---------------------------------------
	binning = (200,-5,5)
	hist_A	 = clone_sample_and_draw( sampMan, sample, plot_var , full_sel_A  , binning )
	c1	= ROOT.TCanvas('c1','c1')
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
	c1.SetGridx()
	c1.SetLogy()
	ROOT.gStyle.SetOptStat(0)
	c1.Divide(2,2,0,0)
	
	c1.cd(1)
	ROOT.gPad.SetTickx(2)
	ROOT.gPad.SetGridx()
	h1.Draw()
	
	c1.cd(2)
	ROOT.gPad.SetGridx()
	ROOT.gPad.SetTickx(2)
	ROOT.gPad.SetTicky(2)
	h2.GetYaxis().SetLabelOffset(0.01)
	h2.Draw()
	
	c1.cd(3)
	ROOT.gPad.SetGridx()
	h3.Draw()
	
	c1.cd(4)
	ROOT.gPad.SetGridx()
	ROOT.gPad.SetTicky(2)
	h4.Draw()
	c1.SaveAs("hist_multi_"+suffix+".pdf")


	#make_efake( sampManElG, 'Z+jets', sel_base_el,'EB', 'ph_eta', suffix='noOvLrm', workspace=elefake, overlaprm=0)


def closure( sampMan, sample, sel_base, eta_cut, var="ph_pt", varbins = range(0,200,25),xtitle="photon p_{T}", plot_var='m_lep_ph', suffix='', workspace=None, overlaprm=1) :
	bins = np.array(varbins,'f')
	h_sr_est= ROOT.TH1F("h_sr_est","",len(varbins)-1,bins)
	h_sr_sum= ROOT.TH1F("h_sr_sum",";%s;count" %(xtitle),len(varbins)-1,bins)
	for i in range(len(varbins)-1):
		sel_base1= sel_base + "&& %s[0]>%3f && %s[0]<=%3f" %(var,varbins[i],var,varbins[i+1])
		filename=var+"%s_sel_%s_%i_%i_%i_" %(plot_var,var,i,varbins[i],varbins[i+1])
		sr_est, sr_est_e, sr_sum, sr_sum_e = closureBin(sampMan, sample, sel_base1, eta_cut, plot_var, overlaprm,filename) 
		#sr_est, sr_est_e, sr_sum = closureBin(sampMan, sample, sel_base1, eta_cut, plot_var, overlaprm,filename) 
		h_sr_sum.SetBinContent(i+1,sr_sum)
		h_sr_sum.SetBinError(i+1,sr_sum_e)
		h_sr_est.SetBinContent(i+1,sr_est)
		h_sr_est.SetBinError(i+1,sr_est_e)
	# draw graph and formatting
	c2 = ROOT.TCanvas("c2","multipads",600,700)
	toppad = ROOT.TPad("top","top",0.01,0.35,0.99,0.99)
	botpad = ROOT.TPad("bot","bot",0.01,0.01,0.99,0.34)
	toppad.SetTopMargin(0.08)
	toppad.SetBottomMargin(0.06)
	toppad.SetLeftMargin(0.15)
	toppad.SetRightMargin(0.05)
	botpad.SetTopMargin(0.00)
	botpad.SetBottomMargin(0.3)
	botpad.SetLeftMargin(0.15)
	botpad.SetRightMargin(0.05)
	c2.cd()
	toppad.Draw()
	#toppad.SetLogy(1)
	botpad.Draw()
	toppad.cd()
	hratio = doratio(h_sr_sum,h_sr_est)
	hformat(h_sr_est,ROOT.kBlack)
	hformat(h_sr_sum,ROOT.kRed)
	h_sr_sum.Draw("e")
	h_sr_est.Draw("e same")
	tl = ROOT.TLegend(0.7,0.72,0.88,0.88,"","brNDC")
	tl.AddEntry("h_sr_sum","signal","PL")
	tl.AddEntry("h_sr_est","A/B*D","PL")
	tl.Draw()
	tex = cmsinternal(toppad)
	# bottom pad
	botpad.cd()
	hratio.Draw("e")
	oneline = ratioline(hratio)
	c2.SaveAs("test.pdf")

def closureBin( sampMan, sample, sel_base, eta_cut, plot_var='m_lep_ph', overlaprm=1,filename = "m_lep_ph_default_", doSave=0) :
	if eta_cut=="EB": sel_base += '&&abs(ph_eta[0])<=2'
	passev = 'ph_mediumPassCSEV_n==1'
	failev = 'ph_mediumFailCSEV_n==1'
	ph_selection_sr  = 'met_pt>25 && ' + passev
	ph_selection_B	 = 'met_pt<25 && ' +failev
	ph_selection_A	 = 'met_pt<25 && ' +passev
	ph_selection_D	 = 'met_pt>25 && ' +failev 
	varfail = '[ptSorted_ph_mediumPassEleOlapFailCSEV_idx[0]]'
	varpass = '[ptSorted_ph_mediumPassEleOlapPassCSEV_idx[0]]'

	if overlaprm:
		ph_selection_sr  = 'met_pt>25 && ph_passEleVeto[0]==1'
		ph_selection_B	 = 'met_pt<25 && ph_passEleVeto[0]==0'
		ph_selection_A	 = 'met_pt<25 && ph_passEleVeto[0]==1'
		ph_selection_D	 = 'met_pt>25 && ph_passEleVeto[0]==0'
		varfail=varpass=''

	full_sel_D = ' && '.join( [sel_base, ph_selection_D, ] )
	full_sel_A	 = ' && '.join( [sel_base, ph_selection_A,] )
	full_sel_B	 = ' && '.join( [sel_base, ph_selection_B,] )
	full_sel_sr    = ' && '.join( [sel_base, ph_selection_sr,] )

	binning = (200,0,200)
	c1 =ROOT.TCanvas()
	hist_D	 = clone_sample_and_draw( sampMan, sample, plot_var+varfail , full_sel_D, binning )
	hist_D.Draw('e')
	if (doSave):c1.SaveAs(filename+"D.pdf")
	hist_A	 = clone_sample_and_draw( sampMan, sample, plot_var+varpass , full_sel_A  , binning )
	hist_A.Draw('e')
	if (doSave):c1.SaveAs(filename+"A.pdf")
	hist_B	 = clone_sample_and_draw( sampMan, sample, plot_var+varfail , full_sel_B  , binning )
	hist_B.Draw('e')
	if (doSave):c1.SaveAs(filename+"B.pdf")
	hist_sr  = clone_sample_and_draw( sampMan, sample, plot_var+varpass , full_sel_sr  , binning )
	hist_sr.Draw('e')
	if (doSave): c1.SaveAs(filename+"S.pdf")
	sampMan.outputs[filename+"A"] = hist_A
	sampMan.outputs[filename+"B"] = hist_B
	sampMan.outputs[filename+"S"] = hist_sr
	sampMan.outputs[filename+"D"] = hist_D
	intrange = [80,100]
	a_sum = hist_A.Integral(*intrange)
	b_sum = hist_B.Integral(*intrange)
	sr_sum = hist_sr.Integral(*intrange)
	d_sum = hist_D.Integral(*intrange)
	print  "Region A: ", a_sum
	print  "Region B: ", b_sum
	print  "Region signal: ", sr_sum
	print  "Region D: ", d_sum
	if b_sum>0: sr_est = a_sum/b_sum*d_sum
	else: return 0,0,sr_sum,0
	if (sr_est>0): sr_err = sr_est*math.sqrt(1/a_sum+ 1/b_sum+1/d_sum)
	else: sr_err =0
	return sr_est, sr_err, sr_sum ,math.sqrt(sr_sum)

#	if b_sum>0: r_est = a_sum/b_sum
#	else: r_est = 0 
#	if d_sum>0: r_sr = sr_sum/d_sum
#	else: r_sr = 0 
#	if (r_est>0): r_err_est = math.sqrt(a_sum *b_sum)/(a_sum+b_sum)    # r_est*math.sqrt(1/a_sum+ 1/b_sum)
#	else: r_err_est =0
#	if (r_sr>0):  r_err_sr  = math.sqrt(sr_sum*d_sum)/(sr_sum+d_sum)  # r_sr*math.sqrt(1/sr_sum+ 1/d_sum)
#	else: r_err_sr =0
#	return r_est, r_err_est, r_sr, r_err_sr 

def doratio(h1,h2):
	hratio = h1.Clone("hratio")
	hratio.Divide(h2)
	hratio.SetMarkerStyle(20)
	hratio.SetMarkerSize(1.1)
	hratio.SetStats(0)
	hratio.SetTitle("")
	hratio.GetYaxis().SetTitle("ratio")
	hratio.SetLineColor(ROOT.kBlack)
	hratio.SetLineWidth(2)
	hratio.GetYaxis().SetTitleSize(0.10)
	hratio.GetYaxis().SetTitleOffset(0.6)
	hratio.GetYaxis().SetLabelSize(0.10)
	hratio.GetXaxis().SetLabelSize(0.10)
	hratio.GetXaxis().SetTitleSize(0.10)
	hratio.GetXaxis().SetTitleOffset(1.0)
	hratio.GetYaxis().SetRangeUser(0.5,1.5)
	hratio.GetYaxis().CenterTitle()
	hratio.GetYaxis().SetNdivisions(506, True)
	return hratio

def cmsinternal(pad):
	pad.cd()
	tex = ROOT.TLatex(0.18,0.93,"CMS Internal")
	tex.SetNDC()
	tex.SetTextSize(0.05)
	tex.SetLineWidth(2)
	tex.Draw()
	return tex

def ratioline(hratio):
	left_edge  = hratio.GetXaxis().GetXmin()
	right_edge = hratio.GetXaxis().GetXmax()
	
	oneline = ROOT.TLine(left_edge, 1, right_edge, 1)
	oneline.SetLineStyle(3)
	oneline.SetLineWidth(2)
	oneline.SetLineColor(ROOT.kBlack)
	oneline.Draw()
	return oneline

def hformat(h1,color):
	h1.SetLineColor(color)
	h1.SetMarkerColor(color)
	h1.SetMarkerStyle(20)
	h1.SetMarkerSize(1.1)
	h1.SetStats(0)
	h1.SetLineWidth(2)
	h1.GetYaxis().SetTitleSize(0.05)
	h1.GetYaxis().SetTitleOffset(1.15)
	h1.GetYaxis().SetLabelSize(0.05)
	h1.GetXaxis().SetLabelSize(0.05)
	h1.GetXaxis().SetTitleSize(0.05)
	h1.GetXaxis().SetTitleOffset(0.8)





main()

	




	
	




