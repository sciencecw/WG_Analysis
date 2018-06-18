from core import Filter
import inspect
import sys
import os

def get_remove_filter() :
    """ Define list of regex strings to filter input branches to remove from the output.
        Defining a non-empty list does not apply the filter, 
        you must also supply --enableRemoveFilter on the command line.
        If both filters are used, all branches in keep_filter are used
        except for those in remove_filter """

    return ['']

def get_keep_filter(tag=None) :
    """ Define list of regex strings to filter input branches to retain in the output.  
        Defining a non-empty list does not apply the filter, 
        you must also supply --enableKeepFilter on the command line
        If both filters are used, all branches in keep_filter are used
        except for those in remove_filter """

    mu_basic = ['mu_n','mu_pt', 'mu_eta', 'mu_phi', 'mu_e', 'mu_charge']
    el_basic = ['el_n', 'el_phi', 'el_eta', 'el_e', 'el_pt', 'el_charge']
    el_addtl = ['el_phiOrig', 'el_sc_eta', 'el_etaOrig', 'el_eOrig', 'el_ptOrig', 
                'el_passVIDHEEP', 'el_passVIDHLT', 'el_passVIDTight', 'el_passVIDVeryLoose', 
                'el_passVIDLoose', 'el_passConvVeto', 'el_passVIDMedium']
    ph_basic = ['ph_n', 'ph_phi', 'ph_eta', 'ph_pt','ph_e','ph_hasPixSeed', 'ph_passEleVeto.*', ]
    ph_addtl = ['ph_passVIDLoose', 'ph_passVIDMedium', 'ph_passVIDTight', 
                'ph_sc_phi', 'ph_sc_eta', 'ph_neuIsoCorr', 'ph_phiOrig', 'ph_etaOrig', 'ph_phiWidth', 'ph_ptOrig', 
                'ph_sigmaIEIEFull5x5', 'ph_r9', 'ph_etaWidth', 'ph_eOrig', 'ph_r9Full5x5', 'ph_sigmaIEIE', 
                'ph_chIsoCorr', 'ph_phoIso', 'ph_chIso', 'ph_neuIso', 'ph_hOverE',]

    met_basic = ['met_pt', 'met_phi']
    met_addtl = ['met_UnclusteredEnUp_phi', 'met_UnclusteredEnDown_phi', 'met_UnclusteredEnDown_pt', 'met_UnclusteredEnUp_pt', 
                 'met_PhotonEnUp_phi', 'met_PhotonEnDown_phi', 'met_PhotonEnDown_pt', 'met_PhotonEnUp_pt', 
                 'met_JetEnUp_phi', 'met_JetEnDown_phi', 'met_JetEnUp_pt', 'met_JetEnDown_pt', 
                 'met_ElectronEnUp_phi', 'met_ElectronEnDown_phi', 'met_ElectronEnDown_pt', 'met_ElectronEnUp_pt', 
                 'met_JetResDown_phi', 'met_JetResUp_phi', 'met_JetResDown_pt', 'met_JetResUp_pt', 
                 'met_MuonEnDown_phi', 'met_MuonEnUp_phi', 'met_MuonEnUp_pt', 'met_MuonEnDown_pt', ]

    jet_basic = ['jet_n', 'jet_pt', 'jet_eta', 'jet_phi', 'jet_e',]
    jet_addtl = ['jet_CSVLoose_n', 'jet_CSVMedium_n', 'jet_CSVTight_n']

    event_basic = ['rho', 'pu_n', 'truepu_n', 'vtx_n', 'pdf_id1', 'pdf_id2', 'pdf_scale', 'pdf_x2', 'pdf_x1', 
                   'lumiSection', 'eventNumber', 'runNumber', 'bxNumber', 'isData', 'EventWeights']


    branches_tight = mu_basic + el_basic + ph_basic + met_basic + jet_basic + event_basic + ["gen_.*"]

    if tag == 'tight' : 
        return branches_tight
    else :
        return branches_tight + el_addtl + ph_addtl + met_addtl + jet_addtl

def config_analysis( alg_list, args ) :
    """ Configure analysis modules. Order is preserved """

    # run on the function provided through
    # the args
    for s in inspect.getmembers(sys.modules[__name__]) :
        if s[0] == args['function'] :
            print '*********************************'
            print 'RUN %s' %( args['function'] )
            print '*********************************'
            s[1]( alg_list, args )

    alg_list.append( weight_event( args ) )

def make_final_mumu( alg_list, args) :

    mu_pt = args.get( 'mu_pt', ' > 25 ' )

    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon( mu_pt) )
    alg_list.append( filter_electron( ) )
    alg_list.append( filter_photon( ) )
    alg_list.append( filter_jet( ) )

    filter_trig = filter_trigger()
    filter_trig.cut_bits = ' == 23 | == 31 '
    alg_list.append( filter_trig )
 
    filter_event = Filter('FilterEvent')
    filter_event.cut_mu_n = ' == 2 '

    alg_list.append( filter_event )

    #alg_list.append( Filter( 'MakePhotonCountVars' ) )
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

def make_final_elel( alg_list, args) :

    el_pt = args.get( 'el_pt', ' > 25 ' )

    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon( ) )
    alg_list.append( filter_electron( el_pt ) )
    alg_list.append( filter_photon( ) )
    alg_list.append( filter_jet( ) )

    filter_trig = filter_trigger()
    filter_trig.cut_bits = ' == 60 '
    alg_list.append( filter_trig )

    filter_event = Filter('FilterEvent')
    filter_event.cut_el_n = ' == 2 '

    alg_list.append( filter_event )

    #alg_list.append( Filter( 'MakePhotonCountVars' ) )
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

def make_final_muel( alg_list, args) :

    mu_pt = args.get( 'mu_pt', ' > 25 ' )
    el_pt = args.get( 'el_pt', ' > 25 ' )

    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon( mu_pt ) )
    alg_list.append( filter_electron(el_pt ) )
    alg_list.append( filter_photon( ) )
    alg_list.append( filter_jet( ) )
    alg_list.append( filter_trigger() )

    filter_event = Filter('FilterEvent')
    filter_event.cut_mu_n = ' == 1 '
    filter_event.cut_el_n = ' == 1 '

    alg_list.append( filter_event )

    #alg_list.append( Filter( 'MakePhotonCountVars' ) )
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

def make_final_mu( alg_list, args) :

    mu_pt = args.get( 'mu_pt', ' > 25 ' )
    ph_id_cut = args.get( 'ph_id_cut', 'medium' )

    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon( mu_pt ) )
    alg_list.append( filter_electron( ) )
    alg_list.append( filter_photon( id_cut=ph_id_cut ) )
    alg_list.append( filter_jet( ) )

    filter_trig = filter_trigger()
    filter_trig.cut_bits = ' == 23 | == 31 '
    alg_list.append( filter_trig )

    filter_event = Filter('FilterEvent')
    filter_event.cut_mu_n = ' == 1 '
    filter_event.cut_mu_pt30_n = ' == 1 '

    alg_list.append( filter_event )

    #alg_list.append( Filter( 'MakePhotonCountVars' ) )
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

def make_final_el( alg_list, args) :

    el_pt = args.get( 'el_pt', ' > 25 ' )
    ph_id_cut = args.get( 'ph_id_cut', 'medium' )

    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon( ) )
    alg_list.append( filter_electron( el_pt ) )
    alg_list.append( filter_photon( id_cut=ph_id_cut ) )
    alg_list.append( filter_jet( ) )

    filter_trig = filter_trigger()
    filter_trig.cut_bits = ' == 60 '
    alg_list.append( filter_trig )

    filter_event = Filter('FilterEvent')
    filter_event.cut_el_n = ' == 1 '
    filter_event.cut_el_pt30_n = ' == 1 '

    alg_list.append( filter_event )

    #alg_list.append( Filter( 'MakePhotonCountVars' ) )
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

def make_final_elg( alg_list, args) :

    mu_pt = args.get( 'mu_pt', ' > 10 ' )
    el_pt = args.get( 'el_pt', ' > 10 ' )
    ph_pt = args.get( 'ph_pt', ' > 15 ' )
    phot_vars = args.get( 'phot_vars', 'False' )
    phot_id = args.get( 'phot_id', 'medium' )
    ph_eta = args.get( 'ph_eta', None )
    sec_lep_veto = args.get( 'sec_lep_veto', 'True' )
    unblind = args.get( 'unblind', 'False' )
    eleVeto = args.get('eleVeto', 'None' )
    eleOlap = args.get('eleOlap', 'True' )


    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon(mu_pt ) )
    alg_list.append( filter_electron(el_pt ) )
    alg_list.append( filter_photon( ph_pt, id_cut=phot_id, ieta_cut=ph_eta,ele_veto=eleVeto, ele_olap=eleOlap  ) )
    alg_list.append( filter_jet( ) )

    filter_trig = filter_trigger()
    filter_trig.cut_bits = ' == 60 '
    alg_list.append( filter_trig )


    filter_event = Filter('FilterEvent')
    if eleOlap == 'False' :
        filter_event.cut_el_pt30_n = ' > 0 '
        filter_event.cut_ph_n = ' > 0 '
    else :
        filter_event.cut_el_pt30_n = ' == 1 '
        filter_event.cut_ph_n = ' > 0 '
        if sec_lep_veto != 'False' :
            filter_event.cut_el_n = ' == 1 '
            filter_event.cut_mu_n = ' == 0 '

    alg_list.append( filter_event )

    if phot_vars == 'True' :
        alg_list.append( Filter( 'MakePhotonCountVars' ) )

    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

    if unblind is not 'True' :
        filter_blind = Filter( 'FilterBlind' )
        filter_blind.cut_mt_lep_met_ph = ' < 100 ' 

        filter_blind.add_var( 'isData', args.get('isData', ' == False' ) )
        alg_list.append( filter_blind )

def make_final_mug( alg_list, args) :

    el_pt = args.get( 'el_pt', ' > 10 ' )
    mu_pt = args.get( 'mu_pt', ' > 10 ' )
    ph_pt = args.get( 'ph_pt', ' > 15 ' )
    muphtrig = args.get( 'muphtrig', 'False' )
    phot_vars = args.get( 'phot_vars', ' False ' )
    phot_id = args.get( 'phot_id', 'medium' )
    ph_eta = args.get( 'ph_eta', None )
    sec_lep_veto = args.get( 'sec_lep_veto', 'True' )
    unblind = args.get( 'unblind', 'False' )

    # order should be muon, electron, photon, jet
    alg_list.append( filter_muon( mu_pt ) )
    alg_list.append( filter_electron( el_pt ) )
    alg_list.append( filter_photon( ph_pt, id_cut=phot_id, ieta_cut=ph_eta ) )
    alg_list.append( filter_jet( ) )

    filter_trig = filter_trigger()
    filter_trig.cut_bits = ' == 23 | == 31 '
    alg_list.append( filter_trig )

    filter_event = Filter('FilterEvent')
    filter_event.cut_mu_pt30_n = ' == 1 '
    filter_event.cut_ph_n = ' > 0 '

    if sec_lep_veto is not 'False' :
        filter_event.cut_mu_n = ' == 1 '
        filter_event.cut_el_n = ' == 0 '

    alg_list.append( filter_event)


    if phot_vars == 'True' :
        alg_list.append( Filter( 'MakePhotonCountVars' ) )
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( build_truth(args) )

    if unblind is not 'True' :
        filter_blind = Filter( 'FilterBlind' )
        filter_blind.cut_mt_lep_met_ph = ' < 100 ' 

        filter_blind.add_var( 'isData', args.get('isData', ' == False' ) )
        alg_list.append( filter_blind )

#def make_final_elgjj( alg_list, args) :
#
#    el_pt = args.get( 'el_pt', ' > 25 ' )
#    ph_pt = args.get( 'ph_pt', ' > 15 ' )
#
#    # order should be muon, electron, photon, jet
#    alg_list.append( filter_muon( ) )
#    alg_list.append( filter_electron(el_pt ) )
#    alg_list.append( filter_photon( ph_pt ) )
#    alg_list.append( filter_jet( ) )
#
#    filter_event = Filter('FilterEvent')
#    filter_event.cut_el_n = ' == 1 '
#    filter_event.cut_ph_n = ' == 1 '
#    filter_event.cut_jet_n = ' > 1 '
#    filter_event.cut_trig_Ele27_eta2p1_tight = ' == True '
#
#    alg_list.append( filter_event )
#
#    alg_list.append( Filter( 'BuildEventVars' ) )
#    alg_list.append( Filter( 'BuildTruth' ) )
#
#    filter_blind = Filter( 'FilterBlind' )
#    filter_blind.cut_abs_dijet_m_from_z = ' < 15 ' 
#
#    filter_blind.add_var( 'isData', args.get('isData', ' == False' ) )
#    alg_list.append( filter_blind )
#
#def make_final_mugjj( alg_list, args) :
#
#    mu_pt = args.get( 'mu_pt', ' > 25 ' )
#    ph_pt = args.get( 'ph_pt', ' > 15 ' )
#
#    # order should be muon, electron, photon, jet
#    alg_list.append( filter_muon( mu_pt ) )
#    alg_list.append( filter_electron( ph_pt ) )
#    alg_list.append( filter_photon( ) )
#    alg_list.append( filter_jet( ) )
#
#    filter_event = Filter('FilterEvent')
#    filter_event.cut_mu_n  = ' == 1 '
#    filter_event.cut_ph_n  = ' == 1 '
#    filter_event.cut_jet_n = ' > 1 '
#    filter_event.cut_trig_Mu24_IsoORIsoTk = ' == True '
#
#    alg_list.append( filter_event )
#
#    alg_list.append( Filter( 'BuildEventVars' ) )
#    alg_list.append( Filter( 'BuildTruth' ) )
#
#    filter_blind = Filter( 'FilterBlind' )
#    filter_blind.cut_abs_dijet_m_from_z = ' < 15 ' 
#
#    filter_blind.add_var( 'isData', args.get('isData', ' == False' ) )
#    alg_list.append( filter_blind )
#
def make_nofilt( alg_list, args ) :

    pass_lepton = args.get('pass_lepton', 'False')
    if pass_lepton == 'True' :
        alg_list.append( filter_muon( mu_pt = ' > 0 ' ) )
        alg_list.append( filter_electron( el_pt = ' > 0 ' ) )

    alg_list.append( filter_photon( ph_pt = ' > 0 ', id_cut = 'medium'   )  )
                           
    alg_list.append( Filter( 'BuildEventVars' ) )
    alg_list.append( Filter( 'BuildTruth' ) )


def filter_trigger() : 

    filter_trigger = Filter('FilterTrigger')

    filter_trigger.add_var( 'triggerBits', '23:HLT_IsoMu24,31:HLT_IsoTkMu24,60:HLT_Ele27_eta2p1_WPTight_Gsf' )

    return filter_trigger


def filter_muon( mu_pt = ' > 25 ', do_cutflow=False, apply_corrections=False, do_hists=False, evalPID='tight' ) :

    filt = Filter('FilterMuon')

    if do_cutflow :
        filt.do_cutflow = True
        filt.add_var('evalPID', evalPID )

    filt.cut_pt           = mu_pt
    filt.cut_eta       = ' < 2.5'
    filt.cut_tight       = ' == True '
    filt.add_var( 'triggerMatchBits', '23,31' )

    filt.cut_isPf_loose         = ' == True '
    filt.cut_isGlobalOrTk_loose = ' == True '

    filt.cut_isGlobal_tight   = ' == True '
    filt.cut_isPF_tight       = ' == True '
    filt.cut_abseta_tight     = ' < 2.4'
    filt.cut_chi2_tight       = ' < 10'
    filt.cut_nMuonHits_tight  = ' > 0 ' 
    filt.cut_nStations_tight  = ' > 1'
    filt.cut_nTrkLayers_tight = ' > 5 ' 
    filt.cut_nPixelHits_tight = ' > 0'
    filt.cut_d0_tight         = ' < 0.2'
    filt.cut_z0_tight         = ' < 0.5'
    filt.cut_corriso_tight    = ' < 0.25'
    filt.cut_trkiso_tight     = ' < 0.05 '

    if apply_corrections :
        filt.add_var( 'apply_corrections', 'true' )

        workarea = os.getenv('WorkArea')
        filt.add_var( 'path', '%s/TreeFilter/RecoWgg/data/MuScleFitCorrector_v4_3/MuScleFit_2012_MC_53X_smearReReco.txt' %workarea )

    if do_hists :
        filt.add_hist( 'cut_pt', 100, 0, 500 )
        filt.add_hist( 'cut_abseta', 50, 0, 5 )
        filt.add_hist( 'cut_chi2', 50, 0, 50 )
        filt.add_hist( 'cut_nTrkLayers', 20, 0, 20 )
        filt.add_hist( 'cut_nStations', 5, 0, 5 )
        filt.add_hist( 'cut_nPixelHits', 20, 0, 20 )
        filt.add_hist( 'cut_d0', 100, -0.05, 0.05 )
        filt.add_hist( 'cut_z0', 100, -0.05, 0.05 )
        filt.add_hist( 'cut_trkiso', 50, 0, 0.5 )
        filt.add_hist( 'cut_corriso', 50, 0, 0.5 )

    return filt

def filter_electron( el_pt = ' > 25 ', do_cutflow=False, do_hists=False, apply_corrections=False, evalPID='tight' ) :

    filt = Filter('FilterElectron')

    filt.cut_pt         = el_pt
    #filt.cut_eta        = ' < 2.5'
    filt.cut_abssceta       = ' <2.5 '

    #filt.cut_tight     = ' == True '
    filt.cut_vid_medium     = ' == True '
    filt.cut_muon_dr    = ' > 0.4 '
    filt.add_var( 'triggerMatchBits', '60' )


    filt.cut_sigmaIEIE_barrel_tight        = ' < 0.00998 '
    filt.cut_absdEtaIn_barrel_tight        = ' < 0.00308 '
    filt.cut_absdPhiIn_barrel_tight        = ' < 0.0816 '
    filt.cut_hovere_barrel_tight           = ' < 0.0414 '
    filt.cut_isoRho_barrel_tight           = ' < 0.0588 '
    filt.cut_ooEmooP_barrel_tight          = ' < 0.0129 '
    #filt.cut_d0_barrel_tight               = ' < 0.0111 '
    #filt.cut_z0_barrel_tight               = ' < 0.0466 '
    filt.cut_misshits_barrel_tight         = ' < 1 '
    filt.cut_passConvVeto_barrel_tight     = ' == 1 '

    filt.cut_sigmaIEIE_barrel_medium       = ' < 0.00998 '
    filt.cut_absdEtaIn_barrel_medium       = ' < 0.00311 '
    filt.cut_absdPhiIn_barrel_medium       = ' < 0.103 '
    filt.cut_hovere_barrel_medium          = ' < 0.253 '
    filt.cut_isoRho_barrel_medium          = ' < 0.0695 '
    filt.cut_ooEmooP_barrel_medium         = ' < 0.134 '
    #filt.cut_d0_barrel_medium              = ' < 0.0118 '
    #filt.cut_z0_barrel_medium              = ' < 0.373 '
    filt.cut_misshits_barrel_medium        = ' < 1 '
    filt.cut_passConvVeto_barrel_medium    = ' == 1 '

    filt.cut_sigmaIEIE_barrel_loose        = ' < 0.011 '
    filt.cut_absdEtaIn_barrel_loose        = ' < 0.00477 '
    filt.cut_absdPhiIn_barrel_loose        = ' < 0.222 '
    filt.cut_hovere_barrel_loose           = ' < 0.298 '
    filt.cut_isoRho_barrel_loose           = ' < 0.0994 '
    filt.cut_ooEmooP_barrel_loose          = ' < 0.241 '
    #filt.cut_d0_barrel_loose               = ' < 0.0261 '
    #filt.cut_z0_barrel_loose               = ' < 0.41 '
    filt.cut_misshits_barrel_loose         = ' < 1 '
    filt.cut_passConvVeto_barrel_loose     = ' == 1 '

    filt.cut_sigmaIEIE_barrel_veryloose    = ' < 0.0115 '
    filt.cut_absdEtaIn_barrel_veryloose    = ' < 0.00749 '
    filt.cut_absdPhiIn_barrel_veryloose    = ' < 0.228 '
    filt.cut_hovere_barrel_veryloose       = ' < 0.356 '
    filt.cut_isoRho_barrel_veryloose       = ' < 0.175 '
    filt.cut_ooEmooP_barrel_veryloose      = ' < 0.299 '
    #filt.cut_d0_barrel_veryloose           = ' < 0.0564 '
    #filt.cut_z0_barrel_veryloose           = ' < 0.472 '
    filt.cut_misshits_barrel_veryloose     = ' < 2 '
    filt.cut_passConvVeto_barrel_veryloose = ' == 1 '

    filt.cut_sigmaIEIE_endcap_tight        = ' < 0.0292 '
    filt.cut_absdEtaIn_endcap_tight        = ' < 0.00605 '
    filt.cut_absdPhiIn_endcap_tight        = ' < 0.0394 '
    filt.cut_hovere_endcap_tight           = ' < 0.0641 '
    filt.cut_isoRho_endcap_tight           = ' < 0.0571 '
    filt.cut_ooEmooP_endcap_tight          = ' < 0.0129 '
    #filt.cut_d0_endcap_tight               = ' < 0.0351 '
    #filt.cut_z0_endcap_tight               = ' < 0.417 '
    filt.cut_misshits_endcap_tight         = ' < 1 '
    filt.cut_passConvVeto_endcap_tight     = ' == 1 '

    filt.cut_sigmaIEIE_endcap_medium       = ' < 0.0298 '
    filt.cut_absdEtaIn_endcap_medium       = ' < 0.00609 '
    filt.cut_absdPhiIn_endcap_medium       = ' < 0.045 '
    filt.cut_hovere_endcap_medium          = ' < 0.0878 '
    filt.cut_isoRho_endcap_medium          = ' < 0.0821 '
    filt.cut_ooEmooP_endcap_medium         = ' < 0.13 '
    #filt.cut_d0_endcap_medium              = ' < 0.0739 '
    #filt.cut_z0_endcap_medium              = ' < 0.602 '
    filt.cut_misshits_endcap_medium        = ' <  1 '
    filt.cut_passConvVeto_endcap_medium    = ' == 1 '

    filt.cut_sigmaIEIE_endcap_loose        = ' < 0.0314 '
    filt.cut_absdEtaIn_endcap_loose        = ' < 0.00868 '
    filt.cut_absdPhiIn_endcap_loose        = ' < 0.213 '
    filt.cut_hovere_endcap_loose           = ' < 0.101 '
    filt.cut_isoRho_endcap_loose           = ' < 0.107 '
    filt.cut_ooEmooP_endcap_loose          = ' < 0.14 '
    #filt.cut_d0_endcap_loose               = ' < 0.118 '
    #filt.cut_z0_endcap_loose               = ' < 0.822 '
    filt.cut_misshits_endcap_loose         = ' <  1 '
    filt.cut_passConvVeto_endcap_loose     = ' == 1 '

    filt.cut_sigmaIEIE_endcap_veryloose    = ' < 0.037 '
    filt.cut_absdEtaIn_endcap_veryloose    = ' < 0.00895 '
    filt.cut_absdPhiIn_endcap_veryloose    = ' < 0.213 '
    filt.cut_hovere_endcap_veryloose       = ' < 0.211 '
    filt.cut_isoRho_endcap_veryloose       = ' < 0.159 '
    filt.cut_ooEmooP_endcap_veryloose      = ' < 0.15 '
    #filt.cut_d0_endcap_veryloose           = ' < 0.222 '
    #filt.cut_z0_endcap_veryloose           = ' < 0.921 '
    filt.cut_misshits_endcap_veryloose     = ' < 3 '
    filt.cut_passConvVeto_endcap_veryloose = ' == 1 '


    if do_cutflow :
        filt.do_cutflow = True
        filt.add_var( 'evalPID', evalPID )

    if apply_corrections :
        workarea = os.getenv('WorkArea')
        filt.add_var('applyCorrections', 'true' )
        filt.add_var('correctionFile', '%s/TreeFilter/RecoWgg/data/step2-invMass_SC-loose-Et_20-trigger-noPF-HggRunEtaR9.dat' %workarea )
        filt.add_var('smearingFile', '%s/TreeFilter/RecoWgg/data/outFile-step4-invMass_SC-loose-Et_20-trigger-noPF-HggRunEtaR9-smearEle.dat' %workarea )


    if do_hists :
        filt.add_hist( 'cut_pt', 100, 0, 500 )
        filt.add_hist( 'cut_abseta', 50, 0, 5 )
        filt.add_hist( 'cut_abseta_crack', 50, 0, 5 )
        filt.add_hist( 'cut_absdEtaIn_barrel_tight', 100, -0.1, 0.1 )
        filt.add_hist( 'cut_absdPhiIn_barrel_tight', 100, -0.1, 0.1 )
        filt.add_hist( 'cut_sigmaIEIE_barrel_tight', 100, 0, 0.05 )
        filt.add_hist( 'cut_hovere_barrel_tight', 100, -1, 1 )
        filt.add_hist( 'cut_d0_barrel_tight', 100, -1, 1 )
        filt.add_hist( 'cut_z0_barrel_tight', 100, -1, 1 )
        filt.add_hist( 'cut_ooEmooP_barrel_tight', 100, 0, 1 )
        filt.add_hist( 'cut_pfIso30_barrel_tight', 100, 0, 10 )
        filt.add_hist( 'cut_passConvVeto_barrel_tight', 2, 0, 2 )
        filt.add_hist( 'cut_misshits_barrel_tight', 10, 0, 10 )


    return filt

def filter_photon( ph_pt = ' > 15 ', id_cut='medium', ieta_cut=None, ele_veto='None', ele_olap='True', do_cutflow=False, do_hists=False, evalPID='medium' ) :

    filt = Filter('FilterPhoton')

    filt.cut_pt           = ph_pt
    filt.cut_abseta       = ' < 2.5'
    filt.cut_abseta_crack = ' > 1.44 & < 1.57 '
    filt.invert('cut_abseta_crack')

    filt.cut_muon_dr    = ' > 0.4 '
    if ele_olap == 'True' :
        filt.cut_electron_dr    = ' > 0.4 '

    if ieta_cut is not None :
        if ieta_cut == 'EB' :
            filt.cut_eb = ' == True '
        if ieta_cut == 'EE' :
            filt.cut_ee = ' == True '

    if ele_veto is not 'None' :
        if ele_veto == 'True' :
            filt.cut_CSEV = ' == True '
        elif ele_veto == 'False'  :
            filt.cut_CSEV = ' == False '
    
    if( id_cut is not 'None' ) :
        setattr( filt, 'cut_%s' %id_cut, ' == True ' )

    filt.cut_sigmaIEIE_barrel_loose  = ' < 0.01031 '
    filt.cut_chIsoCorr_barrel_loose  = ' < 1.295 '
    filt.cut_neuIsoCorr_barrel_loose = ' < 10.910 '
    filt.cut_phoIsoCorr_barrel_loose = ' < 3.630 '
    filt.cut_hovere_barrel_loose = ' < 0.0597 '

    filt.cut_sigmaIEIE_endcap_loose  = ' < 0.03013 '
    filt.cut_chIsoCorr_endcap_loose  = ' < 1.011 '
    filt.cut_neuIsoCorr_endcap_loose = ' < 5.931 '
    filt.cut_phoIsoCorr_endcap_loose = ' < 6.641 '
    filt.cut_hovere_endcap_loose = ' < 0.0481 '

    filt.cut_sigmaIEIE_barrel_medium  = ' < 0.01022 '
    filt.cut_chIsoCorr_barrel_medium  = ' < 0.441 '
    filt.cut_neuIsoCorr_barrel_medium = ' < 2.725 '
    filt.cut_phoIsoCorr_barrel_medium = ' < 2.571 '
    filt.cut_hovere_barrel_medium = ' < 0.0396 '

    filt.cut_sigmaIEIE_endcap_medium  = ' < 0.03001 '
    filt.cut_chIsoCorr_endcap_medium  = ' < 0.442 '
    filt.cut_neuIsoCorr_endcap_medium = ' < 1.715 '
    filt.cut_phoIsoCorr_endcap_medium = ' < 3.863 '
    filt.cut_hovere_endcap_medium = ' < 0.0219 '

    filt.cut_sigmaIEIE_barrel_tight  = ' < 0.00994 '
    filt.cut_chIsoCorr_barrel_tight  = ' < 0.202 '
    filt.cut_neuIsoCorr_barrel_tight = ' < 0.264 '
    filt.cut_phoIsoCorr_barrel_tight = ' < 2.362 '
    filt.cut_hovere_barrel_tight = ' < 0.0269 '

    filt.cut_sigmaIEIE_endcap_tight  = ' < 0.0300 '
    filt.cut_chIsoCorr_endcap_tight  = ' < 0.034 '
    filt.cut_neuIsoCorr_endcap_tight = ' < 0.586 '
    filt.cut_phoIsoCorr_endcap_tight = ' < 2.617 '
    filt.cut_hovere_endcap_tight = ' < 0.0213 '

    if do_cutflow : 
        filt.do_cutflow = True
        filt.add_var( 'evalPID', evalPID )

    if do_hists :
        filt.add_hist( 'cut_pt', 100, 0, 500 )
        filt.add_hist( 'cut_abseta', 50, 0, 5 )
        filt.add_hist( 'cut_abseta_crack', 50, 0, 5 )
        filt.add_hist( 'cut_hovere', 50, 0, 0.1 )
        filt.add_hist( 'cut_eveto', 2, 0, 2 )
        filt.add_hist( 'cut_sigmaIEIE_barrel_medium', 50, 0, 0.05 )
        filt.add_hist( 'cut_chIsoCorr_barrel_medium', 50, 0, 5 )
        filt.add_hist( 'cut_neuIsoCorr_barrel_medium', 50, 0, 5 )
        filt.add_hist( 'cut_phoIsoCorr_barrel_medium', 50, 0, 5 )
        filt.add_hist( 'cut_sigmaIEIE_endcap_medium', 50, 0, 0.05 )
        filt.add_hist( 'cut_chIsoCorr_endcap_medium', 50, 0, 5 )
        filt.add_hist( 'cut_neuIsoCorr_endcap_medium', 50, 0, 5 )
        filt.add_hist( 'cut_phoIsoCorr_endcap_medium', 50, 0, 5 )

        filt.add_hist( 'cut_sigmaIEIE_barrel_loose', 50, 0, 0.05 )
        filt.add_hist( 'cut_chIsoCorr_barrel_loose', 50, 0, 5 )
        filt.add_hist( 'cut_neuIsoCorr_barrel_loose', 100, -5, 5 )
        filt.add_hist( 'cut_phoIsoCorr_barrel_loose', 50, 0, 5 )
        filt.add_hist( 'cut_sigmaIEIE_endcap_loose', 50, 0, 0.05 )
        filt.add_hist( 'cut_chIsoCorr_endcap_loose', 50, 0, 5 )
        filt.add_hist( 'cut_neuIsoCorr_endcap_loose', 100, -5, 5 )
        filt.add_hist( 'cut_phoIsoCorr_endcap_loose', 50, 0, 5 )

    return filt

def filter_jet( jet_pt = ' > 30 ', do_hists=False ) :

    filt = Filter( 'FilterJet' )

    filt.cut_pt = jet_pt
    filt.cut_abseta = ' < 4.5 '

    filt.cut_muon_dr    = ' > 0.4 '
    filt.cut_electron_dr    = ' > 0.4 '
    filt.cut_photon_dr    = ' > 0.4 '

    filt.cut_jet_nhf_central_loose = ' < 0.99 '
    filt.cut_jet_nemf_central_loose = ' < 0.99 '
    filt.cut_jet_nconst_central_loose = ' > 1'
    filt.cut_jet_nhf_central_tight = ' < 0.90 '
    filt.cut_jet_nemf_central_tight = ' < 0.90 '
    filt.cut_jet_nconst_central_tight = ' > 1'
    filt.cut_jet_nhf_central_tightlep = ' < 0.9 '
    filt.cut_jet_nemf_central_tightlep = ' < 0.9 '
    filt.cut_jet_nconst_central_tightlep = '> 1 '
    filt.cut_jet_muf_central_tightlep = '< 0.8 '
    filt.cut_jet_chf_central_loose = ' > 0'
    filt.cut_jet_cmult_central_loose = ' > 0 '
    filt.cut_jet_cemf_central_loose = ' < 0.99 '
    filt.cut_jet_chf_central_tight = ' > 0 '
    filt.cut_jet_cmult_central_tight = ' > 0'
    filt.cut_jet_cemf_central_tight = ' < 0.99 '
    filt.cut_jet_chf_central_tightlep = ' > 0'
    filt.cut_jet_cmult_central_tightlep = ' > 0'
    filt.cut_jet_cemf_central_tightlep = ' < 0.90'
    filt.cut_jet_nhf_transition_loose = ' > 0.01 '
    filt.cut_jet_nemf_transition_loose = ' < 0.98 '
    filt.cut_jet_nmult_transition_loose = ' > 2 '
    filt.cut_jet_nhf_transition_tight = ' > 0.01 '
    filt.cut_jet_nemf_transition_tight = ' < 0.98 '
    filt.cut_jet_nmult_transition_tight = ' > 2 '
    filt.cut_jet_nemf_forward_loose = ' < 0.90 '
    filt.cut_jet_nmult_forward_loose = '> 10 '
    filt.cut_jet_nemf_forward_tight = ' < 0.90 '
    filt.cut_jet_nmult_forward_tight = ' > 10 '

    #filt.cut_jet_el_dr = ' > 0.4 '
    #filt.cut_jet_ph_dr = ' > 0.4 '

    filt.cut_jet_CSV_Loose = ' > 0.5426 '
    filt.cut_jet_CSV_Medium = ' > 0.8484 '
    filt.cut_jet_CSV_Tight = ' > 0.9535 '

    if do_hists :
        filt.add_hist( 'cut_pt', 100, 0, 500 )
        filt.add_hist( 'cut_abseta', 50, 0, 5 )

    return filt

def build_truth( args ) :

    truth_filt = Filter('BuildTruth') 

    truth_filt.cut_lep_mother = ' == 24 || == -24 ||  == 11 || == -11 || == 12 || == -12 || == 13 || == -13 || == 14 || == -14 || == 15 || == -15 || == 16 || == -16 '
    #truth_filt.cut_lep_status = ' != 23 '

    truth_filt.cut_ph_pt = ' > 5 '
    truth_filt.cut_ph_IsPromptFinalState = ' == True '

    doFHPFS = args.get( 'doFHPFS', False )
    if doFHPFS == 'true' :
        truth_filt.cut_ph_FromHardProcessFinalState = ' == True '

    return truth_filt 

def weight_event( args ) :

    filt = Filter( 'WeightEvent' )

    filt_str = args.get( 'ApplyNLOWeight', 'false' )

    filt.add_var( 'ApplyNLOWeight', filt_str )

    if 'sampleFile' in args :

        workarea = os.getenv('WorkArea')
        filt.add_var( 'sample_file', args['sampleFile'])
        filt.add_var( 'data_file', '%s/TreeFilter/RecoResonance/data/DataPileupHistogram2016.root' %workarea )
        filt.add_var( 'sample_hist', 'pileup_true' )
        filt.add_var('data_hist', 'pileup')
    else :
        print 'weight_event requires as a command line argument like --moduleArgs " { \'sampleFile\' : \'/path/histograms.root\'} "'

    return filt

