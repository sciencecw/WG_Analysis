#include "RunAnalysis.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <stdlib.h>

#include "BranchDefs.h"
#include "BranchInit.h"

#include "Util.h"

#include "TFile.h"

int main(int argc, char **argv)
{

    //TH1::AddDirectory(kFALSE);
    CmdOptions options = ParseOptions( argc, argv );

    // Parse the text file and form the configuration object
    AnaConfig ana_config = ParseConfig( options.config_file, options );
    std::cout << "Configured " << ana_config.size() 
	      << " analysis modules " << std::endl;

    RunModule runmod;
    ana_config.Run(runmod, options);

    std::cout << "^_^ Finished ^_^" << std::endl;


}

void RunModule::initialize( TChain * chain, TTree * outtree, TFile *outfile,
                            const CmdOptions & options,
			    std::vector<ModuleConfig> &configs ) {

    _input_tree = chain;
    // *************************
    // initialize trees
    // *************************
    InitINTree(chain);
    InitOUTTree( outtree );

    OUT::mu_pt20_n                              = 0;
    OUT::mu_pt30_n                              = 0;
    OUT::mu_passTight                           = 0;
    OUT::mu_passMedium                          = 0;
    OUT::mu_passLoose                           = 0;
    OUT::mu_hasTruthMatchMu                     = 0;
    OUT::mu_truthMatchMu_dr                     = 0;
    OUT::mu_truthMatchMu_pt                     = 0;
    OUT::mu_hasTrigMatch                        = 0;
    OUT::mu_trigMatch_dr                        = 0;

    OUT::el_pt30_n                              = 0;
    OUT::el_passVeryLoose                       = 0;
    OUT::el_passLoose                           = 0;
    OUT::el_passMedium                          = 0;
    OUT::el_passTight                           = 0;
    OUT::el_hasTrigMatch                        = 0;
    OUT::el_trigMatch_dr                      = 0;
    OUT::el_hasTruthMatchEl                  = 0;
    OUT::el_truthMatchEl_dr                     = 0;
    OUT::el_truthMatchEl_pt                     = 0;

    OUT::ph_chIsoCorr                           = 0;
    OUT::ph_neuIsoCorr                          = 0;
    OUT::ph_phoIsoCorr                          = 0;
    OUT::ph_min_el_dr                          = 0;
    OUT::ph_IsEB                                = 0;
    OUT::ph_IsEE                                = 0;
    OUT::ph_passTight                           = 0;
    OUT::ph_passMedium                          = 0;
    OUT::ph_passLoose                           = 0;
    OUT::ph_passLooseNoSIEIE                    = 0;
    OUT::ph_passHOverELoose                     = 0;
    OUT::ph_passHOverEMedium                    = 0;
    OUT::ph_passHOverETight                     = 0;
    OUT::ph_passSIEIELoose                      = 0;
    OUT::ph_passSIEIEMedium                     = 0;
    OUT::ph_passSIEIETight                      = 0;
    OUT::ph_passChIsoCorrLoose                  = 0;
    OUT::ph_passChIsoCorrMedium                 = 0;
    OUT::ph_passChIsoCorrTight                  = 0;
    OUT::ph_passNeuIsoCorrLoose                 = 0;
    OUT::ph_passNeuIsoCorrMedium                = 0;
    OUT::ph_passNeuIsoCorrTight                 = 0;
    OUT::ph_passPhoIsoCorrLoose                 = 0;
    OUT::ph_passPhoIsoCorrMedium                = 0;
    OUT::ph_passPhoIsoCorrTight                 = 0;
    OUT::ph_hasTruthMatchPh                 = 0;
    OUT::ph_truthMatchPh_dr                     = 0;
    OUT::ph_truthMatchPh_pt                     = 0;

    OUT::jet_CSVLoose_n                         = 0;
    OUT::jet_CSVMedium_n                        = 0;
    OUT::jet_CSVTight_n                         = 0;
    OUT::jet_IdLoose                            = 0;
    OUT::jet_IdTight                            = 0;
    OUT::jet_IdTightLep                         = 0;
    
    OUT::m_lep_ph                               = 0;
    OUT::m_lep_ph_comb_leadLep                  = 0;
    OUT::m_lep_ph_comb_sublLep                  = 0;
    OUT::m_lep_met_ph                           = 0;
    OUT::m_mt_lep_met_ph = 0;
    OUT::m_mt_lep_met_ph_forcewmass = 0;
    OUT::mt_w = 0;
    OUT::mt_res = 0;
    OUT::mt_lep_ph  = 0;
    OUT::dphi_lep_ph                            = 0;
    OUT::dr_lep_ph                              = 0;
    OUT::mt_lep_met                             = 0;
    OUT::m_lep_met                              = 0;
    OUT::pt_lep_met                             = 0;
    OUT::dphi_lep_met                           = 0;
    OUT::mt_lep_met_ph                          = 0;
    OUT::mt_lep_met_ph_inv                      = 0;
    OUT::dphi_ph_w                              = 0;
    OUT::pt_lep_met_ph                          = 0;
    OUT::RecoWMass                              = 0;
    OUT::recoM_lep_nu_ph                        = 0;
    OUT::recoMet_eta                            = 0;
    OUT::recoW_pt                               = 0;
    OUT::recoW_eta                              = 0;
    OUT::recoW_phi                              = 0;
    OUT::m_ll                                   = 0;
    OUT::nu_z_solution_success                  = 0;
   
    OUT::leadjet_pt                             = 0;
    OUT::subljet_pt                             = 0;
    OUT::leaddijet_m                            = 0;
    OUT::leaddijet_pt                           = 0;
    OUT::massdijet_m                            = 0;
    OUT::massdijet_pt                           = 0;

    OUT::isBlinded                              = 0;
    OUT::SampleWeight                           = 0;

    OUT::ph_loose_n                            = 0;
    OUT::ph_medium_n                            = 0;
    OUT::ph_tight_n                            = 0;

    OUT::ph_mediumPassPSV_n                     = 0;
    OUT::ph_mediumFailPSV_n                     = 0;
    OUT::ph_mediumPassCSEV_n                    = 0;
    OUT::ph_mediumFailCSEV_n                    = 0;

    OUT::ph_mediumPassEleOlap_n                 = 0;
    OUT::ph_mediumPassEleOlapPassCSEV_n         = 0;
    OUT::ph_mediumPassEleOlapFailCSEV_n         = 0;

    OUT::ph_mediumNoSIEIE_n                     = 0;
    OUT::ph_mediumNoChIso_n                     = 0;
    OUT::ph_mediumNoNeuIso_n                    = 0;
    OUT::ph_mediumNoPhoIso_n                    = 0;

    OUT::ph_mediumNoSIEIENoChIso_n              = 0;
    OUT::ph_mediumNoSIEIENoNeuIso_n             = 0;
    OUT::ph_mediumNoSIEIENoPhoIso_n             = 0;
    OUT::ph_mediumNoChIsoNoPhoIso_n             = 0;
    OUT::ph_mediumNoChIsoNoNeuIso_n             = 0;
    OUT::ph_mediumNoPhoIsoNoNeuIso_n            = 0;

    OUT::ph_mediumNoSIEIEPassPSV_n              = 0;
    OUT::ph_mediumNoChIsoPassPSV_n              = 0;
    OUT::ph_mediumNoNeuIsoPassPSV_n             = 0;
    OUT::ph_mediumNoPhoIsoPassPSV_n             = 0;

    OUT::ph_mediumNoSIEIEFailPSV_n              = 0;
    OUT::ph_mediumNoChIsoFailPSV_n              = 0;
    OUT::ph_mediumNoNeuIsoFailPSV_n             = 0;
    OUT::ph_mediumNoPhoIsoFailPSV_n             = 0;

    OUT::ph_mediumNoSIEIEPassCSEV_n             = 0;
    OUT::ph_mediumNoChIsoPassCSEV_n             = 0;
    OUT::ph_mediumNoNeuIsoPassCSEV_n            = 0;
    OUT::ph_mediumNoPhoIsoPassCSEV_n            = 0;

    OUT::ph_mediumNoSIEIEFailCSEV_n             = 0;
    OUT::ph_mediumNoChIsoFailCSEV_n             = 0;
    OUT::ph_mediumNoNeuIsoFailCSEV_n            = 0;
    OUT::ph_mediumNoPhoIsoFailCSEV_n            = 0;

    OUT::ptSorted_ph_loose_idx                  = 0;
    OUT::ptSorted_ph_medium_idx                 = 0;
    OUT::ptSorted_ph_tight_idx                  = 0;

    OUT::ptSorted_ph_mediumPassPSV_idx          = 0;
    OUT::ptSorted_ph_mediumFailPSV_idx          = 0;
    OUT::ptSorted_ph_mediumPassCSEV_idx         = 0;
    OUT::ptSorted_ph_mediumFailCSEV_idx         = 0;

    OUT::ptSorted_ph_mediumPassEleOlap_idx         = 0;
    OUT::ptSorted_ph_mediumPassEleOlapPassCSEV_idx = 0;
    OUT::ptSorted_ph_mediumPassEleOlapFailCSEV_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIE_idx          = 0;
    OUT::ptSorted_ph_mediumNoChIso_idx          = 0;
    OUT::ptSorted_ph_mediumNoNeuIso_idx         = 0;
    OUT::ptSorted_ph_mediumNoPhoIso_idx         = 0;

    OUT::ptSorted_ph_mediumNoSIEIENoChIso_idx   = 0;
    OUT::ptSorted_ph_mediumNoSIEIENoNeuIso_idx  = 0;
    OUT::ptSorted_ph_mediumNoSIEIENoPhoIso_idx  = 0;
    OUT::ptSorted_ph_mediumNoChIsoNoPhoIso_idx  = 0;
    OUT::ptSorted_ph_mediumNoChIsoNoNeuIso_idx  = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIEPassPSV_idx   = 0;
    OUT::ptSorted_ph_mediumNoChIsoPassPSV_idx   = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoPassPSV_idx  = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoPassPSV_idx  = 0;

    OUT::ptSorted_ph_mediumNoSIEIEFailPSV_idx   = 0;
    OUT::ptSorted_ph_mediumNoChIsoFailPSV_idx   = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoFailPSV_idx  = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoFailPSV_idx  = 0;

    OUT::ptSorted_ph_mediumNoSIEIEPassCSEV_idx  = 0;
    OUT::ptSorted_ph_mediumNoChIsoPassCSEV_idx  = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoPassCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoPassCSEV_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIEFailCSEV_idx  = 0;
    OUT::ptSorted_ph_mediumNoChIsoFailCSEV_idx  = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoFailCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoFailCSEV_idx = 0;

    OUT::trueph_n                               = 0;
    OUT::trueph_pt                              = 0;
    OUT::trueph_eta                             = 0;
    OUT::trueph_phi                             = 0;
    OUT::trueph_motherPID                       = 0;
    OUT::trueph_status                          = 0;

    OUT::truelep_n                              = 0;
    OUT::truelep_pt                             = 0;
    OUT::truelep_eta                            = 0;
    OUT::truelep_phi                            = 0;
    OUT::truelep_e                              = 0;
    OUT::truelep_motherPID                      = 0;
    OUT::truelep_status                         = 0;
    OUT::truelep_PID                             = 0;

    OUT::truenu_n                               = 0;

    OUT::truelepnu_m                            = 0;
    OUT::truelepnuph_m                          = 0;
    OUT::truelepph_dr                           = 0;
    OUT::truemt_lep_met_ph                      = 0;
    OUT::truemt_res                             = 0;
    OUT::truemt_res_l23                             = 0;
    OUT::truemt_res_lO                             = 0;
    OUT::trueht                                 = 0;

    OUT::isWMuDecay                             = 0;
    OUT::isWElDecay                             = 0;
    OUT::isWTauDecay                            = 0;
    OUT::isWTauElDecay                          = 0;
    OUT::isWTauMuDecay                          = 0;
    OUT::WIDStep                                = 0;

    OUT::NLOWeight                              = 1;

    OUT::PUWeight                               = 1;
    OUT::PUWeightUP5                            = 1;
    OUT::PUWeightUP10                           = 1;
    OUT::PUWeightDN5                            = 1;
    OUT::PUWeightDN10                           = 1;

    // *************************
    // Declare Branches
    // *************************

    bool build_truth = false;
    BOOST_FOREACH( ModuleConfig & mod_conf, configs ) {
        if( mod_conf.GetName() == "BuildTruth" )  build_truth = true;
    }
    outtree->Branch("mu_pt20_n", &OUT::mu_pt20_n, "mu_pt20_n/I"  );
    outtree->Branch("mu_pt30_n", &OUT::mu_pt30_n, "mu_pt30_n/I"  );
    outtree->Branch("mu_passTight", &OUT::mu_passTight            );
    outtree->Branch("mu_passMedium", &OUT::mu_passMedium           );
    outtree->Branch("mu_passLoose", &OUT::mu_passLoose            );
    outtree->Branch("mu_hasTruthMatchMu", &OUT::mu_hasTruthMatchMu );
    outtree->Branch("mu_truthMatchMu_dr", &OUT::mu_truthMatchMu_dr   );
    outtree->Branch("mu_truthMatchMu_pt", &OUT::mu_truthMatchMu_pt      );
    outtree->Branch("mu_hasTrigMatch", &OUT::mu_hasTrigMatch );
    outtree->Branch("mu_trigMatch_dr", &OUT::mu_trigMatch_dr);

    outtree->Branch("el_pt30_n", &OUT::el_pt30_n, "el_pt30_n/I"  );
    outtree->Branch("el_passVeryLoose", &OUT::el_passVeryLoose );
    outtree->Branch("el_passLoose", &OUT::el_passLoose );
    outtree->Branch("el_passMedium", &OUT::el_passMedium );
    outtree->Branch("el_passTight", &OUT::el_passTight );
    outtree->Branch("el_hasTrigMatch", &OUT::el_hasTrigMatch );
    outtree->Branch("el_trigMatch_dr", &OUT::el_trigMatch_dr );
    outtree->Branch("el_hasTruthMatchEl", &OUT::el_hasTruthMatchEl );
    outtree->Branch("el_truthMatchEl_dr", &OUT::el_truthMatchEl_dr   );
    outtree->Branch("el_truthMatchEl_pt", &OUT::el_truthMatchEl_pt      );

    outtree->Branch("ph_chIsoCorr", &OUT::ph_chIsoCorr);
    outtree->Branch("ph_neuIsoCorr", &OUT::ph_neuIsoCorr);
    outtree->Branch("ph_phoIsoCorr", &OUT::ph_phoIsoCorr);
    outtree->Branch("ph_min_el_dr", &OUT::ph_min_el_dr);
    outtree->Branch("ph_IsEB", &OUT::ph_IsEB );
    outtree->Branch("ph_IsEE", &OUT::ph_IsEE );
    outtree->Branch("ph_passTight", &OUT::ph_passTight );
    outtree->Branch("ph_passMedium", &OUT::ph_passMedium );
    outtree->Branch("ph_passLoose", &OUT::ph_passLoose );
    outtree->Branch("ph_passLooseNoSIEIE", &OUT::ph_passLooseNoSIEIE );
    outtree->Branch("ph_passHOverELoose", &OUT::ph_passHOverELoose );
    outtree->Branch("ph_passHOverEMedium", &OUT::ph_passHOverEMedium );
    outtree->Branch("ph_passHOverETight", &OUT::ph_passHOverETight );
    outtree->Branch("ph_passSIEIELoose", &OUT::ph_passSIEIELoose );
    outtree->Branch("ph_passSIEIEMedium", &OUT::ph_passSIEIEMedium );
    outtree->Branch("ph_passSIEIETight", &OUT::ph_passSIEIETight );
    outtree->Branch("ph_passChIsoCorrLoose", &OUT::ph_passChIsoCorrLoose );
    outtree->Branch("ph_passChIsoCorrMedium", &OUT::ph_passChIsoCorrMedium );
    outtree->Branch("ph_passChIsoCorrTight", &OUT::ph_passChIsoCorrTight );
    outtree->Branch("ph_passNeuIsoCorrLoose", &OUT::ph_passNeuIsoCorrLoose );
    outtree->Branch("ph_passNeuIsoCorrMedium", &OUT::ph_passNeuIsoCorrMedium );
    outtree->Branch("ph_passNeuIsoCorrTight", &OUT::ph_passNeuIsoCorrTight );
    outtree->Branch("ph_passPhoIsoCorrLoose", &OUT::ph_passPhoIsoCorrLoose );
    outtree->Branch("ph_passPhoIsoCorrMedium", &OUT::ph_passPhoIsoCorrMedium );
    outtree->Branch("ph_passPhoIsoCorrTight", &OUT::ph_passPhoIsoCorrTight );
    outtree->Branch("ph_hasTruthMatchPh", &OUT::ph_hasTruthMatchPh );
    outtree->Branch("ph_truthMatchPh_dr", &OUT::ph_truthMatchPh_dr   );
    outtree->Branch("ph_truthMatchPh_pt", &OUT::ph_truthMatchPh_pt      );


    outtree->Branch("jet_CSVLoose_n", &OUT::jet_CSVLoose_n, "jet_CSVLoose_n/I"  );
    outtree->Branch("jet_CSVMedium_n", &OUT::jet_CSVMedium_n, "jet_CSVMedium_n/I"  );
    outtree->Branch("jet_CSVTight_n", &OUT::jet_CSVTight_n, "jet_CSVTight_n/I"  );
    outtree->Branch("jet_IdLoose", &OUT::jet_IdLoose );
    outtree->Branch("jet_IdTight", &OUT::jet_IdTight );
    outtree->Branch("jet_IdTightLep", &OUT::jet_IdTightLep );

    outtree->Branch("m_lep_ph"        , &OUT::m_lep_ph        , "m_lep_ph/F"  );
    outtree->Branch("m_lep_ph_comb_sublLep"        , &OUT::m_lep_ph_comb_sublLep);
    outtree->Branch("m_lep_ph_comb_leadLep"        , &OUT::m_lep_ph_comb_leadLep);
    outtree->Branch("m_lep_met_ph"        , &OUT::m_lep_met_ph        , "m_lep_met_ph/F"  );
    outtree->Branch("m_mt_lep_met_ph"        , &OUT::m_mt_lep_met_ph        , "m_mtlep_met_ph/F"  );
    outtree->Branch("m_mt_lep_met_ph_forcewmass"        , &OUT::m_mt_lep_met_ph_forcewmass        , "m_mt_lep_met_ph_forcewmass/F"  );
    outtree->Branch("mt_w"        , &OUT::mt_w        , "mt_w/F"  );
    outtree->Branch("mt_res"        , &OUT::mt_res        , "mt_res/F"  );
    outtree->Branch("mt_lep_ph"        , &OUT::mt_lep_ph        , "mt_lep_ph/F"  );
    outtree->Branch("dphi_lep_ph"        , &OUT::dphi_lep_ph        , "dphi_lep_ph/F"  );
    outtree->Branch("dr_lep_ph"        , &OUT::dr_lep_ph        , "dr_lep_ph/F"  );
    outtree->Branch("mt_lep_met"      , &OUT::mt_lep_met      , "mt_lep_met/F" );
    outtree->Branch("m_lep_met"       , &OUT::m_lep_met       , "m_lep_met/F" );
    outtree->Branch("pt_lep_met"      , &OUT::pt_lep_met      , "pt_lep_met/F" );
    outtree->Branch("dphi_lep_met"    , &OUT::dphi_lep_met    , "dphi_lep_met/F" );
    outtree->Branch("mt_lep_met_ph"   , &OUT::mt_lep_met_ph   , "mt_lep_met_ph/F");
    outtree->Branch("mt_lep_met_ph_inv"   , &OUT::mt_lep_met_ph_inv   , "mt_lep_met_ph_inv/F");
    outtree->Branch("dphi_ph_w"   , &OUT::dphi_ph_w   , "dphi_ph_w/F");
    outtree->Branch("pt_lep_met_ph"   , &OUT::pt_lep_met_ph   , "pt_lep_met_ph/F");
    outtree->Branch("RecoWMass"       , &OUT::RecoWMass       , "RecoWMass/F");
    outtree->Branch("recoM_lep_nu_ph" , &OUT::recoM_lep_nu_ph , "recoM_lep_nu_ph/F");
    outtree->Branch("recoMet_eta" , &OUT::recoMet_eta, "recoMet_eta/F");
    outtree->Branch("recoW_pt" , &OUT::recoW_pt, "recoW_pt/F");
    outtree->Branch("recoW_eta" , &OUT::recoW_eta, "recoW_eta/F");
    outtree->Branch("recoW_phi" , &OUT::recoW_phi, "recoW_phi/F");
    outtree->Branch("m_ll" , &OUT::m_ll, "m_ll/F");
    outtree->Branch("nu_z_solution_success" , &OUT::nu_z_solution_success, "nu_z_solution_success/O");


    outtree->Branch("leadjet_pt", &OUT::leadjet_pt, "leadjet_pt/F" );
    outtree->Branch("subljet_pt", &OUT::subljet_pt, "subljet_pt/F" );
    outtree->Branch("leaddijet_m", &OUT::leaddijet_m, "leaddijet_m/F" );
    outtree->Branch("leaddijet_pt", &OUT::leaddijet_pt, "leaddijet_pt/F" );
    outtree->Branch("massdijet_m", &OUT::massdijet_m, "massdijet_m/F" );
    outtree->Branch("massdijet_pt", &OUT::massdijet_pt, "massdijet_pt/F" );

    outtree->Branch("isBlinded" , &OUT::isBlinded, "isBlinded/O");
    outtree->Branch("SampleWeight" , &OUT::SampleWeight, "SampleWeight/F");

#ifdef MODULE_MakePhotonCountVars

    outtree->Branch("ph_loose_n"           , &OUT::ph_loose_n           , "ph_loose_n/I"           );
    outtree->Branch("ph_medium_n"           , &OUT::ph_medium_n           , "ph_medium_n/I"           );
    outtree->Branch("ph_tight_n"           , &OUT::ph_tight_n           , "ph_tight_n/I"           );

    outtree->Branch("ph_mediumPassPSV_n"           , &OUT::ph_mediumPassPSV_n           , "ph_mediumPassPSV_n/I"           );
    outtree->Branch("ph_mediumFailPSV_n"           , &OUT::ph_mediumFailPSV_n           , "ph_mediumFailPSV_n/I"           );
    outtree->Branch("ph_mediumPassCSEV_n"          , &OUT::ph_mediumPassCSEV_n          , "ph_mediumPassCSEV_n/I"          );
    outtree->Branch("ph_mediumFailCSEV_n"          , &OUT::ph_mediumFailCSEV_n          , "ph_mediumFailCSEV_n/I"          );

    outtree->Branch("ph_mediumPassEleOlap_n"          , &OUT::ph_mediumPassEleOlap_n, "ph_mediumPassEleOlap_n/I"          );
    outtree->Branch("ph_mediumPassEleOlapPassCSEV_n"          , &OUT::ph_mediumPassEleOlapPassCSEV_n, "ph_mediumPassEleOlapPassCSEV_n/I"          );
    outtree->Branch("ph_mediumPassEleOlapFailCSEV_n"          , &OUT::ph_mediumPassEleOlapFailCSEV_n, "ph_mediumPassEleOlapFailCSEV_n/I"          );

    outtree->Branch("ph_mediumNoSIEIE_n"           , &OUT::ph_mediumNoSIEIE_n           , "ph_mediumNoSIEIE_n/I"           );
    outtree->Branch("ph_mediumNoChIso_n"           , &OUT::ph_mediumNoChIso_n           , "ph_mediumNoChIso_n/I"           );
    outtree->Branch("ph_mediumNoNeuIso_n"          , &OUT::ph_mediumNoNeuIso_n          , "ph_mediumNoNeuIso_n/I"          );
    outtree->Branch("ph_mediumNoPhoIso_n"          , &OUT::ph_mediumNoPhoIso_n          , "ph_mediumNoPhoIso_n/I"          );

    outtree->Branch("ph_mediumNoSIEIENoChIso_n", &OUT::ph_mediumNoSIEIENoChIso_n, "ph_mediumNoSIEIENoChIso_n/I" );
    outtree->Branch("ph_mediumNoSIEIENoNeuIso_n", &OUT::ph_mediumNoSIEIENoNeuIso_n, "ph_mediumNoSIEIENoNeuIso_n/I" );
    outtree->Branch("ph_mediumNoSIEIENoPhoIso_n", &OUT::ph_mediumNoSIEIENoPhoIso_n, "ph_mediumNoSIEIENoPhoIso_n/I" );
    outtree->Branch("ph_mediumNoChIsoNoPhoIso_n", &OUT::ph_mediumNoChIsoNoPhoIso_n, "ph_mediumNoChIsoNoPhoIso_n/I" );
    outtree->Branch("ph_mediumNoChIsoNoNeuIso_n", &OUT::ph_mediumNoChIsoNoNeuIso_n, "ph_mediumNoChIsoNoNeuIso_n/I" );
    outtree->Branch("ph_mediumNoPhoIsoNoNeuIso_n", &OUT::ph_mediumNoPhoIsoNoNeuIso_n, "ph_mediumNoPhoIsoNoNeuIso_n/I" );


    outtree->Branch("ph_mediumNoSIEIEPassPSV_n"           , &OUT::ph_mediumNoSIEIEPassPSV_n           , "ph_mediumNoSIEIEPassPSV_n/I"           );
    outtree->Branch("ph_mediumNoChIsoPassPSV_n"           , &OUT::ph_mediumNoChIsoPassPSV_n           , "ph_mediumNoChIsoPassPSV_n/I"           );
    outtree->Branch("ph_mediumNoNeuIsoPassPSV_n"          , &OUT::ph_mediumNoNeuIsoPassPSV_n          , "ph_mediumNoNeuIsoPassPSV_n/I"          );
    outtree->Branch("ph_mediumNoPhoIsoPassPSV_n"          , &OUT::ph_mediumNoPhoIsoPassPSV_n          , "ph_mediumNoPhoIsoPassPSV_n/I"          );

    outtree->Branch("ph_mediumNoSIEIEFailPSV_n"           , &OUT::ph_mediumNoSIEIEFailPSV_n           , "ph_mediumNoSIEIEFailPSV_n/I"           );
    outtree->Branch("ph_mediumNoChIsoFailPSV_n"           , &OUT::ph_mediumNoChIsoFailPSV_n           , "ph_mediumNoChIsoFailPSV_n/I"           );
    outtree->Branch("ph_mediumNoNeuIsoFailPSV_n"          , &OUT::ph_mediumNoNeuIsoFailPSV_n          , "ph_mediumNoNeuIsoFailPSV_n/I"          );
    outtree->Branch("ph_mediumNoPhoIsoFailPSV_n"          , &OUT::ph_mediumNoPhoIsoFailPSV_n          , "ph_mediumNoPhoIsoFailPSV_n/I"          );

    outtree->Branch("ph_mediumNoSIEIEPassCSEV_n"           , &OUT::ph_mediumNoSIEIEPassCSEV_n           , "ph_mediumNoSIEIEPassCSEV_n/I"           );
    outtree->Branch("ph_mediumNoChIsoPassCSEV_n"           , &OUT::ph_mediumNoChIsoPassCSEV_n           , "ph_mediumNoChIsoPassCSEV_n/I"           );
    outtree->Branch("ph_mediumNoNeuIsoPassCSEV_n"          , &OUT::ph_mediumNoNeuIsoPassCSEV_n          , "ph_mediumNoNeuIsoPassCSEV_n/I"          );
    outtree->Branch("ph_mediumNoPhoIsoPassCSEV_n"          , &OUT::ph_mediumNoPhoIsoPassCSEV_n          , "ph_mediumNoPhoIsoPassCSEV_n/I"          );

    outtree->Branch("ph_mediumNoSIEIEFailCSEV_n"           , &OUT::ph_mediumNoSIEIEFailCSEV_n           , "ph_mediumNoSIEIEFailCSEV_n/I"           );
    outtree->Branch("ph_mediumNoChIsoFailCSEV_n"           , &OUT::ph_mediumNoChIsoFailCSEV_n           , "ph_mediumNoChIsoFailCSEV_n/I"           );
    outtree->Branch("ph_mediumNoNeuIsoFailCSEV_n"          , &OUT::ph_mediumNoNeuIsoFailCSEV_n          , "ph_mediumNoNeuIsoFailCSEV_n/I"          );
    outtree->Branch("ph_mediumNoPhoIsoFailCSEV_n"          , &OUT::ph_mediumNoPhoIsoFailCSEV_n          , "ph_mediumNoPhoIsoFailCSEV_n/I"          );

    outtree->Branch("ptSorted_ph_loose_idx"  , &OUT::ptSorted_ph_loose_idx );
    outtree->Branch("ptSorted_ph_medium_idx"  , &OUT::ptSorted_ph_medium_idx );
    outtree->Branch("ptSorted_ph_tight_idx"  , &OUT::ptSorted_ph_tight_idx );

    outtree->Branch("ptSorted_ph_mediumPassPSV_idx"   , &OUT::ptSorted_ph_mediumPassPSV_idx );
    outtree->Branch("ptSorted_ph_mediumFailPSV_idx"   , &OUT::ptSorted_ph_mediumFailPSV_idx );
    outtree->Branch("ptSorted_ph_mediumPassCSEV_idx"  , &OUT::ptSorted_ph_mediumPassCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumFailCSEV_idx"  , &OUT::ptSorted_ph_mediumFailCSEV_idx );

    outtree->Branch("ptSorted_ph_mediumPassEleOlap_idx"  , &OUT::ptSorted_ph_mediumPassEleOlap_idx);
    outtree->Branch("ptSorted_ph_mediumPassEleOlapPassCSEV_idx"  , &OUT::ptSorted_ph_mediumPassEleOlapPassCSEV_idx);
    outtree->Branch("ptSorted_ph_mediumPassEleOlapFailCSEV_idx"  , &OUT::ptSorted_ph_mediumPassEleOlapFailCSEV_idx);

    outtree->Branch("ptSorted_ph_mediumNoSIEIE_idx"   , &OUT::ptSorted_ph_mediumNoSIEIE_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIso_idx"   , &OUT::ptSorted_ph_mediumNoChIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoNeuIso_idx"  , &OUT::ptSorted_ph_mediumNoNeuIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoPhoIso_idx"  , &OUT::ptSorted_ph_mediumNoPhoIso_idx );

    outtree->Branch("ptSorted_ph_mediumNoSIEIENoChIso_idx", &OUT::ptSorted_ph_mediumNoSIEIENoChIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoSIEIENoNeuIso_idx", &OUT::ptSorted_ph_mediumNoSIEIENoNeuIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoSIEIENoPhoIso_idx", &OUT::ptSorted_ph_mediumNoSIEIENoPhoIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIsoNoPhoIso_idx", &OUT::ptSorted_ph_mediumNoChIsoNoPhoIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIsoNoNeuIso_idx", &OUT::ptSorted_ph_mediumNoChIsoNoNeuIso_idx );
    outtree->Branch("ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx", &OUT::ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx );

    outtree->Branch("ptSorted_ph_mediumNoSIEIEPassPSV_idx"   , &OUT::ptSorted_ph_mediumNoSIEIEPassPSV_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIsoPassPSV_idx"   , &OUT::ptSorted_ph_mediumNoChIsoPassPSV_idx );
    outtree->Branch("ptSorted_ph_mediumNoNeuIsoPassPSV_idx"  , &OUT::ptSorted_ph_mediumNoNeuIsoPassPSV_idx );
    outtree->Branch("ptSorted_ph_mediumNoPhoIsoPassPSV_idx"  , &OUT::ptSorted_ph_mediumNoPhoIsoPassPSV_idx );

    outtree->Branch("ptSorted_ph_mediumNoSIEIEFailPSV_idx"   , &OUT::ptSorted_ph_mediumNoSIEIEFailPSV_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIsoFailPSV_idx"   , &OUT::ptSorted_ph_mediumNoChIsoFailPSV_idx );
    outtree->Branch("ptSorted_ph_mediumNoNeuIsoFailPSV_idx"  , &OUT::ptSorted_ph_mediumNoNeuIsoFailPSV_idx );
    outtree->Branch("ptSorted_ph_mediumNoPhoIsoFailPSV_idx"  , &OUT::ptSorted_ph_mediumNoPhoIsoFailPSV_idx );

    outtree->Branch("ptSorted_ph_mediumNoSIEIEPassCSEV_idx"   , &OUT::ptSorted_ph_mediumNoSIEIEPassCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIsoPassCSEV_idx"   , &OUT::ptSorted_ph_mediumNoChIsoPassCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumNoNeuIsoPassCSEV_idx"  , &OUT::ptSorted_ph_mediumNoNeuIsoPassCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumNoPhoIsoPassCSEV_idx"  , &OUT::ptSorted_ph_mediumNoPhoIsoPassCSEV_idx );

    outtree->Branch("ptSorted_ph_mediumNoSIEIEFailCSEV_idx"   , &OUT::ptSorted_ph_mediumNoSIEIEFailCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumNoChIsoFailCSEV_idx"   , &OUT::ptSorted_ph_mediumNoChIsoFailCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumNoNeuIsoFailCSEV_idx"  , &OUT::ptSorted_ph_mediumNoNeuIsoFailCSEV_idx );
    outtree->Branch("ptSorted_ph_mediumNoPhoIsoFailCSEV_idx"  , &OUT::ptSorted_ph_mediumNoPhoIsoFailCSEV_idx );

#endif

    if( build_truth ) {
        outtree->Branch("trueph_n"           , &OUT::trueph_n, "trueph_n/I" );
        outtree->Branch("trueph_pt"           , &OUT::trueph_pt                        );
        outtree->Branch("trueph_eta"          , &OUT::trueph_eta                       );
        outtree->Branch("trueph_phi"          , &OUT::trueph_phi                       );
        outtree->Branch("trueph_motherPID"    , &OUT::trueph_motherPID                 );
        outtree->Branch("trueph_status"       , &OUT::trueph_status                    );

        outtree->Branch("truelep_n"           , &OUT::truelep_n, "truelep_n/I" );
        outtree->Branch("truelep_pt"          , &OUT::truelep_pt                       );
        outtree->Branch("truelep_eta"         , &OUT::truelep_eta                      );
        outtree->Branch("truelep_phi"         , &OUT::truelep_phi                      );
        outtree->Branch("truelep_e"           , &OUT::truelep_e                        );
        outtree->Branch("truelep_motherPID"   , &OUT::truelep_motherPID                );
        outtree->Branch("truelep_status"      , &OUT::truelep_status                   );
        outtree->Branch("truelep_PID"         , &OUT::truelep_PID                      );

        outtree->Branch("truenu_n"            , &OUT::truenu_n, "truenu_n/I"           );

        outtree->Branch("truelepnu_m"         , &OUT::truelepnu_m, "truelepnu_m/F"     );
        outtree->Branch("truelepnuph_m"       , &OUT::truelepnuph_m, "truelepnuph_m/F" );
        outtree->Branch("truelepph_dr"        , &OUT::truelepph_dr, "truelepph_dr/F"   );
        outtree->Branch("truemt_lep_met_ph" , &OUT::truemt_lep_met_ph, "truemt_lep_met_ph/F");
        outtree->Branch("truemt_res" , &OUT::truemt_res, "truemt_res/F");
        outtree->Branch("truemt_res_l23" , &OUT::truemt_res_l23, "truemt_res_l23/F");
        outtree->Branch("truemt_res_lO" , &OUT::truemt_res_lO, "truemt_res_lO/F");
        outtree->Branch("trueht" , &OUT::trueht, "trueht/F");

        outtree->Branch("isWMuDecay"          , &OUT::isWMuDecay, "isWMuDecay/O"       );
        outtree->Branch("isWElDecay"          , &OUT::isWElDecay, "isWElDecay/O"       );
        outtree->Branch("isWTauDecay"         , &OUT::isWTauDecay, "isWTauDecay/O"     );
        outtree->Branch("isWTauElDecay"         , &OUT::isWTauElDecay, "isWTauElDecay/O"     );
        outtree->Branch("isWTauMuDecay"         , &OUT::isWTauMuDecay, "isWTauMuDecay/O"     );
        outtree->Branch("WIDStep"             , &OUT::WIDStep    , "WIDStep/I"     );
    }

    outtree->Branch("NLOWeight", &OUT::NLOWeight, "NLOWeight/F" );
    outtree->Branch("PUWeight", &OUT::PUWeight, "PUWeight/F" );
    outtree->Branch("PUWeightUP5", &OUT::PUWeightUP5, "PUWeightUP5/F" );
    outtree->Branch("PUWeightUP10", &OUT::PUWeightUP10, "PUWeightUP10/F" );
    outtree->Branch("PUWeightDN5", &OUT::PUWeightDN5, "PUWeightDN5/F" );
    outtree->Branch("PUWeightDN10", &OUT::PUWeightDN10, "PUWeightDN10/F" );

    BOOST_FOREACH( ModuleConfig & mod_conf, configs ) {
    
        if( mod_conf.GetName() == "FilterBlind" ) { 
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "isData" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::string data = eitr->second;
                std::transform(data.begin(), data.end(), data.begin(), ::tolower);
                if( data=="true") _isData=true;
                else              _isData=false;
            }
        }
        if( mod_conf.GetName() == "FilterMuon" ) { 
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "triggerMatchBits" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::vector<std::string> trigger_bit_list = Tokenize( eitr->second, "," );
                for( std::vector<std::string>::const_iterator bitr = trigger_bit_list.begin(); bitr != trigger_bit_list.end(); ++bitr ) {
                    std::stringstream ss_id( *bitr );
                    int trig_id;
                    ss_id >> trig_id;
                    _muonTrigMatchBits.push_back(trig_id);
                }
            }
            eitr = mod_conf.GetInitData().find( "evalPID" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::string pid = eitr->second;
                if( pid == "tight"     ) _eval_mu_tight       = true;
                if( pid == "medium"    ) _eval_mu_medium      = true;
                if( pid == "loose"    ) _eval_mu_loose      = true;
            }
        }
        if( mod_conf.GetName() == "FilterElectron" ) { 
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "triggerMatchBits" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::vector<std::string> trigger_bit_list = Tokenize( eitr->second, "," );
                for( std::vector<std::string>::const_iterator bitr = trigger_bit_list.begin(); bitr != trigger_bit_list.end(); ++bitr ) {
                    std::stringstream ss_id( *bitr );
                    int trig_id;
                    ss_id >> trig_id;
                    _electronTrigMatchBits.push_back(trig_id);
                }
            }
            eitr = mod_conf.GetInitData().find( "evalPID" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::string pid = eitr->second;
                if( pid == "tight"     ) _eval_el_tight       = true;
                if( pid == "medium"    ) _eval_el_medium      = true;
                if( pid == "loose"     ) _eval_el_loose       = true;
                if( pid == "veryloose" ) _eval_el_veryloose   = true;
            }
        }
        if( mod_conf.GetName() == "FilterPhoton" ) { 
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "evalPID" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::string pid = eitr->second;
                if( pid == "tight"     ) _eval_ph_tight       = true;
                if( pid == "medium"    ) _eval_ph_medium      = true;
                if( pid == "loose"     ) _eval_ph_loose       = true;
            }
        }
        if( mod_conf.GetName() == "FilterTrigger" ) { 
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "triggerBits" );
            if( eitr != mod_conf.GetInitData().end() ) {
                // get the mapping of trigger ID to the name
                // This could be taken from the TrigInfoTree
                std::vector<std::string> trigger_bit_list = Tokenize( eitr->second, "," );
                for( std::vector<std::string>::const_iterator bitr = trigger_bit_list.begin(); bitr != trigger_bit_list.end(); ++bitr ) {
                    std::vector<std::string> name_id_map = Tokenize( *bitr, ":" );
                    std::stringstream ss_id( name_id_map[0] );
                    // convert the ID to an int
                    int trig_id;
                    ss_id >> trig_id;
                    // make an entry in the output map
                    triggerResults[trig_id] = 0;
                    outtree->Branch(name_id_map[1].c_str(), &(triggerResults[trig_id]), (name_id_map[1]+"/O").c_str() );
                }
            }
            eitr = mod_conf.GetInitData().find( "AuxTreeName" );
            if( eitr != mod_conf.GetInitData().end() ) {
                // if the trigger bits were provided, do not overwrite
                if( triggerResults.size() ) {
                    std::cout << "WARNING -- triggerBits were provided.  Will not get trigger info from the tree.  Please provide either triggerBits or AuxTreeName" << std::endl;
                }
                else {
                    const std::string & trig_tree_name = eitr->second;
                    TFile * fptr = _input_tree->GetFile();
                    TTree * trigTree = dynamic_cast<TTree*>(fptr->Get(trig_tree_name.c_str()));
                    if( !trigTree ) {
                        std::cout << "Did not locate tree with name " + trig_tree_name + ". Will not store trigger results!" << std::endl;
                    }
                    else {
                        Int_t   trigger_ids;
                        Char_t  trigger_names;
                        std::cout << trigTree->GetName() << std::endl;
                        trigTree->SetBranchAddress("trigger_ids", &trigger_ids);
                        trigTree->SetBranchAddress("trigger_names", &trigger_names);

                        for( int tidx = 0; tidx < trigTree->GetEntries(); tidx++  ) {
                            trigTree->GetEntry(tidx);
                            triggerResults[trigger_ids] = 0;
                            std::string tn( &trigger_names);
                            outtree->Branch(tn.c_str(), &(triggerResults[trigger_ids]), (tn+"/O").c_str() );
                        }
                    }
                }
            }
        }
       if( mod_conf.GetName() == "FilterMET" ) {
            // working on met filter. Basically copied from trigger code
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "METFilterBits" );
            if( eitr != mod_conf.GetInitData().end() ) {
                // get the mapping of metfilter index to the name
                // This could be taken from the FilterInfoTree
                std::vector<std::string> metfilter_bit_list = Tokenize( eitr->second, "," );
                for( std::vector<std::string>::const_iterator bitr = metfilter_bit_list.begin(); bitr != metfilter_bit_list.end(); ++bitr ) {
                    std::vector<std::string> name_id_map = Tokenize( *bitr, ":" );
                    std::stringstream ss_id( name_id_map[0] );
                    // convert the ID to an int
                    int metfilter_id;
                    ss_id >> metfilter_id;
                    // make an entry in the output map
                    metfilterResults[metfilter_id] = 0;
                    outtree->Branch(name_id_map[1].c_str(), &(metfilterResults[metfilter_id]), (name_id_map[1]+"/O").c_str() );
                }
            }
            eitr = mod_conf.GetInitData().find( "METAuxTreeName" );
            if( eitr != mod_conf.GetInitData().end() ) {
                // if the metfilter bits were provided, do not overwrite
                if( metfilterResults.size() ) {
                    std::cout << "WARNING -- METFilterBits were provided.  Will not get metfilter info from the tree.  Please provide either metfilterBits or METAuxTreeName" << std::endl;
                }
                else {
                    const std::string & metfilter_tree_name = eitr->second;
                    TFile * fptr = _input_tree->GetFile();
                    TTree * metfilterTree = dynamic_cast<TTree*>(fptr->Get(metfilter_tree_name.c_str()));
                    if( !metfilterTree ) {
                        std::cout << "Did not locate tree with name " + metfilter_tree_name + ". Will not store metfilter results!" << std::endl;
                    }
                    else {
                        Int_t   filter_ids;
                        Char_t  filter_names;
                        std::cout << metfilterTree->GetName() << std::endl;
                        metfilterTree->SetBranchAddress("filter_ids",   &filter_ids);
                        metfilterTree->SetBranchAddress("filter_names", &filter_names);

                        for( int tidx = 0; tidx < metfilterTree->GetEntries(); tidx++  ) {
                            metfilterTree->GetEntry(tidx);
                            metfilterResults[filter_ids] = 0;
                            std::string tn( &filter_names);
                            outtree->Branch(tn.c_str(), &(metfilterResults[filter_ids]), (tn+"/O").c_str() );
                        }
                    }
                }
            }
        }
        if( mod_conf.GetName() == "WeightEvent" ) { 
            std::map<std::string, std::string>::const_iterator itr = mod_conf.GetInitData().find( "ApplyNLOWeight" );
            if( itr != mod_conf.GetInitData().end() ) {
                if( itr->second == "true" ) {
                    _needs_nlo_weght = true;
                }
            }
            std::string sample_filename;
            std::string data_filename;
            std::string sample_histname;
            std::string data_histname;
            std::map<std::string, std::string>::const_iterator sfitr = mod_conf.GetInitData().find( "sample_file" );
            std::map<std::string, std::string>::const_iterator dfitr = mod_conf.GetInitData().find( "data_file" );
            std::map<std::string, std::string>::const_iterator shitr = mod_conf.GetInitData().find( "sample_hist" );
            std::map<std::string, std::string>::const_iterator dhitr = mod_conf.GetInitData().find( "data_hist" );

            bool get_weight_hists = true;

            if( sfitr != mod_conf.GetInitData().end() ) {
                sample_filename = sfitr->second;
            }
            else {
                get_weight_hists = false;
            }
            if( dfitr != mod_conf.GetInitData().end() ) {
                data_filename = dfitr->second;
            }
            else {
                get_weight_hists = false;
            }
            if( shitr != mod_conf.GetInitData().end() ) {
                sample_histname = shitr->second;
            }
            else {
                get_weight_hists = false;
            }
            if( dhitr != mod_conf.GetInitData().end() ) {
                data_histname = dhitr->second;
            }
            else {
                get_weight_hists = false;
            }

            if( get_weight_hists ) {

                _puweight_sample_file = TFile::Open( sample_filename.c_str(), "READ" );
                _puweight_data_file = TFile::Open( data_filename.c_str(), "READ" );

                _puweight_sample_hist = dynamic_cast<TH1F*>(_puweight_sample_file->Get( sample_histname.c_str() ) ) ;
                _puweight_data_hist = dynamic_cast<TH1D*>(_puweight_data_file->Get( data_histname.c_str() ) );
                if( !_puweight_sample_hist ) {
                    std::cout << "Could not retrieve histogram " << sample_histname << " from " << sample_filename  << std::endl;
                }
                if( !_puweight_data_hist ) {
                    std::cout << "Could not retrieve histogram " << data_histname << " from " << data_filename << std::endl;
                }
            }
            else {
                std::cout << "WeightEvent::ERROR - Needed histogram does not exist! Will not apply PU weights!" << std::endl;
            }
            
        }
        if( mod_conf.GetName() == "FilterDataQuality" ) { 
            std::map<std::string, std::string>::const_iterator itr = mod_conf.GetInitData().find( "jsonFile" );
            if( itr != mod_conf.GetInitData().end() ) {
                std::string jsonFile = itr->second;

                std::string line;
                std::ifstream infile( jsonFile.c_str() );
                if( infile.is_open() ) {
                    while( getline( infile, line ) ) {
                        std::vector<std::string> key_val_tok = Tokenize( line, ":" );
                        std::vector<std::string> run_number_tok = Tokenize( key_val_tok[0], "\"" );
                        // find first and last brackets
                        std::string::size_type first_bracket = key_val_tok[1].find( "[" );
                        std::string::size_type last_bracket  = key_val_tok[1].rfind( "]" );

                        std::string rm_last = key_val_tok[1].substr( 0, last_bracket );
                        std::string all_vals  = rm_last.substr( first_bracket+1 );

                        int run_number;
                        std::stringstream run_number_ss( run_number_tok[1] );
                        run_number_ss >> run_number;

                        std::cout << "Run = " << run_number_tok[1] << " values = " << all_vals << std::endl;

                        std::vector<std::string> ranges_tok = Tokenize( all_vals, "[" );

                        std::vector<int> full_range;
                        for( std::vector<std::string>::const_iterator itr = ranges_tok.begin(); itr != ranges_tok.end(); ++itr ) {
                            std::string range_str = *itr;
                            std::string::size_type bracket_pos = range_str.find("]");
                            std::string range = range_str.substr( 0, bracket_pos );

                            std::vector<std::string> range_tok = Tokenize( range, "," );
                            if( range_tok.size() != 2 ) {
                                std::cout << "Expected two entries in the range.  String was " << range << std::endl;
                                continue;
                            }

                            int range_begin;
                            int range_end;

                            std::stringstream range_begin_ss( range_tok[0] );
                            std::stringstream range_end_ss  ( range_tok[1] );

                            range_begin_ss >> range_begin;
                            range_end_ss   >> range_end;

                            for( int ls = range_begin ; ls <= range_end; ++ls ) {
                                full_range.push_back(ls);
                            }

                        }
                        _quality_map[run_number] = full_range;
                    }
                }
            }
        }
    }

}

bool RunModule::execute( std::vector<ModuleConfig> & configs ) {

    // In BranchInit
    CopyInputVarsToOutput();

    // loop over configured modules
    bool save_event = true;
    BOOST_FOREACH( ModuleConfig & mod_conf, configs ) {
        save_event &= ApplyModule( mod_conf );
    }

    return save_event;

}

bool RunModule::ApplyModule( ModuleConfig & config ) {

    bool keep_evt = true;

    if( config.GetName() == "FilterMuon" ) {
        FilterMuon( config );
    }
    if( config.GetName() == "FilterElectron" ) {
        FilterElectron( config );
    }
    if( config.GetName() == "FilterPhoton" ) {
        FilterPhoton( config );
    }
    if( config.GetName() == "FilterJet" ) {
        FilterJet( config );
    }
    if( config.GetName() == "BuildEventVars" ) {
        BuildEventVars( config );
    }
    if( config.GetName() == "MakePhotonCountVars" ) {
        MakePhotonCountVars( config );
    }
    if( config.GetName() == "BuildTruth" ) {
        BuildTruth( config );
    }
    if( config.GetName() == "FilterEvent" ) {
        keep_evt &= FilterEvent( config );
    }
    if( config.GetName() == "FilterTrigger" ) {
        keep_evt &= FilterTrigger( config );
    }
    if( config.GetName() == "FilterMET" ){
//        keep_evt &= FilterMET( config );
         FilterMET( config );
    }
    if( config.GetName() == "WeightEvent" ) {
        WeightEvent( config );
    }
    if( config.GetName() == "FilterDataQuality" ) {
        keep_evt &= FilterDataQuality( config );
    }
    if( config.GetName() == "FilterBlind" ) {
        keep_evt &= FilterBlind( config );
    }


    return keep_evt;


}

void RunModule::FilterMuon( ModuleConfig & config ) {

    OUT::mu_n                  = 0;
    OUT::mu_pt20_n             = 0;
    OUT::mu_pt30_n             = 0;
    OUT::mu_passTight          -> clear();
    OUT::mu_passMedium         -> clear();
    OUT::mu_passLoose          -> clear();
    OUT::mu_hasTrigMatch       -> clear();
    OUT::mu_trigMatch_dr       -> clear();
    // truth matching is happeing in BuildTruth
    //OUT::mu_hasTruthMatchMu       -> clear();
    //OUT::mu_truthMatchMu_dr    -> clear();
    //OUT::mu_truthMatchMu_pt    -> clear();

    ClearOutputPrefix("mu_");

    for( int idx = 0; idx < IN::mu_n; ++idx ) {

        float pt = IN::mu_pt->at(idx);
        float eta = IN::mu_eta->at(idx);

        if( !config.PassFloat( "cut_pt", pt   ) ) continue;
        if( !config.PassFloat( "cut_eta", fabs(eta) ) ) continue;

        bool isPfMu = IN::mu_isPf->at(idx);
        bool isGloMu = IN::mu_isGlobal->at(idx);
        bool isTkMu = IN::mu_isTracker->at(idx);
        float chi2 = IN::mu_chi2->at(idx);
        int nHits =  IN::mu_nHits->at(idx);
        int nStations = IN::mu_nMuStations->at(idx);
        int nPixHits = IN::mu_nPixHits->at(idx);
        int nTrkLayers = IN::mu_nTrkLayers->at(idx);
        float d0 = IN::mu_d0->at(idx);
        float z0 = IN::mu_dz->at(idx);
        float pfIso = IN::mu_pfIso->at(idx);
        float tkIso = IN::mu_trkIso->at(idx);
        
        bool pass_tight = true;
        bool pass_loose = true;
        //bool pass_medium = true;

        bool use_eval = _eval_mu_tight || _eval_mu_medium || _eval_mu_loose;

        // loose cuts
        if( !use_eval || _eval_mu_loose ) {
            if( !config.PassBool ( "cut_isPf_loose"   , isPfMu ) ) {
                pass_loose = false;
                if( _eval_mu_loose) continue;
            }
            if( !config.PassBool ( "cut_isGlobalOrTk_loose"   , (isGloMu || isTkMu ) ) ) {
                pass_loose = false;
                if( _eval_mu_loose) continue;
            }
        }

        // loose cuts
        if( !use_eval || _eval_mu_tight ) {
            if( !config.PassBool ( "cut_isGlobal_tight"       , isGloMu     ) ) {
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassBool ( "cut_isPf_tight"       , isPfMu     ) ) {
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_chi2_tight"       , chi2  ) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_nMuonHits_tight"  , nHits ) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_nStations_tight" ,nStations) ){ 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_nPixelHits_tight" , nPixHits) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_nTrkLayers_tight" , nTrkLayers ) ){ 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_d0_tight"         , fabs(d0) ) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_z0_tight"         , fabs(z0) ) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_trkiso_tight"     , tkIso/pt       ) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
            if( !config.PassFloat( "cut_corriso_tight"    , pfIso     ) ) { 
                pass_tight = false;
                if( _eval_mu_tight ) continue;
            }
        }

        if( !config.PassBool( "cut_tight", pass_tight) ) continue;

        TLorentzVector mulv;
        mulv.SetPtEtaPhiE( IN::mu_pt->at(idx), 
                           IN::mu_eta->at(idx),
                           IN::mu_phi->at(idx),
                           IN::mu_e->at(idx)
                           );

        float mindr = 101.0;
        for( int hltidx = 0 ; hltidx < IN::HLTObj_n; ++hltidx ) {

            std::vector<int> passTrigs = IN::HLTObj_passTriggers->at(hltidx);

            bool foundMatch = false;
            for( std::vector<int>::const_iterator bitr = _muonTrigMatchBits.begin(); bitr != _muonTrigMatchBits.end(); ++bitr ) {

                if( std::find( passTrigs.begin(), passTrigs.end(), (*bitr) ) != 
                              passTrigs.end() ) foundMatch = true;
            }

            if( !foundMatch ) continue;

            TLorentzVector hltlv;
            hltlv.SetPtEtaPhiE( IN::HLTObj_pt->at(hltidx),
                                IN::HLTObj_eta->at(hltidx),
                                IN::HLTObj_phi->at(hltidx),
                                IN::HLTObj_e->at(hltidx)
                              );

            float dr = mulv.DeltaR( hltlv );
            if( dr < mindr ) {
                mindr = dr;
            }
        }


        bool matchTrig = false;
        //if( mindr < 0.2 ) {
        // change it to 0.1 according to the muon pog
        if( mindr < 0.1 ){ 
            matchTrig = true;
        }
        OUT::mu_hasTrigMatch->push_back(matchTrig);
        OUT::mu_trigMatch_dr->push_back(mindr);

        std::vector<int> matchPID;
        matchPID.push_back(13);
        matchPID.push_back(-13);

        // truth matching is happening in BuildTruth
        //float truthMinDR = 100.0;
        //TLorentzVector matchLV;
        //bool match = HasTruthMatch( mulv, matchPID, 0.2, truthMinDR, matchLV );
        //OUT::mu_hasTruthMatchMu->push_back( match );
        //OUT::mu_truthMatchMu_dr->push_back( truthMinDR );
        //OUT::mu_truthMatchMu_pt->push_back( matchLV.Pt() );
        CopyPrefixIndexBranchesInToOut( "mu_", idx );

        OUT::mu_passTight -> push_back( pass_tight );
        OUT::mu_passLoose -> push_back( pass_loose );

        OUT::mu_n++;

        if( IN::mu_pt->at(idx) > 20 ) {
            OUT::mu_pt20_n++;
        }
        if( IN::mu_pt->at(idx) > 30 ) {
            OUT::mu_pt30_n++;
        }

    }


}

void RunModule::FilterElectron( ModuleConfig & config ) {

    OUT::el_n          = 0;
    OUT::el_pt30_n          = 0;
    OUT::el_passVeryLoose->clear();
    OUT::el_passLoose->clear();
    OUT::el_passMedium->clear();
    OUT::el_passTight->clear();
    OUT::el_hasTrigMatch->clear();
    OUT::el_trigMatch_dr->clear();

    ClearOutputPrefix("el_");

    for( int idx = 0; idx < IN::el_n; ++idx ) {

        float pt    = IN::el_pt->at(idx);
        float eta   = IN::el_eta->at(idx);

        if( !config.PassFloat( "cut_pt",  pt ) ) continue;
        if( !config.PassFloat( "cut_eta", fabs(eta) ) ) continue;

        float sceta = IN::el_sc_eta->at(idx);
        float dEtaIn = IN::el_dEtaIn->at(idx);
        float dPhiIn = IN::el_dPhiIn->at(idx);
        float sigmaIEIEFull5x5 = IN::el_sigmaIEIEfull5x5->at(idx);
        float d0 = IN::el_d0->at(idx);
        float z0 = IN::el_dz->at(idx);
        float hovere = IN::el_hOverE->at(idx);
        float ooEmooP = IN::el_ooEmooP->at(idx);
        float iso_rho = IN::el_pfIsoRho->at(idx);
        bool passConvVeto = IN::el_passConvVeto->at(idx);
        int misshits = IN::el_expectedMissingInnerHits->at(idx);


        bool use_eval = _eval_el_tight || _eval_el_medium || _eval_el_loose || _eval_el_veryloose;

        bool pass_tight = true;
        bool pass_medium = true;
        bool pass_loose = true;
        bool pass_veryloose = true;

        if( fabs(sceta) < 1.479 ) { // barrel

            // Tight cuts
            if( !use_eval || _eval_el_tight ) {
                if( !config.PassFloat( "cut_absdEtaIn_barrel_tight"    , fabs(dEtaIn)       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_barrel_tight"    , fabs(dPhiIn)       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_tight" , sigmaIEIEFull5x5    ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_d0_barrel_tight"        , d0           ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_z0_barrel_tight"        , z0           ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_tight"    , hovere       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_barrel_tight"    , ooEmooP       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_barrel_tight"   , iso_rho   ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_barrel_tight"   , passConvVeto      ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_barrel_tight"  , misshits     ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
            }
            
            // Medium cuts
            if( !use_eval || _eval_el_medium ) {
                if( !config.PassFloat( "cut_absdEtaIn_barrel_medium"    , fabs(dEtaIn)       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_barrel_medium"    , fabs(dPhiIn)       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_medium" , sigmaIEIEFull5x5    ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_d0_barrel_medium"        , d0           ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_z0_barrel_medium"        , z0           ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_medium"    , hovere       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_barrel_medium"    , ooEmooP       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_barrel_medium"   , iso_rho   ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_barrel_medium"   , passConvVeto      ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_barrel_medium"  , misshits     ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
            }
            
            // Loose cuts
            if( !use_eval || _eval_el_loose ) {
                if( !config.PassFloat( "cut_absdEtaIn_barrel_loose"    , fabs(dEtaIn)       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_barrel_loose"    , fabs(dPhiIn)       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_loose" , sigmaIEIEFull5x5    ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_d0_barrel_loose"        , d0           ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_z0_barrel_loose"        , z0           ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_loose"    , hovere       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_barrel_loose"    , ooEmooP       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_barrel_loose"   , iso_rho ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_barrel_loose"   , passConvVeto      ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_barrel_loose"  , misshits     ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
            }

            // Very Loose cuts
            if( !use_eval || _eval_el_veryloose ) {
                if( !config.PassFloat( "cut_absdEtaIn_barrel_veryloose"    , fabs(dEtaIn)       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_barrel_veryloose"    , fabs(dPhiIn)       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_veryloose" , sigmaIEIEFull5x5    ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_d0_barrel_veryloose"        , d0           ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_z0_barrel_veryloose"        , z0           ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_veryloose"    , hovere       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_barrel_veryloose"    , ooEmooP       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_barrel_veryloose"   , iso_rho   ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_barrel_veryloose"   , passConvVeto      ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_barrel_veryloose"  , misshits     ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
            }


        }
        else { // endcap

            // Tight cuts
            if( !use_eval || _eval_el_tight ) {
                if( !config.PassFloat( "cut_absdEtaIn_endcap_tight"    , fabs(dEtaIn)       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_endcap_tight"    , fabs(dPhiIn)       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_tight" , sigmaIEIEFull5x5    ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_d0_endcap_tight"        , d0           ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_z0_endcap_tight"        , z0           ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_tight"    , hovere       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_endcap_tight"    , ooEmooP       ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_endcap_tight"   , iso_rho   ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }

                if( !config.PassInt( "cut_passConvVeto_endcap_tight"   , passConvVeto      ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_endcap_tight"  , misshits     ) ) {
                    pass_tight=false;
                    if( _eval_el_tight ) continue;
                }
            }
            
            // Medium cuts
            if( !use_eval || _eval_el_medium ) {
                if( !config.PassFloat( "cut_absdEtaIn_endcap_medium"    , fabs(dEtaIn)       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_endcap_medium"    , fabs(dPhiIn)       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_medium" , sigmaIEIEFull5x5    ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_d0_endcap_medium"        , d0           ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_z0_endcap_medium"        , z0           ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_medium"    , hovere       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_endcap_medium"    , ooEmooP       ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_endcap_medium"   , iso_rho   ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_endcap_medium"   , passConvVeto      ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_endcap_medium"  , misshits     ) ) {
                    pass_medium=false;
                    if( _eval_el_medium ) continue;
                }
            }
            
            // Loose cuts
            if( !use_eval || _eval_el_loose ) {
                if( !config.PassFloat( "cut_absdEtaIn_endcap_loose"    , fabs(dEtaIn)       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_endcap_loose"    , fabs(dPhiIn)       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_loose" , sigmaIEIEFull5x5    ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_d0_endcap_loose"        , d0           ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_z0_endcap_loose"        , z0           ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_loose"    , hovere       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_endcap_loose"    , ooEmooP       ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_endcap_loose"   , iso_rho   ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_endcap_loose"   , passConvVeto      ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_endcap_loose"  , misshits     ) ) {
                    pass_loose=false;
                    if( _eval_el_loose ) continue;
                }
            }

            // Very Loose cuts
            if( !use_eval || _eval_el_veryloose ) {
                if( !config.PassFloat( "cut_absdEtaIn_endcap_veryloose"    , fabs(dEtaIn)       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_absdPhiIn_endcap_veryloose"    , fabs(dPhiIn)       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_veryloose" , sigmaIEIEFull5x5    ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_d0_endcap_veryloose"        , d0           ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_z0_endcap_veryloose"        , z0           ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_veryloose"    , hovere       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_ooEmooP_endcap_veryloose"    , ooEmooP       ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassFloat( "cut_isoRho_endcap_veryloose"   , iso_rho   ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassInt( "cut_passConvVeto_endcap_veryloose"   , passConvVeto      ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
                if( !config.PassInt  ( "cut_misshits_endcap_veryloose"  , misshits     ) ) {
                    pass_veryloose=false;
                    if( _eval_el_veryloose ) continue;
                }
            }

        }

        if( !config.PassBool( "cut_tight"     , pass_tight ) ) continue;
        if( !config.PassBool( "cut_medium"    , pass_medium ) ) continue;
        if( !config.PassBool( "cut_vid_medium", IN::el_passVIDMedium->at(idx) ) ) continue;

        TLorentzVector ellv;
        ellv.SetPtEtaPhiE( IN::el_pt->at(idx), 
                           IN::el_eta->at(idx),
                           IN::el_phi->at(idx),
                           IN::el_e->at(idx)
                           );

        float min_mu_dr = 100.0;
        for( int muidx = 0; muidx < OUT::mu_n; ++muidx ) {

            TLorentzVector mulv;
            mulv.SetPtEtaPhiE( OUT::mu_pt->at(muidx),
                               OUT::mu_eta->at(muidx),
                               OUT::mu_phi->at(muidx),
                               OUT::mu_e->at(muidx)
                             );

            float dr = mulv.DeltaR( ellv );

            if( dr < min_mu_dr ) {
                min_mu_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_muon_dr", min_mu_dr ) ) continue;

        float mindr = 101.0;
        for( int hltidx = 0 ; hltidx < IN::HLTObj_n; ++hltidx ) {

            std::vector<int> passTrigs = IN::HLTObj_passTriggers->at(hltidx);

            bool foundMatch = false;
            for( std::vector<int>::const_iterator bitr = _electronTrigMatchBits.begin(); bitr != _electronTrigMatchBits.end(); ++bitr ) {
                if( std::find( passTrigs.begin(), passTrigs.end(), (*bitr) ) != 
                              passTrigs.end() ) foundMatch = true;
            }

            if( !foundMatch ) continue;

            TLorentzVector hltlv;
            hltlv.SetPtEtaPhiE( IN::HLTObj_pt->at(hltidx),
                                IN::HLTObj_eta->at(hltidx),
                                IN::HLTObj_phi->at(hltidx),
                                IN::HLTObj_e->at(hltidx)
                              );

            float dr = ellv.DeltaR( hltlv );
            if( dr < mindr ) {
                mindr = dr;
            }
        }

        bool matchTrig = false;
        if( mindr < 0.2 ) {
            matchTrig = true;
        }
        OUT::el_hasTrigMatch->push_back(matchTrig);
        OUT::el_trigMatch_dr->push_back(mindr);

        OUT::el_passVeryLoose->push_back( pass_veryloose );
        OUT::el_passLoose->push_back( pass_loose );
        OUT::el_passMedium->push_back( pass_medium );
        OUT::el_passTight->push_back( pass_tight );

        CopyPrefixIndexBranchesInToOut( "el_", idx );
        OUT::el_n++;

        if( IN::el_pt->at(idx) > 30 ) {
            OUT::el_pt30_n++;
        }

    }

}

void RunModule::FilterPhoton( ModuleConfig & config ) {

    OUT::ph_n          = 0;
    OUT::ph_chIsoCorr            -> clear();
    OUT::ph_neuIsoCorr           -> clear();
    OUT::ph_phoIsoCorr           -> clear();
    OUT::ph_min_el_dr           -> clear();
    OUT::ph_IsEB                 -> clear();
    OUT::ph_IsEE                 -> clear();
    OUT::ph_passTight            -> clear();
    OUT::ph_passMedium           -> clear();
    OUT::ph_passLoose            -> clear();
    OUT::ph_passLooseNoSIEIE     -> clear();
    OUT::ph_passHOverELoose      -> clear();
    OUT::ph_passHOverEMedium     -> clear();
    OUT::ph_passHOverETight      -> clear();
    OUT::ph_passSIEIELoose       -> clear();
    OUT::ph_passSIEIEMedium      -> clear();
    OUT::ph_passSIEIETight       -> clear();
    OUT::ph_passChIsoCorrLoose   -> clear();
    OUT::ph_passChIsoCorrMedium  -> clear();
    OUT::ph_passChIsoCorrTight   -> clear();
    OUT::ph_passNeuIsoCorrLoose  -> clear();
    OUT::ph_passNeuIsoCorrMedium -> clear();
    OUT::ph_passNeuIsoCorrTight  -> clear();
    OUT::ph_passPhoIsoCorrLoose  -> clear();
    OUT::ph_passPhoIsoCorrMedium -> clear();
    OUT::ph_passPhoIsoCorrTight  -> clear();


    ClearOutputPrefix("ph_");

    for( int idx = 0; idx < IN::ph_n; ++idx ) {

        float pt = IN::ph_pt->at(idx);
        float pt_uncalib = IN::ph_ptOrig->at(idx);
        float eta = IN::ph_eta->at(idx);
        float sceta = IN::ph_sc_eta->at(idx);

        bool iseb = false;
        bool isee = false;
        if( fabs(sceta) < 1.479 ) {
            iseb = true;
        }
        if( fabs(sceta) > 1.57 ) {
            isee = true;
        }

        if( !config.PassFloat( "cut_pt", pt ) ) continue;
        if( !config.PassFloat( "cut_eta", fabs(eta) ) ) continue;
        if( !config.PassBool( "cut_eb", iseb ) ) continue;
        if( !config.PassBool( "cut_ee", isee ) ) continue;

        float sigmaIEIE = IN::ph_sigmaIEIEFull5x5->at(idx);
        float hovere = IN::ph_hOverE->at(idx);

        float pfChIso  = IN::ph_chIso->at(idx);
        float pfNeuIso = IN::ph_neuIso->at(idx);
        float pfPhoIso = IN::ph_phoIso->at(idx);

        float pfChIsoRhoCorr = 0.0;
        float pfNeuIsoRhoCorr = 0.0;
        float pfPhoIsoRhoCorr = 0.0;
        calc_corr_iso( pfChIso, pfPhoIso, pfNeuIso, IN::rho, sceta, pfChIsoRhoCorr, pfPhoIsoRhoCorr, pfNeuIsoRhoCorr);

        float p1_neu = 0;
        float p2_neu = 0;
        float p1_pho = 0;
        // taken from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recommended_Working_points_for_2
        // Updated Dec 2016
        if( iseb ) {
            p1_neu = 0.0148;
            p2_neu = 0.000017;
            p1_pho = 0.0047;
        }
        else {
            p1_neu = 0.0163;
            p2_neu = 0.000014;
            p1_pho = 0.0034;
        }

        float pfChIsoPtRhoCorr  = pfChIsoRhoCorr;
        float pfNeuIsoPtRhoCorr = pfNeuIsoRhoCorr-p1_neu*pt_uncalib-p2_neu*pt_uncalib*pt_uncalib;
        float pfPhoIsoPtRhoCorr = pfPhoIsoRhoCorr-p1_pho*pt_uncalib;

        bool pass_loose         = true;
        bool pass_loose_nosieie = true;
        bool pass_medium        = true;
        bool pass_tight         = true;

        bool pass_sieie_loose   = true;
        bool pass_sieie_medium   = true;
        bool pass_sieie_tight   = true;

        bool pass_chIsoCorr_loose   = true;
        bool pass_chIsoCorr_medium   = true;
        bool pass_chIsoCorr_tight   = true;

        bool pass_neuIsoCorr_loose   = true;
        bool pass_neuIsoCorr_medium   = true;
        bool pass_neuIsoCorr_tight   = true;

        bool pass_phoIsoCorr_loose   = true;
        bool pass_phoIsoCorr_medium   = true;
        bool pass_phoIsoCorr_tight   = true;

        bool pass_hovere_loose   = true;
        bool pass_hovere_medium   = true;
        bool pass_hovere_tight   = true;

        bool use_eval = _eval_ph_tight || _eval_ph_medium || _eval_ph_loose;

        if( fabs(sceta) < 1.479 ) { // barrel

            //loose 
            if( !use_eval || _eval_ph_loose ) {
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_loose"   , sigmaIEIE         ) ) {
                    pass_loose=false;
                    pass_sieie_loose=false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_chIsoCorr_barrel_loose"   , pfChIsoPtRhoCorr  ) ) {
                    pass_loose=false;
                    pass_chIsoCorr_loose = false;
                    pass_loose_nosieie=false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_neuIsoCorr_barrel_loose"  , pfNeuIsoPtRhoCorr ) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_neuIsoCorr_loose = false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_phoIsoCorr_barrel_loose"  , pfPhoIsoPtRhoCorr ) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_phoIsoCorr_loose = false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_loose"  , hovere) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_hovere_loose = false;
                    if( _eval_ph_loose ) continue;
                }

            }

            // medium
            if( !use_eval || _eval_ph_medium ) {
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_medium"  , sigmaIEIE         ) ) {
                    pass_medium=false;
                    pass_sieie_medium=false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_chIsoCorr_barrel_medium"  , pfChIsoPtRhoCorr  ) ) {
                    pass_medium=false;
                    pass_chIsoCorr_medium = false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_neuIsoCorr_barrel_medium" , pfNeuIsoPtRhoCorr ) ) {
                    pass_medium=false;
                    pass_neuIsoCorr_medium = false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_phoIsoCorr_barrel_medium" , pfPhoIsoPtRhoCorr ) ) {
                    pass_medium=false;
                    pass_phoIsoCorr_medium = false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_medium" , hovere) ) {
                    pass_medium=false;
                    pass_hovere_medium = false;
                    if( _eval_ph_medium ) continue;
                }
            }

            // tight
            if( !use_eval || _eval_ph_tight ) {
                if( !config.PassFloat( "cut_sigmaIEIE_barrel_tight"   , sigmaIEIE         ) ) {
                    pass_tight=false;
                    pass_sieie_tight=false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_chIsoCorr_barrel_tight"   , pfChIsoPtRhoCorr  ) ) {
                    pass_tight=false;
                    pass_chIsoCorr_tight = false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_neuIsoCorr_barrel_tight"  , pfNeuIsoPtRhoCorr ) ) {
                    pass_tight=false;
                    pass_neuIsoCorr_tight = false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_phoIsoCorr_barrel_tight"  , pfPhoIsoPtRhoCorr ) ) {
                    pass_tight=false;
                    pass_phoIsoCorr_tight = false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_hovere_barrel_tight"  , hovere) ) {
                    pass_tight=false;
                    pass_hovere_tight = false;
                    if( _eval_ph_tight ) continue;
                }
            }

        }
        else { // endcap
            // loose
            if( !use_eval || _eval_ph_loose ) {
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_loose"   , sigmaIEIE         ) ) {
                    pass_loose=false;
                    pass_sieie_loose=false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_chIsoCorr_endcap_loose"   , pfChIsoPtRhoCorr  ) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_chIsoCorr_loose = false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_neuIsoCorr_endcap_loose"  , pfNeuIsoPtRhoCorr ) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_neuIsoCorr_loose = false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_phoIsoCorr_endcap_loose"  , pfPhoIsoPtRhoCorr ) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_phoIsoCorr_loose = false;
                    if( _eval_ph_loose ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_loose"  , hovere) ) {
                    pass_loose=false;
                    pass_loose_nosieie=false;
                    pass_hovere_loose = false;
                    if( _eval_ph_loose ) continue;
                }
            }

            // medium
            if( !use_eval || _eval_ph_medium ) {
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_medium"  , sigmaIEIE         ) ) {
                    pass_medium=false;
                    pass_sieie_medium=false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_chIsoCorr_endcap_medium"  , pfChIsoPtRhoCorr  ) ) {
                    pass_medium=false;
                    pass_chIsoCorr_medium = false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_neuIsoCorr_endcap_medium" , pfNeuIsoPtRhoCorr ) ) {
                    pass_medium=false;
                    pass_neuIsoCorr_medium = false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_phoIsoCorr_endcap_medium" , pfPhoIsoPtRhoCorr ) ) {
                    pass_medium=false;
                    pass_phoIsoCorr_medium = false;
                    if( _eval_ph_medium ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_medium" , hovere) ) {
                    pass_medium=false;
                    pass_hovere_medium = false;
                    if( _eval_ph_medium ) continue;
                }
            }

            // tight
            if( !use_eval || _eval_ph_tight ) {
                if( !config.PassFloat( "cut_sigmaIEIE_endcap_tight"   , sigmaIEIE         ) ) {
                    pass_tight=false;
                    pass_sieie_tight=false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_chIsoCorr_endcap_tight"   , pfChIsoPtRhoCorr  ) ) {
                    pass_tight=false;
                    pass_chIsoCorr_tight = false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_neuIsoCorr_endcap_tight"  , pfNeuIsoPtRhoCorr ) ) {
                    pass_tight=false;
                    pass_neuIsoCorr_tight = false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_phoIsoCorr_endcap_tight"  , pfPhoIsoPtRhoCorr ) ) {
                    pass_phoIsoCorr_tight = false;
                    pass_tight=false;
                    if( _eval_ph_tight ) continue;
                }
                if( !config.PassFloat( "cut_hovere_endcap_tight"  , hovere) ) {
                    pass_hovere_tight = false;
                    pass_tight=false;
                    if( _eval_ph_tight ) continue;
                }
            }
        }

        if( !config.PassBool( "cut_tight"    , pass_tight     ) ) continue;
        if( !config.PassBool( "cut_medium"   , pass_medium    ) ) continue;
        if( !config.PassBool( "cut_loose"    , pass_loose     ) ) continue;
        if( !config.PassBool( "cut_CSEV"     , IN::ph_passEleVeto->at(idx) ) ) continue;

        TLorentzVector phlv;
        phlv.SetPtEtaPhiE( IN::ph_pt->at(idx), 
                           IN::ph_eta->at(idx),
                           IN::ph_phi->at(idx),
                           IN::ph_e->at(idx) 
                           );

        float min_mu_dr = 100.0;
        for( int muidx = 0; muidx < OUT::mu_n; ++muidx ) {

            TLorentzVector mulv;
            mulv.SetPtEtaPhiE( OUT::mu_pt->at(muidx),
                               OUT::mu_eta->at(muidx),
                               OUT::mu_phi->at(muidx),
                               OUT::mu_e->at(muidx)
                             );

            float dr = mulv.DeltaR( phlv );

            if( dr < min_mu_dr ) {
                min_mu_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_muon_dr", min_mu_dr ) ) continue;

        float min_el_dr = 100.0;
        for( int elidx = 0; elidx < OUT::el_n; ++elidx ) {

            TLorentzVector ellv;
            ellv.SetPtEtaPhiE( OUT::el_pt->at(elidx),
                               OUT::el_eta->at(elidx),
                               OUT::el_phi->at(elidx),
                               OUT::el_e->at(elidx)
                             );

            float dr = ellv.DeltaR( phlv );

            if( dr < min_el_dr ) {
                min_el_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_electron_dr", min_el_dr ) ) continue;


        OUT::ph_chIsoCorr            -> push_back(pfChIsoPtRhoCorr);
        OUT::ph_neuIsoCorr           -> push_back(pfNeuIsoPtRhoCorr);
        OUT::ph_phoIsoCorr           -> push_back(pfPhoIsoPtRhoCorr);

        OUT::ph_min_el_dr            -> push_back( min_el_dr );

        OUT::ph_passTight            -> push_back(pass_tight);
        OUT::ph_passMedium           -> push_back(pass_medium);
        OUT::ph_passLoose            -> push_back(pass_loose);
        OUT::ph_passLooseNoSIEIE     -> push_back(pass_loose_nosieie);
        OUT::ph_passHOverELoose      -> push_back(pass_hovere_loose);
        OUT::ph_passHOverEMedium     -> push_back(pass_hovere_medium);
        OUT::ph_passHOverETight      -> push_back(pass_hovere_tight);
        OUT::ph_passSIEIELoose       -> push_back(pass_sieie_loose);
        OUT::ph_passSIEIEMedium      -> push_back(pass_sieie_medium);
        OUT::ph_passSIEIETight       -> push_back(pass_sieie_tight);
        OUT::ph_passChIsoCorrLoose   -> push_back(pass_chIsoCorr_loose);
        OUT::ph_passChIsoCorrMedium  -> push_back(pass_chIsoCorr_medium);
        OUT::ph_passChIsoCorrTight   -> push_back(pass_chIsoCorr_tight);
        OUT::ph_passNeuIsoCorrLoose  -> push_back(pass_neuIsoCorr_loose);
        OUT::ph_passNeuIsoCorrMedium -> push_back(pass_neuIsoCorr_medium);
        OUT::ph_passNeuIsoCorrTight  -> push_back(pass_neuIsoCorr_tight);
        OUT::ph_passPhoIsoCorrLoose  -> push_back(pass_phoIsoCorr_loose);
        OUT::ph_passPhoIsoCorrMedium -> push_back(pass_phoIsoCorr_medium);
        OUT::ph_passPhoIsoCorrTight  -> push_back(pass_phoIsoCorr_tight);

        OUT::ph_IsEB -> push_back( iseb );
        OUT::ph_IsEE -> push_back( isee );


        CopyPrefixIndexBranchesInToOut( "ph_", idx );
        OUT::ph_n++;

    }
}

void RunModule::calc_corr_iso( float chIso, float phoIso, float neuIso, float rho, float eta, float &chIsoCorr, float &phoIsoCorr, float &neuIsoCorr )  const
{

    // from https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
    // updated Dec 2016
    
    float ea_ch=0.0;
    float ea_pho=0.0;
    float ea_neu=0.0;

    if( fabs( eta ) < 1.0 ) {
        ea_ch = 0.0360;
        ea_neu = 0.0597;
        ea_pho = 0.1210;
    }
    else if( fabs(eta) >= 1.0 && fabs( eta ) < 1.479 ) {
        ea_ch = 0.0377;
        ea_neu = 0.0807;
        ea_pho = 0.1107;
    }
    else if( fabs(eta) >= 1.479 && fabs( eta ) < 2.0 ) {
        ea_ch = 0.0306;
        ea_neu = 0.0629;
        ea_pho = 0.0699;
    }
    else if( fabs(eta) >= 2.0 && fabs( eta ) < 2.2 ) {
        ea_ch = 0.0283;
        ea_neu = 0.0197;
        ea_pho = 0.1056;
    }
    else if( fabs(eta) >= 2.2 && fabs( eta ) < 2.3 ) {
        ea_ch = 0.0254;
        ea_neu = 0.0184;
        ea_pho = 0.1457;
    }
    else if( fabs(eta) >= 2.3 && fabs( eta ) < 2.4 ) {
        ea_ch = 0.0217;
        ea_neu = 0.0284;
        ea_pho = 0.1719;
    }
    else if( fabs(eta) >= 2.4 ) {
        ea_ch = 0.0167;
        ea_neu = 0.0591;
        ea_pho = 0.1998;
    }

    chIsoCorr  = chIso  - rho*ea_ch;
    phoIsoCorr = phoIso - rho*ea_pho;
    neuIsoCorr = neuIso - rho*ea_neu;

    if( chIsoCorr < 0 ) {
        chIsoCorr = 0;
    }
    if( phoIsoCorr < 0 ) {
        phoIsoCorr = 0;
    }
    if( neuIsoCorr < 0 ) {
        neuIsoCorr = 0;
    }

}


void RunModule::FilterJet( ModuleConfig & config ) const {

    OUT::jet_n          = 0;
    OUT::jet_CSVLoose_n          = 0;
    OUT::jet_CSVMedium_n          = 0;
    OUT::jet_CSVTight_n          = 0;
    OUT::jet_IdLoose        -> clear();
    OUT::jet_IdTight        -> clear();
    OUT::jet_IdTightLep     -> clear();
    ClearOutputPrefix("jet_");

    for( int idx = 0; idx < IN::jet_n; ++idx ) {

        float pt = IN::jet_pt->at(idx);
        float eta = IN::jet_eta->at(idx);

        if( !config.PassFloat( "cut_pt", pt ) ) continue;
        if( !config.PassFloat( "cut_eta", fabs(eta)) ) continue;

        float nhf        = IN::jet_nhf->at(idx);
        float chf        = IN::jet_chf->at(idx);
        float muf        = IN::jet_muf->at(idx);
        float cemf       = IN::jet_cemf->at(idx);
        float nemf       = IN::jet_nemf->at(idx);
        int   cmult      = IN::jet_cmult->at(idx);
        int   nmult      = IN::jet_nmult->at(idx);
        int   nconst     = cmult + nmult;

        bool pass_loose    = true;
        bool pass_tight    = true;
        bool pass_tightlep = true;

        if( fabs( eta ) < 2.7 ) {
            if( !config.PassFloat( "cut_jet_nhf_central_loose" , nhf  ) ) pass_loose = false;
            if( !config.PassFloat( "cut_jet_nemf_central_loose", nemf ) ) pass_loose = false;
            if( !config.PassInt  ( "cut_jet_nconst_central_loose", nconst ) ) pass_loose = false;

            if( !config.PassFloat( "cut_jet_nhf_central_tight" , nhf  ) ) pass_tight = false;
            if( !config.PassFloat( "cut_jet_nemf_central_tight", nemf ) ) pass_tight = false;
            if( !config.PassInt  ( "cut_jet_nconst_central_tight", nconst ) ) pass_tight = false;

            if( !config.PassFloat( "cut_jet_nhf_central_tightlep" , nhf  ) ) pass_tightlep = false;
            if( !config.PassFloat( "cut_jet_nemf_central_tightlep", nemf ) ) pass_tightlep = false;
            if( !config.PassInt  ( "cut_jet_nconst_central_tightlep", nconst ) ) pass_tightlep = false;
            if( !config.PassFloat  ( "cut_jet_muf_central_tightlep", muf ) ) pass_tightlep = false;
            if( fabs( eta ) < 2.4 ) {
                if( !config.PassFloat( "cut_jet_chf_central_loose" , chf  ) ) pass_loose = false;
                if( !config.PassInt( "cut_jet_cmult_central_loose" , cmult  ) ) pass_loose = false;
                if( !config.PassFloat( "cut_jet_cemf_central_loose" , cemf ) ) pass_loose = false;

                if( !config.PassFloat( "cut_jet_chf_central_tight" , chf  ) ) pass_tight = false;
                if( !config.PassInt( "cut_jet_cmult_central_tight" , cmult  ) ) pass_tight = false;
                if( !config.PassFloat( "cut_jet_cemf_central_tight" , cemf ) ) pass_tight = false;

                if( !config.PassFloat( "cut_jet_chf_central_tightlep" , chf  ) ) pass_tightlep = false;
                if( !config.PassInt( "cut_jet_cmult_central_tightlep" , cmult  ) ) pass_tightlep = false;
                if( !config.PassFloat( "cut_jet_cemf_central_tightlep" , cemf ) ) pass_tightlep = false;
            }

        }
        else if( fabs( eta ) < 3.0 ) {
            if( !config.PassFloat( "cut_jet_nhf_transition_loose", nhf ) ) pass_loose = false;
            if( !config.PassFloat( "cut_jet_nemf_transition_loose", nemf ) ) pass_loose = false;
            if( !config.PassInt( "cut_jet_nmult_transition_loose", nmult ) ) pass_loose = false;

            if( !config.PassFloat( "cut_jet_nhf_transition_tight", nhf ) ) pass_tight = false;
            if( !config.PassFloat( "cut_jet_nemf_transition_tight", nemf ) ) pass_tight = false;
            if( !config.PassInt( "cut_jet_nmult_transition_tight", nmult ) ) pass_tight = false;

        }
        else {
            if( !config.PassFloat( "cut_jet_nemf_forward_loose", nemf ) ) pass_loose = false;
            if( !config.PassInt( "cut_jet_nmult_forward_loose", nmult ) ) pass_loose = false;

            if( !config.PassFloat( "cut_jet_nemf_forward_tight", nemf ) ) pass_tight = false;
            if( !config.PassInt( "cut_jet_nmult_forward_tight", nmult ) ) pass_tight = false;
        }

        if( !config.PassBool( "cut_tight", pass_tight) ) continue;
        if( !config.PassBool( "cut_loose", pass_loose) ) continue;

        TLorentzVector jetlv;
        jetlv.SetPtEtaPhiE( IN::jet_pt->at(idx), 
                            IN::jet_eta->at(idx),
                            IN::jet_phi->at(idx),
                            IN::jet_e->at(idx) 
                            );

        float min_mu_dr = 100.0;
        for( int muidx = 0; muidx < OUT::mu_n; ++muidx ) {

            TLorentzVector mulv;
            mulv.SetPtEtaPhiE( OUT::mu_pt->at(muidx),
                               OUT::mu_eta->at(muidx),
                               OUT::mu_phi->at(muidx),
                               OUT::mu_e->at(muidx)
                             );

            float dr = mulv.DeltaR( jetlv );

            if( dr < min_mu_dr ) {
                min_mu_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_muon_dr", min_mu_dr ) ) continue;

        float min_el_dr = 100.0;
        for( int elidx = 0; elidx < OUT::el_n; ++elidx ) {

            TLorentzVector ellv;
            ellv.SetPtEtaPhiE( OUT::el_pt->at(elidx),
                               OUT::el_eta->at(elidx),
                               OUT::el_phi->at(elidx),
                               OUT::el_e->at(elidx)
                             );

            float dr = ellv.DeltaR( jetlv );

            if( dr < min_el_dr ) {
                min_el_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_electron_dr", min_el_dr ) ) continue;

        float min_ph_dr = 100.0;
        for( int phidx = 0; phidx < OUT::ph_n; ++phidx ) {

            TLorentzVector phlv;
            phlv.SetPtEtaPhiE( OUT::ph_pt->at(phidx),
                               OUT::ph_eta->at(phidx),
                               OUT::ph_phi->at(phidx),
                               OUT::ph_e->at(phidx)
                             );

            float dr = phlv.DeltaR( jetlv );

            if( dr < min_ph_dr ) {
                min_ph_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_photon_dr", min_ph_dr ) ) continue;

        CopyPrefixIndexBranchesInToOut( "jet_", idx );
        OUT::jet_n++;

        OUT::jet_IdLoose        -> push_back(pass_loose );
        OUT::jet_IdTight        -> push_back(pass_tight );
        OUT::jet_IdTightLep     -> push_back(pass_tightlep );

        float jet_csv = IN::jet_bTagCisvV2->at(idx);

        if( config.PassFloat( "cut_jet_CSV_Loose",  jet_csv ) ) {
            OUT::jet_CSVLoose_n++;
        }
        if( config.PassFloat( "cut_jet_CSV_Medium",  jet_csv ) ) {
            OUT::jet_CSVMedium_n++;
        }
        if( config.PassFloat( "cut_jet_CSV_Tight",  jet_csv ) ) {
            OUT::jet_CSVTight_n++;
        }

    }
}

bool RunModule::FilterEvent( ModuleConfig & config ) const {

    bool keep_event = true;

    if( !config.PassInt( "cut_el_n"     , OUT::el_n        ) ) keep_event=false;
    if( !config.PassInt( "cut_el_pt30_n", OUT::el_pt30_n   ) ) keep_event=false;
    if( !config.PassInt( "cut_mu_n"     , OUT::mu_n        ) ) keep_event=false;
    if( !config.PassInt( "cut_mu_pt30_n", OUT::mu_pt30_n   ) ) keep_event=false;
    if( !config.PassInt( "cut_mu_pt20_n", OUT::mu_pt20_n   ) ) keep_event=false;
    if( !config.PassInt( "cut_ph_n"     , OUT::ph_n   ) ) keep_event=false;
    if( !config.PassInt( "cut_jet_n"    , OUT::jet_n  ) ) keep_event=false;
    
    //if( !config.PassBool( "cut_trig_Ele27_eta2p1_tight", IN::passTrig_HLT_Ele27_eta2p1_WPTight_Gsf) ) keep_event=false;
    //if( !config.PassBool( "cut_trig_Mu27_IsoORIsoTk", (IN::passTrig_HLT_IsoMu27 | IN::passTrig_HLT_IsoTkMu27) ) ) keep_event=false;
    //if( !config.PassBool( "cut_trig_Mu24_IsoORIsoTk", (IN::passTrig_HLT_IsoMu24 | IN::passTrig_HLT_IsoTkMu24) ) ) keep_event=false;
    //if( !config.PassBool( "cut_trig_Mu24_IsoORIsoTk", (IN::passTrig_HLT_IsoMu24 | IN::passTrig_HLT_IsoTkMu24) ) ) keep_event=false;
    //if( !config.PassBool( "cut_trig_Mu17_Photon30", (IN::passTrig_HLT_Mu17_Photon30_CaloIdL_L1ISO) ) ) keep_event=false;

    return keep_event;
    
}

bool RunModule::FilterTrigger( ModuleConfig & config ) {

    std::vector<int> passed_ids;
    // for each configured trigger, get its decision and store the result 
    // in the output branch
    for( std::map<int, bool>::iterator itr = triggerResults.begin();
            itr != triggerResults.end() ; ++itr ) {

        if( std::find( IN::passedTriggers->begin(), IN::passedTriggers->end(), itr->first ) 
                != IN::passedTriggers->end() ) {
            itr->second = true;
            passed_ids.push_back( itr->first );
        }
        else {
            itr->second = false;
        }
    }

    bool keep_event = true;

    if( !config.PassAnyIntVector( "cut_bits", passed_ids ) ) keep_event = false;

    return keep_event;
    
}

bool RunModule::FilterMET( ModuleConfig & config ) {

    std::vector<int> passed_ids;
    // for each configured metfilter, get its decision and store the result
    // in the output branch
    for( std::map<int, bool>::iterator itr = metfilterResults.begin();
            itr != metfilterResults.end() ; ++itr ) {

        if( std::find( IN::passedFilters->begin(), IN::passedFilters->end(), itr->first )
                != IN::passedFilters->end() ) {
            itr->second = true;
            passed_ids.push_back( itr->first );
        }
        else {
            itr->second = false;
        }
    }

    bool keep_event = true;


    if( !config.PassAnyIntVector( "cut_metfilter_bits", passed_ids ) ) keep_event = false;

    return keep_event;

}

void RunModule::BuildEventVars( ModuleConfig & config ) const {


    OUT::m_lep_ph = 0;
    OUT::m_lep_ph_comb_leadLep->clear();
    OUT::m_lep_ph_comb_sublLep->clear();
    OUT::m_lep_met_ph = 0;
    OUT::m_mt_lep_met_ph = 0;
    OUT::m_mt_lep_met_ph_forcewmass = 0;
    OUT::mt_w = 0;
    OUT::mt_res = 0;
    OUT::mt_lep_ph = 0;
    OUT::dphi_lep_ph = 0;
    OUT::dr_lep_ph = 0;
    OUT::m_lep_met = 0;
    OUT::mt_lep_met = 0;
    OUT::pt_lep_met = 0;
    OUT::dphi_lep_met = 0;
    OUT::mt_lep_met_ph = 0;
    OUT::mt_lep_met_ph_inv = 0;
    OUT::RecoWMass = 0;
    OUT::recoM_lep_nu_ph = 0;
    OUT::recoMet_eta= 0;
    OUT::recoW_pt= 0;
    OUT::recoW_eta= 0;
    OUT::recoW_phi= 0;
    OUT::m_ll = 0;

    OUT::leadjet_pt = 0;
    OUT::subljet_pt = 0;
    OUT::leaddijet_m = 0;
    OUT::leaddijet_pt = 0;
    OUT::massdijet_m = 0;
    OUT::massdijet_pt = 0;

    std::vector<TLorentzVector> leptons;
    std::vector<TLorentzVector> photons;

    for( int idx = 0; idx < OUT::mu_n; ++idx ) {
        TLorentzVector tlv;
        tlv.SetPtEtaPhiE( OUT::mu_pt->at(idx), 
                          OUT::mu_eta->at(idx), 
                          OUT::mu_phi->at(idx), 
                          OUT::mu_e->at(idx) );

        leptons.push_back(tlv);
    }

    for( int idx = 0; idx < OUT::el_n; ++idx ) {
        TLorentzVector tlv;
        tlv.SetPtEtaPhiE( OUT::el_pt->at(idx), 
                          OUT::el_eta->at(idx), 
                          OUT::el_phi->at(idx), 
                          OUT::el_e->at(idx) );

        leptons.push_back(tlv);
    }


    for( int idx = 0; idx < OUT::ph_n; ++idx ) {
        TLorentzVector tlv;
        tlv.SetPtEtaPhiE( OUT::ph_pt->at(idx), 
                          OUT::ph_eta->at(idx), 
                          OUT::ph_phi->at(idx), 
                          OUT::ph_e->at(idx) );

        photons.push_back(tlv);
    }

    TLorentzVector metlv;
    metlv.SetPtEtaPhiM( OUT::met_pt, 0.0, OUT::met_phi, 0.0 );
    TLorentzVector metlvOrig( metlv );

    if( leptons.size() > 0 ) {
        OUT::mt_lep_met = Utils::calc_mt( leptons[0], metlvOrig );
        OUT::m_lep_met = (leptons[0]+metlvOrig).M();
        OUT::pt_lep_met = (leptons[0]+metlvOrig).Pt();
        OUT::dphi_lep_met = leptons[0].DeltaPhi( metlvOrig );

        bool success = get_constriained_nu_pz( leptons[0], metlv );
        OUT::nu_z_solution_success = success;

        OUT::RecoWMass = ( leptons[0] + metlv ).M();

        if( leptons.size() > 1 ) {
            OUT::m_ll = (leptons[0] + leptons[1]).M();
        }

    }

    if( photons.size() > 0 ) {

        if( leptons.size() > 0 ) {

            OUT::m_lep_ph = ( leptons[0] + photons[0] ).M();
            OUT::m_lep_met_ph = ( leptons[0] + photons[0] + metlvOrig ).M();
            OUT::dphi_lep_ph = leptons[0].DeltaPhi(photons[0] );
            OUT::dr_lep_ph = leptons[0].DeltaR(photons[0] );
            OUT::recoM_lep_nu_ph = ( leptons[0] + metlv + photons[0] ).M();
            OUT::recoMet_eta = metlv.Eta() ;
            OUT::mt_lep_met_ph = Utils::calc_mt( leptons[0] + metlvOrig, photons[0]);
            OUT::mt_lep_met_ph_inv = Utils::calc_mt( leptons[0] + photons[0], metlvOrig);
            OUT::dphi_ph_w = ( metlvOrig + leptons[0] ).DeltaPhi( photons[0] );
            OUT::pt_lep_met_ph = ( metlvOrig + leptons[0] + photons[0]).Pt();

            TLorentzVector recoW = leptons[0] + metlv;

            OUT::recoW_pt = recoW.Pt() ;
            OUT::recoW_eta = recoW.Eta() ;
            OUT::recoW_phi = recoW.Phi() ;

            float mt = Utils::calc_mt( leptons[0], metlvOrig );

            TLorentzVector wlv;
            wlv.SetXYZM( leptons[0].Px() + metlvOrig.Px(), leptons[0].Py() + metlvOrig.Py(), leptons[0].Pz(), mt );

            TLorentzVector wlv_force;
            wlv_force.SetXYZM( leptons[0].Px() + metlvOrig.Px(), leptons[0].Py() + metlvOrig.Py(), leptons[0].Pz(), _m_w );
            OUT::m_mt_lep_met_ph = ( wlv + photons[0] ).M();
            OUT::m_mt_lep_met_ph_forcewmass = ( wlv_force + photons[0] ).M();
            OUT::mt_w = mt;

            
            TLorentzVector lep_trans; 
            TLorentzVector ph_trans; 
            lep_trans.SetPtEtaPhiM( leptons[0].Pt(), 0.0, leptons[0].Phi(), leptons[0].M() );
            ph_trans.SetPtEtaPhiM( photons[0].Pt(), 0.0, photons[0].Phi(), photons[0].M() );

            OUT::mt_res = ( lep_trans + ph_trans + metlvOrig ).M();
            OUT::mt_lep_ph = ( lep_trans + ph_trans ).M();

            for( std::vector<TLorentzVector>::const_iterator phitr = photons.begin();
                    phitr != photons.end(); ++phitr ) {

                float mass = ( *phitr + leptons[0] ).M();

                OUT::m_lep_ph_comb_leadLep->push_back( mass );
            }

            if( leptons.size() == 2 ) {
                for( std::vector<TLorentzVector>::const_iterator phitr = photons.begin();
                        phitr != photons.end(); ++phitr ) {

                    float mass = ( *phitr + leptons[1] ).M();

                    OUT::m_lep_ph_comb_sublLep->push_back( mass );
                }
            }
        }
    }

    std::vector<std::pair<float, int> > sorted_jets;
    std::vector<TLorentzVector> jet_lvs;

    for( int jeti = 0; jeti < OUT::jet_n; ++jeti ) {
        sorted_jets.push_back( std::make_pair( OUT::jet_pt->at(jeti ), jeti ) );

        TLorentzVector lv;
        lv.SetPtEtaPhiE( OUT::jet_pt->at(jeti),
                        OUT::jet_eta->at(jeti),
                        OUT::jet_phi->at(jeti),
                        OUT::jet_e->at(jeti)
                );
        jet_lvs.push_back( lv );
    }


    std::sort(sorted_jets.rbegin(), sorted_jets.rend());

    if( OUT::jet_n > 0 ) {
        OUT::leadjet_pt = OUT::jet_pt->at(0);

        if( OUT::jet_n > 1 ) {
            OUT::leadjet_pt = OUT::jet_pt->at(sorted_jets[0].second);
            OUT::subljet_pt = OUT::jet_pt->at(sorted_jets[1].second);

            OUT::leaddijet_m  = (jet_lvs[sorted_jets[0].second]+jet_lvs[sorted_jets[1].second]).M();
            OUT::leaddijet_pt = (jet_lvs[sorted_jets[0].second]+jet_lvs[sorted_jets[1].second]).Pt();

            float min_mass = 100000000.;
            int min_idx1 = -1;
            int min_idx2 = -1;
            for( unsigned i = 0 ; i < jet_lvs.size(); ++i ) {
                for( unsigned j = i+1 ; j < jet_lvs.size(); ++j ) {

                    float mass = ( jet_lvs[i] + jet_lvs[j] ).M();
                    float diff = fabs( 91.2 - mass );

                    if( diff < min_mass ) {
                        min_mass = diff;
                        min_idx1 = i;
                        min_idx2 = j;
                    }
                }
            }

            OUT::massdijet_m  = ( jet_lvs[min_idx1] + jet_lvs[min_idx2] ).M();
            OUT::massdijet_pt = ( jet_lvs[min_idx1] + jet_lvs[min_idx2] ).Pt();
        }

    }

}

void RunModule::MakePhotonCountVars( ModuleConfig & config ) const { 

    OUT::ph_loose_n = 0;
    OUT::ph_medium_n = 0;
    OUT::ph_tight_n = 0;

    OUT::ph_mediumPassPSV_n = 0;
    OUT::ph_mediumFailPSV_n = 0;
    OUT::ph_mediumPassCSEV_n = 0;
    OUT::ph_mediumFailCSEV_n = 0;

    OUT::ph_mediumPassEleOlap_n = 0;
    OUT::ph_mediumPassEleOlapPassCSEV_n = 0;
    OUT::ph_mediumPassEleOlapFailCSEV_n = 0;

    OUT::ph_mediumNoSIEIE_n = 0;
    OUT::ph_mediumNoChIso_n = 0;
    OUT::ph_mediumNoNeuIso_n = 0;
    OUT::ph_mediumNoPhoIso_n = 0;

    OUT::ph_mediumNoSIEIENoChIso_n = 0;
    OUT::ph_mediumNoSIEIENoNeuIso_n = 0;
    OUT::ph_mediumNoSIEIENoPhoIso_n = 0;
    OUT::ph_mediumNoChIsoNoPhoIso_n = 0;
    OUT::ph_mediumNoChIsoNoNeuIso_n = 0;
    OUT::ph_mediumNoPhoIsoNoNeuIso_n = 0;

    OUT::ph_mediumNoSIEIEPassPSV_n = 0;
    OUT::ph_mediumNoChIsoPassPSV_n = 0;
    OUT::ph_mediumNoNeuIsoPassPSV_n = 0;
    OUT::ph_mediumNoPhoIsoPassPSV_n = 0;

    OUT::ph_mediumNoSIEIEFailPSV_n = 0;
    OUT::ph_mediumNoChIsoFailPSV_n = 0;
    OUT::ph_mediumNoNeuIsoFailPSV_n = 0;
    OUT::ph_mediumNoPhoIsoFailPSV_n = 0;

    OUT::ph_mediumNoSIEIEPassCSEV_n = 0;
    OUT::ph_mediumNoChIsoPassCSEV_n = 0;
    OUT::ph_mediumNoNeuIsoPassCSEV_n = 0;
    OUT::ph_mediumNoPhoIsoPassCSEV_n = 0;

    OUT::ph_mediumNoSIEIEFailCSEV_n = 0;
    OUT::ph_mediumNoChIsoFailCSEV_n = 0;
    OUT::ph_mediumNoNeuIsoFailCSEV_n = 0;
    OUT::ph_mediumNoPhoIsoFailCSEV_n = 0;

    OUT::ptSorted_ph_loose_idx->clear();
    OUT::ptSorted_ph_medium_idx->clear();
    OUT::ptSorted_ph_tight_idx->clear();

    OUT::ptSorted_ph_mediumPassPSV_idx->clear();
    OUT::ptSorted_ph_mediumFailPSV_idx->clear();
    OUT::ptSorted_ph_mediumPassCSEV_idx->clear();
    OUT::ptSorted_ph_mediumFailCSEV_idx->clear();

    OUT::ptSorted_ph_mediumPassEleOlap_idx->clear();
    OUT::ptSorted_ph_mediumPassEleOlapPassCSEV_idx->clear();
    OUT::ptSorted_ph_mediumPassEleOlapFailCSEV_idx->clear();

    OUT::ptSorted_ph_mediumNoSIEIE_idx->clear();
    OUT::ptSorted_ph_mediumNoChIso_idx->clear();
    OUT::ptSorted_ph_mediumNoNeuIso_idx->clear();
    OUT::ptSorted_ph_mediumNoPhoIso_idx->clear();

    OUT::ptSorted_ph_mediumNoSIEIENoChIso_idx->clear();
    OUT::ptSorted_ph_mediumNoSIEIENoNeuIso_idx->clear();
    OUT::ptSorted_ph_mediumNoSIEIENoPhoIso_idx->clear();
    OUT::ptSorted_ph_mediumNoChIsoNoPhoIso_idx->clear();
    OUT::ptSorted_ph_mediumNoChIsoNoNeuIso_idx->clear();
    OUT::ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx->clear();

    OUT::ptSorted_ph_mediumNoSIEIEPassPSV_idx->clear();
    OUT::ptSorted_ph_mediumNoChIsoPassPSV_idx->clear();
    OUT::ptSorted_ph_mediumNoNeuIsoPassPSV_idx->clear();
    OUT::ptSorted_ph_mediumNoPhoIsoPassPSV_idx->clear();

    OUT::ptSorted_ph_mediumNoSIEIEFailPSV_idx->clear();
    OUT::ptSorted_ph_mediumNoChIsoFailPSV_idx->clear();
    OUT::ptSorted_ph_mediumNoNeuIsoFailPSV_idx->clear();
    OUT::ptSorted_ph_mediumNoPhoIsoFailPSV_idx->clear();

    OUT::ptSorted_ph_mediumNoSIEIEPassCSEV_idx->clear();
    OUT::ptSorted_ph_mediumNoChIsoPassCSEV_idx->clear();
    OUT::ptSorted_ph_mediumNoNeuIsoPassCSEV_idx->clear();
    OUT::ptSorted_ph_mediumNoPhoIsoPassCSEV_idx->clear();

    OUT::ptSorted_ph_mediumNoSIEIEFailCSEV_idx->clear();
    OUT::ptSorted_ph_mediumNoChIsoFailCSEV_idx->clear();
    OUT::ptSorted_ph_mediumNoNeuIsoFailCSEV_idx->clear();
    OUT::ptSorted_ph_mediumNoPhoIsoFailCSEV_idx->clear();


    std::vector<std::pair<float, int> > sorted_ph_loose;
    std::vector<std::pair<float, int> > sorted_ph_medium;
    std::vector<std::pair<float, int> > sorted_ph_tight;

    std::vector<std::pair<float, int> > sorted_ph_mediumPassPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumFailPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumPassCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumFailCSEV;

    std::vector<std::pair<float, int> > sorted_ph_mediumPassEleOlap;
    std::vector<std::pair<float, int> > sorted_ph_mediumPassEleOlapPassCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumPassEleOlapFailCSEV;

    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIE;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoNeuIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoPhoIso;

    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIENoChIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIENoNeuIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIENoPhoIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIsoNoPhoIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIsoNoNeuIso;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoPhoIsoNoNeuIso;

    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIEPassPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIsoPassPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoNeuIsoPassPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoPhoIsoPassPSV;

    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIEFailPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIsoFailPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoNeuIsoFailPSV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoPhoIsoFailPSV;

    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIEPassCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIsoPassCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoNeuIsoPassCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoPhoIsoPassCSEV;

    std::vector<std::pair<float, int> > sorted_ph_mediumNoSIEIEFailCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoChIsoFailCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoNeuIsoFailCSEV;
    std::vector<std::pair<float, int> > sorted_ph_mediumNoPhoIsoFailCSEV;


    std::vector<TLorentzVector> photons;
    std::vector<std::pair<float, int> > sorted_photons;
    for( int idx = 0; idx < OUT::ph_n; ++idx ) {
        TLorentzVector phot;
        phot.SetPtEtaPhiM(  OUT::ph_pt->at(idx), 
                            OUT::ph_eta->at(idx),
                            OUT::ph_phi->at(idx),
                            0.0
                        );
        photons.push_back(phot);

        std::pair<float, int> sort_pair = std::make_pair( phot.Pt(), idx );

        sorted_photons.push_back( sort_pair );

        bool passHOverEMedium = OUT::ph_passHOverEMedium->at(idx);
        bool passSIEIEMedium  = OUT::ph_passSIEIEMedium->at(idx);
        bool passChIsoMedium  = OUT::ph_passChIsoCorrMedium->at(idx);
        bool passNeuIsoMedium = OUT::ph_passNeuIsoCorrMedium->at(idx);
        bool passPhoIsoMedium = OUT::ph_passPhoIsoCorrMedium->at(idx);
        bool passEleOlap      = (OUT::ph_min_el_dr->at(idx) > 0.4);
        bool passLoose        = OUT::ph_passLoose->at(idx);
        bool passMedium        = OUT::ph_passMedium->at(idx);
        bool passTight        = OUT::ph_passTight->at(idx);

        bool eleVeto          = OUT::ph_passEleVeto->at(idx);
        bool hasPixSeed       = OUT::ph_hasPixSeed->at(idx);

        if( passHOverEMedium ) {

            if( passNeuIsoMedium && passPhoIsoMedium ) {
                sorted_ph_mediumNoSIEIENoChIso.push_back( sort_pair );
            }
            if( passChIsoMedium && passPhoIsoMedium ) {
                sorted_ph_mediumNoSIEIENoNeuIso.push_back( sort_pair );
            }
            if( passChIsoMedium && passNeuIsoMedium ) {
                sorted_ph_mediumNoSIEIENoPhoIso.push_back( sort_pair );
            }
            if( passSIEIEMedium && passNeuIsoMedium ) {
                sorted_ph_mediumNoChIsoNoPhoIso.push_back( sort_pair );
            }
            if( passSIEIEMedium && passPhoIsoMedium ) {
                sorted_ph_mediumNoChIsoNoNeuIso.push_back( sort_pair );
            }
            if( passSIEIEMedium && passChIsoMedium ) {
                sorted_ph_mediumNoPhoIsoNoNeuIso.push_back( sort_pair );
            }
            
            if( passChIsoMedium && 
                passNeuIsoMedium &&  
                passPhoIsoMedium    ) {
                sorted_ph_mediumNoSIEIE.push_back( sort_pair );
                if( eleVeto ) {
                    sorted_ph_mediumNoSIEIEPassCSEV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoSIEIEFailCSEV.push_back(sort_pair);
                }
                if( hasPixSeed ) {
                    sorted_ph_mediumNoSIEIEFailPSV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoSIEIEPassPSV.push_back(sort_pair);
                }
            }
            if( passSIEIEMedium && 
                passNeuIsoMedium &&  
                passPhoIsoMedium    ) {
                sorted_ph_mediumNoChIso.push_back( sort_pair );
                if( eleVeto ) {
                    sorted_ph_mediumNoChIsoPassCSEV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoChIsoFailCSEV.push_back(sort_pair);
                }
                if( hasPixSeed ) {
                    sorted_ph_mediumNoChIsoFailPSV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoChIsoPassPSV.push_back(sort_pair);
                }
            }
            if( passSIEIEMedium && 
                passChIsoMedium &&  
                passPhoIsoMedium    ) {
                sorted_ph_mediumNoNeuIso.push_back( sort_pair );
                if( eleVeto ) {
                    sorted_ph_mediumNoNeuIsoPassCSEV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoNeuIsoFailCSEV.push_back(sort_pair);
                }
                if( hasPixSeed ) {
                    sorted_ph_mediumNoNeuIsoFailPSV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoNeuIsoPassPSV.push_back(sort_pair);
                }
            }
            if( passSIEIEMedium && 
                passChIsoMedium &&  
                passNeuIsoMedium    ) {
                sorted_ph_mediumNoPhoIso.push_back( sort_pair );
                if( eleVeto ) {
                    sorted_ph_mediumNoPhoIsoPassCSEV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoPhoIsoFailCSEV.push_back(sort_pair);
                }
                if( hasPixSeed ) {
                    sorted_ph_mediumNoPhoIsoFailPSV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumNoPhoIsoPassPSV.push_back(sort_pair);
                }
            }
            if( passSIEIEMedium && 
                passChIsoMedium &&  
                passNeuIsoMedium &&
                passPhoIsoMedium     ) {
                sorted_ph_medium.push_back( sort_pair );
                if( eleVeto ) {
                    sorted_ph_mediumPassCSEV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumFailCSEV.push_back(sort_pair);
                }
                if( hasPixSeed ) {
                    sorted_ph_mediumFailPSV.push_back(sort_pair);
                }
                else {
                    sorted_ph_mediumPassPSV.push_back(sort_pair);
                }
            }
        }
        if( passLoose )  {
            sorted_ph_loose.push_back( sort_pair );
        }
        if( passMedium )  {
            if( passEleOlap ) {
                sorted_ph_mediumPassEleOlap.push_back( sort_pair );
                if( eleVeto ) {
                    sorted_ph_mediumPassEleOlapPassCSEV.push_back( sort_pair );
                }
                else {
                    sorted_ph_mediumPassEleOlapFailCSEV.push_back( sort_pair );
                }
            }
        }
        if( passTight )  {
            sorted_ph_tight.push_back( sort_pair );
        }
    }

    std::sort(sorted_ph_loose.rbegin(),sorted_ph_loose.rend()) ;
    std::sort(sorted_ph_medium.rbegin(),sorted_ph_medium.rend()) ;
    std::sort(sorted_ph_tight.rbegin(),sorted_ph_tight.rend()) ;

    std::sort(sorted_ph_mediumPassPSV.rbegin(),sorted_ph_mediumPassPSV.rend()) ;
    std::sort(sorted_ph_mediumFailPSV.rbegin(),sorted_ph_mediumFailPSV.rend()) ;
    std::sort(sorted_ph_mediumPassCSEV.rbegin(),sorted_ph_mediumPassCSEV.rend()) ;
    std::sort(sorted_ph_mediumFailCSEV.rbegin(),sorted_ph_mediumFailCSEV.rend()) ;

    std::sort(sorted_ph_mediumPassEleOlap.rbegin(),sorted_ph_mediumPassEleOlap.rend());
    std::sort(sorted_ph_mediumPassEleOlapPassCSEV.rbegin(),sorted_ph_mediumPassEleOlapPassCSEV.rend());
    std::sort(sorted_ph_mediumPassEleOlapFailCSEV.rbegin(),sorted_ph_mediumPassEleOlapFailCSEV.rend());


    std::sort(sorted_ph_mediumNoSIEIE.rbegin(),sorted_ph_mediumNoSIEIE.rend()) ;
    std::sort(sorted_ph_mediumNoChIso.rbegin(),sorted_ph_mediumNoChIso.rend()) ;
    std::sort(sorted_ph_mediumNoNeuIso.rbegin(),sorted_ph_mediumNoNeuIso.rend()) ;
    std::sort(sorted_ph_mediumNoPhoIso.rbegin(),sorted_ph_mediumNoPhoIso.rend()) ;

    std::sort(sorted_ph_mediumNoSIEIENoChIso.rbegin(),sorted_ph_mediumNoSIEIENoChIso.rend());
    std::sort(sorted_ph_mediumNoSIEIENoNeuIso.rbegin(),sorted_ph_mediumNoSIEIENoNeuIso.rend());
    std::sort(sorted_ph_mediumNoSIEIENoPhoIso.rbegin(),sorted_ph_mediumNoSIEIENoPhoIso.rend());
    std::sort(sorted_ph_mediumNoChIsoNoPhoIso.rbegin(),sorted_ph_mediumNoChIsoNoPhoIso.rend());
    std::sort(sorted_ph_mediumNoChIsoNoNeuIso.rbegin(),sorted_ph_mediumNoChIsoNoNeuIso.rend());
    std::sort(sorted_ph_mediumNoPhoIsoNoNeuIso.rbegin(),sorted_ph_mediumNoPhoIsoNoNeuIso.rend());

    std::sort(sorted_ph_mediumNoSIEIEPassPSV.rbegin(),sorted_ph_mediumNoSIEIEPassPSV.rend()) ;
    std::sort(sorted_ph_mediumNoChIsoPassPSV.rbegin(),sorted_ph_mediumNoChIsoPassPSV.rend()) ;
    std::sort(sorted_ph_mediumNoNeuIsoPassPSV.rbegin(),sorted_ph_mediumNoNeuIsoPassPSV.rend()) ;
    std::sort(sorted_ph_mediumNoPhoIsoPassPSV.rbegin(),sorted_ph_mediumNoPhoIsoPassPSV.rend()) ;

    std::sort(sorted_ph_mediumNoSIEIEFailPSV.rbegin(),sorted_ph_mediumNoSIEIEFailPSV.rend()) ;
    std::sort(sorted_ph_mediumNoChIsoFailPSV.rbegin(),sorted_ph_mediumNoChIsoFailPSV.rend()) ;
    std::sort(sorted_ph_mediumNoNeuIsoFailPSV.rbegin(),sorted_ph_mediumNoNeuIsoFailPSV.rend()) ;
    std::sort(sorted_ph_mediumNoPhoIsoFailPSV.rbegin(),sorted_ph_mediumNoPhoIsoFailPSV.rend()) ;

    std::sort(sorted_ph_mediumNoSIEIEPassCSEV.rbegin(),sorted_ph_mediumNoSIEIEPassCSEV.rend()) ;
    std::sort(sorted_ph_mediumNoChIsoPassCSEV.rbegin(),sorted_ph_mediumNoChIsoPassCSEV.rend()) ;
    std::sort(sorted_ph_mediumNoNeuIsoPassCSEV.rbegin(),sorted_ph_mediumNoNeuIsoPassCSEV.rend()) ;
    std::sort(sorted_ph_mediumNoPhoIsoPassCSEV.rbegin(),sorted_ph_mediumNoPhoIsoPassCSEV.rend()) ;

    std::sort(sorted_ph_mediumNoSIEIEFailCSEV.rbegin(),sorted_ph_mediumNoSIEIEFailCSEV.rend()) ;
    std::sort(sorted_ph_mediumNoChIsoFailCSEV.rbegin(),sorted_ph_mediumNoChIsoFailCSEV.rend()) ;
    std::sort(sorted_ph_mediumNoNeuIsoFailCSEV.rbegin(),sorted_ph_mediumNoNeuIsoFailCSEV.rend()) ;
    std::sort(sorted_ph_mediumNoPhoIsoFailCSEV.rbegin(),sorted_ph_mediumNoPhoIsoFailCSEV.rend()) ;


    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_loose.begin(); itr != sorted_ph_loose.end(); ++itr ) {
        OUT::ptSorted_ph_loose_idx->push_back( itr->second );
        OUT::ph_loose_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_medium.begin(); itr != sorted_ph_medium.end(); ++itr ) {
        OUT::ptSorted_ph_medium_idx->push_back( itr->second );
        OUT::ph_medium_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_tight.begin(); itr != sorted_ph_tight.end(); ++itr ) {
        OUT::ptSorted_ph_tight_idx->push_back( itr->second );
        OUT::ph_tight_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumPassPSV.begin(); itr != sorted_ph_mediumPassPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumPassPSV_idx->push_back( itr->second );
        OUT::ph_mediumPassPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumFailPSV.begin(); itr != sorted_ph_mediumFailPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumFailPSV_idx->push_back( itr->second );
        OUT::ph_mediumFailPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumPassCSEV.begin(); itr != sorted_ph_mediumPassCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumPassCSEV_idx->push_back( itr->second );
        OUT::ph_mediumPassCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumFailCSEV.begin(); itr != sorted_ph_mediumFailCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumFailCSEV_idx->push_back( itr->second );
        OUT::ph_mediumFailCSEV_n++;
    }

    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumPassEleOlap.begin(); itr != sorted_ph_mediumPassEleOlap.end(); ++itr ) {
        OUT::ptSorted_ph_mediumPassEleOlap_idx->push_back( itr->second );
        OUT::ph_mediumPassEleOlap_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumPassEleOlapPassCSEV.begin(); itr != sorted_ph_mediumPassEleOlapPassCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumPassEleOlapPassCSEV_idx->push_back( itr->second );
        OUT::ph_mediumPassEleOlapPassCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumPassEleOlapFailCSEV.begin(); itr != sorted_ph_mediumPassEleOlapFailCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumPassEleOlapFailCSEV_idx->push_back( itr->second );
        OUT::ph_mediumPassEleOlapFailCSEV_n++;
    }


    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIE.begin(); itr != sorted_ph_mediumNoSIEIE.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIE_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIE_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIso.begin(); itr != sorted_ph_mediumNoChIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIso_idx->push_back( itr->second );
        OUT::ph_mediumNoChIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoNeuIso.begin(); itr != sorted_ph_mediumNoNeuIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoNeuIso_idx->push_back( itr->second );
        OUT::ph_mediumNoNeuIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoPhoIso.begin(); itr != sorted_ph_mediumNoPhoIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoPhoIso_idx->push_back( itr->second );
        OUT::ph_mediumNoPhoIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIENoChIso.begin(); itr != sorted_ph_mediumNoSIEIENoChIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIENoChIso_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIENoChIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIENoNeuIso.begin(); itr != sorted_ph_mediumNoSIEIENoNeuIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIENoNeuIso_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIENoNeuIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIENoPhoIso.begin(); itr != sorted_ph_mediumNoSIEIENoPhoIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIENoPhoIso_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIENoPhoIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIsoNoPhoIso.begin(); itr != sorted_ph_mediumNoChIsoNoPhoIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIsoNoPhoIso_idx->push_back( itr->second );
        OUT::ph_mediumNoChIsoNoPhoIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIsoNoNeuIso.begin(); itr != sorted_ph_mediumNoChIsoNoNeuIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIsoNoNeuIso_idx->push_back( itr->second );
        OUT::ph_mediumNoChIsoNoNeuIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoPhoIsoNoNeuIso.begin(); itr != sorted_ph_mediumNoPhoIsoNoNeuIso.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx->push_back( itr->second );
        OUT::ph_mediumNoPhoIsoNoNeuIso_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIEPassPSV.begin(); itr != sorted_ph_mediumNoSIEIEPassPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIEPassPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIEPassPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIsoPassPSV.begin(); itr != sorted_ph_mediumNoChIsoPassPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIsoPassPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoChIsoPassPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoNeuIsoPassPSV.begin(); itr != sorted_ph_mediumNoNeuIsoPassPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoNeuIsoPassPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoNeuIsoPassPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoPhoIsoPassPSV.begin(); itr != sorted_ph_mediumNoPhoIsoPassPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoPhoIsoPassPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoPhoIsoPassPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIEFailPSV.begin(); itr != sorted_ph_mediumNoSIEIEFailPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIEFailPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIEFailPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIsoFailPSV.begin(); itr != sorted_ph_mediumNoChIsoFailPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIsoFailPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoChIsoFailPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoNeuIsoFailPSV.begin(); itr != sorted_ph_mediumNoNeuIsoFailPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoNeuIsoFailPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoNeuIsoFailPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoPhoIsoFailPSV.begin(); itr != sorted_ph_mediumNoPhoIsoFailPSV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoPhoIsoFailPSV_idx->push_back( itr->second );
        OUT::ph_mediumNoPhoIsoFailPSV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIEPassCSEV.begin(); itr != sorted_ph_mediumNoSIEIEPassCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIEPassCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIEPassCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIsoPassCSEV.begin(); itr != sorted_ph_mediumNoChIsoPassCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIsoPassCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoChIsoPassCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoNeuIsoPassCSEV.begin(); itr != sorted_ph_mediumNoNeuIsoPassCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoNeuIsoPassCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoNeuIsoPassCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoPhoIsoPassCSEV.begin(); itr != sorted_ph_mediumNoPhoIsoPassCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoPhoIsoPassCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoPhoIsoPassCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoSIEIEFailCSEV.begin(); itr != sorted_ph_mediumNoSIEIEFailCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoSIEIEFailCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoSIEIEFailCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoChIsoFailCSEV.begin(); itr != sorted_ph_mediumNoChIsoFailCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoChIsoFailCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoChIsoFailCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoNeuIsoFailCSEV.begin(); itr != sorted_ph_mediumNoNeuIsoFailCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoNeuIsoFailCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoNeuIsoFailCSEV_n++;
    }
    for( std::vector<std::pair<float, int> >::const_iterator itr = sorted_ph_mediumNoPhoIsoFailCSEV.begin(); itr != sorted_ph_mediumNoPhoIsoFailCSEV.end(); ++itr ) {
        OUT::ptSorted_ph_mediumNoPhoIsoFailCSEV_idx->push_back( itr->second );
        OUT::ph_mediumNoPhoIsoFailCSEV_n++;
    }


}

void RunModule::BuildTruth( ModuleConfig & config ) const {

#ifdef EXISTS_gen_n

    OUT::trueph_n = 0;
    OUT::trueph_pt->clear();
    OUT::trueph_eta->clear();
    OUT::trueph_phi->clear();
    OUT::trueph_motherPID->clear();
    OUT::trueph_status->clear();

    OUT::truelep_n = 0;
    OUT::truelep_pt->clear();
    OUT::truelep_eta->clear();
    OUT::truelep_phi->clear();
    OUT::truelep_e->clear();
    OUT::truelep_motherPID->clear();
    OUT::truelep_status->clear();
    OUT::truelep_PID->clear();

    OUT::truenu_n = 0;

    OUT::truelepnu_m = 0.0;
    OUT::truelepnuph_m = 0.0;
    OUT::truelepph_dr = 0.0;
    OUT::truemt_lep_met_ph = 0;
    OUT::truemt_res = 0;
    OUT::truemt_res_l23 = 0;
    OUT::truemt_res_lO = 0;

    OUT::trueht = 0;

    OUT::isWMuDecay = false;
    OUT::isWElDecay = false;
    OUT::isWTauDecay = false;
    OUT::isWTauElDecay = false;
    OUT::isWTauMuDecay = false;
    OUT::WIDStep = 0;

    OUT::ph_hasTruthMatchPh->clear();
    OUT::ph_truthMatchPh_pt->clear();
    OUT::ph_truthMatchPh_dr->clear();

    OUT::el_hasTruthMatchEl->clear();
    OUT::el_truthMatchEl_pt->clear();
    OUT::el_truthMatchEl_dr->clear();

    OUT::mu_hasTruthMatchMu->clear();
    OUT::mu_truthMatchMu_pt->clear();
    OUT::mu_truthMatchMu_dr->clear();

    std::vector<TLorentzVector> trueleps;
    std::vector<TLorentzVector> nulvs;
    std::vector<TLorentzVector> photlvs;
    std::vector<std::pair<float, unsigned> >  sorted_tleps;
    std::vector<std::pair<float, unsigned> >  sorted_tphots;

    TLorentzVector nu_sum;

    for( int gidx = 0; gidx < IN::gen_n; ++gidx ) {

        int id = IN::gen_PID->at(gidx);
        int absid = abs(id);
        int st = IN::gen_status->at(gidx);

        // calculate HT at truth level
        if( ( absid < 6 || id == 21 )  && st == 23 ) {
            OUT::trueht += IN::gen_pt->at(gidx);
        }

        if( absid >= 11 && absid <= 16) {

            // keep all leptons at this point
            // so they can be used later
            sorted_tleps.push_back( std::make_pair( IN::gen_pt->at(gidx), gidx ) );

        }
        // neutrinos
        if( absid == 12 || absid == 14 || absid == 16 ) {

            TLorentzVector nulv;
            nulv.SetPtEtaPhiM( IN::gen_pt->at(gidx),
                               IN::gen_eta->at(gidx),
                               IN::gen_phi->at(gidx),
                               0.0
                             );

            nu_sum += nulv;
            nulvs.push_back( nulv );

            OUT::truenu_n++;
        }
        // photons
        if( absid == 22 ) {

            bool pass_ph_cuts = true;

            if( !config.PassInt( "cut_ph_mother", abs(IN::gen_motherPID->at(gidx) ) ) ) pass_ph_cuts = false;
            if( !config.PassInt( "cut_ph_status", IN::gen_status->at(gidx) ) ) pass_ph_cuts = false;
            if( !config.PassBool( "cut_ph_IsPromptFinalState", IN::gen_isPromptFinalState->at(gidx) ) ) pass_ph_cuts = false;
            if( !config.PassBool( "cut_ph_FromHardProcessFinalState", IN::gen_fromHardProcessFinalState->at(gidx) ) ) pass_ph_cuts = false;


            if( pass_ph_cuts ) {

                sorted_tphots.push_back( std::make_pair( IN::gen_pt->at(gidx), gidx ) );
            }
        }
    }

    // sort leptons so the highest pT leptons come first
    std::sort(sorted_tleps.rbegin(), sorted_tleps.rend());
    // sort photons so the highest pT leptons come first
    std::sort(sorted_tphots.rbegin(), sorted_tphots.rend());

    // now fill outputs using the sorted objects
    for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
            sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

        unsigned lidx = litr->second;

        int pid = IN::gen_PID->at(lidx);
        // should be an electron or muon
        if( !( abs(pid) == 11 || abs(pid) == 13 ) ) continue;

        if( !config.PassInt( "cut_lep_mother", IN::gen_motherPID->at(lidx) ) ) continue;
        if( !config.PassInt( "cut_lep_status", IN::gen_status->at(lidx) ) ) continue;

        OUT::truelep_n++;
        OUT::truelep_pt->push_back( IN::gen_pt->at( lidx ) );
        OUT::truelep_eta->push_back( IN::gen_eta->at( lidx ) );
        OUT::truelep_phi->push_back( IN::gen_phi->at( lidx ) );
        OUT::truelep_e->push_back( IN::gen_e->at( lidx ) );
        OUT::truelep_motherPID->push_back( IN::gen_motherPID->at( lidx ) );
        OUT::truelep_status->push_back( IN::gen_status->at( lidx ) );
        OUT::truelep_PID->push_back( IN::gen_PID->at( lidx ) );

        TLorentzVector leplv;
        leplv.SetPtEtaPhiE( IN::gen_pt->at(lidx),
                            IN::gen_eta->at(lidx),
                            IN::gen_phi->at(lidx),
                            IN::gen_e->at(lidx)
                          );

        trueleps.push_back( leplv );

    }

    for( std::vector<std::pair<float, unsigned> >::const_iterator pitr = 
            sorted_tphots.begin(); pitr != sorted_tphots.end(); ++pitr ) {

        unsigned pidx = pitr->second;

        OUT::trueph_n++;
        OUT::trueph_pt->push_back( IN::gen_pt->at( pidx ) );
        OUT::trueph_eta->push_back( IN::gen_eta->at( pidx ) );
        OUT::trueph_phi->push_back( IN::gen_phi->at( pidx ) );
        OUT::trueph_motherPID->push_back( IN::gen_motherPID->at( pidx ) );
        OUT::trueph_status->push_back( IN::gen_status->at( pidx ) );

        TLorentzVector phlv;
        phlv.SetPtEtaPhiM( IN::gen_pt->at(pidx),
                           IN::gen_eta->at(pidx),
                           IN::gen_phi->at(pidx),
                           0.0
                         );

        photlvs.push_back( phlv );
    }

    //photon truth match
    for( int phidx = 0; phidx < OUT::ph_n; ++phidx ) {

        TLorentzVector phlv;
        phlv.SetPtEtaPhiE( OUT::ph_pt->at(phidx),
                           OUT::ph_eta->at(phidx),
                           OUT::ph_phi->at(phidx),
                           OUT::ph_e->at(phidx)
                         );


        float mindr = 101.;
        float matchPt = 0;
        for( int tphidx = 0; tphidx < OUT::trueph_n; ++tphidx ) {

            TLorentzVector tphlv;
            tphlv.SetPtEtaPhiM( OUT::trueph_pt->at(tphidx),
                                OUT::trueph_eta->at(tphidx),
                                OUT::trueph_phi->at(tphidx),
                                0.0
                              );


            float dr = phlv.DeltaR( tphlv );

            if( dr < mindr ) {
                mindr = dr;
                matchPt = tphlv.Pt();
            }

        }

        OUT::ph_hasTruthMatchPh->push_back( ( mindr < 100 ) );
        OUT::ph_truthMatchPh_pt->push_back( matchPt );
        OUT::ph_truthMatchPh_dr->push_back( mindr );

    }

    // muon truth match
    for( int muidx = 0; muidx < OUT::mu_n; ++muidx ) {

        TLorentzVector mulv;
        mulv.SetPtEtaPhiE( OUT::mu_pt->at(muidx),
                           OUT::mu_eta->at(muidx),
                           OUT::mu_phi->at(muidx),
                           OUT::mu_e->at(muidx)
                         );


        float mindr = 101.;
        float matchPt = 0;
        for( int tlepidx = 0; tlepidx < OUT::truelep_n; ++tlepidx ) {

            if( abs( OUT::truelep_PID->at(tlepidx) ) != 13 ) continue;

            TLorentzVector tleplv;
            tleplv.SetPtEtaPhiE( OUT::truelep_pt->at(tlepidx),
                                 OUT::truelep_eta->at(tlepidx),
                                 OUT::truelep_phi->at(tlepidx),
                                 OUT::truelep_e->at(tlepidx)
                               );

            float dr = mulv.DeltaR( tleplv );

            if( dr < mindr ) {
                mindr = dr;
                matchPt = tleplv.Pt();
            }

        }

        OUT::mu_hasTruthMatchMu->push_back( ( mindr < 100 ) );
        OUT::mu_truthMatchMu_pt->push_back( matchPt );
        OUT::mu_truthMatchMu_dr->push_back( mindr );

    }

    // electron truth match
    for( int elidx = 0; elidx < OUT::el_n; ++elidx ) {

        TLorentzVector ellv;
        ellv.SetPtEtaPhiE( OUT::el_pt->at(elidx),
                           OUT::el_eta->at(elidx),
                           OUT::el_phi->at(elidx),
                           OUT::el_e->at(elidx)
                         );


        float mindr = 101.;
        float matchPt = 0;
        for( int tlepidx = 0; tlepidx < OUT::truelep_n; ++tlepidx ) {

            if( abs( OUT::truelep_PID->at(tlepidx) ) != 11 ) continue;

            TLorentzVector tleplv;
            tleplv.SetPtEtaPhiE( OUT::truelep_pt->at(tlepidx),
                                 OUT::truelep_eta->at(tlepidx),
                                 OUT::truelep_phi->at(tlepidx),
                                 OUT::truelep_e->at(tlepidx)
                               );

            float dr = ellv.DeltaR( tleplv );

            if( dr < mindr ) {
                mindr = dr;
                matchPt = tleplv.Pt();
            }

        }

        OUT::el_hasTruthMatchEl->push_back( ( mindr < 100 ) );
        OUT::el_truthMatchEl_pt->push_back( matchPt );
        OUT::el_truthMatchEl_dr->push_back( mindr );

    }

    if( trueleps.size() > 0 && OUT::trueph_n>0) {
   
        TLorentzVector lep_trans;
        TLorentzVector ph_trans;
        lep_trans.SetPtEtaPhiM( trueleps[0].Pt(), 0.0, trueleps[0].Phi(), trueleps[0].M() );
        ph_trans.SetPtEtaPhiM( OUT::trueph_pt->at(0), 0.0, OUT::trueph_phi->at(0), 0.0 );

        TLorentzVector metlv;
        metlv.SetPtEtaPhiM( nu_sum.Pt(), 0.0, nu_sum.Phi(), 0.0 );

        OUT::truemt_res = ( lep_trans + ph_trans + metlv ).M();
        //std::cout << "Truelep pt = " << trueleps[0].Pt() << " nu pt = " << nu_sum.Pt() << std::endl;
        OUT::truemt_lep_met_ph = Utils::calc_mt( trueleps[0], nu_sum);

        if( trueleps.size() == 2  ) { 

            TLorentzVector l23_trans;
            TLorentzVector lO_trans;
            for( int lidx = 0; lidx != OUT::truelep_n ; ++lidx ) {

                TLorentzVector leplv = trueleps[lidx];

                if( OUT::truelep_status->at(lidx) == 23 ) {
                    l23_trans.SetPtEtaPhiM( leplv.Pt(), 0.0, leplv.Phi(), leplv.M() );
                }
                else {
                    lO_trans.SetPtEtaPhiM( leplv.Pt(), 0.0, leplv.Phi(), leplv.M() );
                }
            }
            OUT::truemt_res_l23 = ( l23_trans + ph_trans + metlv ).M();
            OUT::truemt_res_lO = ( lO_trans + ph_trans + metlv ).M();
        }
    }


    OUT::truelepnu_m = 0.0;
    OUT::truelepnuph_m = 0.0;
    OUT::truelepph_dr = 0.0;

    if( OUT::truelep_n == 1 ) {
        TLorentzVector wlv( trueleps[0] );

        for( unsigned i = 0 ; i < nulvs.size(); ++i) {
            wlv = wlv + nulvs[i];
        }
        OUT::truelepnu_m = wlv.M();

        if( photlvs.size() > 0 ) {
            OUT::truelepnuph_m = (wlv + photlvs[0]).M();
            OUT::truelepph_dr =  trueleps[0].DeltaR( photlvs[0] );
        }
    }

    // in dermining the W decay mode use the original
    // vector of leptons.  Do not use trueleps which
    // may have additional cuts
    OUT::isWMuDecay = false;
    OUT::isWElDecay = false;
    OUT::isWTauDecay = false;
    OUT::isWTauElDecay = false;
    OUT::isWTauMuDecay = false;
    OUT::WIDStep = 0;
    bool found_w_mother = false;
    for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
            sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

        if( abs(IN::gen_motherPID->at(litr->second)) != 24 ) continue;
        if( IN::gen_status->at(litr->second) != 3 ) continue;

        found_w_mother = true;

        int lep_id = IN::gen_PID->at(litr->second);

        if( abs(lep_id) == 11 || abs(lep_id) == 12 ) {
            OUT::isWElDecay=true;
        }
        if( abs(lep_id) == 13 || abs(lep_id) == 14 ) {
            OUT::isWMuDecay=true;
        }
        if( abs(lep_id) == 15 || abs(lep_id) == 16 ) {
            OUT::isWTauDecay=true;
        }
    }

    if( !found_w_mother ) {
        OUT::WIDStep++;
        for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
                sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

            if( abs(IN::gen_motherPID->at(litr->second)) != 24 ) continue;
            if( IN::gen_status->at(litr->second) != 23 ) continue;

            found_w_mother = true;

            int lep_id = IN::gen_PID->at(litr->second);

            if( abs(lep_id) == 11 || abs(lep_id) == 12 ) {
                OUT::isWElDecay=true;
            }
            if( abs(lep_id) == 13 || abs(lep_id) == 14 ) {
                OUT::isWMuDecay=true;
            }
            if( abs(lep_id) == 15 || abs(lep_id) == 16 ) {
                OUT::isWTauDecay=true;
            }
        }
    }

    if( !found_w_mother ) {
        OUT::WIDStep++;
        for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
                sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

            if( abs(IN::gen_motherPID->at(litr->second)) != 24 ) continue;
            if( IN::gen_status->at(litr->second) != 2 ) continue;

            found_w_mother = true;

            int lep_id = IN::gen_PID->at(litr->second);

            if( abs(lep_id) == 11 || abs(lep_id) == 12 ) {
                OUT::isWElDecay=true;
            }
            if( abs(lep_id) == 13 || abs(lep_id) == 14 ) {
                OUT::isWMuDecay=true;
            }
            if( abs(lep_id) == 15 || abs(lep_id) == 16 ) {
                OUT::isWTauDecay=true;
            }
        }
    }

    if( !found_w_mother ) {
        OUT::WIDStep++;
        for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
                sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

            if( abs(IN::gen_motherPID->at(litr->second)) != 24 ) continue;
            if( IN::gen_status->at(litr->second) != 1 ) continue;

            found_w_mother = true;

            int lep_id = IN::gen_PID->at(litr->second);
            if( abs(lep_id) == 11 || abs(lep_id) == 12 ) {
                OUT::isWElDecay=true;
            }
            if( abs(lep_id) == 13 || abs(lep_id) == 14 ) {
                OUT::isWMuDecay=true;
            }
            if( abs(lep_id) == 15 || abs(lep_id) == 16 ) {
                OUT::isWTauDecay=true;
            }
        }
    }
    if( !found_w_mother ) {
        OUT::WIDStep++;
        int n_el    = 0;
        int n_mu    = 0;
        int n_tau   = 0;
        int n_elnu  = 0;
        int n_munu  = 0;
        int n_taunu = 0;
        for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
                sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

            if( IN::gen_status->at(litr->second) != 1 ) continue;

            int momId = abs(IN::gen_motherPID->at(litr->second));
            int Id    = abs(IN::gen_PID      ->at(litr->second));

            if( ( Id == momId) || ( momId == 15 ) ) {

                if( Id == 11 ) n_el++;
                if( Id == 12 ) n_elnu++;
                if( Id == 13 ) n_mu++;
                if( Id == 14 ) n_munu++;
                if( Id == 15 ) n_tau++;
                if( Id == 16 ) n_taunu++;

            }
        }

        // if anything tau-like is present then its a tau decay
        if( n_taunu > 0 || n_tau > 0 ) {
            found_w_mother = true;
            OUT::isWTauDecay = true;
        }
        // otherwise the final state objects should be consistent with a W decay
        else if( n_munu > 0 && n_mu > 0 ) {
            OUT::isWMuDecay = true;
        }
        // otherwise the final state objects should be consistent with a W decay
        else if( n_elnu > 0 && n_el > 0 ) {
            OUT::isWElDecay = true;
        }
    }

    if( OUT::isWTauDecay ) {
        std::vector<unsigned> found_leptons;
        for( std::vector<std::pair<float, unsigned> >::const_iterator litr = 
                sorted_tleps.begin(); litr != sorted_tleps.end(); ++litr ) {

            if( IN::gen_status->at(litr->second) != 1 ) continue;

            int momId = abs(IN::gen_motherPID->at(litr->second));
            int Id    = abs(IN::gen_PID      ->at(litr->second));

            if( fabs(momId) == 15 && (fabs( Id ) == 11 || fabs(Id) == 13) ) {

                found_leptons.push_back( litr->second );
            }
        }

        if( found_leptons.size() > 0 ) {
            int lep_id = IN::gen_PID->at(found_leptons[0]);
            if( fabs(lep_id) == 11 ) {
                OUT::isWTauElDecay=true;
            }
            if( fabs(lep_id) == 13 ) {
                OUT::isWTauMuDecay =true;
            }
        }
    }
#endif

}

bool RunModule::FilterDataQuality( ModuleConfig & config) const {

    bool pass_quality = false;
    OUT::PassQuality = -1;

    bool isData = IN::isData;
    int run = IN::runNumber;
    int ls  = IN::lumiSection;

    if( _quality_map.size() > 0 && isData ) {

        std::map<int, std::vector<int> >::const_iterator mitr = _quality_map.find( run );
        if( mitr != _quality_map.end() ) {

            std::vector<int>::const_iterator vitr = std::find( mitr->second.begin(), mitr->second.end(), ls );

            if( vitr != mitr->second.end() ) {
                pass_quality = true;
                OUT::PassQuality = 1;
            }
            else {
                pass_quality = false;
                OUT::PassQuality = 0;
            }
        }
        else {
            pass_quality = false;
            OUT::PassQuality = 0;
        }
    }
    else {
        pass_quality = true;
    }

    return pass_quality;

}


bool RunModule::FilterBlind( ModuleConfig & config ) const {

    bool keep_event = true;

    bool pass_blind = true;
    if( OUT::ph_n > 0 ) {
        if( !config.PassFloat( "cut_ph_pt_lead", OUT::ph_pt->at(0)) ) pass_blind=false;
    }
    if( !config.PassFloat( "cut_mt_lep_met_ph", OUT::mt_lep_met_ph) ) pass_blind=false;
    if( !config.PassFloat( "cut_mt_res", OUT::mt_res) ) pass_blind=false;

    if( OUT::jet_n > 1 ) {
        if( !config.PassFloat( "cut_abs_dijet_m_from_z", fabs(OUT::leaddijet_m-91.2)) ) pass_blind=false;
    }

    if( !pass_blind ) {
        OUT::isBlinded=true;
        if( OUT::isData ) keep_event=false;
    }
    else {
        OUT::isBlinded=false;
    }

    return keep_event;

}


bool RunModule::get_constriained_nu_pz( const TLorentzVector lepton, TLorentzVector &metlv ) const {

    float solved_pz = -1;

    bool desc_pos = calc_constrained_nu_momentum( lepton, metlv, solved_pz );
    if( desc_pos ) {
        metlv.SetXYZM( metlv.Px(), metlv.Py(), solved_pz, 0.0 );
    }
    else {
        //std::cout << "DISCRIMINANT IS NEGATIVE" << std::endl;
        // require the discriminant to be zero
        // solve a second quadratic equation to 
        // rescale MET so that there is a solution

        float alpha = ( lepton.Px()*metlv.Px() + lepton.Py()*metlv.Py() )/ metlv.Pt();
        float delta = ( _m_w*_m_w - lepton.M()*lepton.M() );

        float Aval = 4*lepton.Pz()*lepton.Pz() - 4*lepton.E()*lepton.E() +4*alpha*alpha;
        float Bval = 4*alpha*delta;
        float Cval = delta*delta;

        float solution1=-1;
        float solution2=-1;

        bool success2 = solve_quadratic( Aval, Bval, Cval, solution1, solution2 );

        if( !success2 ) {
            std::cout << "SECOND FAILURE" << std::endl;
        }

        float scale1 = solution1/metlv.Pt();
        float scale2 = solution2/metlv.Pt();

        TLorentzVector metlv_sol1;
        TLorentzVector metlv_sol2;
        metlv_sol1.SetPtEtaPhiM( OUT::met_pt*scale1, 0.0, OUT::met_phi, 0.0 );
        metlv_sol2.SetPtEtaPhiM( OUT::met_pt*scale2, 0.0, OUT::met_phi, 0.0 );

        float pz_sol1 = -1;
        float pz_sol2 = -1;
        bool success_sol1 = calc_constrained_nu_momentum( lepton, metlv_sol1, pz_sol1 );
        bool success_sol2 = calc_constrained_nu_momentum( lepton, metlv_sol2, pz_sol2 );

        if( !success_sol1 ) {
            //std::cout << "FAILURE SOLUTION 1" << std::endl;
            metlv.SetPtEtaPhiM(-1, 0, 0, 0);
            return false;
        }

        if( !success_sol2 ) {
            //std::cout << "FAILURE SOLUTION 2" << std::endl;
            metlv.SetPtEtaPhiM(-1, 0, 0, 0);
            return false;
        }

        TVector3 solved_met3v_sol1;
        TVector3 solved_met3v_sol2;
        solved_met3v_sol1.SetXYZ(metlv_sol1.Px(), metlv_sol1.Py(), pz_sol1);
        solved_met3v_sol2.SetXYZ(metlv_sol2.Px(), metlv_sol2.Py(), pz_sol2);
        TLorentzVector solved_metlv_sol1;
        TLorentzVector solved_metlv_sol2;
        solved_metlv_sol1.SetVectM( solved_met3v_sol1 , 0.0 );
        solved_metlv_sol2.SetVectM( solved_met3v_sol2 , 0.0 );

        float wmass_sol1 = ( lepton + solved_metlv_sol1 ).M();
        float wmass_sol2 = ( lepton + solved_metlv_sol2 ).M();

        if( fabs( wmass_sol1 - _m_w ) < fabs( wmass_sol2 - _m_w ) ) {
            solved_pz = pz_sol1;
            metlv = metlv_sol1;
        }
        else {
            solved_pz = pz_sol2;
            metlv = metlv_sol2;
        }
        
    }
    return desc_pos;
}

bool RunModule::calc_constrained_nu_momentum( const TLorentzVector lepton, const TLorentzVector met, float & result ) const {

   float little_a = _m_w*_m_w - lepton.M()*lepton.M() + 2*( lepton.Px()*met.Px() + lepton.Py()*met.Py() );

   float Aval = ( 4*lepton.E()*lepton.E() ) - ( 4*lepton.Pz()*lepton.Pz() );
   float Bval = -4 * little_a * lepton.Pz();

   float Cval = 4*lepton.E()*lepton.E()*met.Pt()*met.Pt() - little_a*little_a;

   float solution1=-1;
   float solution2=-1;
   bool success = solve_quadratic( Aval, Bval, Cval, solution1, solution2 );

   if ( success ) {
      if( fabs(solution1 - lepton.Pz() ) < fabs( solution2 - lepton.Pz() ) ) {
          result = solution1;
      }
      else {
          result = solution2;
      }
   }
   return success;
}

bool RunModule::solve_quadratic( float Aval, float Bval, float Cval, float & solution1, float &solution2 ) const {

   float discriminant = Bval*Bval - 4*Aval*Cval;

   //std::cout << "DISCRIMINANT = " << discriminant << std::endl;

   if ( discriminant >= 0 ) {
      solution1 = ( -1*Bval + sqrt( discriminant ) ) / ( 2 * Aval ) ; 
      solution2 = ( -1*Bval - sqrt( discriminant ) ) / ( 2 * Aval ) ; 
      return true;
   }
   else {
       return false;
   }
}
void RunModule::WeightEvent( ModuleConfig & config ) const {

#ifdef EXISTS_EventWeights

    OUT::NLOWeight = 1.0;
    if( !IN::isData && _needs_nlo_weght ) {
        if( IN::EventWeights->at(0) < 0 ) {
            OUT::NLOWeight = -1.0;
        }
    }

#endif 
#ifdef EXISTS_truepu_n
    if( !_puweight_data_hist || !_puweight_sample_hist ) {
        OUT::PUWeight = 1.0;
        return;
    }
    float puval = -1;
    puval = IN::truepu_n;
    //#ifdef EXISTS_nVtx
    //puval = IN::nVtx;
    //#endif

    OUT::PUWeight = calc_pu_weight( puval );

    OUT::PUWeightUP5 = calc_pu_weight( puval, 1.05 );
    OUT::PUWeightUP10 = calc_pu_weight( puval, 1.10 );

    OUT::PUWeightDN5 = calc_pu_weight( puval, 0.95 );
    OUT::PUWeightDN10 = calc_pu_weight( puval, 0.9 );

#endif
    
}

float RunModule::calc_pu_weight( float puval, float mod ) const {

    float tot_data   = _puweight_data_hist->Integral();
    float tot_sample = _puweight_sample_hist->Integral();

    int bin_sample = _puweight_sample_hist->FindBin(puval);
    int bin_data = _puweight_data_hist->FindBin(mod*puval);

    float val_data = _puweight_data_hist->GetBinContent( bin_data );
    float val_sample = _puweight_sample_hist->GetBinContent( bin_sample );

    float num = val_data*mod/tot_data;
    float den = val_sample/tot_sample;

    float weight = num/den;

    //std::cout << "puval = " << puval << " data_val = " << val_data << " sample_val = " << val_sample << " num = " << num << " den = " << den << " weight = " << weight << std::endl;

    if( weight < 0.00000005 ) {

        int bin_min_sample = _puweight_sample_hist->FindBin(puval-2.5);
        int bin_max_sample = _puweight_sample_hist->FindBin(puval+2.5);
        int bin_min_data = _puweight_data_hist->FindBin(puval*mod-2.5);
        int bin_max_data = _puweight_data_hist->FindBin(puval*mod+2.5);

        val_data = _puweight_data_hist->Integral(bin_min_data, bin_max_data);
        val_sample = _puweight_sample_hist->Integral(bin_min_sample, bin_max_sample);

        num = val_data/tot_data;
        den = val_sample/tot_sample;

        weight = num/den;

    }
    return weight;
}



bool RunModule::HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR ) const {
    
    float minDR = 100.0;
    TLorentzVector matchLV;
    return HasTruthMatch( objlv, matchPID, maxDR, minDR, matchLV );

}

bool RunModule::HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR, float &minDR ) const {
    
    TLorentzVector matchLV;
    return HasTruthMatch( objlv, matchPID, maxDR, minDR, matchLV );

}

bool RunModule::HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR, float &minDR, TLorentzVector &matchLV ) const {
    
    int motherPID = 0;
    int parentage = -1;
    return HasTruthMatch( objlv, matchPID, maxDR, minDR, matchLV, motherPID, parentage );

}

bool RunModule::HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR, float & minDR, TLorentzVector & matchLV, int &matchMotherPID, int &matchParentage  ) const {
   
    minDR = 100.0;
    matchLV.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    bool match=false;
#ifdef EXISTS_Gen_pt

    // store commom formats for all
    // information 
    std::vector<int> match_pid;
    std::vector<TLorentzVector> match_tlv;
    std::vector<int> match_momPID;
    std::vector<int> match_parentage;
    BOOST_FOREACH( int mpid, matchPID ) {

        // use photons
        if( mpid == 22 ) {

            for( int gidx = 0; gidx < OUT::trueph_n; ++gidx ) {

                TLorentzVector mclv;
                mclv.SetPtEtaPhiM( OUT::trueph_pt->at(gidx), 
                                   OUT::trueph_eta->at(gidx), 
                                   OUT::trueph_phi->at(gidx), 
                                   0. );

                match_pid.push_back( mpid );
                match_tlv.push_back( mclv );
                match_momPID.push_back( OUT::trueph_motherPID->at(gidx) );
                match_parentage.push_back( 0 );
            }
        }
        else if( mpid==11 || mpid==13 || mpid==-11 || mpid==-13 ) {
            for( int lidx = 0; lidx < OUT::truelep_n; ++lidx ) {

                // only use this lepton type
                if( mpid != OUT::truelep_PID->at(lidx) ) continue;

                TLorentzVector mclv;
                mclv.SetPtEtaPhiE( OUT::truelep_pt->at(lidx), 
                                   OUT::truelep_eta->at(lidx), 
                                   OUT::truelep_phi->at(lidx), 
                                   OUT::truelep_e->at(lidx) );

                match_pid.push_back( mpid );
                match_tlv.push_back( mclv );
                match_momPID.push_back( OUT::truelep_motherPID->at(lidx) );
                match_parentage.push_back( 0 );
            }
        }
        //else {
        //    for( unsigned lidx = 0; lidx < IN::GJetAk04Pt->size(); ++lidx ) {
        //        TLorentzVector mclv;
        //        mclv.SetPtEtaPhiE( IN::GJetAk04Pt->at(lidx), IN::GJetAk04Eta->at(lidx), IN::GJetAk04Phi->at(lidx), IN::GJetAk04E->at(lidx) );

        //        match_pid.push_back( mpid );
        //        match_tlv.push_back( mclv );
        //        // FIX should be filled with the originator quark/gluon
        //        match_momPID.push_back( 1 );
        //        match_parentage.push_back( 0 );
        //    }
        //}


    }

    for( unsigned idx = 0; idx < match_tlv.size(); ++idx ) {

        float dr = match_tlv[idx].DeltaR( objlv );
        //std::cout << "dr = " << dr << std::endl;
        if( dr < maxDR) {
            match = true;
            matchMotherPID = match_momPID[idx];
            matchParentage = match_parentage[idx]; 
        }
        // store the minimum delta R
        if( dr < minDR ) {
            minDR = dr;
            matchLV = match_tlv[idx];
        }
    }

#endif
    return match;

}


RunModule::RunModule() {
    _m_w = 80.385;
    _isData = false;
    _eval_mu_loose    =false;
    _eval_mu_medium   =false;
    _eval_mu_tight    =false;
    _eval_ph_tight    =false;
    _eval_ph_medium   =false;
    _eval_ph_loose    =false;
    _eval_el_tight    =false;
    _eval_el_medium   =false;
    _eval_el_loose    =false;
    _eval_el_veryloose=false;
    _needs_nlo_weght  =false;
}
