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
    std::cout << "Configured " << ana_config.size() << " analysis modules " << std::endl;

    RunModule runmod;
    ana_config.Run(runmod, options);

    std::cout << "^_^ Finished ^_^" << std::endl;


}

void RunModule::initialize( TChain * chain, TTree * outtree, TFile *outfile,
                            const CmdOptions & options, std::vector<ModuleConfig> &configs ) {

    // *************************
    // initialize trees
    // *************************
    InitINTree(chain);
    InitOUTTree( outtree );
    
    // *************************
    // Set defaults for added output variables
    // *************************
    // Examples :

    OUT::m_lep_ph        = 0;
    OUT::m_lep_met_ph        = 0;
    OUT::dphi_lep_ph        = 0;
    OUT::dr_lep_ph        = 0;
    OUT::mt_lep_met      = 0;
    OUT::m_lep_met      = 0;
    OUT::pt_lep_met      = 0;
    OUT::dphi_lep_met      = 0;
    OUT::mt_lep_met_ph   = 0;
    OUT::mt_lep_met_ph_inv   = 0;
    OUT::dphi_ph_w   = 0;
    OUT::pt_lep_met_ph   = 0;
    OUT::RecoWMass       = 0;
    OUT::recoM_lep_nu_ph = 0;
    OUT::recoMet_eta= 0;
    OUT::recoW_pt= 0;
    OUT::recoW_eta= 0;
    OUT::recoW_phi= 0;
    OUT::nu_z_solution_success = 0;
    OUT::m_ll = 0;
    OUT::isBlinded = 0;
    OUT::SampleWeight = 0;

    OUT::ph_hasTruthMatchPh = 0;
    OUT::ph_truthMatchPh_dr = 0;
    OUT::ph_truthMatchPh_pt = 0;

    OUT::mu_hasTruthMatchMu = 0;
    OUT::mu_truthMatchMu_dr = 0;
    OUT::mu_truthMatchMu_pt = 0;

    OUT::el_hasTruthMatchEl = 0;
    OUT::el_truthMatchEl_dr = 0;
    OUT::el_truthMatchEl_pt = 0;

    OUT::ph_medium_n = 0;

    OUT::ph_mediumPassPSV_n = 0;
    OUT::ph_mediumFailPSV_n = 0;
    OUT::ph_mediumPassCSEV_n = 0;
    OUT::ph_mediumFailCSEV_n = 0;

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

    OUT::ptSorted_ph_loose_idx = 0;
    OUT::ptSorted_ph_medium_idx = 0;
    OUT::ptSorted_ph_tight_idx = 0;

    OUT::ptSorted_ph_mediumPassPSV_idx = 0;
    OUT::ptSorted_ph_mediumFailPSV_idx = 0;
    OUT::ptSorted_ph_mediumPassCSEV_idx = 0;
    OUT::ptSorted_ph_mediumFailCSEV_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIE_idx = 0;
    OUT::ptSorted_ph_mediumNoChIso_idx = 0;
    OUT::ptSorted_ph_mediumNoNeuIso_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIso_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIENoChIso_idx = 0;
    OUT::ptSorted_ph_mediumNoSIEIENoNeuIso_idx = 0;
    OUT::ptSorted_ph_mediumNoSIEIENoPhoIso_idx = 0;
    OUT::ptSorted_ph_mediumNoChIsoNoPhoIso_idx = 0;
    OUT::ptSorted_ph_mediumNoChIsoNoNeuIso_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIEPassPSV_idx = 0;
    OUT::ptSorted_ph_mediumNoChIsoPassPSV_idx = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoPassPSV_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoPassPSV_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIEFailPSV_idx = 0;
    OUT::ptSorted_ph_mediumNoChIsoFailPSV_idx = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoFailPSV_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoFailPSV_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIEPassCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoChIsoPassCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoPassCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoPassCSEV_idx = 0;

    OUT::ptSorted_ph_mediumNoSIEIEFailCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoChIsoFailCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoNeuIsoFailCSEV_idx = 0;
    OUT::ptSorted_ph_mediumNoPhoIsoFailCSEV_idx = 0;

    // *************************
    // Declare Branches
    // *************************

    outtree->Branch("el_pt30_n", &OUT::el_pt30_n, "el_pt30_n/I"  );
    outtree->Branch("mu_pt20_n", &OUT::mu_pt20_n, "mu_pt20_n/I"  );
    outtree->Branch("mu_pt30_n", &OUT::mu_pt30_n, "mu_pt30_n/I"  );

    outtree->Branch("jet_CSVLoose_n", &OUT::jet_CSVLoose_n, "jet_CSVLoose_n/I"  );
    outtree->Branch("jet_CSVMedium_n", &OUT::jet_CSVMedium_n, "jet_CSVMedium_n/I"  );
    outtree->Branch("jet_CSVTight_n", &OUT::jet_CSVTight_n, "jet_CSVTight_n/I"  );

    outtree->Branch("m_lep_ph"        , &OUT::m_lep_ph        , "m_lep_ph/F"  );
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

    outtree->Branch("truemt_lep_met_ph" , &OUT::truemt_lep_met_ph, "truemt_lep_met_ph/F");
    outtree->Branch("truemt_res" , &OUT::truemt_res, "truemt_res/F");

    outtree->Branch("ph_hasTruthMatchPh", &OUT::ph_hasTruthMatchPh );
    outtree->Branch("ph_truthMatchPh_dr", &OUT::ph_truthMatchPh_dr);
    outtree->Branch("ph_truthMatchPh_pt", &OUT::ph_truthMatchPh_pt);

    outtree->Branch("mu_hasTruthMatchMu", &OUT::mu_hasTruthMatchMu );
    outtree->Branch("mu_truthMatchMu_dr", &OUT::mu_truthMatchMu_dr);
    outtree->Branch("mu_truthMatchMu_pt", &OUT::mu_truthMatchMu_pt);

    outtree->Branch("el_hasTruthMatchEl", &OUT::el_hasTruthMatchEl );
    outtree->Branch("el_truthMatchEl_dr", &OUT::el_truthMatchEl_dr);
    outtree->Branch("el_truthMatchEl_pt", &OUT::el_truthMatchEl_pt);

    outtree->Branch("mu_hasTrigMatch", &OUT::mu_hasTrigMatch );
    outtree->Branch("mu_trigMatchMinDr", &OUT::mu_trigMatchMinDr);

    outtree->Branch("el_hasTrigMatch", &OUT::el_hasTrigMatch );
    outtree->Branch("el_trigMatchMinDr", &OUT::el_trigMatchMinDr);

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
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "triggerBits" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::vector<std::string> trigger_bit_list = Tokenize( eitr->second, "," );
                for( std::vector<std::string>::const_iterator bitr = trigger_bit_list.begin(); bitr != trigger_bit_list.end(); ++bitr ) {
                    std::stringstream ss( *bitr );
                    int bit;
                    ss >> bit;
                    _muonTrigMatchBits.push_back(  bit );
                }
            }
        }
        if( mod_conf.GetName() == "FilterElectron" ) { 
            std::map<std::string, std::string>::const_iterator eitr = mod_conf.GetInitData().find( "triggerBits" );
            if( eitr != mod_conf.GetInitData().end() ) {
                std::vector<std::string> trigger_bit_list = Tokenize( eitr->second, "," );
                for( std::vector<std::string>::const_iterator bitr = trigger_bit_list.begin(); bitr != trigger_bit_list.end(); ++bitr ) {
                    std::stringstream ss( *bitr );
                    int bit;
                    ss >> bit;
                    _electronTrigMatchBits.push_back(  bit );
                }
            }
        }
    }
}

bool RunModule::execute( std::vector<ModuleConfig> & configs ) {

    _selectedMuons.clear();
    _selectedElectrons.clear();
    _selectedPhotons.clear();

    // In BranchInit
    CopyInputVarsToOutput();

    // set as a default
    OUT::SampleWeight = 1.0;

    // loop over configured modules
    bool save_event = true;
    BOOST_FOREACH( ModuleConfig & mod_conf, configs ) {
        save_event &= ApplyModule( mod_conf );
    }

    return save_event;

}

bool RunModule::ApplyModule( ModuleConfig & config ) 
{

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
    if( config.GetName() == "FilterBlind" ) {
        keep_evt &= FilterBlind( config );
    }


    return keep_evt;

}

void RunModule::FilterMuon( ModuleConfig & config )
{

    OUT::mu_n          = 0;
    OUT::mu_pt20_n          = 0;
    OUT::mu_pt30_n          = 0;
    OUT::mu_hasTrigMatch->clear();
    OUT::mu_trigMatchMinDr->clear();

    ClearOutputPrefix("mu_");

    for( int idx = 0; idx < IN::mu_n; ++idx ) {

        if( !config.PassFloat( "cut_pt", IN::mu_pt->at(idx) ) ) {
            continue;
        } 
        if( !config.PassFloat( "cut_eta", fabs(IN::mu_eta->at(idx)) ) ) continue;
        if( !config.PassBool( "cut_tight", IN::mu_passTight->at(idx) ) ) continue;

        TLorentzVector mulv;
        mulv.SetPtEtaPhiE( IN::mu_pt->at(idx), 
                           IN::mu_eta->at(idx),
                           IN::mu_phi->at(idx),
                           IN::mu_e->at(idx)
                           );

        _selectedMuons.push_back( mulv );

        float mindr = 101.0;
        for( unsigned hltidx = 0 ; hltidx < IN::HLTObjPt->size(); ++hltidx ) {

            ULong64_t mubits = IN::HLTObjMuTriggers->at(hltidx);

            if( mubits == 0 ) continue;

            TLorentzVector hltlv;
            hltlv.SetPtEtaPhiM( IN::HLTObjPt->at(hltidx),
                                IN::HLTObjEta->at(hltidx),
                                IN::HLTObjPhi->at(hltidx),
                                IN::HLTObjM->at(hltidx)
                              );


            for( std::vector<int>::const_iterator bitr = _muonTrigMatchBits.begin(); bitr != _muonTrigMatchBits.end(); ++bitr ) {

               ULong64_t bitcheck(1);
               bitcheck = bitcheck << ( *bitr );

               if( (mubits & bitcheck) == bitcheck ) {

                   float dr = mulv.DeltaR( hltlv );
                   if( dr < mindr ) {
                       mindr = dr;
                   }
               }
            }
        }

        bool matchTrig = false;
        if( mindr < 0.2 ) {
            matchTrig = true;
        }
        OUT::mu_hasTrigMatch->push_back(matchTrig);
        OUT::mu_trigMatchMinDr->push_back(mindr);

        CopyPrefixIndexBranchesInToOut( "mu_", idx );
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
    OUT::el_hasTrigMatch->clear();
    OUT::el_trigMatchMinDr->clear();

    ClearOutputPrefix("el_");

    for( int idx = 0; idx < IN::el_n; ++idx ) {

        if( !config.PassFloat( "cut_pt", IN::el_pt->at(idx) ) ) continue;
        if( !config.PassFloat( "cut_eta", fabs(IN::el_eta->at(idx)) ) ) continue;
        if( !config.PassBool( "cut_tight", IN::el_passTight->at(idx) ) ) continue;
        if( !config.PassBool( "cut_medium", IN::el_passMedium->at(idx) ) ) continue;
        if( !config.PassBool( "cut_vid_medium", IN::el_passVIDMedium->at(idx) ) ) continue;

        TLorentzVector ellv;
        ellv.SetPtEtaPhiE( IN::el_pt->at(idx), 
                           IN::el_eta->at(idx),
                           IN::el_phi->at(idx),
                           IN::el_e->at(idx)
                           );

        float min_mu_dr = 100.0;
        for( unsigned muidx=0; muidx < _selectedMuons.size(); ++muidx ) {

            float dr = _selectedMuons[muidx].DeltaR( ellv );

            if( dr < min_mu_dr ) {
                min_mu_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_muon_dr", min_mu_dr ) ) continue;

        _selectedElectrons.push_back( ellv );

        float mindr = 101.0;
        for( unsigned hltidx = 0; hltidx < IN::HLTObjPt->size(); ++hltidx ) {

            ULong64_t elbits = IN::HLTObjElTriggers->at(hltidx);

            if( elbits == 0 ) continue;

            TLorentzVector hltlv;
            hltlv.SetPtEtaPhiM( IN::HLTObjPt->at(hltidx),
                                IN::HLTObjEta->at(hltidx),
                                IN::HLTObjPhi->at(hltidx),
                                IN::HLTObjM->at(hltidx)
                              );


            for( std::vector<int>::const_iterator bitr = _electronTrigMatchBits.begin(); bitr != _electronTrigMatchBits.end(); ++bitr ) {

               ULong64_t bitcheck(1);
               bitcheck = bitcheck << ( *bitr );

               if( (elbits & bitcheck) == bitcheck ) {

                   float dr = ellv.DeltaR( hltlv );
                   if( dr < mindr ) {
                       mindr = dr;
                   }
               }
            }
        }

        bool matchTrig = false;
        if( mindr < 0.2 ) {
            matchTrig = true;
        }
        OUT::el_hasTrigMatch->push_back(matchTrig);
        OUT::el_trigMatchMinDr->push_back(mindr);

        CopyPrefixIndexBranchesInToOut( "el_", idx );
        OUT::el_n++;

        if( IN::el_pt->at(idx) > 30 ) {
            OUT::el_pt30_n++;
        }

    }

}

void RunModule::FilterPhoton( ModuleConfig & config ) {

    OUT::ph_n          = 0;
    ClearOutputPrefix("ph_");

    for( int idx = 0; idx < IN::ph_n; ++idx ) {

        if( !config.PassFloat( "cut_pt", IN::ph_pt->at(idx) ) ) continue;
        if( !config.PassFloat( "cut_eta", fabs(IN::ph_eta->at(idx)) ) ) continue;
        if( !config.PassBool( "cut_loose", IN::ph_passLoose->at(idx) ) ) continue;
        if( !config.PassBool( "cut_medium", IN::ph_passMedium->at(idx) ) ) continue;
        if( !config.PassBool( "cut_tight", IN::ph_passTight->at(idx) ) ) continue;
        if( !config.PassBool( "cut_tight", IN::ph_passTight->at(idx) ) ) continue;
        if( !config.PassBool( "cut_eb", IN::ph_IsEB->at(idx) ) ) continue;
        if( !config.PassBool( "cut_ee", IN::ph_IsEE->at(idx) ) ) continue;
        if( !config.PassBool( "cut_CSEV", IN::ph_eleVeto->at(idx) ) ) continue;

        TLorentzVector phlv;
        phlv.SetPtEtaPhiE( IN::ph_pt->at(idx), 
                           IN::ph_eta->at(idx),
                           IN::ph_phi->at(idx),
                           IN::ph_e->at(idx) 
                           );

        float min_mu_dr = 100.0;
        for( unsigned muidx=0; muidx < _selectedMuons.size(); ++muidx ) {

            float dr = _selectedMuons[muidx].DeltaR( phlv );

            if( dr < min_mu_dr ) {
                min_mu_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_muon_dr", min_mu_dr ) ) continue;

        float min_el_dr = 100.0;
        for( unsigned elidx=0; elidx < _selectedElectrons.size(); ++elidx ) {

            float dr = _selectedElectrons[elidx].DeltaR( phlv );

            if( dr < min_el_dr ) {
                min_el_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_electron_dr", min_el_dr ) ) continue;


        _selectedPhotons.push_back( phlv );

        CopyPrefixIndexBranchesInToOut( "ph_", idx );
        OUT::ph_n++;

    }
}

void RunModule::FilterJet( ModuleConfig & config ) const {

    OUT::jet_n          = 0;
    OUT::jet_CSVLoose_n          = 0;
    OUT::jet_CSVMedium_n          = 0;
    OUT::jet_CSVTight_n          = 0;
    ClearOutputPrefix("jet_");

    for( int idx = 0; idx < IN::jet_n; ++idx ) {

        TLorentzVector jetlv;
        jetlv.SetPtEtaPhiE( IN::jet_pt->at(idx), 
                            IN::jet_eta->at(idx),
                            IN::jet_phi->at(idx),
                            IN::jet_e->at(idx) 
                            );

        if( !config.PassFloat( "cut_pt", IN::jet_pt->at(idx) ) ) continue;

        float min_mu_dr = 100.0;
        for( unsigned muidx=0; muidx < _selectedMuons.size(); ++muidx ) {

            float dr = _selectedMuons[muidx].DeltaR( jetlv );

            if( dr < min_mu_dr ) {
                min_mu_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_muon_dr", min_mu_dr ) ) continue;

        float min_el_dr = 100.0;
        for( unsigned elidx=0; elidx < _selectedElectrons.size(); ++elidx ) {

            float dr = _selectedElectrons[elidx].DeltaR( jetlv );

            if( dr < min_el_dr ) {
                min_el_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_electron_dr", min_el_dr ) ) continue;

        float min_ph_dr = 100.0;
        for( unsigned phidx=0; phidx < _selectedPhotons.size(); ++phidx ) {

            float dr = _selectedPhotons[phidx].DeltaR( jetlv );

            if( dr < min_ph_dr ) {
                min_ph_dr = dr;
            }
        }

        if( !config.PassFloat( "cut_photon_dr", min_ph_dr ) ) continue;

        CopyPrefixIndexBranchesInToOut( "jet_", idx );
        OUT::jet_n++;

        float jet_csv = IN::jet_bTagCSVV2->at(idx);

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
    
    if( !config.PassBool( "cut_trig_Ele27_eta2p1_tight", IN::passTrig_HLT_Ele27_eta2p1_WPTight_Gsf) ) keep_event=false;
    if( !config.PassBool( "cut_trig_Mu27_IsoORIsoTk", (IN::passTrig_HLT_IsoMu27 | IN::passTrig_HLT_IsoTkMu27) ) ) keep_event=false;
    if( !config.PassBool( "cut_trig_Mu24_IsoORIsoTk", (IN::passTrig_HLT_IsoMu24 | IN::passTrig_HLT_IsoTkMu24) ) ) keep_event=false;
    if( !config.PassBool( "cut_trig_Mu24_IsoORIsoTk", (IN::passTrig_HLT_IsoMu24 | IN::passTrig_HLT_IsoTkMu24) ) ) keep_event=false;
    if( !config.PassBool( "cut_trig_Mu17_Photon30", (IN::passTrig_HLT_Mu17_Photon30_CaloIdL_L1ISO) ) ) keep_event=false;

    return keep_event;
    
}

void RunModule::BuildEventVars( ModuleConfig & config ) const {


    OUT::m_lep_ph = 0;
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
        bool passLoose        = OUT::ph_passLoose->at(idx);
        bool passTight        = OUT::ph_passTight->at(idx);

        bool eleVeto          = OUT::ph_eleVeto->at(idx);
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

    OUT::ph_hasTruthMatchPh->clear();
    OUT::ph_truthMatchPh_pt->clear();
    OUT::ph_truthMatchPh_dr->clear();

    OUT::el_hasTruthMatchEl->clear();
    OUT::el_truthMatchEl_pt->clear();
    OUT::el_truthMatchEl_dr->clear();

    OUT::mu_hasTruthMatchMu->clear();
    OUT::mu_truthMatchMu_pt->clear();
    OUT::mu_truthMatchMu_dr->clear();

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

            if( abs( OUT::truelep_Id->at(tlepidx) ) != 13 ) continue;

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

            if( abs( OUT::truelep_Id->at(tlepidx) ) != 11 ) continue;

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

    TLorentzVector nu_sum;
    std::vector<TLorentzVector> trueleps;
    std::vector<std::pair<float, unsigned> >  sorted_tleps;

    for( int tlepidx = 0; tlepidx < OUT::truelep_n; ++tlepidx ) {
        int absid = abs(OUT::truelep_Id->at(tlepidx));

        if( absid == 12 || absid == 14 || absid == 16 ) {

            TLorentzVector nulv;
            nulv.SetPtEtaPhiM( OUT::truelep_pt->at(tlepidx),
                               OUT::truelep_eta->at(tlepidx),
                               OUT::truelep_phi->at(tlepidx),
                               0.0
                             );

            nu_sum += nulv;
        }
        if( absid == 11 || absid == 13 ) {
            TLorentzVector leplv;
            leplv.SetPtEtaPhiE( OUT::truelep_pt->at(tlepidx),
                                OUT::truelep_eta->at(tlepidx),
                                OUT::truelep_phi->at(tlepidx),
                                OUT::truelep_e->at(tlepidx)
                              );

            trueleps.push_back( leplv );
            sorted_tleps.push_back( std::make_pair( leplv.Pt(), trueleps.size()-1) );
        }

    }

    OUT::truemt_lep_met_ph = 0;
    OUT::truemt_res = 0;
    if( trueleps.size() > 0 && OUT::trueph_n>0) {
   
        TLorentzVector lep_trans;
        TLorentzVector ph_trans;
        lep_trans.SetPtEtaPhiM( trueleps[sorted_tleps[0].second].Pt(), 0.0, trueleps[sorted_tleps[0].second].Phi(), trueleps[sorted_tleps[0].second].M() );
        ph_trans.SetPtEtaPhiM( OUT::trueph_pt->at(0), 0.0, OUT::trueph_phi->at(0), 0.0 );

        TLorentzVector metlv;
        metlv.SetPtEtaPhiM( nu_sum.Pt(), 0.0, nu_sum.Phi(), 0.0 );

        OUT::truemt_res = ( lep_trans + ph_trans + metlv ).M();
           
        std::sort(sorted_tleps.rbegin(), sorted_tleps.rend());
        OUT::truemt_lep_met_ph = Utils::calc_mt( trueleps[sorted_tleps[0].second], nu_sum);
    }



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
        if( OUT::EvtIsRealData ) keep_event=false;
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

RunModule::RunModule() {
    _m_w = 80.385;
    _isData = false;
}

