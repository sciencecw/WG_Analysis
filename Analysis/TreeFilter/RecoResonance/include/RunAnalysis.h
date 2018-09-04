#ifndef RUNANALYSIS_H
#define RUNANALYSIS_H

#include "AnalysisBase.h"


#include <vector>


#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"

// The RunModule inherits from RunModuleBase (an Abstract Base Class )
// defined in the Core package so that all
// RunModules present a common interface in a Run function
// This allows the code defined in this package
// to be run from the Core package to minimize
// code duplication in each module
class RunModule : public virtual RunModuleBase {

    public :

        RunModule();

        // The run function must exist and be defined exactly as this
        // because it is defined in RunModuleBase 
        // In src/RunModule.cxx all the analysis is defind
	// in this RunModule function
        void initialize( TChain * chain, TTree *outtree, TFile *outfile, 
			 const CmdOptions & options,
			 std::vector<ModuleConfig> & config) ;
        bool execute( std::vector<ModuleConfig> & config ) ;

        void finalize( ) {};

        // The ApplyModule function calls any other module defined below
        // in src/RunModule.cxx.  This funciton is not strictly required
        // but its a consistent way to apply modules
        bool ApplyModule         ( ModuleConfig & config ) ;

        void FilterMuon          ( ModuleConfig & config ) ;
        void FilterElectron      ( ModuleConfig & config ) ;
        void FilterPhoton        ( ModuleConfig & config ) ;
        void FilterJet           ( ModuleConfig & config ) const;
        void BuildEventVars      ( ModuleConfig & config ) const;
        void MakePhotonCountVars ( ModuleConfig & config ) const;
        void BuildTruth          ( ModuleConfig & config ) const;
        bool FilterTrigger       ( ModuleConfig & config ) ;
        bool FilterMET           ( ModuleConfig & config ) ;
        bool FilterEvent         ( ModuleConfig & config ) const;
        void WeightEvent         ( ModuleConfig & config ) const;

        bool FilterDataQuality   ( ModuleConfig & config ) const;
        bool FilterBlind         ( ModuleConfig & config ) const;

        bool HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR ) const;
        bool HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR, float &minDR ) const;
        bool HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR, float &minDR, TLorentzVector &matchLV ) const;
        bool HasTruthMatch( const TLorentzVector & objlv, const std::vector<int> & matchPID, float maxDR, float &minDR, TLorentzVector &matchLV, int &matchMotherPID, int &matchParentage ) const;

        bool get_constriained_nu_pz( const TLorentzVector lepton, TLorentzVector & metlv ) const ;
        bool calc_constrained_nu_momentum( const TLorentzVector lepton, const TLorentzVector metlv, float &result ) const ;
        bool solve_quadratic( float Aval, float Bval, float Cval, float & solution1, float &solution2 ) const; 
        void calc_corr_iso( float chIso, float phoIso, float neuIso, float rho, float eta, float &chIsoCorr, float &phoIsoCorr, float &neuIsoCorr )  const;

        float calc_pu_weight( float puval, float mod=1.0) const;

    private :

        float _m_w;
        float _m_z;
        bool _isData;

        bool _eval_mu_loose    ;
        bool _eval_mu_medium   ;
        bool _eval_mu_tight    ;
        bool _eval_ph_tight    ;
        bool _eval_ph_medium   ;
        bool _eval_ph_loose    ;
        bool _eval_el_tight    ;
        bool _eval_el_medium   ;
        bool _eval_el_loose    ;
        bool _eval_el_veryloose;
        bool _needs_nlo_weght  ;

        std::map<int, std::vector<int> > _quality_map;

        std::vector<int> _muonTrigMatchBits;
        std::vector<int> _electronTrigMatchBits;

        std::map<int, bool> metfilterResults;
        std::map<int, bool> triggerResults;

        TFile * _puweight_sample_file;
        TFile * _puweight_data_file;
        TH1F * _puweight_sample_hist;
        TH1D * _puweight_data_hist;

        TChain * _input_tree;

};

// Ouput namespace 
// Declare any output variables that you'll fill here
namespace OUT {

    Int_t mu_pt20_n;
    Int_t mu_pt30_n;
    std::vector<Bool_t> *mu_passTight;
    std::vector<Bool_t> *mu_passMedium;
    std::vector<Bool_t> *mu_passLoose;
    std::vector<Bool_t> *mu_hasTruthMatchMu;
    std::vector<float>  *mu_truthMatchMu_dr;
    std::vector<float>  *mu_truthMatchMu_pt;
    std::vector<float>  *mu_hasTrigMatch;
    std::vector<float>  *mu_trigMatch_dr;

    Int_t el_pt30_n;
    std::vector<Bool_t> *el_passVeryLoose;
    std::vector<Bool_t> *el_passLoose;
    std::vector<Bool_t> *el_passMedium;
    std::vector<Bool_t> *el_passTight;
    std::vector<Bool_t> *el_hasTrigMatch;
    std::vector<float>  *el_trigMatch_dr;
    std::vector<Bool_t> *el_hasTruthMatchEl;
    std::vector<float>  *el_truthMatchEl_dr;
    std::vector<float>  *el_truthMatchEl_pt;

    std::vector<float> *ph_chIsoCorr;
    std::vector<float> *ph_neuIsoCorr;
    std::vector<float> *ph_phoIsoCorr;
    std::vector<float> *ph_min_el_dr;
    std::vector<Bool_t> *ph_IsEB;
    std::vector<Bool_t> *ph_IsEE;
    std::vector<Bool_t> *ph_passTight;
    std::vector<Bool_t> *ph_passMedium;
    std::vector<Bool_t> *ph_passLoose;
    std::vector<Bool_t> *ph_passLooseNoSIEIE;
    std::vector<Bool_t> *ph_passHOverELoose;
    std::vector<Bool_t> *ph_passHOverEMedium;
    std::vector<Bool_t> *ph_passHOverETight;
    std::vector<Bool_t> *ph_passSIEIELoose;
    std::vector<Bool_t> *ph_passSIEIEMedium;
    std::vector<Bool_t> *ph_passSIEIETight;
    std::vector<Bool_t> *ph_passChIsoCorrLoose;
    std::vector<Bool_t> *ph_passChIsoCorrMedium;
    std::vector<Bool_t> *ph_passChIsoCorrTight;
    std::vector<Bool_t> *ph_passNeuIsoCorrLoose;
    std::vector<Bool_t> *ph_passNeuIsoCorrMedium;
    std::vector<Bool_t> *ph_passNeuIsoCorrTight;
    std::vector<Bool_t> *ph_passPhoIsoCorrLoose;
    std::vector<Bool_t> *ph_passPhoIsoCorrMedium;
    std::vector<Bool_t> *ph_passPhoIsoCorrTight;

    std::vector<Bool_t> *ph_hasTruthMatchPh;
    std::vector<float>  *ph_truthMatchPh_dr;
    std::vector<float>  *ph_truthMatchPh_pt;

    Int_t jet_CSVLoose_n;
    Int_t jet_CSVMedium_n;
    Int_t jet_CSVTight_n;
    std::vector<Bool_t>  *jet_IdLoose;
    std::vector<Bool_t>  *jet_IdTight;
    std::vector<Bool_t>  *jet_IdTightLep;
    
    float m_lep_ph;
    float m_lpp_best;
    float m_lep_ph_ph;
    std::vector<float> *m_lep_ph_comb_leadLep;
    std::vector<float> *m_lep_ph_comb_sublLep;
    float m_lep_met_ph;
    float m_mt_lep_met_ph;
    float m_mt_lep_met_ph_forcewmass;
    float mt_w;
    float mt_res;
    float mt_lep_ph;
    float dphi_lep_ph;
    float dr_lep_ph;
    float mt_lep_met;
    float m_lep_met;
    float pt_lep_met;
    float dphi_lep_met;
    float mt_lep_met_ph;
    float mt_lep_met_ph_inv;
    float dphi_ph_w;
    float pt_lep_met_ph;
    float RecoWMass;
    float recoM_lep_nu_ph;
    float recoMet_eta;
    float recoW_pt;
    float recoW_eta;
    float recoW_phi;
    float m_ll;
    Bool_t nu_z_solution_success;
   
    float leadjet_pt;
    float subljet_pt;
    float leaddijet_m;
    float leaddijet_pt;
    float massdijet_m;
    float massdijet_pt;

    Bool_t isBlinded;
    float SampleWeight;

    Int_t ph_loose_n;
    Int_t ph_medium_n;
    Int_t ph_tight_n;

    Int_t ph_mediumPassPSV_n;
    Int_t ph_mediumFailPSV_n;
    Int_t ph_mediumPassCSEV_n;
    Int_t ph_mediumFailCSEV_n;

    Int_t ph_mediumPassEleOlap_n;
    Int_t ph_mediumPassEleOlapPassCSEV_n;
    Int_t ph_mediumPassEleOlapFailCSEV_n;

    Int_t ph_mediumNoSIEIE_n;
    Int_t ph_mediumNoChIso_n;
    Int_t ph_mediumNoNeuIso_n;
    Int_t ph_mediumNoPhoIso_n;

    Int_t ph_mediumNoSIEIENoChIso_n;
    Int_t ph_mediumNoSIEIENoNeuIso_n;
    Int_t ph_mediumNoSIEIENoPhoIso_n;
    Int_t ph_mediumNoChIsoNoPhoIso_n;
    Int_t ph_mediumNoChIsoNoNeuIso_n;
    Int_t ph_mediumNoPhoIsoNoNeuIso_n;

    Int_t ph_mediumNoSIEIEPassPSV_n;
    Int_t ph_mediumNoChIsoPassPSV_n;
    Int_t ph_mediumNoNeuIsoPassPSV_n;
    Int_t ph_mediumNoPhoIsoPassPSV_n;

    Int_t ph_mediumNoSIEIEFailPSV_n;
    Int_t ph_mediumNoChIsoFailPSV_n;
    Int_t ph_mediumNoNeuIsoFailPSV_n;
    Int_t ph_mediumNoPhoIsoFailPSV_n;

    Int_t ph_mediumNoSIEIEPassCSEV_n;
    Int_t ph_mediumNoChIsoPassCSEV_n;
    Int_t ph_mediumNoNeuIsoPassCSEV_n;
    Int_t ph_mediumNoPhoIsoPassCSEV_n;

    Int_t ph_mediumNoSIEIEFailCSEV_n;
    Int_t ph_mediumNoChIsoFailCSEV_n;
    Int_t ph_mediumNoNeuIsoFailCSEV_n;
    Int_t ph_mediumNoPhoIsoFailCSEV_n;

    std::vector<Int_t> *ptSorted_ph_loose_idx;
    std::vector<Int_t> *ptSorted_ph_medium_idx;
    std::vector<Int_t> *ptSorted_ph_tight_idx;

    std::vector<Int_t> *ptSorted_ph_mediumPassPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumFailPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumPassCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumFailCSEV_idx;

    std::vector<Int_t> *ptSorted_ph_mediumPassEleOlap_idx;
    std::vector<Int_t> *ptSorted_ph_mediumPassEleOlapPassCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumPassEleOlapFailCSEV_idx;

    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIE_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoNeuIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoPhoIso_idx;

    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIENoChIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIENoNeuIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIENoPhoIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIsoNoPhoIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIsoNoNeuIso_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoPhoIsoNoNeuIso_idx;

    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIEPassPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIsoPassPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoNeuIsoPassPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoPhoIsoPassPSV_idx;

    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIEFailPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIsoFailPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoNeuIsoFailPSV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoPhoIsoFailPSV_idx;

    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIEPassCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIsoPassCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoNeuIsoPassCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoPhoIsoPassCSEV_idx;

    std::vector<Int_t> *ptSorted_ph_mediumNoSIEIEFailCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoChIsoFailCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoNeuIsoFailCSEV_idx;
    std::vector<Int_t> *ptSorted_ph_mediumNoPhoIsoFailCSEV_idx;

    Int_t trueph_n;
    std::vector<float> *trueph_pt;
    std::vector<float> *trueph_eta;
    std::vector<float> *trueph_phi;
    std::vector<Int_t> *trueph_motherPID;
    std::vector<Int_t> *trueph_status;
    std::vector<Int_t> *trueph_nMatchingLep;
    Int_t truephIPFS_n;

    Int_t truelep_n;
    std::vector<float> *truelep_pt;
    std::vector<float> *truelep_eta;
    std::vector<float> *truelep_phi;
    std::vector<float> *truelep_e;
    std::vector<Int_t> *truelep_motherPID;
    std::vector<Int_t> *truelep_status;
    std::vector<Int_t> *truelep_PID;

    Int_t truenu_n;

    float truelepnu_m;
    float truelepnuph_m;
    float truelepph_dr;
    float truemt_lep_met_ph;
    float truemt_res;
    float truemt_res_l23;
    float truemt_res_lO;

    float trueht;

    Bool_t isWMuDecay;
    Bool_t isWElDecay;
    Bool_t isWTauDecay;
    Bool_t isWTauElDecay;
    Bool_t isWTauMuDecay;
    Int_t WIDStep;

    Bool_t PassQuality;

    float NLOWeight;

    float PUWeight;
    float PUWeightUP5;
    float PUWeightUP10;
    float PUWeightDN5;
    float PUWeightDN10;

};

#endif
