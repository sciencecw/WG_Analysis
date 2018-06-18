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

    // *************************
    // initialize trees
    // *************************
    InitINTree(chain);
    InitOUTTree( outtree );
    
	_outfile = outfile;
    // *************************
    // Set defaults for added output variables
    // *************************
    // Examples :
    OUT::ph_n              = 0;
    OUT::ph_pt             = 0;
    OUT::ph_eta            = 0;
    OUT::ph_phi            = 0;
    OUT::ph_e              = 0;

    // *************************
    // Declare Branches
    // *************************

    // Examples :
    outtree->Branch("ph_n"        , &OUT::ph_n        );
    outtree->Branch("ph_pt"       , &OUT::ph_pt       );
    outtree->Branch("ph_eta"      , &OUT::ph_eta      );
    outtree->Branch("ph_phi"      , &OUT::ph_phi      );
    outtree->Branch("ph_e"        , &OUT::ph_e        );
    outtree->Branch("ph_sigmaIEIE", &OUT::ph_sigmaIEIE);

	hout = new TH1F("hout","m_lep_ph",100,0,200);
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

bool RunModule::ApplyModule( ModuleConfig & config ) const {

    // This bool is used for filtering
    // If a module implements an event filter
    // update this variable and return it
    // to apply the filter
    bool keep_evt = true;

    // This part is a bit hacked.  For each module that
    // you write below, you have to put a call to that
    // function with a matching name here.
    // The name is used to match the name used
    // in the python configuration.
    // There are fancy ways to do this, but it
    // would require the code to be much more complicated
    //
    // Example :
    if( config.GetName() == "BuildPhoton" ) {
        BuildPhoton( config );
    }

    // If the module applies a filter the filter decision
    // is passed back to here.  There is no requirement
    // that a function returns a bool, but
    // if you want the filter to work you need to do this
    //
    // Example :
    if( config.GetName() == "FilterEvent" ) {
        keep_evt &= FilterEvent( config );
    }

    return keep_evt;

}

// ***********************************
//  Define modules here
//  The modules can do basically anything
//  that you want, fill trees, fill plots, 
//  caclulate an event filter
// ***********************************
//
// Examples :

void RunModule::BuildPhoton( ModuleConfig & config ) const {

    OUT::ph_pt         -> clear();
    OUT::ph_eta        -> clear();
    OUT::ph_phi        -> clear();
    OUT::ph_e          -> clear();
    OUT::ph_sigmaIEIE  -> clear();
    OUT::ph_n          = 0;

    // Check for preprocessor defined variables 
    // (set in the generated c++ code) to avoid
    // compilation failures if a certain branch
    // does not exist.  
    // You don't have to be pendantic about it
    // and add a check for every branch, but in this
    // template it is done this way to ensure compilation
    // regardles of the input file used
	
	//event condition
	if  (IN::ph_mediumPassCSEV_n>=1 && IN::el_n>=1 && fabs(IN::ph_eta[0][IN::ptSorted_ph_mediumPassCSEV_idx[0][0]])<2 && IN::ph_pt[0][IN::ptSorted_ph_mediumPassCSEV_idx[0][0]]<30){
		hout->Fill(IN::m_lep_ph);
	} 
	return;
	if ( !(IN::ph_mediumPassCSEV_n>=1 && IN::el_n>=1 && fabs(IN::m_lep_ph-60)<5 && IN::ph_pt[0][IN::ptSorted_ph_mediumPassCSEV_idx[0][0]]<30&& fabs(IN::ph_eta[0][IN::ptSorted_ph_mediumPassCSEV_idx[0][0]])<2)) return;

	std::cout<<IN::eventNumber<<" ph_pt0" <<(*IN::ph_pt)[0]<<std::endl; 
	//ele 
	std::cout<< "ele pt: ";
	for (float pt: (*IN::el_pt)) std::cout<< pt<<" "; 
	std::cout<< std::endl; 
	std::cout<< "ele eta: ";
	for (float pt: (*IN::el_eta)) std::cout<< pt<<" "; 
	std::cout<< std::endl; 
	std::cout<< "ele phi: ";
	for (float pt: (*IN::el_phi)) std::cout<< pt<<" "; 
	std::cout<< std::endl; 
	// list index passing eOlap && eVeto 
	for (int idx: (*IN::ptSorted_ph_mediumPassCSEV_idx)) std::cout<< idx<< " ";
	std::cout<< std::endl; 

	std::cout<< "ph pt: ";
	//photon info
    for (int i =0 ;i<IN::ph_n;i++){
			std::cout<< (*IN::ph_pt)[i] <<" ";  
	}
	std::cout<< std::endl; 
	std::cout<< "ph eta: ";
    for (int i =0 ;i<IN::ph_n;i++){
			std::cout<< (*IN::ph_eta)[i] <<" ";  
	}
	std::cout<< std::endl; 
	std::cout<< "ph phi: ";
    for (int i =0 ;i<IN::ph_n;i++){
			std::cout<< (*IN::ph_phi)[i] <<" ";  
	}
	std::cout<< std::endl; 

}

// This is an example of a module that applies an
// event filter.  Note that it returns a bool instead
// of a void.  In principle the modules can return any
// type of variable, you just have to handle it
// in the ApplyModule function

bool RunModule::FilterEvent( ModuleConfig & config ) const {

    bool keep_event = true;

    #ifdef EXISTS_nPho
    int nPho = IN::nPho;
    if( !config.PassInt("cut_nPho", nPho ) )
      keep_event = false;

    #endif
    return keep_event;
    
}

void RunModule::finalize() {
	_outfile->cd();
	hout->Write();
}
