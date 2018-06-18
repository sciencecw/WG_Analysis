#include "AnalysisBase.h"
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/uuid/uuid.hpp>           
#include <boost/uuid/uuid_generators.hpp> 
#include <boost/uuid/uuid_io.hpp> 
#include <sys/types.h>
#include <sys/stat.h>
#include <ctime>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "TFile.h"

AnaConfig::AnaConfig()
{
}

void AnaConfig::AddModule(  ModuleConfig module ) 
{
    confs.push_back( module );
}

ModuleConfig::ModuleConfig( const std::string &_name ) :
    name( _name )
{
}

const ModuleConfig AnaConfig::getEntry( unsigned int i ) const {

    if( i > confs.size()-1 ) {
        std::cout << "AnaConfig::getEntry - ERROR : entry exceeds size" << std::endl;
        return ModuleConfig( "" );
    }

    return confs[i];
}

void AnaConfig::Run( RunModuleBase & runmod, const CmdOptions & options ) {

    if( options.files.size() == 0 ) {
        runmod.initialize( 0, 0, 0, options, getEntries() );
        runmod.execute( getEntries() );
        runmod.finalize();
    }

    int total_njobs = 0;
    for( unsigned fidx = 0; fidx < options.files.size(); ++fidx ) {
        total_njobs += options.files[fidx].jobs.size();
    }

    // loop over subjobs
    int jobidx = -1;
    for( unsigned fidx = 0; fidx < options.files.size(); ++fidx ) {
        TChain *chain = new TChain(options.treeName.c_str() );

        std::vector<std::string> copied_files;
        if( options.copyInputFiles ) {
            boost::uuids::random_generator generator;
            BOOST_FOREACH( const std::string &fpath, options.files[fidx].files ) {
                std::vector<std::string> filepath_split = Tokenize( fpath, "/" );
                std::string fname = filepath_split[filepath_split.size()-1];

                boost::uuids::uuid file_uuid = generator();
                std::stringstream file_uuid_ss;
                file_uuid_ss << file_uuid << ".root";

                std::string filename = file_uuid_ss.str();

                copied_files.push_back( filename );

                boost::filesystem::copy_file(fpath, filename);
                chain->Add( filename.c_str() );
                std::cout << "Add input file " << filename << std::endl;
            }
        }
        else {
            BOOST_FOREACH( const std::string &fpath, options.files[fidx].files ) {
              chain->Add( fpath.c_str() );
              std::cout << "Add input file " << fpath << std::endl;
            }
        }

        // Loop over job ranges
        for( unsigned jidx = 0 ; jidx < options.files[fidx].jobs.size(); ++jidx ) {
            jobidx++;

            int jobid = options.files[fidx].jobs[jidx].first;
            int minevt = options.files[fidx].jobs[jidx].second.first;
            int maxevt = options.files[fidx].jobs[jidx].second.second;

            // Transform the job id into a string.  Eg job #0 -> Job_0000
            std::stringstream jobstrss("");
            jobstrss << "Job_" << std::setw(4) << std::setfill('0') << jobid;
            std::string jobstr = jobstrss.str();

            std::string outputDir = options.outputDir;
            //if( total_njobs > 1 || jobid != 0 ) {}
            // update the output directory to the job directory
            // this now happens in all cases.  
            outputDir += "/" + jobstr;
            
            // create the output directory
            std::cout << "mkdir " << outputDir << std::endl;
            int res = boost::filesystem::create_directories(outputDir.c_str());
            if( res != 1 ) {
                std::cout << "WARNING -- Directory creation failed " << res << std::endl;
            }
            std::string filepath = outputDir + "/" + options.outputFile;

            TFile * outfile  = TFile::Open(filepath.c_str(), "RECREATE");

            outfile->cd();
            TTree * outtree  = 0;
            std::vector<std::string> name_tok = Tokenize( chain->GetName(), "/" ); 

            std::string dir_path = "";
            if( name_tok.size() > 1 ) {
                // if the tree is in a directory
                // reproduce the directory structure
                // and put the tree in there
                for( unsigned i = 0; i < name_tok.size()-1; i++ ) {
                    dir_path = dir_path + name_tok[i]+"/";
                    outfile->mkdir(dir_path.c_str());
                }
                outfile->cd(dir_path.c_str());

                std::string name = name_tok[name_tok.size()-1];
                std::cout << " Create output tree in " << dir_path << "/" << name << std::endl;
            
                outtree = new TTree(name.c_str(), name.c_str());
            }
            else {
                outtree = new TTree(chain->GetName(), chain->GetName());
                outtree->SetDirectory(outfile);
            }

            outfile->cd();
            TH1F * hfilter = new TH1F("filter", "filter", 2, 0, 2);
            hfilter->GetXaxis()->SetBinLabel(1, "Total");
            hfilter->GetXaxis()->SetBinLabel(2, "Filter");

            runmod.initialize( chain, outtree, outfile, options, getEntries() );

            if( maxevt == 0 ) {
                maxevt = chain->GetEntries();
            }

            int n_saved = 0;
            int n_evt = 0;
            std::cout << "Will analyze " << maxevt-minevt << " events between " << minevt << " and " << maxevt << std::endl;
            boost::posix_time::ptime time_start = boost::posix_time::microsec_clock::local_time();
            for( int cidx = minevt; cidx < maxevt; cidx++ ) {

                //if( n_evt == 0 ) time_start = boost::posix_time::microsec_clock::local_time();
                if( n_evt % options.nPrint == 0 && n_evt > 0 ) {
                    boost::posix_time::ptime time_now = boost::posix_time::microsec_clock::local_time();
                    boost::posix_time::time_duration deltat = time_now - time_start;
                    std::cout << "Processed " << n_evt << " entries in " << std::fixed << std::setprecision(2) << deltat.total_milliseconds()/1000. << " seconds" << std::endl;
                    //time_start = boost::posix_time::microsec_clock::local_time();
                }

                n_evt++;

                chain->GetEntry(cidx);

                bool save_event = runmod.execute( getEntries() );

                hfilter->Fill(0);
                if( save_event ) {
                    if( !options.disableOutputTree ) outtree->Fill();
                    hfilter->Fill(1);
                    n_saved++;
                }
            }
            std::cout << "Wrote " << n_saved << " events" << std::endl;

            runmod.finalize();

            //bool has_any_cutflows=false;
            BOOST_FOREACH( ModuleConfig & conf, getEntries() ) {
                outfile->cd();

                if( conf.hasCutFlows() ) {
                    //has_any_cutflows = true;
                    outfile->mkdir( conf.GetName().c_str() );
                    outfile->cd(conf.GetName().c_str() );

                    BOOST_FOREACH( const std::string & cutname, conf.getCutFlows()[0].getOrder() ) {
                        std::string befname = cutname+"_before";
                        std::string aftname = cutname+"_after";
                        if( conf.getCutFlows()[0].hasHist(befname) ) {
                            conf.getCutFlows()[0].getHist(befname).Write();
                        }
                        if( conf.getCutFlows()[0].hasHist(aftname) ) {
                            conf.getCutFlows()[0].getHist(aftname).Write();
                        }
                    }
                    conf.WriteCutFlowHists( outfile );
                    conf.PrintCutFlows( );
                }
            }

            outfile->cd(dir_path.c_str());
            outtree->Write();
            hfilter->Write();
            outfile->Close();

            std::vector<std::string> output_files;
            output_files.push_back(filepath);

            if( options.transferToStorage ) {
                std::string storage_dir = options.storagePath;
                std::string eos = "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select";

                if( jobidx == 0 ) { // make the directory the first time
                    std::string mkdir_cmd = eos + " mkdir " + options.storagePath;
                    std::cout << mkdir_cmd << std::endl;
                    system( mkdir_cmd.c_str() );
                    // copy the configuration script to this directory for future reference
                    std::vector<std::string> conf_file_tok = Tokenize( options.config_file, "/" );
                    std::string conf_file_name = conf_file_tok[conf_file_tok.size()-1];

                    std::string cpy_cmd = eos + " cp " + options.config_file + " " + options.storagePath + "/" + conf_file_name ;
                    std::cout << cpy_cmd << std::endl;
                    system( cpy_cmd.c_str() );
                }
                // update the output directory to the job directory
                // this now happens in all cases.  
                storage_dir += "/" + jobstr;
                std::string mkdir_cmd = eos + " mkdir " + storage_dir;
                std::cout << mkdir_cmd << std::endl;
                system( mkdir_cmd.c_str() );

                BOOST_FOREACH( const std::string & path, output_files ) {
                    std::string copy_cmd = eos + " cp " + path + " " + storage_dir + "/" + options.outputFile;
                    std::cout << copy_cmd << std::endl;
                    system( copy_cmd.c_str() );
                }
                BOOST_FOREACH( const std::string & path, output_files ) {
                    std::string rm_cmd = " rm " + path;
                    std::cout << rm_cmd << std::endl;
                    system( rm_cmd.c_str() );
                }
            }
        }
        if( options.copyInputFiles ) {

            BOOST_FOREACH( const std::string &fname, copied_files ) {
                boost::filesystem::remove( fname );
            }
        }
    }
   
}

std::vector<ModuleConfig> & AnaConfig::getEntries() {
  return confs;
}

Cut::Cut( CutType::Op in_op, CutType::Type in_type, CutType::Comp in_comp,
               bool in_val_bool, int in_val_int, float in_val_float ) : 
    op( in_op ),
    type( in_type ),
    comp( in_comp ),
    val_bool( in_val_bool ),
    val_int ( in_val_int ),
    val_float( in_val_float )
{
}

void Cut::Print() const{ 

    std::string op_str = "UNKNOWN";
    if( op == CutType::EQUAL_TO ) op_str = "EQUAL_TO";
    if( op == CutType::GREATER_THAN) op_str = "GREATER_THAN";
    if( op == CutType::LESS_THAN) op_str = "LESS_THAN";
    if( op == CutType::GREATER_THAN_OR_EQUAL_TO) op_str = "GREATER_THAN_OR_EQUAL_TO";
    if( op == CutType::LESS_THAN_OR_EQUAL_TO) op_str = "LESS_THAN_OR_EQUAL_TO";
    if( op == CutType::NOT_EQUAL_TO ) op_str = "NOT_EQUAL_TO";
    if( op == CutType::OTHER) op_str = "OTHER";

    std::string type_str = "UNKNOWN";
    if( type == CutType::BOOL ) type_str = "BOOL";
    if( type == CutType::INT ) type_str = "INT";
    if( type == CutType::FLOAT) type_str = "FLOAT";

    std::cout << "op = " << op_str << " type = " << type_str << " bool val = " << val_bool << " int val = " << val_int << " float val = " << val_float << std::endl;
}

CutConfig::CutConfig( const std::string &cut_str ) {

    // cut string should have the format name [val]
    // find the location of the brackets
    std::size_t posbeg = cut_str.find("[");
    std::size_t posend = cut_str.find("]");

    // throw an error if the formatting is incorrect
    if( posbeg == std::string::npos || posend == std::string::npos ) {

      std::cout << "CutConfig - ERROR : cut string should have the format name [val]"
                << "Instead got " << cut_str << std::endl;
      return;
    }

    // get the cut name values
    std::string name    = cut_str.substr( 0, posbeg);
    std::string val_str = cut_str.substr( posbeg+1, posend-posbeg-1 );

    bool should_invert = false;
    // if this cut should be inverted there is a ! in the string.
    // Remove ! from the name if it is there
    std::size_t posinv = name.find("!");
    if( posinv != std::string::npos ) {
        should_invert = true;
        name = name.substr(0, posinv);
    }

    // remove whitespace from the name
    boost::algorithm::trim(name);

    // handle the case if a double & or | was used in the cut definition
    // just remove one of them
    size_t loc;
    while( loc = val_str.find("&&"), loc != std::string::npos ) {
        val_str.erase(loc, 1);
    }
    while( loc = val_str.find("||"), loc != std::string::npos ) {
        val_str.erase(loc, 1);
    }


    // treat a comma or & as an and operator
    bool has_comma_op = val_str.find(",") != std::string::npos;
    bool has_and_op   = val_str.find("&") != std::string::npos;
    bool has_or_op    = val_str.find("|") != std::string::npos;

    // throw an error if both ANDs and ORs are present in the cut
    // this case is not handled yet ( requires order of operation handling )
    if( ( has_comma_op || has_and_op ) && has_or_op ) {
        std::cout << "CutConfig - ERROR : both an AND operation and an OR operation are present in the cut.  This case is not yet handled" << std::endl;
    }
          
    std::vector<Cut> conf_cut;

    std::vector<std::string> each_cut; 
    // split by the operator.  If no operator is 
    // found, then just use the whole string
    if     ( has_comma_op ) each_cut = Tokenize(val_str, ",");
    else if( has_and_op )   each_cut = Tokenize(val_str, "&");
    else if( has_or_op )    each_cut = Tokenize(val_str, "|");
    else                    each_cut.push_back(val_str);

    CutType::Comp cut_comp = CutType::AND;
    if( has_or_op ) {
        cut_comp = CutType::OR;
    }
    
    BOOST_FOREACH( const std::string & val_orig, each_cut ) {
      
        // copy the value string and set it to lower case (mainly for True/False -> true/false )
        std::string val( val_orig );
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);

        CutType::Op cut_op;
        CutType::Type cut_type;
        bool  cut_val_bool  = 0;
        int   cut_val_int   = 0;
        float cut_val_float = 0;

        // check if it is a boolean cut first
        if( val.find("true") != std::string::npos ) {
            cut_type = CutType::BOOL;
            if( val.find( "!=" ) != std::string::npos ) {
                cut_op = CutType::NOT_EQUAL_TO;
            }
            else {
                cut_op   = CutType::EQUAL_TO;
            }
            cut_val_bool = true;
        }
        else if( val.find("false") != std::string::npos ) {
            cut_type = CutType::BOOL;
            if( val.find( "!=" ) != std::string::npos ) {
                cut_op = CutType::NOT_EQUAL_TO;
            }
            else {
                cut_op   = CutType::EQUAL_TO;
            }
            cut_val_bool = false;
        }
        else {
            cut_op = attempt_logicalop_parse( val, cut_type, cut_val_int, cut_val_float );
        }

        conf_cut.push_back( Cut(cut_op, cut_type, cut_comp, cut_val_bool, cut_val_int, cut_val_float) );
    }

    SetName( name );
    SetCuts( conf_cut );
    SetIsInverted( should_invert );
    SetCompOp( cut_comp );

}

CutType::Op CutConfig::attempt_logicalop_parse( const std::string & val, CutType::Type & type,
                                                  int & cut_val_int, float &  cut_val_float )
{
    
    CutType::Op cut_op;

    // create a string that will have the operator removed
    std::string mod_val;

    std::size_t pos;
    if( (pos = val.find( "==" ) ) != std::string::npos ) {
        cut_op = CutType::EQUAL_TO;
        mod_val = val.substr( pos+2 );
    }
    else if( (pos = val.find( ">=" ) ) != std::string::npos ) {
        cut_op = CutType::GREATER_THAN_OR_EQUAL_TO;
        mod_val = val.substr( pos+2 );
    }
    else if( (pos = val.find( "<=" ) ) != std::string::npos ) {
        cut_op = CutType::LESS_THAN_OR_EQUAL_TO;
        mod_val = val.substr( pos+2 );
    }
    else if( (pos = val.find( ">" ) ) != std::string::npos ) {
        cut_op = CutType::GREATER_THAN;
        mod_val = val.substr( pos+1 );
    }
    else if( (pos = val.find( "<" ) ) != std::string::npos ) {
        cut_op = CutType::LESS_THAN;
        mod_val = val.substr( pos+1 );
    }
    else if( (pos = val.find( "!=" ) ) != std::string::npos ) {
        cut_op = CutType::NOT_EQUAL_TO;
        mod_val = val.substr( pos + 2 );
    }
    // grab single = but warn
    else if( (pos = val.find( "=" ) ) != std::string::npos ) {
        cut_op = CutType::EQUAL_TO;
        mod_val = val.substr( pos+1 );
        std::cout << "CutConfig::attempt_logicalop_parse - WARNING : Interpreting single = as == for cut " << val << std::endl;
    }
    else {
        cut_op = CutType::OTHER;
    }

    boost::algorithm::trim(mod_val);

    // now grab the value
    if( cut_op != CutType::OTHER ) {

        // Check for float by the presence of a decimal point. 
        // otherwise its an int
        // May need to fix this later
        
        if( mod_val.find(".") != std::string::npos ) {
            type = CutType::FLOAT;
            std::stringstream ss( mod_val );
            ss >> cut_val_float;
        }
        else {
            type = CutType::INT;
            std::stringstream ss( mod_val );
            ss >> cut_val_int;
            // Even if this is configured as an
            // int it can be compared to a float
            // so store it as a float too
            cut_val_float = cut_val_int;
        }
            
    }

    return cut_op;
}

bool CutConfig::PassFloat(const std::string &name, const float cutval ) const {

    int n_pass = 0;
    int n_cuts = 0;
    BOOST_FOREACH( const Cut & cut, GetCuts() ) {
        n_cuts++;

        if( cut.type != CutType::FLOAT && cut.type != CutType::INT ) {
          std::cout << "CutConfig::PassFloat - ERROR : Float cut requested for cut " << name << " but cut was not configured as a float " << std::endl;
          continue;
        }

        if( cut.op == CutType::GREATER_THAN_OR_EQUAL_TO ) {
          if( (cutval >= cut.val_float) ) n_pass++;
        }
        if( cut.op == CutType::LESS_THAN_OR_EQUAL_TO ) {
          if( (cutval <= cut.val_float) ) n_pass++;
        }
        if( cut.op == CutType::GREATER_THAN ) {
          if( (cutval > cut.val_float) ) n_pass++;
        }
        if( cut.op == CutType::LESS_THAN ) {
          if( (cutval < cut.val_float) ) n_pass++;
        }
        if( cut.op == CutType::EQUAL_TO ) {
          if( (cutval == cut.val_float) ) n_pass++;
          std::cout << "CutConfig::PassFloat - WARNING : EQUAL_TO operator used for float comparison for cut " << name << std::endl;
        }
        if( cut.op == CutType::NOT_EQUAL_TO ) {
          if( (cutval != cut.val_float) ) n_pass++;
          std::cout << "CutConfig::PassFloat - WARNING : NOT_EQUAL_TO operator used for float comparison for cut " << name << std::endl;
        }

    }

    // Handle AND vs OR
    bool pass_cuts = false;
    if( comp == CutType::AND ) {
        pass_cuts = ( n_cuts == n_pass );
    }
    else if( comp == CutType::OR ) {
        pass_cuts =(  n_pass > 0 );
    }
    else {
        std::cout << "CutConfig::PassFloat - ERROR : Did not understand Comparison type!" << std::endl;
    }


    if( GetIsInverted() ) {
      pass_cuts = !pass_cuts;
    }

    return pass_cuts;
}
        
bool CutConfig::PassInt(const std::string &name, const int cutval ) const {

    int n_pass = 0;
    int n_cuts = 0;
    BOOST_FOREACH( const Cut & cut, GetCuts() ) {
        n_cuts++;

        if( cut.type != CutType::INT ) {
          std::cout << "CutConfig::PassInt - ERROR : Int cut requested for cut " << name << " but cut was not configured as an int" << std::endl;
          continue;
        }

        if( cut.op == CutType::GREATER_THAN_OR_EQUAL_TO ) {
          if( (cutval >= cut.val_int) ) n_pass++;
        }
        if( cut.op == CutType::LESS_THAN_OR_EQUAL_TO ) {
          if( (cutval <= cut.val_int) ) n_pass++;
        }
        if( cut.op == CutType::GREATER_THAN ) {
          if( (cutval > cut.val_int) ) n_pass++;
        }
        if( cut.op == CutType::LESS_THAN ) {
          if( (cutval < cut.val_int) ) n_pass++;
        }
        if( cut.op == CutType::EQUAL_TO ) {
          if( (cutval == cut.val_int) ) n_pass++;
        }
        if( cut.op == CutType::NOT_EQUAL_TO ) {
          if( (cutval != cut.val_int) ) n_pass++;
        }
    }

    // Handle AND vs OR
    bool pass_cuts = false;
    if( comp == CutType::AND ) {
        pass_cuts = ( n_cuts == n_pass );
    }
    else if( comp == CutType::OR ) {
        pass_cuts =(  n_pass > 0 );
    }
    else {
        std::cout << "CutConfig::PassInt - ERROR : Did not understand Comparison type!" << std::endl;
    }

    if( GetIsInverted() ) {
      pass_cuts = !pass_cuts;
    }

    return pass_cuts;
}
        
bool CutConfig::PassBool(const std::string &name, const bool cutval ) const {

    int n_pass = 0;
    int n_cuts = 0;
    BOOST_FOREACH( const Cut & cut, GetCuts() ) {

        n_cuts++;

        if( cut.type != CutType::BOOL ) {
          std::cout << "CutConfig::PassBool - ERROR : Boolean cut requested for cut " << name << " but cut was not configured as a bool " << std::endl;
          continue;
        }

        if( cut.op == CutType::EQUAL_TO ) {
            if( (cutval == cut.val_bool) ) n_pass++;
        }
        if( cut.op == CutType::NOT_EQUAL_TO ) {
            if( (cutval != cut.val_bool) ) n_pass++;
        }
    }

    // Handle AND vs OR
    bool pass_cuts = false;
    if( comp == CutType::AND ) {
        pass_cuts = ( n_cuts == n_pass );
    }
    else if( comp == CutType::OR ) {
        pass_cuts =(  n_pass > 0 );
    }
    else {
        std::cout << "CutConfig::PassBool - ERROR : Did not understand Comparison type!" << std::endl;
    }

    if( GetIsInverted() ) {
        pass_cuts = !pass_cuts;
    }

    return pass_cuts;
}
        

bool CutConfig::PassAnyIntVector(const std::string &name, const std::vector<int> &cutval ) const {

    int n_pass = 0;
    int n_cuts = 0;
    BOOST_FOREACH( const Cut & cut, GetCuts() ) {
        n_cuts++;

        if( cut.type != CutType::INT ) {
          std::cout << "CutConfig::PassAnyIntVector - ERROR : Int cut requested for cut " << name << " but cut was not configured as an int" << std::endl;
          continue;
        }

        if( cut.op == CutType::EQUAL_TO ) {
          if( std::find(cutval.begin(), cutval.end(), cut.val_int) != cutval.end() ) n_pass++;
        }
        else {
            std::cout << "CutConfig::PassAnyIntVector - ERROR : Expect EQUAL_TO requirement for finding int in vector" << std::endl;
            continue;
        }
    }

    // Handle AND vs OR
    bool pass_cuts = false;
    if( comp == CutType::AND ) {
        pass_cuts = ( n_cuts == n_pass );
    }
    else if( comp == CutType::OR ) {
        pass_cuts =(  n_pass > 0 );
    }
    else {
        std::cout << "CutConfig::PassInt - ERROR : Did not understand Comparison type!" << std::endl;
    }

    if( GetIsInverted() ) {
      pass_cuts = !pass_cuts;
    }

    return pass_cuts;
}
        

bool ModuleConfig::PassInt( const std::string & cutname, const int cutval )
{

    if( HasCut( cutname ) ) {
        const CutConfig & cut_conf = GetCut( cutname );
        bool result = cut_conf.PassInt( cutname, cutval );

        if( cutflows.size() ) { // only assume 1 cutflow for now
            cutflows[0].AddCutDecisionInt( cutname, result, cutval );
        }

        return result;
    }
    else {
        //if the cut doesn't exist then pass
        return true;
    }

}

bool ModuleConfig::PassBool( const std::string & cutname, const bool cutval )
{

    if( HasCut( cutname ) ) {
        const CutConfig & cut_conf = GetCut( cutname );
        bool result = cut_conf.PassBool( cutname, cutval );

        if( cutflows.size() ) { // only assume 1 cutflow for now
            cutflows[0].AddCutDecisionBool( cutname, result, cutval );
        }

        return result;
    }
    else {
        //if the cut doesn't exist then pass
        return true;
    }

}

bool ModuleConfig::PassFloat( const std::string & cutname, const float cutval )
{

    if( HasCut( cutname ) ) {
        const CutConfig & cut_conf = GetCut( cutname );
        bool result = cut_conf.PassFloat( cutname, cutval );

        if( cutflows.size() ) { // only assume 1 cutflow for now
            cutflows[0].AddCutDecisionFloat( cutname, result, cutval );
        }

        return result;
    }
    else {
        //if the cut doesn't exist then pass
        return true;
    }
      

}

bool ModuleConfig::PassAnyIntVector( const std::string & cutname, const std::vector<int> &cutval )
{

    if( HasCut( cutname ) ) {
        const CutConfig & cut_conf = GetCut( cutname );
        bool result = cut_conf.PassAnyIntVector( cutname, cutval );

        if( cutflows.size() ) { 
            std::cout << "No cutflows for PassAnyIntVector" << std::endl;
        }

        return result;
    }
    else {
        //if the cut doesn't exist then pass
        return true;
    }

}


bool ModuleConfig::HasCut( const std::string &name ) const {

    bool found=false;
    BOOST_FOREACH( const CutConfig & cut_conf, configs ) {
        if( name == cut_conf.GetName() ) {
            found = true;
            break;
        }
    }

    return found;
}

const CutConfig & ModuleConfig::GetCut( const std::string & name ) const {

    BOOST_FOREACH( const CutConfig & cut_conf, configs ) {
        if( cut_conf.GetName() == name ) {
            return cut_conf;
        }
    }

    // This will throw an exception.  Its a bit hacked, but its the desired functionality
    std::cout << "Requested cut does not exist!  To avoid this error, use HasCut to check for its existance before calling GetCut" << std::endl;
    assert (false);
    return configs[configs.size()];

}

void ModuleConfig::AddCutFlow( const std::string & name ) {

    if( cutflows.size() == 0 ) {
        cutflows.push_back( CutFlowModule( name ) );
    }

}

void ModuleConfig::AddCut( CutConfig config ) {

    configs.push_back(config);

}

void ModuleConfig::AddHist( const std::string &histname, int nbin, float xmin, float xmax ) {

    if( cutflows.size() ) {
        std::string befname = histname + "_before";
        std::string aftname = histname + "_after";

        cutflows[0].createHist(GetName(), befname, nbin, xmin, xmax );
        cutflows[0].createHist(GetName(), aftname, nbin, xmin, xmax );
        //cutflows[0].getHists()[befname] = new TH1F( (GetName() + ":" + befname).c_str(), befname.c_str(), nbin, xmin, xmax );
        //cutflows[0].getHists()[aftname] = new TH1F( (GetName() + ":" + aftname).c_str(), aftname.c_str(), nbin, xmin, xmax );
    }

}

void ModuleConfig::PrintCutFlows() const {

    BOOST_FOREACH( const CutFlowModule & cf, cutflows ) {
        cf.Print();
    }
}

void ModuleConfig::WriteCutFlowHists( TFile * cutflowfile ) const {

    BOOST_FOREACH( const CutFlowModule & cf, cutflows ) {
        cf.WriteCutFlowHist( cutflowfile );
    }
}

CutFlowModule::CutFlowModule( const std::string & _name ) :
  total(0),
  name(_name)
{
}

CutFlowModule::~CutFlowModule() {
}

void CutFlowModule::WriteCutFlowHist( TDirectory * dir) const
{

    unsigned nbins = order.size() + 1;
    TH1F* cuthist = new TH1F( ( GetName() + ":cuthist").c_str(), "cuthist", nbins, 0, nbins );
    cuthist->GetXaxis()->SetBinLabel( 1, "Total" );
    cuthist->SetBinContent( 1, total );
    for( unsigned binnum = 0; binnum < order.size(); ++binnum ) {
        cuthist->GetXaxis()->SetBinLabel( binnum+2, order[binnum].c_str() );
        cuthist->SetBinContent(binnum+2, counts.at(order[binnum]));
    }
    //cuthist->SetDirectory( thisdir );

    cuthist->Write();

}

    

//CutFlowModule::CutFlowModule(const CutFlowModule & copy ) {
//
//    for( std::map<std::string, TH1F*>::const_iterator itr = copy.getHists().begin(); itr != getHists().end(); ++itr ) {
//        getHists().insert( std::make_pair( itr->first, new TH1F( *(itr->second) ) ) );
//    }
//
//    copy.SetOrder(order);
//    copy.SetCounts(counts);
//    copy.SetTotal(total);
//    copy.SetName(name);
//}

void CutFlowModule::createHist( const std::string &basename, const std::string &histname, 
                                int nbin, float xmin, float xmax) 
{

    hists.insert( std::make_pair(histname, 
                                 TH1F( ( basename + ":" + histname).c_str(), histname.c_str(), 
                                 nbin, xmin, xmax ) ) );
}

void CutFlowModule::AddCutDecision( const std::string & cutname, bool pass, float weight ) {

    std::map<std::string, float>::iterator fitr = counts.find( cutname );
    if( fitr == counts.end() ) {
        order.push_back(cutname);
        counts[cutname] = 0.0;
        counts[cutname] += pass*weight;
    }
    else {
        fitr->second += pass*weight;
    }
    // for the first cut, also add a total that
    // doesn't care about the cut decision
    if( std::find(order.begin(), order.end(), cutname ) == order.begin() ) {
        total += weight;
    }

}

void CutFlowModule::AddCutDecisionFloat( const std::string & cutname, bool pass, float val, float weight ) {

    AddCutDecision( cutname, pass, weight );
    // fill hists
    if( hists.size() ) {
        std::string befname = cutname + "_before";
        std::string aftname = cutname + "_after";
        std::map<std::string, TH1F>::iterator hitr = hists.find(befname);
        if( hitr != hists.end() ) {
            hitr->second.Fill( val );
        }
        // fill the after hist only when the cut passes
        if( pass ) {
            std::map<std::string, TH1F>::iterator hitr = hists.find(aftname);
            if( hitr != hists.end() ) {
                hitr->second.Fill( val );
            }
        }
    }
}

void CutFlowModule::AddCutDecisionInt( const std::string & cutname, bool pass, int val, float weight ) {

    AddCutDecision( cutname, pass, weight );
    // fill hists
    if( hists.size() ) {
        std::string befname = cutname + "_before";
        std::string aftname = cutname + "_after";
        std::map<std::string, TH1F>::iterator hitr = hists.find(befname);
        if( hitr != hists.end() ) {
            hitr->second.Fill( val );
        }
        // fill the after hist only when the cut passes
        if( pass ) {
            std::map<std::string, TH1F>::iterator hitr = hists.find(aftname);
            if( hitr != hists.end() ) {
                hitr->second.Fill( val );
            }
        }
    }
}

void CutFlowModule::AddCutDecisionBool( const std::string & cutname, bool pass, bool val, float weight ) {

    AddCutDecision( cutname, pass, weight );
    // fill hists
    if( hists.size() ) {
        std::string befname = cutname + "_before";
        std::string aftname = cutname + "_after";
        std::map<std::string, TH1F>::iterator hitr = hists.find(befname);
        if( hitr != hists.end() ) {
            hitr->second.Fill( val );
        }
        // fill the after hist only when the cut passes
        if( pass ) {
            std::map<std::string, TH1F>::iterator hitr = hists.find(aftname);
            if( hitr != hists.end() ) {
                hitr->second.Fill( val );
            }
        }
    }
}
void CutFlowModule::Print() const {

    // get max width
    int max_width = 0;
    BOOST_FOREACH( const std::string & name, order ) {
        int width = name.size();
        if( width > max_width ) {
          max_width = width;
        }
    }

    // make it a bit bigger
    max_width += 5;
    int line_width = 50;
    if( max_width > line_width ) {
      line_width = max_width + 5;
    }

    std::cout << std::string(line_width, '-') << std::endl;
    std::cout << "Cut flow : " << name << std::endl;
    std::cout << std::string(line_width, '-') << std::endl;
    std::cout << std::setw(max_width) << std::setfill(' ') << std::setiosflags(std::ios::left) << std::fixed << "Total" << " : " << total << std::endl;
    BOOST_FOREACH( const std::string & name, order ) {
        
        std::cout << std::setw(max_width) << std::setfill(' ') << std::setiosflags(std::ios::left) << name << " : " << counts.find(name)->second << std::endl;
    }
    std::cout << std::string(line_width, '-') << std::endl;

}

AnaConfig ParseConfig( const std::string & fname, CmdOptions & options ) {

    AnaConfig ana_config;

    std::ifstream file(fname.c_str());
    std::string line;
    if( file.is_open() ) {
        
        bool read_modules = false;
        while( getline(file, line) ) {

            // __Modules__ line indicates the beginning of modules
            // when this line is encountered, don't parse the line, just
            // indicate that the following lines should be modules

            std::string valid_line = line;

            // remove comments
            std::string::size_type comment_pos = line.find("#");
            if( comment_pos != std::string::npos ) {
                valid_line = line.substr( 0, comment_pos );
            }

            if( valid_line.empty() ) {
                continue;
            }

            if( valid_line.find("__Modules__") != std::string::npos ) {
                read_modules = true;
            }
            else if( read_modules ) { // when the "Modules" line is found this is set to true
                ReadModuleLine( valid_line, ana_config );
            }
            else {
                // by default read header lines
                ReadHeaderLine( valid_line, options );
            }

        }
    }
    else {
        std::cout << "ERROR - Failed to open config file " << fname << std::endl;
    }
    return ana_config;
}


void ReadModuleLine( const std::string & line, AnaConfig & config ) {

    // Split by : character.  The 0th entry is the
    // module name, 1st entry is the configuration
    std::vector<std::string> module_split = Tokenize( line, ":" );
    if( module_split.size() < 2 ) {
        std::cout << "ParseConfig - ERROR : config file entry does not contain a \":\" "
                  << "character.  Please check the config file" << std::endl;
        return;
    }
    std::string module_name = module_split[0];
    boost::algorithm::trim(module_name);
    ModuleConfig this_module(module_name);

    // get the rest of the string ( takes care if multiple ':' are present ) 
    std::string module_config = boost::algorithm::join( std::vector<std::string>( module_split.begin()+1, module_split.end() ), ":" );

    // Split the module config by ';' to separate cut entries
    std::vector<std::string> cut_split = Tokenize( module_config, ";" );

    std::vector<CutConfig> module_cuts;
    BOOST_FOREACH( std::string & cut, cut_split ) {
        ReadCut( cut, this_module );
    }

    config.AddModule( this_module );

}

void ReadCut( std::string &cut, ModuleConfig & module ) {

    if( cut.find_first_not_of(' ') == std::string::npos ) return; //check if cut is only whitespace 

    //trim whitespae
    boost::algorithm::trim(cut);

    // if value is do_cutflow just call AddCutFlow
    if( cut.find("do_cutflow") != std::string::npos ) {
        module.AddCutFlow( module.GetName() );
        return;
    }
    // if value is hist parse the histogram parameters
    else if( cut.find("hist") == 0 ) {
        ParseHistPars( cut, module );
        return;
    }
    else if( cut.find("init_") == 0 ) {
        ParseDataEntry( cut, module );
        return;
    }
    else {
        // otherwise this is a normal cut
        module.AddCut( CutConfig( cut ) );
    }

}

void ParseHistPars( const std::string & cut, ModuleConfig & module ) {

    // example hist entry : hist [cut_d0_barrel,100,-1.000000,1.000000]
    // first strip off the []
    std::string hist_entries = cut.substr( cut.find("[")+1, cut.find("]")-1 );
    std::vector<std::string> hist_split = Tokenize( hist_entries, "," );
    if( hist_split.size() != 4 ) {
        std::cout << "ParseHistPars - ERROR : cannot parse histogram.  Expect 4 entries : name, nbins, xmin, xmax " << hist_entries << std::endl;
        return;
    }
    std::string name = hist_split[0];
    int nbins;
    float xmin;
    float xmax;
    // trim the entries and convert to integer/float
    boost::algorithm::trim(hist_split[1]);
    boost::algorithm::trim(hist_split[2]);
    boost::algorithm::trim(hist_split[3]);
    std::stringstream nbinstr(hist_split[1]);
    std::stringstream xminstr(hist_split[2]);
    std::stringstream xmaxstr(hist_split[3]);
    nbinstr >> nbins;
    xminstr >> xmin;
    xmaxstr >> xmax;
    // add histogram
    module.AddHist( name, nbins, xmin, xmax );

}
void ParseDataEntry( const std::string & cut_str, ModuleConfig & module ) {

    // cut string should have the format name [val]
    // find the location of the brackets
    std::size_t posbeg = cut_str.find("[");
    std::size_t posend = cut_str.find("]");

    // throw an error if the formatting is incorrect
    if( posbeg == std::string::npos || posend == std::string::npos ) {

      std::cout << "CutConfig - ERROR : cut string should have the format name [val]"
                << "Instead got " << cut_str << std::endl;
      return;
    }

    // get the cut name values
    std::string name    = cut_str.substr( 0, posbeg);
    std::string val_str = cut_str.substr( posbeg+1, posend-posbeg-1 );

    // remove data_ from the name
    name = name.substr( 5 );

    // remove whitespace from the name
    boost::algorithm::trim(name);
    boost::algorithm::trim(val_str);

    module.AddInitData( name, val_str );

}

void ReadHeaderLine( const std::string & line, CmdOptions & options ) {

    // reading header information
    // every header entry should have the format
    // name : value
    std::vector<std::string> header_split = Tokenize( line, ":" );
    if( header_split.size() < 2 ) {
        std::cout << "ReadHeaderLine - ERROR : cannot parse header.  Expect a format of name : value " << std::endl;
        return;
    }

    std::string header_key = header_split[0];
    // the remaining part of the string has the header value
    // takes care of the case when multiple ':' are present in the line
    std::string header_val = boost::algorithm::join( std::vector<std::string>( header_split.begin()+1, header_split.end() ) , ":");
    boost::algorithm::trim(header_val);

    if( header_key.find("files") != std::string::npos ) {// read files
        ParseFiles( header_val, options );
    }
    else if( header_key.find("treeName") != std::string::npos ) {
        options.treeName = header_val;
        boost::algorithm::trim(options.treeName);
    }
    else if( header_key.find("outputDir") != std::string::npos ) {
        options.outputDir = header_val;
        boost::algorithm::trim(options.outputDir);
    }
    else if( header_key.find("outputFile") != std::string::npos ) {
        options.outputFile = header_val;
        boost::algorithm::trim(options.outputFile);
    }
    else if( header_key.find("storagePath") != std::string::npos ) {
        options.storagePath = header_val;
        boost::algorithm::trim(options.storagePath);
        options.transferToStorage = true;
    }
    else if( header_key.find("nevt") != std::string::npos ) {
        std::stringstream ss(header_val);
        ss >> options.nevt;
    }
    else if( header_key.find("nPrint") != std::string::npos ) {
        std::stringstream ss(header_val);
        ss >> options.nPrint;
    }
    else if( header_key.find("sample") != std::string::npos ) {
        std::stringstream ss(header_val);
        options.sample = header_val;
        boost::algorithm::trim(options.sample);
    }
    else if( header_key.find("disableOutputTree") != std::string::npos ) {
        if( header_val.find("true") != std::string::npos ) {
            options.disableOutputTree = true;
        }
    }
    else if( header_key.find("copyInputFiles") != std::string::npos ) {
        if( header_val.find("true") != std::string::npos ) {
            options.copyInputFiles = true;
        }
    }

}

void ParseFiles( const std::string & files_val, CmdOptions & options ) {

    // this is a bit complicated
    // jobs are associated to one or more input files
    // the jobs are labeled by a number ID and are associated
    // with an event range
    // The format is like this
    // [file1,file2,file3][0:(0-10),1:(10-20)];[file4,file5][2:(0-10)]
    //
    // first split into individual [files][jobs] entries with a ;
    std::vector<std::string> file_map_entries = Tokenize( files_val, ";");
    BOOST_FOREACH( const std::string & file_map, file_map_entries ) {
     
        // First check that the string begins with [ and ends with ] 
        if( !(file_map[0] == '[' and file_map[file_map.size()-1] == ']' ) ) {
            std::cout << "ParseFiles - ERROR : String should begin with [ and end with ]" << std::endl;
            continue;
        }

        // remove the '[' at the beginning and ']' at the end
        std::string file_map_mod = file_map.substr(1, file_map.size() - 2 );
        
        // split by "][" to get individual entries.  There should be 2
        std::vector<std::string> file_map_split = Tokenize( file_map_mod, "][" );
        if( !(file_map_split.size() == 2) ) {
            std::cout << "ParseFiles - ERROR : File entry should have a list of files and a list of events" << std::endl;
            continue;
        }
        // entries are printed as python tuples separated by a comma
        std::vector<std::string> job_list = Tokenize( file_map_split[1], "," );

        // vector to collect the job info.  Map the job id to an event range
        std::vector< std::pair< int, std::pair< int, int > > > out_job_list;
        int jobidx = -1;
        BOOST_FOREACH( const std::string event_vals, job_list ) {
            jobidx++;
            // the entry is like 0:(0-500) 
            std::vector<std::string> job_evtrange = Tokenize( event_vals, ":" );
            std::string jobidstr = job_evtrange[0];
            // trim jobid and convert to integer
            boost::algorithm::trim(jobidstr);
            std::stringstream jobidss(jobidstr);
            int jobid;
            jobidss >> jobid;
            
            // Strip the ( ) from the entry
            std::string event_vals_mod = job_evtrange[1].substr( 1, job_evtrange[1].size() - 2 );

            // Split the values by "-" there should be two entries
            std::vector<std::string> vals = Tokenize( event_vals_mod, "-");
            if( !(vals.size() == 2) ) {
                std::cout << "ParseFiles - ERROR : Events entry should have a size 2 tuple of integers" << std::endl;
                continue;
            }
            // trim and convert to integers
            boost::algorithm::trim(vals[0]);
            boost::algorithm::trim(vals[1]);
            std::stringstream minstr(vals[0]);
            std::stringstream maxstr(vals[1]);
            int minval;
            int maxval;
            minstr >> minval ;
            maxstr >> maxval ;

            // store the configuration
            std::pair<int, int> evt_range( minval, maxval );
            out_job_list.push_back( std::make_pair( jobid, evt_range ) );
        }

        std::cout << "Add files " << file_map_split[0] << std::endl;
        std::vector<std::string> files = Tokenize(file_map_split[0], ",");
        FileEntry entry;
        entry.files = files;
        entry.jobs = out_job_list;
        options.files.push_back( entry );
    }
}

CmdOptions::CmdOptions() : nevt(-1), nPrint(10000), transferToStorage(false), disableOutputTree(false), copyInputFiles(false)
{
}

CmdOptions ParseOptions( int argc, char **argv ) 
{

    CmdOptions options;
    const struct option longopts[] =
         {
         {"conf_file", required_argument,        0, 'c'},
         {0,0,0,0},
    };

    int iarg=0;

    //turn off getopt error message
    
    while(iarg != -1) {
      iarg = getopt_long(argc, argv, "c:", longopts, 0);
      switch (iarg) {
        case 'c' : 
          {
          options.config_file = optarg;
          break;
          }
      }
    }

    return options;
}

std::vector<std::string> Tokenize(const std::string & text, const std::string &in_sep ) {

    std::vector<std::string> out_tokens;
    if( text == "" ) {
        std::cout << "Tokenize : WARNING - string to split is empty" << std::endl;
        return out_tokens;
    }

    boost::char_separator<char> sep(in_sep.c_str());
    boost::tokenizer< boost::char_separator<char> > tokens(text, sep);
    BOOST_FOREACH (const std::string& t, tokens) {
        out_tokens.push_back(t);
    }

    return out_tokens;
}
        
    
