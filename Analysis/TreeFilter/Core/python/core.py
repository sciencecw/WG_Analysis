import os
import re
import sys
import imp
import math
import ROOT
import array
import copy
import inspect
import subprocess
import multiprocessing
import logging
import time
import datetime
import filecmp
import uuid
from argparse import ArgumentParser
import eos_utilities as eosutil

def ParseArgs() :

    parser = ArgumentParser(description='')
    
    #-----------------------------
    # input files
    #-----------------------------
    parser.add_argument('--files', dest='files', default=None, help='list of input files (comma separated).  If not provided, --filesDir must be provided')
    
    parser.add_argument('--filesDir', dest='filesDir', default=None, help='Directory to search for files.  If not provided, --files must be provided')

    parser.add_argument('--noInputFiles', dest='noInputFiles', default=False, action='store_true', help='Do not expect any files')

    #-----------------------------
    # Output locations
    #-----------------------------
    
    parser.add_argument('--outputDir', dest='outputDir', default=None, help='<Required> Output sample path.  If nproc is set to > 1 a job directory is appended')
    
    parser.add_argument('--outputFile', dest='outputFile', default=None, help='Name of output file.  If none given it will assume the name of one if the input files')
    
    parser.add_argument('--storagePath', dest='storagePath', default=None, help='After the output file is written, transfer it to storage.  Use --outputDir to give the local output from where the file is transferred.  Local files are removed after the copy')

    parser.add_argument('--confFileName', dest='confFileName', default='analysis_config.txt', help='Name of the configuration file.  Default : analysis_config.txt')

    parser.add_argument('--exeName', dest='exeName', default='RunAnalysis.exe', help='Name of the executable file.  Use to avoid overwriting an executable that is in use')

    #-----------------------------
    # Important additional input info
    #-----------------------------
     
    parser.add_argument('--treeName', dest='treeName', default=None, help='<Required> name of tree to process')

    parser.add_argument('--fileKey', dest='fileKey', default='.root', help='key to match root files')
    
    parser.add_argument('--module', dest='module', default=None, help='<Required> name of the module to import')

    #-----------------------------
    # Logging
    #-----------------------------

    parser.add_argument('--loglevel', dest='loglevel', default='INFO', help='Log level, DEBUG, INFO, WARNING, ERROR, CRITICAL. Default is INFO')

    parser.add_argument('--nPrint', dest='nPrint', default=10000, type=int, help='Frequency at which the event number is printed')

    #-----------------------------
    # Configure job splitting
    #-----------------------------
    
    parser.add_argument('--nsplit', dest='nsplit', type=int, default=0, help='Split into nsplit subjobs.  If --nFilesPerJob is set instead, nsplit will be computed.  Jobs are first split by the number of input files.  If there are fewer input files than jobs requested the input files are further split.  Also set nproc > 1 to use multiprocessing.  However multiprocessing is only used by file')

    parser.add_argument('--nJobs', dest='nJobs', type=int, default=0, help='Equivalent to --nsplit')
    
    parser.add_argument('--nFilesPerJob', dest='nFilesPerJob', type=int, default=0, help='Number of files to run in each job.  If --nsplit is defined, then this is ignored')

    parser.add_argument('--totalEvents', dest='totalEvents', type=int, default=0, help='Total number of events to run over.  Only use to run over fewer than all events')

    parser.add_argument('--nproc', dest='nproc', type=int, default=1, help='Number of processors to use.  If set to > 1, use multiprocessing')

    parser.add_argument('--batch', dest='batch', action='store_true', default=False, help='Submit jobs to the batch')

    parser.add_argument('--condor', dest='condor', action='store_true', default=False, help='Submit jobs to condor (UMD)')

    parser.add_argument('--resubmit', dest='resubmit', action='store_true', default=False, help='Only submit jobs whose output does not exist')
    
    #-----------------------------
    # Filter flags
    #-----------------------------

    parser.add_argument('--enableRemoveFilter', dest='enableRemoveFilter', default=False, action='store_true', help='Do not write out branches in the remove filter in the analysis module')
    
    parser.add_argument('--enableKeepFilter', dest='enableKeepFilter', default=False, action='store_true', help='Only write out branches in the remove filter in the analysis module.  Can specify branches to remove within the keep filter by enabling the remove filter')

    parser.add_argument('--removeFilterSelection', dest='removeFilterSelection', default=None, help='Pass this tag to the analysis module to select branches to remove')

    parser.add_argument('--keepFilterSelection', dest='keepFilterSelection', default=None, help='Pass this tag to the analysis module to select branches to keep')

    #-----------------------------
    # Control flags
    #-----------------------------

    parser.add_argument('--noRun', dest='noRun', default=False, action='store_true', help='Just generate and compile, do not run the merging')
    
    parser.add_argument('--noCompile', dest='noCompile', default=False, action='store_true', help='Just run the merging.  Do not recompile')
    
    parser.add_argument('--noCompileWithCheck', dest='noCompileWithCheck', default=False, action='store_true', help='Recompile only if tree structure has changed')
    
    parser.add_argument('--debugCode', dest='debugCode', default=False, action='store_true', help='Place debugging statements in the written code')
    
    parser.add_argument('--writeExpertCode', dest='writeExpertCode', default=False, action='store_true', help='Write additional functions for expert debugging')
    
    parser.add_argument('--disableOutputTree', dest='disableOutputTree', default=False, action='store_true', help='Do not write events to the output tree')

    parser.add_argument('--copyInputFiles', dest='copyInputFiles', default=False, action='store_true', help='Transfer files to local directory. Only works with batch.')
    
    #-----------------------------
    # Occasionally used
    #-----------------------------

    parser.add_argument('--sample', dest='sample', default=None, help='Name of sample.  May be used in cases where the sample name must be known to the c++ code')
    
    parser.add_argument('--moduleArgs', dest='moduleArgs', default=None, help='Arguments to pass to module.  Should be in the form of a dictionary')

    parser.add_argument('--write_file_list', dest='write_file_list', default=False, action='store_true', help='Write a text file containing files found on eos.  use --read_file_list to use this file in the future')

    parser.add_argument('--read_file_list', dest='read_file_list', default=False, action='store_true', help='Read the file list instaed of looking on eos.  Must have used --write_file_list to create the file' )
    
    return parser.parse_args()

class JobConfig :

    def __init__(self) :

        self.jobid        = None
        self.storage_path = None
        self.output_dir   = None
        self.run_command  = None


def config_and_run( options, package_name ) :

    assert options.noInputFiles or options.files is not None or options.filesDir is not None , 'Must provide a file list via --files or a search directory via --filesDir'
    assert options.outputDir is not None, 'Must provide an output directory via --outputDir'
    assert options.noInputFiles or options.treeName is not None, 'Must provide a tree name via --treeName'
    assert options.module is not None, 'Must provide a module via --module'
    
    if options.batch and options.noCompileWithCheck :
        assert False, "Running noCompileWithCheck with batch mode can result in an executable not being created for a batch job!"

    if options.copyInputFiles and not ( options.batch or options.condor ) :
        print "Can only copy input files in batch mode"
        options.copyInputFiles = False

    if options.loglevel.upper() not in ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'] :
        print 'Loglevel must be one of DEBUG, INFO, WARNING, ERROR, CRITICAL.  Default to INFO'
        options.loglevel = 'INFO'

    #initialize the logger 
    logging.basicConfig(format='%(levelname)s:\t%(message)s', level=getattr(logging, options.loglevel.upper(), None) )

    workarea = os.getenv('WorkArea')
    assert workarea is not None, 'Did not locate WorkArea environment variable.  Please set it to the root of the package.'

    input_files = []
    if options.files is not None :
        if options.files.count('root://') > 0 :
            split_files = options.files.split(',')
            for file in split_files :
                dir_base = '/'.join(file.split('/')[0:-1])
                filename = file.split('/')[-1]
                input_files += collect_input_files_eos( dir_base, filename)
        else :
            input_files = options.files.split(',')
            input_files = filter( lambda s: os.path.isfile( s ), input_files )
    elif options.filesDir is not None :
        for eachdir in options.filesDir.split(',') :
            input_files += collect_input_files( eachdir, options.fileKey, write_file_list=options.write_file_list, read_file_list=options.read_file_list )

    # remove any blank entries
    while '' in input_files :
        input_files.remove('')

    if not options.read_file_list :
        input_files = check_and_filter_input_files( input_files, options.treeName )

    if not options.noInputFiles and not input_files :
        print 'Did not locate any input files with the path and file name given!  Check the inputs.'
        if options.files is not None :
            print options.files
        else :
            print options.filesDir
        return False

    # if output file name is not given grab the name from one of the input files
    if options.outputFile is None :
        if input_files :
            options.outputFile = input_files[0].split('/')[-1]

    logging.info('Get Branch mapping from file %s'%input_files[0])
    branches = get_branch_mapping(input_files[0], options.treeName )

    ImportedModule = import_module( options.module )

    branches_to_keep = get_keep_branches( ImportedModule, branches, 
                                               options.enableKeepFilter, 
                                               options.keepFilterSelection,
                                               options.enableRemoveFilter,
                                               options.removeFilterSelection)

    if options.enableKeepFilter :
        print 'Will keep %d branches from output file : ' %len(branches_to_keep)
        print '\n'.join(branches_to_keep)
    elif options.enableRemoveFilter :
        print 'Will remove %d branches from output file : ' %( len(branches) - len(branches_to_keep))
        print '\n'.join( list( set( [ br['name'] for br in branches ] ) - set( branches_to_keep ) ) )

    # -------------------------------
    # gather module arguments
    # put outputDir into the arguments
    # so it can be picked up from the
    # module
    # -------------------------------
    modargs = {}
    if options.moduleArgs is not None :
        modargs = eval( options.moduleArgs )
    if 'outputDir' not in modargs :
        modargs['outputDir'] = options.outputDir

    # -------------------------------
    # populate the alg_list from 
    # the imported module
    # -------------------------------
    alg_list = []
    try :
        ImportedModule.config_analysis(alg_list, modargs)
    except TypeError, e : 
        logging.warning('********************************')
        logging.warning('Could not call config_analysis with two arguments')
        logging.warning('To maintain compatibility with the old method of using a single argument')
        logging.warning('The function will be called in this way.')
        logging.warning('Just in case the exception is,')
        logging.warning(e)
        logging.warning('********************************')

        ImportedModule.config_analysis(alg_list)

    user_include_file_name = '%s/TreeFilter/%s/include/RunAnalysis.h' %( workarea, package_name )
    out_branches = get_out_branches( user_include_file_name, branches, alg_list )

    # get the executable name.  If .exe does not
    # appear as the extension, add it
    if re.match( '.*?\.exe', options.exeName ) is None :
        options.exeName = options.exeName + '.exe'

    if options.noCompileWithCheck : 

        #-------------------------------------------
        # use nCompileWithCheck when you want to 
        # use the existing executable (as if using --noCompile)
        # except when the tree content has changed
        # The check is done by building new compilation files and
        # doing a diff between the new and old files
        # if a diff is found a compilation is triggered
        #-------------------------------------------

        old_brdef_file_name = '%s/TreeFilter/%s/include/BranchDefs.h' %( workarea, package_name )
        old_header_file_name = '%s/TreeFilter/%s/include/BranchInit.h' %( workarea, package_name )
        old_source_file_name = '%s/TreeFilter/%s/src/BranchInit.cxx' %(workarea, package_name )
        old_linkdef_file_name = '%s/TreeFilter/%s/include/LinkDef.h' %(workarea, package_name )

        tmpname = str( uuid.uuid4() )

        new_brdef_file_name = '/tmp/%s/BranchDefs.h'  %tmpname
        new_header_file_name = '/tmp/%s/BranchInit.h' %tmpname
        new_source_file_name = '/tmp/%s/BranchInit.cxx' %tmpname
        new_linkdef_file_name = '/tmp/%s/LinkDef.h' %tmpname

        os.mkdir( '/tmp/%s' %tmpname )

        # Write the c++ files having the branch definitions and 
        # SetBranchAddress calls
        # Write all output branches in the 
        # header file so that the code will compile
        # only the keep branches will be saved, however
        #write_header_files(new_brdef_file_name, new_linkdef_file_name, branches,out_branches, [ br['name'] for br in branches ], alg_list )
        write_header_files(new_brdef_file_name, new_linkdef_file_name, branches,out_branches, branches_to_keep, alg_list )

        write_source_file(new_source_file_name, new_header_file_name, branches, out_branches, branches_to_keep, write_expert_code=options.writeExpertCode )

        if not (filecmp.cmp(  old_brdef_file_name, new_brdef_file_name ) and
                filecmp.cmp(  old_header_file_name, new_header_file_name ) and
                filecmp.cmp( old_source_file_name, new_source_file_name )  and
                filecmp.cmp( old_linkdef_file_name, new_linkdef_file_name ) ) :

            print '--------------------------------------'
            print 'Difference found in processing code'
            print 'Will use new code and will compile'
            print '--------------------------------------'

            # If we need to recompile then a different executable should
            # be made so as to not overwrite an existing one
            # to be safe, use a random string
            options.exeName = '.'.join(options.exeName.split('.')[0:-1]) + str( uuid.uuid4() ) + '.exe'

            comp_success = compile_code( alg_list, branches, out_branches, branches_to_keep, workarea, package_name, options.exeName, options.writeExpertCode)

            if not comp_success :
                return


    elif not options.noCompile :

        comp_success = compile_code( alg_list, branches, out_branches, branches_to_keep, workarea, package_name, options.exeName, options.writeExpertCode)
        if not comp_success :
            return


    # Get the path of the executable.  First try the
    # WorkArea environment variable which will give
    # an absolute path.  If that doesn't exist, try
    # Using the path of this script to get an
    # absolute path
    exe_path = '%s/TreeFilter/%s/%s' %(workarea, package_name, options.exeName)


    # --------------------------
    # Handle file splitting.
    # Be default nsplit = 1.  In this case the default 
    # behavior is only one job running over all
    # events
    # --------------------------
    if options.nJobs > 0 and options.nsplit == 0 :
        options.nsplit = options.nJobs

    file_evt_list = get_file_evt_map( input_files, options.nsplit, options.nFilesPerJob, options.totalEvents, options.treeName )

    #if options.nproc > 1 and options.nproc > len(file_evt_list) :
    #    options.nproc = len(file_evt_list)


    logging.info('********************************')
    logging.info('Will run a total of %d processes, %d at a time' %(len(file_evt_list), options.nproc) )
    logging.info('********************************')

    if options.nproc > 1 and file_evt_list > 1: #multiprocessing!

        command_info_orig = generate_multiprocessing_commands( file_evt_list, alg_list, exe_path, vars(options) )
        #command_info_orig = generate_multiprocessing_commands( file_evt_list, alg_list, exe_path, options )

        if options.resubmit :
            command_info = filter_jobs_for_resubmit( command_info_orig, options.storagePath, options.outputDir, options.outputFile, options.treeName )
        else :
            command_info = command_info_orig

        commands = []
        for jobid, info in command_info :
            print info
            commands.append( (jobid, make_exe_command( info[0], info[1], info[2] )) ) 

        # Stop here if not running
        logging.info('********************************')
        for cmd in commands :
            logging.info( 'Executing ' + cmd[1] )
        logging.info('********************************')

        if options.noRun :
            logging.info('Not runnning.  Stop here')
            return

        pool = multiprocessing.Pool(options.nproc)
        pool.map( os.system, [c[1] for c in commands])

    elif options.batch :

        command_info_orig = generate_multiprocessing_commands( file_evt_list, alg_list, exe_path, vars(options) )
        #command_info_orig = generate_multiprocessing_commands( file_evt_list, alg_list, exe_path, options )

        # storagePath could have been passed as 'None'
        # make 'None' None
        if options.storagePath == 'None' :
            options.storagePath = None

        if options.resubmit :
            command_info = filter_jobs_for_resubmit( command_info_orig, options.storagePath, options.outputDir, options.outputFile, options.treeName )
        else :
            command_info = command_info_orig

        commands = []
        for jobid, info in command_info :
            commands.append( (jobid, make_exe_command( info[0], info[1], info[2] )) ) 

        wrappers = []
        for jobid, cmd in commands :
            # write a wrapper script 
            #if options.copyLocal :
            ##    get the event range
            wrapper = '%s/wrapper_%s.sh' %( options.outputDir, jobid) 
            file = open( wrapper, 'w' )
            #file.write( 'source /afs/cern.ch/user/j/jkunkle/.bashrc \n' )
            #file.write( 'source  /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.28/x86_64-slc6-gcc47-opt/root/bin/thisroot.sh \n' )
            file.write( 'export LD_LIBRARY_PATH=/afs/cern.ch/sw/lcg/contrib/gcc/4.6/x86_64-slc6-gcc46-opt/lib64:$LD_LIBRARY_PATH\n' )
            file.write( cmd + '\n' )
            file.close()
            os.system('chmod 777 %s' %wrapper )
            wrappers.append((jobid,wrapper))

        print 'Will submit %d jobs' %len(wrappers)
        for jobid, wrap in wrappers :
            os.system( 'bsub -q 1nh -o %s/out_%s.txt -e %s/err_%s.txt %s ' %(options.outputDir, jobid, options.outputDir, jobid, wrap) )

    elif options.condor :

        command_info_orig = generate_multiprocessing_commands( file_evt_list, alg_list, exe_path, vars(options) )

        # storagePath could have been passed as 'None'
        # make 'None' None
        if options.storagePath == 'None' :
            options.storagePath = None

        if options.resubmit :
            command_info = filter_jobs_for_resubmit( command_info_orig, options.storagePath, options.outputDir, options.outputFile, options.treeName )
        else :
            command_info = command_info_orig

        job_desc_file = create_job_desc_file( command_info, {})

        condor_command = 'condor_submit %s ' % job_desc_file 
        logging.info('********************************')
        logging.info( 'Executing ' + condor_command )
        logging.info('********************************')
        os.system(condor_command)

    

    else :
        if options.resubmit :
            print 'Must implement resubmit for single jobs'
            return

        output_file = '%s/%s' %(options.outputDir, options.outputFile )
        conf_file  = '%s/%s' %(options.outputDir, options.confFileName )
        write_config( alg_list, conf_file, options.treeName, options.outputDir, options.outputFile, file_evt_list, options.storagePath, options.sample, options.disableOutputTree, options.nPrint, options.copyInputFiles ) 
        command = make_exe_command( exe_path, conf_file, options.outputDir  )

        # Stop here if not running
        if options.noRun :
            logging.info('Not runnning.  Stop here')
            return

        if not os.path.isdir( options.outputDir ) and options.outputDir.count('root://') != -1 :
            os.makedirs( options.outputDir )

        logging.info('********************************')
        logging.info( 'Executing ' + command )
        logging.info('********************************')
        os.system(command)

    print 'Output written to %s' %(options.outputDir)


def collect_input_files( filesDir, filekey='.root', write_file_list=False, read_file_list=False ) :
    input_files = []
    if filesDir.count('root://') > 0 :
        input_files = collect_input_files_eos( filesDir, filekey, write_file_list=write_file_list, read_file_list=read_file_list )
    else :
        input_files = collect_input_files_local( filesDir, filekey )

    return input_files

def collect_input_files_local( filesDir, filekey='.root' ) :
    input_files = []
    for top, dirs, files in os.walk(filesDir) :
        for f in files :
            if f.count(filekey) > 0 and os.path.isfile( top+'/'+f ) :
                input_files.append(top+'/'+f)

    return input_files

def collect_input_files_eos( filesDir, filekey='.root', write_file_list=False, read_file_list=False ) :
    
    logging.info('Getting list of input files from eos in %s' %filesDir)
    if read_file_list :
        logging.info('Will read file list ' )
        tmpfile = '/tmp/filelist.txt'
        eosutil.copy_eos_to_local( '%s/filelist.txt' %filesDir, tmpfile )

        input_files = []
        file = open( tmpfile, 'r' )
        for line in file :
            input_files.append( line.rstrip('\n') )

        file.close()

        return input_files 

    else :

        input_files = []
        for top, dirs, files, sizes in eosutil.walk_eos(filesDir) :
            for f in files :
                if f.count(filekey) > 0 :
                    input_files.append(top+'/'+f)

        if write_file_list :
            logging.info('Write files to file list' )
            ofile = open( '/tmp/filelist.txt', 'w' )
            for line in input_files :
                ofile.write(line+'\n')
            ofile.close()
            eosutil.copy_local_to_eos( '/tmp/filelist.txt', '%s/filelist.txt' %( filesDir ) )


        return input_files


#-----------------------------------------------------------
def compile_code( alg_list, branches, out_branches, branches_to_keep, workarea, package_name, exeName, writeExpertCode) :

    lockfile = '.compile.lock'
    lockname = '%s/TreeFilter/%s/%s' %( workarea, package_name, lockfile )
    while os.path.isfile( lockname ) :
        print 'Waiting for other compilation to finish.  If this is not expected, do ctrl-z; rm %s; fg ' %lockfile
        time.sleep(10)

    file=open(lockname, 'w')
    file.write('lock')
    file.close()

    brdef_file_name = '%s/TreeFilter/%s/include/BranchDefs.h' %( workarea, package_name )
    header_file_name = '%s/TreeFilter/%s/include/BranchInit.h' %( workarea, package_name )
    source_file_name = '%s/TreeFilter/%s/src/BranchInit.cxx' %(workarea, package_name )
    linkdef_file_name = '%s/TreeFilter/%s/include/LinkDef.h' %(workarea, package_name )

    # Write the c++ files having the branch definitions and 
    # SetBranchAddress calls
    # Write all output branches in the 
    # header file so that the code will compile
    # only the keep branches will be saved, however
    #write_header_files(brdef_file_name, linkdef_file_name, branches,out_branches, [ br['name'] for br in branches ], alg_list )
    write_header_files(brdef_file_name, linkdef_file_name, branches,out_branches, branches_to_keep, alg_list )

    write_source_file(source_file_name, header_file_name, branches, out_branches, branches_to_keep, write_expert_code=writeExpertCode )

    # compile
    logging.info('********************************')
    logging.info('  Begin Compilation ' )
    logging.info('********************************')
    proc = subprocess.Popen(['make', 'clean'])
    proc.wait()


    proc = subprocess.Popen(['make', 'EXENAME=%s' %exeName])
    retcode = proc.wait()
    
    # abort if non-zero return code
    if retcode :
        logging.error( 'Compilation failed.  Will not run' )
        os.system('rm %s' %lockname)

        return False

    logging.info('********************************')
    logging.info('  Compilation Finished ')
    logging.info('********************************')

    os.system('rm %s' %lockname)

    return True
#-----------------------------------------------------------
def get_branch_mapping_from_trees( trees ) :
    """ get all branches with their types and size """

    if not isinstance(trees, list) :
        trees = [trees]

    files_all_branches = []

    if not trees :
        return []

    for tree in trees :

        all_branches = get_branch_map_from_tree( tree )

        files_all_branches.append(all_branches)

    # now collect the maximum values from all files
    return collect_files_branches( files_all_branches )

#-----------------------------------------------------------
def get_branch_mapping( files, treename ) :
    """ get all branches with their types and size """

    if not isinstance(files, list) :
        files = [files]

    files_all_branches = []

    if not files :
        return []

    for filename in files :
        file = ROOT.TFile.Open( filename )
        tree = file.Get( treename )

        all_branches = get_branch_map_from_tree( tree )

        files_all_branches.append(all_branches)
        file.Close()

    # now collect the maximum values from all files
    return collect_files_branches( files_all_branches )

def collect_files_branches( files_all_branches ) :

    all_branches_max = []
    if len(files_all_branches) == 1 :
        all_branches_max = files_all_branches[0]
    else :

        # start by grabbing the first entry
        all_branches_max = files_all_branches[0]
        # now for the branches associated to each file, loop
        # over the branch names and check if the branch
        # size in this file is larger than the current max branch
        # size.  If larger, update.
        # this might be slow, but it runs only once per job
        for file_branches in files_all_branches[1:] :
            for brinfo in file_branches :
                for idx, brinfomax in enumerate(list(all_branches_max)) :

                    # get branches of the same name
                    if brinfo['name'] != brinfomax['name'] :
                        continue

                    #print 'Original size ',brinfomax['totSize']
                    #print 'New size ',brinfo['totSize']
                    if brinfo['totSize'] > brinfomax['totSize'] :
                        #print 'Update the size of branch %s ' %(brinfo['name'])
                        #print all_branches_max[idx]
                        #print brinfo
                        all_branches_max[idx] = brinfo

                    # if we've matched names we can move on to the next file
                    break

    return all_branches_max

#-----------------------------------------------------------
def get_branch_map_from_tree( tree ) :

    all_branches = []
    for br in tree.GetListOfBranches() :
        
        name  = br.GetName()
        title = br.GetTitle()

        #I don't know why this happens.  Just fix it
        title = title.replace(']]', '][')
 
        lf = br.GetListOfLeaves()[0]
        type  = lf.GetTypeName()
        size  = lf.GetLen()
        is_range = lf.IsRange()
        # match any number of brackets and get them.  
        # Also get the varId at the end.  
        # If there are no brackets the first group is empty
        check = re.match('\w+((?:\[\w+\])*)\s*/(\w)', title) 
        size_leaf = None
        bracketEntry = None
        varId = None
        leafEntry = None
        arrayStr = ''
        size_entries = []
        totSize = lf.GetNdata()  # keep the total size for later
        if check is None :
            leafEntry = title
        else :
            bracketEntry = check.group(1)

            if len(check.groups()) > 1 :
                varId = check.group(2)

            # construct the leafEntry
            leafEntry = name

            # there could be any number of brackets in the variable
            # If there are no brackets the leafEntry is simpy the name
            # plus the varId.  Otherwise do some more work
            if not bracketEntry == '' :
                # get each bracket entry
                each_bracket = bracketEntry.split('][')
                each_bracket = map( lambda x : x.rstrip(']').lstrip('['), each_bracket )

                # if one or more entries is a variable length array we need
                # to get its maximum size.  Otherwise we already have the size
                for  eb in each_bracket :

                    # construct the leafEntry
                    leafEntry += '[%s]' %eb.rstrip('_')

                    is_int = True
                    try :
                        test = int(eb)
                    except ValueError :
                        is_int = False

                    # If its not an int, its another branch.  Get it and find its max size
                    if not is_int :

                        size_leaf_name = eb.rstrip('_')

                        size_br = tree.GetBranch(size_leaf_name)
                        size = size_br.GetListOfLeaves()[0].GetMaximum()+1 #Add one for buffer
                        size_entries.append(size)
                    else :
                        size_entries.append(int(eb))

                # This should be the format of the array in the declaration
                arrayStr = '[' + ']['.join([str(x) for x in size_entries]) + ']'

        if varId is not None :
            leafEntry += '/%s' %varId

        prop_dic = {'name' : name, 'type' : type, 'totSize' : totSize, 'arrayStr' : arrayStr, 'leafEntry' : leafEntry, 'sizeEntries' : size_entries }
        if name in [d['name'] for d in all_branches] :
            continue

        if is_range :
            # any variable that is a range variable needs to come first
            all_branches.insert(0, prop_dic) 
        else :
            all_branches.append(prop_dic)

    return all_branches

#-----------------------------------------------------------
def check_and_filter_input_files( input_files, treename ) :

    logging.info('Run input file checks')
    filtered_files = []
    if not isinstance(input_files, list) :
        input_files = [input_files]

    for filename in input_files :

        pass_filter = True

        file = ROOT.TFile.Open( filename )
        tree = file.Get( treename )
        
        if tree == None :
            pass_filter = False
            print 'Removed file %s that did not contain the input tree' %filename
            print 'Existing top level objects : '
            for obj in file.GetListOfKeys() :
                print obj.GetName()

        if pass_filter :
            filtered_files.append(filename)

        file.Close()


    logging.info('Found %d input files.  %d Files were removed' %(len(filtered_files), len(input_files) - len(filtered_files) ) )
    return filtered_files

def import_module( module ) :

    ImportedModule=None

    ispath = ( module.count('/') > 0 )
    module_path = None
    if ispath :
        module_path = module
    else :
        #get path of this script
        script_path = os.path.realpath(__file__)
        module_path = os.path.dirname(script_path) + '/' + module

    try :
        ImportedModule = imp.load_source(module.split('.')[0], module_path)
    except IOError :
        print 'Could not import module %s' %module_path

    if not hasattr(ImportedModule, 'config_analysis') :
        print 'ERROR - module does not implement a function called config_analysis'
        sys.exit(-1)

    return ImportedModule

def get_keep_branches( module, branches, enable_keep_filter, keep_tag, enable_remove_filter, remove_tag ) :

    remove_filter = []
    keep_filter = []
    if hasattr(module, 'get_remove_filter') :
        try : 
            remove_filter = module.get_remove_filter(remove_tag)
        except TypeError :
            remove_filter = module.get_remove_filter()

        remove_filter = module.get_remove_filter()
    if hasattr(module, 'get_keep_filter') :
        try : 
            keep_filter = module.get_keep_filter(keep_tag)
        except TypeError :
            keep_filter = module.get_keep_filter()

    # by default keep all branches
    all_branches = [ br['name'] for br in branches ]
    branches_to_keep =  list(all_branches)
    if enable_keep_filter :
        #only keep these branches.  Reset the keep list
        branches_to_keep = []
        for kregex in keep_filter :
            matches = [ re.match( kregex, br['name'] ) for br in branches ]
            successful_matches = filter( lambda x : x is not None, matches )
            branches_to_keep += [ x.group(0) for x in successful_matches]

        branches_to_keep = list( set( branches_to_keep ) )

    if enable_remove_filter :
        for rregex in remove_filter :
            matches = [ re.match( rregex, br['name'] ) for br in branches ]
            successful_matches = filter( lambda x : x is not None, matches )
            branches_to_remove = [ x.group(0) for x in successful_matches]
            branches_to_keep = list( set(branches_to_keep) - set(branches_to_remove) )
            
    return branches_to_keep

def get_out_branches( include_name, branches, alg_list ) :
    """ Read include file to get the output branch names """

    ofile = open( include_name , 'r' )
    in_out_namespace = False

    active_defs = get_active_defines( branches, alg_list )

    out_branches = []

    in_defines = []
    remove_strings = ['extern', 'unsigned', 'int', 'bool', 'float', 'string', 'std::string', 'Int_t', 'Bool_t', 'Float_t', 'Double_t', '}', '', '*', ';' ]
    in_mulitline_comment = False
    for line in ofile :

        linestr = line
        # ignore precompiler includes
        if line.count('#') : 
            if in_out_namespace :
                if line.count( '#ifdef' ) :
                    in_defines.append( line.split(' ')[1] )
                if line.count( '#endif' ) :
                    in_defines = in_defines[:-1]
                continue
            else :
                continue

        # ignore single line comments
        if line.count('//') :
            continue

        if line.count('*/') : 
            linestr = line[line.find('*/')+2:]
            in_multiline_comment = False

        if in_mulitline_comment :
            continue

        #ignore multline comments
        if line.count('/*') : 
            in_mulitline_comment = True
            linestr = line[:line.find('/*')]


        # find OUT namespace 
        if linestr.count('namespace') and linestr.count('OUT') :
            in_out_namespace = True
            continue

        if in_out_namespace :

            if linestr.count('}') :
                break

            full_line = linestr.lstrip(' ')

            split_line = linestr.split()
            for rmstr in remove_strings :
                if rmstr in split_line :
                    split_line.remove( rmstr )

            to_rm = []
            for entry in split_line :
                if entry.count('vector') :
                    to_rm.append(entry)

            for rm in to_rm :
                split_line.remove( rm )



            if len( split_line ) > 1 :
                print 'Cannot handle line ', line
                print split_line
                print 'Please check the format of the line and fix the code to read it properly'
                raw_input('cont')

            if len( split_line) :

                save_entry = True
                if in_defines :
                    for define in in_defines :
                        if define not in active_defs :
                            save_entry = False
                            break
                if save_entry :
                    out_branches.append( {'full_entry' : full_line, 'branch' : split_line[0] })
            
    # clean extra characters from names
    cleaned_branches = []
    for brdic in out_branches :
        cleaned = brdic['branch'].replace(';', '').replace('*','')
        if cleaned == '' :
            continue
        cleaned_branches.append( {'full_entry' : brdic['full_entry'], 'branch' : cleaned } )


    return cleaned_branches 

def make_exe_command( exe_path, conf_file, output_loc ) :

    command = [exe_path,
               ' --conf_file %s' %conf_file,
               ' | tee %s/stdout' %(output_loc),
              ]

    return ' '.join(command)

def generate_multiprocessing_commands( file_evt_list, alg_list, exe_path, kwargs ) :

    commands = []
    # if the number of files is less than 
    # nproc, split the list so that it can be
    # used for multiprocessing
    file_evt_list_mod = []
    if kwargs['nproc'] > 1 :
        if len(file_evt_list) < kwargs['nproc'] : 
            for entry in file_evt_list :
                file = entry[0]
                for evtrange in entry[1] :
                    file_evt_list_mod.append( (file, [evtrange]) )
    else :
        for entry in file_evt_list :
            file = entry[0]
            for evtrange in entry[1] :
                file_evt_list_mod.append( (file, [evtrange]) )

    if not file_evt_list_mod :
        file_evt_list_mod = file_evt_list

    for idx, file_split in enumerate(file_evt_list_mod) :
        jobid = 'Job_%04d' %idx
        outputDir = kwargs['outputDir']
        conf_file = '%s/%s_%s.%s' %(outputDir, kwargs['confFileName'].split('.')[0], jobid, '.'.join(kwargs['confFileName'].split('.')[1:]))
        if not os.path.isdir( outputDir ) :
            os.makedirs( outputDir )
        
        write_config(alg_list, conf_file, kwargs['treeName'], outputDir, kwargs['outputFile'], [file_split], kwargs['storagePath'], kwargs['sample'], kwargs['disableOutputTree'], kwargs['nPrint'], kwargs['copyInputFiles'], idx )
        #commands.append( ( jobid, make_exe_command( exe_path, conf_file, outputDir+'/'+jobid ) ) )
        #write_config(alg_list, conf_file, options.treeName, outputDir, options.outputFile, [file_split], options.storagePath, options.sample, options.disableOutputTree, idx )
        #commands.append( ( jobid, make_exe_command( exe_path, conf_file, outputDir+'/'+jobid ) ) )
        commands.append( ( jobid, (exe_path, conf_file, outputDir+'/'+jobid) ) )

    return commands

def filter_jobs_for_resubmit( orig_commands, storagePath, outputDir, outputFile, treeName=None ) :

    commands = []
        
    if storagePath is not None :
        for jobid, cmd_info in orig_commands:
            file_path = '%s/%s/' %(storagePath, jobid )
            exists = False
            for top, dirs, files, sizes in eosutil.walk_eos(file_path) :
                if outputFile in files :
                    exists = True

            if not exists :
                commands.append( ( jobid, cmd_info ) )
    else :
        for jobid, cmd_info in orig_commands :
            file_path = '%s/%s/' %(outputDir, jobid )
            exists = False
            for top, dirs, files in os.walk(file_path) :
                if outputFile in files :
                    exists = True

            if not exists :
                commands.append( ( jobid, cmd_info ) )


    return commands

def get_file_evt_map( input_files, nsplit, nFilesPerJob, totalEvents, treeName ) :

    if not input_files :
        return []

    nFilesPerSplit = 0
    if nFilesPerJob > 0 :
        nFilesPerSplit = nFilesPerJob
    elif nFilesPerJob == 0 and nsplit > 0 :
        # split by the number of files
        nFilesPerSplit = int(math.ceil(float(len(input_files))/nsplit))

    # by default put all the files into one job
    if nFilesPerSplit == 0 :
        nFilesPerSplit = len(input_files)
    # split into sub-lists based on the number of jobs.  last job may have fewer files

    split_files = [input_files[i:i+nFilesPerSplit] for i in range(0, len(input_files), nFilesPerSplit)]
    if nsplit > 0 and nFilesPerJob > 0 :
        print 'Both --nsplit and --nFilesPerJob were provided.  This could result in only a part of the sample being analyzed'
        split_files = split_files[0:nsplit]

    # now split by events.  If no further splitting is requested
    # then each file runs over the full set of events
    files_evtsplit = []

    # determine how many times to split each file.  Do this by 
    # adding a split to each file until the number
    # of splits is achieved
    files_nsplit = [ 1 ]*len(split_files)
    n_addtl_split = nsplit - len(split_files) 
    files_idx = 0
    for sp in range(0, n_addtl_split ) :
        files_nsplit[files_idx] += 1

        # increment the index
        # if we're at the end of the list go back to the beginning
        files_idx += 1
        if files_idx >= len(files_nsplit) :
            files_idx = 0

    #get the total number of events for each file
    files_nevt = [0]*len(split_files)
    for idx, files in enumerate(split_files) :
        for file in files :
            tmp = ROOT.TFile.Open(file)
            tree = tmp.Get(treeName)
            files_nevt[idx] += tree.GetEntries()
            tmp.Close()

    # for each file get the range to use
    for files, nsplit, filetotevt in zip(split_files, files_nsplit, files_nevt) :

        if totalEvents > 0 :
            totevt = totalEvents
        else :
            totevt = filetotevt
        
        splitlist = []
        split_base = int(totevt)/int(nsplit)
        prev_val = 0
        for splitidx in range(0, nsplit) :
            # at the last entry the upper limit is the number of events
            # this prevents any missed events
            if splitidx == nsplit-1 :
                splitlist.append( (prev_val, totevt) )
            else :
                splitlist.append( (prev_val, prev_val+split_base) )

            prev_val = prev_val + split_base

        files_evtsplit.append(splitlist)

    assert len(files_evtsplit) == len(split_files), 'ERROR - size mismatch.  This should not happen!'

    split_files_evt_match = []
    for flist, evtsplit in zip(split_files, files_evtsplit) :
        split_files_evt_match.append( (flist, evtsplit ))

    return split_files_evt_match

def write_config( alg_list, filename, treeName, outputDir, outputFile, files_list, storage_path=None, sample=None, disableOutputTree=False, nPrint=None, copyInputFiles=False, start_idx=0 ) :

    base_dir = os.path.split( filename )[0]
    if not os.path.isdir(base_dir ) :
        os.makedirs( base_dir)

    cfile = open( filename, 'w')
    cfile.write( '# %s\n' %datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S') )
    # first write header information
    jobid = start_idx
    file_line = ''
    for flist, evtlist in files_list :
        file_line += flist.__str__().replace('\'','').replace(' ', '')
        file_line += '['
        for min, max in evtlist :
            file_line += '%d:(%d-%d),' %( jobid, min, max )
            jobid += 1
        file_line.rstrip(',')
        file_line += '];'

    cfile.write( 'files : %s\n' %file_line )
    cfile.write( 'treeName : %s\n' %( treeName ) )
    cfile.write( 'outputDir : %s\n' %( outputDir ) )
    cfile.write( 'outputFile : %s\n' %( outputFile ) )
    if storage_path is not None :
        cfile.write( 'storagePath : %s\n' %storage_path )
    if nPrint is not None :
        cfile.write( 'nPrint : %d\n' %nPrint )
    if sample is not None :
        cfile.write( 'sample : %s\n' %sample )
    if disableOutputTree :
        cfile.write( 'disableOutputTree : true\n' )
    if copyInputFiles :
        cfile.write( 'copyInputFiles : true\n' )

    cfile.write( '__Modules__\n' )

    for alg in alg_list :
        conf_string = alg.name + ' : '
        for name, val in inspect.getmembers(alg) :

            if name[0] == '_' : continue

            if name[0:3] == 'cut' :
                # check if this cut is inverted
                inv_str = ''
                if alg.is_inverted(name) :
                    inv_str = '!'
                conf_string += '%s %s[%s] ; ' %(name, inv_str, val)

        for name, val in alg.vars.iteritems() :
            conf_string += 'init_%s [%s] ; ' %( name, val )

        if alg.do_cutflow :
            conf_string += 'do_cutflow [] ; '

        if alg.hists :
        
            for hist in alg.hists :
                conf_string+='hist [%s,%d,%f,%f];' %(hist['name'], hist['nbin'], hist['xmin'], hist['xmax'])

        cfile.write(conf_string + '\n')

    cfile.close()

#-----------------------------------------------------------------------------------
def create_job_desc_file(command_info, kwargs) :

    priority              = kwargs.get('priority', 0)

    # the header of the job description file
    desc_entries = [
                    '#Use only the vanilla universe',
                    'universe = vanilla',
                    '# This is the executable to run.  If a script,',
                    '#   be sure to mark it "#!<path to interp>" on the first line.',
                    '# Filename for stdout, otherwise it is lost',
                    'output = stdout.txt',
                    'error = stderr.txt',
                    '# Copy the submittor environment variables.  Usually required.',
                    'getenv = True',
                    '# Copy output files when done.  REQUIRED to run in a protected directory',
                    'when_to_transfer_output = ON_EXIT_OR_EVICT',
                    'priority=%d' %priority
                    ]

    for job_index, args in command_info :

        initialdir = os.path.dirname( args[1] ) + '/' + job_index

        if not os.path.isdir( initialdir ) :
            os.mkdir( initialdir )
 
        desc_entries += [
                         'Executable = %s' %args[0],
                         'Initialdir = %s' %initialdir,
                         '# This is the argument line to the Executable'
                        ]
        # assemble the argument command
        arg_command = ['arguments = "',
                       '--conf_file %s' %args[1],
                      ]

        arg_command += ['"']

        desc_entries +=  [ ' '.join(arg_command) ]

        desc_entries += [
                          '# Queue job',
                          'queue'
                        ]

    #Determine the job description file name
    desc_name = 'job_desc.txt'

    desc_file_name = os.path.dirname( args[1] ) +'/'+desc_name

    descf = open(desc_file_name, 'w')

    descf.write('\n'.join(desc_entries))
    descf.close()

    return desc_file_name

##-----------------------------------------------------------------------------------
# The original function from submit_jobs_penn
##-----------------------------------------------------------------------------------
#def create_job_desc_file(command_info, kwargs) :
#
#    # input options
#    job_file_list         = kwargs.get('job_file_list')
#    formatted_file_list   = kwargs.get('formatted_file_list')
#    output                = kwargs.get('output')
#    xrootd_mirror         = kwargs.get('xrootd_mirror', None)
#    script                = kwargs.get('script')
#    module_args           = kwargs.get('module_args', '')
#    input_command         = kwargs.get('input_command')
#    additional_input_path = kwargs.get('additional_input_path', None)
#    #transfer_input_files  = kwargs.get('transfer_input_files', None)
#    transfer_output_files = kwargs.get('transfer_output_files', False)
#    priority              = kwargs.get('priority', 0)
#    starting_job_id       = kwargs.get('starting_job_id', 0)
#    filelist_name         = kwargs.get('filelist_name', None)
#    filelist_sep          = kwargs.get('filelist_sep', '\n')
#
#    # The user has requested to transfer the input files
#    # If the input files are on the xrootd space, we have to create
#    # a script that copies them manually and use that later
#    orig_script = script
#    files_output = output
#    if xrootd_mirror is not None :
#        files_output = xrootd_mirror
#
#    orig_script_name = os.path.basename(orig_script)
#    
#    newname = make_transfer_script( files_output, orig_script_name, kwargs )
#    script = files_output+'/' + newname
#
#    if transfer_input_files :
#        #Generate or get file transfer queue
#        #Generate in 10 minute intervals
#        curr_gmtime = time.gmtime()
#        curr_min = (int(curr_gmtime[4])/10)*10
#        queue_name = time.strftime('%Y_%m_%d_%H',curr_gmtime)
#        queue_name += '_%d' %curr_min
#        user = getpass.getuser()
#        file_queue_dir = '/disk/userdata00/atlas_data2/%s/condor_interface_queues/%s' %(user,queue_name)
#        if not os.path.isdir( file_queue_dir ) :
#            os.makedirs( file_queue_dir )
#            os.system('touch %s/queue.txt' %file_queue_dir )
#
#
#
#    # the header of the job description file
#    desc_entries = [
#                    '#Use only the vanilla universe',
#                    'universe = vanilla',
#                    '# This is the executable to run.  If a script,',
#                    '#   be sure to mark it "#!<path to interp>" on the first line.',
#                    'Executable = %s' %script,
#                    '# Filename for stdout, otherwise it is lost',
#                    'output = stdout.txt',
#                    'error = stderr.txt',
#                    '# Copy the submittor environment variables.  Usually required.',
#                    'getenv = True',
#                    '# Copy output files when done.  REQUIRED to run in a protected directory',
#                    'when_to_transfer_output = ON_EXIT_OR_EVICT',
#                    'priority=%d' %priority
#                    ]
#
#    
#    #if xrootd_mirror is None :
#    #    desc_entries += [ 'should_transfer_files = YES' ]
#
#    for job_index, files_arg in formatted_file_list.iteritems() :
#
#        outid = 'Job_%04d' %(job_index+starting_job_id)
#        initialdir = '%s/%s' %(files_output, outid)
#
#        job_module_args = module_args
#
#        desc_entries += [
#                         'Initialdir = %s' %initialdir,
#                         '# This is the argument line to the Executable'
#                        ]
#        # assemble the argument command
#
#        arg_command = ['arguments = "',
#                       '--script %s' %orig_script,
#                       '--files %s' %files_arg,
#                       '--input-command \' %s \'' %input_command,
#                       '--module-args \' %s \'' %job_module_args,
#                      ]
#
#        # if transferring files, pass the 
#        # queue to the job.  otherwise tell the script to 
#        # not transfer files
#        #if transfer_input_files :
#        #    arg_command += ['--file-queue-dir %s' %file_queue_dir ]
#        #else :
#        arg_command += [ '--no-transfer-input' ]
#
#        # if transferring output files, pass the output
#        if transfer_output_files :
#            new_output = output+'/Job_%04d' %(job_index+starting_job_id)
#            arg_command += [ '--output %s '%new_output ]
#
#        if additional_input_path is not None :
#            arg_command += [ ' --additional-input %s' %(additional_input_path % outid) ]
#
#        if filelist_name is not None :
#            arg_command += [ ' --filelist-name %s --filelist-sep %s ' %( filelist_name, filelist_sep )]
#
#        # End of arguments.  Add additional arguments above
#        arg_command += ['"']
#
#        desc_entries +=  [ ' '.join(arg_command) ]
#
#        desc_entries += [
#                          '# Queue job',
#                          'queue'
#                        ]
#
#    #Determine the job description file name
#    desc_name = 'job_desc.txt'
#
#    desc_file_name = files_output+'/'+desc_name
#
#    descf = open(desc_file_name, 'w')
#
#    descf.write('\n'.join(desc_entries))
#    descf.close()
#
#    print 'cp %s %s' %(desc_file_name, output) 
#    os.system('cp %s %s' %(desc_file_name, output) )
#
#    return desc_file_name

def write_header_files( brdefname, linkdefname, branches, out_branches=[], keep_branches=[], alg_list=[] ) :

    branch_header = open(brdefname, 'w')
    branch_header.write('#ifndef BRANCHDEFS_H\n')
    branch_header.write('#define BRANCHDEFS_H\n')
    branch_header.write('#include "TTree.h"\n')
    branch_header.write('#include <vector>\n')

    branch_header.write('//This next bit of code defines whether each branch\n'
                        '//exists or not.  This can be used to solve the problem\n'
                        '//when a branch in available in some ntuples, but not others.\n'
                        '//If this happens, the code will not compile because the\n'
                        '//branch is not written in the header file.  To fix this problem\n'
                        '//simply surround the offending code with #ifdef EXISTS_MYVAR ... #endif\n'
                        '//and if the variable does not exist the preprocessor will ignore that code\n')

    active_defs = get_active_defines( branches, alg_list )
    for define in active_defs :
        branch_header.write('#define %s\n'%define)

    branch_header.write('//Define variables as extern below and declare them in the .cxx file to avoid multiple definitions\n')
    branch_header.write('namespace IN {\n');

    for conf in branches :

        name = conf['name']

        # declare the variable differently for a vector
        if conf['type'].count('vector') :
            modtype = conf['type'].replace('vector', 'std::vector')
            modtype = modtype.replace('string', 'std::string')
            file_entry = ' extern %s\t\t\t\t*%s;\n' %(modtype, name)
        else :
            file_entry = ' extern %s\t\t\t\t%s%s;\n' %(conf['type'], name, conf['arrayStr'])

        branch_header.write(file_entry)


    branch_header.write('};\n');

    branch_header.write('namespace OUT {\n');

    for conf in branches :

        name = conf['name']

        if name not in keep_branches :
            continue

        if conf['type'].count('vector') :
            modtype = conf['type'].replace('vector', 'std::vector')
            modtype = modtype.replace('string', 'std::string')
            file_entry = ' extern %s\t\t\t\t*%s;\n' %(modtype, name)
        else :
            file_entry = ' extern %s\t\t\t\t%s%s;\n' %(conf['type'], name, conf['arrayStr'])

        branch_header.write(file_entry)

    for brdic in out_branches :
        file_entry = ' extern %s ' %brdic['full_entry'] 
        branch_header.write(file_entry)

    branch_header.write('};\n');
    branch_header.write('#endif\n')

    branch_header.close()

    # this file is used to create root dictionaries
    # for vector<vector<>> types.  It will always be
    # created and if no such types exist it will
    # still be compiled, but will do nothing

    # first collect a list of all the types that need to be added
    link_types = []
    for conf in branches :
        if conf['type'].count('vector') > 0 :
            modtype = conf['type'].replace('vector','std::vector')
            link_types.append( modtype )

    # remove duplicates
    link_types = set( link_types )

    link_header = open(linkdefname, 'w')
    link_header.write('#ifdef __CINT__\n')
    link_header.write('#pragma link off all globals;  \n')
    link_header.write('#pragma link off all classes;  \n')
    link_header.write('#pragma link off all functions;\n')
    link_header.write('#pragma link C++ nestedclasses;\n\n')

    for type in link_types :
        link_header.write('#pragma link C++ class %s+;\n' %type)

    link_header.write('#endif\n\n')

    link_header.close()
    

def write_source_file(source_file_name, header_file_name, branches, out_branches=[], keep_branches=[], debugCode=False, write_expert_code=False) :

    branch_header = open(header_file_name, 'w')
    branch_header.write('#ifndef BRANCHINIT_H\n')
    branch_header.write('#define BRANCHINIT_H\n')
    branch_header.write('#include "TTree.h"\n')
    branch_header.write('#include "TChain.h"\n')
    branch_header.write('void InitINTree( TChain * tree );\n')
    branch_header.write('void InitOUTTree( TTree * tree );\n')
    branch_header.write('void CopyInputVarsToOutput(std::string prefix = std::string() );\n')
    branch_header.write('void CopyPrefixBranchesInToOut( const std::string & prefix );\n' )
    branch_header.write('void CopyPrefixIndexBranchesInToOut( const std::string & prefix, unsigned index, bool quiet=false, std::string veto = std::string() );\n')
    branch_header.write('void CopyMapVarsToOut( const std::map<std::string, float> & results);\n')
    branch_header.write('void CopyVectorMapVarsToOut( const std::map<std::string, std::vector< float> > & results);\n')
    branch_header.write('void CopyIntMapVarsToOut( const std::map<std::string, int> & results);\n')
    branch_header.write('void CopyIntVectorMapVarsToOut( const std::map<std::string, std::vector< int > > & results);\n')
    branch_header.write('void ClearOutputPrefix ( const std::string & prefix );\n')
    if write_expert_code:
        branch_header.write('void CheckVectorSize ( const std::string & prefix, unsigned expected );\n')

    for br in branches :
        name = br['name']
        if name not in keep_branches :
            continue
        branch_header.write('void Copy%sInToOut( std::string prefix = std::string() ); \n' %name)

        if br['type'].count('vector') :
            branch_header.write('void Copy%sInToOutIndex( unsigned index, std::string prefix = std::string(), bool quiet=false, std::string veto = std::string() ); \n' %name)
            branch_header.write('void ClearOutput%s( std::string prefix ); \n' %name)
            if write_expert_code:
                branch_header.write('void Check%sVectorSize( std::string prefix, unsigned expected  ); \n' %name)


    branch_header.write('#endif\n')
    branch_header.close()

    print 'WRITE ', source_file_name

    branch_setting = open(source_file_name, 'w')
    branch_setting.write('#include <algorithm>\n')
    branch_setting.write('#include <iostream>\n')
    branch_setting.write('#include "TTree.h"\n')
    branch_setting.write('#include "TChain.h"\n')
    branch_setting.write('#include "include/BranchInit.h"\n' )
    branch_setting.write('#include "include/BranchDefs.h"\n\n' )

    branch_setting.write('namespace IN {\n');
    for conf in branches :

        name = conf['name']

        # declare the variable differently for a vector
        if conf['type'].count('vector') :
            modtype = conf['type'].replace('vector', 'std::vector')
            modtype = modtype.replace('string', 'std::string')
            file_entry = ' %s\t\t\t\t*%s;\n' %(modtype, name)
        else :
            file_entry = ' %s\t\t\t\t%s%s;\n' %(conf['type'], name, conf['arrayStr'])

        branch_setting.write(file_entry)

    branch_setting.write('};\n');

    branch_setting.write('namespace OUT {\n');

    for conf in branches :

        name = conf['name']

        if name not in keep_branches :
            continue

        if conf['type'].count('vector') :
            modtype = conf['type'].replace('vector', 'std::vector')
            modtype = modtype.replace('string', 'std::string')
            file_entry = ' %s\t\t\t\t*%s;\n' %(modtype, name)
        else :
            file_entry = ' %s\t\t\t\t%s%s;\n' %(conf['type'], name, conf['arrayStr'])

        branch_setting.write(file_entry)

    branch_setting.write('};\n');

    branch_setting.write('void InitINTree( TChain * tree) {\n\n')
    counter = -1
    
    for conf in branches :

        name = conf['name']

        counter += 1
        # might need to fix this for variable length arrays
        set_line = '  tree->SetBranchAddress("%s", &IN::%s);\n' %(name, name)
        #if type.count('vector') :
        #else :
        #    set_line = '  tree->SetBranchAddress("%s", &IN::%s);\n' %(name, name)

        branch_setting.write(set_line)

    branch_setting.write('};\n\n')

    branch_setting.write('void InitOUTTree( TTree * tree ) {\n')

    for conf in branches :

        name = conf['name']

        if name not in keep_branches :
            continue

        if conf['type'].count('vector') :
            set_line = '  tree->Branch("%s", &OUT::%s);\n' %(name, name)
        else :
            set_line = '  tree->Branch("%s", &OUT::%s, "%s");\n' %(name, name, conf['leafEntry'])

        branch_setting.write(set_line)

    branch_setting.write('}\n')

    branch_setting.write('void CopyInputVarsToOutput( std::string prefix) {\n')

    for conf in branches :

        name = conf['name']
        if name not in keep_branches :
            continue

        branch_setting.write('    Copy%sInToOut( prefix ); \n' %name )

    branch_setting.write('}\n\n')

    branch_setting.write('// The next set of functions allows one to copy \n')
    branch_setting.write('// input variables to the outputs based on a key\n')
    branch_setting.write('// A copy function is generated for each pair of variables\n')
    branch_setting.write('// The copy function holds the name of the function to compare\n')
    branch_setting.write('// to the input key.  If the variables are vectors, a second function\n')
    branch_setting.write('// is generated that allows one to copy all variables matching a given key\n')
    branch_setting.write('// at a certain index and pushes that back on the output variable\n\n')

    branch_setting.write('void CopyPrefixBranchesInToOut( const std::string & prefix ) {\n' )
    branch_setting.write('// Just call each copy function with the prefix \n\n')

    for conf in branches :
        name = conf['name']
        if name not in keep_branches :
            continue

        branch_setting.write( '    Copy%sInToOut( prefix );\n' %name)

    branch_setting.write('}; \n\n' )
    
    branch_setting.write('void CopyPrefixIndexBranchesInToOut( const std::string & prefix, unsigned index, bool quiet, std::string veto ) { \n\n')

    branch_setting.write('// Just call each copy function with the prefix \n\n')

    for conf in branches :
        name = conf['name']
        if name not in keep_branches :
            continue

        if conf['type'].count('vector') :
            branch_setting.write( '    Copy%sInToOutIndex( index, prefix, quiet, veto );\n' %name)

    branch_setting.write('}; \n\n' )

    branch_setting.write('void CopyMapVarsToOut( const std::map<std::string, float> & results ) {\n')

    valid_branches = []
    for conf in branches :
        if conf['type'].count('vector') :
            continue
        if conf['arrayStr'] :
            continue

        name = conf['name']
        if name not in keep_branches :
            continue


        valid_branches.append( name )

    for brdic in out_branches :

        if not brdic['full_entry'].count('vector') :
            valid_branches.append( brdic['branch'] )

    branch_setting.write( '    std::vector<std::string> filled_branches;\n' )
    branch_setting.write( '    std::map<std::string, float>::const_iterator mitr; \n'  )
    for br in valid_branches :
        branch_setting.write( '    mitr = results.find( "%s" ); \n' %br )
        branch_setting.write( '    if( mitr != results.end() ) { \n' )
        branch_setting.write( '        OUT::%s = mitr->second;\n' %br )
        branch_setting.write( '        filled_branches.push_back( "%s" );\n' %(br) )
        branch_setting.write( '    }\n' )
    branch_setting.write( '    if( filled_branches.size() < results.size() ) { \n' )
    branch_setting.write( '        for( std::map<std::string, float>::const_iterator mitr = results.begin(); mitr != results.end(); ++mitr ) { \n' )
    branch_setting.write( '            if( std::find( filled_branches.begin(), filled_branches.end(), mitr->first ) == filled_branches.end()  ) { \n' )
    branch_setting.write( '                std::cout << "Could not fill branch " << mitr->first << std::endl; \n '  )
    branch_setting.write( '            }\n' )
    branch_setting.write( '        }\n' )
    branch_setting.write( '    }\n' )
    branch_setting.write( '}\n' )

    branch_setting.write('void CopyVectorMapVarsToOut( const std::map<std::string, std::vector<float> > & results ) {\n')

    valid_branches = []
    for conf in branches :
        nvec = conf['type'].count('vector')
        if not nvec==1 :
            continue

        if not conf['type'].count('float') :
            continue

        name = conf['name']
        if name not in keep_branches :
            continue

        valid_branches.append( name )

    for brdic in out_branches :

        if brdic['full_entry'].count('vector') and brdic['full_entry'].count('float') :
            valid_branches.append( brdic['branch'] )

    branch_setting.write( '    std::vector<std::string> filled_branches;\n' )
    branch_setting.write( '    std::map<std::string, std::vector<float> >::const_iterator mitr; \n' )
    for br in valid_branches :
        branch_setting.write( '    mitr = results.find( "%s" ); \n' %br )
        branch_setting.write( '    if( mitr != results.end() ) { \n' )
        branch_setting.write( '        OUT::%s->assign(mitr->second.begin(), mitr->second.end() );\n' %br )
        branch_setting.write( '        filled_branches.push_back( "%s" );\n' %(br) )
        branch_setting.write( '    }\n' )
    branch_setting.write( '    if( filled_branches.size() < results.size() ) { \n' )
    branch_setting.write( '        for( std::map<std::string, std::vector<float> >::const_iterator mitr = results.begin(); mitr != results.end(); ++mitr ) { \n' )
    branch_setting.write( '            if( std::find( filled_branches.begin(), filled_branches.end(), mitr->first ) == filled_branches.end()  ) { \n' )
    branch_setting.write( '                std::cout << "Could not fill branch " << mitr->first << std::endl; \n '  )
    branch_setting.write( '            }\n' )
    branch_setting.write( '        }\n' )
    branch_setting.write( '    }\n' )
    branch_setting.write( '}\n' )

    branch_setting.write('void CopyIntMapVarsToOut( const std::map<std::string, int> & results ) {\n')

    valid_branches = []
    for conf in branches :
        if conf['type'].count('vector') :
            continue
        if conf['arrayStr'] :
            continue

        name = conf['name']
        if name not in keep_branches :
            continue


        valid_branches.append( name )

    for brdic in out_branches :

        if not brdic['full_entry'].count('vector') :
            valid_branches.append( brdic['branch'] )

    branch_setting.write( '    std::vector<std::string> filled_branches;\n' )
    branch_setting.write( '    std::map<std::string, int>::const_iterator mitr; \n'  )
    for br in valid_branches :
        branch_setting.write( '    mitr = results.find( "%s" ); \n' %br )
        branch_setting.write( '    if( mitr != results.end() ) { \n' )
        branch_setting.write( '        OUT::%s = mitr->second;\n' %br )
        branch_setting.write( '        filled_branches.push_back( "%s" );\n' %(br) )
        branch_setting.write( '    }\n' )
    branch_setting.write( '    if( filled_branches.size() < results.size() ) { \n' )
    branch_setting.write( '        for( std::map<std::string, int>::const_iterator mitr = results.begin(); mitr != results.end(); ++mitr ) { \n' )
    branch_setting.write( '            if( std::find( filled_branches.begin(), filled_branches.end(), mitr->first ) == filled_branches.end()  ) { \n' )
    branch_setting.write( '                std::cout << "Could not fill branch " << mitr->first << std::endl; \n '  )
    branch_setting.write( '            }\n' )
    branch_setting.write( '        }\n' )
    branch_setting.write( '    }\n' )
    branch_setting.write( '}\n' )

    branch_setting.write('void CopyIntVectorMapVarsToOut( const std::map<std::string, std::vector<int> > & results ) {\n')

    valid_branches = []
    for conf in branches :
        nvec = conf['type'].count('vector')
        if not nvec==1 :
            continue

        if not conf['type'].count('int') :
            continue

        name = conf['name']
        if name not in keep_branches :
            continue

        valid_branches.append( name )

    for brdic in out_branches :

        if brdic['full_entry'].count('vector') and brdic['full_entry'].count('int') :
            valid_branches.append( brdic['branch'] )

    branch_setting.write( '    std::vector<std::string> filled_branches;\n' )
    branch_setting.write( '    std::map<std::string, std::vector<int> >::const_iterator mitr; \n' )
    for br in valid_branches :
        branch_setting.write( '    mitr = results.find( "%s" ); \n' %br )
        branch_setting.write( '    if( mitr != results.end() ) { \n' )
        branch_setting.write( '        OUT::%s->assign(mitr->second.begin(), mitr->second.end() );\n' %br )
        branch_setting.write( '        filled_branches.push_back( "%s" );\n' %(br) )
        branch_setting.write( '    }\n' )
    branch_setting.write( '    if( filled_branches.size() < results.size() ) { \n' )
    branch_setting.write( '        for( std::map<std::string, std::vector<int> >::const_iterator mitr = results.begin(); mitr != results.end(); ++mitr ) { \n' )
    branch_setting.write( '            if( std::find( filled_branches.begin(), filled_branches.end(), mitr->first ) == filled_branches.end()  ) { \n' )
    branch_setting.write( '                std::cout << "Could not fill branch " << mitr->first << std::endl; \n '  )
    branch_setting.write( '            }\n' )
    branch_setting.write( '        }\n' )
    branch_setting.write( '    }\n' )
    branch_setting.write( '}\n' )

    branch_setting.write('void ClearOutputPrefix ( const std::string & prefix ) {\n')
    for conf in branches :
        name = conf['name']
        if name not in keep_branches :
            continue

        if conf['type'].count('vector') :
            branch_setting.write( '    ClearOutput%s( prefix );\n' %name)

    branch_setting.write('}; \n\n' )

    if write_expert_code :
        branch_setting.write('void CheckVectorSize( const std::string & prefix, unsigned expected ) { \n\n')

        branch_setting.write('// Just call each function with the prefix \n\n')

        for conf in branches :
            name = conf['name']
            if name not in keep_branches :
                continue

            if conf['type'].count('vector') :
                branch_setting.write( '    Check%sVectorSize( prefix, expected);\n' %name)

        branch_setting.write('}; \n\n' )
    
    for conf in branches :
        name = conf['name']
        if name not in keep_branches :
            continue

        # determine how to copy the variable based on
        # its type
        set_line = ''
        if conf['type'].count('vector') :
            modtype = conf['type'].replace('vector', 'std::vector')
            modtype = modtype.replace('string', 'std::string')
            set_line = '  *OUT::%s = %s(*IN::%s);\n' %(name, modtype, name)
        else :
            if conf['totSize'] > 1 :
                #set_line = '  memcpy(OUT::%s, IN::%s, %d);\n' %(name, name, conf['totSize'])
                if conf['totSize'] < 200 :
                    set_line = '  std::copy(IN::%s, IN::%s+%d, OUT::%s);\n ' %(name, name, conf['totSize'], name)
                    if debugCode :
                        set_line = ' std::cout << "Attempt to copy variable %s" << std::endl;\n ' %name + set_line
                else :
                    sizes = conf['sizeEntries']
                    set_line = ''
                    for array_val in generate_array_loop( sizes ) :
                        set_line += '  OUT::%s%s = IN::%s%s;\n ' %(name, array_val, name, array_val)


                    # code crashes when copying large arrays the only solution I've found is to do it manually

            else :
                set_line = '  OUT::%s = IN::%s;\n' %(name, name)

        branch_setting.write('void Copy%sInToOut( std::string prefix ) { \n\n' %name)
        branch_setting.write('    std::string my_name = "%s";\n' %name)
        branch_setting.write('    std::size_t pos = my_name.find( prefix ); \n')
        branch_setting.write('    // if the filter is given only continue if its matched at the beginning \n' )
        branch_setting.write('    if( prefix != "" &&  pos != 0 ) return; \n' )
        branch_setting.write(set_line )
        branch_setting.write('}; \n\n ')

        if conf['type'].count('vector') :
            branch_setting.write('void Copy%sInToOutIndex( unsigned index, std::string  prefix, bool quiet, std::string veto ) { \n\n' %name)
            branch_setting.write('    std::string my_name = "%s";\n' %name)
            branch_setting.write('    std::size_t pos = my_name.find( prefix ); \n')
            branch_setting.write('    std::size_t pos_veto = std::string::npos; \n')
            branch_setting.write('    if( !veto.empty() ) pos_veto = my_name.find( veto ); \n')
            branch_setting.write('    // if the filter is given only continue if its matched at the beginning \n' )
            branch_setting.write('    if( prefix != "" &&  pos != 0 ) return; \n' )
            branch_setting.write('    if( pos_veto == 0 ) return; \n' )
            branch_setting.write('    if( index >= IN::%s->size() ) {\n ' %name)
            branch_setting.write('        if( !quiet ) {\n ' )
            branch_setting.write('        std::cout << "Vector size exceeded for branch IN::%s" << std::endl;\n ' %name)
            branch_setting.write('        }; \n')
            branch_setting.write('        return; \n ')
            branch_setting.write('    }; \n\n ')
            branch_setting.write('    //std::cout << "Copy varaible %s" << " at index " << index << ", prefix = " << prefix << std::endl; \n ' %name)
            branch_setting.write('    OUT::%s->push_back( IN::%s->at(index) ); \n ' %(name, name))
            branch_setting.write('}; \n\n ')

            branch_setting.write('void ClearOutput%s( std::string  prefix ) { \n\n' %name)
            branch_setting.write('    std::string my_name = "%s";\n' %name)
            branch_setting.write('    std::size_t pos = my_name.find( prefix ); \n')
            branch_setting.write('    // if the filter is given only continue if its matched at the beginning \n' )
            branch_setting.write('    if( prefix != "" &&  pos != 0 ) return; \n' )
            branch_setting.write('    //std::cout << "Clear varaible %s, prefix = " << prefix << std::endl; \n ' %name)
            branch_setting.write('    OUT::%s->clear(); \n ' %name)
            branch_setting.write('}; \n\n ')

            if write_expert_code :
                branch_setting.write('void Check%sVectorSize( std::string  prefix, unsigned expected ) { \n\n' %name)
                branch_setting.write('    std::string my_name = "%s";\n' %name)
                branch_setting.write('    std::size_t pos = my_name.find( prefix ); \n')
                branch_setting.write('    // if the filter is given only continue if its matched at the beginning \n' )
                branch_setting.write('    if( prefix != "" &&  pos != 0 ) return; \n' )
                branch_setting.write('    if( IN::%s->size() != expected ) { \n'  %name)
                branch_setting.write('        std::cout << "Vector size for branch %s is not as expected.  Got " << IN::%s->size() << ", expected " << expected << std::endl; \n'%(name, name) )
                branch_setting.write('    } \n' )
                branch_setting.write('}; \n\n ')



    branch_setting.close()

def get_active_defines( branches, alg_list ) :
    #-----------------------
    # Define a preprocessor variable
    # for each input variable to
    # allow the code to check if
    # the variable exists
    #-----------------------
    defs = []
    for conf in branches :
        name = conf['name']
        defs.append('EXISTS_%s\n'%name)

    #-----------------------
    # Define a preprocessor variable
    # for each module that is scheduled
    # to allow the code to check
    # if that module is used
    #-----------------------
    all_filters = [f.name for f in alg_list]
    all_filters = list( set( all_filters ) ) 
    for filter in all_filters :
        defs.append( 'MODULE_%s\n' %filter )
    
    return defs

def generate_array_loop( array_sizes ) :

    def _increment( array_val, array_sizes, depth ) :

        if depth == -1 :
            return array_val
        if array_val[depth] == ( array_sizes[depth] - 1 ) : #at the last value, reset to 0 and increment the previous
            array_val[depth] = 0
            _increment( array_val, array_sizes, depth-1)
        else :
            array_val[depth] += 1

        return array_val


    n_depth = len(array_sizes)
    array_val = [0]*n_depth

    total_size = reduce(lambda x, y: x*y, array_sizes)
    for idx in range(0, total_size) :
        yield '[' + ']['.join([str(x) for x in array_val]) + ']'

        _increment( array_val, array_sizes, n_depth-1 );



class Filter :

    def __init__(self, name) :
        self.invert_list = []
        self.name = name
        self.do_cutflow = False
        self.hists = []
        self.vars = {}

    def invert(self, name) :
        self.invert_list.append(name)

    def is_inverted(self, name) :
        return name in self.invert_list

    def add_var(self, name, val) :
        self.vars[name] = val

    def add_hist(self, name, nbin, xmin, xmax ) :
        self.hists.append( {'name':name, 'nbin':nbin, 'xmin':xmin, 'xmax':xmax } )


## for python only code
#
#class AnalysisConfig :
#
#    def __init__(self) :
#        self.modules = []
#        self.files = []
#        self.keep_branches = []
#
#    def add_module(self, module) :
#        self.modules.append(module)
#
#    def get_modules(self) :
#        return self.modules
#
#    def Run( self, runmod, options ) :
#
#        if not hasattr(runmod, 'Run') :
#            print 'AnalysisConfig.Run - ERROR : Run module does not implement a function called Run'
#
#        runmod.in_tree = ROOT.TChain( options.treeName, options.treeName )
#        for file in self.files :
#            runmod.in_tree.Add( file )
#
#        runmod.outfile = ROOT.TFile.Open( '%s/%s' %(options.outputDir, options.outputFile ), 'RECREATE' )
#        runmod.outfile.cd()
#
#        # make the output tree
#        # if the input tree is within a directory, put the output tree
#        # in a directory of the same name
#        name_split = runmod.in_tree.GetName().split('/')
#        dir_path = ''
#        if len(name_split) > 1 :
#            for dir in name_split[:-1] :
#                dir_path = dir_path + dir + '/'
#                print 'TFile.mkdir %s' %dir_path
#                runmod.outfile.mkdir( dir_path )
#
#            runmod.outfile.cd( dir_path )
#
#            tree_name = name_split[-1]
#
#            runmod.out_tree = ROOT.TTree( tree_name, tree_name )
#        else :
#            runmod.out_tree = ROOT.TTRee( runmod.in_tree.GetName(), runmod.in_tree.GetName() )
#            runmod.out_tree.SetDirectory( runmod.outfile )
#
#        self.init_out_tree( runmod )
#
#        print 'OUT branches'
#        for br in runmod.out_tree.GetListOfBranches() :
#            print br.GetName()
#
#        runmod.Run( self.get_modules(), options, 0, 0 )
#
#        runmod.outfile.Write()
#        runmod.outfile.Close()
#
#    def init_out_tree( self, module ) :
#
#        all_branches = get_branch_mapping( module.in_tree )
#
#        for br in all_branches :
#            brname = br['name']
#
#            if brname not in self.keep_branches :
#                continue
#
#            if br['type'].count('vector') :
#                print 'FIX ME'
#                setattr(module.out_tree, brname, std.vector(int)())
#                out_tree.Branch(brname, getattr(module.out_tree, brname))
#            else :
#                array_type = ''
#                br_type = br['type']
#                if br_type == 'int' or br_type == 'Int_t' :
#                    array_type = 'i'
#                elif br_type == 'float' or br_type == 'Float_t' :
#                    array_type = 'f'
#                elif br_type == 'double' or br_type == 'Double_t' :
#                    array_type = 'd'
#                elif br_type == 'Long64_t' :
#                    array_type = 'L'
#                elif br_type == 'Bool_t' :
#                    array_type = 'B'
#                else :
#                    print 'No conversion for type %s' %br_type
#                    continue
#               
#                    #    root2PythonDictionary = { 
#                    #            #root:python 
#                    #            'C':'c',#char 
#                    #            'B':'b',#signed char 
#                    #            'b':'B',#unsigned char 
#                    #            'S':'h',#signed short 
#                    #            's':'H',#unsigned short 
#                    #            'I':'i',#signed int 
#                    #            'i':'I',#unsigned int 
#                    #            'L':'l',#signed long 
#                    #            'l':'L',#unsigned long 
#                    #            'F':'f',#float 
#                    #            'D':'d'#double 
#                    #    } 
#                setattr(module, brname, array.array(array_type, br['totSize']*[0]) )
#                print brname
#                print getattr(module, brname)
#                print br['leafEntry']
#                module.out_tree.Branch(brname, getattr(module, brname), br['leafEntry'])
#
#class ModuleConfig : 
#    
#    def __init__( self, alg ) :
#
#        self.name = alg.name
#
#        self.cuts = {}
#
#        for member in alg.__dict__ : 
#            if not member.count('cut_') :
#                continue
#
#            self.cuts[member] = getattr(alg, member)
#
#            
#    def pass_cut( self, cut, val ) :
#
#        # if cut hasn't been configured, then pass it
#        if cut not in self.cuts :
#            return True
#
#        return eval( 'val' + self.cuts[cut] )
#
#    def Print( self ) :
#        print 'Module name is %s' %self.name
#        print 'Cuts are : '
#        for name, val in self.cuts.iteritems() :
#            print '%s %s' %(name, val)
#
#def ConfigureAnalysis( runmod, options ) :
#
#    # read input files
#    assert options.files is not None or options.filesDir is not None, 'Must provide a list of files via --files or a directory containing files via --filesDir'
#
#    ana_config = AnalysisConfig()
#    if options.files is not None :
#        ana_config.files = options.files.split(',')
#    elif options.filesDir is not None :
#        ana_config.files = collect_input_files( options.filesDir )
#        
#    if not ana_config.files :
#        print 'No files were found'
#        return 
#
#    # get the configuration
#    ImportedModule = import_module( options.module )
#
#    alg_list = []
#    ImportedModule.config_analysis( alg_list )
#
#    for alg in alg_list :
#        ana_config.add_module(ModuleConfig(alg))
#
#    keep_filter   = []
#    remove_filter = []
#
#    branches = get_branch_mapping( ana_config.files, options.treeName )
#
#    ana_config.keep_branches = get_keep_branches( ImportedModule, branches, options.enableKeepFilter, options.enableRemoveFilter )
#
#    return ana_config
#    
#def copy_event( module ) :
#
#    for br in module.out_tree.GetListOfBranches() :
#        brname = br.GetName()
#
#        #module.nEle[0] = module.in_tree.nEle
#        print brname
#        print module.out_tree.GetBranch(brname).GetTitle()
#        brtitle = module.out_tree.GetBranch(brname).GetTitle()
#        is_array = brtitle.count('[') and brtitle.count(']')
#        #print getattr(module, brname)
#        #varlen = len(getattr(module, brname))
#        varlen = module.in_tree.GetBranch(brname).GetLeaf(brname).GetNdata()
#        print varlen
#        #print 'leaf type', module.in_tree.GetBranch(brname).GetLeaf(brname).GetTypeName()
#        #print getattr(module, brname)
#        
#        if is_array :
#            #setattr(module, brname, getattr(module.in_tree, brname) )
#            for idx in xrange(0, len(getattr(module.in_tree, brname ))) :
#                exec('module.%s[%d] = module.in_tree.%s[%d]' %(brname, idx, brname, idx) )
#        else :
#            exec('module.%s[0] = module.in_tree.%s' %(brname, brname ) )
#
#        #if varlen == 1 :
#        #    exec('module.%s[0] = module.in_tree.%s' %(brname, brname ) )
#        #else :
#        #    setattr(module, brname, getattr(module.in_tree, brname) )
#        #    #for idx in xrange(0, len(getattr(module.in_tree, brname ))) :
#        #    #    exec('module.%s[%d] = copy.copymodule.in_tree.%s[%d]' %(brname, idx, brname, idx) )
#        ##setattr(module, brname+'[0]', getattr(module.in_tree, brname) )
