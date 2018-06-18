# setup CMSSW (without explicitly having a CMSSW release)
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_25/src
eval `scramv1 runtime -csh`
cd -

setenv WorkArea ${PWD}/Analysis
if( $?PYTHONPATH == 0 ) then
   setenv PYTHONPATH ${PWD}/Analysis/TreeFilter/Core/python:${PWD}/Analysis/Util/python:${PWD}/Plotting
else 
   setenv PYTHONPATH ${PWD}/Analysis/TreeFilter/Core/python:${PWD}/Analysis/Util/python:${PWD}/Plotting:$PYTHONPATH
endif


