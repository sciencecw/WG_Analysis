import os
import scheduler_base
from scheduler_base import JobConf

from argparse import ArgumentParser
p = ArgumentParser()
p.add_argument( '--run', dest='run', default=False, action='store_true', help='Run filtering' )
p.add_argument( '--check', dest='check', default=False, action='store_true', help='Run check of completion' )
p.add_argument( '--resubmit', dest='resubmit', default=False, action='store_true', help='Only submit missing output' )
p.add_argument( '--batch', dest='batch', default=False, action='store_true', help='Run on batch, not locally ' )
options = p.parse_args()

if not options.check :
    options.run = True
else :
    options.run = False

options.local = ( not options.batch )


base = '/data/users/jkunkle/Resonances'

# ----------------------------
# The suffix that appears on the 
# end of the output jobs
# ----------------------------

args_pholap           = { 'function' : 'filter_photon', 'pt_cut' : ' > 10 ', 'nph_cut' : ' == 0 ' }
args_PtMaxInc         = { 'function' : 'filter_photon', 'pt_cut' : ' > 180 ', 'nph_cut' : ' == 0' }
args_PtMax130         = { 'function' : 'filter_photon', 'pt_cut' : ' > 550 ', 'nph_cut' : ' == 0' }
args_PtMin130         = { 'function' : 'filter_photon', 'pt_cut' : ' > 180 ', 'nph_cut' : ' > 0' }
args_PtMin500         = { 'function' : 'filter_photon', 'pt_cut' : ' > 550 ', 'nph_cut' : ' > 0 ' }
args_HT               = { 'function' : 'filter_genht' , 'trueht_cut' :  ' < 100 ' }
args_mtMax600         = { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 600 ' }
args_mtMax1300        = { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 1300 ' }
args_mtMax500         = { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 500 ' }
args_mtMax400         = { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 400 ' }
args_mtMin400Max1200  = { 'function' : 'filter_mtres' , 'mtres_cut' : ' >= 400 && < 1200 ' }
args_mtMin1200        = { 'function' : 'filter_mtres' , 'mtres_cut' : ' >= 1200 ' }

jobs = [
    #JobConf(base, 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', tags=['HT'], suffix='TrueHTOlap'),
    #JobConf(base, 'WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8' , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'  , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'  , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'               , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'          , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'            , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8', tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'   , tags=['pholap']  , suffix = 'PhOlap'  ),
    #JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'                , tags=['PtMaxInc'], suffix = 'PhCutMax'),
    #JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'                 , tags=['PtMaxInc'], suffix = 'PhCutMax'),
    #JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'        , tags=['PtMax130'], suffix = 'PhCutMax'),
    #JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'         , tags=['PtMax130'], suffix = 'PhCutMax'),
    #JobConf(base, 'WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'        , tags=['PtMin500'], suffix = 'PhCutMin'),
    #JobConf(base, 'WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'         , tags=['PtMin500'], suffix = 'PhCutMin'),

    JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8PhCutMax', tags=['PtMin130'], suffix = 'PhCutMin'),
    JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-madgraphMLM-pythia8PhCutMax' , tags=['PtMin130'], suffix = 'PhCutMin'),
    JobConf(base, 'WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8TrueHTOlap'     , tags=['pholap']  , suffix = 'PhOlap'  ),

    ##JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8PhCutMax'                , tags=['mtMax600']       , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8PhCutMaxPhCutMin', tags=['mtMax1300']      , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8PhCutMax'                 , tags=['mtMax500']       , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-madgraphMLM-pythia8PhCutMaxPhCutMin' , tags=['mtMax1300']      , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'                        , tags=['mtMax400']       , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'                         , tags=['mtMax400']       , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'                , tags=['mtMin400Max1200'], suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_PtG-130_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'                 , tags=['mtMin400Max1200'], suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'                , tags=['mtMin1200']      , suffix = 'MTResCut'),
    ##JobConf(base, 'WGToLNuG_PtG-500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'                 , tags=['mtMin1200']      , suffix = 'MTResCut'),
]

options.nFilesPerJob = 1
options.nproc = 2
options.treename='UMDNTuple/EventTree'
options.exename='RunAnalysis'
options.copyInputFiles=True
options.enableKeepFilter=False

#input_dirs = ['SingleLepNoPhID_el_2017_02_03', 'SingleLepNoPhID_mu_2017_02_03']
#input_dirs = ['SingleLepNoPhID_mu_2017_04_12', 'SingleLepNoPhID_el_2017_04_12','LepGamma_mug_2017_04_12', 'LepGamma_elg_2017_04_12', 'LepLep_mumu_2017_04_12', 'LepLep_elel_2017_04_12', 'LepLep_muel_2017_04_12']
input_dirs =[ 
              
              #'LepGammaNoEleOlap_elg_2018_03_28',
              #'SingleLepNoPhID_mu_2018_03_28', 'SingleLepNoPhID_el_2018_03_28', 
              #'LepLep_mumu_2018_03_28', 'LepLep_muel_2018_03_28', 'LepLep_elel_2018_03_28',
              #'LepGamma_mug_2018_03_28','LepGamma_elg_2018_03_28', 
              #'LepGammaNoPhId_mug_2018_03_28',
              #'LepGammaNoPhId_elg_2018_03_28', 
              'LepGammaNoEleOlap_elg_2018_04_10', 
]

module = 'Conf.py'

configs = []

for input_dir in input_dirs : 

    configs.append( {
                     'module' : module,
                     'args' : {},
                     'input' : input_dir, 
                     'output' : base + '/' + input_dir,
                     'tag' : 'olap', 
                     'args_tag_pholap'           : { 'function' : 'filter_photon', 'pt_cut' : ' > 10 ', 'nph_cut' : ' == 0 ' },
                     'args_tag_PtMaxInc'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 180 ', 'nph_cut' : ' == 0' },
                     'args_tag_PtMax130'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 550 ', 'nph_cut' : ' == 0' },
                     'args_tag_PtMin130'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 180 ', 'nph_cut' : ' > 0' },
                     'args_tag_PtMin500'         : { 'function' : 'filter_photon', 'pt_cut' : ' > 550 ', 'nph_cut' : ' > 0 ' },
                     'args_tag_HT'               : { 'function' : 'filter_genht' , 'trueht_cut' :  ' < 100 ' },
                     'args_tag_mtMax600'         : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 600 ' },
                     'args_tag_mtMax1300'        : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 1300 ' },
                     'args_tag_mtMax500'         : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 500 ' },
                     'args_tag_mtMax400'         : { 'function' : 'filter_mtres' , 'mtres_cut' : ' < 400 ' },
                     'args_tag_mtMin400Max1200'  : { 'function' : 'filter_mtres' , 'mtres_cut' : ' >= 400 && < 1200 ' },
                     'args_tag_mtMin1200'        : { 'function' : 'filter_mtres' , 'mtres_cut' : ' >= 1200 ' },
    })
                     

scheduler_base.RunJobs( jobs, configs, options)

