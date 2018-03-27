OUTPUTMAIN=~/Plots/Resonance/Plots_2018_01_22

DIR_MUG=/data/users/jkunkle/Resonances/LepGamma_mug_2018_01_21
DIR_ELG=/data/users/jkunkle/Resonances/LepGamma_elg_2018_01_21
DIR_MUGNoId=/data/users/jkunkle/Resonances/LepGammaNoPhId_mug_2018_01_21
DIR_ELGNoId=/data/users/jkunkle/Resonances/LepGammaNoPhId_elg_2018_01_21
DIR_MU=/data/users/jkunkle/Resonances/SingleLepNoPhID_mu_2018_01_21
DIR_EL=/data/users/jkunkle/Resonances/SingleLepNoPhID_el_2018_01_21
DIR_MUMU=/data/users/jkunkle/Resonances/LepLep_mumu_2018_01_21
DIR_ELEL=/data/users/jkunkle/Resonances/LepLep_elel_2018_01_21

DIR_SIGTruth=/data/users/jkunkle/Resonances/SigNoFilt_2018_01_21


all : plots workspaces limits

plots :
	#python MakePlots.py --baseDirMu ${DIR_MU} --baseDirEl ${DIR_EL} --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG}  --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --baseDirMuMu ${DIR_MUMU} --baseDirElEl ${DIR_ELEL} --baseDirNoFilt ${DIR_SIGTruth} --makeZCR --makeWCR --makeSR --makeJetBkg --makeEleBkg --outputDir ${OUTPUTMAIN}/Plots
#	python MakePlots.py --baseDirMu ${DIR_MU} --baseDirEl ${DIR_EL} --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG}  --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --baseDirMuMu ${DIR_MUMU} --baseDirElEl ${DIR_ELEL} --baseDirNoFilt ${DIR_SIGTruth} --makeJetBkg --outputDir ${OUTPUTMAIN}/Plots
	python MakePlots.py --baseDirMu ${DIR_MU} --baseDirEl ${DIR_EL} --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG}  --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --baseDirMuMu ${DIR_MUMU} --baseDirElEl ${DIR_ELEL} --baseDirNoFilt ${DIR_SIGTruth} --makeSR --outputDir ${OUTPUTMAIN}/Plots

wjets : 

	#python MakeWJetsWS.py --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --doClosure  --outputDir ${OUTPUTMAIN}/WJetsWS
	python MakeWJetsWS.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --outputDir ${OUTPUTMAIN}/WJetsWS

signal : 

	python MakeSignalWS.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --outputDir ${OUTPUTMAIN}/SignalWS

mcbkg : 

	python MakeMCBkgWS.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --outputDir ${OUTPUTMAIN}/MCBkgWS

efake : 

	python MakeEFakeWS.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --outputDir ${OUTPUTMAIN}/EFakeWS


workspaces :
	#python MakeWorkspaces.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --doSignal --doWGamma --doTop --doZGamma --doWJets --doEleFake --doClosure --outputDir ${OUTPUTMAIN}/Workspaces
	#python MakeWorkspaces.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId}  --doSignal --doWGamma --doTop --doZGamma  --outputDir ${OUTPUTMAIN}/Workspaces
	#python MakeWorkspaces.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --doWJets --doEleFake --doClosure  --outputDir ${OUTPUTMAIN}/Workspaces
	#python MakeWorkspaces.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --doWGamma --outputDir ${OUTPUTMAIN}/Workspaces
	python MakeWorkspaces.py --baseDirMuG ${DIR_MUG} --baseDirElG ${DIR_ELG} --baseDirMuGNoId ${DIR_MUGNoId} --baseDirElGNoId ${DIR_ELGNoId} --doSignal --outputDir ${OUTPUTMAIN}/Workspaces

fits :
	python MakeFits.py --baseDir ${OUTPUTMAIN}/Workspaces/ --combineDir /home/jkunkle/WG_Analysis/CMSSW_8_1_0 --outputDir ${OUTPUTMAIN} --doVarOptimization   --useHistTemplate  --noRunCombine
	#python MakeFits.py --baseDir ${OUTPUTMAIN}/Workspaces/ --combineDir ~/usercode/Plotting/LimitSetting/CMSSW_7_4_7/ --outputDir ${OUTPUTMAIN} --doVarOptimization  --limitMethod MaxLikelihoodFit
	#python MakeFits.py --baseDir ${OUTPUTMAIN}/Workspaces/ --combineDir ~/usercode/Plotting/LimitSetting/CMSSW_7_4_7/ --outputDir ${OUTPUTMAIN} --doVarOptimization  --limitMethod HybridNew
	#python MakeFits.py --baseDir ${OUTPUTMAIN}/Workspaces/ --combineDir /home/jkunkle/WG_Analysis/CMSSW_8_1_0 --outputDir ${OUTPUTMAIN} --doVarOptimization  --limitMethod Asymptotic   --useToySignal --useToyBackground

signalfits : 
	python MakeResSignalFits.py --baseDirMu  ${DIR_MUG} --baseDirEl ${DIR_ELG} --outputDir ${OUTPUTMAIN}/SignalFits
