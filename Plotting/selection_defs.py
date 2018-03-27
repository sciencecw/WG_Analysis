def get_base_selection( channel ) :

    if channel == 'mu' :
        return 'mu_pt30_n==1 && mu_n==1 && el_n==0 '  # require 1 muon with pt > 30 and 1 muon with pt > 10 (second lepton veto) and 0 electrons with pt > 10 (second lepton veto)
    if channel == 'el' :
        return 'el_pt30_n==1 && el_n==1 && mu_n==0 ' # require 1 electron with pt > 30 and 1 electron with pt > 10 (second lepton veto) and 0 muonswith pt > 10 (second lepton veto)
    if channel == 'mumu' :
        return 'mu_pt30_n>=1 && mu_n==2  '  # require 1 muon with pt > 30 and 2 muons with pt > 10 
    if channel == 'elel' :
        return 'el_pt30_n>=1 && el_n==2 ' # require 1 electron with pt > 30 and 2 electrons with pt > 10
    if channel == 'muel' :
        return 'mu_pt30_n>=1 && el_n==1  '  # require 1 muon with pt > 30 and 2 muons with pt > 10 

def get_weight_str( ) :

    #return ' ( NLOWeight * PUWeight + isData ) '
    return '1'

def get_phid_selection( sel1, sel2='' ) :

    if sel1 == 'all' :
        return 'ph_n'
    if sel1 == 'medium' :
        return 'ph_medium_n' 
    if sel1 == 'chIso' :
        if sel2 == 'sigmaIEIE' :
            return 'ph_mediumNoSIEIENoChIso_n'
        else :
            return 'ph_mediumNoChIso_n'
    if sel1 == 'sigmaIEIE' :
        if sel2 == 'chIso' :
            return 'ph_mediumNoSIEIENoChIso_n'
        else :
            return 'ph_mediumNoSIEIE_n'

    assert( 'get_phid_selection -- Could not parse selection vars!' )

def get_phid_idx( sel1, sel2='' ) :

    if sel1 == 'all' :
        return '0'
    if sel1 == 'medium' :
        return 'ptSorted_ph_medium_idx[0]' 
    if sel1 == 'chIso' :
        if sel2 == 'sigmaIEIE' :
            return 'ptSorted_ph_mediumNoSIEIENoChIso_idx[0]'
        else :
            return 'ptSorted_ph_mediumNoChIso_idx[0]'
    if sel1 == 'sigmaIEIE' :
        if sel2 == 'chIso' :
            return 'ptSorted_ph_mediumNoSIEIENoChIso_idx[0]'
        else :
            return 'ptSorted_ph_mediumNoSIEIE_idx[0]'

    assert( 'get_phid_idx -- Could not parse selection vars!' )

def get_phid_cut_var( ivar ) :

    if ivar == 'chIso' :
        return 'ph_chIsoCorr'
    if ivar == 'sigmaIEIE' :
        return 'ph_sigmaIEIE'

    assert( 'get_phid_cut_var -- Could not parse selection vars!' )



