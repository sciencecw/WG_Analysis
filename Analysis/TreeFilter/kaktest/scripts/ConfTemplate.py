from core import Filter

def get_remove_filter() :
    """ Define list of regex strings to filter input branches to remove from the output.
        Defining a non-empty list does not apply the filter, 
        you must also supply --enableRemoveFilter on the command line.
        If both filters are used, all branches in keep_filter are used
        except for those in remove_filter """

    return []

def get_keep_filter() :
    """ Define list of regex strings to filter input branches to retain in the output.  
        Defining a non-empty list does not apply the filter, 
        you must also supply --enableKeepFilter on the command line
        If both filters are used, all branches in keep_filter are used
        except for those in remove_filter """

    return []

def config_analysis( alg_list, args ) :
    """ Configure analysis modules. Order is preserved """

    
    # for complicated configurations, define a function
    # that returns the Filter object and append it to the
    # alg list.  Otherwise you can directly append 
    # a Filter object to the list
    # There is no restriction on the naming or inputs to these funtions
    alg_list.append( build_photon( do_cutflow=True, do_hists=True ) )

    filter_event = Filter('FilterEvent')
    filter_event.cut_nPho = ' > 0 '

    alg_list.append( filter_event )

def build_photon( do_cutflow=False, do_hists=False ) :

    # Define the filter object
    # The name must mathch that checked in the ApplyModule
    # function in src/RunAnalysis.cxx
    filt = Filter('BuildPhoton')

    # filter objects have a do_cutflow option
    # enabling this flag keeps a count, for each cut, of
    # the number of times a PassInt, PassFloat, or PassBool
    # is called and how many times it passes
    filt.do_cutflow = do_cutflow

    # define cuts
    # a few examples are given below
    # complicated cuts, such as
    # vetoing the crack region in eta
    # require some special handling
    filt.cut_pt           = ' > 15 '
    filt.cut_abseta       = ' < 2.5'
    filt.cut_abseta_crack = ' > 1.479 & < 1.566 '
    filt.invert('cut_abseta_crack')

    filt.cut_sigmaIEIE = ' < 0.011 '

    if do_hists :
        # The Filter object also retains a list of 
        # histograms.  Each time a PassInt, PassFloat, or PassBool
        # is called for a given cut the histogram is filled
        # with the value input.
        # The histogram name
        # (the first argument to add_hist) must
        # match the cut name for it to work.
        # by default two histograms are made for
        # each cut -- a before and after.
        # The before histogram is filled on each call
        # while the after histogram is filled only
        # if the cut passes
        filt.add_hist( 'cut_pt', 100, 0, 500 )
        filt.add_hist( 'cut_abseta', 50, 0, 5 )
        filt.add_hist( 'cut_abseta_crack', 50, 0, 5 )
        filt.add_hist( 'cut_sigmaIEIE', 50, 0, 0.05 )

    return filt

