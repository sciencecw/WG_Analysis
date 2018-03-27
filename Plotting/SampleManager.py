import os
import sys
import math
import re
import imp
import ROOT
import copy
import getpass
import uuid
import itertools
import eos_utilities
import random
from array import array
import time
import analysis_utils

from uncertainties import ufloat
from uncertainties import umath
import pickle
import core
import subprocess
import multiprocessing
import collections

ROOT.gROOT.SetBatch(False)

ROOT.gStyle.SetPalette(1)

class Sample :
    """ Store information about one sample """

    def __init__ ( self, name, **kwargs ) :

        self.name = name
        self.hist = None
        self.loop_hists = []

        # list of open files. Only used if extracting histograms
        self.ofiles = []

        # chain is the TChain for this sample.  
        # Will be created/overwritten if Addfiles is callled.  
        # Should not be filled for a sample group
        self.chain = kwargs.get('chain', None)

        # isActive determines if the sample will be drawn, default=True
        self.init_isActive = kwargs.get('isActive', True) 
        self.isActive      = self.init_isActive
        self.isActiveReq   = self.init_isActive

        # if isData is true the sample will be drawn as points with an error bar, default=False
        self.isData   = kwargs.get('isData', False)

        # if isSignal is true the sample will be drawn as a line and not stacked, default=False
        self.isSignal = kwargs.get('isSignal', False)
        self.lineStyle = kwargs.get('sigLineStyle', 7 )
        self.lineWidth = kwargs.get('sigLineWidth', 2 )

        # if the sample is to be used as a ratio
        self.isRatio = kwargs.get('isRatio', False)

        # if true for any of the samples that are drawn
        # a shaded error band is drawn on the top of the stack
        # using errors from all samples
        self.displayErrBand = kwargs.get('displayErrBand', False)

        # color determines the histogram or line color, default=black
        self.color = kwargs.get('color', ROOT.kBlack)

        # Provide a legend name.  by default the legend name will be the sample name
        self.legendName = kwargs.get('legendName', name)
        if self.legendName is None :
            self.legendName = name

        # if drawRatio is true, the ratio between this and the stack will be
        # drawn in the ratio box, default=False
        self.drawRatio = kwargs.get('drawRatio', False)

        # scale is the weight applied to this sample, default=1.0
        self.scale         = kwargs.get('scale', 1.0)

        # hold the cross section info
        self.cross_section = kwargs.get('cross_section', 0.0)
        self.total_events = kwargs.get('total_events', -1)

        self.averageSamples = kwargs.get('averageSamples', False )

        # groupedSamples is filled if this is a grouping of other samples.
        # if a sample is grouped, it is drawn as the sum of all entries
        # in this list.  Be careful to set isActive to False for the
        # sub samples or they will be drawn as well.  Grouped samples
        # can be further grouped default=[] 
        self.groupedSamples = kwargs.get('groupedSamples', [])

        # set a sample to temporary to delete it
        # when clear_all is called
        self.temporary = False

        self.failed_draw = False

        self.list_of_branches = []

    def SetHist( self, hist=None ) :
        if hist is not None :
            self.hist = hist
        self.InitHist()

    def InitHist(self) :
        self.hist.SetLineColor( self.color )
        self.hist.SetMarkerColor( self.color )
        self.hist.SetTitle('')
        self.hist.Scale( self.scale )
        print 'Scale %s by %f' %( self.name, self.scale )
        if self.isData :
            self.hist.SetMarkerStyle( 20 )
            self.hist.SetMarkerSize( 1.15 )
            self.hist.SetStats(0)
        if self.isSignal :
            self.hist.SetLineWidth(self.lineWidth)
        if self.isRatio :
            self.hist.SetMarkerStyle( 10 )
            self.hist.SetMarkerSize( 1.1 )
            #self.hist.SetNdivisions(509, True )


    def AddFiles( self, files, treeName=None, readHists=False ) :
        """ Add one or more files and grab the tree named treeName """

        if not isinstance(files, list) :
            files = [files]

        if treeName is not None :
            self.chain = ROOT.TChain(treeName, self.name)
            for file in files :
                self.chain.Add(file)

            self.chain.SetBranchStatus('*', 0 )

        if readHists :
            for file in files :
                self.ofiles.append( ROOT.TFile.Open( file ) )

    def AddGroupSamples( self, samples ) :
        """ Add subsamples to this sample """

        if not isinstance( samples, list) :
            samples = [samples]

        self.groupedSamples += samples

    def IsGroupedSample(self) :
        return ( len( self.groupedSamples ) > 0 )


    def enable_parsed_branches( self, brstr ) :
        if self.chain is not None :
            for br in self.chain.GetListOfBranches() :
                if brstr.count( br.GetName() ) > 0 and self.chain.GetBranchStatus( br.GetName() ) == 0 :
                    self.chain.SetBranchStatus( br.GetName(), 1)

    def get_list_of_branches(self) :
        if self.list_of_branches :
            return self.list_of_branches
        else :
            branches = []
            if self.chain is not None :
                for br in self.chain.GetListOfBranches() :
                    branches.append(br.GetName())

            self.list_of_branches = branches

            return branches

    def reset_status( self ) :
        """ Reset sample to the initial status """
        self.isActive = self.isActiveReq

    def getLineStyle( self ) :
        return self.lineStyle

class LegendConfig :

    def __init__(self) : 
        pass

class LabelConfig :

    def __init__(self) : 
        pass

class HistConfig :

    def __init__(self) : 
        pass

class DrawConfig :
    """ Store and process all informaiton necessary for making a histogram """

    used_names = []

    def __init__(self, var, selection, histpars, samples=None, hist_config={}, legend_config={}, label_config={}, replace_selection_for_sample={}) :

        self.var = var
        self.selection = selection

        if not isinstance( self.var, list ) :
            self.var = [self.var]
        if not isinstance( self.selection, list ) :
            self.selection = [self.selection]

        self.samples = samples
    
        if self.samples is not None and not isinstance( self.samples, list ) :
            self.samples = [self.samples]

        self.histpars = histpars
        self.hist_config = hist_config
        self.legend_config = legend_config
        self.label_config = label_config
        self.replace_selection_for_sample = replace_selection_for_sample

        self.modified_selection = None
        self.stack_save_params = {}
        self.stack_dump_params = {}

        self.compare_hists = False
        self.no_auto_draw = False

        self.hist_configs = collections.OrderedDict()

        self.samp_man_id = None

    def doRatio(self) :
        return self.hist_config.get('doratio', False)

    def get_ylabel(self ) :
        #FIX
        ylabel = self.hist_config.get('ylabel', None) 
        if ylabel is None :
            if isinstance( self.histpars, tuple ) :
                bin_width = ( self.histpars[2] - self.histpars[1] )/self.histpars[0]
                bin_width_f = ( self.histpars[2] - self.histpars[1] )/float(self.histpars[0])
            else :
                bin_width = 1
                bin_width_f = 1
                

            if math.fabs(bin_width_f - bin_width) != 0 :
                ylabel = 'Events / %.1f GeV' %bin_width_f
            else :
                ylabel = 'Events / %d GeV' %bin_width

        return ylabel

    def get_xlabel(self) :

        return self.hist_config.get('xlabel', '')

    def get_rlabel(self) :
        rlabel = self.hist_config.get('rlabel', None) 
        if rlabel is None :
            rlabel = 'Data / MC'

        return rlabel

    def get_tick_x_format(self) :
        return self.hist_config.get( 'ticks_x', None )
    def get_tick_y_format(self) :
        return self.hist_config.get( 'ticks_y', None )
    

    def get_legend_entries(self) :

        legend_entries = self.legend_config.get('legend_entries', [])
        if len( legend_entries)  != len(self.samples) :
            legend_entries = self.samples

        return legend_entries

    def get_ymin( self ) :
        return self.hist_config.get('ymin', None)
    def get_ymax( self ) :
        return self.hist_config.get('ymax', None)
    def get_ymax_scale( self ) :
        return self.hist_config.get('ymax_scale', None)
    def get_rmin( self ) :
        return self.hist_config.get('rmin', 0 )
    def get_rmax( self ) :
        return self.hist_config.get('rmax', 2 )
    def get_logy( self ) :
        return self.hist_config.get('logy', False )
    def get_normalize( self ) :
        return self.hist_config.get('normalize', False )

    def save_stack( self, filename, dirname, canname ) :

        self.stack_save_params['filename'] = filename
        self.stack_save_params['dirname'] = dirname
        self.stack_save_params['canname'] = canname

    def dump_stack( self, filename, dirname ) :

        self.stack_dump_params['filename'] = filename
        self.stack_dump_params['dirname'] = dirname

    def add_label( self, config ) :
        if 'extra_label' not in self.label_config :
            self.label_config


    def get_labels( self ) :

        labels=[]
        
        stattext = self.label_config.get('statsLabel', None )
        if stattext is not None :
            statlabel  = ROOT.TLatex()
            statlabel.SetTextFont( 42 )
            statlabel   .SetNDC()
            statlabel  .SetTextSize(0.045)
            statlabel.SetText(0.16, 0.93, stattext)

            labels.append(statlabel)

        labelStyle = self.label_config.get('labelStyle', None)
        labelLoc = self.label_config.get('labelLoc', None)
        if labelStyle is None:

            text_x = 0.18
            text_y = 0.93

            if labelLoc == 'topright' :
                text_x = 0.75
                text_y = 0.93
            cmslabel = ROOT.TLatex()
            cmslabel.SetNDC()
            cmslabel.SetTextSize( 0.04 )
            cmslabel.SetText(text_x, text_y, 'CMS Internal')

            labels.append(cmslabel)

        elif labelStyle.count('fancy') :
            extText = 'Internal'
            if labelStyle.count('prelim') :
                extText = 'Preliminary'

            labeltext = '19.4 fb^{-1} (8 TeV)'
            if labelStyle.count('13') :
                labeltext = '36 fb^{-1} (13 TeV)'

            rootslabel = ROOT.TLatex()
            cmslabel = ROOT.TLatex()
            extlabel = ROOT.TLatex()

            rootslabel.SetTextFont(42)
            cmslabel.SetTextFont( 61 )
            extlabel.SetTextFont(52)

            extlabel  .SetNDC()
            rootslabel .SetNDC()
            cmslabel   .SetNDC()

            if not labelStyle.count( 'paper') :
                extlabel  .SetTextSize(0.045)
            rootslabel .SetTextSize(0.045)
            cmslabel  .SetTextSize(0.055)

            extlabel.SetText( 0.25, 0.93, extText )
            #rootslabel.SetText(0.65, 0.93, '#font[132]{#sqrt{s} = 8 TeV, L = 19.4 fb^{-1} }' )

            rootslabel.SetText(0.73, 0.93, labeltext  )

            if not labelStyle.count('paper') :
                labels.append(extlabel)
            labels.append(cmslabel)
            labels.append(rootslabel)

        extra_label = self.label_config.get( 'extra_label', None )
        if extra_label is not None :

            if not isinstance( extra_label, list) :
                extra_label = [extra_label]

            for lab in extra_label :
                extra_label_str = '#font[132]{'+lab+'}'
                extra_label_loc = self.label_config.get( 'extra_label_loc', None )
                labels.append(self.place_extra_label( lab, extra_label_loc ))

        return labels

    #--------------------------------
    def place_extra_label(self, text, location=None) :

        label = ROOT.TLatex()
        label.SetNDC()
        label.SetTextSize( 0.045 )
        xval = 0.6
        yval = 0.7
        if location is None : 
            print 'Please give a location for the label'
        elif isinstance(location, tuple) : 
            xval = location[0]
            yval = location[1]
        elif location == 'TopLeft' :
            xval = 0.15
            yval = 0.85

        elif location == 'BottomLeft' :
            xval = 0.25
            yval = 0.25
        else :
            xval = 0.6
            yval = 0.7

        label.SetText(xval, yval, text)
        return label

    def get_unique_name( self, var ) :

        outname = ''
        
        basename = var

        # first remove brackets
        pos_begin = 100
        while pos_begin >= 0 :
            pos_begin = basename.rfind( '[' )
            if pos_begin < 0 :
                break
            pos_end = basename.find( ']', pos_begin )
            basename = basename[:pos_begin] + basename[pos_end+1:]
        if basename.count('+') :
            basename = basename.split('[')[0]
            basename = basename.replace('+', '_')
        if basename.count('(') :
            basename = basename.replace('(', '_')
        if basename.count(')') :
            basename = basename.replace(')', '_')
        if basename.count(':') :
            basename = basename.replace(':', '_')

        if basename in self.used_names :
            for i in range( 0, 1000000  ) :
                outname = '%s_%d' %(basename, i)
                if outname not in self.used_names :
                    self.used_names.append(outname)
                    break
        else :
            outname = basename
            self.used_names.append(basename)

        return outname


    def get_var_val(self, sample, treename) :
        mod_var = self.var
        all_branches = sample.get_list_of_branches()
        for br in all_branches :
            if mod_var.count( br ) and not mod_var.count(treename+'.'+br)  :
                mod_var = mod_var.replace( br, treename+'.'+br )

        return mod_var

    def compile_selection_string( self, sample, treename ) :
        eval_str = self.get_eval_selection_string( sample, treename )
        self.compiled_selection_string = compile( eval_str, '<string>', 'eval')

    def compile_var_string( self, sample, treename ) :
        var_str = self.get_var_val( sample, treename )
        self.compiled_var_str = compile( var_str, '<string>', 'eval')

    def get_compiled_selection_string( self ) :
        return self.compiled_selection_string

    def get_compiled_var( self ) :
        return self.compiled_var_str

    def get_names( self ) :
        return self.hist_configs.keys()

    def get_hist_type( self ) :

        if self.var[0].count(':') == 2 :
            return 'TH3F'
        elif self.var[0].count(':') == 1 :
            return 'TH2F'
        else :
            return 'TH1F'

    def get_hist_declarations( self ) :

        hist_decs = []

        for name in self.hist_configs.keys() :

            if type( self.histpars ) is tuple : 
                if self.var[0].count(':') == 1 : # 2-D histogram
                    if len(self.histpars) == 2 and type( self.histpars[0] ) is list and type(self.histpars[1]) is list : #both axes are variably binned
                        text = 'double %sxarr[%d] = {'%(name, len(self.histpars[0])) + ','.join( [str(x) for x in self.histpars[0]] ) + '}; \n '
                        text += 'double %syarr[%d] = {'%(name, len(self.histpars[1])) + ','.join( [str(y) for y in self.histpars[1]] ) + '}; \n '
                        text += r' hist_%s = new TH2F( "%s", "", %d, %sxarr, %d, %syarr );' %( name, name, len(self.histpars[0])-1, name, len(self.histpars[1])-1, name ) 
                        hist_decs.append(text)
                    else :
                        if len(self.histpars) != 6 :
                            print 'varable expression requests a 2-d histogram, please provide 6 hist parameters, nbinsx, xmin, xmax, nbinsy, ymin, ymax'
                            return
                        text = r' hist_%s = new TH2F( "%s", "", %d, %f, %f, %d, %f, %f );' %( name, name, self.histpars[0], self.histpars[1], self.histpars[2], self.histpars[3], self.histpars[4], self.histpars[5]  ) 
                        hist_decs.append(text)
                elif self.var[0].count(':') == 2 and not self.var[0].count('::') : # make a 3-d histogram
                    if len(self.histpars) != 9 :
                        print 'varable expression requests a 3-d histogram, please provide 9 hist parameters, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax'
                        return
                    text = r' hist_%s = new TH3F( "%s", "", %d, %f, %f, %d, %f, %f, %d, %f, %f );' %( name, name, self.histpars[0], self.histpars[1], self.histpars[2], self.histpars[3], self.histpars[4], self.histpars[5], self.histpars[6], self.histpars[7], self.histpars[8]  ) 
                    hist_decs.append(text)
                else : # 1-d histogram
                    text = r' hist_%s = new TH1F( "%s", "", %d, %f, %f );' %( name, name, self.histpars[0], self.histpars[1], self.histpars[2] ) 
                    hist_decs.append(text)

            elif type( self.histpars ) is list : # variable rebinning
                text = 'double %sxarr[%d] = {'%(name, len(self.histpars)) + ','.join( [str(x) for x in self.histpars] ) + '}; \n '
                text += r' hist_%s = new TH1F( "%s", "", %d, %sxarr );' %( name, name, len(self.histpars)-1, name ) 
                hist_decs.append(text)

            else :
                print 'No histogram parameters were passed'

        return hist_decs

    def create_hist_configs( self, branches=None ) :

        if branches is None :
            if len( self.var ) == 1 : # make one name for each selection
                if self.samples and ( len(self.selection) == len( self.samples ) ) :
                    if len(self.hist_config.get('colors', [])) != len( self.selection) :
                        self.hist_config['colors'] = [ROOT.kBlack]*len(self.selection)
                    for samp, sel, color, leg_entry in zip(self.samples, self.selection, self.hist_config['colors'], self.get_legend_entries() ) :
                        name = self.get_unique_name( self.var[0] ) 
                        self.hist_configs[name] = {'var' : self.var[0], 'selection' : sel, 'sample' : samp, 'color' : color, 'legend_entry' : leg_entry} 

                else :
                    for sel in self.selection :
                        name = self.get_unique_name( self.var[0] ) 
                        self.hist_configs[name] = {'var' : self.var[0], 'selection' : sel} 
            else : 
                if self.samples and ( len(self.selection) == len( self.samples ) and len( self.var ) == len( self.samples) ) :
                    if len(self.hist_config.get('colors', [])) != len( self.selection) :
                        self.hist_config['colors'] = [ROOT.kBlack]*len(self.selection)
                for var, samp, sel, color, leg_entry in zip( self.var, self.samples, self.selection, self.hist_config['colors'], self.get_legend_entries() ) :
                    name = self.get_unique_name( self.var[0] ) 
                    self.hist_configs[name] = {'var' : var, 'selection' : sel, 'sample' : samp, 'color' : color, 'legend_entry' : leg_entry} 

                print 'Case when multiple vars is used is not implemented'
        else :

            if len( self.var ) == 1 : # make one name for each selection
                if self.samples and ( len(self.selection) == len( self.samples ) ) :
                    if len(self.hist_config.get('colors', [])) != len( self.selection) :
                        self.hist_config['colors'] = [ROOT.kBlack]*len(self.selection)
                    for samp, sel, color, leg_entry in zip(self.samples, self.selection, self.hist_config['colors'], self.get_legend_entries() ) :
                        name = self.get_unique_name( self.var[0] ) 
                        var = self.get_cpp_var_str( self.var[0], branches ) 
                        selection = self.get_cpp_selection_str( sel, branches )
                        self.hist_configs[name] = {'var' : self.var[0], 'selection' : sel, 'cppvar':var, 'cppselection':selection, 'sample' : samp, 'color' : color, 'legend_entry' : leg_entry} 
                else :

                    for sel in self.selection :
                        name = self.get_unique_name( self.var[0] ) 
                        var = self.get_cpp_var_str( self.var[0], branches ) 
                        selection = self.get_cpp_selection_str( sel, branches ) 
                        self.hist_configs[name] = {'var' : self.var[0], 'selection' : sel, 'cppvar':var, 'cppselection':selection} 
    
            else : #unclear if this case exists, don't implement for now
                print 'Case when multiple vars is used is not implemented'


    def get_cpp_selection_strs( self ) :
        return [v['cppselection'] for v in self.hist_configs.values()]

    def get_cpp_selection_str( self, selection, branches ) :
        """ Create selection string in c++ 

            branches is a list of branches in the tree
            such that we can search the selection string
            for occurences of these branches
        """

        # for each branch, find all occurances of that branch
        # and store the branch name and end of the name
        # indexed by the start index
        # Then make the replacement for the largest
        # string.  This avoids problems when one branch
        # is a substring of another and the incorrect
        # replacement is made
        matched_locations = {}
        for br in branches :
            for selitr in re.finditer( br['name'], selection ) :
                start = selitr.start()
                end   = selitr.end()
                br_range = (start, end)
                matched_locations.setdefault( br_range , [] ).append( br['name'] )

        max_ranges = {}
        ranges = list( matched_locations.keys() )
        for idx1, range1 in enumerate(ranges) :
            has_sub_range = False
            for idx2, range2 in enumerate(ranges) :
                if idx1 == idx2 : 
                    continue
                if range1[0] == range2[0] : #they start at the same place
                    if range1[1] <= range2[1] :
                        has_sub_range = True
                if range1[1] == range2[1] : # they end at the same place 
                    if range1[0] >= range2[0] :
                        has_sub_range = True
                if range1[0] > range2[0] and range1[1] < range2[1] :
                    has_sub_range = True

            if not has_sub_range :
                max_ranges[range1] = matched_locations[range1]

        modified_selection = str(selection)
        # start from the end of the string (max start)
        # so that the found indices don't change
        # when the string is modified
        all_start = max_ranges.keys()
        all_start.sort()
        for max_range in reversed( all_start ) :
            max_br = max_ranges[max_range]
            if len( max_br ) != 1 :
                print 'Matched multiple branches with the max range!  This must be fixed!'
                sys.exit(1)

            modified_selection = modified_selection[:max_range[0]] + 'IN::' + max_br[0] + modified_selection[max_range[1]:]

        #for br in branches :
        #    if modified_selection.count(br['name']) and not modified_selection.count( 'IN::'+br['name']):
        #        modified_selection = modified_selection.replace( br['name'], 'IN::'+br['name'])
        # a bit hacked
        start_bracket = 100
        while start_bracket > 0 :
            start_bracket = modified_selection.rfind('[')
            if start_bracket < 0 :
                break
            end_bracket = modified_selection.find( ']', start_bracket )
            
            # put an at->( in place of the first bracket
            modified_selection = modified_selection[:start_bracket] + '->at(' + modified_selection[start_bracket+1:end_bracket] + ')' + modified_selection[end_bracket+1:]

        return modified_selection

    def get_cpp_var_strs(self ) :
        return [v['cppvar'] for v in self.hist_configs.values()]

    def get_cpp_var_str(self, var, branches) :

        # for a 1-D histogram the var can undergo the
        # same modification as the selection str
        # for 2-D and 3-D histograms the colon separating
        # the variables should become a comma so that
        # the Fill command works.  However the ordering
        # of the var string and the Fill command is reversed
        # eg for a 3-D histogram the var is varZ:varY:varX
        # but the fill is Fill( varX, varY, varZ)

        input_var = var
        if var.count( ':' ) :
            split_var = var.split( ':' )
            split_var.reverse()
            input_var = ','.join( split_var )

        var_str = self.get_cpp_selection_str( input_var, branches )
        return var_str

        ## for each branch, find all occurances of that branch
        ## and store the branch name and end of the name
        ## indexed by the start index
        ## Then make the replacement for the largest
        ## string.  This avoids problems when one branch
        ## is a substring of another and the incorrect
        ## replacement is made
        #matched_locations = {}
        #for br in branches :
        #    for selitr in re.finditer( br['name'], var ) :
        #        start = selitr.start()
        #        matched_locations.setdefault( start , {} )
        #        matched_locations[start][selitr.end()-start] = br['name']

        #modified_var = str(var)
        ## start from the end of the string (max start)
        ## so that the found indices don't change
        ## when the string is modified
        #all_start = matched_locations.keys()
        #all_start.sort()
        #for start in reversed( all_start ) :
        #    matches = matched_locations[start]
        #    max_match  = max( matches.keys() )
        #    max_br = matches[max_match]

        #    modified_var = modified_var[:start] + 'IN::' + max_br + modified_var[start+max_match:]

        ##modified_var = var
        ##for br in branches :
        ##    if modified_var.count(br['name']) and not modified_var.count( 'IN::'+br['name']):
        ##        modified_var= modified_var.replace( br['name'], 'IN::'+br['name'])

        ## a bit hacked
        #for i in range(0, 10) :
        #    modified_var= modified_var.replace('[%s]'%i, '->at(%s)'%i )

        #return modified_var

    
    def get_eval_selection_string(self, sample, treename) :

        if self.modified_selection is not None :
            return self.modified_selection
        else :

            self.modified_selection = self.selection

            # append treename to all identified branches
            all_branches = sample.get_list_of_branches()
            for br in all_branches :
                if self.modified_selection.count(br) and not self.modified_selection.count( treename+'.'+br):
                    self.modified_selection = self.modified_selection.replace( br, treename+'.'+br)

            if self.modified_selection.count('&&') :
                self.modified_selection = self.modified_selection.replace( '&&', 'and')

            if self.modified_selection.count('||') :
                self.modified_selection = self.modified_selection.replace( '||', 'or')

            if self.modified_selection.count('fabs') :
                self.modified_selection = self.modified_selection.replace( 'fabs', 'math.fabs')

            return self.modified_selection

    def get_selection_string( self, name ) :
        res = self.replace_selection_for_sample.get( name, None )
        if res is not None :
            return res
        else :
            return self.selection[0]

    def init_hist( self, name ) :

        hist = None
        histname = str(uuid.uuid4())

        if type( self.histpars ) is tuple :
            if self.var[0].count(':') == 1 : 
                if len(self.histpars) == 2 and type( self.histpars[0] ) is list and type(self.histpars[1]) is list :
                    hist = ROOT.TH2F( histname, '', len(self.histpars[0])-1, array('f', self.histpars[0]), len(self.histpars[1])-1, array('f', self.histpars[1]) )
                else :
                    if len(self.histpars) != 6 :
                        print 'varable expression requests a 2-d histogram, please provide 6 hist parameters, nbinsx, xmin, xmax, nbinsy, ymin, ymax'
                        return
                    hist = ROOT.TH2F( histname, '', self.histpars[0], self.histpars[1], self.histpars[2], self.histpars[3], self.histpars[4], self.histpars[5])
            elif self.var[0].count(':') == 2 and not self.var.count('::') : # make a 3-d histogram
                if len(self.histpars) != 9 :
                    print 'varable expression requests a 3-d histogram, please provide 6 hist parameters, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax'
                    return
                hist= ROOT.TH3F( histname, '',self.histpars[0], self.histpars[1], self.histpars[2], self.histpars[3], self.histpars[4], self.histpars[5], self.histpars[6], self.histpars[7], self.histpars[8] )
            else : # 1-d histogram

                hist= ROOT.TH1D( histname, '', self.histpars[0], self.histpars[1], self.histpars[2])

        elif type( self.histpars ) is list :
            hist = ROOT.TH1D( histname, '', len(self.histpars)-1, array('f', self.histpars))
        else :
            print 'No histogram parameters were passed'

        if hist is not None :
            hist.SetTitle( name )
            hist.Sumw2()

        return hist


    
class SampleManager :
    """ Manage input samples and drawn histograms """

    def __init__(self, base_path, treeName=None, mcweight=1.0, treeNameModel='events', filename='ntuple.root', base_path_model=None, xsFile=None, lumi=None, readHists=False, quiet=False) :

        #
        # This plotting module assumes that root files are 
        # organized under sample directories.  The directories
        # that are read are configured through the module passed
        # to ReadSamples
        #
        # All drawn objects are kept in SampleManager in variables
        # starting with curr_ 
        #

        #
        # path to directory containing samples
        # the samples that are read are configured through
        # the input module
        #
        self.base_path       = base_path 

        # the name of the tree to read
        self.treeName        = treeName

        # Name of the file.  This can be overwritten in the configuration module
        self.fileName        = filename

        #the name of the tree to read for models
        self.treeNameModel   = treeNameModel

        #
        # path to directory containing samples for models
        # the samples that are read are configured through
        # the input module
        #
        self.base_path_model = base_path_model

        
        # weight that is applied to all
        # samples not labeled as Data
        self.mcweight        = mcweight

        # store all samples
        self.samples               = []

        # store model samples
        self.modelSamples          = []

        # store the order that the samples were added
        # in the configuration module.  The samples
        # are stacked in this order
        self.stack_order                 = []

        # if the cross section file is given, open it
        # and grab the cross section map out of it
        self.weightMap = analysis_utils.read_xsfile( xsFile, lumi, print_values=True )

        self.curr_hists            = {}
        self.curr_canvases         = {}
        self.curr_stack            = None
        self.curr_legend           = None

        self.legendLimits          = None

        self.legendLoc             = 'Nominal'
        self.legendCompress        = 1.0
        self.legendWiden           = 1.0
        self.legendTranslateX      = 0.0
        self.legendTranslateY      = 0.0
        self.entryWidth            = 0.06
        # Save any plot decorations such as labels here
        # This guarantees that the objects stay in memory
        self.curr_decorations      = []

        # When variable binning is used the binning is stored
        # so that it can be used in future calls
        self.binning = None

        # read histograms instead of trees
        self.readHists = readHists

        self.quiet = quiet

        self.transient_data = {}

        self.samples_conf=None

        self.draw_commands=[]

        self.collect_commands=False

        # make a unique ID
        # this is used to indicate
        # to which SampleManger
        # a DrawConfig was associated to
        self.id = str(uuid.uuid4())

        # keep track if a sample group has been added
        self.added_sample_group=False

            
    #--------------------------------
    def create_sample( self, name, **kwargs ) :

        if name in self.get_sample_names() :
            print 'Sample with name %s already exists' %name
            return None

        new_sample = Sample( name=name )
        new_sample.hist = kwargs.pop('hist', None)
        for arg, val in kwargs.iteritems() :
            if hasattr( new_sample, arg ) :
                setattr( new_sample, arg, val )

        if new_sample.hist is not None :
            new_sample.SetHist( )

        self.samples.append(new_sample)
        return new_sample

    #--------------------------------
    def clone_sample( self, oldname, newname, **kwargs ) :

        if isinstance( oldname, Sample ) :
            oldsample = oldname
        else :
            _oldsample = self.get_samples( name=oldname )

            if _oldsample :
                oldsample = _oldsample[0]
            else :
                print 'Could not locate old sample'
                print 'new name = ', newname
                print 'old name = ', oldname
                return None

        newsample = copy.copy( oldsample )
        newsample.name = newname
        
        histval = kwargs.pop('hist', None)
        for arg, val in kwargs.iteritems() :
            if hasattr( newsample, arg ) :
                setattr( newsample, arg, val )

        if histval is not None :
            newsample.SetHist( histval )
        else :
            newsample.hist=None

        self.samples.append(newsample)

        return newsample

    #--------------------------------
    def create_ratio_sample( self, name, num_sample, den_sample, color=ROOT.kBlack ) :

        if name in self.get_sample_names() :
            print 'Sample %s already exists!  Will not create!' %name
            return None

        if not isinstance( num_sample, Sample ) :
            num_sample_list = self.get_samples( name=num_sample )
            if not num_sample_list :
                print 'create_ratio_sample - ERROR : Numerator sample does not exist, %s' %num_sample
                return None
            num_sample = num_sample_list[0]

        if not isinstance( den_sample, Sample ) :
            den_sample_list = self.get_samples( name=den_sample )
            if not den_sample_list :
                print 'create_ratio_sample - ERROR : Denominator sample does not exist'
                return None
            den_sample = den_sample_list[0]

        ratio_hist = num_sample.hist.Clone( name )
        ratio_hist.Divide( den_sample.hist )

        ratio_hist.SetMarkerStyle(20)
        ratio_hist.SetMarkerSize(1.1)
        #ratio_hist.SetNdivisions(509, True )
        ratio_hist.SetStats(0)
        ratio_hist.SetTitle('')

        return self.create_sample( name=name, isRatio=True, hist=ratio_hist, temporary=True, color=color )

    #--------------------------------
    def config_legend(self, **kwargs ) :

        config = {}

        config['legendLoc']        = kwargs.pop('legendLoc'        , 'Nominal')
        config['legendCompress']   = kwargs.pop('legendCompress'   , 1.0)
        config['legendWiden']      = kwargs.pop('legendWiden'      , 1.0)
        config['legendTranslateX'] = kwargs.pop('legendTranslateX' , 0.0)
        config['legendTranslateY'] = kwargs.pop('legendTranslateY' , 0.0)
        config['entryWidth']       = kwargs.pop('entryWidth' , 0.06)

        for key, val in kwargs.iteritems() :
            config[key] = val
        
        return config

    #--------------------------------
    def apply_lenged_conf(self, config ) :

        if config is None :
            return

        for key, val in config.iteritems() :
            if hasattr( self, key ) :
                setattr(self, key, val)

    #--------------------------------
    def get_samples(self, **kwargs) :

        # if no arguments provided, return all samples
        if not kwargs :
            return self.samples

        # collect results for each argument provided
        # then require and AND of samples matching each
        # criteria
        each_results = []
        for arg, val in kwargs.iteritems() :
            val_list = val
            if not isinstance(val_list, list) :
                val_list = [val_list]
            if hasattr( Sample(''), arg ) :
                each_results.append([ samp  for samp in self.samples if getattr( samp, arg ) in val_list])

        common_results = list( reduce( lambda x,y : set(x) & set(y), each_results ) ) 
        if not common_results :
            print 'WARNING : Found zero samples matching criteria!  Sample matching criteria were : '
            #assert( '' )
            print kwargs
            #for s in self.get_samples() :
            #    print s.name
            return []

        return common_results
        ## default, return all samples
        #if names is None :
        #    return self.samples

        ## empty list provided, return empty list
        #if not names :
        #    return []

        #if not isinstance(names, list) :
        #    names = [names]

        #return filter( lambda x : x.name in names, self.samples )

    #--------------------------------
    def get_model_samples(self, names=[]) :
        if not isinstance(names, list) :
            names = [names]

        if names :
            return filter( lambda x : x.name in names, self.modelSamples )
        else :
           return self.modelSamples

    #--------------------------------
    def get_signal_samples(self) :
        return filter( lambda x : x.isSignal and x.isActive , self.samples )

    #--------------------------------
    def get_sample_names(self) :
        return [x.name for x in self.samples]
    
    #--------------------------------
    def get_model_sample_names(self) :
        return [x.name for x in self.modelSamples]


    #--------------------------------
    def add_decoration(self, obj) :
        self.curr_decorations.append(obj)

    #--------------------------------
    def activate_sample(self, samp_name) :
        
        name = samp_name
        if not isinstance( samp_name, str ) :
            name = samp_name.name

        sel_samps = self.get_samples(name=name)
        if not sel_samps :
            print 'No sample with name %s' %name
        elif len(sel_samps) > 1 :
            print 'Located multiple samples with name %s' %name
        else :
            print 'Activate sample %s' %sel_samps[0].name
            sel_samps[0].isActiveReq=True
            sel_samps[0].isActive=True

    #--------------------------------
    def deactivate_sample(self, samp_name) :
        sel_samps = self.get_samples(name=samp_name)
        if not sel_samps :
            print 'No sample with name %s' %samp_name
        elif len(sel_samps) > 1 :
            print 'Located multiple samples with name %s' %samp_name
        else :
            print 'Deactivate sample %s' %sel_samps[0].name
            sel_samps[0].isActiveReq=False
            sel_samps[0].isActive=False

    #--------------------------------
    def deactivate_all_samples(self) :
        sel_samps = self.get_samples()
        for samp in sel_samps :
            samp.isActive=False
            samp.isActiveReq=False

    #--------------------------------
    def clear_all(self) :
        """ clear all objects """

        self.clear_hists()

        for can in self.curr_canvases.values() :
            if can != None :
                can.Close()
        self.curr_canvases         = {}

        #if self.curr_stack != None :
        ##    self.curr_stack.Delete()
        #    self.curr_stack.Clear()
        #    self.curr_stack.Delete()
        if self.curr_legend != None :
            #self.curr_legend.Delete()
            self.curr_legend = None

        self.curr_decorations      = []

        self.legendLoc='Nominal'
        self.legendCompress=1.0
        self.legendWiden=1.0
        self.legendTranslateX=0.0
        self.legendTranslateY=0.0
        self.entryWidth = 0.052

        self.transient_data= {}
        self.stored_command=''

        for samp in self.samples :
            samp.reset_status()


    #--------------------------------
    def clear_hists(self) :
        rm_samples = []
        for idx, samp in enumerate(self.samples) :
            samp.hist=None
            if samp.temporary :
                print 'removing sample %s' %samp.name
                rm_samples.append(samp)

        for samp in rm_samples:
            if samp.hist is not None :
                samp.hist.Delete()
            self.samples.remove(samp)

        for samp in self.modelSamples :
            samp.hist=None

    #--------------------------------
    def clear_samples(self) :
        for samp in self.samples :
            if samp.chain is not None :
                for fileobj in samp.chain.GetListOfFiles() :
                    file = ROOT.TFile(fileobj.GetTitle())
                    file.Close()
    
        self.samples = []

    #--------------------------------
    def get_grouped_sample_names(self) :
        names = []
        for samp in self.samples :
            if samp.IsGroupedSample() :
                names.append(samp.name)

        return names
            

    #--------------------------------

    def get_sample_branch_minmax( self, samp, branch ) :

        print 'Get %s, %s' %( samp, branch)

        if isinstance( samp, str ) :
            samp = self.get_samples(name=samp)[0]

        if samp.groupedSamples :
            min_val_tot = None
            max_val_tot = None
            for subsampname in samp.groupedSamples :
                min_val, max_val = self.get_sample_branch_minmax(subsampname, branch )
                if min_val_tot is None or min_val_tot < min_val :
                    min_val_tot = min_val
                if max_val_tot is None or max_val_tot > max_val :
                    max_val_tot = max_val
            return (min_val_tot, max_val_tot )

        else :

            samp.chain.SetBranchStatus(branch, 1)
            max_val = samp.chain.GetMaximum(branch)
            min_val = samp.chain.GetMinimum(branch)

            return (min_val, max_val)

    #--------------------------------
    def add_temp_sample(self, samp) :
        samp.temporary = True
        self.samples.append(samp)

    #--------------------------------
    def get_stack_order( self ) :

        out_order = []

        samps = self.get_samples( name=self.stack_order )

        for samp in samps :
            if not samp.isData and not samp.isSignal and samp.isActive :
                out_order.append( samp.name )

        return out_order

    #---------------------------------------
    def GetLowestGroupedSamples( self, sample ) :
        lowest = []
        for subsamp in self.get_samples(name=sample.groupedSamples ) :
            if subsamp.IsGroupedSample() :
                lowest += self.GetLowestGroupedSamples( subsamp )
            else :
                lowest.append(subsamp)
        return lowest


    #---------------------------------------
    def start_command_collection( self ) :
        self.collect_commands = True
        for sample in self.samples :
            sample.loop_hists=[]
        

    #---------------------------------------
    def add_draw_config( self, varexp, selection, histpars, hist_config={}, label_config={}, legend_config={}, replace_selection_for_sample={}  ) :
        self.draw_commands.append( DrawConfig( varexp, selection, histpars, hist_config=hist_config, label_config=label_config, legend_config=legend_config, replace_selection_for_sample=replace_selection_for_sample ) )

        self.draw_commands[-1].samp_man_id = self.id

        return self.draw_commands[-1]

    #---------------------------------------
    def add_compare_config( self, varexp, selection, samples, histpars, hist_config={}, label_config={}, legend_config={}, replace_selection_for_sample={}  ) :
        self.draw_commands.append( DrawConfig( varexp, selection, histpars, samples=samples, hist_config=hist_config, label_config=label_config, legend_config=legend_config, replace_selection_for_sample=replace_selection_for_sample ) )

        self.draw_commands[-1].compare_hists=True

        self.draw_commands[-1].samp_man_id = self.id

        return self.draw_commands[-1]

    #---------------------------------------
    def config_and_queue_hist( self, samp, var, sel, binning ) :

        self.activate_sample( samp )
        draw_conf = DrawConfig( var, sel, binning, samples=samp )
        self.queue_draw( draw_conf)
        return draw_conf

    #---------------------------------------
    def queue_draw( self, draw_config ) :
        self.draw_commands.append( draw_config )

        self.draw_commands[-1].samp_man_id = self.id

        self.draw_commands[-1].no_auto_draw=True

    #---------------------------------------
    def add_save_stack( self, filename, outputDir, canname=None) :
        if canname is None :
            canname = 'base'
        self.draw_commands[-1].save_stack( filename, outputDir, canname )

    #---------------------------------------
    def add_dump_stack( self, filename, outputDir) :
        self.draw_commands[-1].dump_stack( filename, outputDir )

    #---------------------------------------
    def run_commands( self, nsplit=0, nFilesPerJob=5, treeName='' ) :

        self.collect_commands = False

        #------------------------------------------------------------
        # 1) Write the code that loads braches, loops, and fills hists
        # 2) compile and run the code
        # 3) Loop over the draw configs and make and save canvases
        #------------------------------------------------------------

        workarea = os.getenv('WorkArea')
        user = getpass.getuser()

        #compile_base = '%s/../Plotting/compiled_code' %workarea
        compile_base = '/tmp/%s/compiled_code/%s' %(user,self.id)
        os.makedirs( '%s/include' %compile_base)
        os.makedirs( '%s/src' %compile_base)
        os.makedirs( '%s/obj' %compile_base)
        os.system( 'cp %s/../Plotting/compiled_code/Makefile %s '%( workarea, compile_base ))

        brdef_file_name = '%s/include/BranchDefs.h'  %( compile_base )
        header_file_name = '%s/include/BranchInit.h' %( compile_base )
        source_file_name = '%s/src/BranchInit.cxx'   %( compile_base )
        linkdef_file_name = '%s/include/LinkDef.h'   %( compile_base )

        runsrc_file_name = '%s/src/RunAnalysis.cxx' %compile_base
        runinc_file_name = '%s/include/RunAnalysis.h' %compile_base

        # collect a complete list of branches used in all draw commands
        # normally this'll be fine, but its possible to
        # run into the situation where some trees may not have certain
        # branches and it'll break
        all_sample_branches = set()
        all_sample_chains = []
        for sample in self.samples :
            if sample.chain is not None :
                all_sample_chains.append( sample.chain )

        all_sample_branches = core.get_branch_mapping_from_trees( all_sample_chains )

        draw_branches = []
        draw_names = []
        for draw_config in self.draw_commands:
            draw_strs = ''
            for var in draw_config.var :
                draw_strs += var
            for sel in draw_config.selection :
                draw_strs += sel
            for br in all_sample_branches :
                if draw_strs.count(br['name']) :
                    if br['name'] not in draw_names :
                        draw_branches.append( br )
                        draw_names.append( br['name'] )

        n_tot = 0
        for draw_command in self.draw_commands :
            draw_command.create_hist_configs( draw_branches )
            n_tot += len(draw_command.hist_configs)

        print 'Will create %d histograms!' %n_tot

        if n_tot == 0 :
            print 'No histograms were scheduled.  Aborting!'
            return

        output_loc = '/tmp/%s/drawn_histograms/%s' %(user,self.id)
        os.makedirs( output_loc )

        # create the source code file
        self.write_source_code( self.draw_commands, runsrc_file_name, draw_branches )

        # create the header code file
        self.write_header_code( self.draw_commands, runinc_file_name )

        # Write the c++ files having the branch definitions and 
        # SetBranchAddress calls
        core.write_header_files(brdef_file_name, linkdef_file_name, draw_branches )

        core.write_source_file(source_file_name, header_file_name, draw_branches )

        # compile
        exename= 'RunAnalysis.exe'
        os.system( 'cd %s ; rm %s ; make EXENAME=%s; cd - '%(compile_base, exename, exename) )

        all_samples = []
        for sample in self.samples :
            if sample.isActive :
                if sample.IsGroupedSample() :
                    for subsamp in self.GetLowestGroupedSamples(sample) :
                        all_samples.append(subsamp)
                else :
                    all_samples.append(sample)

        # remove output directory
        if os.path.isdir( output_loc ) :
            os.system( 'rm -rf %s' %output_loc )

        output_file = 'hist.root'

        configs = []
        for sample in all_samples :
            config_name = '%s/configs/config_%s.txt' %(compile_base, sample.name)
            entries = sample.chain.GetEntries()
            #file_evt_map = core.get_file_evt_map( [f.GetTitle() for f in sample.chain.GetListOfFiles()], nsplit=nsplit, nFilesPerJob=nFilesPerJob, totalEvents=None, treeName='ggNtuplizer/EventTree')

            #exe_path = '%s/RunAnalysis' %compile_base

            #options = {}
            #options['nproc']             = 6
            #options['outputDir']         = '%s/%s' %(output_loc, sample.name)
            #options['confFileName']      = 'config_%s.txt' %sample.name
            #options['treeName']          = sample.chain.GetName()
            #options['outputFile']        = output_file
            #options['storagePath']       = None
            #options['sample']            = sample.name
            #options['disableOutputTree'] = True
            #options['nPrint']            = 10000

            #commands_orig = core.generate_multiprocessing_commands( file_evt_map, [], exe_path, options )

            #configs += commands_orig 

            ##core.write_config([], config_name, sample.chain.GetName(), output_loc, '%s.root'%sample.name, file_evt_map, sample=sample.name, disableOutputTree=True )
            ##configs.append((entries, config_name))
            file_evt_map = core.get_file_evt_map( [f.GetTitle() for f in sample.chain.GetListOfFiles()], nsplit=nsplit, nFilesPerJob=nFilesPerJob, totalEvents=None, treeName=treeName)
            core.write_config([], config_name, sample.chain.GetName(), output_loc, '%s.root'%sample.name, file_evt_map, sample=sample.name, disableOutputTree=True )
            configs.append((entries, config_name))

        
        #configs.sort(reverse=True)

        run_cmds = ['%s/%s --conf_file %s' %( compile_base,exename, c[1] ) for c in configs ]
        print run_cmds

        nproc = 6
        if len(run_cmds) < nproc :
            nproc = len(run_cmds)

        
        p=multiprocessing.Pool(nproc)
        p.map(os.system, run_cmds)

        file_map = {}
        # ------------------------
        # For older style without multiprocessing
        #for file in os.listdir( '%s/Job_0000' %output_loc ) :
        #    file_map[file] = []

        #for top, dirs, files in os.walk( output_loc ) :
        #    for file in files :
        #        if file in file_map.keys() :
        #            file_map[file].append( top+'/'+file ) 

        #comb_dir = output_loc+'/COMB'
        #os.makedirs( comb_dir )
        #for base, files in file_map.iteritems() :
        #    if len(files) == 1 :
        #        cp_cmd = 'cp %s %s/' %(files[0], comb_dir )
        #        print cp_cmd
        #        os.system( cp_cmd )
        #    else :
        #        os.system( 'hadd %s/%s %s' %( comb_dir, base, ' '.join( files ) ) )
        # ------------------------
        for sample in all_samples :
            os.system( 'hadd %s/%s/%s %s/%s/Job*/%s' %( output_loc, sample.name, output_file, output_loc, sample.name, output_file ) )

        self.output_loc = output_loc

        # Now get the histograms and draw
        for draw_config in self.draw_commands:
            if draw_config.no_auto_draw : 
                continue
            if draw_config.compare_hists :
                self.CompareFromHistFiles( draw_config )
            else :
                self.DrawFromHistFiles( draw_config )

    def DrawFromHistFiles(self,  draw_config ) :

        self.clear_all()

        all_samples = []
        for sample in self.get_samples(isActive=True):
            if sample.IsGroupedSample() :
                for subsamp in self.GetLowestGroupedSamples(sample) :
                    all_samples.append( subsamp)
            else :
                all_samples.append( sample )

        for sample in all_samples :
            # ------------------------
            # For older style without multiprocessing
            #filename = '%s/%s.root' %( self.output_loc, sample.name )
            # ------------------------
            filename = '%s/%s/hist.root' %( self.output_loc, sample.name )
            for name  in draw_config.hist_configs.keys() :
                self.load_hist_from_file_cache( sample, name, filename )

        # handle grouped samples
        for sample in self.samples :
            if sample.IsGroupedSample() and sample.isActive :
                self.group_sample(sample, isModel=False)
                
        if isinstance( draw_config.histpars, tuple) and len(draw_config.histpars) == 4 :
            if isinstance( draw_config.histpars[3], list ) :
                self.variable_rebinning(binning=draw_config.histpars[3]) 
            else :
                self.variable_rebinning(threshold=draw_config.histpars[3]) 

        self.MakeStack( draw_config )

        self.DrawCanvas(self.curr_stack, draw_config, datahists=['Data'], sighists=self.get_signal_samples())

        if draw_config.stack_dump_params :
            self.DumpStack( draw_config.stack_dump_params['dirname'], draw_config.stack_dump_params['filename'] )
        if draw_config.stack_save_params :
            self.SaveStack( draw_config.stack_save_params['filename'], draw_config.stack_save_params['dirname'], draw_config.stack_save_params['canname'] )
                

    def CompareFromHistFiles(self, draw_config ) :

        self.clear_all()

        ##-------------------
        ## list of samples may have duplicates
        ## get_samples does not pay attention to
        ## duplicates, go one-by-one
        ##-------------------
        #in_samples = []
        #for rsamp in draw_config.samples :
        #    in_samples += self.get_samples(name=rsamp) 

        #-----------------------
        # To handle the case when the
        # same sample is requested multiple
        # times, create a new sample for each
        # hist config
        #-----------------------
        created_samples = []
        for name, conf in draw_config.hist_configs.iteritems() :

            newsamp = self.clone_sample( oldname=conf['sample'], newname=name, temporary=True )
            print 'Create %s' %name

            if newsamp.IsGroupedSample() :
                for subsamp in self.GetLowestGroupedSamples(newsamp) :
                    # ------------------------
                    # For older style without multiprocessing
                    #filename = '%s/%s.root' %( self.output_loc, subsamp.name )
                    # ------------------------
                    filename = '%s/%s/hist.root' %( self.output_loc, subsamp.name )
                    self.load_hist_from_file_cache( subsamp, name, filename, debug=True )
                    subsamp.hist.Draw()

                self.group_sample(newsamp, isModel=False)
                newsamp.hist.Draw()
                    
            else :
                filename = '%s/%s.root' %( self.output_loc, newsamp.name )
                self.load_hist_from_file_cache( newsamp , name, filename )

            created_samples.append(newsamp)

        if not created_samples :
            print 'No hists were created from samples %s' %(', '.join(reqsamples) )
            return created_samples

        if isinstance( draw_config.histpars, tuple) and len(draw_config.histpars) == 4 :
            if isinstance( draw_config.histpars[3], list ) :
                self.variable_rebinning(binning=draw_config.histpars[3], samples=created_samples) 
            else :
                self.variable_rebinning(threshold=draw_config.histpars[3], samples=created_samples) 

        if draw_config.doRatio() :
            self.create_top_canvas_for_ratio('same')
        else :
            self.create_standard_canvas('same')

        self.curr_canvases['same'].cd()

        self.DrawSameCanvas( self.curr_canvases['same'], created_samples, draw_config )

        if draw_config.doRatio() :
            #rname = created_samples[0].name + '_ratio'
            for samp, hc in zip(created_samples[1:], draw_config.hist_configs.values()[1:]) :

                color = hc['color']
                rcolor = color

                if len( created_samples ) == 2 :
                    rcolor = ROOT.kBlack

                rname = 'ratio%s' %samp.name
                rsamp = self.create_ratio_sample( rname, num_sample = created_samples[0], den_sample=samp, color=rcolor)

                rsamp.legend_entry = hc.get('legend_entry', None )


        # make the legend
        step = len(created_samples)
        self.curr_legend = self.create_standard_legend(step, draw_config=draw_config )

        self.create_same_legend( draw_config.get_legend_entries() , created_samples )

        self.DrawCanvas(self.curr_canvases['same'], draw_config)

        if draw_config.stack_dump_params :
            self.DumpStack( draw_config.stack_dump_params['dirname'], draw_config.stack_dump_params['filename'], doRatio=draw_config.doRatio() )
        if draw_config.stack_save_params :
            self.SaveStack( draw_config.stack_save_params['filename'], draw_config.stack_save_params['dirname'], draw_config.stack_save_params['canname'] )
            
    def load_samples( self, draw_config ) :

        if draw_config.samp_man_id is not None and draw_config.samp_man_id != self.id :
            print 'Provided DrawConfig was not made with this SampleManager.  Exiting!'
            return []

        created_samples = []

        if not draw_config.hist_configs :
            print draw_config.samples[0].name
            matched_samples = self.get_samples( name=draw_config.samples[0].name)
            print matched_samples
            if matched_samples :
                if matched_samples[0].hist is not None :
                    created_samples.append( matched_samples[0] )

        else :
            for name, conf in draw_config.hist_configs.iteritems() :

            
                newsamp = self.clone_sample( oldname=conf['sample'], newname=name, temporary=True )
                print 'Create %s' %name

                if newsamp.IsGroupedSample() :
                    for subsamp in self.GetLowestGroupedSamples(newsamp) :
                        # ------------------------
                        # For older style without multiprocessing
                        #filename = '%s/%s.root' %( self.output_loc, subsamp.name )
                        # ------------------------
                        filename = '%s/%s/hist.root' %( self.output_loc, subsamp.name )
                        self.load_hist_from_file_cache( subsamp, name, filename, debug=True )
                        subsamp.hist.Draw()

                    self.group_sample(newsamp, isModel=False)
                        
                else :
                    # ------------------------
                    # For older style without multiprocessing
                    # filename = '%s/%s.root' %( self.output_loc, conf['sample'].name )
                    # ------------------------
                    filename = '%s/%s/hist.root' %( self.output_loc, conf['sample'].name )
                    self.load_hist_from_file_cache( newsamp , name, filename )

                created_samples.append(newsamp)

        if not created_samples :
            print 'No hists were created'
            return created_samples

        if isinstance( draw_config.histpars, tuple) and len(draw_config.histpars) == 4 :
            if isinstance( draw_config.histpars[3], list ) :
                self.variable_rebinning(binning=draw_config.histpars[3], samples=created_samples) 
            else :
                self.variable_rebinning(threshold=draw_config.histpars[3], samples=created_samples) 

        return created_samples



    def load_hist_from_file_cache( self, sample, name, filename, debug=False ) :

        if debug :
            print 'Load hist %s from %s into sample %s' %( name, filename, sample.name )
        sample.hist = None
        if not hasattr(sample, 'file' ) :
            sample.file = ROOT.TFile.Open( filename, 'READ')
        sample.hist = sample.file.Get(name).Clone()
        sample.hist.SetDirectory(0)
        sample.hist.Sumw2()
        self.format_hist( sample )
            

    def write_source_code( self, draw_commands, file, branches ) :

        text = ''

        text += r'#include "include/RunAnalysis.h"' + '\n'
        text += r'#include <iostream>' + '\n'
        text += r'#include <iomanip>' + '\n'
        text += r'#include <fstream>' + '\n'
        text += r'#include <sstream>' + '\n'
        text += r'#include <boost/foreach.hpp>' + '\n'
        text += r'#include <boost/algorithm/string.hpp>' + '\n'
        text += r'#include <sys/types.h>' + '\n'
        text += r'#include <sys/stat.h>' + '\n'
        text += r'#include <math.h>' + '\n'
        text += r'#include <stdlib.h>' + '\n'
        text += r'#include "include/BranchDefs.h"' + '\n'
        text += r'#include "include/BranchInit.h"' + '\n'
        text += r'#include "Core/Util.h"' + '\n'
        text += r'#include "TFile.h"' + '\n'
        text += r'int main(int argc, char **argv)' + '\n'
        text += r'{' + '\n'
        text += r'    CmdOptions options = ParseOptions( argc, argv );' + '\n'
        text += r'    AnaConfig ana_config = ParseConfig( options.config_file, options );' + '\n'
        text += r'    RunModule runmod;' + '\n'
        text += r'    ana_config.Run(runmod, options);' + '\n'
        text += r'    std::cout << "^_^ Finished ^_^" << std::endl;' + '\n'
        text += r'}' + '\n'
        text += r'void RunModule::initialize( TChain * chain, TTree * outtree, TFile *outfile,' + '\n'
        text += r'                            const CmdOptions & options, std::vector<ModuleConfig> &configs ) {' + '\n'
        text += r'    f = outfile; '+ '\n'
        text += r'    f->cd(); '+ '\n'
        text += r'    InitINTree(chain);' + '\n'

        for draw_config in draw_commands :
            for hist_str in draw_config.get_hist_declarations() :
               text += r' %s' %( hist_str ) + '\n\n';
               #text += r' hist_%s = new TH1F( "%s", "", %d, %f, %f );' %( draw_config.name, draw_config.name, draw_config.histpars[0], draw_config.histpars[1], draw_config.histpars[2] ) + '\n\n';
        text += r'}' + '\n'
        text += r'bool RunModule::execute( std::vector<ModuleConfig> & configs ) {' + '\n'
        for draw_config in self.draw_commands :
            for name in draw_config.get_names() :
                text += '    Draw%s(  ); \n' %name
            #text += '    Draw%s(  ); \n' %draw_config.name
        text += r'    return false;' + '\n'
        text += r'}' + '\n\n'

        text += r'void RunModule::finalize(  ) {' + '\n'
        for draw_config in self.draw_commands :
            for name in draw_config.get_names() :
                text += '    hist_%s->Write(); \n' %name
        text += r'}' + '\n\n'

        for draw_config in self.draw_commands :
            for name, config in draw_config.hist_configs.iteritems() :

                text += 'void RunModule::Draw%s( ) const { \n' %name
                first_replace = True
                for samp, rselection in draw_config.replace_selection_for_sample.iteritems() :
                    if first_replace :
                        text += '    if( curr_sample == %s ) { \n ' %samp
                        first_replace=False
                    else :
                        text += '    else if( curr_sample == %s ) { \n ' %samp

                    text += '        weight = %s; \n ' %rselection
                    text += '        if( weight != 0 ) { \n ' 
                    text += '        hist_%s->Fill(%s, weight); \n ' %(name, config['cppvar'])
                    text += '        } \n ' 
                    text += '    } \n ' 


                # just check if the replacement was done
                if first_replace : # no replacement
                    text += '    // Original selection : %s \n ' %config['selection']
                    text += '    float weight = %s; \n ' %config['cppselection']
                    text += '        if( weight != 0 ) { \n ' 
                    text += '        hist_%s->Fill(%s, weight); \n '  %(name, config['cppvar'])
                    text += '        } \n ' 
                else :
                    text += '    else { \n'
                    text += '        // Original selection : %s \n ' %config['selection']
                    text += '        float weight = %s; \n ' %config['cppselection']
                    text += '        if( weight != 0 ) { \n ' 
                    text += '        hist_%s->Fill(%s, weight); \n '  %(name, config['cppvar'])
                    text += '        } \n ' 
                    text += '    }\n '
                text += '}\n'


        ofile = open( file, 'w' )
        ofile.write(text)
        ofile.close()


    def write_header_code( self, draw_commands, file ) :

        text = ''

        text += '#ifndef RUNANALYSIS_H' + '\n'
        text += '#define RUNANALYSIS_H' + '\n'
        text += '#include "../../../Analysis/TreeFilter/Core/Core/AnalysisBase.h"' + '\n'
        text += '#include <string>' + '\n'
        text += '#include <vector>' + '\n'
        text += '#include "TTree.h"' + '\n'
        text += '#include "TChain.h"' + '\n'
        text += '#include "TH1F.h"' + '\n'
        text += '#include "TH2F.h"' + '\n'
        text += '#include "TH3F.h"' + '\n'
        text += '#include "TLorentzVector.h"' + '\n'
        text += 'class RunModule : public virtual RunModuleBase {' + '\n'
        text += '    public :' + '\n'
        text += '        RunModule() {}' + '\n'
        text += '        void initialize( TChain * chain, TTree *outtree, TFile *outfile, const CmdOptions & options, std::vector<ModuleConfig> & configs) ;' + '\n'
        text += '        bool execute( std::vector<ModuleConfig> & config ) ;' + '\n'
        text += '        void finalize( ) ;' + '\n'

        for draw_config in draw_commands :
            for name in draw_config.get_names() :
                text += '        void Draw%s ( ) const;' %name + '\n'

        for draw_config in draw_commands :
            for name in draw_config.get_names() :
                text += '        %s * hist_%s; '%(draw_config.get_hist_type(), name) + '\n'

        text += '            TFile * f;\n '

        text += '};' + '\n'
        text += 'namespace OUT {' + '\n'
        text += '};' + '\n'
        text += '#endif' + '\n'

        ofile = open( file, 'w' )
        ofile.write(text)
        ofile.close()


    #---------------------------------------
    def create_queued_hists( self, sample ) :


        for draw_config in self.draw_commands :
            sample.loop_hists.append( draw_config.init_hist(sample.name) )
            

        nentries = sample.chain.GetEntries()

        sample.chain.SetBranchStatus('*', 1)

        for draw_config in self.draw_commands :
            draw_config.compile_selection_string(sample, treename='sample.chain')
            draw_config.compile_var_string(sample, treename='sample.chain')

        for entry in sample.chain :

            for draw_config, hist in zip(self.draw_commands, sample.loop_hists ) :

                try :
                    weight = eval(draw_config.get_compiled_selection_string() )
                except :
                    print 'Failed to evaluate draw command.  Please check command and fix'
                    print draw_config.get_eval_selection_string(sample, treename='sample.chain')
                    raise
                if weight != 0 :
                    try :
                        hist.Fill( eval(draw_config.get_compiled_var( ) ), weight )
                    except :
                        print 'Failed to eval var.  Please check and fix'
                        print draw_config.get_var_val(sample, treename='sample.chain')
                        raise

        return len(sample.loop_hists)

    #---------------------------------------
    #def wait_on_draws(self ) :

    #    while self.curent_draws :
    #        to_rm = []
    #        for dr in self.curent_draws :

    #        self.curent_draws = [!(x.ready) for self.current

    #---------------------------------------
    def ListBranches(self, key=None ) :
        """ List all available branches.  If key is provided only show those that match the key """
    
        # grab list from 0th sample.  This may not work in some cases
        for br in self.samples[0].chain.GetListOfBranches() :
            if key is None :
                print br.GetName()
            else :
                if br.GetName().count(key) :
                    print br.GetName()

    #---------------------------------------
    def list_hists( self ) :
        all_hists = []
        for samp in self.get_samples() :
            for ofile in samp.ofiles :
                if ofile is not None :
                    for objkey in ofile.GetListOfKeys() :
                        obj = ofile.Get( objkey.GetName() )
                        class_name = obj.ClassName()
                        if class_name.count( 'TH1' ) :
                            all_hists.append(obj.GetName() )

        uniq_hists = list( set( all_hists ) )
        uniq_hists.sort()

        for h in uniq_hists :
            print h

    #---------------------------------------
    def SaveStack( self, filename, outputDir=None, canname=None, write_command=False, command_file='commands.txt'  ) :
        """ Save current plot to filename.  Must supply --outputDir  """
        
        if outputDir is None :
            print 'No output directory provided.  Will not save.'
        else :
            
            # write the command to a file if requested
            if write_command :
                if not os.path.isdir( outputDir ) :
                    os.makedirs( outputDir )
                cmdfile = open( outputDir +'/' +command_file, 'a' )
                cmdfile.write( '%s : %s \n' %( filename, self.transient_data.get( 'command', 'NO COMMAND STORED') ) )
                cmdfile.close()

            if self.collect_commands :
                self.add_save_stack( filename, outputDir, canname )
                return
    
            if not os.path.isdir( outputDir ) :
                print 'Creating directory %s' %outputDir
                os.makedirs(outputDir)

            filenamesplit = filename.split('.')
            if len( filenamesplit ) > 1 :
                filenamestrip = '.'.join( filenamesplit[:-1] )
            else :
                filenamestrip = filenamesplit[0]

            histnameeps = outputDir + '/' + filenamestrip+'.eps'
            #histnameeps = outputDir + '/' + filenamestrip+'.C'
            if not (filename.count( '.pdf' ) or filename.count('.png') ):
                histnamepdf = outputDir + '/' + filenamestrip+'.pdf'
            else :
                histnamepdf = outputDir + '/' +filename
    
            if len( self.curr_canvases ) == 0 :
                print 'No canvases to save'
            elif len( self.curr_canvases ) == 1  :
                self.curr_canvases.values()[0].SaveAs(histnamepdf)
                self.curr_canvases.values()[0].SaveAs(histnameeps)
            else :
                if canname is not None :
                    if canname not in self.curr_canvases :
                        print 'provided can name does not exist'
                    else :
                        self.curr_canvases[canname].SaveAs(histnamepdf)
                        self.curr_canvases[canname].SaveAs(histnameeps)

                else :
    
                    print 'Multiple canvases available.  Select which to save'
                    keys = self.curr_canvases.keys() 
                    for idx, key in enumerate(keys) :
                        print '%s (%d)' %(key, idx)
                    selidx = int(raw_input('enter number 0 - %d' %( len(keys)-1 )))
                    selkey = keys[selidx]
                    self.curr_canvases[selkey].SaveAs(histnamepdf)
                    self.curr_canvases[selkey].SaveAs(histnameeps)

    
    #---------------------------------------
    def DumpStack( self, outputDir=None, txtname=None, doRatio=None, details=False ) :
    
        if self.collect_commands :
            self.add_dump_stack( txtname, outputDir )
            return

        if doRatio is None :
            if self.draw_commands :
                doRatio = self.draw_commands[-1].doRatio()

        # store the signal and stack entries
        stack_entries  = {}
        signal_entries = {}
        ratio_entries  = {}
        detail_entries = {}
    
        # get samples with the MC stack, data, and signal samples
        samp_list = self.get_samples(name=self.get_stack_order()) + self.get_samples(isData=True) + self.get_samples(isSignal=True)
    
        # get the integrals
        for samp in samp_list :
            if samp.hist == None :
                continue
            err = ROOT.Double()
            integral = samp.hist.IntegralAndError( 1, samp.hist.GetNbinsX(), err )
            if samp.isSignal : 
                signal_entries[samp.name] = ufloat(integral, err)
            else :
                stack_entries[samp.name] = ufloat(integral, err )
        
        #collect the list to be printed 
        order = list(self.get_stack_order())
        if 'Data' in stack_entries :
            order.insert(0, 'Data')
        
        # get the sum over the full stack
        bkg_sum = ufloat(0.0, 0.0)
        for name, vals in stack_entries.iteritems() :
            if name != 'Data' :
                bkg_sum += vals
    
        sig_sum = ufloat(0.0, 0.0)
        for name, vals in signal_entries.iteritems() :
            sig_sum += vals
    
        latex_lines = []
        latex_lines.append( r'\begin{tabular}{| l | c |} ' )
        latex_lines.append( r'Sample & Events \\ \hline ' )
    
        lines = []
        for nm in order :
            if nm in stack_entries :
                lines.append('%s : \t %s' %( nm, stack_entries[nm] ))
                latex_lines.append( '%s & %s ' %( nm, stack_entries[nm] ) + r'\\')
    
        for sig in signal_entries :
            lines.append('%s : \t %s' %( sig, signal_entries[sig] ))
            latex_lines.append( '%s & %s ' %( sig, signal_entries[sig] )  + r'\\')
    
        lines.append('Stack Sum : \t %s' %(bkg_sum))
        latex_lines.append('Stack Sum & %s ' %(bkg_sum) + r'\\')
    
        for sig in signal_entries :
            den = umath.sqrt(signal_entries[sig] + bkg_sum )
            if den != 0 :
                lines.append('S/sqrt(S+B) (S=%s,B=All Bkg) : %s' %( sig, (signal_entries[sig]/den )) )
                latex_lines.append('S/sqrt(S+B) (S=%s) & %s ' %( sig, (signal_entries[sig]/den ) ) + r'\\')
            else :
                lines.append('S/sqrt(S+B) (S=%s,B=All Bkg) : nan' %( sig ) )
                latex_lines.append('S/sqrt(S+B) (S=%s) & nan ' %( sig) + r'\\')
    
        for sig in signal_entries :
            for st in stack_entries :
                den = umath.sqrt(signal_entries[sig] + stack_entries[st] )
                if den.n != 0 :
                    lines.append('S/sqrt(S+B) (S=%s,B=%s) : %s' %( sig, st,  (signal_entries[sig]/ den)  ))
                    #latex_lines.append('S/sqrt(S+B) (S=%s,B=%s) & %.2f' %( sig, st,  signal_entries[sig][0]/ den ) + r'\\')
                else :
                    lines.append('S/sqrt(S+B) (S=%s,B=%s) : NAN +- NAN' %( sig, st  ))
                    #latex_lines.append('S/sqrt(S+B) (S=%s,B=%s) & NAN' %( sig, st  ) + r'\\')
    
        if doRatio is not None and doRatio :
            rsamps = self.get_samples( isRatio=True )
            if rsamps :
                for rsamp in rsamps :
                    legend_entry = None
                    if hasattr( rsamp, 'legend_entry' ) :
                        legend_entry = rsamp.legend_entry
                    ratio_entries[rsamp.name] = {'legend_entry' : legend_entry, 'bins' : [] }
                    for _bin in range( 0, rsamp.hist.GetNbinsX() ) :
                        bin = _bin + 1
                        lines.append('%s, bin %d : %.3f += %.4f ' %( rsamp.name, bin, rsamp.hist.GetBinContent(bin), rsamp.hist.GetBinError(bin) ) )
                        ratio_entries[rsamp.name]['bins'].append(
                                                           {'bin' : bin, 'val' : rsamp.hist.GetBinContent(bin), 'err' : rsamp.hist.GetBinError(bin), 
                                                           'min' : rsamp.hist.GetXaxis().GetBinLowEdge(bin), 
                                                           'max' : rsamp.hist.GetXaxis().GetBinUpEdge(bin) }  )
        
        if details :
            detail_entries['detail'] = {}

            active_samps = self.get_samples(isActive=True)

            for samp in active_samps :
                if samp.hist is None :
                    continue
                detail_entries['detail'][samp.name] = {'bins' : {} }
                for bin in range(1, samp.hist.GetNbinsX()+1) :
                    val = samp.hist.GetBinContent(bin)
                    err = samp.hist.GetBinError(bin)

                    min = samp.hist.GetXaxis().GetBinLowEdge(bin)
                    max = samp.hist.GetXaxis().GetBinUpEdge(bin)

                    detail_entries['detail'][samp.name]['bins'][str(bin)] = {}
                    detail_entries['detail'][samp.name]['bins'][str(bin)]['val'] = ufloat( val, err )
                    detail_entries['detail'][samp.name]['bins'][str(bin)]['min'] = min
                    detail_entries['detail'][samp.name]['bins'][str(bin)]['max'] = max

                    lines.append('%s, bin %d (%s-%s) : %.3f += %.4f ' %( samp.name, bin, str(min), str(max), val, err ) )

        for line in lines :
            print line
    
        latex_lines.append( r'\end{tabular}' )

        if txtname is not None and outputDir is not None  :

            if txtname.count('.txt') == 0 :
                latexname = txtname + '.tex'
                picname = txtname + '.pickle'
                txtname += '.txt'
            else :
                latexname = txtname.rstrip('txt') + 'tex'
                picname = txtname.rstrip('txt') + 'pickle'

            if not os.path.isdir(outputDir ) :
                os.makedirs( outputDir )

            ##txtfile = open( outputDir + '/' + txtname, 'w')
            #txtfile = open( txtname, 'w')
            #for line in lines :
            #    txtfile.write( line + '\n' )
            #txtfile.close()
            #os.system( 'mv %s %s' %( txtname, outputDir ) )

            ##latexfile = open(outputDir + '/' + latexname, 'w')
            #latexfile = open(latexname, 'w')
            #for line in latex_lines  :
            #    latexfile.write( line + '\n' )
            #latexfile.close()
            #os.system( 'mv %s %s' %( latexname, outputDir ) )

            # write a pickle file
            stack_entries.update(signal_entries)
            stack_entries.update(ratio_entries)
            stack_entries.update(detail_entries)
            stack_entries['All Bkg'] = bkg_sum
            stack_entries['Total Expected'] = bkg_sum+sig_sum

            picfile = open( outputDir + '/' + picname, 'w' )
            pickle.dump( stack_entries, picfile )
            picfile.close()

        return
    
    #---------------------------------------
    def DumpRoc( self, outputDir=None, txtname=None, inDirs='' ) :
    
        output = []
        for title, entries in self.transient_data.iteritems() :
            output.append( title )
            print output[-1]
            for entry in entries :
                output.append('Cutval=%(CutVal)f, nSig=%(nSig)f, nBkg=%(nBkg)f, sigEff=%(sigEff)f, bkgEff=%(bkgEff)f, S/sqrt(S+B)=%(SoverRootSplusB)f ' %entry )
                print output[-1]
    
        if txtname is not None and outputDir is not None  :
    
            outdir = outputDir + '/' + inDirs
    
            if not os.path.isdir( outdir ) :
                print 'Making directory : ', outdir
                os.makedirs( outdir )
    
            if txtname.count('.txt') == 0 :
                txtname += '.txt'
    
            txtfile = open( outdir + '/' + txtname, 'w' )
            for out in output :
                txtfile.write( out + '\n' )
            txtfile.close()
    
    
    def ReloadSamples(self ) :

        #for samp in self.samples :
        #    if samp.chain is not None :
        #        samp.chain.Delete()

        self.samples = []

        self.ReadSamples(self.samples_conf )

    #--------------------------------
    def AddSample(self, name, path=None, filekey=None, isData=False, scale=None, isSignal=False, sigLineStyle=7, sigLineWidth=2, drawRatio=False, plotColor=ROOT.kBlack, lineColor=None, isActive=True, required=False, displayErrBand=False, useXSFile=False, XSName=None, legend_name=None) :
        """ Create an entry for this sample """

        if self.added_sample_group :
            print 'WARNING -- Please add all samples before adding a sample group, othwerise you may see funny behavior'

        # get all root files under this sample
        input_files = []

        # accept path as string (-> single entry list) or list
        if not isinstance(path, list) :
            path = [path]

        # collect files from each path 
        #
        # keep a list of paths that have
        # been traveresed and print a warning
        # if we're reading inconsistent paths
        subpaths_used = []
        paths_used    = []

        input_files = self.collect_input_files( self.base_path, path, paths_used, subpaths_used, filekey=filekey) 

        if input_files :
            #
            # Print a warning if we might be getting the wrong files
            #
            n_unique_pathlenghts = len( set( [ len( upath.split('/') ) for upath in paths_used ] ) )
            if n_unique_pathlenghts > 1 :
                print 'Found ntuples in parent directories, these may be duplicates'

            thisscale = 1.0
            # multiply by command line MC weight only for MC
            #if self.mcweight is not None and not isData :
            if self.mcweight is not None  :
                thisscale *= self.mcweight

            # multply by scale provided to this function
            if scale is not None :
                thisscale *= scale

            thisxs = 0
            totevt = -1
            if useXSFile :
                xsname = name
                if XSName is not None :
                    xsname = XSName
                if xsname in self.weightMap  :
                    thisscale *= self.weightMap[xsname]['scale']
                    thisxs = self.weightMap[xsname]['cross_section']
                    totevt = self.weightMap[xsname]['n_evt']
                    print 'Update scale for %s' %name

            thisSample = Sample(name, isActive=isActive, isData=isData, isSignal=isSignal, sigLineStyle=sigLineStyle, sigLineWidth=sigLineWidth, displayErrBand=displayErrBand, color=plotColor, drawRatio=drawRatio, scale=thisscale, cross_section=thisxs, total_events=totevt, legendName=legend_name)
            thisSample.AddFiles( input_files, self.treeName, self.readHists )

            self.samples.append(thisSample)

            # keep the order that this sample was added
            self.stack_order.append(name)

        if not input_files and required :
            print '***********************************************' 
            print 'Sample, %s does not exist and is required by this module' %name
            print self.base_path
            print '***********************************************' 
            sys.exit(-1)

	print_prefix = "AddSample: Reading %s (%s) " %(name, path )
        print_prefix = print_prefix.ljust(60)
        if not input_files :
            print print_prefix + " [ \033[1;31mFailed\033[0m  ]"
        else :
            if len(set(path) & set(subpaths_used)) != len(path) :
                print print_prefix + " [ \033[1;33mPartial\033[0m ]"
                print path
                print subpaths_used
            else :
                print print_prefix + " [ \033[1;32mSuccess\033[0m ]" 

    #--------------------------------
    def collect_input_files( self, base_path, path_list, paths_used, subpaths_used, filekey=None ) :

        input_files = []

        # if the necessary inputs are not provided return an empty list
        if base_path is None or path_list is None or not path_list :
            return input_files

        for subpath in path_list :
            fullpath = base_path + '/' + subpath
            # if files have been provided, read them directly
            if os.path.isfile(fullpath) :
                input_files.append(fullpath)
            else : #otherwise search directories for the needed files
                if base_path.count( 'root://eoscms' ) :
                    for top, dirs, files, sizes in eos_utilities.walk_eos( base_path + '/' + subpath ) :
                        for file in files :
                            if filekey is not None :
                                if file.count(filekey) == 0 : 
                                    continue
                            elif file.count(self.fileName) == 0 :
                                continue
                            paths_used.append(top)
                            subpaths_used.append(subpath)
                            input_files.append(top+'/'+file)
                            
                else : # local directories
                    for top, dirs, files in os.walk( base_path +'/' + subpath , followlinks=True) :
                        for file in files :
                            if filekey is not None :
                                if file.count(filekey) == 0 : 
                                    continue
                            elif file.count(self.fileName) == 0 :
                                continue
                            paths_used.append(top)
                            subpaths_used.append(subpath)
                            input_files.append(top+'/'+file)

        return input_files

    #--------------------------------
    def AddSampleGroup(self, name, input_samples=[], isData=False, scale=None, isSignal=False, sigLineStyle=7,drawRatio=False, plotColor=ROOT.kBlack, lineColor=None, legend_name=None, isActive=True, displayErrBand=False, averageSamples=False ) :
        """Make a new sample from any number of samples that have already been added via AddSample

           For example if a process is made of a number of individual samples that each have their
           own weight, first add those samples using AddSample with their own scale ( be sure
           to also give isActive=True or the individual samples will be drawn).  Then call
           Group Samples with the list of input samples and the new name.  
        """

        self.added_sample_group=True
        #check if input samples actually exist
        available_samples = []
        for samp in input_samples :
            if samp not in self.get_sample_names() :
                print 'WARNING - Child sample, %s, does not exist!' %samp
            else :
                available_samples.append(samp)

        # if no input samples exist, exit
        if not available_samples :
            return

        # keep the order that this sample was added
        self.stack_order.append(name)

        thisscale = 1.0
        # multply by scale provided to this function (MCweight is applied to input samples)
        if scale is not None :
            thisscale *= scale
        
        print 'Grouping %s' %name
        thisSample = Sample(name, isActive=isActive, isData=isData, isSignal=isSignal, sigLineStyle=sigLineStyle, displayErrBand=displayErrBand, color=plotColor, drawRatio=drawRatio, scale=thisscale, legendName=legend_name, averageSamples=averageSamples)

        for samp in available_samples :

            is_a_grouped_sample = ( name in self.get_grouped_sample_names() )

            if is_a_grouped_sample :
                group_samples = self.get_samples(name=name)[0].groupedSamples
                thisSample.AddGroupSamples( group_samples )
            else :
                thisSample.AddGroupSamples( samp )

        self.samples.append(thisSample)
                
    def AddModelSampleGroup(self, name, input_samples=[], isData=False, scale=None, isSignal=False, drawRatio=False, plotColor=ROOT.kBlack, lineColor=None, legend_name=None, isActive=True) :
        """Make a new sample from any number of samples that have already been added via AddSample

           For example if a process is made of a number of individual samples that each have their
           own weight, first add those samples using AddSample with their own scale ( be sure
           to also give isActive=True or the individual samples will be drawn).  Then call
           Group Samples with the list of input samples and the new name.  
        """

        thisscale = 1.0
        # multply by scale provided to this function (MCweight is applied to input samples)
        if scale is not None :
            thisscale *= scale
        
        print 'Grouping %s' %name
        thisSample = Sample(name, isActive=isActive, isData=isData, isSignal=isSignal, color=plotColor, drawRatio=drawRatio, scale=thisscale, legendName=legend_name)

        for samp in input_samples :
            is_a_grouped_sample = ( name in self.get_grouped_sample_names() )

            if is_a_grouped_sample :
                group_samples = self.get_samples(name=name).groupedSamples
                thisSample.AddGroupSamples( group_samples )
            else :
                thisSample.AddGroupSamples( samp )

        self.modelSamples.append(thisSample)


    def AddModelSample(self, name, legend_name=None, path=None, scale=1.0 , filekey=None, plotColor=ROOT.kBlack) :
        input_files = []

        if not isinstance(path, list) :
            path = [path]

        # keep a list of paths that have
        # been traveresed and print a warning
        # if we're reading inconsistent paths
        subpaths_used = []
        paths_used    = []

        input_files = self.collect_input_files( self.base_path_model, path, paths_used, subpaths_used, filekey=filekey) 
        
        if input_files :
            #
            # Print a warning if we might be getting the wrong files
            #
            n_unique_pathlenghts = len( set( [ len( upath.split('/') ) for upath in paths_used ] ) )
            if n_unique_pathlenghts > 1 :
                print 'Found ntuples in parent directories, these may be duplicates'

            thisscale = 1.0
            if scale is not None :
                thisscale *= scale

            thisSample = Sample(name, color=plotColor, scale=thisscale, legend_name=legend_name)
            thisSample.AddFiles( self.treeNameModel, input_files )
            self.modelSamples.append(thisSample)

	print_prefix = "AddModelSample: Reading %s (%s) " %(name, path )
        print_prefix = print_prefix.ljust(60)
        if not input_files :
            print print_prefix + " [ \033[1;31mFailed\033[0m  ]"
        else :
            if len(set(path) & set(subpaths_used)) != len(path) :
                print print_prefix + " [ \033[1;33mPartial\033[0m ]"
            else :
                print print_prefix + " [ \033[1;32mSuccess\033[0m ]" 


    def ReadSamples(self, conf, expected=[], failOnMissing=False ) :

        self.samples_conf = conf

        self.expected_samples = expected
        self.fail_on_missing = failOnMissing

        ImportedModule=None

        ispath = ( conf.count('/') > 0 )
        module_path = None
        if ispath :
            module_path = conf
        else :
            #get path of this script
            script_path = os.path.realpath(__file__)
            module_path = os.path.dirname(script_path) + '/' + conf

        try :
            ImportedModule = imp.load_source(conf.split('.')[0], module_path)
        except IOError :
            print 'Could not import module %s' %module_path
        
        if hasattr(ImportedModule, 'config_samples') :
            print '-------------------------------------'
            print 'BEGIN READING SAMPLES'
            print '-------------------------------------'
            ImportedModule.config_samples(self)
        else :
            print 'ERROR - samplesConf does not implement a function called config_samples '
            sys.exit(-1)
        
        if hasattr(ImportedModule, 'print_examples') :
            ImportedModule.print_examples()
        else :
            print 'WARNING - samplesConf does not implement a function called print_examples '


    def DrawHist(self, histpath, rebin=None, varRebinThresh=None, doratio=False, subtract_bkg=False, ylabel=None, xlabel=None, rlabel=None, logy=False, ymin=None, ymax=None, rmin=None, rmax=None, xmin=None, xmax=None, normalize=False, ticks_x=None, ticks_y=None, label_config={}, legend_config={}) :

        self.clear_all()
        for sample in self.samples :
            self.get_hist( sample, histpath )

        if varRebinThresh is not None :
            self.variable_rebinning(varRebinThresh) 

        if rebin is not None :
            for samp in self.get_samples( isActive=True ) :
                if samp.hist is not None :
                    samp.hist.Rebin(rebin)

        draw_config = DrawConfig( histpath, None, None, hist_config={'doratio' : doratio, 'xlabel' : xlabel, 'ylabel' : ylabel, 'ymin' : ymin, 'ymax' : ymax, 'logy' : logy, 'rmin' : rmin, 'rmax' : rmax, 'normalize' : normalize, 'ticks_x' : ticks_x, 'ticks_y' : ticks_y} , label_config=label_config, legend_config=legend_config)

        # make a stack if there are samples to be 
        # stacked.  If no stacked samples exist
        # then just histograms are being drawn
        if self.get_stack_order()  :
            self.MakeStack(draw_config )

        if ylabel is None :
            binwidth = self.get_samples(isActive=True)[0].hist.GetBinWidth(1)
            ylabel = 'Events / %.1f GeV' %binwidth
        if rlabel is None :
            rlabel = 'Data / MC'

        datahists = ['Data']
        if subtract_bkg :
            new_data = self.clone_sample( 'Data', 'DataSub', temporary=True )
            data_hist = self.get_samples( name='Data' )[0].hist
            bkg_hist = self.get_samples( name='Background' )[0].hist
            new_bkg = self.clone_sample( 'Background', 'BackgroundSub', temporary=True )

            new_data.hist = data_hist.Clone( 'DataSub' )
            new_bkg.hist = bkg_hist.Clone( 'BackgroundSub' )
            for bin in range( 1, data_hist.GetNbinsX() + 1 ) :
                data_val = data_hist.GetBinContent( bin )
                data_err = data_hist.GetBinError( bin )
                
                bkg_val = bkg_hist.GetBinContent( bin )
                bkg_err = bkg_hist.GetBinError( bin )

                new_bkg.hist.SetBinContent( bin, 0. )

                new_val = data_val - bkg_val
                #new_err = math.sqrt( data_err*data_err + bkg_err*bkg_err )
                new_err = data_err

                new_data.hist.SetBinContent( bin, new_val )
                new_data.hist.SetBinError( bin, new_err )

            datahists = ['DataSub']

        errSamp = None
        err_samp_list = [ s for s in self.get_samples( isActive=True ) if s.displayErrBand ]
        if subtract_bkg :
            err_samp_list = self.get_samples( name='DataSub' )
            #err_samp_list = []

        for samp in err_samp_list :
            bkg_name = '__AllStack__'
            if subtract_bkg :
                bkg_name = 'BackgroundSub'

            all_stack_hist = self.get_samples( name=bkg_name )[0].hist
            if errSamp is None :
                errSamp = self.clone_sample( samp.name, 'err_band', temporary=True )
                errSamp.hist = all_stack_hist.Clone('err_band')
                for bin in range(1, errSamp.hist.GetNbinsX() + 1 ) :
                    errSamp.hist.SetBinContent(bin, 0)
                    errSamp.hist.SetBinError(bin, 0)
            for bin in range(1, all_stack_hist.GetNbinsX()+1 ) :
                # set the bin content to the stack sum
                # unless the background is subtracted
                # at which point it should be at zero
                if not subtract_bkg :
                    errSamp.hist.SetBinContent( bin, all_stack_hist.GetBinContent( bin ) )
                curr_err = errSamp.hist.GetBinError( bin )
                this_err = samp.hist.GetBinError(bin)
                errSamp.hist.SetBinError( bin, math.sqrt( curr_err*curr_err + this_err*this_err )  )

        
        # for plotting a TGraphErrors is needed
        if errSamp is not None :
            errSamp.graph = ROOT.TGraphErrors( errSamp.hist )

        # calculate a chi-2 with respect to the errors
        if errSamp is not None :
            chi2sum = 0
            dataSamp = self.get_samples( name=datahists[0] )[0]
            nbins = errSamp.hist.GetNbinsX()
            for binnum in range( 1,  nbins +1 ) :

                exp_bkg = errSamp.hist.GetBinContent( binnum )
                exp_err = errSamp.hist.GetBinError( binnum )
                data_val = dataSamp.hist.GetBinContent( binnum )
                data_err = dataSamp.hist.GetBinError( binnum )

                tot_err = math.sqrt( data_err*data_err + exp_err*exp_err )
        
                if not tot_err == 0 :
                    chi2sum += math.pow( (exp_bkg-data_val) / tot_err, 2 )

            pzero = ROOT.TMath.Prob( chi2sum, nbins )
            print pzero

            kstest = dataSamp.hist.KolmogorovTest( errSamp.hist  )
            print kstest

            label_style = label_config.get('labelStyle', '')
            if not label_style.count( 'nostats' ) :
                draw_config.label_config['statsLabel'] = '#chi^{2} = %.1f, p0 = %.3f ' %(chi2sum , pzero )

        topcan = self.curr_stack

        self.DrawCanvas(topcan, draw_config, datahists=datahists, sighists=self.get_signal_samples(), errhists=['err_band']  )

        if xmin is not None and xmax is not None :
            for prim in self.curr_canvases['top'].GetListOfPrimitives() :
                if isinstance( prim, ROOT.TH1F ) :
                    prim.GetXaxis().SetRangeUser( xmin, xmax )

        if not self.get_stack_order()  :
            # if histograms are being drawn then
            # most of the work happens in DrawCanvas
            # but we still need to create a legend

            # first find the active samples
            active_samps = self.get_samples( isActive=True )
            self.curr_legend = self.create_standard_legend( len( active_samps ), draw_config=draw_config)

            legendOrder = legend_config.get( 'legendOrder', [] )

            if legendOrder :
                active_samps = []

                for entry in legendOrder :
                    active_samps.append( self.get_samples( name = entry )[0] )
                    
                
            self.create_same_legend(legend_entries=legendOrder, created_samples=active_samps)

            self.curr_legend.Draw()

    def Draw(self, varexp, selection, histpars, hist_config={}, label_config={}, legend_config={}, treeHist=None, treeSelection=None, generate_data_from_sample=None, replace_selection_for_sample={} , useModel=False ) :

        
        if self.collect_commands :
            self.add_draw_config( varexp, selection, histpars, hist_config=hist_config, label_config=label_config, legend_config=legend_config, replace_selection_for_sample=replace_selection_for_sample  )
            return

        config = DrawConfig( varexp, selection, histpars, hist_config=hist_config, label_config=label_config, legend_config=legend_config, replace_selection_for_sample=replace_selection_for_sample  )

        
        command = 'samples.Draw(\'%s\',  \'%s\', %s )' %( varexp, selection, str( histpars ) )
        self.transient_data['command'] = command 

        self.draw_and_configure( config, generate_data_from_sample=generate_data_from_sample, useModel=useModel, treeHist=treeHist, treeSelection=treeSelection )

    def draw_and_configure( self, draw_config, generate_data_from_sample=None, useModel=False, treeHist=None, treeSelection=None ) :

        self.clear_all()

        #move to somewhere else
        #self.apply_lenged_conf( legendConfig )

        res = self.draw_active_samples( draw_config )
        if not res :
            return

        if generate_data_from_sample is not None :
            samp_list = self.get_samples( name=generate_data_from_sample )
            if samp_list : 
                rand = ROOT.TRandom3()
                rand.SetSeed( int( time.time() ) )
                nbins = samp_list[0].hist.GetNbinsX()
                for bin in range( 1, nbins+1 ) :
                    newval = rand.Poisson( samp_list[0].hist.GetBinContent(bin) )
                    samp_list[0].hist.SetBinContent( bin, newval ) 
                    if newval > 0 :
                        samp_list[0].hist.SetBinError( bin, math.sqrt(newval) ) 


        if useModel :
            for sample in self.modelSamples :
                self.create_hist_new( draw_config, sample, isModel=True )

            # Model is created, replace the sample in self.samples with the
            # sample having the same name in self.modelSamples
            for samp in self.modelSamples :
                if samp.name in self.get_sample_names() :
                    self.get_samples(name=name).hist = samp.hist
                    self.get_samples(name=name).legendName = samp.legendName

        if isinstance( draw_config.histpars, tuple) and len(draw_config.histpars) == 4 :
            if isinstance( draw_config.histpars[3], list ) :
                self.variable_rebinning(binning=draw_config.histpars[3]) 
            else :
                self.variable_rebinning(threshold=draw_config.histpars[3]) 

        self.MakeStack(draw_config, useModel, treeHist, treeSelection )

        self.DrawCanvas(self.curr_stack, draw_config, datahists=['Data'], sighists=self.get_signal_samples())

    def Draw3DProjections(self, varexp, selection, histpars=None, x_by_y_bin_vals={}, doratio=False, ylabel=None, xlabel=None, rlabel=None, logy=False, ymin=None, ymax=None, ymax_scale=None, rmin=None, rmax=None, showBackgroundTotal=False, backgroundLabel='AllBkg', removeFromBkg=[], addToBkg=[], useModel=False, treeHist=None, treeSelection=None, labelStyle=None, extra_label=None, extra_label_loc=None, generate_data_from_sample=None, replace_selection_for_sample={}, legendConfig=None  ) :

        command = 'samples.Draw(\'%s\',  \'%s\', %s )' %( varexp, selection, str( histpars ) )

        if not x_by_y_bin_vals :
            print 'Must give a dictionary that maps y bins to x bins'

        self.clear_all()

        self.transient_data['command'] = command 

        self.apply_lenged_conf( legendConfig )

        res = self.draw_active_samples( draw_config )

        if not res :
            return

        if generate_data_from_sample is not None :
            samp_list = self.get_samples( name=generate_data_from_sample )
            if samp_list : 
                rand = ROOT.TRandom3()
                rand.SetSeed( int( time.time() ) )
                nbins = samp_list[0].hist.GetNbinsX()
                for bin in range( 1, nbins+1 ) :
                    newval = rand.Poisson( samp_list[0].hist.GetBinContent(bin) )
                    samp_list[0].hist.SetBinContent( bin, newval ) 
                    if newval > 0 :
                        samp_list[0].hist.SetBinError( bin, math.sqrt(newval) ) 


        if useModel :
            for sample in self.modelSamples :
                self.create_hist( sample, treeHist, treeSelection, histpars, isModel=True )

            # Model is created, replace the sample in self.samples with the
            # sample having the same name in self.modelSamples
            for samp in self.modelSamples :
                if samp.name in self.get_sample_names() :
                    self.get_samples(name=name).hist = samp.hist
                    self.get_samples(name=name).legendName = samp.legendName

        for sample in self.get_samples() :
            if sample.hist is not None :
                sample.main_hist = sample.hist.Clone( )


            else :
                sample.main_hist = None

        for (xmin, xmax), yvals_raw in x_by_y_bin_vals.iteritems() :

            if not isinstance( yvals_raw[0], tuple) :
                yvals = []
                for idx, yval_raw_min in enumerate( yvals_raw[:-1] ) :
                    yval_raw_max = yvals_raw[idx+1]
                    yvals.append( (yval_raw_min, yval_raw_max) )
            else :
                yvals = yvals_raw
            
            for ymin, ymax in yvals :
                for sample in self.get_samples() :
                    if hasattr(sample, 'main_hist' ) and sample.main_hist is not None :
                
                        if xmin is None :
                            xbin_min = 1
                        else :
                            xbin_min = sample.main_hist.GetXaxis().FindBin( xmin )

                        if xmax is None :
                            xbin_max = sample.main_hist.GetNbinsX() 
                        else :
                            xbin_max = sample.main_hist.GetXaxis().FindBin( xmax )

                        if ymin is None :
                            ybin_min = 1
                        else :
                            ybin_min = sample.main_hist.GetYaxis().FindBin( ymin )

                        if ymax is None :
                            ybin_max = sample.main_hist.GetNbinsY() 
                        else :
                            ybin_max = sample.main_hist.GetYaxis().FindBin( ymax )

                        print 'xmin = %d, xmax = %s, ymin = %f, ymax = %f, xbinmin = %d, xbinmax = %d, ybinmin = %d, yminmax = %d' %( xmin, xmax, ymin, ymax, xbin_min, xbin_max, ybin_min, ybin_max) 
                        sample.hist = sample.main_hist.ProjectionZ( str( uuid.uuid4()), xbin_min, xbin_max, ybin_min, ybin_max )
                    else :
                        sample.hist = None
        
                
                self.MakeStack(varexp, doratio, showBackgroundTotal, backgroundLabel, removeFromBkg, addToBkg, useModel, treeHist, treeSelection )

                if ylabel is None :
                    bin_width = ( histpars[2] - histpars[1] )/histpars[0]
                    bin_width_f = ( histpars[2] - histpars[1] )/float(histpars[0])
                    if math.fabs(bin_width_f - bin_width) != 0 :
                        ylabel = 'Events / %.1f GeV' %bin_width_f
                    else :
                        ylabel = 'Events / %d GeV' %bin_width
                if rlabel is None :
                    rlabel = 'Data / MC'
                    
                self.DrawCanvas(self.curr_stack, ylabel=ylabel, xlabel=xlabel, rlabel=rlabel, logy=logy, ymin=ymin, ymax=ymax, ymax_scale=ymax_scale, rmin=rmin, rmax=rmax, datahists=['Data'], sighists=self.get_signal_samples(), doratio=doratio, labelStyle=labelStyle, extra_label=extra_label, extra_label_loc=extra_label_loc )

                yield (xmin, xmax, ymin, ymax)


    def DrawSamples(self, varexp, selection, samples, histpars=None, normalize=False, doratio=False, useTreeModel=False, treeHist=None, treeSelection=None ) :
        if not isinstance( samples, list ) :
            samples = [samples]

        self.MakeSameCanvas(samples, varexp, selection, histpars, doratio)
        self.DrawSameCanvas(normalize, doratio)

    def MakeStack(self, draw_config, useModel=False, treeHist=None, treeSelection=None ) :

        # Get info for summed sample
        bkg_name = '__AllStack__'
        # get all stacked histograms and add them
        stack_samples = self.get_samples( name=self.get_stack_order() )

        if stack_samples :
            sum_hist = stack_samples[0].hist.Clone(bkg_name)
            for samp in stack_samples[1:] :
                sum_hist.Add(samp.hist)

            self.create_sample( bkg_name, isActive=False, hist=sum_hist, temporary=True )

        doratio = draw_config.doRatio()

        if doratio :
            # when stacking, the ratio is made with respect to the data.  Find the sample that
            # is labeled as data.  Throw an error if one data sample is not found
            data_samples = self.get_samples(isData=True)
            if not data_samples :
                print 'MakeStack - ERROR : No data samples found!'

            self.create_ratio_sample( 'ratio', num_sample=data_samples[0], den_sample='__AllStack__' )

            # make ratio histograms for signal samples
            signal_samples = self.get_samples( isSignal=True, drawRatio=True )
            for sample in signal_samples :
                ratio_samp = self.create_ratio_sample( sample.name + '_ratio', num_sample=data_samples[0], den_sample=sample )
                ratio_samp.hist.SetLineColor( sample.color )
                ratio_samp.isSignal = True

        #make the stack and fill
        self.curr_stack = (ROOT.THStack(str(uuid.uuid4()), ''))

        # reverse so that the stack is in the correct order
        orderd_samples = []
        for sampname in self.get_stack_order():
            samplist = self.get_samples(name=sampname, isActive=True )
            if samplist :
                orderd_samples.append(samplist[0])
        for samp in reversed(orderd_samples) :              
            samp.hist.SetFillColor( samp.color )
            samp.hist.SetLineColor( ROOT.kBlack )
            samp.hist.SetLineWidth( 1 )

            self.curr_stack.Add(samp.hist, 'HIST')

        # additional formatting
        data_samp = self.get_samples(name='Data')

        # make the legend
        # In placing the legend move the bottom down 0.05 for each entry
        # calculate the step usa
        drawn_samples = []
        for sname in self.get_stack_order():
            drawn_samples+=self.get_samples( name=sname, failed_draw=False, isActive=True )
        step = len(drawn_samples)

        self.curr_legend = self.create_standard_legend( step, draw_config=draw_config)

        legendTextSize = draw_config.legend_config.get('legendTextSize', 0.04 )

        # format the entries
        tmp_legend_entries = []
        legend_entries = []

        if data_samp and data_samp[0].isActive :
            tmp_legend_entries.append(  (data_samp[0].hist, data_samp[0].legendName, 'PE') )

        for samp in drawn_samples :
            tmp_legend_entries.append( ( samp.hist, samp.legendName,  'F') )

        print '********************NOT FILLING SIGNAL ENTRY IN LEGEND**********************'
        #for samp in self.get_signal_samples() :
        #    if samp.isActive :
        #        tmp_legend_entries.append( ( samp.hist, samp.legendName, 'L') )

        if self.legendLoc=='Double' :
            legend_entries = [None]*len(tmp_legend_entries)
            if len(legend_entries)%2 == 0 :
                n_first_col = len(legend_entries)/2 
            else :
                n_first_col = (len(legend_entries)/2) + 1
            n_2nd_col = 0
            for idx in range(0, len(tmp_legend_entries) ) :
                if idx%2 == 0 : 
                    if idx < n_first_col :
                        newidx = idx*2
                        legend_entries[newidx] = tmp_legend_entries[idx]
                    else :
                        n_2nd_col+=1
                        newidx = n_2nd_col*2-1
                        legend_entries[newidx] = tmp_legend_entries[idx]
                else : 
                    if idx < n_first_col :
                        newidx = idx*2;
                        legend_entries[newidx] = tmp_legend_entries[idx]
                    else :
                        n_2nd_col+=1
                        newidx = n_2nd_col*2-1
                        legend_entries[newidx] = tmp_legend_entries[idx]

        else :
            legend_entries = tmp_legend_entries
        for le in legend_entries :
            entry = self.curr_legend.AddEntry(le[0], le[1], le[2])
            if legendTextSize is not None :
                entry.SetTextSize(legendTextSize)


    #----------------------------------------------------
    def MakeSameCanvas(self, draw_config, useStoredBinning=False, preserve_hists=False, useModel=False, treeHist=None, treeSelection=None) :

        if not preserve_hists :
            self.clear_all()

        created_samples = []
        for hist_name, hist_config in draw_config.hist_configs.iteritems() :

            samp = hist_config['sample']
            selection = hist_config['selection']

            # In this case the same sample may be drawn multiple
            # times.  to avoid any conflicts, add new samples
            # and draw into those

            newname = hist_name

            newsamp = self.clone_sample( oldname=samp, newname=newname, temporary=True )

            if useModel :
                self.create_hist( newsamp, treeHist, treeSelection, draw_config.histpars, isModel=True)
            else :
                self.create_hist( newsamp, hist_config['var'], hist_config['selection'], draw_config.histpars)

            created_samples.append( newsamp )

        if not created_samples :
            print 'No hists were created'
            return created_samples

        if isinstance( draw_config.histpars, tuple) and len(draw_config.histpars) == 4 :
            if isinstance( draw_config.histpars[3], list ) :
                self.variable_rebinning(binning=draw_config.histpars[3], samples=created_samples, useStoredBinning=useStoredBinning) 
            else :
                self.variable_rebinning(threshold=draw_config.histpars[3], samples=created_samples, useStoredBinning=useStoredBinning) 

        if draw_config.doRatio() :
            self.create_top_canvas_for_ratio('same')
        else :
            self.create_standard_canvas('same')

        self.curr_canvases['same'].cd()

        self.DrawSameCanvas( self.curr_canvases['same'], created_samples, draw_config )

        if draw_config.doRatio() :
            #rname = created_samples[0].name + '_ratio'
            for hist_name in draw_config.hist_configs.keys()[1:] :
                hist_config = draw_config.hist_configs[hist_name]

                samp = hist_config['sample']
                color = hist_config['color']

                rcolor = color
                if len( draw_config.hist_configs ) == 2 :
                    rcolor = ROOT.kBlack

                rname = 'ratio%s' %samp
                if rname in self.get_sample_names() :
                    for i in range(0, 100 ) :
                        rname = 'ratio%s_%d' %(samp, i)
                        if rname not in self.get_sample_names() :
                            break
                rsamp = self.create_ratio_sample( rname, num_sample = draw_config.hist_configs.keys()[0], den_sample=hist_name, color=rcolor)
                rsamp.legend_entry = hist_config.get('legend_entry', None )

        return created_samples


    def get_hist( self, sample, histpath ) :
        sampname = sample.name
        print 'Getting hist for %s' %sampname

        # check that this histogram hasn't been drawn
        if sample.hist is not None :
            print 'Histogram already extracted for %s' %sampname
            return

        # Draw the histogram.  Use histpars as the bin limits if given
        if sample.IsGroupedSample() :
            for subsampname in sample.groupedSamples :
                subsamp = self.get_samples( name=subsampname )[0]
                print 'Extract grouped hist %s' %subsampname
                if subsampname in self.get_sample_names() :
                    self.get_hist( subsamp, histpath )

            
            self.group_sample( sample )
            return

        else :

            thishist = sample.ofiles[0].Get(histpath)
            if thishist == None :
                print 'Could not get hist!'
                sample.isActive=False
            else :
                for ofile in sample.ofiles[1:] :
                    thishist.Add( ofile.Get(histpath) )

                if sample.hist is not None :
                    sample.hist.Delete()

                sample.hist = thishist
                if sample.hist is not None :
                    self.format_hist( sample )
    

    def create_hist( self, sample, varexp, selection, histpars, isModel=False ) :

        if isinstance( sample, str) :
            slist = self.get_samples( name=sample )
            if not slist :
                print 'Could not retrieve sample, %s' %sample
                return False
            if len(slist) > 1 :
                print 'Located multiple samples with name %s' %sample
                return False
            sample = slist[0]

        sampname = sample.name
    
        if not self.quiet : print 'Creating hist for %s' %sampname
        if not self.quiet : print '%s : %s ' %( varexp, selection )
        if not self.quiet : print histpars

        ## check that this histogram hasn't been drawn
        #if sample.hist is not None :
        #    print 'Histogram already drawn for %s' %sampname
        #    return

        histname = str(uuid.uuid4())

        full_selection = selection

        # enable branches for all variables matched in the varexp and selection
        sample.enable_parsed_branches( varexp+selection ) 

        sample.hist = None
        if type( histpars ) is tuple :
            if varexp.count(':') == 1 : 
                if len(histpars) == 2 and type( histpars[0] ) is list and type(histpars[1]) is list :
                    sample.hist = ROOT.TH2F( histname, '', len(histpars[0])-1, array('f', histpars[0]), len(histpars[1])-1, array('f', histpars[1]) )
                else :
                    if len(histpars) != 6 :
                        print 'varable expression, %s,  requests a 2-d histogram, please provide 6 hist parameters, nbinsx, xmin, xmax, nbinsy, ymin, ymax' %varexp
                        return
                    sample.hist = ROOT.TH2F( histname, '', histpars[0], histpars[1], histpars[2], histpars[3], histpars[4], histpars[5])
            elif varexp.count(':') == 2 and not varexp.count('::') : # make a 3-d histogram
                if len(histpars) != 9 :
                    print 'varable expression, %s,  requests a 3-d histogram, please provide 9 hist parameters, nbinsx, xmin, xmax, nbinsy, ymin, ymax, nbinsz, zmin, zmax' %varexp
                    return
                sample.hist= ROOT.TH3F( histname, '',histpars[0], histpars[1], histpars[2], histpars[3], histpars[4], histpars[5], histpars[6], histpars[7], histpars[8] )
            else : # 1-d histogram

                sample.hist= ROOT.TH1F( histname, '', histpars[0], histpars[1], histpars[2])

        elif type( histpars ) is list :
            sample.hist = ROOT.TH1F( histname, '', len(histpars)-1, array('f', histpars))
        else :
            print 'No histogram parameters were passed'

        if sample.hist is not None :
            sample.hist.SetTitle( sampname )
            sample.hist.Sumw2()
            ROOT.SetOwnership(sample.hist, False )

        # Draw the histogram.  Use histpars as the bin limits if given
        if sample.IsGroupedSample() :
            for subsampname in sample.groupedSamples :
                subsamp = self.get_samples( name=subsampname )[0]
                
                if not self.quiet : print 'Draw grouped hist %s' %subsampname

                if isModel and subsampname in [s.name for s in self.get_model_samples()] :
                    self.create_hist( subsamp, varexp, selection, histpars, isModel=isModel )
                elif subsampname in self.get_sample_names() :
                    self.create_hist( subsamp, varexp, selection, histpars, isModel=isModel )

            sample.failed_draw=False
            for subsampname in sample.groupedSamples :
                subsamp = self.get_samples( name=subsampname )[0]
                if subsamp.failed_draw :
                    sample.failed_draw=True


            self.group_sample( sample, isModel=isModel )

            return True

        else :
            if sample.chain is not None :
                #self.draw_hist( sample, varexp, histname, full_selection, draw_opt='goff' )
                res = sample.chain.Draw(varexp + ' >> ' + histname, full_selection, 'goff' )
                if res < 0 :
                    sample.failed_draw=True
                else :
                    sample.failed_draw=False
            else :
                sample.failed_draw=True

            if sample.hist is not None :
                self.format_hist( sample )

            return True

        # Group draw parallelization
        # wait for draws to finish
        #self.wait_on_draws()

    #def draw_hist( self, sample, varexp, histname, selection, draw_opt='' ) :

    def create_hist_new( self, draw_config, sample=None, isModel=False ) :

        if sample is None :
            sample = draw_config.samples[0]

        if isinstance( sample, str) :
            slist = self.get_samples( name=sample )
            if not slist :
                print 'Could not retrieve sample, %s' %sample
            if len(slist) > 1 :
                print 'Located multiple samples with name %s' %sample
            sample = slist[0]

        if not self.quiet : print 'Creating hist for %s' %sample.name

        # enable branches for all variables matched in the varexp and selection

        sample.hist = draw_config.init_hist(sample.name)
        selection = draw_config.get_selection_string( sample.name )
        varexp    = draw_config.var[0]

        if not self.quiet : print selection

        sample.enable_parsed_branches( varexp+selection ) 

        # Draw the histogram.  Use histpars as the bin limits if given
        if sample.IsGroupedSample() :
            for subsampname in sample.groupedSamples :
                subsamp = self.get_samples( name=subsampname )[0]
                
                if not self.quiet : print 'Draw grouped hist %s' %subsampname

                if isModel and subsampname in [s.name for s in self.get_model_samples()] :
                    self.create_hist_new( draw_config, subsamp, isModel=isModel )
                elif subsampname in self.get_sample_names() :
                    self.create_hist_new( draw_config, subsamp, isModel=isModel )

            sample.failed_draw=False
            for subsampname in sample.groupedSamples :
                subsamp = self.get_samples( name=subsampname )[0]
                if subsamp.failed_draw :
                    sample.failed_draw=True


            self.group_sample( sample, isModel=isModel )

            return

        else :
            if sample.chain is not None :
                res = sample.chain.Draw(varexp + ' >> ' + sample.hist.GetName(), selection , 'goff' )
                if res < 0 :
                    sample.failed_draw=True
                else :
                    sample.failed_draw=False
            else :
                sample.failed_draw=True

            if sample.hist is not None :
                self.format_hist( sample )

        # Group draw parallelization
        # wait for draws to finish
        #self.wait_on_draws()


    def format_hist( self, sample ) :
        self.AddOverflow( sample.hist )
        sample.SetHist()

    def AddOverflow(self,  hist ) :

        # account for overflow and underflow
        # go 3-2-1 because TH3 inherits from TH2 etc..
        if isinstance( hist, ROOT.TH3 ) :
            nbinsx = hist.GetNbinsX()
            nbinsy = hist.GetNbinsY()
            nbinsz = hist.GetNbinsZ()

            #print 'Content before'
            #sum = 0
            #for xbin in range( 0, nbinsx+2 ) :
            #    for ybin in range( 0, nbinsy+2 ) :
            #        for zbin in range( 0, nbinsz+2 ) :
            #            sum += hist.GetBinContent( xbin, ybin, zbin ) 
            #            print 'Content (%d, %d, %d) = %d' %( xbin, ybin, zbin, hist.GetBinContent( xbin, ybin, zbin ) )
            #print 'SUM = ', sum
            for xbin in range( 1, nbinsx+1 ) :
                for ybin in range( 1, nbinsy+1 ) :
                    for zbin in range( 1, nbinsz+1 ) :
                        if xbin == 1 or xbin == nbinsx or ybin == 1 or ybin == nbinsy or zbin == 1 or zbin == nbinsz  :
                            # only combine on edge bins

                            bin = [ xbin, ybin, zbin ]

                            all_mod_bins = {}
                            # get singly modified bins (TH1, TH2, TH3)
                            self.get_single_mod_bins( bin, all_mod_bins )
                            # get doubly modified bins (TH2, TH3)
                            self.get_multi_mod_bins( bin, all_mod_bins )
                            # get triply modified bins (TH3)
                            self.get_multi_mod_bins( bin, all_mod_bins )

                            #get the full list of bins that are overflow
                            all_bins = []
                            #print all_mod_bins
                            for vals in all_mod_bins.values() :
                                for val in vals :
                                    if ( (val[0]==nbinsx+1 or val[0]==0 or val[0]==xbin) and 
                                         (val[1]==nbinsy+1 or val[1]==0 or val[1]==ybin) and 
                                         (val[2]==nbinsz+1 or val[2]==0 or val[2]==zbin) ) :
                                        all_bins.append(val)

                            # unique the list
                            unique_bins = [list(x) for x in set( tuple(x)  for x in all_bins) ]

                            for ovbin in unique_bins :
                                #print 'Combine bin %s with bin %s' %( str(bin), str(ovbin) ) 
                                self.combine_overflow_bin( hist, bin, ovbin )
            #print 'Content after'
            #sum = 0
            #for xbin in range( 0, nbinsx+2 ) :
            #    for ybin in range( 0, nbinsy+2 ) :
            #        for zbin in range( 0, nbinsz+2 ) :
            #            sum += hist.GetBinContent( xbin, ybin, zbin )
            #            print 'Content (%d, %d, %d) = %d' %( xbin, ybin, zbin, hist.GetBinContent( xbin, ybin, zbin ) )
            #print 'SUM = ', sum

        elif isinstance( hist, ROOT.TH2 ) :

            nbinsx = hist.GetNbinsX()
            nbinsy = hist.GetNbinsY()
            #print 'Content before'
            #sum = 0
            #for xbin in range( 0, nbinsx+2 ) :
            #    for ybin in range( 0, nbinsy+2 ) :
            #        sum += hist.GetBinContent( xbin, ybin )
            #        print 'Content (%d, %d ) = %d' %( xbin, ybin,  hist.GetBinContent( xbin, ybin ) )
            #print 'SUM = ', sum

            for xbin in range( 1, nbinsx+1 ) :
                for ybin in range( 1, nbinsy+1 ) :
                    if xbin == 1 or xbin == nbinsx  or ybin == 1 or ybin == nbinsy :
                        bin = [ xbin, ybin ]

                        all_mod_bins = {}
                        # get singly modified bins (TH1, TH2, TH3)
                        self.get_single_mod_bins( bin, all_mod_bins )
                        # get doubly modified bins (TH2, TH3)
                        self.get_multi_mod_bins( bin, all_mod_bins )

                        #print all_mod_bins

                        #get the full list of bins that are overflow
                        all_bins = []
                        for vals in all_mod_bins.values() :
                            for val in vals :
                                if ( (val[0]==nbinsx+1 or val[0]==0 or val[0]==xbin) and 
                                    (val[1]==nbinsy+1 or val[1]==0 or val[1]==ybin) ) :
                                    all_bins.append(val)

                        # unique the list
                        unique_bins = [list(x) for x in set( tuple(x)  for x in all_bins) ]

                        for ovbin in unique_bins :
                            self.combine_overflow_bin( hist, bin, ovbin )

            #print 'Content after'
            #sum = 0
            #for xbin in range( 0, nbinsx+2 ) :
            #    for ybin in range( 0, nbinsy+2 ) :
            #        sum += hist.GetBinContent( xbin, ybin )
            #        print 'Content (%d, %d) = %d' %( xbin, ybin,  hist.GetBinContent( xbin, ybin ) )
            #print 'SUM = ', sum
        elif isinstance( hist, ROOT.TH1 ) :

            nbinsx = hist.GetNbinsX()
            #overflow
            self.combine_overflow_bin( hist, nbinsx, nbinsx+1 )
            #underflow
            self.combine_overflow_bin( hist, 1, 0 )


    def get_single_mod_bins( self, bin, bin_collection ) :
        # for each entry in the bin, 
        # make a modified bin
        # for that value up and that value down
        for binentry in range(0, len( bin ) ):
            binup = list(bin)
            binup [binentry]+=1
            bindn = list(bin)
            bindn [binentry]-=1
            bin_collection.setdefault( (binentry,), []).append(binup)
            bin_collection.setdefault( (binentry,), []).append(bindn)

    def get_multi_mod_bins( self, bin, bin_collection ) :
        if not bin_collection :
            print 'Must have already defined bins to run multi'
            return

        # loop over the already defined bins
        for modbinval in bin_collection.keys() :
            modbins = bin_collection[modbinval]
            # grab entries that have fewer than the number of axes modified
            if len(modbinval) < len(bin) :
                for modbin in modbins :
                    # get an entry that was not modified before
                    # and modify it
                    for binentry in range(0, len( bin ) ) :
                        if binentry in modbinval :
                            continue
                        binup = list( modbin )
                        binup[binentry]+=1
                        bindn = list( modbin )
                        bindn[binentry]-=1
                        bin_collection.setdefault( modbinval+(binentry,), []).append(binup)
                        bin_collection.setdefault( modbinval+(binentry,), []).append(bindn)


    def combine_overflow_bin( self, hist, bin, overbin ) :

        if not isinstance( overbin, tuple ) :
            if isinstance( overbin, list ) :
                overbin = tuple( overbin )
            else :
                overbin = (overbin,)
        if not isinstance( bin, tuple ) :
            if isinstance( bin, list) :
                bin = tuple( bin )
            else :
                bin = (bin,)

        over_val = hist.GetBinContent( *overbin )

        if over_val != 0 :

            orig_val = hist.GetBinContent( *bin )
            orig_err = hist.GetBinError( *bin )
            over_err = hist.GetBinError( *overbin )

            new_val = orig_val + over_val

            if orig_val == 0 :
                new_err = over_err
            else :
                new_err = math.sqrt( orig_err*orig_err + over_err*over_err )

            hist.SetBinContent( *( bin + (new_val,) ) )
            hist.SetBinError( *( bin + (new_err,) ) )
            hist.SetBinContent( *( overbin + (0,) ) )
            hist.SetBinError( *( overbin + (0,) ) )


    def get_active_samples( self, histpath ) :

        for sample in self.samples :
            if sample.isActive :
                self.get_hist( sample, histpath )

    def draw_active_samples( self, draw_config ) :

        failed_samples = []
        success_samples = []
        for sample in self.samples :
            if sample.isActive :
                self.create_hist_new( draw_config, sample )
                if sample.failed_draw :
                    failed_samples.append( sample.name )
                else :
                    success_samples.append( sample.name )

        for samp in failed_samples :
            print 'Failed to draw sample %s' %samp

        if not success_samples :
            return False

        return True
            

    def variable_rebinning(self, threshold=None, binning=None, samples=[], useStoredBinning=False) :

        if not samples:
            samples = self.get_samples(name=self.get_stack_order())

        # variable r
        if binning is not None :
            for samp in self.get_samples() :
                if samp.hist is not None :
                    samp.SetHist(self.do_variable_rebinning(samp, binning))
            return

        elif threshold is not None :

            if useStoredBinning :
                binning = self.binning
            else :
                all_stack_hists = []
                for samp in samples :
                    all_stack_hists.append(samp.hist)

                binning = self.make_variable_binning( all_stack_hists, threshold)

                # store binning for future use
                self.binning = binning

            for samp in self.get_samples() :
                if samp.hist is not None :
                    samp.SetHist(self.do_variable_rebinning(samp, binning))
            return

        else :
            print 'variable_rebinning : Must provide a rebinning threshold, or a binning scheme'
            return


    def group_sample(self, sample, isModel=False) :

        if not sample.IsGroupedSample() :
            print 'Trying to group a sample that is not a grouped sample'
            return
        
        subsamp_names = sample.groupedSamples
        if not self.quiet : print 'RUN GROUPING FOR %s' %sample.name
        if not self.quiet : print subsamp_names
        
        if isModel :
            model_subsamps = self.get_model_samples(subsamp_names)
            sample.hist = model_subsamps[0].hist.Clone()
            for msamp in model_subsamps[1:] :

                if sample.averageSamples : 
                    self.addHistsWeightedAvg( sample.hist, samp.hist )

                sample.hist.Add( msamp.hist )
                sample.hist.Draw()
            #sample.hist.Scale(sample.scale)
            #self.modelSamples.append(sample)
        else :
            subsamps = self.get_samples(name=subsamp_names)

            for samp in subsamps :
                if samp.IsGroupedSample() :
                    self.group_sample( samp, isModel=False )

            valid_samps = [s for s in subsamps if s.hist is not None]

            sample.hist = valid_samps[0].hist.Clone()
            for samp in valid_samps[1:] :

                if sample.averageSamples : 
                    self.addHistsWeightedAvg( sample.hist, samp.hist )

                sample.hist.Add( samp.hist )
                #sample.hist.Draw()
            #sample.hist.Scale(sample.scale)

        sample.InitHist()

    def get_list_from_tree(self, vars, selection, sample ) :

        output = []

        if not isinstance( vars, list ) :
            vars = [vars]

        if sample.IsGroupedSample() :
            for subsampname in sample.groupedSamples :
                subsamp = self.get_samples( name=subsampname )[0]
                
                if not self.quiet : print 'running on grouped hist %s' %subsampname

                if subsampname in self.get_sample_names() :
                    output += self.get_list_from_tree( vars, selection, subsamp )

            return output

        else :

            sample.copied_tree = sample.chain.CopyTree( selection )

            nentries = sample.copied_tree.GetEntries()

            for i in xrange( 0, nentries ) :
                sample.copied_tree.GetEntry( i )
                evt_entries = []
                for var in vars :

                    newvar = var
                    #if an index is requested, get the name and index separately
                    res = re.match( '(\w+)\[(\d)\]', var )
                    if res is not None :
                        newvar = res.group(1)
                        idx = int(res.group(2) )
                        vec = getattr( sample.copied_tree, newvar )
                        evt_entries.append( vec[idx] )
                    else :
                        evt_entries.append( getattr( sample.copied_tree, var ) )

                evt_entries.append( sample.scale )

                output.append( tuple( evt_entries ) )

        return output

    def set_canvas_default_formatting(self, topcan, doratio, logy=False, ylabel=None, ymin=None, ymax=None ) :

        for prim in topcan.GetListOfPrimitives() :
            if isinstance(prim, ROOT.TH1F) :
                if ymin is not None and ymax is not None :
                    prim.GetYaxis().SetRangeUser( ymin, ymax )
                prim.SetTitle('')
                offset = 1.15
                if logy :
                    offset = 1.1
                if doratio == True or doratio == 1 : # canvas sizes differ for ratio, so title, label sizes are different
                    prim.GetYaxis().SetTitleSize(0.06)
                    prim.GetYaxis().SetTitleOffset(offset)
                    prim.GetYaxis().SetLabelSize(0.06)
                    prim.GetXaxis().SetLabelSize(0.06)
                    prim.GetXaxis().SetTitleSize(0.06)
                elif doratio == 2 : 
                    prim.GetYaxis().SetTitleSize(0.06)
                    prim.GetYaxis().SetTitleOffset(offset)
                    prim.GetYaxis().SetLabelSize(0.06)
                    prim.GetXaxis().SetLabelSize(0.06)
                    prim.GetXaxis().SetTitleSize(0.06)
                else :
                    prim.GetYaxis().SetTitleSize(0.05)
                    prim.GetYaxis().SetTitleOffset( offset )
                    prim.GetYaxis().SetLabelSize(0.05)
                    prim.GetXaxis().SetLabelSize(0.05)
                    prim.GetXaxis().SetTitleSize(0.05)
                    prim.GetXaxis().SetTitleSize(0.05)

    def set_stack_default_formatting(self, topcan, doratio, logy=False ) :
        if topcan.GetHists() != None :
            if topcan.GetHists().GetSize() > 0 :
                offset = 1.15
                if logy :
                    offset = 1.1
                if doratio : # canvas sizes differ for ratio, so title, label sizes are different
                    topcan.GetHistogram().GetYaxis().SetTitleSize(0.06)
                    topcan.GetHistogram().GetYaxis().SetTitleOffset(offset)
                    topcan.GetHistogram().GetYaxis().SetLabelSize(0.06)
                    topcan.GetHistogram().GetXaxis().SetLabelSize(0.06)
                    topcan.GetHistogram().GetXaxis().SetTitleSize(0.06)
                else :
                    topcan.GetHistogram().GetYaxis().SetTitleSize(0.05)
                    topcan.GetHistogram().GetYaxis().SetTitleOffset(offset)
                    topcan.GetHistogram().GetYaxis().SetLabelSize(0.05)
                    topcan.GetHistogram().GetXaxis().SetLabelSize(0.05)
                    topcan.GetHistogram().GetXaxis().SetTitleSize(0.05)

    def set_ratio_default_formatting(self, canvas, ratiosamps, draw_config ) :

            canvas.cd()
                
            doratio = draw_config.doRatio()
            rlabel = draw_config.get_rlabel()
            rmin   = draw_config.get_rmin()
            rmax   = draw_config.get_rmax()



            for idx, ratiosamp in enumerate( ratiosamps ) :
                drawopt = 'same'
                if idx == 0 :
                    drawopt = ''
                if ratiosamp.isSignal :
                    drawopt += 'HIST'

                if doratio == True or doratio == 1 :
                    ratiosamp.hist.GetYaxis().SetTitleSize(0.10)
                    ratiosamp.hist.GetYaxis().SetTitleOffset(0.6)
                    ratiosamp.hist.GetYaxis().SetLabelSize(0.12)
                    ratiosamp.hist.GetXaxis().SetLabelSize(0.12)
                    ratiosamp.hist.GetXaxis().SetTitleSize(0.12)
                    ratiosamp.hist.GetXaxis().SetTitleOffset(1.0)
                elif doratio==2 :
                    ratiosamp.hist.GetYaxis().SetTitleSize(0.06)
                    ratiosamp.hist.GetYaxis().SetTitleOffset(0.8)
                    ratiosamp.hist.GetYaxis().SetLabelSize(0.06)
                    ratiosamp.hist.GetXaxis().SetLabelSize(0.06)
                    ratiosamp.hist.GetXaxis().SetTitleSize(0.06)
                    ratiosamp.hist.GetXaxis().SetTitleOffset(1.0)

                ratiosamp.hist.SetStats(0)
                ratiosamp.hist.SetMarkerStyle(20)
                ratiosamp.hist.SetMarkerSize(1.1)
                ratiosamp.hist.SetTitle('')
                ratiosamp.hist.GetYaxis().CenterTitle()
                ratiosamp.hist.GetYaxis().SetNdivisions(506, True)
                if rlabel is not None :
                    ratiosamp.hist.GetYaxis().SetTitle(rlabel)
                if rmin is not None and rmax is not None :
                    ratiosamp.hist.GetYaxis().SetRangeUser(rmin, rmax)

            #left_edge  = ratiosamps[0].hist.GetXaxis().GetXmin()
            #right_edge = ratiosamps[0].hist.GetXaxis().GetXmax()

            #canvas.cd()

            #oneline = ROOT.TLine(left_edge, 1, right_edge, 1)
            #oneline.SetLineStyle(3)
            #oneline.SetLineWidth(2)
            #oneline.SetLineColor(ROOT.kBlack)
            #oneline.Draw()
            #self.add_decoration(oneline)

    def calc_yaxis_limits(self, draw_config ) :

        ymin       = draw_config.get_ymin()
        ymax       = draw_config.get_ymax()
        ymax_scale = draw_config.get_ymax_scale()
        
        calcymax = 0
        calcymin = 0.5
        for samp in ( self.get_samples(isActive=True ) + self.get_samples( name='__AllStack__' ) )  :
            if samp.hist == None :
                continue
            max = samp.hist.GetMaximum()
            min = samp.hist.GetMaximum()

            if max > calcymax :
                calcymax = max
            if min < ymin :
                calcymin = min

        if ymax is None :
            ymax = calcymax 

        if ymin is None :
            ymin = calcymin 

        if ymax_scale is not None :
            ymax *= ymax_scale
        else :
            ymax *= 1.2

        return (ymin, ymax)

    def create_standard_canvas(self, name='base') :

        xsize = 650 
        ysize = 500
        self.curr_canvases[name] = ROOT.TCanvas(name, name, xsize, ysize)

        self.curr_canvases['top'] = self.curr_canvases[name]
        #self.curr_canvases[name].SetTopMargin(0.08)
        #self.curr_canvases[name].SetBottomMargin(0.13)
        #self.curr_canvases[name].SetLeftMargin(0.13)
        #self.curr_canvases[name].SetTitle('')
        self.curr_canvases[name].SetTopMargin(0.08)
        self.curr_canvases[name].SetBottomMargin(0.13)
        self.curr_canvases[name].SetLeftMargin(0.15)
        self.curr_canvases[name].SetRightMargin(0.05)
        self.curr_canvases[name].SetTitle('')

    def create_standard_ratio_canvas(self) :

        xsize = 620 
        ysize = 620
        self.curr_canvases['base'] = ROOT.TCanvas('basecan', 'basecan', xsize, ysize)

        self.curr_canvases['bottom'] = ROOT.TPad('bottompad', 'bottompad', 0.01, 0.01, 0.99, 0.34)
        self.curr_canvases['top'] = ROOT.TPad('toppad', 'toppad', 0.01, 0.35, 0.99, 0.99)
        self.curr_canvases['top'].SetTopMargin(0.08)
        self.curr_canvases['top'].SetBottomMargin(0.06)
        self.curr_canvases['top'].SetLeftMargin(0.15)
        self.curr_canvases['top'].SetRightMargin(0.05)
        #self.curr_canvases['bottom'].SetTopMargin(0.05)
        self.curr_canvases['bottom'].SetTopMargin(0.00)
        self.curr_canvases['bottom'].SetBottomMargin(0.3)
        self.curr_canvases['bottom'].SetLeftMargin(0.15)
        self.curr_canvases['bottom'].SetRightMargin(0.05)
        self.curr_canvases['base'].cd()
        self.curr_canvases['bottom'].Draw()
        self.curr_canvases['top'].Draw()

    def create_large_ratio_canvas(self) :

        xsize = 600 
        ysize = 850
        self.curr_canvases['base'] = ROOT.TCanvas('basecan', 'basecan', xsize, ysize)

        self.curr_canvases['bottom'] = ROOT.TPad('bottompad', 'bottompad', 0.01, 0.01, 0.99, 0.48)
        self.curr_canvases['top'] = ROOT.TPad('toppad', 'toppad', 0.01, 0.49, 0.99, 0.99)
        self.curr_canvases['top'].SetTopMargin(0.08)
        self.curr_canvases['top'].SetBottomMargin(0.02)
        self.curr_canvases['top'].SetLeftMargin(0.15)
        self.curr_canvases['top'].SetRightMargin(0.05)
        #self.curr_canvases['bottom'].SetTopMargin(0.05)
        self.curr_canvases['bottom'].SetTopMargin(0.03)  # non-zero to allow space for numbers on the axis that go above the plot
        self.curr_canvases['bottom'].SetBottomMargin(0.2)
        self.curr_canvases['bottom'].SetLeftMargin(0.15)
        self.curr_canvases['bottom'].SetRightMargin(0.05)
        self.curr_canvases['base'].cd()
        self.curr_canvases['bottom'].Draw()
        self.curr_canvases['top'].Draw()

    def create_top_canvas_for_ratio(self, name='can') :

        xsize = 650 
        ysize = 500
        self.curr_canvases[name] = ROOT.TCanvas(name, name, xsize, ysize)

        self.curr_canvases[name].SetTopMargin(0.08)
        self.curr_canvases[name].SetBottomMargin(0.08)
        self.curr_canvases[name].SetLeftMargin(0.15)
        self.curr_canvases[name].SetRightMargin(0.05)
        self.curr_canvases[name].SetTitle('')


    def DrawCanvas(self, topcan, draw_config, datahists=[], sighists=[], errhists=[] ) :

        doratio=draw_config.doRatio()
        if doratio == True or doratio == 1 :
            self.create_standard_ratio_canvas()
        elif doratio == 2 :
            self.create_large_ratio_canvas()
        else :
            self.create_standard_canvas() 

        (ymin, ymax) = self.calc_yaxis_limits( draw_config )
        
        self.curr_canvases['top'].cd()

        ylabel = draw_config.get_ylabel()
        ticksx = draw_config.get_tick_x_format()
        ticksy = draw_config.get_tick_y_format()

        # in what cases is topcan 
        # an already filled TCanvas?
        if isinstance(topcan, ROOT.TCanvas ) :
            for prim in topcan.GetListOfPrimitives() :
                if isinstance(prim, ROOT.TH1F) :
                    if ylabel is not None :
                        prim.GetYaxis().SetTitle(ylabel)
                    if ticksx is not None :
                        prim.GetXaxis().SetNdivisions( ticksx[0],ticksx[1],ticksx[2] )
                    if ticksy is not None :
                        prim.GetYaxis().SetNdivisions( ticksy[0],ticksy[1],ticksy[2] )



            topcan.DrawClonePad()

        elif isinstance(topcan, ROOT.THStack ) :
            if topcan.GetHists() == None :
                print 'Top stack has no hists! Exiting'
                return
            topcan.Draw()
            topcan.SetMinimum(ymin)
            topcan.SetMaximum(ymax)
            self.set_stack_default_formatting( topcan, doratio, logy=draw_config.get_logy())

        normalize = draw_config.get_normalize()

        # draw the data
        for dsamp in self.get_samples( name=datahists, isActive=True ):
            dsamp.hist.SetMarkerStyle(20)
            if normalize :
                dsamp.hist.DrawNormalized('PE same')
            else :
                dsamp.hist.Draw('PE same')


        # draw the signals
        legendTextSize = draw_config.legend_config.get('legendTextSize', 0.04 )
        print legendTextSize
        if sighists :
            sigsamps = self.get_samples(name=sighists)
            for samp in sighists : 
                print samp.isActive
                if samp.isActive :
                    print 'Draw Signal hist ', samp.name
                    #samp.hist.SetLineWidth(3)
                    samp.hist.SetLineStyle(samp.getLineStyle())
                    samp.hist.SetStats(0)
                    if normalize :
                        samp.hist.DrawNormalized('HIST same')
                    else :
                        samp.hist.Draw('HIST same')
                    entry = self.curr_legend.AddEntry( samp.hist, samp.legendName, 'L')
                    entry.SetTextSize(legendTextSize)

        if errhists :
            errsamps = self.get_samples(name=errhists )
            for samp in errsamps :
                samp.graph.SetFillStyle(3004)
                samp.graph.SetFillColor(ROOT.kBlack)
                samp.graph.Draw('2same')

                entry = self.curr_legend.AddEntry( samp.graph, 'Total uncertainty', 'F')
                entry.SetTextSize(legendTextSize)


        if doratio :
            self.curr_canvases['bottom'].cd()
            ratiosamps =  self.get_samples( isRatio=True )
            self.set_ratio_default_formatting( self.curr_canvases['bottom'], ratiosamps, draw_config )

            for idx, samp in enumerate(ratiosamps) :
                drawopt = 'same'
                if idx == 0 :
                    drawopt = ''
                if samp.isSignal :
                    drawopt += 'HIST'

                samp.hist.Draw(drawopt)

        self.curr_canvases['top'].cd()

        # if topcan is None then a new
        # canvas has been created above
        # and the formatting should be set
        if topcan is None :
            self.set_canvas_default_formatting(self.curr_canvases['top'], doratio, logy=draw_config.get_logy(), ylabel=ylabel, ymin=ymin, ymax=ymax )

        # draw the legend
        if self.curr_legend is not None :
            self.curr_legend.Draw()

        # draw the plot status label
        labels = draw_config.get_labels()
        for lab in labels :
            lab.Draw()
            self.curr_decorations.append( lab )

        if doratio :

            self.curr_canvases['bottom'].cd()

            left_edge  = ratiosamps[0].hist.GetXaxis().GetXmin()
            right_edge = ratiosamps[0].hist.GetXaxis().GetXmax()

            oneline = ROOT.TLine(left_edge, 1, right_edge, 1)
            oneline.SetLineStyle(3)
            oneline.SetLineWidth(2)
            oneline.SetLineColor(ROOT.kBlack)
            oneline.Draw()
            self.add_decoration(oneline)

        xlabel = draw_config.get_xlabel()
        if xlabel is not None :
            if doratio :
                ratiosamp = self.get_samples(isRatio=True)[0]
                ratiosamp.hist.GetXaxis().SetTitle(xlabel)
            else :
                if isinstance(topcan, ROOT.THStack ) :
                    topcan.GetHistogram().GetXaxis().SetTitle(xlabel)

                for samp in self.get_samples() :
                    if samp.hist :
                        samp.hist.GetXaxis().SetTitle(xlabel)

        if isinstance(topcan, ROOT.THStack ) :
            #print topcan.GetHistogram().GetNdivisions()
            #raw_input('cont')
            if ylabel is not None :
                topcan.GetHistogram().GetYaxis().SetTitle(ylabel)
            if ticksx is not None : 
                topcan.GetHistogram().GetXaxis().SetNdivisions(ticksx[0], ticksx[1], ticksx[2],True)
            if ticksy is not None : 
                topcan.GetHistogram().GetYaxis().SetNdivisions(ticksy[0], ticksy[1], ticksy[2], True)


            wh = self.curr_canvases['top'].GetWh()
            ww = self.curr_canvases['top'].GetWw()

            tl = 15.

            tick_frac_y = tl / ( ww )
            tick_frac_x = tl / ( wh )


            topcan.GetHistogram().GetXaxis().SetTickLength(tick_frac_x)
            topcan.GetHistogram().GetYaxis().SetTickLength(tick_frac_y)

        if draw_config.get_logy():
            self.curr_canvases['top'].SetLogy()

    def DrawSameCanvas(self, canvas, samples, draw_config, drawHist=False ) :

        canvas.cd()

        if not drawHist :
            drawHist = [0]*len(samples)

        ymin = draw_config.get_ymin()
        ymax = draw_config.get_ymax()
        ymax_scale = draw_config.get_ymax_scale()
        normalize = draw_config.get_normalize()

        calcymax = 0
        calcymin = 0.5
        for samp in samples :
            max = samp.hist.GetMaximum()
            min = samp.hist.GetMaximum()
            if normalize and samp.hist.Integral() != 0 :
                max = max / samp.hist.Integral()
                min = min / samp.hist.Integral()
            if max > calcymax :
                calcymax = max
            if min < calcymin :
                calcymin = min

        if ymax is None :
            ymax = calcymax
        if ymin is None :
            ymin = calcymin
        if ymax_scale is not None :
            ymax *= ymax_scale
        else :
            ymax *= 1.2
        ymin *= 0.8

        first = True
        for hist_name, hist_config in draw_config.hist_configs.iteritems() :

            draw_samp = self.get_samples( name=hist_name  )
            if draw_samp :
                draw_samp = draw_samp[0]
            else :
                draw_samp = None
                print 'WARNING Did not get a sample associated with the hist_config'

            drawcmd = 'same'
            if first :
                drawcmd = ''
                first = False
            if draw_samp is not None and draw_samp.isSignal  :
                drawcmd+='hist'

            
            if draw_samp is not None :

                draw_samp.hist.GetYaxis().SetTitle( draw_config.get_ylabel() )
                if not draw_config.doRatio()  :
                    draw_samp.hist.GetXaxis().SetTitle( draw_config.get_xlabel() )

                draw_samp.hist.SetLineColor( hist_config['color'] )
                draw_samp.hist.SetLineWidth( 2 )
                draw_samp.hist.SetMarkerSize( 1.1 )
                draw_samp.hist.SetMarkerStyle( 20 )
                draw_samp.hist.SetMarkerColor(hist_config['color'])
                draw_samp.hist.SetStats(0)

                if normalize and draw_samp.hist.Integral() != 0  :
                    draw_samp.hist.Scale(1.0/draw_samp.hist.Integral())

                draw_samp.hist.GetYaxis().SetRangeUser(ymin, ymax)

                #draw_samp.hist.Draw(drawcmd+'goff')
                draw_samp.hist.Draw(drawcmd)

    def CompareSelections( self, varexp, selections, reqsamples, histpars, hist_config={}, label_config={}, legend_config={}, same=False, useModel=False, treeHist=None, treeSelection=None ) :
        assert len(selections) == len(reqsamples), 'selections and samples must have same length'

        if 'colors' in hist_config :
            if len(hist_config['colors']) != len( selections ) :
                print 'Size of colors input does not match size of vars input!'

                hist_config['colors'] = [ self.get_samples(name=s)[0].color for s in reqsamples ]

        if self.collect_commands :
            self.add_compare_config( varexp, selections, reqsamples, histpars, hist_config=hist_config, label_config=label_config, legend_config=legend_config)
            return

        if not same :
            self.clear_all()

        config = DrawConfig( varexp, selections, histpars, samples=reqsamples, hist_config=hist_config, label_config=label_config, legend_config=legend_config )
        config.create_hist_configs()

        self.draw_commands.append(config)
        
        created_samples = self.MakeSameCanvas(config, preserve_hists=True, useModel=useModel, treeHist=treeHist, treeSelection=treeSelection )
        if not created_samples :
            print 'No histograms were created'
            return

        # make the legend
        step = len(created_samples)
        self.curr_legend = self.create_standard_legend(step, draw_config=config )

        legend_entries = config.get_legend_entries()
        self.create_same_legend( legend_entries , created_samples )

        #self.DrawCanvas(self.curr_canvases['same'], ylabel=ylabel, xlabel=xlabel, rlabel=rlabel, doratio=doratio, labelStyle=labelStyle, rmin=rmin, rmax=rmax, ymax=ymax, ymin=ymin, logy=logy, extra_label=extra_label, extra_label_loc=extra_label_loc)
        self.DrawCanvas(self.curr_canvases['same'], config )

    def create_same_legend(self,  legend_entries, created_samples ) :

        # check for an input legend_entries
        if not legend_entries : 
            legend_entries = [s.legendName for s in created_samples]

        for idx, samp in enumerate(created_samples) :
            drawopt = 'PL'
            if samp.isSignal :
                drawopt = 'L'
            legname = legend_entries[idx]
            self.curr_legend.AddEntry(samp.hist, legname,  drawopt)
            self.curr_legend.SetMargin(0.2)

    def Draw2D( self, varexp, selections, sample_names, histpars=None, drawopts='', xlabel=None, ylabel=None) :

        self.clear_hists()

        if not isinstance(sample_names, list) :
            sample_names = [sample_names]
        if not isinstance(selections, list) :
            selections = [selections]

        if len(sample_names) != len(selections) :
            print 'Length of samples does not match length of selections'

        created_samples=[]
        for idx, (samp_name, selection) in enumerate(zip(sample_names, selections)) :

            if not self.get_samples(name=samp_name) :
                print 'No sample with name ', samp_name
            samp = self.get_samples(name=samp_name)[0]
            newname = '%s_%d' %(samp.name, idx)
            newsamp = self.clone_sample( oldname=samp.name, newname=newname, temporary=True )
            newsamp.hist = None

            self.create_hist( newsamp, varexp, selection, histpars)

            if xlabel is not None :
                samp.hist.GetXaxis().SetTitle( xlabel )
            if ylabel is not None :
                samp.hist.GetYaxis().SetTitle( ylabel )
            
            created_samples.append(newsamp)

        for idx, samp in enumerate(created_samples) :

            self.curr_canvases['base%d'%idx] = ROOT.TCanvas('basecan%d'%idx, '')
            self.curr_canvases['base%d'%idx].SetRightMargin(0.12)
            self.curr_canvases['base%d'%idx].cd()

            samp.hist.Draw(drawopts)

    def CompareVars( self, varexps, selections, sample_names, histpars=None, hist_config={}, label_config={}, legend_config={}, same=False ) :

        self.clear_all()

        if not isinstance( varexps, list ) : varexps = [varexps]
        if not isinstance( selections, list ) : selections= [selections]
        if not isinstance( sample_names, list ) : sample_names= [sample_names]

        if len(selections) < len(varexps) and len(selections)==1 :
            selections = [selections[0]]*len(varexps)


        if 'colors' in hist_config :
            if len(hist_config['colors']) != len( selections ) :
                print 'Size of colors input does not match size of vars input!'

                hist_config['colors'] = [ self.get_samples(name=s)[0].color for s in reqsamples ]

        config = DrawConfig( varexps, selections, histpars, samples=sample_names, hist_config=hist_config, label_config=label_config, legend_config=legend_config )
        config.create_hist_configs()
        self.draw_commands.append(config)

        created_samples = self.MakeSameCanvas(config, preserve_hists=True )
        if not created_samples :
            print 'No histograms were created'
            return

        #samples = []
        #for sn in sample_names :
        #    samples.append( self.get_samples( name=sn )[0] )

        #if len(samples) < len(varexps) and len(samples)==1 :
        #    samples = [samples[0]]*len(varexps)

        #created_hists = []
        #for var, selection, sample in zip(varexps, selections, samples) :

        #    newname = '%s_%s' %(sample.name, var)
        #    newsamp = self.clone_sample( oldname=sample.name, newname = newname, temporary=True )
        #    created_hists.append(newsamp.name)
        #    self.create_hist( newsamp, var, selection, histpars)
        #    self.curr_hists[newsamp.name] = newsamp.hist.Clone(var)

        #if doratio :
        #    self.create_ratio_sample( 'ratio', num_sample=created_hists[0], den_sample=created_hists[1] )


        #self.curr_canvases['same'] = ROOT.TCanvas('same', '')
        #for idx, name in enumerate( created_hists ) :
        #    self.curr_hists[name].SetMarkerColor( colors[idx] )
        #    self.curr_hists[name].SetLineColor( colors[idx] )

        #    if idx == 0 :
        #        if normalize :
        #            self.curr_hists[name].DrawNormalized('hist')
        #        else :
        #            self.curr_hists[name].Draw('hist')
        #    else :
        #        if normalize :
        #            self.curr_hists[name].DrawNormalized('samehist')
        #        else :
        #            self.curr_hists[name].Draw('samehist')

        ##self.DrawCanvas(self.curr_canvases['same'])

        ## make the legend
        ## In placing the legend move the bottom down 0.05 for each entry
        ##step = len(varexps)
        ##self.curr_legend = self.create_standard_legend(step, doratio)

        ##for var, lab in zip( varexps, labels ) :
        ##    histname = '%s_%s' %(sample.name, var)
        ##    self.curr_legend.AddEntry(self.curr_hists[histname], lab, 'L' )
        ##    
        #self.DrawCanvas(self.curr_canvases['same'], ylabel=ylabel, xlabel=xlabel, doratio=doratio)

        # make the legend
        step = len(created_samples)
        self.curr_legend = self.create_standard_legend(step, draw_config=config )

        legend_entries = config.get_legend_entries()
        self.create_same_legend( legend_entries , created_samples )

        #self.DrawCanvas(self.curr_canvases['same'], ylabel=ylabel, xlabel=xlabel, rlabel=rlabel, doratio=doratio, labelStyle=labelStyle, rmin=rmin, rmax=rmax, ymax=ymax, ymin=ymin, logy=logy, extra_label=extra_label, extra_label_loc=extra_label_loc)
        self.DrawCanvas(self.curr_canvases['same'], config )

    #---------------------------------------
    def CompareRatios( self, varexp, selections, samples, histpars=None, normalize=False, doratio=False, colors = [], legend_entries=[], xlabel=None, ylabel=None, rlabel=None, rmin=None, rmax=None, ymin=None, ymax=None, logy=False ) :

        self.clear_all()

        assert  len(samples) == len(selections), 'sample and selection must have the same length'

        if len( selections )%2 != 0 :
            print 'Expect an even number of selecitions, even entries should be numerators, odd entries are denominators '

        if len(colors) != (len( samples )/2) :
            if colors :
                print 'Size of colors input does not match size of vars input!'
    
            colors = [ self.get_samples(name=s)[0].color for s in samples ]


        num_samps = []
        den_samps = []
        for idx, (selection, samplename) in enumerate(zip(selections, samples)) :

            sample = self.get_samples(name=samplename)[0]

            if idx%2 == 0 :
                newname = '%snum%d' %(sample.name, idx)
            else :
                newname = '%sden%d' %(sample.name, idx)

            newsamp = self.clone_sample( oldname=sample.name, newname=newname, temporary=True )
            self.create_hist( newsamp, varexp, selection, histpars)

            if normalize :
                newsamp.hist.Scale( 1.0/newsamp.hist.Integral() )
            
            if idx%2 == 0 :
                num_samps.append( newsamp )
            else :
                den_samps.append( newsamp )

        ratio_samps = []
        for num, den in zip( num_samps, den_samps ) :
            rsamp = self.create_ratio_sample( num.name+'ratio', num_sample=num, den_sample=den )
            rsamp.isRatio=False #trick it into not counting these as ratios
            ratio_samps.append(rsamp)
            
        if doratio :
            self.create_top_canvas_for_ratio('same')
        else :
            self.create_standard_canvas('same')

        self.DrawSameCanvas(self.curr_canvases['same'], ratio_samps, draw_config)

        
        if not legend_entries :
            legend_entries = [s.name for s in ratio_samps]

        self.curr_legend = self.create_standard_legend( len(ratio_samps) )
        for idx, samp in enumerate(ratio_samps) :
            drawopt = 'PL'
            if samp.isSignal :
                drawopt = 'L'
            legname = legend_entries[idx]
            self.curr_legend.AddEntry(samp.hist, legname,  drawopt)
            self.curr_legend.SetMargin(0.2)

            samp.hist.SetLineColor( colors[idx] )
            samp.hist.SetMarkerColor( colors[idx] )

        self.DrawCanvas(self.curr_canvases['same'], ylabel=ylabel, xlabel=xlabel, rlabel=rlabel, doratio=doratio, rmin=rmin, rmax=rmax, ymax=ymax, ymin=ymin, logy=logy)
    

    #---------------------------------------
    def DoTAndP( self, varexp, num_selection, den_selection, sample, histpars=None, binning=[], xlabel=None, ylabel=None, ymin=None, ymax=None, label=None, normalize=False, colors=[], legend_entries=[] ) :
    
        self.clear_all()

        if not isinstance(sample, list) :
            sample = [sample]
        if not isinstance(num_selection, list) :
            num_selection= [num_selection]
        if not isinstance(den_selection, list) :
            den_selection= [den_selection]
    
        assert  len(sample) == len(num_selection), 'sample and numerator selection must have the same length'
        assert  len(sample) == len(den_selection), 'sample and denominator selection must have the same length'
    
        if len(colors) != len( sample ) :
            if colors :
                print 'Size of colors input does not match size of vars input!'
    
            colors = [ self.get_samples(name=s)[0].color for s in sample ]
    
        num_hists = {}
        den_hists = {}
    
    
        # create a new sample manager to hold histograms
        #sman = SampleManager('', '' )
        for idx, (num, den, insamp) in enumerate( zip(num_selection, den_selection, sample) ) :
            use_stored_first = False
            if idx > 0 :
                use_stored_first = True
    
            samp = insamp + str(idx)

            num_name = '%s_num' %samp
            den_name = '%s_den' %samp

            self.clone_sample(oldname=insamp, newname=num_name, temporary=True )
            self.clone_sample(oldname=insamp, newname=den_name, temporary=True )
    
            created_num_samp = []
            created_den_samp = []
            created_num_samp += self.MakeSameCanvas([num_name], varexp, num, histpars, doratio=False, useStoredBinning=use_stored_first, preserve_hists=True )
            num_samp = created_num_samp[0]
            #sman.samples.append(num_samp)

            created_den_samp += self.MakeSameCanvas([den_name], varexp, den, histpars, doratio=False, useStoredBinning=True, preserve_hists=True)
            den_samp = created_den_samp[0]
            #sman.samples.append( den_samp )
            #self.clear_hists()

            if normalize :
                num_samp.hist.Scale( 1.0/num_samp.hist.Integral() )
                den_samp.hist.Scale( 1.0/den_samp.hist.Integral() )

            rsamp = self.create_ratio_sample( '%s_ratio' %samp, num_sample=num_samp, den_sample=den_samp )
            if xlabel is not None :
                rsamp.hist.GetXaxis().SetTitle( xlabel )
                #rsamp.hist.GetXaxis().SetTitleSize( 0.05 )
                #rsamp.hist.GetXaxis().SetTitleOffset( 1.1 )
            if ylabel is not None :
                rsamp.hist.GetYaxis().SetTitle( ylabel )
                #rsamp.hist.GetYaxis().SetTitleSize( 0.05 )
                #rsamp.hist.GetYaxis().SetTitleOffset( 1.25 )
    
        plot_ymax = 0.0
        plot_ymin = 100.0
    
        for ratio in self.get_samples(isRatio=True) :
            ratio.hist.SetTitle('')
            max = ratio.hist.GetMaximum()
            min = ratio.hist.GetMinimum()
            if max > plot_ymax :
                plot_ymax = max
            if min < plot_ymin :
                plot_ymin = min
        if ymin is not None :
            plot_ymin = ymin
        else :
            plot_ymin *= 0.9
        if ymax is not None :
            plot_ymax = ymax
        else :
            plot_ymax *= 1.1
    
    
        for ratio in self.get_samples(isRatio=True) :
            ratio.hist.GetYaxis().SetRangeUser( plot_ymin, plot_ymax )
    
        self.create_standard_canvas( )
        self.curr_canvases['base'].SetBottomMargin(0.13)
        self.curr_canvases['base'].SetLeftMargin(0.13)
        self.curr_canvases['base'].SetTitle('')
    
        for idx, ( ratio, color ) in enumerate( zip(self.get_samples(isRatio=True), colors ) ) :
            drawcmd = ''
            if idx > 0 :
                drawcmd = 'same'
            if varexp.count(':') :
                drawcmd += 'colz'
            ratio.hist.SetLineColor( color )
            ratio.hist.SetMarkerColor( color )
            ratio.hist.SetMarkerStyle( 20 )
            ratio.hist.SetMarkerSize( 1.2 )
            ratio.hist.Draw(drawcmd)
            self.add_decoration( ratio.hist )
    
        leg_xmin = 0.6
        leg_xmax = 0.8
        leg_ymin = 0.68
        leg_ymax = 0.88
    
        if varexp.count('eta') :
            leg_xmin -= 0.2
            leg_xmax -= 0.2
            leg_ymin -= 0.1
            leg_ymax -= 0.1
    
        leg = ROOT.TLegend( leg_xmin, leg_ymin, leg_xmax,  leg_ymax)
        leg.SetBorderSize(0)
        leg.SetFillColor(ROOT.kWhite)
    
        # check for an input legend_entries
        if not legend_entries : 
            legend_entries = sample
    
        for idx, (inname, legentry) in enumerate( zip(sample, legend_entries) ) :
            name = inname+str(idx)
            ratio = self.get_samples(name='%s_ratio' %name )[0]
            leg.AddEntry( ratio.hist, legentry, 'LPE' )
    
        self.curr_canvases['base'].cd()
        leg.Draw()
        self.curr_legend = leg
    
        if label is not None :
            lab = ROOT.TLatex()
            lab.SetNDC()
    
            lab.SetText( 0.2, 0.83, label )
            lab.SetTextFont(42)
    
            lab.Draw()
    
            if varexp.count('eta') :
                xsize = lab.GetXsize()
                newxmin = 0.5-( xsize/2 )
                lab.SetX(newxmin)
            self.curr_decorations .append(lab)

    def MakeRocCurve( self, varexps, selections, signal, background , histpars=None, colors=[], markers=[], legend_entries=[], debug=False, less_than=False, doSoverB=False, extra_label=None, extra_label_loc=None, ymin=None, ymax=None, legendConfig=None ) :

        if not isinstance( selections, list ) :
            selections = [selections]
        if not isinstance( varexps, list ) :
            varexps = [varexps]
        if not isinstance( signal, list ) :
            signal = [signal]
        if not isinstance( background, list ) :
            background = [background]
        if not isinstance( histpars, list ) :
            histpars = [histpars]
        if not isinstance( less_than, list ) :
            less_than = [less_than]

        self.clear_all()

        self.apply_lenged_conf( legendConfig )

        sig_samps = []
        bkg_samps = []
        for sig in signal :
            sig_samps.append(self.get_samples(name=sig)[0])
        for bkg in background :
            bkg_samps.append(self.get_samples(name=bkg)[0])

        if not colors :
            colors = [ s.color for s in sig_samps ]

        if not markers :
            markers = [3]*len(signal)

        if len(histpars) != len(signal) :
            histpars = histpars*len(signal)

        assert len(sig_samps) == len(bkg_samps), 'Did not retrieve all signal or background samples'

        created_hists = []

        self.create_standard_canvas()
        self.curr_canvases['base'].cd()

        for idx, (varexp, selection, histpar, sig_samp, bkg_samp, color, marker, lt) in enumerate( zip( varexps, selections, histpars, sig_samps, bkg_samps, colors, markers, less_than) ) :

            newsigsamp = copy.copy(sig_samp)
            newsigsamp.hist = None
            newsigsamp.name = '%s%d' %(sig_samp.name, idx)
            newbkgsamp = copy.copy(bkg_samp)
            newbkgsamp.hist = None
            newbkgsamp.name = '%s%d' %(bkg_samp.name, idx)

            self.create_hist( newsigsamp, varexp, selection, histpar)
            self.create_hist( newbkgsamp, varexp, selection, histpar)

            sig_eff = []
            bkg_eff = []
            nsig = []
            nbkg = []
            sig_tot = newsigsamp.hist.Integral()
            bkg_tot = newbkgsamp.hist.Integral()
            for bin in range( 1, histpar[0]+1 ) :
                min = bin
                max = histpar[0]
                if lt :
                    min = 0 
                    max = bin
                sig_eff.append( newsigsamp.hist.Integral( min, max )/sig_tot )
                bkg_eff.append( newbkgsamp.hist.Integral( min, max )/bkg_tot )
                nsig.append( newsigsamp.hist.Integral( min, max ))
                nbkg.append( newbkgsamp.hist.Integral( min, max ))

            name='roc%d' %idx

            nbins = histpar[0]
            # if there are only two bins, this is a boolean and the TGraph should only have one point
            if nbins == 2 :
                nbins = 1

            self.curr_hists[name] =  ROOT.TGraph(nbins)
            self.curr_hists[name].SetName(name)

            if ymin is not None :
                self.curr_hists[name].SetMinimum( ymin )
            if ymax is not None :
                self.curr_hists[name].SetMaximum( ymax )

            if doSoverB :
                for bin, (sigeff, sig, bkg) in enumerate(zip( sig_eff, nsig, nbkg) ):
                    if sig == 0 and bkg == 0 :
                        continue
                    if histpar[0] <= 2 and sigeff==1 :
                        continue;
                    if debug : print 'Place point at %f, bin %d, nSig=%f, nBkg=%f, %f, %f ' %( newsigsamp.hist.GetBinLowEdge( bin+1 ), bin+1, sig, bkg, sigeff, sig/math.sqrt(sig+bkg))
                    self.curr_hists[name].SetPoint( bin, sigeff, sig/math.sqrt(sig+bkg) );

            else :
                bin=0
                for (sig, bkg) in zip( sig_eff, bkg_eff ):
                    if histpar[0] <= 2 and sig==1 :
                        continue;
                    bin+=1
                    if debug : print 'Place point at %f, bin %d : %f, %f ' %( newsigsamp.hist.GetBinLowEdge( bin ), bin, sig, 1.0-bkg)
                    self.curr_hists[name].SetPoint( bin, sig, 1.0-bkg );

            #store numerical results

            for bin, (sigeff, bkgeff, sig, bkg) in enumerate(zip( sig_eff, bkg_eff, nsig, nbkg) ):
                den = sig+bkg
                if den != 0 :
                    self.transient_data.setdefault( 'Sig: %s, Bkg: %s' %(sig_samp.name, bkg_samp.name), [] ).append( { 'CutVal' : newsigsamp.hist.GetBinLowEdge( bin+1 ), 'nSig' : sig, 'nBkg' : bkg, 'bkgEff' : bkgeff, 'sigEff': sigeff, 'SoverRootSplusB' : sig/math.sqrt(den) } )


            self.curr_hists[name].SetLineColor(color)
            self.curr_hists[name].SetMarkerColor(color)
            self.curr_hists[name].SetMarkerStyle(marker)
            self.curr_hists[name].SetMarkerSize(1.2)

            drawcmd = ''
            if idx == 0 :
                #drawcmd = 'ACP'
                drawcmd = 'AP'
                if histpar[0]>2 : drawcmd = 'C'+drawcmd
            else :
                #drawcmd = 'CPsame'
                drawcmd = 'Psame'
                if histpar[0]>2 : drawcmd = 'C'+drawcmd

            self.curr_hists[name].Draw(drawcmd)
            self.curr_hists[name].GetXaxis().SetTitle('Signal Efficiency')
            if doSoverB :
                self.curr_hists[name].GetYaxis().SetTitle('S / #sqrt{ S + B }')
            else :
                self.curr_hists[name].GetYaxis().SetTitle('Background Rejection')

            self.curr_hists[name].GetYaxis().SetTitleSize(0.05)
            self.curr_hists[name].GetXaxis().SetTitleSize(0.05)
            self.curr_hists[name].GetYaxis().SetLabelSize(0.05)
            self.curr_hists[name].GetXaxis().SetLabelSize(0.05)
            self.curr_hists[name].SetTitle()
            created_hists.append(self.curr_hists[name])

        if legend_entries :

            if len(legend_entries) != len(self.curr_hists) :
                print 'Length of lengend entries does not match length of created histograms'
            else :

                self.curr_legend = self.create_standard_legend( len(sig_samps) )
                for hist, entry in zip( created_hists, legend_entries) :
                    self.curr_legend.AddEntry(hist, entry, 'PE')

                self.curr_legend.Draw()

        if extra_label is not None :
            self.place_extra_label( extra_label, extra_label_loc )


    # ------------------------------------------------------------
    #   Do variable rebinning for a stack plot
    # ------------------------------------------------------------
    def make_variable_binning(self, stacklist,threshold=3):

        if not isinstance(stacklist, list) :
            stacklist = [stacklist]
    
        # sum stack list
        sum = stacklist[0].Clone('sum')
        sum.Reset()
        for h in stacklist:
            sum.Add(h)
    
        # make binning
        bins=[]
        axis=sum.GetXaxis()
        bins.append(axis.GetXmin())
        count=0
        for b in range(1, sum.GetNbinsX()+1):
            # this special case is to not extend the first 
            # filled bin to the edge of the histogram
            if sum.GetBinContent(b)>0 and count==0 and len(bins)==1:
                bins.append(axis.GetBinLowEdge(b))
            count+=sum.GetBinContent(b)
            if count>threshold:
                bins.append(axis.GetBinUpEdge(b))
                count=0
        if count!=0:
            bins.append(axis.GetXmax())
        print bins,count
        return bins
    
    # ------------------------------------------------------------
    #   Do variable rebinning for a stack plot
    # ------------------------------------------------------------
    def do_variable_rebinning(self, samp,bins):

        if isinstance( samp, Sample ) :
            oldhist = samp.hist
        if isinstance( samp, ROOT.TH1 ) :
            oldhist = samp
        newhist=ROOT.TH1F(oldhist.GetName()+"_rebin",
        oldhist.GetTitle()+";"+oldhist.GetXaxis().GetTitle()+";"+oldhist.GetYaxis().GetTitle(),len(bins)-1,array('d',bins))
        a=oldhist.GetXaxis()
        newa=newhist.GetXaxis()
        for b in range(1, oldhist.GetNbinsX()+1):
            newb=newa.FindBin(a.GetBinCenter(b))
            val=newhist.GetBinContent(newb)
            err=newhist.GetBinError(newb)
            ratio_bin_widths=newa.GetBinWidth(newb)/a.GetBinWidth(b)
            val=val+oldhist.GetBinContent(b)/ratio_bin_widths
            err=math.sqrt(err*err+oldhist.GetBinError(b)/ratio_bin_widths*oldhist.GetBinError(b)/ratio_bin_widths)
            newhist.SetBinContent(newb,val)
            newhist.SetBinError(newb,err)
    
        return newhist

    # ----------------------------------------------------------------------------
    def create_standard_legend(self, nentries,draw_config=None ) :

        legend_config = {}
        if draw_config is not None :
            legend_config = draw_config.legend_config

        legendTranslateX = legend_config.get('legendTranslateX', self.legendTranslateX )
        legendTranslateY = legend_config.get('legendTranslateY', self.legendTranslateY )
        
        legendWiden      = legend_config.get('legendWiden'   , self.legendWiden    )
        legendCompress   = legend_config.get('legendCompress', self.legendCompress )

        legendLoc        = legend_config.get('legendLoc', self.legendLoc )

        entryWidth       = legend_config.get('entryWidth', self.entryWidth )
        

        if legendLoc == 'TopLeft' :
            legend_limits = { 'x1' : 0.2+legendTranslateX, 'y1' : 0.88-legendCompress*entryWidth*nentries+legendTranslateY, 'x2' : 0.5*legendWiden+legendTranslateX, 'y2' : 0.88+legendTranslateY }
        elif legendLoc == 'Double' :
            legend_limits = { 'x1' : 0.15+legendTranslateX, 'y1' : 0.90-legendCompress*entryWidth*nentries+legendTranslateY, 'x2' : 0.65*legendWiden+legendTranslateX, 'y2' : 0.85+legendTranslateY }
        else :
            legend_limits = { 'x1' : 0.9-0.25*legendWiden+legendTranslateX, 'y1' : 0.90-legendCompress*entryWidth*nentries+legendTranslateY, 'x2' : 0.90+legendTranslateX, 'y2' : 0.90+legendTranslateY }

        # modify for different canvas size
        if draw_config.doRatio():
            legend_limits['y1'] = legend_limits['y1']*0.90
            #legend_limits['y2'] = legend_limits['y2']*1.05

        # grab stored legend limits
        if self.legendLimits :
            legend_limits = self.legendLimits

        leg = ROOT.TLegend(legend_limits['x1'], legend_limits['y1'],
                           legend_limits['x2'], legend_limits['y2'])
        leg.SetFillColor(ROOT.kWhite)
        leg.SetBorderSize(0)

        if legendLoc == 'Double' :
            leg.SetNColumns(2)
        

        return leg
    
    # ----------------------------------------------------------------------------
    def store_current_legend_placement(self) :
    
        self.legendLimits = {}
    
        leg = self.curr_legend
        self.legendLimits['x1'] = leg.GetX1NDC()
        self.legendLimits['y1'] = leg.GetY1NDC() 
        self.legendLimits['x2'] = leg.GetX2NDC()
        self.legendLimits['y2'] = leg.GetY2NDC()

        print self.legendLimits
        
    # ----------------------------------------------------------------------------
    def outputExists(self, name, dir) :
        exists = False
        if os.path.isdir(dir) :
            for file in os.listdir( dir ) :
                if file.count(name) :
                    exists=True

        return exists

    # ----------------------------------------------------------------------------
    def addHistsWeightedAvg( self, dest_hist, add_hist ) :

        def _add_generic_bin( dest_hist, add_hist, ibin ) :
            dest_val = dest_hist.GetBinContent( ibin )
            dest_err = dest_hist.GetBinError( ibin )

            add_val  = add_hist.GetBinContent( ibin )
            add_err  = add_hist.GetBinError( ibin )

            if dest_val == 0 :
                avg = ufloat( add_val, add_err )
            elif add_val == 0 :
                avg = ufloat( dest_val, dest_err )
            else :
                dest_uf = ufloat( dest_val, dest_err )
                add_uf  = ufloat( add_val , add_err )

                dest_w = 1./( dest_err * dest_err )
                add_w   = 1./( add_err * add_err )

                avg = (dest_w * dest_uf + add_w * add_uf) / ( dest_w + add_w )

                print 'bin = ', ibin
                print 'dest = ${:.2ufL}$' .format( dest_uf )
                print 'add = ${:.2ufL}$' .format( add_uf )
                print 'avg = ${:.2ufL}$' .format( avg )

            dest_hist.SetBinContent( ibin, avg.n )
            dest_hist.SetBinError( ibin, avg.s )



        if isinstance( dest_hist, ROOT.TH3 ) :
            if dest_hist.GetNbinsX() != add_hist.GetNbinsX() or dest_hist.GetNbinsY() != add_hist.GetNbinsY() or dest_hist.GetNbinsZ() != add_hist.GetNbinsZ() :
                print 'addHistsWeightedAvg: ERROR Bins must have the same lengths!  Exiting!'
                return

            for xbin in range( 1, dest_hist.GetNbinsX()+1 ) :
                for ybin in range( 1, dest_hist.GetNbinsY()+1 ) :
                    for zbin in range( 1, dest_hist.GetNbinsZ()+1 ) :
                        ibin = dest_hist.GetBin( xbin, ybin, zbin )
                        _add_generic_bin( dest_hist, add_hist, ibin )

        elif isinstance( dest_hist, ROOT.TH2 ) :
            if dest_hist.GetNbinsX() != add_hist.GetNbinsX() or dest_hist.GetNbinsY() != add_hist.GetNbinsY() :
                print 'addHistsWeightedAvg: ERROR Bins must have the same lengths!  Exiting!'
                return

            for xbin in range( 1, dest_hist.GetNbinsX()+1 ) :
                for ybin in range( 1, dest_hist.GetNbinsY()+1 ) :
                    ibin = dest_hist.GetBin( xbin, ybin )
                    _add_generic_bin( dest_hist, add_hist, ibin )

        elif isinstance( dest_hist, ROOT.TH1 ) :
            if dest_hist.GetNbinsX() != add_hist.GetNbinsX() :
                print 'addHistsWeightedAvg: ERROR Bins must have the same lengths!  Exiting!'
                return

            for xbin in range( 1, dest_hist.GetNbinsX()+1 ) :
                ibin = dest_hist.GetBin( xbin )
                _add_generic_bin( dest_hist, add_hist, ibin )

