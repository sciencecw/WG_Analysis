import os
import pdb
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
from functools import wraps

from uncertainties import ufloat
from uncertainties import umath
import pickle
import core
import subprocess
import multiprocessing
import collections

ROOT.gROOT.SetBatch(False)

ROOT.gStyle.SetPalette(1)


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

    def get_doratio(self) :
        """ return doratio key value set in hist_config """
        return self.hist_config.get('doratio', False)

    def get_reverseratio(self) :
        """ return reverseratio key value set in hist_config """
        return self.hist_config.get('reverseratio', False)

    def get_binomunc(self) :
        """ return binomunc key value set in hist_config, default is False """
        return self.hist_config.get('binomunc', False)

    def get_drawhist(self) :
        return self.hist_config.get('drawhist', False)

    def get_ylabel(self) :
        """ return y axis label set in hist_config
            if not set, returns Events/ x GeV with bin width calculated from histpars"""
        ylabel = self.hist_config.get('ylabel', None) 
        if ylabel is None :
            if isinstance( self.histpars, tuple ) :
                bin_width = ( self.histpars[2] - self.histpars[1] )/self.histpars[0]
                bin_width_f = ( self.histpars[2] - self.histpars[1] )/float(self.histpars[0])
            else :
                bin_width = 1
                bin_width_f = 1
                

            xunit = self.hist_config.get('xunit', 'GeV') 
            if math.fabs(bin_width_f - bin_width) != 0 :
                ylabel = 'Events / %.1f %s' %(bin_width_f,xunit)
            else :
                ylabel = 'Events / %d %s' %(bin_width, xunit)

        return ylabel

    def get_xlabel(self) :

        xunit = self.hist_config.get('xunit', 'GeV') 
        return self.hist_config.get('xlabel',"" ) + '[%s]' %xunit

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

    def get_weight( self ) :
        """ defaults to empty string: no weights
            string: Apply weights
        """
        return self.hist_config.get('weight', "")

    def get_unblind( self ) :
        """ defaults to False: skipping Data
            1, True: show Data
            string: Apply blinding selection
            This setting only has effect on Data
        """
        return self.hist_config.get('unblind', False)

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
            statlabel  .SetNDC()
            statlabel  .SetTextSize(0.045)
            statlabel.SetText(0.16, 0.93, stattext)

            labels.append(statlabel)

        labelStyle = self.label_config.get('labelStyle', None)
        labelLoc = self.label_config.get('labelLoc', None)
        if labelStyle is None:

            text_dx = self.label_config.get("dx",0)
            text_x = text_dx +0.18
            text_dy = self.label_config.get("dy",0)
            text_y = text_dy +0.88

            if labelLoc == 'topright' :
                text_x = 0.75
                text_y = 0.93
            cmslabel = ROOT.TLatex()
            cmslabel.SetNDC()
            cmslabel.SetTextSize( 0.05 )
            cmslabel.SetText(text_x, text_y, 'CMS')
            wiplabel = ROOT.TLatex()
            wiplabel.SetNDC()
            wiplabel.SetTextSize( 0.05 )
            wiplabel.SetText(text_x+0.07, text_y, 'Simulation Work in Progress')
            wiplabel.SetTextFont(52)
            labeltext = '36 fb^{-1} (13 TeV)'
            rootslabel = ROOT.TLatex()
            rootslabel.SetText(text_dx+0.75, text_dy+0.88, labeltext  )
            rootslabel.SetTextFont(42)
            rootslabel .SetNDC()
            rootslabel .SetTextSize(0.045)

            labels.append(rootslabel)
            labels.append(cmslabel)
            labels.append(wiplabel)

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



    ####################################################


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
        """  Fills hist_configs ordered dictionary
        Each sample gets one dictionary entry with var, selection, sample, color ,legend entry according to order of input
        Argument:
            branches: ??
        """

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
        """
        Initialize histogram
        Parameters:
            name: histogram name
        Returns:
            TH1-3, depending on length of histpars
        """

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
