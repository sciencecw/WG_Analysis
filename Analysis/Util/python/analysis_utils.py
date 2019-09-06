import os

class Printer() :
    """ print tabulated data with variable column width """

    def __init__ (self) :
        self.entries = []

    def AddLine( self, line ) :
        self.entries.append( line )

    def Print(self) :
        message = self.GetMessage()
        for msg in message:
            print msg
        return

    def GetMessage(self):
        # get number of columns
        max_cols = max( [len(line) for line in self.entries ] )
        colwidths = [0]*max_cols

        # collect max length of each column
        for line in self.entries :
            for idx, col in enumerate(line) :
                if len(col) > colwidths[idx] :
                    colwidths[idx] = len(col)

        #for line in self.entries :
        #    print ' '.join( [col.ljust(colwidths[idx]) for idx, col in enumerate(line) ] )
        return   [ ' '.join( [col.ljust(colwidths[idx]) for idx, col in enumerate(line) ] )  for line in self.entries ]

def read_xsfile( xsfile, lumi, print_values=False ) :
    if print_values:
        print '-------------------------------------'
        print ' LOAD CROSS SECTION INFO'
        print '-------------------------------------'
    weightMap = {}
    if xsfile is None :
        return weightMap
    if lumi is None :
        print 'Cannot calculate weights without a luminosity'
        return weightMap
    if not os.path.isfile( xsfile ) :
        print 'Could not locate cross section file.  No values will be loaded.'
        raise RuntimeError ## assume we always need one
        return weightMap

    ofile = open( xsfile )
    xsdict = eval( ofile.read() )
    xs_printer = Printer()

    for name, values in xsdict :

        values["lumi"] = lumi ## add lumi value for scale_calc() to unpack)
        scale = scale_calc(**values)
        xs_printer.AddLine( ['Sample %s ' %name, 'cross section : %.4g' %values['cross_section'], 'pb  N Events : %d' %(values['n_evt']), 'Scale : %.4g' %scale ] )

        weightMap[name] = {}  
        weightMap[name]['scale'] = scale
        weightMap[name]['cross_section'] = values['cross_section']
        weightMap[name]['n_evt'] = values['n_evt']
        weightMap[name]['k_factor'] = values['k_factor']
        weightMap[name]['gen_eff'] = values['gen_eff']

    if print_values :
        xs_printer.Print()

    return weightMap, xs_printer


def scale_calc(cross_section, lumi,  n_evt, gen_eff = 1.0, k_factor = 1.0, **kwargs):
    """ calculate scaling of MC samples for a given lumi
        n_evt: total number of events in MC sample before selection
        * ignores extra keyword arguments
    """
    return nevents_calc(cross_section, lumi, gen_eff,  k_factor) / n_evt


def nevents_calc(cross_section, lumi, gen_eff = 1.0, k_factor = 1.0 , **kwargs):
    """ calculate expectant normalization of MC events for a given lumi """
    return cross_section * lumi * gen_eff * k_factor



