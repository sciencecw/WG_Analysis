import ROOT
from ROOT import RooFit
from uncertainties import ufloat
import uuid
import re
from collections import namedtuple, OrderedDict
from functools import wraps
from DrawConfig import DrawConfig
from pprint import pprint


ROOT.gStyle.SetPalette(ROOT.kBird) 

def f_Obsolete(f):
        @wraps(f)
        def f_wrapper(*args, **kws):
                print f.__name__,": This method is obsolete"
                return f(*args,**kws)
        return f_wrapper

xname = "Reco mass m(#gamma,e)"
class FitManager : 
    """ Aim to collect all fitting machinery here """

    ParamDCB = ["dcb_mass"   ,"dcb_sigma" ,
               "dcb_alpha1" ,"dcb_power1",
               "dcb_alpha2" ,"dcb_power2"]
    setuparray = {
                    "cb":
                        [("cb_mass"  ,"Mass"  ,90,80,100),
                         ("cb_sigma" ,"Sigma" ,1,0.1,100),
                         ("cb_alpha" ,"Alpha" ,-1,-10,10),
                         ("cb_power" ,"Power" ,2,0,10),
                        ],
                    "dcbp": # reparametrize IC
                        [
                        ("x"         ,20,200),
                        ("dcb_mass"  ,91,87,95),
                        ("dcb_sigma" ,3,1,10),
                        ("dcb_alpha1",1,0.5,3),
                        ("dcb_power1",1,1,50),
                        ("dcb_alpha2",1,0.5,3),
                        ("dcb_power2",2,0.5,8),
                        ],
                    "dcb":
                        [
                        ("x"         ,20,200),
                        ("dcb_mass"  ,91,87,95),
                        ("dcb_sigma" ,3,1,10),
                        ("dcb_alpha1",1,0.2,3),
                        ("dcb_power1",1,1,100),
                        ("dcb_alpha2",1,0.2,3),
                        ("dcb_power2",2,0.5,8),
                        ],
                    "expo":
                        [
                        ('x',),#no redefinition
                        ('expo_c',0,-0.1,0.01)
                        #('expo_c',0,), # const
                        ],
                    "gauzz": # gaussian modelling Z peak
                        [
                        ("x",20,200),
                        ('gauzz_mean',90,80,100),
                        ('gauzz_sig' ,3,0.5,10),
                        ],
                    "gauszg": # model Z gamma
                        [
                        ('x',),#no redefinition
                        ('gauszg_mean', 85,60,90),
                        ('gauszg_sig' , 6,5,15),
                        ],
                    "gaus":
                        [
                        ('x',),#no redefinition
                        ('gaus_mean',100,-100,200),
                        ('gaus_sig' , 60,50,200),
                        ],
                    "simul":
                        [
                        ("Nsig",100000,10000,10000000),
                        ("Nbkg", 10000, 1000,10000000),
                        ],
                    "simulzg":
                        [
                        ("Nz"  ,100000,  100,10000000),
                        ("Nzg", 100000,  100,10000000),
                        ]
                 }
    LineDefs = [ 
            [RooFit.LineColor(ROOT.kRed)],
            [RooFit.LineColor(ROOT.kBlue),   RooFit.LineStyle(ROOT.kDashed),],
            [RooFit.LineColor(ROOT.kCyan+1), RooFit.LineStyle(ROOT.kDashed),],
               ]
    def __init__(self, fname, hist=None, xvardata = None, #,xvarfit=None,
                    label="", sample_params={}) :

        # initialize object container
        self.defs = OrderedDict()          # store RooRealVar for fit functions
        self.func = None
        self.func_pdf = None
        self.fit_params = {}
        self.fitrange = None
        self.loaded = 0
        self.curr_decorations = []

        # copy input arguments
        self.func_name = fname
        self.label = label

        self.xvardata = xvardata  # data plot range
        #self.xvarfit = xvarfit    # fit range
        # if x range is not given, then make one
        if self.xvardata == None:
                self.xvardata = ROOT.RooRealVar("x",xname,0,200,"GeV")
        if isinstance(self.xvardata,tuple):
                self.xvardata = ROOT.RooRealVar("x",xname,*self.xvardata)
        #if self.xvarfit == None:
        #        self.xvarfit = self.xvardata 
        self.pdfplotrange = False

        # toggles and extra parameters
        self.sample_params = sample_params

        self.canvas = None
        self.curr_canvases = {}
        self.objs = []
        self.addhist(hist) #default hist
        self.wk = None

    def addhist(self,hist,name="datahist", bkgd = False):
        if hist:
            assert isinstance(hist,ROOT.TH1)
            his = hist.Clone()
            # histogram style
            his.SetMarkerStyle( 20 )
            his.SetMarkerSize( 1.0 )
            if bkgd == False:
                his.SetLineColor( ROOT.kBlack )
                his.SetMarkerColor( ROOT.kBlack )
                self.hist = his
                # make datahist
                self.datahist = ROOT.RooDataHist( '%s%s' %(self.label,name), 'data', 
                            ROOT.RooArgList(self.xvardata), his )

                self.datahistlist['data'] = self.datahist
                ROOT.SetOwnership( self.datahist, False )
            else:
                his.SetLineColor( ROOT.kBlue-10 )
                his.SetMarkerColor( ROOT.kBlue-10 )
                self.histlist[bkgd] = his
                self.datahistlist[bkgd] = ROOT.RooDataHist( '%s%s' %(self.label,name), 'bkgd', 
                                        ROOT.RooArgList(self.xvardata), his )
                ROOT.SetOwnership(self.datahistlist[bkgd], False )
        else: 
            self.hist = None #ROOT.TH1F()
            self.datahist = None
            self.histlist={}
            self.datahistlist={}

    @f_Obsolete
    def MakeROOTObj( self, root_obj, *args ) :
        """ Generic function for making ROOT objects."""  

        try :
            thisobj = getattr(ROOT, root_obj)( *args )
            ROOT.SetOwnership( thisobj, False )
            return thisobj

        except TypeError :
            print '***********************************************************'
            print 'FitManager.MakeROOTObj -- Failed to create a %s.  Please check the arguments :'%root_obj
            print args
            print 'Exception is below'
            print '***********************************************************'
            raise


    @f_Obsolete ## FIXME
    def Integral( self ) :

        err = ROOT.Double()
        val = self.hist.IntegralAndError( self.hist.FindBin( self.xvar.getMin() ), self.hist.FindBin( self.xvar.getMax() ), err )

        return ufloat( val, err )

    def MakeIntegral(self,intrange=None, pdfname = None):
        """ make pdf integral that returns in range [0,1] """
        if intrange is None:
            intrange = self.fitrange
        if pdfname is None:
            pdf = self.func_pdf
        else:
            pdf = self.wk.pdf(pdfname)
        self.integral = pdf.createIntegral(self.xvardata,RooFit.NormSet(self.xvardata),RooFit.Range(*intrange))
        return self.integral

    def get_integral_value(self, intrange=None, pdfname = None):
        """ calculate pdf integral return [0,1] """
        ### FIXME
        pass

    def add_vars( self, arg_list) :
        if self.func_name == 'dijet' : 
            for i in range( 1, self.func_norders+1 ) :
                short_name = 'power' 
                if i > 1 :
                    short_name = 'logcoef%d' %(i-1)

                long_name = '%s_order%d_%s' %( self.func_name, i, self.label )
                this_def = self.get_vals( self.func_name, i )
                var = ROOT.RooRealVar( long_name, long_name, *this_def)
                ROOT.SetOwnership(var, False)
                arg_list.add( var  ) 
                self.fit_params[short_name] = var
        return

    @f_Obsolete
    def add_vars_old(self):
        pass
        # These need to be modified to match the format above
        #if self.func_name == 'atlas' : 
        #    def_num_power = self.get_vals( self.func_name+'_num_power', 1 )
        #    def_den_power = self.get_vals( self.func_name+'_den_power', 1 )

        #    var_num_power = ROOT.RooRealVar( 'num_power', 'num_power', def_num_power[0], def_num_power[1], def_num_power[2] )
        #    var_den_power = ROOT.RooRealVar( 'den_power', 'den_power', def_den_power[0], def_den_power[1], def_den_power[2] )
        #    ROOT.SetOwnership(var_num_power, False)
        #    ROOT.SetOwnership(var_den_power, False)
        #    arg_list.add( var_num_power )
        #    arg_list.add( var_den_power )
        #    for i in range( 1, self.func_norders+1 ) :
        #        def_den_logcoef = self.get_vals( self.func_name+'_den_logcoef', i )
        #        var_den_locoef  = ROOT.RooRealVar( 'den_logcoef_order%d' %i, 'den_logcoef_order%d' %i, def_den_logcoef[0], def_den_logcoef[1], def_den_logcoef[2] )
        #        ROOT.SetOwnership(var_den_locoef, False)
        #        #var.SetName( '%s_order%d' %( self.func_name, i ) )
        #        arg_list.add( var_den_locoef  ) 

        #if self.func_name == 'power' : 
        #    for i in range( 1, self.func_norders+1 ) :
        #        this_def_coef = self.get_vals( self.func_name+'_coef', i )
        #        this_def_pow  = self.get_vals( self.func_name+'_pow', i )

        #        var_coef = ROOT.RooRealVar( 'coef%d' %i, 'coef%d'%i, this_def_coef[0], this_def_coef[1], this_def_coef[2] )
        #        ROOT.SetOwnership(var_coef, False)
        #        var_pow = ROOT.RooRealVar( 'pow%d' %i, 'pow%d'%i, this_def_pow[0], this_def_pow[1], this_def_pow[2] )
        #        ROOT.SetOwnership(var_pow, False)

        #        arg_list.add( var_coef )
        #        arg_list.add( var_pow )


    def get_fit_function( self, forceUseRooFit=False ) :

        function = ''

        if self.func_name == 'dijet' : 
            order_entries = []

            log_str = 'TMath::Log10(@0/13000)'
            for i in range( 1, self.func_norders) :
                order_entries.append('@'+str(i+1)+'*' + '*'.join( [log_str]*i))
            
            function = 'TMath::Power( @0/13000., @1 + ' + '+'.join( order_entries) + ')'


        if self.func_name == 'atlas' : 

            function = 'TMath::Power( (1-(@0/13000.)), @1 ) / ( TMath::Power( @0/13000. , @2+ '
            order_entries = []
            for i in range( 0, self.func_norders ) :
                order_entries.append( '@%d*TMath::Log10( @0/13000.)'  %( i+3 ) )

            function += ' + '.join( order_entries )
            function += ' ) )'

        if self.func_name == 'power' : 

            order_entries = []
            for i in range( 0, self.func_norders ) :
                order_entries.append( '@%d*TMath::Power( @0, @%d )'  %( i*2+1, i*2+2 ) )

            function = '+'.join( order_entries )

        mod_function = function.replace( '@0', 'x' )
        mod_function = '[0]*' + mod_function
        for i in range( 0, 9 ) :
            mod_function = mod_function.replace( '@%d' %i, '[%d]' %i )
        return mod_function

    @f_Obsolete
    def fit_histogram( self, workspace=None ) :
        self.run_fit()
        self.calculate_func_pdf()
        return self.get_results( workspace )

    def create_standard_ratio_canvas(self) :
        """ ported from SampleManager: setup canvas with ratio inset"""
        xsize = 800 
        ysize = 750
        self.curr_canvases['base'] = ROOT.TCanvas('basecan', 'basecan', xsize, ysize)

        self.curr_canvases['bottom'] = ROOT.TPad('bottompad', 'bottompad', 0.01, 0.01, 0.99, 0.34)
        self.curr_canvases['top'] = ROOT.TPad('toppad', 'toppad', 0.01, 0.34, 0.99, 0.99)
        self.curr_canvases['top'].SetTopMargin(0.08)
        # so that the ratio plot touches main plot, ie no gaps in btwn
        self.curr_canvases['top'].SetBottomMargin(0.0)
        self.curr_canvases['top'].SetLeftMargin(0.1)
        self.curr_canvases['top'].SetRightMargin(0.05)
        self.curr_canvases['bottom'].SetTopMargin(0.00) #no gaps
        self.curr_canvases['bottom'].SetBottomMargin(0.3)
        self.curr_canvases['bottom'].SetLeftMargin(0.1)
        self.curr_canvases['bottom'].SetRightMargin(0.05)
        self.curr_canvases['base'].cd()
        self.curr_canvases['bottom'].Draw()
        self.curr_canvases['top'].Draw()

    def fitrangehelper(self,fitrange):
        """ Convert fit range to tuples """
        if not fitrange: return
        if isinstance(fitrange,tuple) and len(fitrange)==2:
                return fitrange
        #if isinstance(fitrange,ROOT.RooRealVar):
        #        print "updating fit range"
        #        self.xvarfit = fitrange
        #xmin = self.xvarfit.getMin()
        #xmax = self.xvarfit.getMax()
        xmin = fitrange.getMin()
        xmax = fitrange.getMax()
        return (xmin,xmax)

    def setup_fit(self ,fitrange = None, dofit=False, **setupargs) : 
        """ fitter using RooFit package """
        # get fit function string if it is simple
        func_str = self.get_fit_function() 
        fitrange = self.fitrangehelper(fitrange)
    
        #arg_list = ROOT.RooArgList()
        #ROOT.SetOwnership(arg_list, False)
        
        #arg_list.add( self.xvardata )
        # adding extra variables for dijet functions
        # arglist is passed to take new variables
        #self.add_vars( arg_list)

        if self.func_name == 'bwxcb' : self.init_bwxcb()
        elif self.func_name == 'cb'  : self.init_cb()
        elif self.func_name == 'dcb' : self.init_dcb(reparam=False, **setupargs)
        elif self.func_name == 'dcbp' : self.init_dcb(**setupargs)
        elif self.func_name == 'dcbexpo' : self.init_dcbexpo(**setupargs)
        else :
             # generic function maker
             self.func_pdf = self.MakeROOTObj('RooGenericPdf', '%s_%s'
                    %(self.func_name, self.label), self.func_name, func_str, arg_list)

        # do fit if requested
        if dofit:  
            self.run_fit(fitrange)
        return self.func_pdf

    def run_fit(self, fitrange = None):
        """ run RooFit Fitter """
        self.fitrange = self.fitrangehelper(fitrange)
        self.xvardata.setRange("runfit",*self.fitrange)
        self.fitresult = self.func_pdf.fitTo( self.datahist,
                                         RooFit.Range("runfit"),
                                         RooFit.SumCoefRange("runfit"),
                                         #ROOT.RooFit.Range(*self.fitrange),
                                         #ROOT.RooFit.Extended(),
                                         RooFit.SumW2Error(True), 
                                         ROOT.RooCmdArg( 'Strategy', 3 ) ,
                                         RooFit.Save(ROOT.kTRUE))
        return

    def run_fit_chi2(self, fitrange = None):
        """ run RooFit Chi2 Fitter
            - To be used with weighted data
        """
        self.fitrange = self.fitrangehelper(fitrange)
#        self.option = ROOT.RooLinkedList()
 #       self.option.Add(ROOT.RooFit.Range(*self.fitrange))
        self.fitresult = self.func_pdf.chi2FitTo( self.datahist,ROOT.RooLinkedList())
        #ROOT.RooFit.Range(*self.fitrange))
                     #ROOT.RooFit.SumW2Error(True), ROOT.RooCmdArg( 'Strategy', 1) ,
                     #ROOT.RooFit.SumW2Error(True), ROOT.RooCmdArg( 'Strategy', 3 ) ,
                     #ROOT.RooFit.Save(ROOT.kTRUE))
        return

    def run_simplefit(self):
        """ simple fitter """
        fitrange = self.fitrangehelper(fitrange)
        func_str = self.get_fit_function() 

        self.func = self.MakeROOTObj( 'TF1', 'tf1_%s' %self.label, func_str, xmin, xmax )

        if self.func_name == 'dijet' : 

            self.func.SetParameter( 0, 0.0000001)
            param = 1
            for i in range( 1, self.func_norders+1 ) :
                this_def = self.get_vals( self.func_name, i )
                self.func.SetParameter( param, this_def[0] )
                param += 1
        self.hist.Fit( self.func, 'R' )
        return self.func


    def get_parameters(self):
        return self.defs.items()
    
    def get_parameter_values(self):
        plist = {}
        for name, parm in self.get_parameters():
            plist[name] = ufloat(parm.getVal(),parm.getError())
        return plist


    @f_Obsolete
    def runfit_oldcomments(self):
            pass
            #nll = self.func_pdf.createNLL(self.datahist) ;
            #m = ROOT.RooMinimizer(nll) 
            #m.setStrategy(2)
            #raw_input('cont2')
            ##// Activate verbose logging of MINUIT parameter space stepping
            #m.setVerbose(ROOT.kTRUE) 
            #raw_input('cont3')
            ##// Call MIGRAD to minimize the likelihood
            #m.simplex() 
            #raw_input('cont3.5')
            #m.migrad() 
            #raw_input('cont4')
            ##// Print values of all parameters, that reflect values (and error estimates)
            ##// that are back propagated from MINUIT
            #self.func_pdf.getParameters(ROOT.RooArgSet(self.xvar)).Print("s") 
            #raw_input('cont5')
            ##// Disable verbose logging
            #m.setVerbose(ROOT.kFALSE) 
            #raw_input('cont6')
            ##// Run HESSE to calculate errors from d2L/dp2
            #m.hesse() 
            #raw_input('cont7')
            ##// Print value (and error) of sigma_g2 parameter, that reflects
            ##// value and error back propagated from MINUIT
            ##cb_cut.Print()
            ##raw_input('cont10')
            #self.cb_sigma.Print()
            ##raw_input('cont11')
            ##cb_power.Print()
            ##raw_input('cont12')
            ##cb_m0.Print()
            #raw_input('cont13')
            ##// Run MINOS on sigma_g2 parameter only
            #m.minos(ROOT.RooArgSet(self.cb_sigma)) 
            #raw_input('cont14')
            #self.cb_sigma.Print()
            #raw_input('cont15')
            ##// Print value (and error) of sigma_g2 parameter, that reflects
            ##// value and error back propagated from MINUIT
            ##sigma_g2.Print() ;

    @f_Obsolete
    def calculate_func_pdf( self ) :

        if self.func_pdf is not None : 
            print 'The PDF Function already exists.  It will be overwritten'
        
        if self.func_name == 'dijet' :
            for i in range( 1, self.func_norders+1 ) :
                fitted_result = self.func.GetParameter(i)
                fitted_error = self.func.GetParError(i)

                self.defs['dijet'][i] = ( fitted_result, fitted_result - fitted_error, fitted_result + fitted_error )

            arg_list = self.MakeROOTObj( 'RooArgList' )
            arg_list.add( self.xvar )
            self.add_vars( arg_list )

            for i in range( 1, arg_list.getSize() )  :
                fitted_error = self.func.GetParError(i)
                arg_list[i].setError( fitted_error )

            func_str = self.get_fit_function( forceUseRooFit=True) 

            self.func_pdf = self.MakeROOTObj( 'RooGenericPdf', '%s_%s' %( self.func_name, self.label), self.func_name, func_str, arg_list)


    def save_fit( self, sampMan=None, workspace = None, logy=False, stats_pos='right' ) :

        if sampMan is not None :

            for name, param in self.fit_params.iteritems() : 
                param.SetName( name) 

            can = ROOT.TCanvas( str(uuid.uuid4()), '' )
            frame = self.xvar.frame() 
            self.datahist.plotOn(frame)
            self.func_pdf.plotOn( frame )
            if stats_pos == 'left' : 
                self.func_pdf.paramOn(frame, ROOT.RooFit.ShowConstants(True), ROOT.RooFit.Layout(0.1,0.5,0.9), ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(3)));
            if stats_pos == 'right' :
                self.func_pdf.paramOn(frame, ROOT.RooFit.ShowConstants(True), ROOT.RooFit.Layout(0.5,0.9,0.9), ROOT.RooFit.Format("NEU",ROOT.RooFit.AutoPrecision(2)));
            frame.Draw()
            if logy :
                ymax = frame.GetMaximum()
                frame.SetMinimum( 0.001 )
                frame.SetMaximum( ymax*10 )
                can.SetLogy()
                can.SetLogx()

            sampMan.outputs[self.label] = can

    def get_results( self, workspace = None) :


        results = {}
        for param in self.fit_params.values() :
            results[param.GetName()] = ufloat( param.getValV(), param.getErrorHi() )


        #power_res = ufloat( power.getValV(), power.getErrorHi() )
        #log_res   = ufloat( logcoef.getValV(), logcoef.getErrorHi())
        #int_res   = ufloat( integral, math.sqrt( integral ) )


        results['integral'] = self.Integral( )
        integral_var = ROOT.RooRealVar('dijet_%s_norm' %( self.label ), 'normalization', results['integral'].n )
        integral_var.setError( results['integral'].s )

        if workspace is not None :
            getattr( workspace , 'import' ) ( self.datahist )
            getattr( workspace , 'import' ) ( self.func_pdf )
            getattr( workspace , 'import' ) ( integral_var )


        return results


    def get_correlations(self):
        ROOT.gStyle.SetOptStat(0)
        if not self.fitresult: return
        self.curr_canvases["corr"] = ROOT.TCanvas("c_corr","correlation",400,400)
        self.curr_canvases["corr"].SetLeftMargin(0.15)
        self.curr_canvases["corr"].SetRightMargin(0.15)
        self.hcorr = self.fitresult.correlationHist()
        self.hcorr.Draw("COLZ")
        return self.curr_canvases["corr"]
        

    def draw(self,title=" ",yrange=None,subplot="",component=False,**kw):
        """ 
            component: True: take default component list from self.component
                       list: list of component names to be drawn
                       default is false, not drawing components
            subplot  : pull or residue
        """
        # Make canvas
        dology = kw.get("logy",0)
        if subplot:
            self.create_standard_ratio_canvas()
            self.curr_canvases["top"].SetLogy(dology)
            self.curr_canvases["top"].cd()
        elif not subplot:
            self.canvas = ROOT.TCanvas("cfit%s"%self.label,"Fitter",800,500) #FIXME: change canvas
            self.canvas.SetLogy(dology)
            self.canvas.cd()

        # make frame
        self.frame = self.xvardata.frame(RooFit.Title(title)) 
        # plot data histogram
        if component: 
            hlist = [h for h in self.datahistlist.values() if h!= self.datahist]
            for h in hlist: h.plotOn(self.frame,RooFit.DataError(ROOT.RooAbsData.SumW2),
                    RooFit.DrawOption("B"),RooFit.DataError(ROOT.RooAbsData.None),
                        RooFit.XErrorSize(0),RooFit.FillColor(ROOT.kBlue-10))
        self.datahist.plotOn(self.frame,RooFit.DataError(ROOT.RooAbsData.SumW2))

        # plotting fitted function
        if self.func_pdf:
            #plotparm   = [self.frame,RooFit.NormRange('myrange')]
            plotparm   = [self.frame,]
            if  self.pdfplotrange and self.fitrange:
                plotparm.append(RooFit.Range(*self.fitrange))
                plotparm.append(RooFit.NormRange(*self.fitrange))
            if component == True:
                component = self.components
            if isinstance(component,list):
                for i,comp in enumerate(map(RooFit.Components,component)):
                    self.func_pdf.plotOn(*(plotparm+FitManager.LineDefs[i+1]+[comp,]))
            self.func_pdf.plotOn(*(plotparm+FitManager.LineDefs[0]))
            pmlayout = kw.get("paramlayout",(0.65,0.9,0.8))
            ## toggle off parameters with None in layout
            if  pmlayout: 
                self.func_pdf.paramOn(self.frame,
                            RooFit.Layout(*pmlayout)).Draw()
        # set y axis range
        if dology:
           self.frame.SetMinimum(1e-5)
        if isinstance(yrange,tuple):
            if yrange[0]!=None: self.frame.SetMinimum(yrange[0])
            if yrange[1]!=None: self.frame.SetMaximum(yrange[1])
        self.frame.Draw()
        ROOT.gPad.SetTicks(1,1)
        # draw the plot status label
        self.draw_label(label_config={"dx":-0.05,"dy":-0.03})
        # draw pull or residual plot
        if subplot:
            self.curr_canvases["bottom"].cd()
            if subplot == "pull": self.makepull(self.xvardata)
            if subplot == "resd": self.makeresd(self.xvardata)
            self.subframe.Draw()
            #if subplot == "pull": self.subframe.GetYaxis().SetRangeUser(-6,6)
            self.ratio_formatting()
            ROOT.gPad.SetTicks(1,1)
            chi  =self.getchisquare()
            # print chi-sq
            self.latex = ROOT.TLatex()
            self.latex.SetTextSize(0.1) 
            self.latex.SetTextColor(ROOT.kRed) 
            self.latex.DrawLatexNDC(0.12,0.9,"#chi^{2}/ndof = %.4g" %chi)
            self.subframe.addObject(self.latex)

            # draw line 
            x=self.subframe.GetXaxis()
            xlow =  x.GetXmin()
            xhigh = x.GetXmax()
            self.line = ROOT.TLine(xlow,0,xhigh,0)
            self.line.SetLineStyle(3)
            self.line.SetLineWidth(2)
            self.line.Draw()
            return self.curr_canvases["base"]
        return self.canvas

    def getchisquare(self,dof = None):
        if dof is None:
            return self.frame.chiSquare(6) #FIXME: dof

    def draw_label(self,**kw):
        draw_config = DrawConfig("","","",**kw) 
        labels = draw_config.get_labels()
        for lab in labels :
            lab.Draw()
            self.curr_decorations.append( lab )

    def ratio_formatting(self):
        """ resize ratio label sizes as they are plotted in smaller pad """
        xAxs = self.subframe.GetXaxis()
        yAxs = self.subframe.GetYaxis()
        xAxs.SetTitleSize(0.085)
        xAxs.SetLabelSize(0.08)
        xAxs.SetTitleOffset(0.9)
        yAxs.SetTitleSize(0.08)
        yAxs.SetLabelSize(0.08)
        yAxs.SetTitleOffset(0.45)
        return

    def makepull(self,x):
        """ make pull histogram: helper function to FitManager.draw() """
        hpull = self.frame.pullHist()
        self.subframe = x.frame(RooFit.Title(" ")) 
        self.subframe.addPlotable(hpull,"P")
        self.subframe.SetYTitle("Pull")
        return

    def makeresd(self,x):
        """ make residual histogram: helper function to FitManager.draw() """
        hresd = self.frame.residHist()
        self.subframe = x.frame(RooFit.Title(" ")) 
        self.subframe.addPlotable(hresd,"P")
        self.subframe.SetYTitle("Residual")
        return

    @f_Obsolete
    def get_defaults( self, sample, var, ieta ) :
    
        self.set_vals('dijet', 1, ( -10.5, -20, 0 ) )
        self.set_vals('dijet', 2, (-2.03, -10, 0) )
        self.set_vals('dijet', 3, ( 0.0, -10, 10) )
        self.set_vals('dijet', 4, (0.0, -10 ,10 ) )
        self.set_vals('dijet', 5, (0.0, -10 ,10 ) )
    
        self.set_vals('power_coef', 1, ( 1000, 0, 10000000 ) )
        self.set_vals('power_coef', 2, (1000, 0, 100000000) )
        self.set_vals('power_coef', 3, ( 0.0, 0, 10) )
        self.set_vals('power_coef', 4, (0.0, 0 ,10 ) )
        self.set_vals('power_coef', 5, (0.0, 0 ,10 ) )
    
        self.set_vals('power_pow', 1, ( -9.9, -100, 100 ) )
        self.set_vals('power_pow', 2, (-0.85, -10, 10) )
        self.set_vals('power_pow', 3, ( 0.0, -10, 10) )
        self.set_vals('power_pow', 4, (0.0, -10 ,10 ) )
        self.set_vals('power_pow', 5, (0.0, -10 ,10 ) )
    
        self.set_vals('atlas_num_power', 1, ( -9.9, -100, 100 ) )
        self.set_vals('atlas_den_power', 1, ( -9.9, -100, 100 ) )
        self.set_vals('atlas_den_logcoef', 1, ( -9.9, -100, 100 ) )
        self.set_vals('atlas_den_logcoef', 2, ( -9.9, -100, 100 ) )

        #self.set_vals('cb_sigma', 450, ( 26.41, 15., 35. ) )
        #self.set_vals('cb_power', 450, ( 2.15, 1., 4. ) )
        #self.set_vals('cb_mass',  450, ( -18.1, -30., -10. ) )
        
        #self.set_vals('cb_sigma', 450, ( 28., 1., 100. ) )
        #self.set_vals('cb_power', 450, ( 30., 0., 130. ) )
        #self.set_vals('cb_mass',  450, ( -18, -100, 0 ) )
    
        self.set_vals('cb_sigma', 500, ( 28., 1., 100. ) )
        self.set_vals('cb_power', 500, ( 2.15, 0., 10. ) )
        self.set_vals('cb_mass',  500, ( -18, -100, 0 ) )
        
        self.set_vals('cb_sigma', 450, ( 28., 1., 100. ) )
        self.set_vals('cb_power', 450, ( 2.15, 0., 10. ) )
        self.set_vals('cb_mass',  450, ( -18, -100, 0 ) )
        
        self.set_vals('cb_sigma', 400, ( 28., 1., 100. ) )
        self.set_vals('cb_power', 400, ( 2.15, 0., 10. ) )
        self.set_vals('cb_mass',  400, ( -18, -100, 0 ) )

        self.set_vals('cb_sigma', 90, ( 28., 1., 100. ) )
        self.set_vals('cb_power', 90, ( 2.15, 0., 10. ) )
        self.set_vals('cb_mass',  90, ( -18, -100, 0 ) )
        
        self.set_vals('cb_sigma', 0, ( 2.8, 1., 100. ) )
        self.set_vals('cb_power', 0, ( 2.15, 0., 10. ) )
        self.set_vals('cb_mass',  0, ( 90, 0, 100 ) )
        self.set_vals('cb_alpha', 0, ( -1, -10, 10,"" ) )

    


    def init_bwxcb(self):
        mass =  self.sample_params['mass']
        width = self.sample_params['width']

        bw_width = mass*width
        if bw_width < 2  :
            bw_width = 2

        bw_m = self.MakeROOTObj( 'RooRealVar', 'bw_mass_%s' %self.label, 'Resonance  Mass', mass, xmin, xmax, 'GeV' )
        bw_w = self.MakeROOTObj('RooRealVar', 'bw_width_%s' %self.label, 'Breit-Wigner width',bw_width, 0, 200,'GeV')
        bw_m.setConstant()
        #bw_w.setConstant()
        bw = self.MakeROOTObj('RooBreitWigner','bw_%s' %self.label, 'A Breit-Wigner Distribution', self.xvar, bw_m,bw_w)

        #------------------------------
        # crystal ball, has four parameters
        #------------------------------
        sigma_vals = self.defs['cb_sigma'][mass]
        power_vals = self.defs['cb_power'][mass]
        mass_vals  = self.defs['cb_mass'][mass]
        cb_cut   = self.MakeROOTObj('RooRealVar','cb_cut_%s' %self.label, 'Cut'  , 0.5, 0.5, 0.50 , '')
        cb_sigma = self.MakeROOTObj('RooRealVar','cb_sigma_%s' %self.label, 'Width', sigma_vals[0], sigma_vals[1], sigma_vals[2], 'GeV')
        cb_power = self.MakeROOTObj('RooRealVar','cb_power_%s' %self.label, 'Power', power_vals[0], power_vals[1], power_vals[2], '')
        cb_m0    = self.MakeROOTObj('RooRealVar','cb_mass_%s' %self.label, 'mass' , mass_vals[0], mass_vals[1], mass_vals[2],'GeV')

        #cb_cut.setConstant()
        #cb_power.setConstant()
        #cb_m0.setConstant()

        cb_cut.setError( 0.05 )
        cb_sigma.setError( 0.5 )
        cb_power.setError( 1. )
        cb_m0.setError( 1. )

        cb = self.MakeROOTObj('RooCBShape','cb_%s' %self.label, 'A  Crystal Ball Lineshape', self.xvar, cb_m0, cb_sigma,cb_cut,cb_power)

        self.func_pdf = self.MakeROOTObj('RooFFTConvPdf','sig_model_%s' %self.label,'Convolution', self.xvar, bw, cb)

        self.fit_params['bw_mass'] = bw_m
        self.fit_params['bw_width'] = bw_w
        #self.fit_params['bw'] = bw

        self.fit_params['cb_cut'] = cb_cut
        self.fit_params['cb_sigma'] = cb_sigma
        self.fit_params['cb_power'] = cb_power
        self.fit_params['cb_mass'] = cb_m0
        #self.fit_params['cb'] = cb

    def make_factory_string(self, classname, pdfname, valsar, norange=False):
        """ input string for RooFactory custom pdf
        
            valsar   : list of tuples of the format
                       (variable name, range)
            classname: pdf class name
            pdfname  : pdf instance name
        """
        factstr = []
        for v in valsar:
            name = v[0] 
            if not norange and v[1:]:
                vrange = ",".join(map(str,v[1:]))
                factstr.append("%s[%s]" %(name,vrange))
            else:
                factstr.append(name)

        argstr = ",".join(factstr)
        return "%s::%s(%s)" %(classname,pdfname,argstr)

    def make_expression_string(self, varname, varexp):
        """ make factory expression with name and formula of expression """
        print varexp
        varlist = re.split('[+\-*/() ]',varexp) # minus needs escape
        varlist = [v for v in varlist if v]
        varlist = ','.join(list(set(varlist)))
        return "expr::%s('%s',%s)" %(varname, varexp, varlist)


    def init_cb(self,icond = "cb"):
        #------------------------------
        # crystal ball, has four parameters
        #------------------------------
        valsar = FitManager.setuparray[icond] #[Setup("cb",100)] #[0][1] #FIXME
        cbvals = []
        for v in valsar:
            name = v[0] #+self.label
            vtmp = ROOT.RooRealVar(*v)
            self.defs[name] = vtmp # store in instance dict
            cbvals.append(vtmp)    # values given in 

        # assume that error is relatively stable
        self.defs["cb_mass"].setError( 1. )
        self.defs["cb_sigma"].setError( 0.5 )
        self.defs["cb_alpha"].setError( 0.05 )
        self.defs["cb_power"].setError( 1. )
                                                                                                    # x mean sigma alpha n
        cb = self.MakeROOTObj('RooCBShape','cb_%s' %self.label, 'A  Crystal Ball Lineshape', self.xvardata,*cbvals) 

        self.func_pdf = cb 

    def init_dcb(self,icond="dcbp", pdflabel = "", 
            makenewfactory = True, reparam = True, seterr = True):
        """ reparam: reparametrize by rooFactory expr
            makenewfactory: replace existing factory object
        """
        #------------------------------
        # double crystal ball reparametrized
        #------------------------------
        if pdflabel:
            pdflabel = "_"+pdflabel
            makenewfactory = False
        if isinstance(icond, list):
            valsar = icond
        elif isinstance(icond,str):
            valsar = FitManager.setuparray[icond] 
        factstr = self.make_factory_string("DoubleCB","doublecb"+pdflabel,valsar) 
        print factstr
        if not self.loaded:
                self.loaded = ROOT.gROOT.ProcessLineSync(".x DoubleCB.cxx+") 

        # make factory
        if makenewfactory or not self.wk: self.wk = ROOT.RooWorkspace("doublecb")
        # make ordinary double crystal ball
        self.wk.factory(factstr) 
        if reparam:
            # define reparametrization
            #self.wk.factory("expr::scaled_power1('dcb_power1/dcb_alpha1', 
            #dcb_power1, dcb_alpha1)")
            #self.wk.factory("expr::scaled_power2('dcb_power2/dcb_alpha2',
            #dcb_power2, dcb_alpha2)")
            exprstr1 = self.make_expression_string('scaled_power1',
                    'dcb_power1/dcb_alpha1')
            exprstr2 = self.make_expression_string('scaled_power2',
                    'dcb_power2/dcb_alpha2')
            self.wk.factory(exprstr1)
            self.wk.factory(exprstr2)
            self.wk.factory("DoubleCB::dcbp%s(x,dcb_mass,dcb_sigma,"
                               "dcb_alpha1,scaled_power1, dcb_alpha2, scaled_power2)" %pdflabel)
            self.func_pdf = self.wk.pdf("dcbp"+pdflabel)
        else:
            self.func_pdf = self.wk.pdf("doublecb"+pdflabel)
        self.retrieve_param(valsar)
        if seterr: self.dcbseterr()
        # retrieve parameters
        #for v in valsar:
        #    name = v[0] 
        #    if name!="x": self.defs[name] = self.wk.var(name)
        self.xvardata = self.wk.var("x") #update xvardata
        self.xvardata.SetTitle(xname)
        self.xvardata.setUnit("GeV")
        # explicitly require fit range in case of factory function
        self.pdfplotrange = True 

    def dcbseterr(self):
        # initial step sizes: set them to be reasonably small
        # ballparking only and need no change for the same category of fits
        pprint(self.defs)
        self.defs["dcb_mass"]  .setError( 0.1 )
        self.defs["dcb_sigma"] .setError( 0.1 )
        self.defs["dcb_alpha1"].setError( 0.01 )
        self.defs["dcb_power1"].setError( 0.1 )
        self.defs["dcb_alpha2"].setError( 0.01 )
        self.defs["dcb_power2"].setError( 0.1 )

    @f_Obsolete
    def init_dcb_old(self,icond="dcb"):
        #------------------------------
        # double crystal ball
        #------------------------------
        valsar = FitManager.setuparray[icond] ##FIXME
        factstr = self.make_factory_string("DoubleCB","pdf",valsar) 
        print factstr
        if not self.loaded:
                self.loaded = ROOT.gROOT.ProcessLineSync(".x DoubleCB.cxx+") 

        # make factory
        self.wk = ROOT.RooWorkspace("doublecb")
        self.wk.factory(factstr) 
        self.func_pdf = self.wk.pdf("pdf")
        self.retrieve_param(valsar)
        # retrieve parameters
        #for v in valsar:
        #    name = v[0] 
        #    if name!="x": self.defs[name] = self.wk.var(name)
        self.xvardata = self.wk.var("x") #update xvardata
        self.xvardata.SetTitle(xname)
        self.xvardata.setUnit("GeV")
        self.dcbseterr()
        # explicitly require fit range in case of factory function
        self.pdfplotrange = True 

    def retrieve_param(self,valsar):
        # retrieve parameters
        for v in valsar:
            name = v[0] 
            if isinstance(v, str): 
                name = v
            if name!="x":
                ## if variable is already registered
                #if name in self.defs: 
                #    print "already registered: ",name, self.defs[name], self.wk.var(name)
                self.defs[name] = self.wk.var(name)
        self.xvardata = self.wk.var("x") #update xv

    #def init_dcbexpo(self, iconddcb = "", sig = "dcbp", bkgd = "gaus", ext='simul', icondzg=None):
    def init_dcbexpo(self, sig = "dcbp", bkgd = "gaus", ext='simul', ic = {}):
        """ icond: initial condition for dcb """
        self.wk = ROOT.RooWorkspace("dcb_bkgd") 
        ## step 1: Double CB shape
        name1 = sig
        if   sig == "dcbp":
            valsar1 = FitManager.ParamDCB #NOTE contain no initial values
            self.init_dcb(icond = ic.get(sig), makenewfactory = False, seterr = False)
        elif sig == "gauzz":
            valsar1 = ic.get(sig,FitManager.setuparray[sig])
            factstr1 = self.make_factory_string("Gaussian"   ,sig,valsar1) 
            self.wk.factory(factstr1) 
        ## step 2: prepare background
        name2   = "bkgd"
        valsar2 = FitManager.setuparray[bkgd]
        if   "expo" in bkgd:
            factstr2 = self.make_factory_string("Exponential",name2,valsar2) 
        elif "gaus" in bkgd:
            factstr2 = self.make_factory_string("Gaussian"   ,name2,valsar2) 

        ## step 2a: prepare optional z gamma gaussian

        #if bkgd2 is not None:
        #    valsar = self.make_factory_string("Gaussian", name3, valsar3)
        ## step 3: extended pdf
        ## SUM for pdf instead of functions
        valsar3 = FitManager.setuparray[ext]
        factstr3 = "SUM::"+ext+"(%s[%g,%g,%g]"%valsar3[0]+ "*%s," %name1+\
                              "%s[%g,%g,%g]"%valsar3[1]+ "*%s)" %name2 
        ## fill component list
        self.components=[name2,'dcbp']

        #make factory
        self.wk.factory(factstr2) 
        self.wk.factory(factstr3) 
        ## print current factory setup
        self.wk.Print()
        self.retrieve_param(valsar1)
        self.retrieve_param(valsar2)
        self.retrieve_param(valsar3)
        self.func_pdf = self.wk.pdf(ext)
        print "DEFINED VARIABLES"
        pprint(self.defs.items())
        self.pdfplotrange = False


