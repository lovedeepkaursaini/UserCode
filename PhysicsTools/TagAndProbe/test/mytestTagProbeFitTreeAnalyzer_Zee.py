import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring("rfio:/castor/cern.ch/user/l/lovedeep/VJets_TnP/mrgd_148952-149063_149181-149442_Nov3.root"),
    InputDirectoryName = cms.string("IdToHLT"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("IsoToId_EffiEle15_ID95.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_gsfEle_pt = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        probe_sc_et = cms.vstring("Probe e_{T}", "0", "1000", "GeV"),
        probe_sc_eta = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        probe_sc_abseta = cms.vstring("Probe |#eta|", "0", "2.5", ""),
        probe_sc_phi = cms.vstring("Probe #phi", "-4.0", "4.0", ""),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        probe_passing = cms.vstring("Probe Passing", "dummy[pass=1,fail=0]")
#        passing = cms.vstring("passing", "dummy[pass=1,fail=0]")
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        gen_exp = cms.vstring(
     "CBExGaussShape::signalRes(mass, mean[1.99996e-01, -10., 10.], sigma[8.269e-04],alpha[3.8296e-04], n[6.7489e+00], sigma_2[2.5849e+00], frac[6.5704e-01])",      
    "ZGeneratorLineShape::signalPhy(mass)",
    "RooExponential::backgroundPass(mass, cPass[-0.025, -5, 0])",
    "RooExponential::backgroundFail(mass, cFail[-0.025, -5, 0])",
    "FCONV::signal(mass, signalPhy, signalRes)",
    "efficiency[0.9,0,1]",
    "signalFractionInPassing[0.9]"
        ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory

        eta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("probe_passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                probe_sc_eta = cms.vdouble(-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5)
            ),
            #first string is the default followed by binRegExp - PDFname pairs
            BinToPDFmap = cms.vstring("gen_exp")
        ),
        phi = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("probe_passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                probe_sc_phi = cms.vdouble(-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0)
            ),
            #first string is the default followed by binRegExp - PDFname pairs
            BinToPDFmap = cms.vstring("gen_exp")
        ),
        et = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("probe_passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                probe_sc_et = cms.vdouble(20,25,30,35,40,45,50,55,60,70,80,100,120)
            ),
            #first string is the default followed by binRegExp - PDFname pairs
            BinToPDFmap = cms.vstring("gen_exp")
        ),
        abseta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("probe_passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                probe_sc_abseta = cms.vdouble(0.0, 1.5,2.5)
            ),
            BinToPDFmap = cms.vstring("gen_exp")
        ),
        et_abseta = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("probe_passing","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                probe_sc_et = cms.vdouble(20,25,30,35,40,45,50,55,60,70,80,100,120),
                probe_sc_abseta = cms.vdouble(0.0, 1.5, 2.5)
            ),
            BinToPDFmap = cms.vstring("gen_exp")
        ),
    )
)

process.fit = cms.Path(process.TagProbeFitTreeAnalyzer)
