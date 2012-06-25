import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################


isMC = False
#InputFileName = "TnPTree_DoubleElecRun2011A.root"
InputFileName ="/hdfs/store/user/anil79/TnPTree_DYToEEFall11POWHEG_May30-newformcMakeTnPtree/1/newformcMakeTnPtree-FE34DF89-152A-E111-BE45-0015178C4D94.root"
#InputFileName = "/hdfs/store/user/anil79/mrgdTightIDTnPTree_DYJetsSummer12.root"
#InputFileName = "/hdfs/store/user/anil79/mrgdLooseIDTnPTree_DYJetsSummer12.root"
#InputFileName = "/hdfs/store/user/anil79/mrgdVetoIDTnPTree_DYJetsSummer12.root"
OutputFilePrefix = "/scratch/anil79/lvefficiency-MCDYJETS-MEDIUM"



################################################
HLTDef = "probe_passingHLT"
PDFName = "pdfSignalPlusBackground"

if isMC:
    InputFileName = "TP_Spring11_DYToEE.root"
    PDFName = ""
    OutputFilePrefix = "efficiency-Spring11-DYToEE-pt20-"
################################################

#specifies the binning of parameters
EfficiencyBins = cms.PSet(
    probe_sc_et = cms.vdouble( 10,15,20,30,40,50, 200 ),
    probe_sc_abseta = cms.vdouble( 0,0.8,1.4442,1.556,2.0,2.5 ),
#    run = cms.vdouble(160431, 167676, 173198, 176929, 178151, 180252),
 #   run = cms.vdouble(190000, 200000),
)
## for super clusters
EfficiencyBinsSC = cms.PSet(
    probe_et = cms.vdouble( 10,15,20,30,40,50, 200 ),
    probe_abseta = cms.vdouble( 0,0.8,1.4442,1.556,2.0,2.5 ),
#    probe_eta = cms.vdouble( -1.44, 1.44 ),
#    run = cms.vdouble(160431, 167676, 173198, 176929, 178151, 180252),
    #run = cms.vdouble(190000, 200000),
)

#### For data: except for HLT step
EfficiencyBinningSpecification = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(EfficiencyBins),
    #first string is the default followed by binRegExp - PDFname pairs
    BinToPDFmap = cms.vstring(PDFName)
)


#### For super clusters
EfficiencyBinningSpecificationSC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsSC),
    BinToPDFmap = cms.vstring(PDFName)
)
EfficiencyBinningSpecificationSCMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsSC,mcTrue = cms.vstring("true")),
    BinToPDFmap = cms.vstring()  
)


#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
    probe_sc_et = cms.vdouble( 10,15,20,30,40,50, 200 ),
    probe_sc_abseta = cms.vdouble( 0,0.8,1.4442,1.556,2.0,2.5 ),
#    run = cms.vdouble(160431, 167676, 173198, 176929, 178151, 180252),
#    run = cms.vdouble(190000, 200000),
    mcTrue = cms.vstring("true")
    ),
    BinToPDFmap = cms.vstring()  
)

#### For HLT step: just do cut & count
EfficiencyBinningSpecificationHLT = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBins),
    BinToPDFmap = cms.vstring()  
)

##########################################################################################
############################################################################################
if isMC:
    mcTruthModules = cms.PSet(
        MCtruth_WPMedium = cms.PSet(
        EfficiencyBinningSpecificationMC,   
        EfficiencyCategoryAndState = cms.vstring("probe_isWPMedium","pass"),
        ),
        )
else:
    mcTruthModules = cms.PSet()
##########################################################################################
##########################################################################################


        

############################################################################################
############################################################################################
####### GsfElectron->Id / selection efficiency 
############################################################################################
############################################################################################

process.GsfElectronToId = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(InputFileName),
    InputDirectoryName = cms.string("GsfElectronToId"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string(OutputFilePrefix+"GsfElectronToId.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_sc_et = cms.vstring("Probe E_{T}", "0", "200", "GeV/c"),
        probe_sc_abseta = cms.vstring("Probe |#eta|", "0.", "2.5", ""),
#                    run = cms.vstring("Run number", "160431", "180252", ""),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        probe_isWPMedium = cms.vstring("probe_isWPMedium", "dummy[pass=1,fail=0]"),
        ),
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring(
"BreitWigner::physicsShape(mass,m0[91.1876,80,100],fwhm[2.495])",  #replace with floating + constriant later?
            "CBShape::resolution(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
            "FCONV::signal(mass,physicsShape,resolution)",
#            "Exponential::backgroundPass(mass, lp[0,-5,5])",
#            "Exponential::backgroundFail(mass, lf[0,-5,5])",
    "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
    "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
            "efficiency[0.9,0,1]",
            "signalFractionInPassing[0.9]" 

##     "CBExGaussShape::signalRes(mass, mean[2.0946e-01], sigma[8.5695e-04],alpha[3.8296e-04], n[6.7489e+00], sigma_2[2.5849e+00], frac[6.5704e-01])",  
### the signal function goes here
#     "CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",  ### signal resolution for "pass" sample
 #    "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",  ### signal resolution for "fail" sample     
  #  "ZGeneratorLineShape::signalPhy(mass)", ### NLO line shape
   # "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
    #"RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
#    "FCONV::signalPass(mass, signalPhy, signalResPass)",
 #   "FCONV::signalFail(mass, signalPhy, signalResFail)",     
  #  "efficiency[0.9,0,1]",
#    "signalFractionInPassing[1.0]"     
   # "Gaussian::signal(mass, mean[91.2, 89.0, 93.0], sigma[2.3, 0.5, 10.0])",
#    "RooExponential::backgroundPass(mass, cPass[-0.02,-5,0])",
 #   "RooExponential::backgroundFail(mass, cFail[-0.02,-5,0])",
  #  "efficiency[0.9,0,1]",
#    "signalFractionInPassing[0.9]"
            #"BreitWigner::signal(mass,m0[91.1876,80,100],fwhm[2.495])",  #replace with floating + constriant later?
#"BreitWigner::physicsShape(mass,m0[91.1876,80,100],fwhm[2.495])",  #replace with floating + constriant later?
 #           "CBShape::resolution(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
  #          "FCONV::signal(mass,physicsShape,resolution)",
#            "Exponential::backgroundPass(mass, lp[0,-5,5])",
 #           "Exponential::backgroundFail(mass, lf[0,-5,5])",
  #          "efficiency[0.9,0,1]",
#            "signalFractionInPassing[0.9]" 

       ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
    mcTruthModules,
##     #the name of the parameter set becomes the name of the directory
    WPMedium = cms.PSet(
    EfficiencyBinningSpecification,   
    EfficiencyCategoryAndState = cms.vstring("probe_isWPMedium","pass"),
    ),
    )
)

process.GsfElectronToId.Variables.WeightVariable = cms.string("PUweight")
process.GsfElectronToId.Efficiencies.WPMedium.BinToPDFmap= cms.vstring()


############################################################################################
############################################################################################
####### SC->GsfElectron efficiency 
############################################################################################
############################################################################################
if isMC:
    SCmcTruthModules = cms.PSet(
        MCtruth_efficiency = cms.PSet(
        EfficiencyBinningSpecificationSCMC,
        EfficiencyCategoryAndState = cms.vstring( "probe_passingGsf", "pass" ),
        ),
    )
else:
    SCmcTruthModules = cms.PSet()    


process.SCToGsfElectron = process.GsfElectronToId.clone()
process.SCToGsfElectron.InputDirectoryName = cms.string("SuperClusterToGsfElectron")
process.SCToGsfElectron.OutputFileName = cms.string(OutputFilePrefix+"SCToGsfElectron.root")
process.SCToGsfElectron.Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_et = cms.vstring("Probe E_{T}", "0", "200", "GeV/c"),
        probe_abseta = cms.vstring("Probe |#eta|", "0.", "2.5", ""),#                        run = cms.vstring("Run number", "160431", "180252", ""),
    )
process.SCToGsfElectron.Categories = cms.PSet(
    mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),           
    probe_passingGsf = cms.vstring("probe_passingGsf", "dummy[pass=1,fail=0]"),                        
    )
process.SCToGsfElectron.Efficiencies = cms.PSet(
    SCmcTruthModules,
    efficiency = cms.PSet(
    EfficiencyBinningSpecificationSC,
    EfficiencyCategoryAndState = cms.vstring( "probe_passingGsf", "pass" ),
    ),
)
#fit nai ho reha
process.SCToGsfElectron.Efficiencies.efficiency.BinToPDFmap = cms.vstring()
process.SCToGsfElectron.Variables.WeightVariable = cms.string("PUweight")

############################################################################################
############################################################################################
####### HLT efficiency 
############################################################################################
############################################################################################


if isMC:
    HLTmcTruthModules = cms.PSet(
        MCtruth_efficiency = cms.PSet(
        EfficiencyBinningSpecificationMC,
        EfficiencyCategoryAndState = cms.vstring( HLTDef, "pass" ),
        ),    
    )
else:
    HLTmcTruthModules = cms.PSet()


EfficienciesPset = cms.PSet(
    HLTmcTruthModules,
    efficiency = cms.PSet(
    EfficiencyBinningSpecificationHLT,
    EfficiencyCategoryAndState = cms.vstring( HLTDef, "pass" ),
    ),
)

########
process.WPMediumToHLTEle = process.GsfElectronToId.clone()
process.WPMediumToHLTEle.InputDirectoryName = cms.string("WPMediumToHLTEle")
process.WPMediumToHLTEle.OutputFileName = cms.string(OutputFilePrefix+"WPMediumToHLTEle.root")
process.WPMediumToHLTEle.Categories = cms.PSet(
    mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),           
    probe_passingHLT = cms.vstring("probe_passingHLT", "dummy[pass=1,fail=0]"), 
    )
process.WPMediumToHLTEle.Efficiencies = EfficienciesPset
process.WPMediumToHLTEle.Efficiencies.efficiency.BinToPDFmap = cms.vstring()
########
process.WPMediumToHLTEle17 = process.WPMediumToHLTEle.clone()
process.WPMediumToHLTEle17.InputDirectoryName = cms.string("WPMediumToHLTEle17")
process.WPMediumToHLTEle17.OutputFileName = cms.string(OutputFilePrefix+"WPMediumToHLTEle17.root")
process.WPMediumToHLTEle17.Variables.WeightVariable = cms.string("PUweight")

########
process.WPMediumToHLTEle8NotEle17 = process.WPMediumToHLTEle.clone()
process.WPMediumToHLTEle8NotEle17.InputDirectoryName = cms.string("WPMediumToHLTEle8NotEle17")
process.WPMediumToHLTEle8NotEle17.OutputFileName = cms.string(OutputFilePrefix+"WPMediumToHLTEle8NotEle17.root")
process.WPMediumToHLTEle8NotEle17.Variables.WeightVariable = cms.string("PUweight")

process.fit = cms.Path(
    process.GsfElectronToId   
 #  process.SCToGsfElectron  
    # process.WPMediumToHLTEle17  
#     process.WPMediumToHLTEle8NotEle17 
    )
