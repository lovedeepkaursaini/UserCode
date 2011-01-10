import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.TagProbeFitTreeAnalyzer = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    #InputFileNames = cms.vstring("rfio:/castor/cern.ch/user/l/lovedeep/VJets_TnP/mrgd_136035to149181-149442_Nov4ReReco_Jan8_v2.root"),
    InputFileNames = cms.vstring("rfio:/castor/cern.ch/user/l/lovedeep/VJets_TnP/mrgd_DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Jan8_v3.root"),
    #InputFileNames = cms.vstring("/uscms_data/d2/kalanand/allTPtrees_35invpb.root"), 
    InputDirectoryName = cms.string("SCToGsf"),
    InputTreeName = cms.string("fitter_tree"),
    OutputFileName = cms.string("testEfficiency.root"),
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
        probe_sc_et = cms.vstring("Probe SC E_{T}", "0", "1000", "GeV/c"),
        probe_sc_abseta = cms.vstring("Probe |#eta|", "-2.5", "2.5", ""),
        probe_nJets = cms.vstring("Jet Mult.", "0", "7", ""),
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        #mcTrue = cms.vstring("MC true", "dummy[false=1,true=0]"),
        probe_passingALL = cms.vstring("isMuon", "dummy[pass=1,fail=0]")
    ),

    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
#    PDFs = cms.PSet(
#        gaussPlusLinear = cms.vstring(
#     "CBExGaussShape::signalRes(mass, mean[2.0946e-01, -5., 5.], sigma[8.5695e-04],alpha[3.8296e-04], n[6.7489e+00], sigma_2[2.5849e+00], frac[6.5704e-01])",  ### the signal function goes here     
#     "CBExGaussShape::signalRes(mass, mean[2.0946e-01, -5., 5.], sigma[8.5695e-04,-1.,1.],alpha[3.8296e-04,-1.,1.], n[6.7489e+00], sigma_2[2.5849e+00,-1.,1.], frac[6.5704e-01])",  ### the signal function goes here     
#    "ZGeneratorLineShape::signalPhy(mass)",
#    "RooExponential::backgroundPass(mass, cPass[-0.02, -5, 0])",
#    "RooExponential::backgroundFail(mass, cFail[-0.02, -5, 0])",
#    "FCONV::signal(mass, signalPhy, signalRes)",
#    "efficiency[0.9,0,1]",
#    "signalFractionInPassing[1.0]"
    #"Gaussian::signal(mass, mean[91.2, 89.0, 93.0], sigma[2.3, 0.5, 10.0])",
    #"RooExponential::backgroundPass(mass, cPass[0,-10,10])",
    #"RooExponential::backgroundFail(mass, cFail[0,-10,10])",
    #"efficiency[0.9,0,1]",
    #"signalFractionInPassing[0.9]"
#        ),
#    ),

    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
           "Voigtian::signal1(mass, mean1[90,80,100], width[2.495,-5,5],sigma1[2,-3,3])",
           "Voigtian::signal2(mass, mean2[90,80,100], width[2.495,-5,5],sigma2[4,-10,10])",
           "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
           "Exponential::backgroundPass(mass, lp[-0.1,-1,1])",
           "Exponential::backgroundFail(mass, lf[-0.1,-1,1])",
           "efficiency[0.9,0,1]",
           "signalFractionInPassing[0.9]"
        ),
    ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        #the name of the parameter set becomes the name of the directory
        nJets = cms.PSet(
            EfficiencyCategoryAndState = cms.vstring("probe_passingALL","pass"),
            UnbinnedVariables = cms.vstring("mass"),
            BinnedVariables = cms.PSet(
                 probe_nJets = cms.vdouble( 0,1,2 ),
            ),
            BinToPDFmap = cms.vstring("gaussPlusLinear")
 ),
#         eta_nJets = cms.PSet(
 #           EfficiencyCategoryAndState = cms.vstring("probe_passingALL","pass"),
  #          UnbinnedVariables = cms.vstring("mass"),
   #3         BinnedVariables = cms.PSet(
    #3             probe_sc_abseta = cms.vdouble( 0.,1.5,2.5 ),
     #            probe_nJets = cms.vdouble( 0,1,2,3,4 ),
      #      ),
       #     BinToPDFmap = cms.vstring("gaussPlusLinear")
       # ),
    )
)


stages = ["SCToGsf"]
#stages = ["SCToGsf","GsfToIso","IsoToId","IdToHLT"]
d={}
for stage in stages:
    setattr(process,stage,process.TagProbeFitTreeAnalyzer.clone())
    getattr(process,stage).InputDirectoryName=stage
    getattr(process,stage).OutputFileName=stage+".root"
seq = getattr(process, stages[0])
for i in range(1,len(stages)):
    seq = seq+getattr(process,stages[i])
process.myseq = cms.Sequence(seq)
process.fit = cms.Path(process.myseq)
