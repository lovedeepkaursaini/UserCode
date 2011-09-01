# Auto generated configuration file
# using: 
# Revision: 1.284.2.4 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: MY/MENLOPS/sherpa_7TeV_ewk_Wlepton4jetsincl_cff.py -s GEN --conditions START3X_V25::All --datatier GEN-SIM-RAW --eventcontent RAWSIM --customise MY/MENLOPS/sherpa_custom_cff.py -n 100 --no_exec --fileout outputGEN.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)


# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.284.2.4 $'),
    annotation = cms.untracked.string('MY/MENLOPS/sherpa_7TeV_ewk_Wlepton4jetsincl_cff.py nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('outputGEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)
# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START3X_V25::All'

#process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")
#process.RandomNumberGeneratorService.generator.initialSeed = 12345XXXX

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
  #   moduleSeeds = cms.PSet(
   #     caloRecHits = cms.untracked.uint32(754321),
    #     VtxSmeared = cms.untracked.uint32(223458),
     #    muonCSCDigis = cms.untracked.uint32(525432),
      #   muonDTDigis = cms.untracked.uint32(67673876),
       #  famosSimHits = cms.untracked.uint32(235791312),
        # MuonSimHits = cms.untracked.uint32(834032),
#         famosPileUp = cms.untracked.uint32(918273),
 #        muonRPCDigis = cms.untracked.uint32(524964),
  #       siTrackerGaussianSmearingRecHits = cms.untracked.uint32(34680),
   #      generator= cms.untracked.uint32(123445)
    # ),
    generator = cms.PSet(
    initialSeed = cms.untracked.uint32(123445),
    engineName = cms.untracked.string('TRandom3')
  ),
    VtxSmeared = cms.PSet(
    initialSeed = cms.untracked.uint32(123445),
    engineName = cms.untracked.string('TRandom3')
  ),

     sourceSeed = cms.untracked.uint32(1234)
)

#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#sourceSeed = cms.untracked.uint32(1234)
#generator1 = cms.PSet(
#    initialSeed = cms.untracked.uint32(123445),
#    engineName = cms.untracked.string('TRandom3')
#  )
#)


process.generator = cms.EDFilter("SherpaGeneratorFilter",
    resultDir = cms.untracked.string('Result'),
    SherpaParameters = cms.PSet(
        parameterSets = cms.vstring('Run'),
        Run = cms.vstring('(run){', 
            ' EVENTS = 1000;', 
            ' EVENT_MODE = HepMC;', 
            ' PH_CSS_KT_RESOLUTION=8;', 
            ' # avoid comix re-init after runcard modification', 
            ' WRITE_MAPPING_FILE 3;', 
            '}(run)', 
            '(beam){', 
            ' BEAM_1 = 2212; BEAM_ENERGY_1 = 3500.;', 
            ' BEAM_2 = 2212; BEAM_ENERGY_2 = 3500.;', 
            '}(beam)', 
            '(processes){', 
            ' Process 93 93 -> 90 91 93{4};', 
            ' NLO_QCD_Part BVIRS {2};', 
            ' Loop_Generator Internal;', 
            ' ME_Generator Amegic {2,3};', 
            ' Order_EW 2;', 
            ' CKKW sqr(20./E_CMS);', 
            ' Enhance_Function PPerp(p[4]) {3};', 
            ' Enhance_Function PPerp(p[4])+PPerp(p[5]) {4};', 
            ' Enhance_Function PPerp(p[4])+PPerp(p[5])+PPerp(p[6]) {5};', 
            ' Enhance_Function PPerp(p[4])+PPerp(p[5])+PPerp(p[6])+PPerp(p[7]) {6};', 
            ' Integration_Error 0.02 {5,6};', 
            ' End process;', 
            '}(processes)', 
            '(selector){', 
            ' Mass  90 91 2. E_CMS;', 
            '}(selector)', 
            '(me){', 
            ' ME_SIGNAL_GENERATOR = Comix Amegic Internal;', 
            ' EVENT_GENERATION_MODE = Weighted;', 
            ' NLO_Mode 2;', 
            '}(me)', 
            '(mi){', 
            ' MI_HANDLER = Amisic  # None or Amisic', 
            '}(mi)')
    ),
    libDir = cms.untracked.string('/afs/hep.wisc.edu/user/kaur/CMSSW_4_1_7/src/MY/MENLOPS/normal/SherpaRun'),
    filterEfficiency = cms.untracked.double(1.0),
    crossSection = cms.untracked.double(-1),
    maxEventsToPrint = cms.untracked.int32(0)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.ProductionFilterSequence * getattr(process,path)._seq 

# customisation of the process


# Automatic addition of the customisation function from MY.MENLOPS.sherpa_custom_cff

def customise(process):

	process.genParticles.abortOnUnknownPDGCode = False

	return(process)


process = customise(process)


# End of customisation functions
