
#
#Anil Singh
#NCU, Taiwan
#

import FWCore.ParameterSet.Config as cms

## PARAMETERIZED FIT: NOT SEXY ENOUGH...TOO MUNDANE
bwConvCB = cms.PSet(
    fittingFunction = cms.vstring(
        "BreitWigner::physicsShape(mass,m0[91.1876,80,100],fwhm[2.495])",
        "CBShape::resolution(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "FCONV::signal(mass,physicsShape,resolution)",
        "Exponential::backgroundPass(mass, lp[0,-5,5])",
        "Exponential::backgroundFail(mass, lf[0,-5,5])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


## we need a separate PDF distribution for each bin.
## I create a mechanism below which may not be very
## elegant (because of an undelying C++ ugly class)
## but it works.

#    "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
 #   "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
#    "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
 #   "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",

##WARNING: if ur python is not good:  "DONT TOUCH THIS THING !!!"
surrogate = cms.PSet(
    fittingFunction   = cms.vstring(
    "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
    "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
    "ZGeneratorLineShape::signalPhy1(mass,0,30,2,2.5,1)",
    "ZGeneratorLineShape::signalPhy2(mass,20,30,2,2.5,0)",
    "Exponential::backgroundPass(mass, lp[0,-5,5])",
    "Exponential::backgroundFail(mass, lf[0,-5,5])",
    "FCONV::signalPass(mass, signalPhy1, signalRes1)",
    "FCONV::signalFail(mass, signalPhy2, signalRes2)",
    "efficiency[0.9,0,1]",
    "signalFractionInPassing[0.9]"
    ),
    )


def SetSignalShape(pdfPset,ptLow,ptHi,etaLow,etaHi): 
    """
    WARNING: THIS FUNCTION IS EXTREMELY DEPENDENT
    ON HOW THE SURROGATE WAS DEFINED. SO THE
    DEFINITION OF \"surrorgate\" ABOVE SHOULD NOT BE
    ALTERED. THE CODE:
    
    pdfPset = Surrogate PSet
    token2 = {1,2,3,4}, 1= (20<M<30), 2 = (30<M<40), 3=(40<M<50) so on...
    token3 = {0,1}, 0 =  Endcap, 1 = Barrel
    
    """
    
    ##make string for passing signal shape
    passing = "ZGeneratorLineShape::signalPhy1(mass,"+ptLow+","+ptHi+","+etaLow+","+etaHi+",1)"
    ##make string for failing signal shape
    failing = "ZGeneratorLineShape::signalPhy2(mass,"+ptLow+","+ptHi+","+etaLow+","+etaHi+",0)"
        
    ##find index for the passing signal shape
    indexPass = pdfPset.fittingFunction.index("ZGeneratorLineShape::signalPhy1(mass,0,30,2,2.5,1)")
    ##find index for the failing signal shape
    indexFail = pdfPset.fittingFunction.index("ZGeneratorLineShape::signalPhy2(mass,20,30,2,2.5,0)")
    
    ##replace old pass string with new one
    pdfPset.fittingFunction[indexPass] = passing
    ###replace old fail string with new one
    pdfPset.fittingFunction[indexFail] = failing
    print passing, '\n', failing

genConvCBEx = surrogate.clone()
SetSignalShape(genConvCBEx, "40", "50","0","0.08")

#rint genConvCBEx
    

igenConvCBEx = cms.PSet(
    fittingFunction   = cms.vstring(
        "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "ZGeneratorLineShape::signalPhy1(mass,1,1)",
        "ZGeneratorLineShape::signalPhy2(mass,0,1)",
        "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
        "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
        "FCONV::signalPass(mass, signalPhy1, signalRes1)",
        "FCONV::signalFail(mass, signalPhy2, signalRes2)",     
        "efficiency[0.9,0,1]", 
        "signalFractionInPassing[0.9]"
        ),
    )


genConvCBExHiPt= cms.PSet(
    fittingFunction   = cms.vstring(
        "CBShape::signalRes1(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "CBShape::signalRes2(mass,0,res[2.61,.001,10],alp[25,0.001,50],n[1.42,0.001,50])",
        "ZGeneratorLineShape::signalPhy1(mass,1,3)",
        "ZGeneratorLineShape::signalPhy2(mass,0,3)",
        "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
        "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
        "FCONV::signalPass(mass, signalPhy1, signalRes1)",
        "FCONV::signalFail(mass, signalPhy2, signalRes2)",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


recConvCBEx = cms.PSet(
    fittingFunction   = cms.vstring(
        "CBExGaussShape::signalResPass(mass, meanP[0.], sigmaP[8.5695e-04, 0., 3.],alphaP[3.8296e-04], nP[6.7489e+00], sigmaP_2[2.5849e+00], fracP[6.5704e-01])",
        "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -5., 5.], sigmaF[8.5695e-04, 0., 5.],alphaF[3.8296e-04], nF[6.7489e+00], sigmaF_2[2.5849e+00], fracF[6.5704e-01])",
        "ZGeneratorLineShape::passSignalPhy(mass, /afs/hep.wisc.edu/home/kaur/Effi/CMSSW_4_4_2/src/PhysicsTools/TagAndProbe/test/ScToPf.root, hMassPassHi)",
        "ZGeneratorLineShape::failSignalPhy(mass, ScToPf.root, hMassPassLow)",
        "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[90.0])",
        "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[90.0])",
        "FCONV::signalPass(mass, passSignalPhy, signalResPass)",
        "FCONV::signalFail(mass, failSignalPhy, signalResFail)",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


gausPlusLinear = cms.PSet(
    fittingFunction = cms.vstring(
        "Voigtian::signal1(mass, mean1[90,80,100], width[2.495,-5,5],sigma1[2,-3,3])",
        "Voigtian::signal2(mass, mean2[90,80,100], width[2.495,-5,5],sigma2[4,-10,10])",
        "SUM::signal(vFrac[0.8,0,1]*signal1, signal2)",
        "Exponential::backgroundPass(mass, lp[-0.1,-1,1])",
        "Exponential::backgroundFail(mass, lf[-0.1,-1,1])",
        "efficiency[0.9,0,1]",
        "signalFractionInPassing[0.9]"
        ),
    )


ptBin=cms.vdouble(20,30)
print ptBin[0]
