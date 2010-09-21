#!/usr/bin/python

#
#Lovedeep,Anil
#Panjab University
#Chandigarh
#

### We need to create a small private library of functions, before we start
##  doing the things


########################################################################
def genCode(allnames):
    '''
    This generate the additional python commands to
    run the standard TagProbe cfg files over the
    HLT paths of our choice.
    '''
    hltPath='''
\nprocess.hltPath_sequence = cms.Sequence('''
    treePath='''
\nprocess.treePath_sequence= cms.Sequence('''
    s=''
    codelines = []
    for trig in allnames:
        line=trig.replace('_','')
        path="Path"+line
        tree="Tree"+line
        l4='''
process.''' + path +''' =process.PassingHLT.clone()
process.''' + path +'''.hltTag =  cms.untracked.InputTag("'''+trig+'''", "", "HLT")

process.''' + tree +''' =process.IdToHLT.clone()
process.''' + tree + '''.flags = cms.PSet( probe_passing = cms.InputTag("'''+path+'''"))
        '''
        codelines.append(l4)
        hltPath += ' process.'+path+' +'
        treePath += ' process.'+tree+' +'
    hltPath=hltPath[:-1]
    hltPath=hltPath+")"
    treePath=treePath[:-1]
    treePath=treePath+")\n"
    for item in codelines:
        s=s+item
    return s+hltPath+treePath
########################################################################




#########################################################################
def genFitCode(allnames):
    '''
    This generate code for performing fit, to determine
    the efficiency for trigger paths of our choice.
    '''
    fitPath='''
\nprocess.fit_hlt = cms.Sequence('''
    s=''
    codelines = []
    for trig in allnames:
        line=trig.replace('_','')
        name="Fit_"+line
        tree="Tree"+line
        l4='''
process.''' + name + ''' =process.TagProbeFitTreeAnalyzer.clone()
process.''' + name + '''.InputDirectoryName = cms.string("'''+tree+'''")
process.''' + name + '''.OutputFileName =  cms.string("'''+"testEff_"+trig+".root"+'''")
        '''
        codelines.append(l4)
        fitPath += ' process.'+name+' +'
    fitPath=fitPath[:-1]
    fitPath=fitPath+")\n"
    for item in codelines:
        s=s+item
    return s+fitPath
#########################################################################




#########################################################################################
# Create a list of trig names from command arguments, or read them from the text file
########################################################################################

import sys

help = '''
    = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    Usage:  ./TnPJob.py -names HLT_NAME1 HLT_NAME1 HLT_NAME1 
    Usage:  ./TnPJob.py -file trigNames.txt 
    = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    '''


trigNames = []
arg = sys.argv
if not (len(arg)>1):
    print help
    sys.exit()
if arg[1]=="-names":
    trigNames=trigNames+arg[2:]
    
if arg[1]=="-file":
    filename =open(arg[2],'r')
    trigNames=filename.readlines()

if not (arg[1]=="-names" or arg[1]=="-file"):
    print help
    sys.exit()

## This bussiness is necessary for cleaning off the trigger names
    
dummy = []
for trigName in trigNames:
    trigName=trigName.replace(' ','')
    if(trigName[-1])=='\n':
        trigName=trigName[:-1]
        dummy.append(trigName)
    else:
        dummy.append(trigName)


## cleaned strings (free of '\n' and  leading-trailing whitespaces) 
trigNames = dummy

################################################################
# Generate the tree maker Code
##################################################################

standardFile = ''
inputFile = open('Electron_TagProbeTreeProducer_cfg.py','r')
for line in inputFile:
    standardFile=standardFile+line

pathString = standardFile[standardFile.find("process.tagAndProbe"):]
standardFile = standardFile.replace(pathString,'')
surgicalString = pathString.replace("process.tree_sequence","process.tree_sequence + process.hltPath_sequence +process.treePath_sequence")

code= standardFile+genCode(trigNames)
code = code+'\n'+surgicalString
theOutputFile=open('NewElectron_TagProbeTreeProducer_cfg.py','w')
theOutputFile.write(code)






################################################################
# Generate the tree fitter Code
##################################################################

fitterFile = ''
inputFile = open('mytestTagProbeFitTreeAnalyzer_Zee.py','r')
for line in inputFile:
    fitterFile=fitterFile+line

pathString = fitterFile[fitterFile.find("process.fit"):]
fitterFile = fitterFile.replace(pathString,'')
surgicalFitString = pathString.replace("process.TagProbeFitTreeAnalyzer","process.TagProbeFitTreeAnalyzer + process.fit_hlt")
Fitcode= fitterFile+genFitCode(trigNames)
Fitcode = Fitcode+'\n'+surgicalFitString
OutputFile=open('NewtestTagProbeFitTreeAnalyzer_Zee.py','w')
OutputFile.write(Fitcode)


