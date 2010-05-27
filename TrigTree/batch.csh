#!/bin/bash
cd /afs/cern.ch/user/l/lovedeep/scratch0/TRIG/CMSSW_3_6_1/src/UserCode/TrigTree/
# do root in batch mode and quit when finished
$ROOTSYS/bin/root -b <<EOF
{
      gROOT->LoadMacro("/afs/cern.ch/user/l/lovedeep/scratch0/TRIG/CMSSW_3_6_1/src/UserCode/TrigTree/trigEffi.C+");
      trigEffi t;
      t.Loop("mc");
}
EOF
