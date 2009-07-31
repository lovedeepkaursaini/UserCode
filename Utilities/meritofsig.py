#!/usr/bin/env python
import math
from ROOT import *

def Effi(file, fout, histname, xmin, xmax):
   integral=0
   file1=TFile(file)
   h1=file1.Get(histname)
   integralfull1=h1.Integral()
   integral1=GetInt(h1,xmin,xmax,integral)
   eff=integral1/integralfull1
   fout.write(str(xmax)+"   "+str(eff)+'\n')

def Intg(Sigfile,Bakfiles,fout,histname,xmin,xmax):
   integral=0
   file1=TFile(Sigfile)
   h1=file1.Get(histname)
   integralfull1=h1.Integral()
   integral1=GetInt(h1,xmin,xmax,integral)

   inte2=float()
   full2=float()
   for bak in Bakfiles:
      file2=TFile(bak)
      h2=file2.Get(histname)
      integralfull2=h2.Integral()
      integral2=GetInt(h2,xmin,xmax,integral)
      inte2 = inte2+integral2
      full2 = full2+integralfull2

   ms=integral1/math.sqrt(integral1+inte2)
   print xmin,xmax,ms
   Seff=integral1/integralfull1
   Beff=inte2/full2
   fout.write(str(xmax)+"   "+str(ms)+'\n')

def GetInt(h1,xmin,xmax,integral):
   axis=h1.GetXaxis()
   bmin=axis.FindBin(xmin)
   bmax= axis.FindBin(xmax)
   integral=h1.Integral(bmin,bmax)
   integral -=h1.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin))/axis.GetBinWidth(bmin)
   integral -=h1.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax)/axis.GetBinWidth(bmax)
   return integral


# for isoComb
#for i in range(1,50):
 #   Intg(0,0.04*i)  

# for Sigma eta eta
Bakfiles=[]
Bakfiles.append("Bak_QCD.root")
Bakfiles.append("Bak_ttbar.root")
Bakfiles.append("Bak_WJets.root")
fout = open("Ms.txt",'w')
fouteff = open("Seff.txt",'w')

outfiles = []
outfiles.append("Bak_QCD.txt")
outfiles.append("Bak_ttbar.txt")
outfiles.append("Bak_WJets.txt")

for i in range(1,25):
   #Intg("Signal.root",Bakfiles,fout,"Sigee_zdg_pat_1",0,0.001*i)
   #Effi("Signal.root",fouteff,"Sigee_zdg_pat_1",0,0.001*i)

   #Intg("Signal.root",Bakfiles,fout,"isoComb_zdg_pat_1",0,0.04*i)
   #Effi("Signal.root",fouteff,"isoComb_zdg_pat_1",0,0.04*i)

   Intg("Signal.root",Bakfiles,fout,"dEta_zdg_pat_1",-0.0008*i,0.0008*i)
   Effi("Signal.root",fouteff,"dEta_zdg_pat_1",-0.0008*i,0.0008*i)

for bak in range (0,len(Bakfiles)):
    outf = open(outfiles[bak],'w')
    for i in range(1,25):
      # Effi(Bakfiles[bak],outf,"Sigee_zdg_pat_1",0,0.001*i)

       Effi(Bakfiles[bak],outf,"dEta_zdg_pat_1",-0.0008*i,0.0008*i)
       #Effi(Bakfiles[bak],outf,"isoComb_zdg_pat_1",0,0.04*i)




# for dEtaIn
#for i in range(1,50):
 #   Intg(-0.0008*i,0.0008*i)  
   
 
