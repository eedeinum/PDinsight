#!/usr/bin/python3

from math import log,pi,sqrt,sin
import sys
import operator
from copy import deepcopy
from collections import defaultdict
from statistics import mean
import random as rand

VERY_MUCH = sys.maxsize
VERSION = "2.0b"

try:
  from numpy import arange
except:
  def arange(xStart,xStop,xStep):
    aList = [xStart]
    while aList[-1] < xStop:
      aList.append(aList[-1]+xStep)
    return aList


## This script computes the requirements for obtaining a given target effective wall permeability (Peff) through diffusive symplastic transport. Depending on the computation modes used, either the density or the maximum particle size that can pass the model PD is given as an input parameter. Directly computing Peff from a combination of parameters is also possible. Multiple values of input parameters can be evaluated in a single run. See parameter file for explanation of the different modes.
## This script belongs to the manuscript "From plasmodesma geometry to effective symplasmic permeability through biophysical modelling" Deinum et al 2019, eLife;  doi: 10.7554/eLife.49000. When using values from this script in a publication, please cite Deinum et al 2019. 
## When using singleInterface or bootstrapInterface modes, please cite Deinum 2022: More Insights from Ultrastructural and Functional Plasmodesmata Data Using PDinsight, Methods Mol Biol . 2022:2457:443-456.  doi: 10.1007/978-1-0716-2132-5_30. 
## In this script, PDs are modeled as a set of concentric cylinders. The space available for symplastic transport, the cytoplasmic sleeve, is confined between desmotubule (DT, radius Rdt) and outer radius (Rn, Rn2, neck radius and/or Rc, central radius). 
### If Rn == Rc, the parameter Lneck (length of the neck region) does not have any effect.
### We found that Fih, a correction factor for the fact that the cell wall is only permeable at discrete spots (the PDs), does not depend on Lcell, the cell length. This parameter, therefore, also has no effect.

## usage: [python] ./PDinsight.py [parameterfile] 
## if no parameter file is given, "parameters.txt" from the current directory is used; program crashes if parameter file cannot be found
## after reading, all parameters are printed to "allParameters.txt"  or fileTag.pars, if fileTag is set

class Parameters():
  versionNumber=VERSION
  parExtraInfo = defaultdict(lambda:None)
  info = { "grid": 'Available options: "triangular", "square", "hexagonal" or "hex" and "random"', "gridList": 'Available options: "triangular", "square", "hexagonal" or "hex" and "random"', "pitList": "Available options: 1, 2, 3, 4, 5, 6, 7, 12, 19", "twinningList": "Available options: 1, 2, 3, 4, 5, 6, 7, 12, 19", "xMax" : "Maximum particle size: xMax = (Rn - Rdt)/ 2 (unobstructed sleeve model). In other words: Rn = Rdt + 2*xMax.", "xMaxList" : "Maximum particle size: xMax = (Rn - Rdt)/ 2 (unobstructed sleeve model). In other words: Rn = Rdt + 2*xMax.", "Rc": "If Rc < Rn, a straight channel with radius Rn is assumed. Use 0 to only study straight channels." , "RcList": "If Rc < Rn, a straight channel with radius Rn is assumed. Use 0 to only study straight channels.", "DTformula": "If true, DT= DTintercept + DTvar1*DTcoef1 + ...", "DTcoefList": "model: DT= DTintercept + DTvar1*DTcoef1 + ...", "DTvarList": "model: DT= DTintercept + DTvar1*DTcoef1 + ...", "DTintercept": "model: DT= DTintercept + DTvar1*DTcoef1 + ...", "randomSeed": "Positive integers are used as seed. For negative numbers, the random generator is initialized without seed.", "noiseModelList": "Available options: AU (additive, uniform), AG = AN (additive, normally distributed), M (multiplicative; lognormal distributed)" ,"boostrapSampleSize": "Determined automatically if 0. Normally, there is no need to set this by hand."}
  for i in info.keys():
    parExtraInfo[i] = info[i]
  parDescriptions = {"RnList": "Neck radius (RnList; nm)", "Rn2List": "2nd neck radius (Rn2List; nm)","Rn": "Neck radius (Rn; nm)", "Rdt": "Desmotubule radius (Rdt; nm)", "RdtList": "Desmotubule radius (RdtList; nm)", "Rc": "Central region radius (Rc; nm)", "RcList": "Central region radius (RcList; nm)" ,"Lpd": "PD length (Lpd; nm)" ,"LpdList": "PD length (LpdList; nm)" ,"Lneck": "Neck length (Lneck; nm)" ,"LneckList": "Neck length (LneckList; nm)" ,"Lneck2List": "2nd neck length (Lneck2List; nm)" ,"xMax": "Max particle size (xMax; nm)" ,"xMaxList": "Max particle size (xMaxList; nm)" ,"grid": "PD grid type (grid)" ,"gridList": "PD grid type (gridList)" , "densList": "PD density (densList; PDs/um2)" , "diff": "Particle diffusion constant (diff; um2/s)" , "dInc": "PD density increment (dInc; PDs/um2)" , "x": "Particle radius (x; nm)" , "xStart": "Starting max particle radius Rn-dens graph (xStart; nm)" , "xStep": "Max particle radius increment Rn-dens graph (xStep; nm)" , "xInc": "Max particle radius increment computeAperture (xInc; nm)" , "Lcell": "Cell length = irrelevant (Lcell; um)" , "dPit": "PD-PD distance inside pit (dPit; nm)" , "dPitList": "PD-PD distance inside pit (dPitList; nm)" , "pitDens": "Pit density (pitDens; PDpits/um2)" , "pitList": "PDs per pit (pitList)" , "twinningList": "PDs per pit (twinningList)" , "x": "Particle radius (x; nm)" , "peList": "Target effective permeability (peList; um/s)" , "fileTag": "output file tag (fileTag)" , "sensIncList": "Stepsize for numerical differentiation (a.u.)", "DTvarList": "Model variables in DTmodel (DTvarList)","DTcoefList": "Coefficients for DTmodel variables (DTcoefList)","DTintercept": "Intercept of DTmodel (DTintercept)", "bootstrapSamples": "resampling n (bootstrapSamples)", "randomSeed": "random seed for bootstrapping", "noiseVarList":"Variables with noise (noiseVarList)","noiseMagnitudeList":"noiseMagnitudeList","noiseModelList":"noiseModelList (AU, AN, M)","noData":"token for missing data (noData)","bootstrapSingleSampleSize": "number of PDs per resampled set"}
  optionDescriptions = {"compSubNano": "Also compute values for sub-nano channel model (ignored for asymmetricPDs)","computeTwinning": "compute impact on Peff due to twinning, starting from isolated PDs", "computeClusterIncrease": "consider (repeated) twinning from starting densities in densList", "printRn": "print calculated Rn in output", "doNotCombine": "no: compute for all combinations of parameters in lists; yes: only compute for same position in list","asymmetricPDs": "Allow different neck length/radius on either side", "DTformula": "yes: Use linear model for DT radius", "addNoise": "Add noise when resampling"} 
  requiredPars = {"sensitivityAnalysis": ["RnList", "RdtList", "RcList","LpdList","densList","Lneck","sensIncList"], "computeVals": ["LpdList","LneckList","Lneck2List","densList","RnList","Rn2List","RcList","RdtList","pitList","dPitList","xMaxList","computeClusterIncrease","doNotCombine","compSubNano","asymmetricPDs","grid"] , "computeUnitVals": ["LpdList","LneckList","Lneck2List","RnList","Rn2List","RcList","RdtList","doNotCombine","asymmetricPDs"] , "computeDens": ["xMaxList","peList","Lpd","Rdt","dInc","printRn","grid"], "computeAperture": ["peList","Rdt","densList","Lpd","grid","xInc","compSubNano"] , "computeRnDensityGraph": ["peList","Lpd","Rdt","xStart","xMax","xStep","grid","dInc"], "all": ["fileTag","x","diff","outputType"],\
          "singleInterface": ["LpdList","LneckList","Lneck2List", "densList","RnList","Rn2List","RcList","RdtList","pitList","dPitList","asymmetricPDs","grid","PDmodel","doNotCombine","DTformula","DTvarList","DTcoefList","DTintercept","noData"],
          "bootstrapInterface": ["LpdList","LneckList","Lneck2List", "densList","RnList","Rn2List","RcList","RdtList","pitList","dPitList","asymmetricPDs","grid","PDmodel","bootstrapSamples","randomSeed","doNotCombine","DTformula","DTvarList","DTcoefList","DTintercept","addNoise","noiseVarList","noiseMagnitudeList","noiseModelList","noData","bootstrapSingleSampleSize"] 
          } ## point of debugging: may be missing a value or two; TODO update for singleInterface & bootstrapInterface. ## maybe also add to singleInterface: ,"addNoise","noiseVarList","noiseMagnitudeList","noiseModelList"
  optionList = ["compSubNano", "computeTwinning", "computeClusterIncrease", "printRn", "doNotCombine","asymmetricPDs","DTformula","addNoise"]
  MCoptionList = ["outputType", "PDmodel"]
  simpleFloats = ["Rdt", "x", "Lpd", "dPit", "Lneck", "xStep", "xStart", "xMax", "xInc","DTintercept" ]
  factorFloats = {"diff":1e6,"Lcell":1e3, "dInc":1e-6, "pitDens":1e-6}
  listFloats = [ "peList", "xMaxList", "twinningList", "pitList", "dPitList", "LpdList", "RnList", "RcList", "RdtList", "sensIncList" ,"LneckList","Lneck2List","Rn2List","DTcoefList","noiseMagnitudeList"]
  listFactorFloats = {"densList":1e-6} 
  simpleInts = ["bootstrapSamples","randomSeed","noData","bootstrapSingleSampleSize"]
  simpleBools = ["compSubNano", "computeFih_subNano", "computeFih_pitField_xMax", "computeFih_pitField_dens", "computeTwinning", "computeClusterIncrease", "computeRnDensityGraph","computeDens","computeAperture","computeVals", "printRn","sensitivityAnalysis", "computeUnitVals","singleInterface","bootstrapInterface","doNotCombine","asymmetricPDs","bootstrapInterface","singleInterface","addNoise","DTformula"]
  simpleStrings = ["grid", "fileTag","dataSep"]
  listStrings = ["gridList","DTvarList","noiseModelList","noiseVarList"]
  def __init__ (self):
    ### Default parameters ###  
    ## switches for computing various different quantities
    self.compSubNano=False # switch: compare with sub-nano channel model
    self.computeFih_subNano=False # switch: compute Fih values for default and sub-nano channel model
    self.computeFih_pitField_xMax=False # switch: compute Fih for pit fields as a function of maximum particle radius (xMax / alpha_bar) 
    self.computeFih_pitField_dens=False # switch: compute Fih for pit fields as a function of PD density (rho)
    self.computeTwinning=False # switch: compute impact on Peff due to twinning, starting from isolated PDs that follow the distribution of parameter grid. 
    self.computeClusterIncrease=False # switch: consider (repeated) twinning from the starting densities in densList rather than homogeneous (high) densities. Useful for computing requirements for a sudden increase of Peff. 
    self.computeRnDensityGraph=False # switch: compute Rn, rho graphs that yield target Peff in peList. Rn is increased by steps of xStep; Rn runs from Rdt+xStart to Rdt + xMax; incompatible with compSubNano, because correction factors can't be used for large Rn
    self.computeDens=False # switch: calculate required densities for Peff targets in peList and alpha_bar values in xMaxList ; output is written to file "requiredDensityTable.tsv"
    self.computeAperture=False ## calculate required alpha_bar and Rn values for Peff targets in peList and densities in densList; if computeClusterIncrease == True: calculate required alpha_bar and Rn values for Peff targets in peList, starting densities in densList and the n-fold increase in clusters from twinningList ; output is written to file "requiredApertureTable.tsv"
    self.computeVals=False # calculate Peff values from all given combinations of variable parameters: densList x xMaxList [[ x RcList]] x ....
    self.computeUnitVals=False # calculate unit permeability (Punit) values: for a density of 1 PD/um2, not considering inhomogeneity factor Fih. Calculate Punit from all given combinations of variable parameters: RnList x RdtList x LpdList [x RcList (insofar Rn <= Rc)] ; currently incompatible with compSubNano

    self.doNotCombine=False # switch: affects computeUnitVals and computeVals. Instead of looping over all combinations, equal length lists (or single values) of the applicable variables are required, and values are computed per i-th set of parameters. 
    self.printRn=True  # Print Rn values to table. Under the unobstructed sleeve model: Rn = Rdt + 2 xMax (so, typically, Rn is redundant, but may be useful to print, e.g., for plotting). ; For computeVals if compSubNano is True, Rn values are always printed, as they are non-trivial for the sub-nano model.
    self.sensitivityAnalysis=False  # Perform sensitivity analysis: Numerically compute derivatives of Peff and other quantities (to Rn, Rc, Rdt, Lpd, Lneck, rho=dens ) around a starting point given by parameters.
    self.asymmetricPDs=False  # Allow for different Lneck and Rn on either side of the PD. Required inputs: Rn1List, Rn2List, Lneck1List,Lneck2List (all equal length). Currenly only works for computeVals and computeUnitVals with doNotCombine==True
    ## new options for singleInterface and bootstrapInterface
    self.singleInterface=False # compute Punit per PD + Peff based on average Punit for complete data; averages of parameters; 
    self.bootstrapInterface=False # compute Peff + confidence intervals based on resampled data
    self.bootstrapDiagnosticGraphs=False # Not implemented yet.
    self.DTformula=False # infer Rdt from formula (linear statistical model) ; currently only for singleInterface and bootstrapInterface
    self.addNoise=False # flexible mode for adding measurement noise when resampling during bootstrap

    self.fileTag = "" # if exists, file names are appended with _fileTag
    self.outputType = "tsv" # output type: tsv or csv
    self.colSep="\t"
    self.fileExt="tsv"
    self.typeDict={ "tsv":["\t","tsv"], "TSV":["\t","tsv"], "csv":[",","csv"], "CSV":[",","csv"]} 
    self.outputTypeOptions=["tsv","csv"]
    self.PDmodelOptions=["1cyl","2cyl","3cyl"]

    ## Parameters
    ## internally: length for all parameters in nm
    self.Rdt = 8. # Desmotubule (DT) radius
    self.RdtList = [ ] # List of desmotubule radius (Rdt).
    self.x = 0.5 # particle radius # CF radius (approximate, nm) Called alpha in manuscript
    self.Lpd=100. # PD length. Usually similar to cell wall thickness. 
    self.dPit=120. # centre-to-centre distance between PDs within a pit field
    self.Lneck=25. # neck region length
    self.LneckList=[] # neck region length ; only used in computeVals/computeUnitVals
    self.Lneck2List=[] # neck region length Lneck of "right" neck; only used in computeVals/computeUnitVals if asymmetricPDs == True and doNotCombine == True
    self.diff = 1.62e8 # particle diffusion constant (in nm2/s!) 1.62e8 is used for carboxyfluorescein
    self.Lcell=1e4 # cell length (in nm!) ; does not affect Peff, Fih, etc. The model is only valid if Lcell > spacing between PDs (or pit fields). 
    self.grid="hexagonal" # distribution of PDs on regular grid. Choose from: "triangular", "square", "hexagonal" or "hex" and "random"
    self.gridList=[] # distribution of PDs on regular grid. Choose from: "triangular", "square", "hexagonal" or "hex" and "random"

    self.peList = [1., 6., 8.5, 25. ]              # list of target permeabilities (Peff)
    self.xMaxList = [2. , 2.5, 3., 3.4, 3.5, 4.]   # list of maximum particle sizes (also alpha_bar): corresponding densities will be computed for obtaining target permeabilities
    self.densList=[0.5e-5, 1e-5, 1.3e-5]      # list of target densities: corresponding xMax and Rn (neck radius) values will be computed for obtaining target permeabilities
    self.twinningList=[1,2,3,4,5,6,7,12,19]       # fold-increases for density through (repeated) twinning. Only used if computeClusterIncrease == True. The list should include 1 to compute the reference values before twinning ; So far, only used in computing required xMax / Rn given (starting) densities in densList; Allowed numbers: 1, 2,3,4,5, 6, 7, 12, 19. 
    self.RcList = [ 0. ] # List of Central region radius (Rc). Program uses Rc = Rn if Rc < Rn. 

    ## parameters used in computing correction factor Fih 
    self.pitList=[1,2,3,4,5,6,7, 12,19] ## Number of PDs per pit field. Allowed numbers: 1, 2,3,4,5, 6, 7, 12, 19. ; Only used for computing Fih in pit fields and twinning effect.  
    #self.dPitList = [80, 100,120, 130, 140, 150] ## Center-to-center distance between PDs within pit field
    self.dPitList = [100,120] ## Center-to-center distance between PDs within pit field
    self.LpdList=[100.,200.,500.] # PD length ; the list is only used for computing Fih in pit fields and twinning effect.
    self.pitDens=6e-6 # pit density; only used for computing Fih in pit fields as a function of xMax (computeFih_pitField_xMax)
    self.RnList=[12.,14.] # Neck radius Rn ; the list is only used for computing Fih in pit fields and twinning effect and computeUnitVals.
    self.Rn2List=[] # Neck radius Rn of "right" neck; only used in computeVals/computeUnitVals if asymmetricPDs == True and doNotCombine == True

    ## parameters used in computing Rn, Density graphs
    self.xStep = 0.1
    self.xStart = 1.
    self.xMax = 36.

    ## new parameters for singleInterface and bootstrapInterface
    self.dataSep=","
    self.noData=-1
    self.bootstrapSamples=10000
    self.bootstrapSingleSampleSize=0
    self.CIlevel=0.95
    self.consistentInterfaceLabels="False" # often, it is too hard to identify the orientation of the PD, making neck1 and neck2 arbitrary among PDs in the same sample. --> bootstrapping strategy. 
    self.PDmodel="3cyl" # (in future use this for all inputs? compute from data? maybe sometimes a 3cyl input and a 2cyl model actually is best...)
    self.randomSeed=-1 # initialize random if < 0
    self.DTvarList=[] # parameter required if DTformula == True
    self.DTcoefList=[] # parameter required if DTformula == True
    self.DTintercept=self.Rdt # parameter required if DTformula == True
    self.noiseModelList=[] # currently: A (additive), M (multiplicative) ; all noise is uniform (U). Future: also Guassian/Normal (G or N)?
    self.noiseVarList=[]  # noise may be applied to ....
    self.noiseMagnitudeList=[] # for uniform noise: xi ~ U(-m,m) ; future: for Gaussian noise: xi ~ N(0,m)
    self.downSample="avg"

    ## numerical parameters 
    self.dInc = 1e-7 # step size for finding target density (followed by linear interpolation)
    self.xInc = 0.01 # step size for finding target PD aperture (alpha_bar, Rn) (followed by linear interpolation)
    self.sensIncList = [ 0.001, 0.0001 ] # step sizes for numerically computing derivatives in sensitivity analysis

  def read(self,filename):
    infile = open(filename,'r')
    print ("Opened", filename, "for reading.")
    print ("Note that internally, all lengths are in nm.")
    for line in infile:
      if line[0] == '#':
        continue
      sp = line.split('#')[0].split()
      if len(sp) == 0:
        continue
      if sp[0] in self.simpleFloats:
        setattr(self,sp[0],float(sp[1]))
      elif sp[0] in self.factorFloats:
        setattr(self,sp[0],float(sp[1])*self.factorFloats[sp[0]])
      elif sp[0] in self.listFloats:
        setattr(self,sp[0],[float(s) for s in sp[1:]])
      elif sp[0] in self.listFactorFloats:
        setattr(self,sp[0],[float(s)*self.listFactorFloats[sp[0]] for s in sp[1:]])
      elif sp[0] in self.simpleInts:
        setattr(self,sp[0],int(sp[1]))
      elif sp[0] in self.simpleBools:
        setattr(self,sp[0],bool(int(sp[1])))
      elif sp[0] in self.listStrings:
        setattr(self,sp[0],[s for s in sp[1:]])
      elif sp[0] == "outputType":
        try:
          self.colSep = self.typeDict[sp[1]][0]
          self.fileExt = self.typeDict[sp[1]][1]
          self.outputType = sp[1]
        except:
          print("File type not supported: ", sp[1], ". Using tsv instead") 
          self.colSep = '\t'
          self.fileExt = "tsv"
          self.outputType = "tsv"
      elif sp[0] == "PDmodel":
        if sp[1] in self.PDmodelOptions:  
          setattr(self,sp[0],sp[1])
        else:
          print("PDmodel not recorgnized", sp[1], "using 1cyl (straight channels) instead") 
          setattr(self,sp[0],"1cyl")
      else: ## assume single string
        try:
          setattr(self,sp[0],sp[1])
        except:
          print ("Parameter not recognized: ", sp[0])
      print(sp[0],getattr(self,sp[0]))
    self.completePars()

  def guiFormat(self):
    self.completePars()
    for p in self.factorFloats:
      newVal = getattr(self,p)/self.factorFloats[p]
      setattr(self,p,newVal)
    for p in self.listFactorFloats:
      newVals = getattr(self,p)
      for i in range(len(newVals)):
        newVals[i]/=self.listFactorFloats[p]
      setattr(self,p,newVals)
      #outfile.write('\t'.join([p]+[str(s/self.listFactorFloats[p]) for s in getattr(self,p)])+"\n") 

  def completePars(self):
    if len(self.gridList) == 0:
      self.gridList.append(self.grid)
    if len(self.RdtList) == 0:
      self.RdtList.append(self.Rdt)
    if len(self.LneckList) == 0:
      self.LneckList.append(self.Lneck)
            
  def printAll(self,filename="allParameters.txt"): 
    outfile = open(filename,'w')
    outfile.write("## All parameters are written to this file, possibly including ones that are not used in the computation.\n")
    for p in self.simpleBools:
      outfile.write(p+'\t'+str(int(getattr(self,p)))+"\n") 
    for p in self.simpleInts:
      outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
    for p in self.simpleFloats:
      outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
    for p in self.factorFloats:
      outfile.write(p+'\t'+str(getattr(self,p)/self.factorFloats[p])+"\n") 
    for p in self.simpleStrings:
      outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
    for p in self.listFloats:
      outfile.write('\t'.join([p]+[str(s) for s in getattr(self,p)])+"\n") 
    for p in self.listFactorFloats:
      outfile.write('\t'.join([p]+[str(s/self.listFactorFloats[p]) for s in getattr(self,p)])+"\n") 
    outfile.write("outputType\t"+self.outputType+'\n')
    outfile.write("PDmodel\t"+self.PDmodel+'\n')
    outfile.close()
            
  def writeSingle(self,p,outfile,guiForm=True): 
    if p in self.simpleBools:
      outfile.write(p+'\t'+str(int(getattr(self,p)))+"\n") 
    elif p in self.simpleInts:
      outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
    elif p in self.simpleFloats:
      outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
    elif p in self.factorFloats:
      if guiForm:
        outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
      else:
        outfile.write(p+'\t'+str(getattr(self,p)/self.factorFloats[p])+"\n") 
    elif p in self.simpleStrings + ["outputType","PDmodel"]:
      outfile.write(p+'\t'+str(getattr(self,p))+"\n") 
    elif p in self.listFloats:
      outfile.write('\t'.join([p]+[str(s) for s in getattr(self,p)])+"\n") 
    elif p in self.listFactorFloats:
      if guiForm:
        outfile.write('\t'.join([p]+[str(s) for s in getattr(self,p)])+"\n") 
      else:
        outfile.write('\t'.join([p]+[str(s/self.listFactorFloats[p]) for s in getattr(self,p)])+"\n") 
    elif p in self.listStrings:
      outfile.write('\t'.join([p]+getattr(self,p))+'\n')


  def validate(self):
    if self.singleInterface or self.bootstrapInterface:
      if self.DTformula:
        if len(getattr(self,"DTvarList")) != len(getattr(self,"DTcoefList")):
          print ("Cannot interpret statistical model: unequal list lengths", "DTvarList", "DTcoefList")
          return False
    if self.bootstrapInterface:
      if self.addNoise:
        if len(getattr(self,"noiseModelList")) != len(getattr(self,"noiseMagnitudeList")):
          print ("Cannot interpret noise model: unequal list lengths", "noiseModelList", "noiseMagnitudeList")
          return False
        if len(getattr(self,"noiseVarList")) != len(getattr(self,"noiseMagnitudeList")):
          print ("Cannot interpret noise model: unequal list lengths", "noiseVarList", "noiseMagnitudeList")
          return False
    if self.doNotCombine == True:
      requiredLists = []
      if self.computeUnitVals:
        testList = deepcopy(self.requiredPars["computeUnitVals"])
        if not self.asymmetricPDs:
          testList.remove("Rn2List")
          testList.remove("Lneck2List")
      if self.computeVals:
        testList = deepcopy(self.requiredPars["computeVals"])
        if self.asymmetricPDs:
          testList.remove("xMaxList")
        else:
          testList.remove("RnList")
          testList.remove("Rn2List")
          testList.remove("Lneck2List")
      if self.singleInterface:
        testList = deepcopy(self.requiredPars["singleInterface"])
        testList.remove("DTvarList")
        testList.remove("DTcoefList")
        testList.remove("DTintercept")
        #testList.remove("DTvarList","DTcoefList","DTintercept","noiseVarList", "noiseMagnitudeList","noiseModelList")
      if self.bootstrapInterface:
        testList = deepcopy(self.requiredPars["bootstrapInterface"])
        testList.remove("DTvarList")
        testList.remove("DTcoefList")
        testList.remove("DTintercept")
        testList.remove("noiseMagnitudeList")
        testList.remove("noiseModelList")
      for r in testList:
        if r[-4:] == "List":
          requiredLists.append(r)
      print ("Required list parameters: ", requiredLists)
      lenList = []
      for r in requiredLists: 
        try:
          lenList.append(len(getattr(self,r)))
        except:
          lenList.append(0)
          print (r)
      lenList.sort()
      print(lenList)
      lenSet = set(lenList)
      if len(lenSet) > 2:
        print ("Too many different list lengths for option \"doNotCombine\". Lengths:", lenSet)
        return False
      elif len(lenSet) ==2 and lenList[0] != 1:
        print ("Too many different list lengths for option \"doNotCombine\". Lengths: ", lenSet)
        return False
      elif len(lenSet) == 1 and lenList[0] == 0:
        print ("All lists empty. Nothing can be computed.")
        return False
    return True

class PD():
    #self.Rdt = 8. # Desmotubule (DT) radius
    #self.Lpd=100. # PD length. Usually similar to cell wall thickness. 
    #self.Lneck=25. # neck region length
    #self.Lneck2=27. # 2nd neck region length
    #self.Rn = 14. # neck region radius
    #self.Rn2 = 17. # 2nd neck region radius
    #self.Rc = 23. # central region radius
    def __init__(self,rdt=0,lpd=0,rn=0,ln=0,rc=0,rn2=0,ln2=0):
    #def __init__(self,lpd=0,rn=0,ln=0,rc=0,rn2=0,ln2=0):
      self.Rdt = rdt # Desmotubule (DT) radius ## removed DT from class PD, because it is unlikely that this can be measured for each PD independenlty.
      self.Lpd=lpd # PD length. Usually similar to cell wall thickness. 
      self.Rn = rn # neck region radius
      self.Lneck=ln # neck region length
      self.Rc = rc # central region radius
      self.Rn2 = rn2 # 2nd neck region radius
      self.Lneck2=ln2 # 2nd neck region length
## add stuff.
      #self.Punit = {tag:-1 for tag in ["1cyl", "2cyl", "3cyl"]} # available PD models ## becomes more complicated w/o Rdt value!
      self.Punit={}
      self.isComplete={tag:0 for tag in ["1cyl", "2cyl", "3cyl"]} # available PD models
    def copy(self):
      return PD(self.Rdt,self.Lpd, self.Rn,self.Lneck,self.Rc,self.Rn2,self.Lneck2)
    def checkComplete(self):
      #if self.Rdt > 0 and self.Lpd > 0 and self.Rn > 0 :
      if self.Lpd > 0 and self.Rn > 0 :
        self.isComplete["1cyl"] = 1
        if self.Lneck > 0 and self.Rc > 0:
          self.isComplete["2cyl"] = 1
          if self.Rn2 > 0 and self.Lneck2 >0:
            self.isComplete["3cyl"] = 1
        return 1
      return 0
    def makeComplete(self,PDdata,pars):
      stopList = ["1cyl","2cyl","3cyl"]
      parLists = [["Rn","Lpd"],["Rn","Rc","Lneck"],["Rn2","Lneck2"]]
      for i in range(len(stopList)):
        for par in parLists[i]:
          if getattr(self,par) == 0:
            setattr(self,par,PDdata.samplePar(par))
            while getattr(self,par) == 0:  ## NOTE: inefficient with high % missing values
              setattr(self,par,PDdata.samplePar(par))
        if pars.PDmodel == stopList[i]:
          if pars.DTformula :
            self.calcRdt(pars)
          else:
            setattr(self,"Rdt",PDdata.samplePar("Rdt"))
          return #self
    def print(self):
      print(self.Rdt,self.Lpd,self.Rn,self.Lneck,self.Rc,self.Rn2,self.Lneck2)
      #print("Lpd:",self.Lpd,"Rn:",self.Rn,"Lneck:",self.Lneck,"Rc:",self.Rc,"Rn2:",self.Rn2,"Lneck2:",self.Lneck2)
    def calcRdt(self,pars):
      varPart = 0
      for i in range(len(pars.DTvarList)):
        varPart += getattr(self,pars.DTvarList[i])*pars.DTcoefList[i]
      #print (pars.DTintercept + varPart)
      self.Rdt = pars.DTintercept + varPart
      return self.Rdt
    def addNoise(self,pars): # possibly has to be moved to a different position. 
      if pars.addNoise:
        #print (pars.noiseVarList, pars.noiseModelList, pars.noiseMagnitudeList)
        tempPD = self.copy()
        #tempPD.print()
        for attempt in range(30): ## max number of attempts to find a valid PD
          #print (attempt)
          for i in range(len(pars.noiseVarList)):
            p=pars.noiseVarList[i]
            s=pars.noiseMagnitudeList[i]
            nm = pars.noiseModelList[i]
            #try: 
            val=getattr(self,p)
            if nm == "AU": ## additive uniform
              val += rand.uniform(-s,s)
            if nm == "M": ## log normal: log (X) ~N(0,s)
              val *= rand.lognormvariate(0,s)
            if nm == "AN" or nm == "AG": ## additive normal/gaussian
              val += rand.gauss(0,s)
            setattr(tempPD,p,max(val,0))  ## possibly do more elaborate validatation, e.g. Rdt <= Rn, Rc, Rn2; Lneck < 1/2 Lpd or Lneck + Lneck2 < Lpd
            #except:
              #pass
          #tempPD.print()
          if tempPD.realityCheck(pars):
            self = tempPD 
            #print ("completed at Attempt", attempt)
            break
          ## else: try again.
          print ("Failed attempt for adding noise: ", attempt)
          #self.print()
      return
    def realityCheck(self,pars):
      if self.Rn < self.Rdt:
        return False
      if pars.PDmodel == "3cyl":
        if self.Lneck + self.Lneck2 > self.Lpd :
          return False
        if self.Rc < self.Rdt:
          return False
      if pars.PDmodel == "2cyl":
        if self.Lneck * 2 > self.Lpd:
          return False
        if self.Rc < self.Rdt:
          return False
      return True
    def simplify(self,pars):
      if pars.PDmodel == "3cyl":
        return self
      if pars.PDmodel == "2cyl":
        if pars.downSample == "avg":
          if self.Rn2 > 0 or self.Lneck2 > 0:
            rPD = self.copy()
            if self.Rn2 > 0:
              rPD.Rn = 0.5*(self.Rn + self.Rn2)
            if self.Lneck2 > 0:
              rPD.Lneck = 0.5*(self.Lneck + self.Lneck2)
            return rPD
          return self
        elif pars.downSample == "sample":
          if self.Rn2 > 0:
            if self.Lneck2 > 0:
              if rand.uniform() > 0.5:
                rPD = self.copy()
                rPD.Rn = self.Rn2
                rPD.Lneck = self.Lneck2
                return rPD
          return self
      if pars.PDmodel == "1cyl":
        ## do not use central radius to calculate avg radius or so. 
        print ("simplification to 1cyl model not yet implemented. Using Rn.")
        return self
      print ("PD cannot be simplified: ", pars.PDmodel)
      return self

class PDdata():
  def __init__(self):
    self.RdtAvg=0 ## removed DT from class PD, because it is unlikely that this can be measured for each PD independenlty.
    self.RnAvg=0
    self.Rn2Avg=0
    self.RcAvg=0
    self.LpdAvg=0
    self.LneckAvg=0
    self.Lneck2Avg=0

  def read(self,pars):
    #self.allLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList","Lneck2List","Rn2List"]
    self.allLists=["LpdList", "RnList", "RcList", "LneckList","Lneck2List","Rn2List"]
    self.PDdataDict={tag:[0,0,0,getattr(pars,tag)] for tag in self.allLists}
    ## [nAll,nReal,Avg, ...List]
    for key in self.PDdataDict:
      ll = self.PDdataDict[key]
      ll[0] = len(ll[3])
      for i in range(ll[0]):
        if ll[3][i] == pars.noData:
          ll[3][i] = 0
        if ll[3][i] > 0:
          ll[2] += ll[3][i]
          ll[1] += 1
      print (key)
      if ll[0] > 0:
        ll[2] /= ll[1]
      print (ll)
      print (self.PDdataDict[key])
  
  def processData(self,pars):
    for key in self.PDdataDict:
      setattr(self,key[:-4]+"Avg",self.PDdataDict[key][2])

  def fillRdt(self,pars,val=-1):
    if val>0:
      self.RdtAvg=val
      self.PDdataDict["RdtList"] = [1,1,val,[val]]

  def startup(self,pars,RdtValue=-1):
    self.read(pars)
    stop= self.checkData(pars)
    if stop:
      return stop
    self.processData(pars)
    self.fillRdt(pars,RdtValue)

  def samplePar(self,tag):
    if tag in self.constPars:
      return getattr(self,tag+"Avg")
    else:
      return rand.choice(self.PDdataDict[tag+"List"][3])


class PDdataStructured(PDdata):
  def __init__(self):
    #self.RdtAvg=0 ## removed DT from class PD, because it is unlikely that this can be measured for each PD independenlty.
    #self.RnAvg=0
    #self.Rn2Avg=0
    #self.RcAvg=0
    #self.LpdAvg=0
    #self.LneckAvg=0
    #self.Lneck2Avg=0
    super().__init__()
    self.complete1=[]
    self.nComplete1=0
    self.complete2=[]
    self.nComplete2=0
    self.complete3=[]
    self.nComplete3=0
    self.partial=[]
    self.nPartial=0
    self.hasMissingData = False

  #def read(self,pars):
    ##self.allLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList","Lneck2List","Rn2List"]
    #self.allLists=["LpdList", "RnList", "RcList", "LneckList","Lneck2List","Rn2List"]
    #self.PDdataDict={tag:[0,0,0,getattr(pars,tag)] for tag in self.allLists}
    ### [nAll,nReal,Avg, ...List]
    #for key in self.PDdataDict:
      #ll = self.PDdataDict[key]
      #ll[0] = len(ll[3])
      #for i in range(ll[0]):
        #if ll[3][i] == pars.noData:
          #ll[3][i] = 0
        #if ll[3][i] > 0:
          #ll[2] += ll[3][i]
          #ll[1] += 1
      #print (key)
      #if ll[0] > 0:
        #ll[2] /= ll[1]
      #print (ll)
      #print (self.PDdataDict[key])

  def checkData(self,pars):
    if pars.PDmodel == "1cyl":
      #requiredLists=["LpdList", "RnList", "RdtList", "LneckList"]
      requiredLists=["LpdList", "RnList", "LneckList"]
    if pars.PDmodel == "2cyl":
      #requiredLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList"]
      requiredLists=["LpdList", "RnList", "RcList", "LneckList"]
    if pars.PDmodel == "3cyl":
      #requiredLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList","Lneck2List","Rn2List"]
      requiredLists=["LpdList", "RnList", "RcList", "LneckList","Lneck2List","Rn2List"]
    #self.nPDs = max(len(pars.LpdList),len(pars.RnList),len(pars.RcList),len(pars.RdtList),len(pars.LneckList),len(pars.Lneck2List),len(pars.Rn2List))
    #self.nPDs = max([len(getattr(pars,tag)) for tag in requiredLists])
    self.nPDs = max([self.PDdataDict[tag][0] for tag in requiredLists])
    bad=False
    self.constPars=[]
    self.listPars=[]
    for key in requiredLists:
      if self.PDdataDict[key][0] == 1:
        self.constPars.append(key)
      elif self.PDdataDict[key][0] == self.nPDs:
        self.listPars.append(key)
      else:
        print ("PD data validation: wrong number of items for ", key, ": should be 1 or ", str(self.nPDs), "\n")
        bad=True
    if len(pars.RdtList) == 0 and not pars.DTformula:
      bad=True
    if bad:
      return 1
    return 0
  
  def processData(self,pars):
  ## get data in good format ## class PD.
    #def __init__(self,rdt=0,lpd=0,rn=0,ln=0,rc=0,rn2=0,ln2=0)
    #def __init__(self,lpd=0,rn=0,ln=0,rc=0,rn2=0,ln2=0)
    ## possibly add to class: Punit (-1 if not computed), per PD model; bool isComplete (or list)
  ## identify complete data
    constPD = PD()
    for c in self.constPars:
      setattr(constPD,c[:-4],getattr(pars,c)[0]) ## remove trailing end "List"
      #setattr(self,c[:-4]+"Avg",getattr(pars,c)[0])  ## superfluous.
    for i in range(self.nPDs):
      newPD = deepcopy(constPD)
      for p in self.listPars:
        setattr(newPD,p[:-4],self.PDdataDict[p][3][i])
      if not newPD.checkComplete():
        self.partial.append(newPD)
        self.hasMissingData = True
      else:
        if newPD.isComplete["3cyl"]:
          self.complete3.append(newPD)
        elif newPD.isComplete["2cyl"]:  ## question: include PDs in multiple lists, or just one?
          self.complete2.append(newPD)
        else:
          self.complete1.append(newPD)
    if pars.PDmodel == "3cyl":
      if len(self.complete2) + len(self.complete1) > 0:
        self.hasMissingData = True
    if pars.PDmodel == "2cyl":
      if len(self.complete1) > 0:
        self.hasMissingData = True
    self.nComplete3 = len(self.complete3)
    self.nComplete2 = len(self.complete2)
    self.nComplete1 = len(self.complete1)
    self.nPartial = len(self.partial)
    super().processData(pars)
    #for key in self.PDdataDict:
      #setattr(self,key[:-4]+"Avg",self.PDdataDict[key][2])
    #self.nPDsClass={tag:len(getattr(self,tag)) for tag in ["partial"]+pars.PDmodelOptions} ## TODO if necessary: fix.
    ## HIER: calculate parameter avg for complete PDs for each model; couple to generating compact resampling sets per parameter.
      #constPD.print()    
      #newPD.print()    
    #print(self.partial)
    #print(self.complete1)
    #print(self.complete2)
    #print(self.complete3)

  #def startup(self,pars):
    #self.read(pars)
    #stop= self.checkData(pars)
    #if stop:
      #return stop
    #self.processData(pars)

  def fillRdt(self,pars,val=-1):
    super().fillRdt(pars,val)
    if pars.DTformula:
      self.PDdataDict["RdtList"] = [0,0,0,[]]
      for pd in self.complete1 + self.complete2 + self.complete3 + self.partial:
        try:
          #pd.calcRdt(pars)
          #self.PDdataDict["RdtList"][3].append(pd.Rdt)
          self.PDdataDict["RdtList"][3].append(pd.calcRdt(pars)) ## calcRdt fills pd.Rdt and returns Rdt. 
          self.PDdataDict["RdtList"][0] += 1
          self.PDdataDict["RdtList"][1] += 1
          self.PDdataDict["RdtList"][2] += pd.Rdt
        except:
          pass
      self.PDdataDict["RdtList"][2] /= self.PDdataDict["RdtList"][0] 
      self.RdtAvg = self.PDdataDict["RdtList"][2] 
    elif val > 0:
      for pd in self.complete1 + self.complete2 + self.complete3 + self.partial:
        pd.Rdt = val
    return


  def samplePD(self,model,pars,complete=True):
    if complete:
      #if model == "1cyl":
        #relevant = self.complete3 + self.complete2 + self.complete1
      #if model == "2cyl":
        #relevant = self.complete3 + self.complete2
      #if model == "3cyl":
        #relevant = self.complete3 
      #if pars.DTformula:
        #return deepcopy(rand.choice(relevant)).simplify(pars)
      #return rand.choice(relevant).simplify(pars)
      if pars.DTformula:
        return deepcopy(rand.choice(self.getCompleteList(model))).simplify(pars)
      return rand.choice(self.getCompleteList(model)).simplify(pars)
    else:
      mustSimplify = False
      ni = rand.randrange(self.nPDs)
      if ni < self.nComplete1:
        # 1-cyl model never has to be simplified
        if pars.DTformula:
          return deepcopy(self.complete1[ni])
        return self.complete1[ni]
      #if model == "1cyl":
        #mustSimplify = True
      ni -= self.nComplete1
      if ni < self.nComplete2:
        if pars.DTformula:# or mustSimplify:
          return deepcopy(self.complete2[ni]).simplify(pars)
        return self.complete2[ni].simplify(pars)
      #if model == "2cyl":
        #mustSimplify = True
      ni -= self.nComplete2
      if ni < self.nComplete3:
        if pars.DTformula or mustSimplify:
          return deepcopy(self.complete3[ni]).simplify(pars)
        return self.complete3[ni].simplify(pars)
      ni -= self.nComplete3
      return deepcopy(self.partial[ni]).makeComplete(self,pars)


  def getCompleteList(self,model):
    if model == "1cyl":
      return (self.complete3 + self.complete2 + self.complete1)
    if model == "2cyl":
      return (self.complete3 + self.complete2)
    if model == "3cyl":
      return (self.complete3)
    else:
      print ("Unknown PDmodel: ", model)
      return [] 
    
  def avgPar(self,tag):
    return getattr(self,tag+"Avg")
  def makeAvgPD_all(self):
    return PD(self.RdtAvg,self.LpdAvg,self.RnAvg,self.LneckAvg,self.RcAvg,self.Rn2Avg,self.Lneck2Avg)
    #return PD(self.LpdAvg,self.RnAvg,self.LneckAvg,self.RcAvg,self.Rn2Avg,self.Lneck2Avg)
  def makeAvgPD_complete(self,model):
    nn = 0
    rdt = 0.
    lpd = 0.
    rn = 0.
    lneck= 0.
    rc = 0.
    rn2= 0.
    lneck2 = 0.
    for i in self.complete3:
      rdt += i.Rdt
      lpd += i.Lpd
      rn += i.Rn
      lneck += i.Lneck
      rc += i.Rc
      rn2 += i.Rn2
      lneck2 += i.Lneck2
      nn += 1
    if model == "3cyl":
      return PD(rdt/nn,lpd/nn,rn/nn,lneck/nn,rc/nn,rn2/nn,lneck2/nn)
      #return PD(lpd/nn,rn/nn,lneck/nn,rc/nn,rn2/nn,lneck2/nn)
    if nn > 0:
      rn2/=nn
      lneck2/=nn
    for i in self.complete2:
      rdt += i.Rdt
      lpd += i.Lpd
      rn += i.Rn
      lneck += i.Lneck
      rc += i.Rc
      # assume not specified for 2cyl only models
      #rn2 += i.Rn2
      #lneck2 += i.Lneck2
      nn += 1
    if model == "2cyl":
      return PD(rdt/nn,lpd/nn,rn/nn,lneck/nn,rc/nn,rn2,lneck2)
      #return PD(lpd/nn,rn/nn,lneck/nn,rc/nn,rn2,lneck2)
    if nn > 0:
      lneck /= nn
      rc /= nn
    for i in self.complete1:
      rdt += i.Rdt
      lpd += i.Lpd
      rn += i.Rn
      nn += 1
    return PD(rdt/nn,lpd/nn,rn/nn,lneck,rc,rn2,lneck2)
    #return PD(lpd/nn,rn/nn,lneck,rc,rn2,lneck2)

class PDdataUnstructured(PDdata):
  def __init__(self):
    super().__init__()
    #self.RdtAvg=0 ## removed DT from class PD, because it is unlikely that this can be measured for each PD independenlty. ; later returned because of calcRdt. 
    #self.RnAvg=0
    #self.Rn2Avg=0
    #self.RcAvg=0
    #self.LpdAvg=0
    #self.LneckAvg=0
    #self.Lneck2Avg=0

  def read(self,pars):
    super().read(pars)
    # remove missing data zeros
    for key in self.PDdataDict:
      ll = self.PDdataDict[key]
      changed = False
      while ll[0] > ll[1]:
        ll[3].remove(0)
        ll[0] -= 1
        changed = True
      if changed:
        print (self.PDdataDict[key])
    ##self.allLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList","Lneck2List","Rn2List"]
    #self.allLists=["LpdList", "RnList", "RcList", "LneckList","Lneck2List","Rn2List"]
    #self.PDdataDict={tag:[0,0,0,getattr(pars,tag)] for tag in self.allLists}
    ### [nAll,nReal,Avg, ...List]
    #for key in self.PDdataDict:
      #ll = self.PDdataDict[key]
      #ll[0] = len(ll[3])
      #for i in range(ll[0]):
        #if ll[3][i] == pars.noData:
          #ll[3][i] = 0
        #if ll[3][i] > 0:
          #ll[2] += ll[3][i]
          #ll[1] += 1
      #print (key)
      #if ll[0] > 0:
        #ll[2] /= ll[1]
      #print (ll)
      #print (self.PDdataDict[key])

  def checkData(self,pars):
    if pars.PDmodel == "1cyl":
      #requiredLists=["LpdList", "RnList", "RdtList", "LneckList"]
      requiredLists=["LpdList", "RnList", "LneckList"]
    if pars.PDmodel == "2cyl":
      #requiredLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList"]
      requiredLists=["LpdList", "RnList", "RcList", "LneckList"]
    if pars.PDmodel == "3cyl":
      #requiredLists=["LpdList", "RnList", "RcList", "RdtList", "LneckList","Lneck2List","Rn2List"]
      requiredLists=["LpdList", "RnList", "RcList", "LneckList","Lneck2List","Rn2List"]
    self.nPDs = max([self.PDdataDict[tag][0] for tag in requiredLists])
    bad=False
    self.constPars=[]
    self.listPars=[]
    for key in requiredLists:
      if self.PDdataDict[key][0] == 1:
        self.constPars.append(key)
      elif self.PDdataDict[key][0] > 1:
        self.listPars.append(key)
      else:
        print ("PD data validation: no items for ", key, str(self.nPDs), "\n")
        bad=True
    if bad:
      return 1
    return 0
  
  #def processData(self,pars):
  ### get data in good format ## class PD.
    ##def __init__(self,rdt=0,lpd=0,rn=0,ln=0,rc=0,rn2=0,ln2=0)
    ##def __init__(self,lpd=0,rn=0,ln=0,rc=0,rn2=0,ln2=0)
    ### possibly add to class: Punit (-1 if not computed), per PD model; bool isComplete (or list)
    #for key in self.PDdataDict:
      #setattr(self,key[:-4]+"Avg",self.PDdataDict[key][2])

  def fillRdt(self,pars,val=-1):
    super().fillRdt(pars,val)
    if pars.DTformula:
      self.RdtAvg = self.makeAvgPD_all().calcRdt(pars) 

  #def startup(self,pars):
    #self.read(pars)
    #stop= self.checkData(pars)
    #if stop:
      #return stop
    #self.processData(pars)
  
  #def samplePar(self,tag):
    #if tag in self.constPars:
      #return getattr(self,tag+"Avg")
    #else:
      #return rand.choice(self.PDdataDict[tag+"List"][3])
  def samplePD(self,model,pars,complete=True):
    pd = PD()
    if model == "3cyl":
      parList=["Lpd", "Rn", "Rc", "Lneck","Lneck2","Rn2"]
    if model == "2cyl":
      parList=["Lpd", "Rn", "Rc", "Lneck"]
    if model == "1cyl":
      parList=["Lpd", "Rn"]
    for p in parList:
      setattr(pd,p,self.samplePar(p))
    if pars.DTformula:
      pd.calcRdt(pars)
    else:
      pd.Rdt = self.RdtAvg
    return pd
    
  def avgPar(self,tag):
    return getattr(self,tag+"Avg")
  def makeAvgPD_all(self):
    return PD(self.RdtAvg,self.LpdAvg,self.RnAvg,self.LneckAvg,self.RcAvg,self.Rn2Avg,self.Lneck2Avg)
    #return PD(self.LpdAvg,self.RnAvg,self.LneckAvg,self.RcAvg,self.Rn2Avg,self.Lneck2Avg)


def main(parfile):

  pars = Parameters()

  pars.read(parfile)
  if pars.fileTag:
    #pars.printAll(pars.fileTag+".pars")
    pars.printAll(pars.fileTag+"_pars.txt")
  else:
    pars.printAll()

  ## actual computations
  ## choose functions to compute by commenting / uncommenting the relevant lines below
  ## computed quantities:
  ## alpha_bar: maximum particle size that fits through a model PD; also xMax in this script
  ## Fih (also f_ih): The correction factor for wall permeability at discrete spots. This is a number between 0 and 1
  ## Peff: Effective wall permeability (for diffusive transport through PDs). For straight channels: Peff = Fih * Available PD orifice area / total wall area 

  ## get density conversion factors for comparison with sub-nano channel model (if required)  
  if pars.compSubNano:
    fName="densityConversionFactorTable."+pars.fileExt
    if pars.fileTag:
      fName="densityConversionFactorTable_"+pars.fileTag+"."+pars.fileExt
    outfile=open(fName, 'w')
    sys.stdout = outfile
    cfDict=printConversionFactors(pars.Rdt,pars.x,pars.xMaxList,pars.Lpd,pars.Lneck)
    outfile.close()
    sys.stdout = sys.__stdout__
  else:
    cfDict = {}

  ## Calculate values of Peff, based on input parameters, varying alpha_bar (from xMaxList), density (from densList), Rc insofar > Rn (from RcList), PD/pit (from pitList) dPit (from dPitList)
  ## if compSubNano, only values for full length sub-nano structure are computed. 
  if pars.computeVals:
    computePeffValues(pars)

  ## Calculate values of Punit, based on input parameters, varying R_nr (from RnList), Rc insofar > Rn (from RcList); Fih=1, density = 1 PD/um2
  if pars.computeUnitVals:
    computePunitValues(pars)

  ## Calculate values of Peff for all data using all complete PDs. Also compute Punit for all complete PDs and Peff based on averages of all data (as a reference). 
  if pars.singleInterface:
    analyzeSingleInterface(pars)

  ## Calculate value of Peff using all complete data + confidence intervals based on resampling on all complete data and/or all data.
  if pars.bootstrapInterface:
    performBootstrap(pars)

  ## calculate required densities for Peff targets in peList and alpha_bar values in xMaxList
  ## uses cfDict, so this must be initialized before calling the function below.
  if pars.computeDens:
    fName="requiredDensityTable."+pars.fileExt
    if pars.fileTag:
      fName="requiredDensityTable_"+pars.fileTag+"."+pars.fileExt
    outfile=open(fName, 'w')
    sys.stdout = outfile
    printTargetDensitiesStraightPD(pars.Rdt,pars.x,pars.xMaxList,pars.peList,pars.diff,pars.Lpd,pars.Lcell,pars.grid,pars.Lneck,pars.dInc,cfDict,pars.colSep,False,pars.printRn)
    outfile.close()
    sys.stdout = sys.__stdout__

  if pars.computeAperture:
    fName="requiredApertureTable."+pars.fileExt
    if pars.fileTag:
      fName="requiredApertureTable_"+pars.fileTag+"."+pars.fileExt
    outfile=open(fName, 'w')
    sys.stdout = outfile
    if pars.computeClusterIncrease:
      ## calculate required alpha_bar and Rn values for Peff targets in peList, starting densities in densList and the n-fold increase in clusters from twinningList
      printTargetAperturesStraightPD_twinning(pars.Rdt,pars.x,pars.densList,pars.peList,pars.diff,pars.Lpd,pars.Lcell,pars.grid,pars.Lneck,pars.xInc,pars.twinningList,pars.dPit,pars.compSubNano,pars.colSep)

    else:
      ## calculate required alpha_bar and Rn values for Peff targets in peList and densities in densList
      printTargetAperturesStraightPD(pars.Rdt,pars.x,pars.densList,pars.peList,pars.diff,pars.Lpd,pars.Lcell,pars.grid,pars.Lneck,pars.xInc,pars.compSubNano,pars.colSep)
    outfile.close()
    sys.stdout = sys.__stdout__

  ## function for computing Rn, rho graphs that meet a certain target permeability   
  if pars.computeRnDensityGraph:
    ### compSubNano needs redo: don't use conversion factors but real Fih.; invert curve using printTargetAperturesStraightPD TODO 
    #if compSubNano:
      #cfDictG=printConversionFactors(Rdt,x,arange(xStart,xMax,xStep),Lpd,Lneck,False) # do not print output
    #else:
      #cfDictG = {}
    for i in range(len(pars.peList)):
      fName="Rn_dens_Peff"+str(pars.peList[i])+"_l"+str(int(pars.Lpd))+"."+pars.fileExt
      if pars.fileTag:
        fName="Rn_dens_Peff"+str(pars.peList[i])+"_l"+str(int(pars.Lpd))+"_"+pars.fileTag+"."+pars.fileExt
      outfile=open(fName, 'w')
      sys.stdout = outfile
      printTargetDensitiesStraightPD(pars.Rdt,pars.x,arange(pars.xStart,pars.xMax,pars.xStep),[pars.peList[i]],pars.diff,pars.Lpd,pars.Lcell,pars.grid,pars.Lneck,pars.dInc,{},pars.colSep,True,True) # 2x True at the end: no blank lines within file for easier plotting ; include Rn in output.
      outfile.close()
      sys.stdout = sys.__stdout__

  ## function used for verifying Fih for sub-nano channel model
  # density used is the first number of densList
  if pars.computeFih_subNano:
    for Lpd in pars.LpdList:
      for dens in pars.densList:
        for grid in pars.gridList:
            fName="Fih_subNano_dens"+str(dens*1e6)+"_l"+str(int(Lpd))+"_"+grid+"."+pars.fileExt
            if pars.fileTag:
              fName="Fih_subNano_dens"+str(dens*1e6)+"_l"+str(int(Lpd))+"_"+grid+"_"+pars.fileTag+"."+pars.fileExt
            outfile=open(fName, 'w')
            sys.stdout = outfile
            compFih(0.00001,pars.diff,Lpd,pars.Lcell,dens,grid,pars.Rdt)
            outfile.close()
            sys.stdout = sys.__stdout__
  ## conclusion: Fih for the sub-nano channel model is computed reasonably for alpha_bar up to ~15 nm (density 10 or 13 PDs/um2). Results with larger alpha_bar should not be used.

  ## function used for computing Fih for pit fields
  ## output written to multipe files (tsv); the first line includes the parameters Rn, dPit, Lpd
  for dPit in pars.dPitList:
    for Lpd in pars.LpdList:
      for Rn in pars.RnList:
        fNameTail="_pit_Rn"+str(int(Rn))+"_d"+str(int(dPit))+"_l"+str(int(Lpd))+"."+pars.fileExt
        if pars.fileTag:
          fNameTail="_pit_Rn"+str(int(Rn))+"_d"+str(int(dPit))+"_l"+str(int(Lpd))+"_"+pars.fileTag+"."+pars.fileExt
        if  pars.computeFih_pitField_xMax:
            outfile=open("Fih_xMax"+fNameTail, 'w')
            outfile.write("#Rn = " + str(Rn) + ", dPit = " + str(dPit) + ", Lpd = " + str(Lpd) + '\n')
            sys.stdout = outfile
            ## compute Fih as a function of maximum particle size (alpha_bar)
            compFih_pitField(0.00001,pars.diff,Lpd,pars.Lcell,pars.pitDens,pars.grid,pars.Rdt,pars.pitList,dPit) ## particle size has no impact on Fih, except that Fih is not defined when particles are too large to pass (alpha > alpha_bar; equivalently x > xMax)
            outfile.close()
        if  pars.computeFih_pitField_dens:
            outfile=open("Fih_dens"+fNameTail, 'w')
            outfile.write("#Rn = " + str(Rn) + ", dPit = " + str(dPit) + ", Lpd = " + str(Lpd) + '\n')
            sys.stdout = outfile
            ## compute Fih as a function of rho (density)
            compFih_pitField2(pars.x,pars.diff,Lpd,pars.Lcell,pars.grid,pars.Rdt,pars.pitList,dPit,Rn)
            outfile.close()
        if  pars.computeTwinning:
            outfile=open("twinning"+fNameTail, 'w')
            outfile.write("#Rn = " + str(Rn) + ", dPit = " + str(dPit) + ", Lpd = " + str(Lpd) + '\n')
            sys.stdout = outfile
            ## compute the impact of (repeated) PD twinning (resulting in pit fields with increasing numbers of PDs) as a function of rho (density)
            compTwinningEffect_pitField(pars.x,pars.diff,Lpd,pars.Lcell,pars.grid,pars.Rdt,pars.pitList,dPit,Rn)
            outfile.close()
        sys.stdout = sys.__stdout__


  ## sensitivity analysis
  if pars.sensitivityAnalysis:
    fName="sensitivity_"+pars.fileTag+"."+pars.fileExt
    outfile=open(fName, 'w')
    sys.stdout = outfile
    performSensitivityAnalysis(pars)
    outfile.close()
    sys.stdout = sys.__stdout__
  ## conclusion: Fih for the sub-nano channel model is computed reasonably for alpha_bar up to ~15 nm (density 10 or 13 PDs/um2). Results with larger alpha_bar should not be used.
  return

def performBootstrap(pars):
  fName="bootstrapInterface_"+pars.fileTag+"."+pars.fileExt
  outfile=open(fName, 'w')
  sys.stdout = outfile
  if pars.DTformula:
    #RdtList = []
    pars.RdtList = [pars.RdtList[0]] ## better done in a check step. 
  if pars.randomSeed >= 0:
    rand.seed(pars.randomSeed)
  else:
    rand.seed()
  if pars.doNotCombine:
    dataPDs={rdt:PDdataStructured() for rdt in pars.RdtList}
  else:
    dataPDs={rdt:PDdataUnstructured() for rdt in pars.RdtList}
  for rdt in pars.RdtList:
    dataPDs[rdt].startup(pars,rdt)
  if pars.doNotCombine:
    print ("\n## Resampling complete PDs")
  else:
    print ("\n## Unstructured Resampling")
  PeffDict,RdtList = bootstrapResample(dataPDs,pars,True)
  bootstrapReport(dataPDs,PeffDict,RdtList,pars)
  if pars.doNotCombine and dataPDs[pars.RdtList[0]].hasMissingData:
    print ("\n## Resampling all data")
    PeffDict,RdtList = bootstrapResample(dataPDs,pars,False)
    bootstrapReport(dataPDs,PeffDict,RdtList,pars)
  outfile.close()
  sys.stdout = sys.__stdout__
  return 

def bootstrapResample(dataPDs,pars,completeOnly):
  ab=[0,0]
  getAB(ab,pars.grid)
  distDict=getDistDict() 
  PeffDict={rdt:[] for rdt in pars.RdtList}
  RdtList = []
  if pars.bootstrapSingleSampleSize <= 0:
    samSize =  dataPDs[pars.RdtList[0]].nPDs
  else:
    samSize =  pars.bootstrapSingleSampleSize 
  print ("#Sample size:", samSize)
  for rdt in pars.RdtList:
    for n in range(pars.bootstrapSamples):
      resampledPDs = []
      ## resample interface
      for nn in range(samSize):
        resampledPDs.append(dataPDs[rdt].samplePD(pars.PDmodel,pars,completeOnly)) 
        if pars.addNoise:
          resampledPDs[-1].addNoise(pars)
      ## compute Peff (compact format)
      PeffDict[rdt].append(computePeffInterface(resampledPDs,pars,ab[0],ab[1],distDict))
      ## if calculated, collect Rdt values. Assume all can be piled over all samples
      if pars.DTformula:
        RdtList += [pd.calcRdt(pars) for pd in resampledPDs]
  return PeffDict,RdtList

def bootstrapReport(dataPDs,PeffDict,RdtList,pars):
  ## compute confidence intervals
  if pars.DTformula:
    RdtList.sort() 
    print ("\nRdt Values: (percentile, value)")
    percentiles = [ (1.-pars.CIlevel)*0.5, 0.25, 0.5, 0.75, pars.CIlevel+(1.-pars.CIlevel)*0.5] 
    nVal = len(RdtList)
    for p in percentiles:
      print ( p*100, RdtList[int(p*nVal)])
  for rdt in pars.RdtList:
    PeffDict[rdt].sort()
  ## TODO: optional: implement more refined methods for calculating n% CI
  lowID = int((1.-pars.CIlevel)*0.5*pars.bootstrapSamples) ## TODO: check docs to verify this is rounding down
  medianID = int(pars.bootstrapSamples/2)
  highID = int((pars.CIlevel+(1.-pars.CIlevel)*0.5)*pars.bootstrapSamples) ## TODO: check docs to verify this is rounding down
  for rdt in pars.RdtList:
    if pars.DTformula:
      print ("#\nRdt from statistical model")
    else:
      print ("#\nRdt", rdt)
    printSingleInterfaceHeader(pars)
    #print ("\n## Densities used in calculation of Peff (PDs/um2):", pars.colSep.join([str(i/pars.listFactorFloats["densList"]) for i in pars.densList])) 
    print ("Median bootstrap Peff values n=",pars.bootstrapSamples,":", pars.colSep.join([str(i) for i in PeffDict[rdt][medianID] for rdt in pars.RdtList]))
    print ("Bootstrap Peff Lower bound (CL: ",(1.-pars.CIlevel)*0.5,"): ", pars.colSep.join([str(i) for i in PeffDict[rdt][lowID] for rdt in pars.RdtList]))
    print ("Bootstrap Peff Upper bound (CL: ",pars.CIlevel+(1.-pars.CIlevel)*0.5,"): ", pars.colSep.join([str(i) for i in PeffDict[rdt][highID] for rdt in pars.RdtList]))
  return

def analyzeSingleInterface(pars):
  fName="singleInterface_"+pars.fileTag+"."+pars.fileExt
  outfile=open(fName, 'w')
  sys.stdout = outfile
  ab=[0,0]
  getAB(ab,pars.grid)
  distDict=getDistDict()
  if pars.DTformula:
    RdtList = []
    pars.RdtList = [pars.RdtList[0]] ## better done in a check step. 
  PeffDict={rdt:[] for rdt in pars.RdtList}
  if pars.doNotCombine:
    #dataPDs=PDdataStructured()
    dataPDs={rdt:PDdataStructured() for rdt in pars.RdtList}
  else:
    dataPDs={rdt:PDdataUnstructured() for rdt in pars.RdtList}
  for rdt in pars.RdtList:
    dataPDs[rdt].startup(pars,rdt)
  #if pars.doNotCombine:
    #dataPDs=PDdataStructured()
  #else:
    #dataPDs=PDdataUnstructured()
  #dataPDs.startup(pars)
  #printSingleInterfaceHeader(pars)
  for rdt in pars.RdtList:
    if pars.doNotCombine:
      completePDs=dataPDs[rdt].getCompleteList(pars.PDmodel)
      PeffAll = computePeffInterface(completePDs,pars,ab[0],ab[1],distDict,True)
      ## compute Punit for all complete PDs, incl relative contribution to Peff
      # done in above function
      ## compute Peff based on all complete PDs
      #print ("\n## Densities used in calculation of Peff (PDs/um2):", pars.colSep.join([str(i/pars.listFactorFloats["densList"]) for i in pars.densList])) 
  #if pars.doNotCombine:
      printSingleInterfaceHeader(pars)
      print ("Peff based on all the complete PDs:", pars.colSep.join([str(i) for i in PeffAll])) 
      dataPDs[rdt].makeAvgPD_complete(pars.PDmodel).print()
      print ("Peff based on averages of complete PDs:", pars.colSep.join([str(i) for i in computePeffInterface([dataPDs[rdt].makeAvgPD_complete(pars.PDmodel)],pars,ab[0],ab[1],distDict)]))
  #else: (applied only to printing of data)
    dataPDs[rdt].makeAvgPD_all().print()
    print ("Peff based on averages of all data:", pars.colSep.join([str(i) for i in computePeffInterface([dataPDs[rdt].makeAvgPD_all()],pars,ab[0],ab[1],distDict)]))
    #Peff_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Rpit,nPit):
    #Peff_pitField_asym(x,Rn,Rn2,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Lneck2,Rpit,nPit):
  ## compute Peff based on per parameter averages (complete data sets only + all data)
  ## report other useful blah.
  outfile.close()
  sys.stdout = sys.__stdout__
  return

def printSingleInterfaceHeader(pars):
  #HIER
  #if pars.DTformula:
    #print ("\n##Densities used in calculation of Peff (PDs/um2_nm):", pars.colSep.join([str(i/pars.listFactorFloats["densList"]) for i in pars.densList ])) 
  #else:
    #print ("\n##Densities_Rdt used in calculation of Peff (PDs/um2_nm):", pars.colSep.join([str(i/pars.listFactorFloats["densList"])+"_"+str(j) for i in pars.densList for j in pars.RdtList])) 
    print ("\n##Densities_clusterSize_dPit used in calculation of Peff (PDs/um2_._nm):", pars.colSep.join([str(i/pars.listFactorFloats["densList"])+"_"+str(j)+"_"+str(k) for i in pars.densList for j in pars.pitList for k in pars.dPitList])) 

def computePeffInterface(PDlist,pars,A,B,distDict,verbose=False):
  PunitList=[] ## because Rdt no longer in class PD
  #if pars.DTformula:
    #PunitList=[[]] ## because Rdt no longer in class PD
  #else:
    #PunitList=[[] for rdt in pars.RdtList]
  ## TODO: RnList, LpdList not necessary as lists, counters also possible
  RnList=[]
  LpdList=[]
  if len(PDlist) == 0:
    return 0
  for pd in PDlist:
    PunitList.append(calcPunit(pd,pars)) 
    #if pars.DTformula:
      #PunitList[0].append(calcPunit(pd,pd.calcRdt(pars),pars)) ## !!!!! TODO issue here with "add noise"
    #else:
      #for i in range(len(pars.RdtList)):
        #PunitList[i].append(calcPunit(pd,pars.RdtList[i],pars))
    RnList.append(calcFihRn(pd,pars.PDmodel)) 
    LpdList.append(pd.Lpd)
  #PunitSum=[sum(PunitList[i]) for i in range(len(PunitList))]
  PunitSum=sum(PunitList)
  Peff = []
  for dens in pars.densList:
    for pit in pars.pitList:
      for dPit in pars.dPitList:
        Rpit = getRpitDefault(distDict,mean(RnList),dPit,pit) 
        if pars.computeClusterIncrease:
          Peff.append( mean(PunitList) * dens*1e6 * Fih_pitField(pars.x,mean(RnList),pars.diff,mean(LpdList),pars.Lcell,dens,A,B,Rpit,pit)) # avg Punit * dens * Fih(avgRn) ## because Rdt no longer in class PD ## probably undesired behaviour: total dens == dens * pit. 
        else:
          Peff.append( mean(PunitList) * dens/pit*1e6 * Fih_pitField(pars.x,mean(RnList),pars.diff,mean(LpdList),pars.Lcell,dens/pit,A,B,Rpit,pit)) # avg Punit * dens * Fih(avgRn) ## because Rdt no longer in class PD ## probably desired behaviour: total dens = dens
    #for i in range(len(PunitList)):
      #Peff.append( mean(PunitList[i]) * dens*1e6 * Fih_pitField(pars.x,mean(RnList),pars.diff,mean(LpdList),pars.Lcell,dens,A,B,Rpit,pars.pitList[0])) # avg Punit * dens * Fih(avgRn)
    ## remember: density is /nm2
    #Peff_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Rpit,nPit):
    #Peff_pitField_asym(x,Rn,Rn2,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Lneck2,Rpit,nPit):
    ## conclusion: write Fih function.
  if verbose:
    # output PunitList
    print("\n## Punit for all individual PDs")
    print("ID", "Punit", "fraction")
    if pars.DTformula:
      print("## DTformula: Rdt = ",pars.DTintercept, " ".join([" ".join(["+",str(pars.DTcoefList[i]), "*",pars.DTvarList[i]]) for i in range(len(pars.DTvarList))]))
      # old; resurrected: Rdt no longer in class PD
      for i in range(len(PDlist)):
        print(str(i), str(PunitList[i]), str(PunitList[i]/PunitSum), PDlist[i].Rdt) # TODO: add PD info
    else:
      print("## Rdt: ", PDlist[0].Rdt)
      # old; resurrected: Rdt no longer in class PD
      for i in range(len(PDlist)):
        print(str(i), str(PunitList[i]), str(PunitList[i]/PunitSum)) # TODO: add PD info
    ## removed: Rdt  back in class PD
    #for j in range(len(PunitList)):
      #if pars.DTformula:
        #print("## DTformula: Rdt = ",pars.DTintercept, " ".join([" ".join(["+",str(pars.DTcoefList[i]), "*",pars.DTvarList[i]]) for i in range(len(pars.DTvarList))]))
        #for i in range(len(PDlist)):
          #print(str(i), str(PunitList[j][i]), str(PunitList[j][i]/PunitSum[j]), str(PDlist[i].calcRdt(pars))) # TODO: add PD info
      #else:
        #print("## Rdt: ", pars.RdtList[j])
        #for i in range(len(PDlist)):
          #print(str(i), str(PunitList[j][i]), str(PunitList[j][i]/PunitSum[j])) # TODO: add PD info
  return Peff

## TODO: make function of PD
def calcFihRn(pd,model):
  # returns: consistent with other data: just avg 2 necks if given. 
  if model == "3cyl":
    if pd.Rn2 == 0: ## other option: only use Rn2 if > Rdt
      return pd.Rn
    else:
      return (pd.Rn+pd.Rn2)*0.5
  return pd.Rn

## TODO: make function of PD
def calcBestRn(pd,model):
  # returns: best weighted avg based on model and all available data
  ## actually better: 1/x avg: Lpd / (sum (Lx/Rx)).
  ## or more consistent with other data: just avg 2 necks if given. 
  if model == "1cyl":
    if pd.Rn2 == 0:
      if pd.Rc == 0:
        return pd.Rn
      else:
        #return (pd.Rn*pd.Lneck*2 + pd.Rc*(pd.Lpd-2*pd.Lneck))/pd.Lpd
        return pd.Lpd/(pd.Lneck*2/pd.Rn  + (pd.Lpd-2*pd.Lneck)/pd.Rc)
    else: # assume Rn2 defined means Rc also defined.
      #return (pd.Rn*pd.Lneck + pd.Rn2*pd.Lneck2 + pd.Rc*(pd.Lpd-pd.Lneck-pd.Lneck2))/pd.Lpd
      return pd.Lpd/(pd.Lneck/pd.Rn + pd.Lneck2/pd.Rn2 + (pd.Lpd-pd.Lneck-pd.Lneck2)/pd.Rc)
  elif model == "2cyl":
    if pd.Rn2 == 0:
      return pd.Rn
    else:
      #return (pd.Rn*pd.Lneck + pd.Rn2*pd.Lneck2)/(pd.Lneck + pd.Lneck2)
      return (pd.Lneck + pd.Lneck2)/(pd.Lneck/pd.Rn + pd.Lneck2/pd.Rn2)
  elif model == "3cyl":
    #return (pd.Rn*pd.Lneck + pd.Rn2*pd.Lneck2)/(pd.Lneck + pd.Lneck2)
    return (pd.Lneck + pd.Lneck2)/(pd.Lneck/pd.Rn + pd.Lneck2/pd.Rn2)


## TODO: make function of PD
#def calcPunit(pd,rdt,pars):
def calcPunit(pd,pars):
  ## this goes wrong with multiple values for pars.whatever (used: x, D) )
#  ## Rdt removed from class PD for reasons of experimental difficulties. Entered separately. 
#  ## option: add value to pd class -> efficiency (if many complete records)
  #tag=pars.PDmodel+"_"+str(rdt)
  tag=pars.PDmodel
  if tag in pd.Punit:
    return pd.Punit[tag] 
  if pars.PDmodel == "1cyl":
    #if pd.Punit["1cyl"] >= 0: 
      #return pd.Punit["1cyl"] 
    pe = Punit(pars.x,pd.Rn,pars.diff,pd.Lpd,pd.Rdt,pd.Rn,pd.Lpd/4.) # Lneck is arbitrary here, so do not assume it is defined.
    #pe = Punit(pars.x,pd.Rn,pars.diff,pd.Lpd,rdt,pd.Rn,pd.Lpd/4.) # Lneck is arbitrary here, so do not assume it is defined.
    #pd.Punit["1cyl"]  = pe
    pd.Punit[tag]  = pe
    return pe
  elif pars.PDmodel == "2cyl":
    #if pd.Punit["2cyl"] >= 0: 
      #return pd.Punit["2cyl"] 
    pe = Punit(pars.x,pd.Rn,pars.diff,pd.Lpd,pd.Rdt,pd.Rc,pd.Lneck)
    #pe = Punit(pars.x,pd.Rn,pars.diff,pd.Lpd,rdt,pd.Rc,pd.Lneck)
    #pd.Punit["2cyl"]  = pe
    pd.Punit[tag]  = pe
    return pe
  elif pars.PDmodel == "3cyl":
    #if pd.Punit["3cyl"] >= 0: 
      #return pd.Punit["3cyl"] 
    pe = PunitAsym(pars.x,pd.Rn,pd.Rn2,pars.diff,pd.Lpd,pd.Rdt,pd.Rc,pd.Lneck,pd.Lneck2) 
    #pe = PunitAsym(pars.x,pd.Rn,pd.Rn2,pars.diff,pd.Lpd,rdt,pd.Rc,pd.Lneck,pd.Lneck2) 
    #pd.Punit["3cyl"]  = pe
    pd.Punit[tag]  = pe
    return pe
  else:
    print ("PD model", pars.PDmodel, "not defined.")
    return 0


def performSensitivityAnalysis(pars):
  ab=[0,0]
  getAB(ab,pars.grid)
  print("## Elasticity := parameterValue / baseValue * partial derivative (=delta-BaseValue/delta-ParameterValue)")
  for base_Rn in pars.RnList:
    for base_Rdt in pars.RdtList:
      for base_Rc in pars.RcList:
        if base_Rc < base_Rn:
          if base_Rc < 0.0000001:
            base_Rc = base_Rn
          else:
            continue
        for base_Lpd in pars.LpdList:
          for base_dens in pars.densList:
            base_Lneck = pars.Lneck
            base_x = pars.x
            print("## Parameter elasticities around alpha=",base_x,"Rn=", base_Rn, ",Rdt=", base_Rdt, ",Rc=", base_Rc, ",Lpd=", base_Lpd, ",Lneck=", base_Lneck, ",Rho=", base_dens*1e6 ) 
            if base_Rn == base_Rc:
              print ("#Peff"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rn+Rc","Rn","Rdt","Rc","Lpd","Lneck","Rho"] for inc in pars.sensIncList))
            else:
              print ("#Peff"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rn","Rdt","Rc","Lpd","Lneck","Rho"] for inc in pars.sensIncList))
            outline = []
            baseVal =Peff(pars.x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc,base_Lneck)
            outline.append(baseVal)
            for inc in pars.sensIncList:
              outline.append(base_x/baseVal*(Peff(base_x+inc,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc,base_Lneck)-baseVal)/inc)
            if base_Rn == base_Rc:
              for inc in pars.sensIncList:
                outline.append(base_Rn/baseVal*(Peff(base_x,base_Rn+inc,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc+inc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rn/baseVal*(Peff(base_x,base_Rn+inc,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rdt/baseVal*(Peff(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt+inc,base_Rc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rc/baseVal*(Peff(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc+inc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Lpd/baseVal*(Peff(base_x,base_Rn,pars.diff,base_Lpd+inc,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Lneck/baseVal*(Peff(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt,base_Rc,base_Lneck+inc)-baseVal)/inc)
            for inc in pars.sensIncList:
              incDens=1e-6*inc
              outline.append(base_dens/baseVal*(Peff(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens+incDens,ab[0],ab[1],base_Rdt,base_Rc,base_Lneck)-baseVal)/incDens)
            print(pars.colSep.join(str(i) for i in outline))

            if base_Rn == base_Rc:
              print ("#Punit"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rn+Rc","Rn","Rdt","Rc","Lpd","Lneck"] for inc in pars.sensIncList))
            else:
              print ("#Punit"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rn","Rdt","Rc","Lpd","Lneck"] for inc in pars.sensIncList))
            outline = []
            baseVal =Punit(base_x,base_Rn,pars.diff,base_Lpd,base_Rdt,base_Rc,base_Lneck)
            #Punit(x,Rn,D,Lpd,Rdt,Rc,Lneck):
            outline.append(baseVal)
            for inc in pars.sensIncList:
              outline.append(base_x/baseVal*(Punit(base_x+inc,base_Rn,pars.diff,base_Lpd,base_Rdt,base_Rc,base_Lneck)-baseVal)/inc)
            if base_Rn == base_Rc:
              for inc in pars.sensIncList:
                outline.append(base_Rn/baseVal*(Punit(base_x,base_Rn+inc,pars.diff,base_Lpd,base_Rdt,base_Rc+inc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rn/baseVal*(Punit(base_x,base_Rn+inc,pars.diff,base_Lpd,base_Rdt,base_Rc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rdt/baseVal*(Punit(base_x,base_Rn,pars.diff,base_Lpd,base_Rdt+inc,base_Rc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rc/baseVal*(Punit(base_x,base_Rn,pars.diff,base_Lpd,base_Rdt,base_Rc+inc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Lpd/baseVal*(Punit(base_x,base_Rn,pars.diff,base_Lpd+inc,base_Rdt,base_Rc,base_Lneck)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Lneck/baseVal*(Punit(base_x,base_Rn,pars.diff,base_Lpd,base_Rdt,base_Rc,base_Lneck+inc)-baseVal)/inc)
            print(pars.colSep.join(str(i) for i in outline))

            print ("#Hindrance(neck)"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rn","Rdt"] for inc in pars.sensIncList))
            outline = []
            #hLam(Rx,aa,Rdt):
            baseVal =hLam(base_Rn,base_x,base_Rdt)
            outline.append(baseVal)
            for inc in pars.sensIncList:
              outline.append(base_x/baseVal*(hLam(base_Rn,base_x+inc,base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rn/baseVal*(hLam(base_Rn+inc,base_x,base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rdt/baseVal*(hLam(base_Rn,base_x,base_Rdt+inc)-baseVal)/inc)
            print(pars.colSep.join(str(i) for i in outline))

            print ("#Hindrance(central_region)"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rc","Rdt"] for inc in pars.sensIncList))
            outline = []
            #hLam(Rx,aa,Rdt):
            baseVal =hLam(base_Rc,base_x,base_Rdt)
            outline.append(baseVal)
            for inc in pars.sensIncList:
              outline.append(base_x/baseVal*(hLam(base_Rc,base_x+inc,base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rc/baseVal*(hLam(base_Rc+inc,base_x,base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rdt/baseVal*(hLam(base_Rc,base_x,base_Rdt+inc)-baseVal)/inc)
            print(pars.colSep.join(str(i) for i in outline))

            print ("#Fih"+pars.colSep+ pars.colSep.join(par+'_'+str(inc) for par in ["alpha","Rn","Rdt","Lpd","Rho"] for inc in pars.sensIncList))
            outline = []
            #getFih_sleeve(x,xMax,D,Lpd,Lcell,dens,A,B,Rdt):
            baseVal =getFih_sleeve(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt)
            outline.append(baseVal)
            for inc in pars.sensIncList:
              outline.append(base_x/baseVal*(getFih_sleeve(base_x+inc,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rn/baseVal*(getFih_sleeve(base_x,base_Rn+inc,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Rdt/baseVal*(getFih_sleeve(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt+inc)-baseVal)/inc)
            for inc in pars.sensIncList:
              outline.append(base_Lpd/baseVal*(getFih_sleeve(base_x,base_Rn,pars.diff,base_Lpd+inc,pars.Lcell,base_dens,ab[0],ab[1],base_Rdt)-baseVal)/inc)
            for inc in pars.sensIncList:
              incDens=1e-6*inc
              outline.append(base_dens/baseVal*(getFih_sleeve(base_x,base_Rn,pars.diff,base_Lpd,pars.Lcell,base_dens+incDens,ab[0],ab[1],base_Rdt)-baseVal)/incDens)
            print(pars.colSep.join(str(i) for i in outline))
  return


def fdt(Rdt,x) :
  # the same as n_c in main text
  return ((Rdt+2.*x)**2-Rdt**2)/x**2

def lam(Rx,aa,Rdt) :
  # relative particle size for the default (unobstruced sleeve) model
  return 2.*aa/(Rx-Rdt)


def hSlit(x):
  # hindrance factor for slit pore
  if x >= 1:
    return 0
  return 1 + 9./16.*x*log(x) - 1.19358*x + 0.4285*x**3 - 0.3192*x**4 + 0.08428*x**5 

def hLam(Rx,aa,Rdt):
  # hindrance based on (relative) particle size; geometry with DT
  return hSlit(lam(Rx,aa,Rdt))


def hCyl(x) :
  # hindrance factor for cylindrical channel
  if  x >= 1:
    return 0
  if  x < 0.95: 
    return 1 + 9./8.*x*log(x) - 1.56034*x + 0.528155*x**2+ 1.91521*x**3 - 2.81903*x**4 + 0.270788*x**5 + 1.10115*x**6 - 0.435933*x**7
  return (1.-x)**2*(0.984*sqrt((1.-x)/x)**5)


def hLamCyl(Rx,aa):
  # hindrance based on (relative) particle size; cylinder geometry 
  return hCyl(aa/Rx)



def gdtH(Rdt,x,xMax):
  # the same as ttA_n,dt / ttA_n,circle in main text
  if x<xMax :
    return (hLam(Rdt+2.*xMax,x,Rdt)*((Rdt+2.*xMax )**2-(Rdt)**2)/(hLamCyl(xMax,x)*(xMax)**2))/fdt(Rdt,xMax)
  return VERY_MUCH

def printConversionFactors(Rdt,x,xMaxList,Lpd,Lneck,printOut=True,colSep='\t'):
  if printOut:
    print("## Required density conversion factors for sub-nano model and derivates")
    print("alpha"+colSep+"alpha_bar"+colSep+"sub-nano"+colSep+"sub-nano_neckOnly"+colSep+"spoke_entrance")
  data={}
  for xMax in xMaxList:
    a = gdtH(Rdt,x,xMax)*fdt(Rdt,xMax)/9. # assume 9 sub-nano channels per PD. 
    fNeck= (2.*(Lneck+x)/Lpd)
    fGate= (2.*(1.+x)/Lpd)
    if printOut:
      print(x, xMax, a, a*(fNeck +(1.-fNeck)/a), a*(fGate+(1.-fGate)/a))
    data[xMax] =[ a, a*(fNeck +(1.-fNeck)/a), a*(fGate+(1.-fGate)/a) ]
  if printOut:
    print("")
  return data

def getAB(ab,gridType):
    #after Berezhkovskii et at 2006, DOI: 10.1063/1.2161196
    # these are the values a and b used in function fSig, which are based on numerical computation for different layouts of traps 
    AAs = 1.75; BBs = 2.02 # 2006 square
    AAh = 1.37; BBh = 2.59 # 2006 hex
    AAt = 1.62; BBt = 1.36 # 2006 triangular
    AAr = 0.34; BBr = -0.58 # 2006 random
    AAc = 1.37; BBc = 0.37 # 2006 cylinder ## not used in this script: cylinders cannot produce a space filling packing of a plane.
    if gridType == "square":
        ab[0] = AAs
        ab[1] = BBs
    elif gridType == "hex" or gridType == "hexagonal":
        ab[0] = AAh
        ab[1] = BBh
    elif gridType == "triangular":
        ab[0] = AAt
        ab[1] = BBt
    elif gridType == "random":
        ab[0] = AAr
        ab[1] = BBr
    return

def fSig(s,a,b):
  #after Berezhkovskii et at 2006, DOI: 10.1063/1.2161196
  return (1+a*sqrt(s)-b*s**2)/(1-s)**2


def ttA(Rx,aa,Rdt):
  # available cross section area, hindrance corrected, with DT
  return hLam(Rx,aa,Rdt)*pi*(Rx**2-Rdt**2)

def ttA_cyl(Rx,aa):
  # available cross section area, hindrance corrected, no DT
  return hLamCyl(Rx,aa)*pi*(Rx**2)

def tNot_MakTrap(x,Rn,D,ltt,lc,dens,A,B):
  # diffusion time tau|| from calculation of the f_ih correction factor
  return 1./(2.*D)*((lc-2.*x)**2+ (ltt+2.*x)**2 +2.*D*((lc-2.*x)/(dens) + (ltt+2.*x)*hLamCyl(Rn,x)*Rn**2*pi) /(4.*D*hLamCyl(Rn,x)*Rn*fSig(dens*(Rn-x)**2*pi,A,B)) +(lc-2.*x)*(ltt+2.*x)*(hLamCyl(Rn,x)*Rn**2*pi*dens+ 1./(hLamCyl(Rn,x)*Rn**2*pi*dens)))

def Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B):
  # effective permeability based on tau||. Used in calculation of the f_ih correction factor
  return Lcell/(2.*tNot_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B) - Lcell**2/D)

def tNot_MakTrap_multiCyl(x,Rn,D,ltt,lc,dens,A,B,xMax,nCyl):
  # diffusion time tau|| from calculation of the f_ih correction factor
  return 1./(2.*D)*((lc-2.*x)**2+ (ltt+2.*x)**2 +2.*D*((lc-2.*x)/(dens) + (ltt+2.*x)*(nCyl*hLamCyl(xMax,x))*xMax**2*pi) /(4.*D*hLamCyl(Rn,x)*Rn*(nCyl*(xMax-x)**2/(Rn-x)**2)*fSig(dens*(Rn-x)**2*pi,A,B)) +(lc-2.*x)*(ltt+2.*x)*(nCyl*hLamCyl(xMax,x)*xMax**2*pi*dens+ 1./(nCyl*hLamCyl(xMax,x)*xMax**2*pi*dens)))

def Peff_MakTrap_multiCyl(x,Rn,D,Lpd,Lcell,dens,A,B,xMax,nCyl):
  # effective permeability based on tau||. Used in calculation of the f_ih correction factor
  return Lcell/(2.*tNot_MakTrap_multiCyl(x,Rn,D,Lpd,Lcell,dens,A,B,xMax,nCyl) - Lcell**2/D)

def tNot_MakTrap_pitField(x,Rn,D,ltt,lc,dens,A,B,Rpit,nCyl):
  # diffusion time tau|| from calculation of the f_ih correction factor
  #return 1./(2.*D)*((lc-2.*x)**2+ (ltt+2.*x)**2 +2.*D*((lc-2.*x)/(dens) + (ltt+2.*x)*(nCyl*hLamCyl(xMax,x))*xMax**2*pi) /(4.*D*hLamCyl(Rn,x)*Rn*(nCyl*(xMax-x)**2/(Rn-x)**2)*fSig(dens*(Rn-x)**2*pi,A,B)) +(lc-2.*x)*(ltt+2.*x)*(nCyl*hLamCyl(xMax,x)*xMax**2*pi*dens+ 1./(nCyl*hLamCyl(xMax,x)*xMax**2*pi*dens))) ### multicyl
  return 1./(2.*D)*((lc-2.*x)**2+ (ltt+2.*x)**2 +2.*D*((lc-2.*x)/(dens) + (ltt+2.*x)*nCyl*hLamCyl(Rn,x)*Rn**2*pi) /(4.*D*hLamCyl(Rpit,x)*Rpit*(nCyl*(Rn-x)**2/(Rpit-x)**2)*fSig(dens*(Rpit-x)**2*pi,A,B)) +(lc-2.*x)*(ltt+2.*x)*(nCyl*hLamCyl(Rn,x)*Rn**2*pi*dens+ 1./(nCyl*hLamCyl(Rn,x)*Rn**2*pi*dens)))

def Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,nCyl):
  # effective permeability based on tau||. Used in calculation of the f_ih correction factor
  return Lcell/(2.*tNot_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,nCyl) - Lcell**2/D)

def Peff_naive(x,Rn,D,Lpd,dens,Rdt,Rc,Lneck):
  # effective permeability without correction factor f_ih, including desmotubule and Rc >= Rn
  #return dens*D*ttA(Rn,x,Rdt)*ttA(Rc,x,Rdt)/(2.*(Lneck+x)*ttA(Rc,x,Rdt)+(Lpd-2.*(Lneck+x))*ttA(Rn,x,Rdt) )
  if Rdt + 2.*x >= Rn:
    return 0.
  return dens*D*ttA(Rn,x,Rdt)*ttA(Rc,x,Rdt)/(2.*(Lneck+x)*ttA(Rc,x,Rdt)+(Lpd-2.*(Lneck+x))*ttA(Rn,x,Rdt) )

def Punit(x,Rn,D,Lpd,Rdt,Rc,Lneck):
  # effective permeability without correction factor f_ih, assuming unit density (rho = 1 PD/um2), including desmotubule and Rc >= Rn
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  return 1e-3*Peff_naive(x,Rn,D,Lpd,1e-6,Rdt,Rc,Lneck)

def Peff_naive_asym(x,Rn,Rn2,D,Lpd,dens,Rdt,Rc,Lneck,Lneck2):
  # effective permeability without correction factor f_ih, including desmotubule and Rc >= Rn
  #return dens*D*ttA(Rn,x,Rdt)*ttA(Rc,x,Rdt)/(2.*(Lneck+x)*ttA(Rc,x,Rdt)+(Lpd-2.*(Lneck+x))*ttA(Rn,x,Rdt) )
  if Rdt + 2*x >= min(Rn,Rn2):
    return 0.
  return dens*D/((Lneck+x)/ttA(Rn,x,Rdt)+ (Lneck2+x)/ttA(Rn2,x,Rdt)+ (Lpd-(Lneck+Lneck2+2.*x))/ttA(Rc,x,Rdt) )

def PunitAsym(x,Rn,Rn2,D,Lpd,Rdt,Rc,Lneck,Lneck2):
  # effective permeability without correction factor f_ih, assuming unit density (rho = 1 PD/um2), including desmotubule and Rc >= Rn
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  return 1e-3*Peff_naive_asym(x,Rn,Rn2,D,Lpd,1e-6,Rdt,Rc,Lneck,Lneck2)
  
  

def Peff_naive_cyl(x,Rn,D,Lpd,dens):
  # effective permeability without correction factor f_ih, straight cylindrical channel (used for calculating f_ih)
  if x >= Rn:
    return 0.
  return dens*D*ttA_cyl(Rn,x)/Lpd

def Peff(x,Rn,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck):
  # effective permeability including f_ih. f_ih =  Peff_MakTrap(..)/ Peff_naive_cyl(..)
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  return 1e-3* Peff_naive(x,Rn,D,Lpd,dens,Rdt,Rc,Lneck) *  Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,Rn,D,Lpd,dens)

def Peff_asym(x,Rn,Rn2,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Lneck2):
  # effective permeability including f_ih. f_ih =  Peff_MakTrap(..)/ Peff_naive_cyl(..)
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  # using average Rn to compute Fih, i.e., using that Fih is only weakly dependent on Rn. Assumes small differences in Rn only. 
  RnAvg = 0.5*(Rn+Rn2)
  return 1e-3* Peff_naive_asym(x,Rn,Rn2,D,Lpd,dens,Rdt,Rc,Lneck,Lneck2) *  Peff_MakTrap(x,RnAvg,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,RnAvg,D,Lpd,dens)

def Peff_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Rpit,nPit):
  # effective permeability including f_ih for pitfields. f_ih =  Peff_MakTrap_pitField(..)/ Peff_naive_cyl(..)
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  if nPit == 1:
    return Peff(x,Rn,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck)
  # Check whether Fih value is lower than with uniformly distributed PDs ( = no clustering). If not, computation invalid. 
  FihP =Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,nPit) /  Peff_naive_cyl(x,Rn,D,Lpd,dens)
  if FihP > nPit*Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,Rn,D,Lpd,dens):
    ## line below was for verification purposes: if Fih(pit) > Fih(single), then Rpit > 0.25*dPit (based on a triangular grid or square)
    print(Rpit, sqrt(1./dens*16./9.*sqrt(3.))/2. ,dens, sqrt(1/dens))
    ## return value -1 signals that the result of the computation is invalid
    return -1
  return 1e-3* Peff_naive(x,Rn,D,Lpd,dens,Rdt,Rc,Lneck) * FihP 

def Peff_pitField_asym(x,Rn,Rn2,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Lneck2,Rpit,nPit):
  # effective permeability including f_ih for pitfields. f_ih =  Peff_MakTrap_pitField(..)/ Peff_naive_cyl(..)
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  if nPit == 1:
    return Peff_asym(x,Rn,Rn2,D,Lpd,Lcell,dens,A,B,Rdt,Rc,Lneck,Lneck2)
  # Check whether Fih value is lower than with uniformly distributed PDs ( = no clustering). If not, computation invalid. 
  RnAvg = 0.5*(Rn+Rn2)
  FihP =Peff_MakTrap_pitField(x,RnAvg,D,Lpd,Lcell,dens,A,B,Rpit,nPit) /  Peff_naive_cyl(x,RnAvg,D,Lpd,dens)
  if FihP > nPit*Peff_MakTrap(x,RnAvg,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,RnAvg,D,Lpd,dens):
    ## line below was for verification purposes: if Fih(pit) > Fih(single), then Rpit > 0.25*dPit (based on a triangular grid or square)
    print(Rpit, sqrt(1./dens*16./9.*sqrt(3.))/2. ,dens, sqrt(1/dens))
    ## return value -1 signals that the result of the computation is invalid
    return -1
  return 1e-3* Peff_naive_asym(x,Rn,Rn2,D,Lpd,dens,Rdt,Rc,Lneck,Lneck2) * FihP 

def Fih_pitField(x,RnAvg,D,Lpd,Lcell,dens,A,B,Rpit,nPit):
    if nPit == 1:
      # from Peff_asym (equiv: Peff)
      return Peff_MakTrap(x,RnAvg,D,Lpd,Lcell,dens,A,B) / Peff_naive_cyl(x,RnAvg,D,Lpd,dens)
    #from Peff_pitField_asym above
    return Peff_MakTrap_pitField(x,RnAvg,D,Lpd,Lcell,dens,A,B,Rpit,nPit) /  Peff_naive_cyl(x,RnAvg,D,Lpd,dens)


def Rn_multiCyl(xMax,Rdt,nCyl):
  ## calculates Rn for sub-nano model. Minimum requirement: 1 nm spacers between cylindrical sub-nano channels. This requirement may cause Rn > Rdt + 2*xMax. 
  rNaive = Rdt + 2.*xMax
  r1 = xMax + (xMax + 0.5)/(sin(pi/nCyl))
  return max(rNaive,r1)

def  Peff_multiCyl(x,xMax,D,Lpd,Lcell,dens,A,B,Rdt,nCyl):
  ## effective permeability for sub-nano model including f_ih. f_ih =  Peff_MakTrap(..)/ Peff_naive_cyl(..)
  # The factor 1e-3 converts Peff to um/s (not nm/s)

  ## Simple model: use naively calculated Rn for calculating f_ih. 
  #Rn=Rdt+2.*xMax ## for simplicity, ignore overlap and reduced trapping efficiency due to the parts that separate the channels. 
  #return 1e-3* nCyl * Peff_naive_cyl(x,xMax,D,Lpd,dens) *  Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,Rn,D,Lpd,dens)

  ## Alternative model, somewhat more precise by using Rn based on 1 nm spacers.
  rNaive = Rdt + 2.*xMax
  r1 = xMax + (xMax + 0.5)/(sin(pi/nCyl)) 
  Rn = max(rNaive,r1)
  return 1e-3* nCyl * Peff_naive_cyl(x,xMax,D,Lpd,dens) *  (Peff_MakTrap_multiCyl(x,Rn,D,Lpd,Lcell,dens,A,B,xMax,nCyl) / (nCyl* Peff_naive_cyl(x,xMax,D,Lpd,dens)))

def Peff_multiCyl_pitField(x,xMax,D,Lpd,Lcell,dens,A,B,Rdt,nCyl,Rpit,nPit):
  # effective permeability including f_ih for sub-nano channel model in pitfields. f_ih =  Peff_MakTrap_pitField(..)/ Peff_naive_cyl(..)
  # f_ih is computed the same as for the default model in pit fields, but using 
  # The factor 1e-3 converts Peff to um/s (not nm/s)
  if nPit == 1:
    return Peff_multiCyl(x,xMax,D,Lpd,Lcell,dens,A,B,Rdt,nCyl)
  Rn = Rdt + 2.*xMax # naive computation; may not fit with 1 nm spacers. 
  r1 = xMax + (xMax + 0.5)/(sin(pi/nCyl)) 
  if r1 > Rn:
    Rpit = Rpit - Rn  + r1
    Rn = r1
  # Check whether Fih value is lower than with uniformly distributed PDs ( = no clustering). If not, computation invalid. 
  FihP =Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,nPit) / Peff_naive_cyl(x,Rn,D,Lpd,dens) 
  if FihP > nPit  *  (Peff_MakTrap_multiCyl(x,Rn,D,Lpd,Lcell,dens,A,B,xMax,nCyl) / (nCyl* Peff_naive_cyl(x,xMax,D,Lpd,dens))):
    # return value -1 signals that the result of the computation is invalid
    return -1
  return 1e-3* nCyl * Peff_naive_cyl(x,xMax,D,Lpd,dens) *  FihP
  #return 1e-3* nCyl * Peff_naive_cyl(x,xMax,D,Lpd,dens) *  Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,nPit) / Peff_naive_cyl(x,Rn,D,Lpd,dens)
  #return 1e-3* nCyl*nPit * Peff_naive_cyl(x,xMax,D,Lpd,dens) *  (Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,nPit) /( nPit*  Peff_naive_cyl(x,Rn,D,Lpd,dens)))

def compFih(x,D,Lpd,Lcell,dens,grid,Rdt):
  nCyl=9
  ab=[0,0]
  getAB(ab,grid)
  xTop = 50.
  xMax = 2.*x
  print("alpha_bar"+pars.colSep+"Fih_default"+pars.colSep+"Fih_sub-nano")
  while xMax < xTop:
    printFih(x,xMax,D,Lpd,Lcell,dens,ab[0],ab[1],Rdt,nCyl)
    xMax += 0.005
  return
    
def getDistDict():
  return {1:0,2:0.5, 3:sqrt(3.)/3,4:sqrt(2.)/2.,5:1.,6:2.*sqrt(3.)/3,7:1.,12:sqrt(13.)/3., 19:2.} # actual Pit Radius: dPit * distDict[x] + Rn

def getRpitDefault(distDict,Rn,dPit,nPit):
  return Rn + distDict[nPit]*dPit

def compFih_pitField(x,D,Lpd,Lcell,dens,grid,Rdt,pitList,dPit):
  #distDict={1:0,2:0.5, 3:sqrt(3.)/3,4:sqrt(2.)/2.,6:2.*sqrt(3.)/3,7:1.,12:sqrt(13.)/3., 19:2.} # actual Pit Radius: dPit * distDict[x] + Rn
  distDict=getDistDict()
  ab=[0,0]
  getAB(ab,grid)
  abRandom=[0,0]
  getAB(abRandom,"random")
  xTop = 50.
  xMax = 2.*x
  print("alpha_bar"+pars.colSep + pars.colSep.join(["Fih_" + str(s) for s in pitList]))
  while xMax < xTop:
    printFih_pitField(x,xMax,D,Lpd,Lcell,dens,ab[0],ab[1],Rdt,dPit,distDict,pitList,False,abRandom)
    xMax += 0.005
  return

def compFih_pitField2(x,D,Lpd,Lcell,grid,Rdt,pitList,dPit,Rn):
  #distDict={1:0,2:0.5, 3:sqrt(3.)/3,4:sqrt(2.)/2.,6:2.*sqrt(3.)/3,7:1.,12:sqrt(13.)/3., 19:2.} # actual Pit Radius: dPit * distDict[x] + Rn
  distDict=getDistDict()
  print("# #PDs/field"+pars.colSep+"PD_area/Pit_area"+pars.colSep+"pit_radius")
  for pj in pitList:
    print(pars.colSep.join(["#", str(pj), str(Rn**2*pj/(Rn+dPit*distDict[pj])**2), str(Rn+dPit*distDict[pj])]))
  ab=[0,0]
  getAB(ab,grid)
  abRandom=[0,0]
  getAB(abRandom,"random")
  dTop=3e-5
  xMax = 0.5*(Rn-Rdt)
  dens=0.1e-6
  print("density"+pars.colSep + pars.colSep.join(["Fih_" + str(s) for s in pitList]) + pars.colSep+"Fih_random")
  while dens < dTop:
    printFih_pitField(x,xMax,D,Lpd,Lcell,dens,ab[0],ab[1],Rdt,dPit,distDict,pitList,True,abRandom)
    dens += 0.1e-6
  return

def compTwinningEffect_pitField(x,D,Lpd,Lcell,grid,Rdt,pitList,dPit,Rn):
    ## computes the impact of increasing the numbers of PDs per pit through (repeated) twinning. 
  #distDict={1:0,2:0.5, 3:sqrt(3.)/3,4:sqrt(2.)/2.,6:2.*sqrt(3.)/3,7:1.,12:sqrt(13.)/3., 19:2.} # actual Pit Radius: dPit * distDict[x] + Rn
  distDict=getDistDict()
  print("# #PDs/field"+pars.colSep+"PD_area/Pit_area"+pars.colSep+"pit_radius")
  for pj in pitList:
    print(pars.colSep.join(["#", str(pj), str(Rn**2*pj/(Rn+dPit*distDict[pj])**2), str(Rn+dPit*distDict[pj])]))
  valid=[1 for pj in pitList]
  ab=[0,0]
  getAB(ab,grid)
  abRandom=[0,0]
  getAB(abRandom,"random")
  dTop=3e-5
  xMax = 0.5*(Rn-Rdt)
  dens=0.1e-6
  print("density"+pars.colSep + pars.colSep.join(["Fih_" + str(s) for s in pitList]) + pars.colSep+"Fih_random")
  while dens < dTop:
    printFih_pitField(x,xMax,D,Lpd,Lcell,dens,ab[0],ab[1],Rdt,dPit,distDict,pitList,True,abRandom,True,valid)
    dens += 0.1e-6
  return

def printFih_pitField(x,xMax,D,Lpd,Lcell,dens,A,B,Rdt,dPit,distDict,pitList,printDens=False,abRandom=[],twinning=False,valid=[]):
  Rn = Rdt + 2.*xMax
  if printDens:
    outList= [dens*1e6] ## convert to PDs/um2
  else:
    outList= [xMax]
  #for pj in pitList:
  for i in range(len(pitList)):
    pj = pitList[i]
    Rpit = Rn + distDict[pj]*dPit
    if pj == 1:
      outList.append( Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,Rn,D,Lpd,dens)  ) 
    else:
      if twinning:
        outList.append( Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens,A,B,Rpit,pj) /  (Peff_naive_cyl(x,Rn,D,Lpd,dens)  ) )
      else:
        outList.append( Peff_MakTrap_pitField(x,Rn,D,Lpd,Lcell,dens/pj,A,B,Rpit,pj) /  (pj*Peff_naive_cyl(x,Rn,D,Lpd,dens/pj)  ) )
      if valid:
        # assumes that pitList[1] == 1. 
        if not valid[i]:
          outList[-1]="x"
        else:
          if outList[-1] > pj* outList[1]:
            valid[i] = 0
            outList[-1]="x"
      ## functions below were used for testing only: They calculate Fih based on a single level of large clusters with radius Rpit and density the cluster density. From a theoretical perspective, this should provide a lower bound to the more elaborate formulas used above. It does. 
      ## Below: only consider the pit radius. This is used for testing agreement with theoretical results. The functions below should yield a lower bound to Fih as computed with the above, more elaborate method. Indeed, they yield somewhat lower values when densities are reasonably low. When densities get so high that actually, the pit fields are barely separated, the values get much lower. This is in the regime that the twinning calculation is no longer valid. 
      #if twinning:
        #outList.append(pj* Peff_MakTrap(x,Rpit,D,Lpd,Lcell,dens,A,B) /  (Peff_naive_cyl(x,Rpit,D,Lpd,dens)  ) )
      #else:
        #outList.append(Peff_MakTrap(x,Rpit,D,Lpd,Lcell,dens/pj,A,B) /  (Peff_naive_cyl(x,Rpit,D,Lpd,dens/pj)  ) )
  if abRandom:
    outList.append( Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,abRandom[0],abRandom[1]) /  Peff_naive_cyl(x,Rn,D,Lpd,dens)  ) 
      
  #print xMax, Peff_MakTrap(x,rNaive,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,rNaive,D,Lpd,dens),   (Peff_MakTrap_multiCyl(x,Rn,D,Lpd,Lcell,dens,A,B,xMax,nCyl) / (nCyl* Peff_naive_cyl(x,xMax,D,Lpd,dens)))
  print(pars.colSep.join([str(s) for s in outList]))
  return
  

def printFih(x,xMax,D,Lpd,Lcell,dens,A,B,Rdt,nCyl):
  ## prints xmax and correction factor Fih for default model (unobstructed sleeve) and sub-nano channel model with nCyl cylindrical channels per PD  
  rNaive = Rdt + 2.*xMax
  r1 = xMax + (xMax + 0.5)/(sin(pi/nCyl)) 
  Rn = max(rNaive,r1)
  print(xMax, Peff_MakTrap(x,rNaive,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,rNaive,D,Lpd,dens),   (Peff_MakTrap_multiCyl(x,Rn,D,Lpd,Lcell,dens,A,B,xMax,nCyl) / (nCyl* Peff_naive_cyl(x,xMax,D,Lpd,dens))))
  return

def getFih_sleeve(x,Rn,D,Lpd,Lcell,dens,A,B,Rdt):
  ## prints xmax and correction factor Fih for default model (unobstructed sleeve) 
  #rNaive = Rdt + 2.*xMax
  #return Peff_MakTrap(x,rNaive,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,rNaive,D,Lpd,dens)
  return Peff_MakTrap(x,Rn,D,Lpd,Lcell,dens,A,B) /  Peff_naive_cyl(x,Rn,D,Lpd,dens)
  
  

def getRn_multiCyl(xMax,Rdt,nCyl):
  # returns: naive Rn, Rn with (at least) 1 nm spacers, Rn with (at least) 0 nm spacers
  if xMax < 0:
    return [ "--", "--", "--"]
  rNaive = Rdt + 2.*xMax
  r1 = xMax + (xMax + 0.5)/(sin(pi/nCyl)) 
  r0 = xMax + (xMax + 0.)/(sin(pi/nCyl)) 
  return [rNaive, max(rNaive,r1), max(rNaive, r0)]

def linearInterpolate(yTarget,x1,x2,y1,y2):
  ## returns the x-value that best matches yTarget based on linear interpolation between (x,y) = (x1,y1) and (x2,y2)  
  return (yTarget - y1)*(x2-x1)/(y2-y1) + x1


def  computePeffValues(pars):
  fName="calculatedPeff."+pars.fileExt
  if pars.fileTag:
    fName="calculatedPeff_"+pars.fileTag+"."+pars.fileExt
  outfile=open(fName, 'w')
  sys.stdout = outfile
  print("## Calculated values of Peff, given all parameters. ")
  print("## Computations with computeClusterIncrease ==",pars.computeClusterIncrease, "(False: Rho = total PD density; True: Rho = PD cluster density, total PD density = Rho * n_cluster)")
  #if pars.doNotCombine and pars.asymmetricPDs:
  if pars.asymmetricPDs:
      print("# alpha"+pars.colSep+"Rho"+pars.colSep+"n_cluster"+pars.colSep+"alpha_bar"+pars.colSep+"Rn"+pars.colSep+"Rn2"+pars.colSep+"Rc"+pars.colSep+"Rdt"+pars.colSep+"Lpd"+pars.colSep+"Lneck"+pars.colSep+"Lneck2"+pars.colSep+"dPit"+pars.colSep+"Peff_default")
  else:
    if pars.compSubNano:
        print("# alpha"+pars.colSep+"Rho"+pars.colSep+"n_cluster"+pars.colSep+"alpha_bar"+pars.colSep+"Rn"+pars.colSep+"Rc"+pars.colSep+"Rdt"+pars.colSep+"Lpd"+pars.colSep+"Lneck"+pars.colSep+"dPit"+pars.colSep+"Peff_default"+pars.colSep+"Rn_subNano"+pars.colSep+"Peff_sub-nano")
    else:
      if pars.printRn:
        print("# alpha"+pars.colSep+"Rho"+pars.colSep+"n_cluster"+pars.colSep+"alpha_bar"+pars.colSep+"Rn"+pars.colSep+"Rc"+pars.colSep+"Rdt"+pars.colSep+"Lpd"+pars.colSep+"Lneck"+pars.colSep+"dPit"+pars.colSep+"Peff_default")
      else:
        print("# alpha"+pars.colSep+"Rho"+pars.colSep+"n_cluster"+pars.colSep+"alpha_bar"+pars.colSep+"Rc"+pars.colSep+"Rdt"+pars.colSep+"Lpd"+pars.colSep+"Lneck"+pars.colSep+"dPit"+pars.colSep+"Peff_default")
  
  #Rdt = pars.Rdt
  ab=[0,0]
  getAB(ab,pars.grid)
  distDict=getDistDict()
  if pars.doNotCombine:
    if pars.asymmetricPDs: 
      loopList=["LpdList","LneckList","Lneck2List","densList","RnList","Rn2List","RcList","RdtList","pitList","dPitList"]
      ## NOTE: uses Rn etc as input, not xMax. Confusing with symmetric case
    else:
      loopList=["LpdList","LneckList","densList","xMaxList","RcList","RdtList","pitList","dPitList"]
      ## NOTE: uses xMax as input, not Rn etc. Confusing with asymmetric case, but compatible with compSubNano
    class localPars(object):
      Lpd = 0
      Lneck = 0
      Lneck2 = 0
      dens = 0
      xMax = 0
      Rn = 0
      Rn2 = 0
      Rc = 0
      pit = 0
      dPit = 0
    #lenDict={x[0]:len(x[1]) for x in [["LpdList",pars.LpdList], ["densList",pars.densList],["xMaxList",pars.xMaxList],["RcList",pars.RcList],["pitList",pars.pitList],["dPitList",pars.dPitList]]}
    lenDict={x[0]:len(x[1]) for x in [[j,getattr(pars,j)] for j in loopList]}
    lMax = max(lenDict.items(),key=operator.itemgetter(1))[1]
    singleList = []
    multiList = []
    for i in lenDict.items():
      if i[1] == 1:
        singleList.append(i[0])
      elif i[1] == lMax:
        multiList.append(i[0])
      else:
        print ("Wrong number of values in ", i[0], ". Should be 1 or ", lMax)
        return
    for i in range(lMax):
      for j in loopList:
        if j in singleList:
          setattr(localPars,j[:-4],getattr(pars,j)[0])
        if j in multiList:
          setattr(localPars,j[:-4],getattr(pars,j)[i])
      if pars.asymmetricPDs:
        xMax = 0.5*(min(localPars.Rn,localPars.Rn2)-localPars.Rdt)
        if xMax > pars.x:
          Rpit = 0.5*(localPars.Rn+localPars.Rn2) + distDict[localPars.pit]*localPars.dPit
          if pars.computeClusterIncrease:
            pe = Peff_pitField_asym(pars.x,localPars.Rn,localPars.Rn2,pars.diff,localPars.Lpd,pars.Lcell,localPars.dens,ab[0],ab[1],localPars.Rdt,localPars.Rc,localPars.Lneck,localPars.Lneck2,Rpit,localPars.pit)
          else: 
            pe = Peff_pitField_asym(pars.x,localPars.Rn,localPars.Rn2,pars.diff,localPars.Lpd,pars.Lcell,localPars.dens/localPars.pit,ab[0],ab[1],localPars.Rdt,localPars.Rc,localPars.Lneck,localPars.Lneck2,Rpit,localPars.pit)
        else:
          pe = 0
        print(pars.colSep.join(str(s) for s in [pars.x, localPars.dens*1e6, localPars.pit, xMax, localPars.Rn, localPars.Rn2, localPars.Rc, localPars.Rdt, localPars.Lpd, localPars.Lneck, localPars.Lneck2, localPars.dPit, pe]))
      else:
        Rn = localPars.xMax * 2.+ localPars.Rdt     
        Rpit = Rn + distDict[localPars.pit]*localPars.dPit
        if pars.computeClusterIncrease:
          pe = Peff_pitField(pars.x,Rn,pars.diff,localPars.Lpd,pars.Lcell,localPars.dens,ab[0],ab[1],localPars.Rdt,localPars.Rc,localPars.Lneck,Rpit,localPars.pit)
        else: 
          pe = Peff_pitField(pars.x,Rn,pars.diff,localPars.Lpd,pars.Lcell,localPars.dens/localPars.pit,ab[0],ab[1],localPars.Rdt,localPars.Rc,localPars.Lneck,Rpit,localPars.pit)
        if pars.compSubNano:
          ## NOTE: computed for straight channel; Rc does not take any effect here. 
          if pars.computeClusterIncrease:
            peSN =  Peff_multiCyl_pitField(pars.x,localPars.xMax,pars.diff,localPars.Lpd,pars.Lcell,localPars.dens,ab[0],ab[1],localPars.Rdt,9,Rpit,localPars.pit)
          else: 
            peSN =  Peff_multiCyl_pitField(pars.x,localPars.xMax,pars.diff,localPars.Lpd,pars.Lcell,localPars.dens/localPars.pit,ab[0],ab[1],localPars.Rdt,9,Rpit,localPars.pit)
          print(pars.colSep.join(str(s) for s in [pars.x, localPars.dens*1e6, localPars.pit, localPars.xMax, Rn, localPars.Rc, localPars.Rdt, localPars.Lpd, localPars.Lneck, localPars.dPit, pe,Rn_multiCyl(xMax,localPars.Rdt,9), peSN]))
        else:
          if pars.printRn:
            print(pars.colSep.join(str(s) for s in [pars.x, localPars.dens*1e6, localPars.pit, localPars.xMax, Rn, localPars.Rc, localPars.Rdt, localPars.Lpd, localPars.Lneck, localPars.dPit, pe]))
          else:
            print(pars.colSep.join(str(s) for s in [pars.x, localPars.dens*1e6, localPars.pit, localPars.xMax, localPars.Rc, localPars.Rdt, localPars.Lpd, localPars.Lneck, localPars.dPit, pe]))
  else:
    if pars.compSubNano:
      for Lpd in pars.LpdList:
        for Lneck in pars.LneckList:
          for dd in pars.densList:
            for xMax in pars.xMaxList:
              gotit = False
              for Rdt in pars.RdtList:
                for Rc in pars.RcList: ## TODO
                  Rn = xMax *2.+Rdt
                  if Rc < Rn:
                    if gotit:
                      continue
                    gotit=True
                    Rc = Rn
                  for nPit in pars.pitList:
                    for dPit in pars.dPitList:
                      Rpit = Rn + distDict[nPit]*dPit
                      if pars.computeClusterIncrease:
                        pe = Peff_pitField(pars.x,Rn,pars.diff,Lpd,pars.Lcell,dd,ab[0],ab[1],Rdt,Rc,Lneck,Rpit,nPit)
                      else: 
                        pe = Peff_pitField(pars.x,Rn,pars.diff,Lpd,pars.Lcell,dd/nPit,ab[0],ab[1],Rdt,Rc,Lneck,Rpit,nPit)
                      if pars.compSubNano:
                        ## NOTE: computed for straight channel; Rc does not take any effect here. 
                        if pars.computeClusterIncrease:
                          peSN =  Peff_multiCyl_pitField(pars.x,xMax,pars.diff,Lpd,pars.Lcell,dd,ab[0],ab[1],Rdt,9,Rpit,nPit)
                        else: 
                          peSN =  Peff_multiCyl_pitField(pars.x,xMax,pars.diff,Lpd,pars.Lcell,dd/nPit,ab[0],ab[1],Rdt,9,Rpit,nPit)
                        print(pars.colSep.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rn, Rc, Rdt, Lpd, Lneck, dPit, pe,Rn_multiCyl(xMax,Rdt,9), peSN]))
                      else:
                        if pars.printRn:
                          print(pars.colSep.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rn, Rc, Rdt, Lpd, Lneck, dPit, pe]))
                        else:
                          print(pars.colSep.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rc, Rdt, Lpd, Lneck, dPit, pe]))
                gotit = False
          print("")
    else:
      for Lpd in pars.LpdList:
        for Lneck in pars.LneckList:
          for dd in pars.densList:
            for Rn in pars.RnList:
              gotit = False
              for Rdt in pars.RdtList:
                xMax = 0.5*(Rn-Rdt)
                for Rc in pars.RcList: 
                  if Rc < Rn:
                    if gotit:
                      continue
                    gotit=True
                    Rc = Rn
                  for nPit in pars.pitList:
                    for dPit in pars.dPitList:
                      if pars.asymmetricPDs:
                        for Rn2 in pars.Rn2List:
                          for Lneck2 in pars.Lneck2List:
                            Rpit = 0.5*(Rn+Rn2) + distDict[nPit]*dPit
                            if pars.computeClusterIncrease:
                              pe = Peff_pitField_asym(pars.x,Rn,Rn2,pars.diff,Lpd,pars.Lcell,dd,ab[0],ab[1],Rdt,Rc,Lneck,Lneck2,Rpit,nPit)
                            else: 
                              pe = Peff_pitField_asym(pars.x,Rn,Rn2,pars.diff,Lpd,pars.Lcell,dd/nPit,ab[0],ab[1],Rdt,Rc,Lneck,Lneck2,Rpit,nPit)
                            print(pars.colSep.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rn, Rn2, Rc, Rdt, Lpd, Lneck,Lneck2, dPit, pe]))
                      else:
                        Rpit = Rn + distDict[nPit]*dPit
                        if pars.computeClusterIncrease:
                          pe = Peff_pitField(pars.x,Rn,pars.diff,Lpd,pars.Lcell,dd,ab[0],ab[1],Rdt,Rc,Lneck,Rpit,nPit)
                        else: 
                          pe = Peff_pitField(pars.x,Rn,pars.diff,Lpd,pars.Lcell,dd/nPit,ab[0],ab[1],Rdt,Rc,Lneck,Rpit,nPit)
                        if pars.printRn:
                          print(pars.colSep.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rn, Rc, Rdt, Lpd, Lneck, dPit, pe]))
                        else:
                          print(pars.colSep.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rc, Rdt, Lpd, Lneck, dPit, pe]))
                gotit = False
          print("")

  outfile.close()
  sys.stdout = sys.__stdout__
  return

def  computePunitValues(pars):
  fName="calculatedPunit."+pars.fileExt
  if pars.fileTag:
    fName="calculatedPunit_"+pars.fileTag+"."+pars.fileExt
  outfile=open(fName, 'w')
  sys.stdout = outfile
  print("## Calculated values of Punit, given Rn, Rc, Rdt, Lpd, Lneck. ")
  if pars.asymmetricPDs:
    print("# alpha"+pars.colSep+"alpha_bar"+pars.colSep+"Rn"+pars.colSep+"Rn2"+pars.colSep+"Rc"+pars.colSep+"Rdt"+pars.colSep+"Lpd"+pars.colSep+"Lneck"+pars.colSep+"Lneck2"+pars.colSep+"Punit_default")
  else:
    if pars.compSubNano:
        print ("# Not implemented")
        outfile.close()
        sys.stdout = sys.__stdout__
        return
    else:
        print("# alpha"+pars.colSep+"alpha_bar"+pars.colSep+"Rn"+pars.colSep+"Rc"+pars.colSep+"Rdt"+pars.colSep+"Lpd"+pars.colSep+"Lneck"+pars.colSep+"Punit_default")
  if pars.doNotCombine:
    if pars.asymmetricPDs:
      loopList=["LpdList","RdtList","RnList","Rn2List","RcList","LneckList","Lneck2List"]
    else:
      loopList=["LpdList","RdtList","RnList","RcList","LneckList"]
    class localPars(object):
      Lpd = 0
      Rdt = 0
      Rn = 0
      Rn2 = 0
      Rc = 0
      Lneck= 0
      Lneck2 = 0
    lenDict={x[0]:len(x[1]) for x in [[j,getattr(pars,j)] for j in loopList]}
    lMax = max(lenDict.items(),key=operator.itemgetter(1))[1]
    singleList = []
    multiList = []
    for i in lenDict.items():
      if i[1] == 1:
        singleList.append(i[0])
      elif i[1] == lMax:
        multiList.append(i[0])
      else:
        print ("Wrong number of values in ", i[0], ". Should be 1 or ", lMax)
        return
    for i in range(lMax):
      for j in loopList:
        if j in singleList:
          setattr(localPars,j[:-4],getattr(pars,j)[0])
        if j in multiList:
          setattr(localPars,j[:-4],getattr(pars,j)[i])
      if localPars.Rc < localPars.Rn:
          continue
      if pars.asymmetricPDs:
        xMax = 0.5*(min(localPars.Rn,localPars.Rn2)-localPars.Rdt)
        if xMax > pars.x:
          pe = PunitAsym(pars.x,localPars.Rn,localPars.Rn2,pars.diff,localPars.Lpd,localPars.Rdt,localPars.Rc,localPars.Lneck,localPars.Lneck2)
        else:
          pe = 0
        if pars.compSubNano:
          continue
        else:
          print(pars.colSep.join(str(s) for s in [pars.x,xMax, localPars.Rn, localPars.Rn2, localPars.Rc, localPars.Rdt, localPars.Lpd, localPars.Lneck, localPars.Lneck2, pe]))
      else:
        xMax = 0.5*(localPars.Rn-localPars.Rdt)
        if xMax > pars.x:
          pe = Punit(pars.x,localPars.Rn,pars.diff,localPars.Lpd,localPars.Rdt,localPars.Rc,localPars.Lneck)
        else:
          pe = 0
        if pars.compSubNano:
          continue
        else:
          print(pars.colSep.join(str(s) for s in [pars.x,xMax, localPars.Rn, localPars.Rc, localPars.Rdt, localPars.Lpd, localPars.Lneck, pe]))
    print("")

  else:
    for Lpd in pars.LpdList:
      for Rdt in pars.RdtList:
        for Rn in pars.RnList:
          gotit = False
          xMax = 0.5*(Rn-Rdt)
          for Rc in pars.RcList: 
            #Rn = xMax *2.+Rdt
            if Rc < Rn:
              if gotit:
                continue
              gotit=True
              Rc = Rn
            for Lneck in pars.LneckList:
              if pars.asymmetricPDs:
                for Rn2 in pars.Rn2List:
                  for Lneck2 in pars.Lneck2List:
                    pe = PunitAsym(pars.x,Rn,Rn2,pars.diff,Lpd,Rdt,Rc,Lneck,Lneck2)
                    if pars.compSubNano:
                      continue
                    else:
                      xMax = 0.5*(min(Rn,Rn2)-Rdt)
                      print(pars.colSep.join(str(s) for s in [pars.x,xMax, Rn, Rn2, Rc, Rdt, Lpd, Lneck,Lneck2,pe]))
              else:
                pe = Punit(pars.x,Rn,pars.diff,Lpd,Rdt,Rc,Lneck)
                if pars.compSubNano:
                  ## NOTE: computed for straight channel; Rc does not take any effect here. 
                  # peSN =  Peff_multiCyl_pitField(pars.x,xMax,pars.diff,Lpd,pars.Lcell,dd,ab[0],ab[1],Rdt,9,Rpit,nPit)
                      #print('\t'.join(str(s) for s in [pars.x, dd*1e6, nPit, xMax, Rn, Rc, Lpd, dPit, pe,Rn_multiCyl(xMax,Rdt,9), peSN]))
                  continue
                else:
                  print(pars.colSep.join(str(s) for s in [pars.x,xMax, Rn, Rc, Rdt, Lpd, Lneck,pe]))
      print("")
  outfile.close()
  sys.stdout = sys.__stdout__
  return

def  printTargetDensitiesStraightPD(Rdt,x,xMaxList,peList,diff,Lpd,Lcell,grid,Lneck,dInc,cfDict,colSep,compact=False,printRn=False):
  print("## Required PD densities for target Peff, given alpha_bar")
  if cfDict:
    if printRn:
      print("alpha"+colSep+"alpha_bar"+colSep+"Rn"+colSep+"Peff"+colSep+"Rho_default"+colSep+"deviation_from_Peff"+colSep+"Rho_sub-nano"+colSep+"Rho_sub-nano_neckOnly"+colSep+"Rho_spoke_entrance")
    else:
      print("alpha"+colSep+"alpha_bar"+colSep+"Peff"+colSep+"Rho_default"+colSep+"deviation_from_Peff"+colSep+"Rho_sub-nano"+colSep+"Rho_sub-nano_neckOnly"+colSep+"Rho_spoke_entrance")
  else:
    if printRn:
      print("alpha"+colSep+"alpha_bar"+colSep+"Rn"+colSep+"Peff"+colSep+"Rho_default"+colSep+"deviation_from_Peff")
    else:
      print("alpha"+colSep+"alpha_bar"+colSep+"Peff"+colSep+"Rho_default"+colSep+"deviation_from_Peff")

  ab=[0,0]
  getAB(ab,grid)
  dStart = dInc
  for xMax in xMaxList:
    dd = dStart
    pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
    peTargetOld = peList[0]
    for peTarget in peList:
      if peTarget < peTargetOld:
        dd = dStart
        pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
      peTargetOld = peTarget
      while pe < peTarget:
        ddPrev = dd          
        dd += dInc 
        pePrev = pe
        pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
      ddBest=linearInterpolate(peTarget,ddPrev, dd,pePrev,pe) 
      peBest=Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,ddBest,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
      if cfDict:
        if printRn:
          print(colSep.join(str(s) for s in [x, xMax, 2.*xMax +Rdt, peTarget, ddBest*1e6,  peBest-peTarget]+[ddBest*i  *1e6 for i in cfDict[xMax]]))
        else:
          print(colSep.join(str(s) for s in [x, xMax, peTarget, ddBest*1e6,  peBest-peTarget]+[ddBest*i  *1e6 for i in cfDict[xMax]]))
      else:
        if printRn:
          print(colSep.join(str(s) for s in [x, xMax, 2.*xMax+Rdt, peTarget, ddBest*1e6,  peBest-peTarget]))
        else:
          print(colSep.join(str(s) for s in [x, xMax, peTarget, ddBest*1e6,  peBest-peTarget]))
    if not compact:
      print("")
  return

def  printTargetAperturesStraightPD(Rdt,x,densList,peList,diff,Lpd,Lcell,grid,Lneck,xInc,compSubNano,colSep):
  print("## Required alpha_bar and corresponding Rn for target Peff, given PD density")
  if compSubNano:
    print("alpha"+colSep+"Rho"+colSep+"Peff"+colSep+"alpha_bar_default"+colSep+"Rn_default"+colSep+"deviation_from_Peff"+colSep+"alpha_bar_sub-nano"+colSep+"Rn_sub-nano_naive"+colSep+"Rn_sub-nano_spacers1nm"+colSep+"Rn_sub-nano_spacers0nm") 
  else:
    print("alpha"+colSep+"Rho"+colSep+"Peff"+colSep+"alpha_bar_default"+colSep+"Rn_default"+colSep+"deviation_from_Peff")
  ab=[0,0]
  getAB(ab,grid)
  xStart = x+xInc
  for dd in densList:
    xMax=xStart
    pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
    peTargetOld = peList[0]
    for peTarget in peList:
      if peTarget < peTargetOld:
        xMax=xStart
        pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
      peTargetOld = peTarget
      while pe < peTarget:
        xMaxPrev=xMax
        xMax += xInc 
        pePrev = pe
        pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
      xMaxBest=linearInterpolate(peTarget,xMaxPrev, xMax,pePrev,pe) 
      peBest=Peff(x,xMaxBest*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMaxBest*2.+Rdt,Lneck)
      outList = [ x, dd*1e6,  peTarget, xMaxBest,Rdt+2*xMaxBest, peBest -peTarget]
      if compSubNano:
        xSN = max(xMax -7,xStart)
        peSN =  Peff_multiCyl(x,xSN,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,9)
        print(peSN,peTarget)
        while peSN < peTarget:
          peSNprev = peSN
          xSNprev = xSN
          xSN += xInc 
          peSN =  Peff_multiCyl(x,xSN,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,9)
        xSNbest=linearInterpolate(peTarget,xSNprev, xSN,peSNprev,peSN) 
        outList.append(xSNbest)
        outList += getRn_multiCyl(xSNbest,Rdt,9)
      print(colSep.join([str(s) for s in outList]))
    print("")
  return 

def  printTargetAperturesStraightPD_twinning(Rdt,x,densList,peList,diff,Lpd,Lcell,grid,Lneck,xInc,twinningList,dPit,compSubNano,colSep):
  print("## Required alpha_bar and corresponding Rn for target Peff, given PD density")
  if compSubNano:
    print("alpha"+colSep+"Rho"+colSep+"n_cluster"+colSep+"Peff"+colSep+"alpha_bar_default"+colSep+"Rn_default"+colSep+"deviation_from_Peff"+colSep+"alpha_bar_sub-nano"+colSep+"Rn_sub-nano_naive"+colSep+"Rn_sub-nano_spacers1nm"+colSep+"Rn_sub-nano_spacers0nm") 
  else:
    print("alpha"+colSep+"Rho"+colSep+"n_cluster"+colSep+"Peff"+colSep+"alpha_bar_default"+colSep+"Rn_default"+colSep+"deviation_from_Peff")
  ab=[0,0]
  getAB(ab,grid)
  distDict=getDistDict()
  xStart = x+xInc
  for dd in densList:
    xMax=xStart
    pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
    for peTarget in peList:
        for pj in twinningList:
          xMax=xStart
          if pj > 1: ## also compute values with homogeneous density increase
            pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd*pj,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
            while pe < peTarget:
              xMaxPrev=xMax
              xMax += xInc 
              pePrev = pe
              pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd*pj,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
            xMaxBest=linearInterpolate(peTarget,xMaxPrev, xMax,pePrev,pe) 
            peBest=Peff(x,xMaxBest*2.+Rdt,diff,Lpd,Lcell,dd*pj,ab[0],ab[1],Rdt,xMaxBest*2.+Rdt,Lneck)
            outList = [ x, dd*1e6*pj, 1, peTarget, xMaxBest,Rdt+2*xMaxBest, peBest -peTarget]
            if compSubNano:
              xSN = max(xMax -1,xStart)
              peSN =  Peff_multiCyl(x,xSN,diff,Lpd,Lcell,dd*pj,ab[0],ab[1],Rdt,9)
              while peSN < peTarget:
                peSNprev = peSN
                xSNprev = xSN
                xSN += xInc 
                peSN =  Peff_multiCyl(x,xSN,diff,Lpd,Lcell,dd*pj,ab[0],ab[1],Rdt,9)
              xSNbest=linearInterpolate(peTarget,xSNprev, xSN,peSNprev,peSN) 
              outList.append(xSNbest)
              outList += getRn_multiCyl(xSNbest,Rdt,9)
            print(colSep.join([str(s) for s in outList]))

          Rpit = xMax*2.+Rdt + distDict[pj]*dPit
          xMax = max(xMax -1,xStart)
          pe = Peff_pitField(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck,Rpit,pj)
          #pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
          while pe < peTarget:
            xMaxPrev=xMax
            xMax += xInc 
            pePrev = pe
            pe = Peff_pitField(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck,Rpit,pj)
            #pe =  Peff(x,xMax*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMax*2.+Rdt,Lneck)
            if pe < 0: 
              print("Computation invalid:")
              break
          if pe < 0:
            outList = [ x, dd*1e6*pj, pj, peTarget, "--","--","--"]
          else:
            xMaxBest=linearInterpolate(peTarget,xMaxPrev, xMax,pePrev,pe) 
            #peBest=Peff(x,xMaxBest*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMaxBest*2.+Rdt,Lneck)
            peBest=Peff_pitField(x,xMaxBest*2.+Rdt,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,xMaxBest*2.+Rdt,Lneck,Rpit,pj)
            outList = [ x, dd*1e6*pj, pj, peTarget, xMaxBest,Rdt+2*xMaxBest, peBest -peTarget]
          if compSubNano:
            #xSN = max(xMax -1,xStart)
            xSN = xStart
            peSN =  Peff_multiCyl_pitField(x,xSN,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,9,Rpit,pj)
            while peSN < peTarget:
              peSNprev = peSN
              xSNprev = xSN
              xSN += xInc 
              peSN =  Peff_multiCyl_pitField(x,xSN,diff,Lpd,Lcell,dd,ab[0],ab[1],Rdt,9,Rpit,pj)
              if peSN < 0: 
                break
            if peSN < 0:
              outList.append ("--")
              xSNbest = -1
            else:
              xSNbest=linearInterpolate(peTarget,xSNprev, xSN,peSNprev,peSN) 
              outList.append(xSNbest)
            outList += getRn_multiCyl(xSNbest,Rdt,9)
          print(colSep.join([str(s) for s in outList]))
    print("")
  return 

if __name__ == '__main__':
  try:
    parfile = sys.argv[1]
  except:
    parfile = "parameters.txt"
  main(parfile)
