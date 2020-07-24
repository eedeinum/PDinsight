#!/usr/bin/python3

import tkinter as tk
#from tkinter import *
from tkinter import messagebox
from tkinter.ttk import *
import PDinsight as pd

class PDgui(tk.Frame):

  def fill(self):
    # choose mode. 
    self.rModeLabel = tk.Label(self,text="Choose mode (load default parameters to change)")
    self.rModeLabel.grid(row=0,column=0,sticky="W",padx=5,pady=1)
    self.rMode = tk.ttk.Combobox(self)
    self.rMode["values"] = ( "computeDens", "computeAperture", "computeVals","computeUnitVals","sensitivityAnalysis", "computeRnDensityGraph")
    self.rMode.current(3)
    self.rMode.grid(row=0,column=1,sticky="W",padx=5,pady=1)
    ## later get value using self.rMode.get()
    self.loadButton = tk.Button(self,text="Load default parameters", command=self.loadDefaults)
    self.loadButton.grid(row=0,column=2,sticky="W",padx=5,pady=1)
    

    # finally: buttons for running & just exporting parameters
    self.parLabel = tk.Label(self,text="Parameter file name: ")
    self.parFileBox = tk.Entry(self,textvariable=self.parFileName)
    self.runButton = tk.Button(self,text="Run PDinsight", command=self.runClicked,state=tk.DISABLED)
    self.expButton = tk.Button(self,text="Export parameter file",command=self.parClicked,state=tk.DISABLED)
    self.infoButton = tk.Button(self,text="Info",command=self.infoClicked)
    self.quitButton = tk.Button(self,text="Quit",command=self.quit)
    self.parLabel.grid(row=1,column=0,sticky="W",columnspan=1,padx=5,pady=1)
    self.parFileBox.grid(row=1,column=1,sticky="W",columnspan=1,padx=5,pady=1)
    self.runButton.grid(row=2,column=0,sticky="W",columnspan=1,padx=5,pady=1)
    self.expButton.grid(row=2,column=1,sticky="W",columnspan=1,padx=5,pady=1)
    self.infoButton.grid(row=2,column=2,sticky="W",columnspan=1,padx=5,pady=1)
    self.quitButton.grid(row=2,column=3,sticky="E",columnspan=1,padx=5,pady=1)
    self.optLabel= tk.Label(self,text="Available options:") 
    self.optLabel.grid(row=3,column=1,sticky="W",columnspan=2,padx=5,pady=1)

  def loadDefaults(self):
    self.runType = self.rMode.get()
    print(self.runType)
    self.getOptions(True)
    self.runButton.config(state=tk.NORMAL)
    self.expButton.config(state=tk.NORMAL)


  def getOptions(self,useDefault=False):
    self.needPars = set(self.defaultPars.requiredPars[self.runType] + self.defaultPars.requiredPars["all"])
    for m in self.rMode["values"]:  ## first set all modes to false to avoid conflicts
      setattr(self.newPars,m,False)
    setattr(self.newPars,self.runType,True)
    self.options = []
    for i in self.needPars:
      if i in self.defaultPars.optionList:
        self.options.append(i)
        setattr(self.newPars,i,getattr(self.defaultPars,i))
    for i in self.options: # removal must appear in a separate step, or fails
      self.needPars.remove(i)
    self.radioList={i:[tk.IntVar(),getattr(self.defaultPars,i)] for i in self.options}
    self.radioList["outputType"] = [tk.IntVar(),self.defaultPars.outputOptions.index(getattr(self.defaultPars,"outputType"))]
    for i in self.options+["outputType"]  :
      self.radioList[i][0].set(self.radioList[i][1])
    try:
      for r in self.radioButtons:
        r.destroy()
    except:
      pass
    nrow=4
    self.radioButtons = [tk.Label(self, text="outputType", justify = tk.LEFT, padx = 20)]
    self.radioButtons[-1].grid(row=nrow,column=0,sticky="W",columnspan=1,padx=5,pady=1) 
    for i in range(len(self.defaultPars.outputOptions)):
      self.radioButtons.append(tk.Radiobutton(self,text=self.defaultPars.outputOptions[i],variable=self.radioList["outputType"][0],value=i,command=self.updatePars))
      self.radioButtons[-1].grid(row=nrow,column=1,sticky="W",columnspan=1,padx=5+65*i,pady=1) 
    nrow=5
    for r in self.radioList:
      if r == "outputType":
        continue
      self.radioButtons.append(tk.Label(self, text=r, justify = tk.LEFT, padx = 20))
      self.radioButtons[-1].grid(row=nrow,column=0,sticky="W",columnspan=1,padx=5,pady=1) 
      #self.radioButtons.append(tk.Radiobutton(self,text="yes",variable=self.radioList[r][0],value=1,command=lambda: self.updatePars([r]))) ## issue: r is evaluated at button press
      self.radioButtons.append(tk.Radiobutton(self,text="yes",variable=self.radioList[r][0],value=1,command=self.updatePars))
      self.radioButtons[-1].grid(row=nrow,column=1,sticky="W",columnspan=1,padx=5,pady=1) 
      #self.radioButtons.append(tk.Radiobutton(self,text="no",variable=self.radioList[r][0],value=0,command=lambda: self.updatePars([r]))) ## issue: r is evaluated at button press
      self.radioButtons.append(tk.Radiobutton(self,text="no",variable=self.radioList[r][0],value=0,command=self.updatePars)) 
      self.radioButtons[-1].grid(row=nrow,column=1,sticky="W",columnspan=1,padx=70,pady=1) 
      self.radioButtons.append(tk.Label(self,text=self.defaultPars.optionDescriptions[r]))
      self.radioButtons[-1].grid(row=nrow,column=1,sticky="W",columnspan=3,padx=130,pady=1) 
      nrow +=1
    self.fixPars() # conditionally alters needPars based on options 
    self.drawPars(useDefault)

  def drawPars(self,useDefault=False):
      ## TODO: make fancy descriptions; below is a silly temporary solution
      #self.descriptions = {p:p for p in self.needPars}
      self.descriptions = {p:self.defaultPars.parDescriptions[p] for p in self.needPars}
      self.extraInfo = {p:self.defaultPars.parExtraInfo[p] for p in self.needPars}
      try:
          for p in self.parBoxes:
            self.parBoxes[p][0].destroy()
            if len(self.parBoxes[p]) == 4:
                self.parBoxes[p][3].destroy()
            else:
                self.parBoxes[p][2].destroy()
          self.parLabel.destroy()
          self.parLabel2.destroy()
      except Exception as err:
          #print (err)
          pass
      self.parLabel = tk.Label(self,text="Required parameters:")
      self.parLabel.grid(row=10,column=1,sticky="W",columnspan=1,padx=5,pady=1)
      self.parLabel2 = tk.Label(self,text="Multiple values allowed if name ends with \"List\"")
      self.parLabel2.grid(row=10,column=2,sticky="W",columnspan=2,padx=5,pady=1)
      self.parBoxes = {}
      nrow=11
      self.parLists = {c:[] for c in self.parCats}
      for p in self.needPars:
        found =0
        for c in self.parCats[:-1]:
          if p in self.fullParLists[c]:
            self.parLists[c].append(p)
            found =1
        if found == 0:
          self.parLists["Misc"].append(p)
      for c in self.parCats:
        self.parBoxes[c] = [tk.Label(self,text=c+" parameters",justify=tk.LEFT),"",tk.Label(self,text="dummy")]
        self.parBoxes[c][0].grid(row=nrow,column=0,sticky="W",columnspan=1,padx=5,pady=1)
        nrow+=1
        for p in self.parLists[c]:
          #vp = (self.register(lambda: self.getValidateCommand(p)), '%P')  
          vp = self.getValidateCommand(p)  
          #vpr = (self.register(vp), '%P','%s')
          vpr = (self.register(vp), '%P')
          if self.extraInfo[p]:
            self.parBoxes[p]=[tk.Label(self,text=self.descriptions[p],justify=tk.LEFT,padx=20,cursor="question_arrow"),tk.StringVar(self,self.formatPar(p,useDefault)),self.extraInfo[p]] 
            k=1
          else:
            self.parBoxes[p]=[tk.Label(self,text=self.descriptions[p],justify=tk.LEFT,padx=20),tk.StringVar(self,self.formatPar(p,useDefault))] 
            k=0
          self.parBoxes[p].append(tk.Entry(self,textvariable=self.parBoxes[p][1],justify=tk.LEFT,width=100,validatecommand=vpr,validate="all"))
          self.parBoxes[p][0].grid(row=nrow,column=0,sticky="W",columnspan=1,padx=5,pady=1)
          self.parBoxes[p][2+k].grid(row=nrow,column=1,sticky="W",columnspan=4,padx=5,pady=1)
          if self.extraInfo[p]:
            self.parBoxes[p][0].bind("<Button-1>", lambda event, p=p: messagebox.showinfo(p,self.parBoxes[p][2]))
          nrow+=1

        # valid percent substitutions (from the Tk entry man page)
        # note: you only have to register the ones you need; 
        #
        # %d = Type of action (1=insert, 0=delete, -1 for others)
        # %i = index of char string to be inserted/deleted, or -1
        # %P = value of the entry if the edit is allowed
        # %s = value of entry prior to editing
        # %S = the text string being inserted or deleted, if any
        # %v = the type of validation that is currently set
        # %V = the type of validation that triggered the callback
        #      (key, focusin, focusout, forced)
        # %W = the tk name of the widget

  def getValidateCommand(self,p):
    if p in ["pitList","twinningList"]:
      return self.verifyPit
    if p == "grid":
      return self.verifyGrid
    if p == "gridList":
      return self.verifyListGrid
    if p in self.defaultPars.simpleFloats + list(self.defaultPars.factorFloats.keys()):
      return self.verifyFloat
    if p in  self.defaultPars.listFloats + list(self.defaultPars.listFactorFloats.keys()):
      return self.verifyListFloat
    if p in self.defaultPars.simpleStrings:
      return self.verifyString
    ## Below option does not occur (yet), as the only listStrings is "gridList"
    #if p in self.defaultPars.listStrings:
      #return self.verifyListString

  def testPartial(self,tList,t):
    L=len(t)
    for tt in tList:
      try:
        if t == tt[:L]:
          return True
      except:
        pass
    return False

  def testPartialFloat(self,t):
    endFloat = ["E","e","-","+"]
    if t[-1] in endFloat:
      return True
    endFloat2 = ["e-", "e+","E-", "E+"]
    if t[-2:] in endFloat2:
      return True
    return False

  def verifyPit(self,P): ## actually, this is a list!
    pitOK = [1,2,3,4,5,6,7,12,19]
    sp = P.split()
    for i in sp:
      try:
        if int(i) not in pitOK: 
          print ("Invalid entry: ", P, ", failing at: ", i)
          print ("Valid options: ", ' '.join([str(p) for p in pitOK]))
          return False
      except ValueError:
        print ("Invalid entry: ", P, ", failing at: ", i)
        print ("Valid options: ", ' '.join([str(p) for p in pitOK]))
        return False
    return True

  def verifyGrid(self,P): 
    gridOK = ["triangular", "square", "hexagonal", "hex", "random"]
    if not P:
      return True
    if P not in gridOK: 
      if self.testPartial(gridOK,P):
        return True
      else:
        print ("Invalid entry: ", P)
        print ("Valid options: ", ' '.join(gridOK))
        return False
    return True

  def verifyListGrid(self,P): 
    gridOK = ["triangular", "square", "hexagonal", "hex", "random"]
    sp = P.split()
    for i in sp:
      if i not in gridOK: 
        if self.testPartial(gridOK,P):
          continue
        else:
          print ("Invalid entry: ", P, ", failing at: ", i)
          print ("Valid options: ", ' '.join(gridOK))
          return False
    return True

  def verifyFloat(self,P):
    if not P:
      return True
    try:
      f = float(P)
    except ValueError:
      if not self.testPartialFloat(P):
        print ("Invalid entry: ", P)
        print ("Valid entries: floating point number")
        return False
    return True

  def verifyListFloat(self,P):
    sp = P.split()
    for i in sp:
      try:
        f = float(i)
      except ValueError:
        if not self.testPartialFloat(i):
          print ("Invalid entry: ", P, ", failing at: ", i)
          print ("Valid entries: floating point numbers (white space separated)")
          return False
    return True

  def verifyString(self,P):
    sp = P.split()
    if len(sp) > 1:
      print ("Invalid entry: ", P)
      print ("Valid entries: any string without white space")
      return False
    return True

  def updateParameter(self,p):
    setattr(self.newPars,p,processPar(self.parBoxes[p][1]))

  def processPar(self,p):
    return 1

  def formatPar(self,p,useDefault=False):
    if useDefault:
      source=self.defaultPars
    else:
      source=self.newPars
    if p in self.listPars:
        s = getattr(source,p)
        return "\t".join([str(i) for i in s])
    else:
        return str(getattr(source,p))

  def updatePars(self,checkList=None):
    try:
      self.readFromGui()
    except:
      pass
    if checkList:
        for i in checkList:
          setattr(self.newPars,i,self.radioList[i][0].get())
        setattr(self.newPars,"outputType",self.defaultPars.outputOptions[self.radioList["outputType"][0].get()])
        self.fixPars(checkList)
    else:
        for i in self.options:
          setattr(self.newPars,i,self.radioList[i][0].get())
        setattr(self.newPars,"outputType",self.defaultPars.outputOptions[self.radioList["outputType"][0].get()])
        self.fixPars()
    self.drawPars()

  def fixPars(self,checkList=None):
    if checkList:
      optionList = checkList
    else:
      optionList = self.options
    for i in optionList:
      if i == "asymmetricPDs":
        if getattr(self.newPars,i) == False:
          try:
            self.needPars.remove("Lneck2List")
            self.needPars.remove("Rn2List")
          except:
            pass
          if self.newPars.computeVals == True:
            self.needPars.add("xMaxList")
            try:
              self.needPars.remove("RnList")
            except:
              pass
        elif getattr(self.newPars,i) == True:
          # use set.add -> automatically checks if option is already there.
          self.needPars.add("Lneck2List")
          self.needPars.add("Rn2List")
          if self.newPars.computeVals == True:
            try:
              self.needPars.remove("xMaxList")
            except:
              pass
            self.needPars.add("RnList")
      ## TODO: add more    

  def runClicked(self):
    fName = self.parFileName.get()
    messagebox.showinfo("Running", "Running with parameter file "+ fName)
    self.writePars(fName)
    pd.main(fName)
    messagebox.showinfo("Finished", "Finished running with parameter file "+ fName)

  def parClicked(self):
    fName = self.parFileName.get()
    messagebox.showinfo("Export", "Exporting values to parameter file "+ fName)
    self.writePars(fName)

  def infoClicked(self):
    messagebox.showinfo("Usage information", "\
Welcome to PDinsight.\n\n\
To get started, choose a mode and click 'Load default parameters'.\n\
To change mode, click 'Load default parameters' again after selecting the new modes. Note that this erases any user input.\n\n\
'Run PDinsight' writes a parameter file to the current directory and runs PDinsight. Additional output files are also written to the current directory.\n\
'Export parameter file' only writes a parameter file to the current directory.\n\n\
For certain parameters, additional information is available by clicking the description. If available, the mouse cursor turns into question mark.\n\n\
Note that the clicking buttons only works if all dialog windows are closed.")

  def readFromGui(self):
    for p in self.needPars:
      try:
          if p in self.listPars:
            setattr(self.newPars,p,self.parBoxes[p][2].get().split())
          else:
            setattr(self.newPars,p,self.parBoxes[p][2].get())
          ## missing above: check for valid entries. 
      except:
        pass
    for p in self.options:
      try:
        setattr(self.newPars,p,self.radioList[p][0].get())
      except:
        pass

  def writePars(self,fName):  
    outfile = open(fName,'w')
    self.readFromGui()
    self.newPars.validate()
    for p in self.needPars:
      self.newPars.writeSingle(p,outfile,True)
    for p in self.options:
      self.newPars.writeSingle(p,outfile,True)
    for m in self.rMode["values"]:  
      self.newPars.writeSingle(m,outfile,True)
    self.newPars.writeSingle("outputType",outfile,True)
    outfile.close()

  def __init__(self,master=None):
    tk.Frame.__init__(self,master)
    self.defaultPars = pd.Parameters()
    self.defaultPars.guiFormat()
    self.newPars = pd.Parameters()
    self.newPars.guiFormat()
    self.defaultPars.fileTag = "gui"
    self.parFileName=tk.StringVar(self,"parameters.txt")
    self.listPars= self.defaultPars.listFloats + list(self.defaultPars.listFactorFloats.keys()) + self.defaultPars.listStrings
    self.parCats = ["Particle","Channel","Distribution","Misc"]
    self.fullParLists = {"Particle": ["x","diff"], "Channel": ["Rn","RnList","Rn2List","Rdt","RdtList","Rc","RcList","Lneck","LneckList","Lneck2List","Lpd","LpdList","xMax","xMaxList"], "Distribution": ["densList","grid","gridList","dPit","dPitList","pitList","twinningList"], }
    self.pack()
    self.fill()

def main():
  gui = tk.Tk()
  gui.title("PD insight version 13")
  guiRun  = PDgui(master=gui)
  guiRun.mainloop()
  #gui.destroy()

if __name__ == '__main__':
  main()
