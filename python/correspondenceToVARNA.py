# fastaToModelToVARNA reads the alignment of a single sequence in a FASTA file,
# reads the alignment to a motif group, and writes out a VARNA file to
# describe the basepair diagram to make

import sys
import os
import re
import string

from CorrespondenceUtilities import readcorrespondencesfromfile
from CorrespondenceUtilities import alignmentrowshtml
from CorrespondenceUtilities import alignmentheaderhtml
from CorrespondenceUtilities import keyforsortbynumber
from CorrespondenceUtilities import positionkeyforsortbynumber
from CorrespondenceUtilities import columnkeyforsortbynumber

def correspondenceToVARNA(libDirectory,motifID,correspondenceFile,outputFileName):
  # read correspondences from the fasta file to the model
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit = readcorrespondencesfromfile(correspondenceFile)

  print "Read alignment to model from " + correspondenceFile

#  print(SequenceToModel)
  sequenceNumber = 1             # some FASTA files may have multiple sequences, need to say which

  FASTAdict = {}
  for sp in SequenceToModel:
    fields = sp.split("_")
    FASTAdict[fields[3]] = fields[4]
  FASTASequence = ""
  for sp in range(0,len(SequenceToModel)):
    FASTASequence += FASTAdict[str(sp+1)]

  ModelToSequence = {}
  for sp in SequenceToModel:
    ModelToSequence[SequenceToModel[sp]] = sp     # will use latest correspondence, for multiple insertions

  FN = libDirectory + "\\" + motifID + "_correspondences.txt"

  interactionsFile = libDirectory + "\\" + motifID + "_interactions.txt"

  with open(interactionsFile,'r') as fh:
    interactions = fh.readlines()

  # read correspondences for the given motif group; there are many such correspondences
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModelDummy, ModelHasName, ModelHasScore, ModelInteriorEdit, ModelFullEdit, ModelCutoffValue, ModelCutoffScore, ModelDeficit = readcorrespondencesfromfile(FN)

#  print(GroupToModel)

  allBPs = ""

  for interaction in interactions:
    fields = interaction.split(" ")
    a = fields[0]
    b = fields[1]
    interactionType = fields[2]
    if "c" in interactionType or "t" in interactionType:
      try:
        aa = ModelToSequence[GroupToModel[motifID+"_Column_"+a]]
        bb = ModelToSequence[GroupToModel[motifID+"_Column_"+b]]

  #      viroid_loop_4... has sequence position to Model as well, just need to work it backward in one case
  #      is there a many to one problem, though?

        print(a,b,interactionType,aa,bb)

        aafields = aa.split('_')
        bbfields = bb.split('_')
        allBPs += "("+aafields[3]+","+bbfields[3]+")"
        orientation = interactionType[0].replace("t","trans").replace("c","cis")
        firstedge = string.lower(interactionType[1]).replace("w","wc")
        secondedge = string.lower(interactionType[2]).replace("w","wc")
        allBPs += ':edge5="'+firstedge+'",edge3="'+secondedge+'",stericity="'+orientation+'";'
      except:
        print("Missing interaction " + interaction)

  option = {}
  option["sequenceDBN"] = FASTASequence
  option["structureDBN"] = "."*len(FASTASequence)
  option["auxBPs"] = allBPs
  if "Loop_1_" in outputFileName:
    option["rotation"] = "270.0"
  else:
    option["rotation"] = "90.0"
  option["eror"] = "True"
  #option["drawBackbone"] = "false"
  #option["annotations"] = "27:anchor=1,type=B"

  VARNA = "java -cp C:\\Users\\zirbel\\Documents\\JAR3D\\viroid\\VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd "
  for a in option:
    VARNA += " -" + a + " " + option[a]
  VARNA += " -o " + outputFileName

  print(VARNA)

  return VARNA

if __name__ == "__main__":
  fastatomodelalignment(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
