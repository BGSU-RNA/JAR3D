# GroupToModelDiagnostic.py reads a text file listing correspondences between
# 3D motifs, motif groups, JAR3D models, and fasta files

import sys
import os
import re
import string

from CorrespondenceUtilities import readcorrespondencesfromfile
from CorrespondenceUtilities import alignmentrowshtml
from CorrespondenceUtilities import alignmentheaderhtml
from CorrespondenceUtilities import keyforsortbynumber
from CorrespondenceUtilities import positionkeyforsortbynumber

def fastatomodelalignment(libDirectory,motifID,alignmentfile,outputfile):
  # read correspondences from the fasta file to the model
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit = readcorrespondencesfromfile(alignmentfile)

  print "Read alignment to model from " + alignmentfile

  FN = libDirectory + "\\" + motifID + "_correspondences.txt"

  # read correspondences for the given motif group; there are many such correspondences
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModelDummy, ModelHasName, ModelHasScore, ModelInteriorEdit, ModelFullEdit, ModelCutoffValue, ModelCutoffScore, ModelDeficit = readcorrespondencesfromfile(FN)

#  print HasScore
#  print ModelHasName

  HasName.update(ModelHasName)
  HasScore.update(ModelHasScore)

  print "Read model correspondences from " + FN
              
  # Loop through instances from 3D and from the sequence alignment and put in an alignment to display
  DisplayColor = {}
  aligdata = {}                                      # new dictionary
  
  for a in InstanceToGroup.iterkeys():
    m = re.search("(.+Instance_[0-9]+)",a)
    aligdata[m.group(1)] = []                        # initialize this key with empty list
    DisplayColor[m.group(1)] = 'blue'                # default display color

  for a in SequenceToModel.iterkeys():
    m = re.search("(Sequence_[0-9]+)",a)
    aligdata[m.group(1)] = []                        # initialize this key with empty list
    DisplayColor[m.group(1)] = 'black'               # default display color
    
  for a in aligdata.iterkeys():
    for j in range(0,len(ModelToColumn)):
      aligdata[a].append('')                         # initialize with blank
                                                     # sorting by key should keep insertions in order
  for a in sorted(InstanceToGroup.iterkeys(), key=positionkeyforsortbynumber):
    m = re.search("(.+Instance_[0-9]+)",a)
    t = int(ModelToColumn[GroupToModel[InstanceToGroup[a]]]) # map position in group to the correct column in the model and in the alignment
    aligdata[m.group(1)][t-1] += a[len(a)-1]         # last character of the key is the base for this position
  
  for a in sorted(SequenceToModel.iterkeys(), key=positionkeyforsortbynumber):
    m = re.search("(Sequence_[0-9]+)",a)
    t = int(ModelToColumn[SequenceToModel[a]])
    aligdata[m.group(1)][t-1] += a[len(a)-1]
      
  f = open(outputfile,"w")
  f.write("<html><title>Alignment to "+motifID+"</title>\n")    
  f.write("<h1>Alignment of " + alignmentfile +" to "+motifID+"</h1>\n")    
  f.write("<a href=\"http://rna.bgsu.edu/rna3dhub/motif/view/" + motifID + "\" target=\"_blank\">Motif atlas entry for " + motifID + "</a><br>")
  f.write("The correspondence between sequences from 3D structures and the motif group is shown in blue, JAR3D alignments of sequences to the motif group are shown in black, and sequences which are too long or too short to be aligned are indicated by : characters.")
  f.write("<table>")
  f.write(alignmentheaderhtml(ModelToColumn, GroupToModel)+'\n')
  f.write(alignmentrowshtml(DisplayColor, aligdata, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit))
  f.write("</table>")
 
  InteractionsFile = libDirectory + "\\" + motifID + "_interactions.txt"
  
  f.write('<br><b>Conserved interactions between motif group positions in ' + motifID + ':</b>')
  f.write('<pre>')
  with open(InteractionsFile,"r") as mf:
      for line in mf.readlines():
        f.write(line)
  f.write("</pre>")

  ModelFile = libDirectory + "\\" + motifID + "_model.txt"
  f.write('<b>JAR3D SCFG/MRF model for ' + motifID + ':</b>')
  f.write('<pre>')
  with open(ModelFile,"r") as mf:
      for line in mf.readlines():
        f.write(line)
  f.write("</pre>")
  f.write("</html>")
  f.close()
  
  print "Wrote html file with alignment of 3D instances and sequences for " + motifID

  return aligdata
  
if __name__ == "__main__":
  fastatomodelalignment(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])  
