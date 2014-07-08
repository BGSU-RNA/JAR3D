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

def fastatomodelalignment(motifID,libDirectory,fastafile,outputfile):
  # read correspondences from the fasta file to the model
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore = readcorrespondencesfromfile(fastafile)

  print "Read alignment to model from " + fastafile

  FN = libDirectory + "\\" + motifID + "_correspondences.txt"

  # read correspondences for the given motif group; there are many such correspondences
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModelDummy, ModelHasName, ModelHasScore = readcorrespondencesfromfile(FN)

  HasName.update(ModelHasName)
  HasScore.update(ModelHasScore)

  print HasScore
  
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
  f.write("<h1>Alignment of " + fastafile +" to "+motifID+"</h1>\n")    
  f.write("The sequences of instances from 3D structures are listed first and in blue<br>")
  f.write("<a href=\"http://rna.bgsu.edu/rna3dhub/motif/view/" + motifID + "\">Motif atlas entry for " + motifID + "</a>  ")
  f.write("<table>")    
  f.write(alignmentheaderhtml(ModelToColumn)+'\n')
  f.write(alignmentrowshtml(DisplayColor,aligdata,HasName,HasScore))
  f.write("</table>")
 
  f.write('<br>SCFG/MRF modelf for ' + motifID)
  f.write('<pre>')
  ModelFile = libDirectory + "\\" + motifID + "_model.txt"
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
