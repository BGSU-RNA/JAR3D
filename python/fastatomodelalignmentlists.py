# fastatomodelalignmentlists.py reads text files listing correspondences between 3D instances, their motif group, JAR3D models, and fasta files

import sys
import os
import re
import string

from CorrespondenceUtilities import readcorrespondencesfromfile
from CorrespondenceUtilities import alignmentrowshtml
from CorrespondenceUtilities import alignmentheaderhtml
from CorrespondenceUtilities import keyforsortbynumber
from CorrespondenceUtilities import positionkeyforsortbynumber

def fastatomodelalignmentlists(libDirectory,motifID,alignmentfile,outputfile):
  # read correspondences from the fasta file to the model
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit = readcorrespondencesfromfile(alignmentfile)

  print "Read alignment to model from " + alignmentfile

  FN = libDirectory + "\\" + motifID + "_correspondences.txt"

  # read correspondences for the given motif group; there are many such correspondences
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModelDummy, ModelHasName, ModelHasScore, ModelInteriorEdit, ModelFullEdit, ModelCutoffValue, ModelCutoffScore, ModelDeficit = readcorrespondencesfromfile(FN)

#  print HasScore
#  print ModelHasName

#  HasName.update(ModelHasName)
#  HasScore.update(ModelHasScore)

  print "Read model correspondences from " + FN

  ColumnHeader = [''] * len(ModelToColumn)
  for a in ModelToColumn.iterkeys():
    ColumnHeader[int(ModelToColumn[a])-1] = a

  PositionNumber = [''] * (len(ModelToColumn)+1)
  for a in GroupToModel.iterkeys():
    colnum = ModelToColumn[GroupToModel[a]]
    m = re.search("Position_([0-9]+)$",a)
    if m is not None:
      PositionNumber[int(colnum)] = m.group(1)

  HeaderColumnNumber = ['Column number'];
  for i in range(0,len(ModelToColumn)):
    HeaderColumnNumber.append(str(i+1))

  HeaderNodeNumber = ['Node number'];
  for i in range(0,len(ModelToColumn)):
    m = re.search("Node_([0-9]+)",ColumnHeader[i])
    a = m.group(1)
    HeaderNodeNumber.append(a)

  HeaderInsertion = ['Insertion positions indicated by I'];
  for i in range(0,len(ModelToColumn)):
    if re.search("Insertion",ColumnHeader[i]):
      HeaderInsertion.append('I')
    else:
      HeaderInsertion.append('')

  HeaderPosition = ['Position in motif group'];
  for i in range(0,len(ModelToColumn)):
    HeaderPosition.append(PositionNumber[i+1])

  HeaderPosition.append('In acceptance region')
  HeaderPosition.append('Cutoff score')
  HeaderPosition.append('Full edit distance')
  HeaderPosition.append('Interior edit distance')
  HeaderPosition.append('Alignment score deficit')

  print HeaderColumnNumber
  print HeaderNodeNumber
  print HeaderInsertion
  print HeaderPosition

  # Loop through instances from 3D and from the sequence alignment and put in an alignment to display
  InstanceList = {}                                  # dictionary of lists for alignment of 3D instances to the motif group

  for a in InstanceToGroup.iterkeys():
    m = re.search("(.+Instance_[0-9]+)",a)
    InstanceList[m.group(1)] = []
    for j in range(0,1+len(ModelToColumn)):
      InstanceList[m.group(1)].append('')                         # initialize with blank

  for a in sorted(InstanceToGroup.iterkeys(), key=positionkeyforsortbynumber):
    m = re.search("(.+Instance_[0-9]+)",a)
    InstanceList[m.group(1)][0] = ModelHasName[m.group(1)]
    t = int(ModelToColumn[GroupToModel[InstanceToGroup[a]]]) # map position in group to the correct column in the model and in the alignment
    InstanceList[m.group(1)][t] += a[len(a)-1]         # last character of the key is the base for this position

  print InstanceList

  SequenceList = {}                                  # disctionary of lists for alignment of new sequences to the motif group

  for a in SequenceToModel.iterkeys():
    m = re.search("(Sequence_[0-9]+)",a)
    SequenceList[m.group(1)] = []
    for j in range(0,1+len(ModelToColumn)):
      SequenceList[m.group(1)].append('')                         # initialize with blank

  for a in sorted(SequenceToModel.iterkeys(), key=positionkeyforsortbynumber):
    m = re.search("(Sequence_[0-9]+)",a)
    SequenceList[m.group(1)][0] = HasName[m.group(1)]
    t = int(ModelToColumn[SequenceToModel[a]]) # map position in group to the correct column in the model and in the alignment
    SequenceList[m.group(1)][t] += a[len(a)-1]               # last character of the key is the base for this position

  for a in SequenceList.iterkeys():
    SequenceList[a].append(HasCutoffValue.get(a,''))
    x = HasCutoffScore.get(a,'')
    CutoffScore = "%.0f" % float(x) if '.' in x else ''
    SequenceList[a].append(CutoffScore)
    SequenceList[a].append(HasFullEdit.get(a,''))
    SequenceList[a].append(HasInteriorEdit.get(a,''))
    x = HasAlignmentScoreDeficit.get(a,'')
    Deficit = "%.2f" % float(x) if '.' in x else ''
    SequenceList[a].append(Deficit)


  print SequenceList

  
if __name__ == "__main__":
  fastatomodelalignmentlists(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])  
