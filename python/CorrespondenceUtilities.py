# CorrespondenceReader.py reads a text file listing correspondences between
# 3D motifs, motif groups, JAR3D models, and fasta files

import sys
import os
import re
import string

# python CorrespondenceDiagnostic.py "C:\\Documents and Settings\\zirbel\\My Documents\\My Dropbox\\BGSURNA\\Motifs\\diagnostic\\IL\\0.6\\bp_models\\"

def keyforsortbynumber(tag):
  m = re.search("([0-9]+$)",tag)
  n = m.group(1)
  while (len(n) < 20):
    n = "0" + n
  m = re.search("(Sequence_)",tag)
  if m is None:
    a = "1"
  else:
    a = "2"
  m = re.search("(^Sequence_)",tag)           # sequence from a fasta file, not from 3D structure
  if m is not None:
    n = "1" + n                                # put such sequences second
  return n + a

def positionkeyforsortbynumber(tag):
  m = re.search("Position_([0-9]+)[_-]",tag)
  n = m.group(1)
  while (len(n) < 20):
    n = "0" + n

  m = re.search("Insertion_([0-9]+)_",tag)
  if m is None:
    nn = ""
  else:
    nn = m.group(1)
  while (len(nn) < 20):
    nn = "0" + nn

  return n+nn

def columnkeyforsortbynumber(tag):
  m = re.search("Column_([0-9]+)[_-]",tag)
  n = m.group(1)
  while (len(n) < 20):
    n = "0" + n

  m = re.search("Insertion_([0-9]+)_",tag)
  if m is None:
    nn = ""
  else:
    nn = m.group(1)
  while (len(nn) < 20):
    nn = "0" + nn

  return n+nn

def readcorrespondencesfromtext(lines):

  InstanceToGroup = {}          # instance of motif to conserved group position
  InstanceToPDB = {}            # instance of motif to NTs in PDBs
  InstanceToSequence = {}       # instance of motif to position in fasta file
  GroupToModel = {}             # positions in motif group to nodes in JAR3D model
  ModelToColumn = {}            # nodes in JAR3D model to display columns
  HasName = {}                  # organism name in FASTA header
  SequenceToModel = {}          # sequence position to node in JAR3D model
  HasScore = {}                 # score of sequence against JAR3D model
  HasInteriorEdit = {}          # minimum interior edit distance to 3D instances from the motif group
  HasFullEdit = {}              # minimum full edit distance to 3D instances from the motif group
  HasCutoffValue = {}           # cutoff value 'true' or 'false'
  HasCutoffScore = {}           # cutoff score, 100 is perfect, 0 is accepted, negative means reject
  HasAlignmentScoreDeficit = {} # alignment score deficit, how far below the best score among 3D instances in this group

  for line in lines:
    if re.search("corresponds_to_group",line):
        m = re.match("(.*) (.*) (.*)",line)
        InstanceToGroup[m.group(1)] = m.group(3)
    elif re.search("corresponds_to_PDB",line):
        m = re.match("(.*) (.*) (.*)",line)
        InstanceToPDB[m.group(1)] = m.group(3)
    elif re.search("corresponds_to_JAR3D",line):
        m = re.match("(.*) (.*) (.*)",line)
        GroupToModel[m.group(1)] = m.group(3)
    elif re.search("corresponds_to_sequence",line):
        m = re.match("(.*) (.*) (.*)",line)
        InstanceToSequence[m.group(1)] = m.group(3)
    elif re.search("appears_in_column",line):
        m = re.match("(.*) (.*) (.*)",line)
        ModelToColumn[m.group(1)] = m.group(3)
    elif re.search("aligns_to_JAR3D",line):
        m = re.match("(.*) (.*) (.*)",line)
        SequenceToModel[m.group(1)] = m.group(3)
    elif re.search("has_name",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasName[m.group(1)] = m.group(3)
    elif re.search("has_score",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasScore[m.group(1)] = m.group(3)
    elif re.search("has_minimum_interior_edit_distance",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasInteriorEdit[m.group(1)] = m.group(3)
    elif re.search("has_minimum_full_edit_distance",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasFullEdit[m.group(1)] = m.group(3)
    elif re.search("has_cutoff_value",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasCutoffValue[m.group(1)] = m.group(3)
    elif re.search("has_cutoff_score",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasCutoffScore[m.group(1)] = m.group(3)
    elif re.search("has_alignment_score_deficit",line):
        m = re.match("(.*) (.*) (.*)",line)
        HasAlignmentScoreDeficit[m.group(1)] = m.group(3)

  return InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit

def readcorrespondencesfromfile(filenamewithpath):

    InstanceToGroup = {}          # instance of motif to conserved group position
    InstanceToPDB = {}            # instance of motif to NTs in PDBs
    InstanceToSequence = {}       # instance of motif to position in fasta file
    GroupToModel = {}             # positions in motif group to nodes in JAR3D model
    ModelToColumn = {}            # nodes in JAR3D model to display columns
    HasName = {}                  # header lines from FASTA file
    SequenceToModel = {}          # sequence position to node in JAR3D model
    HasScore = {}                 # score of sequence against JAR3D model
    HasInteriorEdit = {}          # minimum interior edit distance to 3D instances from the motif group
    HasFullEdit = {}              # minimum full edit distance to 3D instances from the motif group
    HasCutoffValue = {}           # cutoff value 'true' or 'false'
    HasCutoffScore = {}           # cutoff score, 100 is perfect, 0 is accepted, negative means reject
    HasAlignmentScoreDeficit = {} # alignment score deficit, how far below the best score among 3D instances in this group

    with open(filenamewithpath,"r") as f:
      lines = f.readlines()

    InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit = readcorrespondencesfromtext(lines)

    return InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit

def alignmentheaderhtml(ModelToColumn,GroupToModel):

  ColumnHeader = [''] * len(ModelToColumn)
  for a in ModelToColumn.keys():
    ColumnHeader[int(ModelToColumn[a])-1] = a

  PositionNumber = [''] * (len(ModelToColumn)+1)
  for a in GroupToModel.keys():
    colnum = ModelToColumn[GroupToModel[a]]
    m = re.search("Column_([0-9]+)$",a)
    if m is not None:
      PositionNumber[int(colnum)] = m.group(1)

  t = '<tr valign=top><th align="left">Column number:<br>Node in JAR3D model:<br>Insertion positions indicated by I:<br>Position in motif group:</th>'
  for i in range(0,len(ModelToColumn)):
    t = t + '<th>'
    t = t + str(i+1) + '<br>'
    m = re.search("Node_([0-9]+)",ColumnHeader[i])
    a = m.group(1)
    t = t + a + '<br>'
    if re.search("Insertion",ColumnHeader[i]):
      t = t + 'I'
    t = t + '<br>' + PositionNumber[i+1] + '<br>'
    t = t + '</th>'
  t = t + '<th><br>In<br>acceptance<br>region</th>'
  t = t + '<th><br><br>Cutoff<br>score</th>'
  t = t + '<th><br>Full<br>edit<br>distance</th>'
  t = t + '<th><br>Interior<br>edit<br>distance</th>'
  t = t + '<th><br>Alignment<br>score<br>deficit</th>'
  t = t + '</tr>'

  return t

def alignmentrowshtml(DisplayColor,aligdata,HasName,HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit):

  t = ""

  for a in sorted(aligdata.keys(),key = keyforsortbynumber):
    if HasScore[a] == ".":
      Score = 0
      DisplayScore = ""
    else:
      Score = float(HasScore[a])
      DisplayScore = "%0.4f" % Score

    if Score > -1000:                                   # an arbitrary minimum score
      t = t + '<tr><td><font color = "' + DisplayColor[a]+ '">'+a+' '+HasName[a]+'</td>'
      for i in range(0,len(aligdata[a])):
        t = t + '<td><font color = "' + DisplayColor[a] + '">'+aligdata[a][i]+'</td>'
      t = t + '<td align="right">' + HasCutoffValue.get(a,'') + '</td>'
      x = HasCutoffScore.get(a,'')
      CutoffScore = "%.0f" % float(x) if '.' in x else ''
      t = t + '<td align="right">' + CutoffScore + '</td>'
      t = t + '<td align="right">' + HasFullEdit.get(a,'') + '</td>'
      t = t + '<td align="right">' + HasInteriorEdit.get(a,'') + '</td>'     # '' is the default
      x = HasAlignmentScoreDeficit.get(a,'')
      Deficit = "%.2f" % float(x) if '.' in x else ''
      t = t + '<td align="right">' + Deficit + '</td>'     # '' is the default
    else:
      t = t + '<tr><td><font color = "' + DisplayColor[a]+ '">'+a+' '+HasName[a]+' (no alignment possible)</td>\n'
      for i in range(0,len(aligdata[a])):
        t = t + '<td></td>'
      t = t + '<td align="right">' + HasCutoffValue.get(a,'') + '</td>'
      t = t + '<td></td>'
      t = t + '<td align="right">' + HasFullEdit.get(a,'') + '</td>'
      t = t + '<td align="right">' + HasInteriorEdit.get(a,'') + '</td>'     # '' is the default
      t = t + '<td></td>'     # '' is the default
    t = t + '</tr>\n'

  return t

def alignsequencesandinstancesfromtext(MotifCorrespondenceText,SequenceCorrespondenceText):

  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn = readcorrespondencesfromtext(MotifCorrespondenceText)[:5]
  NotInstanceToGroup, NotInstanceToPDB, NotInstanceToSequence, NotGroupToModel, NotModelToColumn, SequenceToModel = readcorrespondencesfromtext(SequenceCorrespondenceText)[:6]

  motifalig = {}

  for a in InstanceToGroup.keys():
    m = re.search("(Instance_[0-9]+)",a)
    motifalig[m.group(1)] = [''] * len(ModelToColumn)     # start empty

  for a in sorted(InstanceToGroup.keys()):
    m = re.search("(Instance_[0-9]+)",a)
    t = int(ModelToColumn[GroupToModel[InstanceToGroup[a]]])
    motifalig[m.group(1)][t-1] += a[len(a)-1]

  sequencealig = {}

  for a in SequenceToModel.keys():
    m = re.search("(Sequence_[0-9]+)",a)
    sequencealig[m.group(1)] = [''] * len(ModelToColumn)  # start empty

  for a in sorted(SequenceToModel.keys()):
    m = re.search("(Sequence_[0-9]+)",a)
    t = int(ModelToColumn[SequenceToModel[a]])
    sequencealig[m.group(1)][t-1] += a[len(a)-1]

  header = {}              # new dictionary
  header['columnname'] = [''] * len(ModelToColumn)
  header['nodes'] = [''] * len(ModelToColumn)
  header['positions'] = [''] * len(ModelToColumn)
  header['insertions'] = [''] * len(ModelToColumn)

  for a in ModelToColumn.keys():
    header['columnname'][int(ModelToColumn[a])-1] = a

  for i in range(0,len(ModelToColumn)):
    m = re.search("Node_([0-9]+)",header['columnname'][i])
    a = m.group(1)
    header['nodes'][i] = a
    if re.search("Insertion",header['columnname'][i]):
      header['insertions'][i] = 'Insertion'

  for a in GroupToModel.keys():
    m = re.search("Column_([0-9]+)$",a)
    if m is not None:
      colnum = ModelToColumn[GroupToModel[a]]
      header['positions'][int(colnum)] = m.group(1)

  return header, motifalig, sequencealig

def alignsequencesandinstancesfromfiles(MotifCorrespondenceFile,SequenceCorrespondenceFile):

  with open(MotifCorrespondenceFile,"r") as f:
    MotifLines = f.readlines()

  print("Read motif correspondence file " + MotifCorrespondenceFile)

  with open(SequenceCorrespondenceFile,"r") as f:
    SequenceLines = f.readlines()

  print("Read sequence correspondence file " + SequenceCorrespondenceFile)

  header, motifalig, sequencealig = alignsequencesandinstancesfromtext(MotifLines,SequenceLines)

  return header, motifalig, sequencealig

# Use JAR3D to create a file of sequence-to-model correspondences
# Next, specify a motif correspondence file and a file of sequence to model correspondences

def alignsequencesandinstancesfromfileshtml(MotifCorrespondenceFile,SequenceCorrespondenceFile):

    InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, NotSequenceToModel = readcorrespondencesfromfile(MotifCorrespondenceFile)

    NotInstanceToGroup, NotInstanceToPDB, NotInstanceToSequence, NotGroupToModel, NotModelToColumn, SequenceToModel = readcorrespondencesfromfile(SequenceCorrespondenceFile)

    DisplayColor = {}

    for i in InstanceToPDB.keys():
      a = re.search("(Instance_[0-9]+)",i)
      DisplayColor[a.group(1)] = 'black'            # default display color

    for i in SequenceToModel.keys():
      a = re.search("(Sequence_[0-9]+)",i)
      DisplayColor[a.group(1)] = 'black'            # default display color

    # Loop through instances from 3D and from the sequence alignment and put in an alignment to display

    aligdata = {}                                      # new dictionary

    for a in InstanceToGroup.keys():
      m = re.search("(Instance_[0-9]+)",a)
      aligdata[m.group(1)] = []

    for a in SequenceToModel.keys():
      m = re.search("(Sequence_[0-9]+)",a)
      aligdata[m.group(1)] = []

    for a in aligdata.keys():
      for j in range(0,len(ModelToColumn)):
        aligdata[a].append('')                        # clumsy but effective
                                                      # sorting by key should keep insertions in order
    for a in sorted(InstanceToGroup.keys()):
      m = re.search("(Instance_[0-9]+)",a)
      t = int(ModelToColumn[GroupToModel[InstanceToGroup[a]]])
      aligdata[m.group(1)][t-1] += a[len(a)-1]

    for a in sorted(SequenceToModel.keys()):
      m = re.search("(Sequence_[0-9]+)",a)
      t = int(ModelToColumn[SequenceToModel[a]])
      aligdata[m.group(1)][t-1] += a[len(a)-1]

    t = '<table>'
    t = t + alignmentheaderhtml(ModelToColumn)+'\n'
    t = t + alignmentrowshtml(DisplayColor,aligdata)
    t = t + '<table>'

    return t


# python CorrespondenceUtilities.py "C:\\Documents and Settings\\zirbel\\My Documents\\My Dropbox\\BGSURNA\\Motifs\\lib\\IL\\0.6\\bp_models\\IL_01239.1_correspondences.txt" "C:\\Documents and Settings\\zirbel\\My Documents\\My Dropbox\\BGSURNA\\Motifs\\diagnostic\\IL\\0.6\\bp_models\\IL_01239.1_diagnostics.txt"

if __name__ == "__main__":
#  alignsequencesandinstancesfromfiles(sys.argv[1],sys.argv[2])

  header, motifalig, sequencealig = alignsequencesandinstancesfromfiles(sys.argv[1],sys.argv[2])
  print("Result from alignsequencesandinstancesfromfiles: ")
  print(header)
  print(motifalig)
  print(sequencealig)


