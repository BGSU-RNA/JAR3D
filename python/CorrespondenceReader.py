# CorrespondenceReader.py reads a text file listing correspondences between
# 3D motifs, motif groups, JAR3D models, and fasta files

import re
import string

InstanceToGroup = {}          # instance of motif to conserved group position
InstanceToPDB = {}            # instance of motif to NTs in PDBs
InstanceToSequence = {}       # instance of motif to position in fasta file
GroupToModel = {}             # positions in motif group to nodes in JAR3D model
SequenceToModel = {}          # sequence position to node in JAR3D model
ModelToColumn = {}            # nodes in JAR3D model to display columns
FN = "C:\Documents and Settings\zirbel\My Documents\My Dropbox\BGSURNA\Motifs\Correspondences\IL_018_Correspondences.txt"

with open(FN,"r") as f:
    for line in f.readlines():
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
        elif re.search("aligns_to_JAR3D",line):
            m = re.match("(.*) (.*) (.*)",line)
            SequenceToModel[m.group(1)] = m.group(3)
        elif re.search("appears_in_column",line):
            m = re.match("(.*) (.*) (.*)",line)
            ModelToColumn[m.group(1)] = m.group(3)

for nt in InstanceToPDB.iterkeys():
  if GroupToModel[InstanceToGroup[nt]] != SequenceToModel[InstanceToSequence[nt]]:
    print nt, GroupToModel[InstanceToGroup[nt]], SequenceToModel[InstanceToSequence[nt]]
      
      
 # Loop through instances from 3D and from the sequence alignment and put in an alignment to display

 
aligdata = {}                                      # new dictionary

for a in InstanceToGroup.iterkeys():
  m = re.search("(Instance_[0-9]+)",a)
  aligdata[m.group(1)] = []

for a in SequenceToModel.iterkeys():
  m = re.search("(Sequence_[0-9]+)",a)
  aligdata[m.group(1)] = []
  
for a in aligdata.iterkeys():
  for j in range(0,len(ModelToColumn)):
    aligdata[a].append('')
                                                  # sorting by key should keep insertions in order
for a in sorted(InstanceToGroup.iterkeys()):
  m = re.search("(Instance_[0-9]+)",a)
  t = int(ModelToColumn[GroupToModel[InstanceToGroup[a]]])
  aligdata[m.group(1)][t-1] += a[len(a)-1]

for a in sorted(SequenceToModel.iterkeys()):
  m = re.search("(Sequence_[0-9]+)",a)
  t = int(ModelToColumn[SequenceToModel[a]])
  aligdata[m.group(1)][t-1] += a[len(a)-1]

#for a,b in aligdata.iteritems():
#  for i in range(0,len(b)-1):
#    print '<td>'+aligdata[a][i]+'</td>',
#  print
  
f = open("testalign.html","w")
f.write("<html>"+FN+" alignment\n<table>")    
for a in sorted(aligdata.iterkeys()):
  f.write('<tr><td>'+a+'</td>')
  for i in range(0,len(aligdata[a])-1):
    f.write('<td>'+aligdata[a][i]+'</td>')
  f.write('</tr>\n')
f.write("</table></html>")
f.close()

print "Wrote html file with alignment of 3D instances and sequences"