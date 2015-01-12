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

def onemodeldiagnostic(motifID,libDirectory,diagDirectory,prevHTML,nextHTML):
  n = 1
  if n > 0:

    FN = diagDirectory + "\\" + motifID + "_diagnostics.txt"

    # read correspondences for the given motif group; there are many such correspondences
    InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit = readcorrespondencesfromfile(FN)
    
    print "Read diagnostics from " + FN
    
    # default display color, indexed by instance; each instance will be displayed in one row
    DisplayColor = {}
    
    # loop through instances from the motif group and set the color that it will be displayed
    for i in InstanceToPDB.iterkeys():
      a = re.search("(.+Instance_[0-9]+)",i)
      DisplayColor[a.group(1)] = 'blue'            # default display color
      
    # loop through sequences from the motif group and set the default display color in a dictionary
    for i in SequenceToModel.iterkeys():
      a = re.search("(.+Sequence_[0-9]+)",i)
      DisplayColor[a.group(1)] = 'black'            # default display color

    MisAlign = 0

#    print GroupToModel
#    print SequenceToModel
#    print InstanceToPDB
          
    for nt in sorted(InstanceToPDB.iterkeys()):
      if GroupToModel[InstanceToGroup[nt]] != SequenceToModel[InstanceToSequence[nt]]:
        print nt + ' belongs to ' + GroupToModel[InstanceToGroup[nt]] + ' but was aligned to ' + SequenceToModel[InstanceToSequence[nt]]
        MisAlign += 0.5
#        a = re.search("(.+Instance_[0-9]+)",nt)
#        DisplayColor[a.group(1)] = 'red'
        a = re.search("(.+Sequence_[0-9]+)",InstanceToSequence[nt])
        DisplayColor[a.group(1)] = 'red'
                
     # Loop through instances from 3D and from the sequence alignment and put in an alignment to display
     
    aligdata = {}                                      # new dictionary
    
    for a in InstanceToGroup.iterkeys():
      m = re.search("(.+Instance_[0-9]+)",a)
      aligdata[m.group(1)] = []                        # initialize this key with empty list
    
    for a in SequenceToModel.iterkeys():
      m = re.search("(.+Sequence_[0-9]+)",a)
      aligdata[m.group(1)] = []                        # initialize this key with empty list
      
    for a in aligdata.iterkeys():
      for j in range(0,len(ModelToColumn)):
        aligdata[a].append('')                         # initialize with blank
                                                       # sorting by key should keep insertions in order
    for a in sorted(InstanceToGroup.iterkeys(), key=positionkeyforsortbynumber):
      m = re.search("(.+Instance_[0-9]+)",a)
      t = int(ModelToColumn[GroupToModel[InstanceToGroup[a]]]) # map position in group to the correct column in the model and in the alignment
      aligdata[m.group(1)][t-1] += a[len(a)-1]         # last character of the key is the base for this position
    
    for a in sorted(SequenceToModel.iterkeys(), key=positionkeyforsortbynumber):
      m = re.search("(.+Sequence_[0-9]+)",a)
      t = int(ModelToColumn[SequenceToModel[a]])
      aligdata[m.group(1)][t-1] += a[len(a)-1]
    
#    for a,b in aligdata.iteritems():
#      for i in range(0,len(b)-1):
#        print '<td>'+aligdata[a][i]+'</td>',
#      print
      
    f = open(diagDirectory+"\\"+motifID+"_GroupToModelDiagnostic.html","w")
    f.write("<html><title>"+motifID+" alignment</title>\n")    
    f.write("<h1>Alignment of "+motifID+" sequences from 3D to the JAR3D model</h1>\n")    
    f.write("<a href=\"" + prevHTML + "\">Previous group</a> | ")
    f.write("<a href=\"" + nextHTML + "\">Next group</a> | ")
    f.write("<a href=\"GroupToModelDiagnostic.html\">List of all groups</a> | ")
    f.write("<a href=\"http://rna.bgsu.edu/rna3dhub/motif/view/" + motifID + "\" target=\"_blank\">Motif atlas entry for " + motifID + "</a>  ")
    f.write("<br>The correspondence between sequences from 3D structures and the motif group is shown in blue and the JAR3D alignment of the sequences to the motif group is shown in black.  Occasionally the two disagree, in which case the JAR3D alignment is shown in red.")
    f.write("<table>")    
    f.write(alignmentheaderhtml(ModelToColumn,GroupToModel)+'\n')
    f.write(alignmentrowshtml(DisplayColor,aligdata,HasName,HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit))
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

    FASTAFile = libDirectory + "\\" + motifID + ".fasta"
    
    f.write('<b>Sequences of instances from ' + motifID + ':<b>')
    f.write('<pre>')
    with open(FASTAFile,"r") as mf:
        for line in mf.readlines():
          f.write(line)
    f.write("</pre>")

    f.write("</html>")
    f.close()
    
    print "Wrote html file with alignment of 3D instances and sequences for " + motifID
  
    return aligdata, MisAlign

def allmodelsdiagnostic(directory):

  print "Starting all models diagnostic in directory " + directory

  libDirectory = directory + "\\" + "lib"                       # replace \\ with os. separator character
  diagDirectory = directory + "\\" + "diagnostic"
  #  diagDirectory = re.sub('bp_models','interactions',libDirectory)  
 
  dirList = os.listdir(libDirectory) 
#  print dirList

  motifIDList = []
  for fn in dirList:
    if re.search("_correspondences.txt",fn):
      motifIDList += [re.sub("_correspondences.txt","",fn)]

#  print motifIDList

  prevHTML = ""

  # write the header of the overall html file

  tf = open(diagDirectory+"\\" + "GroupToModelDiagnostic.html","w")
  tf.write("<html><body><h1>Comparison of motif group and JAR3D alignment of sequences each motif group</h1>")

  # iterate through files, call onemodeldiagnostic(motifID,path)

  for i in range(0,len(motifIDList)):
#  for i in range(0,2):
    motifID = motifIDList[i]
    currentHTML = motifID + "_GroupToModelDiagnostic.html"

    print "Starting motif group " + motifID

#    m = re.search("(.*)(_diagnostics)(.*)",fn)

    # just in case you want to look at what is in the _correspondences.txt file for this motif group
    if 0 > 1:
      InteractionFile = libDirectory + "\\" + motifID + "_correspondences.txt"
      with open(InteractionFile,"r") as intf:
        for line in intf.readlines():
          print line,
        
    if i < len(motifIDList)-1:
      nextHTML = motifIDList[i+1] + "_GroupToModelDiagnostic.html"
    else:
      nextHTML = ""

    aligdata, MisAlign = onemodeldiagnostic(motifID,libDirectory,diagDirectory,prevHTML,nextHTML)
    if int(MisAlign) == 1:
      tf.write('<a href="' + currentHTML + '">' + motifID + '</a><font color = "red"> had ' + str(int(MisAlign)) + ' misalignment<font color = "black"><br>')
    elif MisAlign > 0:
      tf.write('<a href="' + currentHTML + '">' + motifID + '</a><font color = "red"> had ' + str(int(MisAlign)) + ' misalignments<font color = "black"><br>')
    else:
      tf.write('<a href="' + currentHTML + '">' + motifID + '</a> had ' + str(MisAlign) + ' misalignments<br>')

    prevHTML = currentHTML
        
  tf.close()

  return dirList
  
if __name__ == "__main__":
  allmodelsdiagnostic(sys.argv[1])  
