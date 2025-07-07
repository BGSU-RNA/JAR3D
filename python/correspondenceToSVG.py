# fastaToModelToSVG reads the alignment of a single sequence in a FASTA file,
# reads the alignment to a motif group, and writes out an SVG file to
# describe the basepair diagram to make

import sys
import os
import re
import string
import numpy as np
import math

from CorrespondenceUtilities import readcorrespondencesfromfile
from CorrespondenceUtilities import alignmentrowshtml
from CorrespondenceUtilities import alignmentheaderhtml
from CorrespondenceUtilities import keyforsortbynumber
from CorrespondenceUtilities import positionkeyforsortbynumber
from CorrespondenceUtilities import columnkeyforsortbynumber

def correspondenceToInteractionList(libDirectory,motifID,correspondenceFile):
  # read correspondences from the fasta file to the model
  InstanceToGroup, InstanceToPDB, InstanceToSequence, GroupToModel, ModelToColumn, SequenceToModel, HasName, HasScore, HasInteriorEdit, HasFullEdit, HasCutoffValue, HasCutoffScore, HasAlignmentScoreDeficit = readcorrespondencesfromfile(correspondenceFile)

#  print "Read alignment to model from " + correspondenceFile

  sequenceNumber = 1             # some FASTA files may have multiple sequences, need to say which

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

  allInteractions = []

  # record all interactions in terms of the sequence positions
  for interaction in interactions:
    fields = interaction.split(" ")
    a = fields[0]  # motif group positions of the two interacting bases
    b = fields[1]
    interactionType = fields[2]

    try:
      aa = ModelToSequence[GroupToModel[motifID+"_Column_"+a]]
      bb = ModelToSequence[GroupToModel[motifID+"_Column_"+b]]

      aafields = aa.split('_')
      bbfields = bb.split('_')

      newInteraction = {}
      newInteraction["firstSequencePosition"] = int(aafields[3])-1   # JAR3D counts from 1, Python counts from 0
      newInteraction["secondSequencePosition"] = int(bbfields[3])-1
      newInteraction["type"] = interactionType

      allInteractions.append(newInteraction)
    except:
      print("Deleted interaction " + interaction)


  return allInteractions

def drawLWSymbol(orientation,edge,c,color,rotation=[]):

  Hoogsteen = 2*np.array(([1,-1,-1,1],[1,1,-1,-1]))   # corners of square for Hoogsteen symbol
  Sugar     = 3*np.array(([1,-0.5,-0.5],[0,math.sqrt(3)/2,-math.sqrt(3)/2]))  # corners of triangle
  tWradius = 2
  cWradius = 2.5

  SVGtext = ""
  if orientation == "c":
    if edge == "W":
      SVGtext += '<circle cx="%0.8f" cy="%0.8f" r="%0.8f" stroke="none" stroke-width="1.0" fill="%s"/>\n' % (c[0],c[1],cWradius,color)
    elif edge == "H":
      square = np.column_stack((c,c,c,c)) + np.matmul(rotation,Hoogsteen)
      SVGtext += '<path d="M%0.8f %0.8f L%0.8f %0.8f L%0.8f %0.8f L%0.8f %0.8f Z" ' % tuple(np.reshape(square,8,order='F'))
      SVGtext += 'style="fill:%s; stroke:none; stroke-width:1.0;"/>\n' % color
    elif edge == "S":
      triangle = np.column_stack((c,c,c)) + np.matmul(rotation,Sugar)
      SVGtext += '<path d="M%0.8f %0.8f L%0.8f %0.8f L%0.8f %0.8f Z" ' % tuple(np.reshape(triangle,6,order='F'))
      SVGtext += 'style="fill:%s; stroke:none; stroke-width:1.0;"/>\n' % color
  elif orientation == "t":
    if edge == "W":
      SVGtext += '<circle cx="%0.8f" cy="%0.8f" r="%0.8f" stroke="%s" stroke-width="1.0" fill="rgb(100%%, 100%%, 100%%)"/>\n' % (c[0],c[1],tWradius,color)
    elif edge == "H":
      square = np.column_stack((c,c,c,c)) + np.matmul(rotation,Hoogsteen)
      SVGtext += '<path d="M%0.8f %0.8f L%0.8f %0.8f L%0.8f %0.8f L%0.8f %0.8f Z" ' % tuple(np.reshape(square,8,order='F'))
      SVGtext += 'style="fill:rgb(100%%, 100%%, 100%%); stroke:%s; stroke-width:1.0;"/>\n' % color
    elif edge == "S":
      triangle = np.column_stack((c,c,c)) + np.matmul(rotation,Sugar)
      SVGtext += '<path d="M%0.8f %0.8f L%0.8f %0.8f L%0.8f %0.8f Z" ' % tuple(np.reshape(triangle,6,order='F'))
      SVGtext += 'style="fill:rgb(100%%, 100%%, 100%%); stroke:%s; stroke-width:1.0;"/>\n' % color

  return SVGtext

def somethingBetween(SVG,i,j):
  for a in range(0,len(SVG)):
    if not a == i and not a == j and abs(np.linalg.norm(SVG[i]["p"]-SVG[j]["p"]) - np.linalg.norm(SVG[i]["p"]-SVG[a]["p"]) - np.linalg.norm(SVG[a]["p"]-SVG[j]["p"])) < 0.01:
      return True

  return False

def interactionListToSVG(sequence,interactionList,numbers=[]):

  leftX = 10.0
  deltaX = 25.0
  topStrandY = 20.0
  bottomStrandY = 50.0
  arcRadiusY = 8
  arcShift = np.array([0,arcRadiusY])
  radius = 5.0               # for circles around letters
  textShift = 2.75           # shift capital letter text down to get into centers of circles
  numberShiftTop = 6.6       # how far to put numbers above/below letters
  numberShiftBottom = 12    # how far to put numbers above/below letters
  GCseparation = 2.5         # how far apart the GC cWW basepair lines should be
  bulge = (bottomStrandY - topStrandY) / 4.0        # how far out to move a bulged base
  instrand = 0 * (bottomStrandY - topStrandY) / 4.5        # how far in to move an in-strand paired base
  scarlet = "rgb(187,0,0)"
  lightScarlet = "rgb(255,62,150)"
  pink = "rgb(255,192,203)"

  SVG = {}
  SVGtext = ""

  currentX = leftX

  stacking = np.zeros((len(sequence),len(sequence)))
  pairs = np.zeros((len(sequence),len(sequence)))

  # track all stacking and pairing interactions
  for interaction in interactionList:
    t = interaction["type"]
    i = interaction["firstSequencePosition"]
    j = interaction["secondSequencePosition"]
    if "s" in t:
      stacking[i][j] = 1
      stacking[j][i] = 1
    elif "c" in t or "t" in t:
      pairs[i][j] = 1
      pairs[j][i] = 1

  both = stacking + pairs

  # assign x and y coordinates to the bases
  if "*" in sequence:
    s = sequence.find("*")    # location of the star
    topNum = s
    bottomNum = len(sequence)-s-1
    width = deltaX * max(topNum-1,bottomNum-1)
    deltaXTop = width / (topNum-1)
    deltaXBottom = width / (bottomNum-1)
    y = topStrandY
    upper = range(0,s)
    lower = range((s+1),len(sequence))

    movedIn = False
    for i in range(0,len(sequence)):
      SVG[i] = {}
      if i < s:
        if np.sum(stacking[i,...]) == 0 and np.sum(pairs[i,...]) == 0:  # base all by itself, move out
          SVG[i]["p"] = np.array([currentX,topStrandY-bulge])
        elif not movedIn and i > 0 and i < s-1 and np.sum(pairs[i,upper]) > 0 and np.sum(pairs[i,lower]) > 0 and np.sum(both[0:i,i:(s-1)]) > 1:
          SVG[i]["p"] = np.array([currentX,topStrandY+instrand])
          movedIn = True
        elif not movedIn and i < s-1 and pairs[0,i] > 0 and (i>1 or np.sum(both[0,upper]) > 1):
          SVG[i]["p"] = np.array([currentX,topStrandY+instrand])
          movedIn = True
        elif not movedIn and i > 0 and pairs[i,s-1] > 0 and (s-1-i>1 or np.sum(both[s-1,upper]) > 1):
          SVG[i]["p"] = np.array([currentX,topStrandY+instrand])
          movedIn = True
        else:
          SVG[i]["p"] = np.array([currentX,topStrandY])
          movedIn = False
        currentX += deltaXTop
      elif i == s:
        currentX -= deltaXTop
        movedIn = False
        SVG[i]["p"] = np.array([float('inf'),float('inf')])
      elif i > s:
        if np.sum(stacking[i,...]) == 0 and np.sum(pairs[i,...]) == 0:  # base all by itself
          SVG[i]["p"] = np.array([currentX,bottomStrandY+bulge])
        elif i == s+1 and pairs[i,len(sequence)-1] > 0:
          SVG[i]["p"] = np.array([currentX,bottomStrandY-instrand])
          movedIn = True
        elif not movedIn and i > s+1 and i+1 < len(sequence) and np.sum(pairs[i,lower]) > 0 and np.sum(pairs[i,upper]) and np.sum(both[(s+1):i,i:(len(sequence)-1)]) > 1:     # another interaction crosses over i
          SVG[i]["p"] = np.array([currentX,bottomStrandY-instrand])
          movedIn = True
        elif not movedIn and i > s+1 and i+1 < len(sequence) and pairs[i,len(sequence)-1] > 0 and (len(sequence)-1-i > 1 or np.sum(both[len(sequence)-1,lower]) > 1):
          SVG[i]["p"] = np.array([currentX,bottomStrandY-instrand])
          movedIn = True
        elif not movedIn and i > s+1 and pairs[i,s+1] > 0 and (i-(s+1)>1 or np.sum(both[s+1,lower]) > 1):
          SVG[i]["p"] = np.array([currentX,bottomStrandY-instrand])
          movedIn = True
        else:
          SVG[i]["p"] = np.array([currentX,bottomStrandY])
          movedIn = False
        currentX -= deltaXBottom
  else:                                  # hairpin x,y layout
    i = 0
    while i < len(sequence)-1 and stacking[i,i+1] > 0:
      SVG[i] = {}
      SVG[i]["p"] = np.array([currentX,topStrandY])
      currentX += deltaX
      i += 1

    i = len(sequence)-1
    while i < len(sequence)-1 and stacking[i,i+1] > 0:
      SVG[i] = {}
      SVG[i]["p"] = np.array([currentX,bottomStrandY])
      currentX += deltaX
      i -= 1

    print(SVG)

  # draw the stacking interactions
  for interaction in interactionList:
    t = interaction["type"]
    if t[0] == "s":
      i = interaction["firstSequencePosition"]
      j = interaction["secondSequencePosition"]
      p1 = SVG[i]["p"]
      p2 = SVG[j]["p"]
      x1 = p1[0]
      y1 = p1[1]
      x2 = p2[0]
      y2 = p2[1]

      if abs(y1-y2) < 0.01 and abs(i-j) > 1 and somethingBetween(SVG,i,j):
        if y1 > (topStrandY + bottomStrandY)/2:
          SVGtext += '<path d="M%0.8f,%0.8f A%0.8f,%0.8f 0 0,0 %0.8f,%0.8f" fill="none" stroke="rgb(102,102,102)" stroke-width="2.0" />\n' % (x1,y1,(x2-x1)/2,arcRadiusY/2,x2,y2)
        else:
          SVGtext += '<path d="M%0.8f,%0.8f A%0.8f,%0.8f 0 0,1 %0.8f,%0.8f" fill="none" stroke="rgb(102,102,102)" stroke-width="2.0" />\n' % (x2,y2,(x2-x1)/2,arcRadiusY/2,x1,y1)
      else:
        SVGtext += '<line x1="%0.8f" y1="%0.8f" x2="%0.8f" y2="%0.8f" stroke="rgb(102,102,102)" stroke-width="2.0" />\n' % (x1,y1,x2,y2)

  # draw out the basepairing interactions
  for interaction in interactionList:
    t = interaction["type"]
    if t[0] == "c" or t[0] == "t":
      i = interaction["firstSequencePosition"]
      j = interaction["secondSequencePosition"]
      p1 = SVG[i]["p"]
      p2 = SVG[j]["p"]
      x1 = p1[0]
      y1 = p1[1]
      x2 = p2[0]
      y2 = p2[1]
      unit = (p2-p1)/np.linalg.norm(p2-p1)  # unit vector from 1 to 2
      perp = np.array([-unit[1],unit[0]])   # unit vector perpendicular to the interaction
      rotation = np.column_stack((unit,perp))
      pair = sequence[i] + sequence[j]

      for a in range(min(i,j)+1,max(i,j)):
        for b in range(0,min(i,j)):
          pairs[a][b] *= 2              # mark this as having been crossed
          pairs[b][a] *= 2              # mark this as having been crossed
        for b in range(max(i,j)+1,len(sequence)):
          pairs[a][b] *= 2              # mark this as having been crossed
          pairs[b][a] *= 2              # mark this as having been crossed

      if pairs[i][j] == 1:
        color = scarlet
      elif pairs[i][j] == 2:
        color = pink
      else:
        color = lightScarlet            # if already crossed, draw lighter

      # if the interaction is on the same horizontal level, use an elliptical arc
      if abs(y1-y2) < 0.01 and abs(i-j) > 1 and somethingBetween(SVG,i,j):
        if y1 > (topStrandY + bottomStrandY)/2:
          sweepFlag = 0
          arcShift = np.array([0,-arcRadiusY])
          SVGtext += '<path d="M%0.8f,%0.8f A%0.8f,%0.8f 0 0,0 %0.8f,%0.8f" fill="none" stroke="%s" stroke-width="1.0" />\n' % (x1,y1,(x2-x1)/2,arcRadiusY,x2,y2,color)
        else:
          sweepFlag = 1
          arcShift = np.array([0,arcRadiusY])
          SVGtext += '<path d="M%0.8f,%0.8f A%0.8f,%0.8f 0 0,1 %0.8f,%0.8f" fill="none" stroke="%s" stroke-width="1.0" />\n' % (x2,y2,(x2-x1)/2,arcRadiusY,x1,y1,color)

        # the following code works for everything but GC cWW and AU cWW.  We can do those later.
        if t[1] == t[2]:
          SVGtext += drawLWSymbol(t[0],t[1],(p1+p2)/2+arcShift,color,rotation)   # position relative to midpoint of interaction line
        else:
          SVGtext += drawLWSymbol(t[0],t[1],(p1+p2)/2-3*unit+arcShift,color,-rotation)   # position relative to midpoint of interaction line
          SVGtext += drawLWSymbol(t[0],t[2],(p1+p2)/2+3*unit+arcShift,color,rotation)
      elif t == "cWW" and (pair == "CG" or pair == "GC"):
        # draw a double line
        a = p1 + 0.5*GCseparation*perp
        b = p2 + 0.5*GCseparation*perp
        SVGtext += '<line x1="%0.8f" y1="%0.8f" x2="%0.8f" y2="%0.8f" stroke="%s" stroke-width="1.0" />\n' % (a[0],a[1],b[0],b[1],color)
        a = p1 - 0.5*GCseparation*perp
        b = p2 - 0.5*GCseparation*perp
        SVGtext += '<line x1="%0.8f" y1="%0.8f" x2="%0.8f" y2="%0.8f" stroke="%s" stroke-width="1.0" />\n' % (a[0],a[1],b[0],b[1],color)
      else:
        SVGtext += '<line x1="%0.8f" y1="%0.8f" x2="%0.8f" y2="%0.8f" stroke="%s" stroke-width="1.0" />\n' % (x1,y1,x2,y2,color)
        if t[1] == t[2]:
          if not (t == "cWW" and (pair == 'AU' or pair == 'UA')):
            SVGtext += drawLWSymbol(t[0],t[1],(p1+p2)/2,color,rotation)
        else:
          SVGtext += drawLWSymbol(t[0],t[1],(p1+p2)/2-3*unit,color,-rotation)   # position relative to midpoint of interaction line
          SVGtext += drawLWSymbol(t[0],t[2],(p1+p2)/2+3*unit,color,rotation)

  for i in range(0,len(sequence)-1):
    if not i == s and not i+1 == s:
      p1 = SVG[i]["p"]
      p2 = SVG[i+1]["p"]
      x1 = p1[0]
      y1 = p1[1]
      x2 = p2[0]
      y2 = p2[1]
      SVGtext += '<line x1="%0.8f" y1="%0.8f" x2="%0.8f" y2="%0.8f" stroke="rgb(0%%, 0%%, 0%%)" stroke-width="0.5" />\n' % (x1,y1,x2,y2)

  # draw out the circles then letters then numbers
  for i in range(0,len(sequence)):
    if not sequence[i] == "*":
      p = SVG[i]["p"]
      x = p[0]
      y = p[1]
      c = sequence[i]
      SVGtext += '<circle cx="%0.8f" cy="%0.8f" r="%0.8f" stroke="none" stroke-width="1.0" fill="rgb(100%%, 100%%, 100%%)"/>\n' % (x,y,radius)
      SVGtext += '<circle cx="%0.8f" cy="%0.8f" r="%0.8f" stroke="rgb(0%%, 0%%, 0%%)" stroke-width="1.0" fill="none"/>\n' % (x,y,radius)
      SVGtext += '<text x="%0.8f" y="%0.8f" text-anchor="middle" font-family="Verdana" font-size="7.5" fill="rgb(0%% 0%%, 0%%)" >%s</text>\n' % (x,y+textShift,c)
    if "*" in sequence:
      if i < s:
        SVGtext += '<text x="%0.8f" y="%0.8f" text-anchor="middle" font-family="Verdana" font-size="7" fill="rgb(0%% 0%%, 0%%)" >%d</text>\n' % (x,y-numberShiftTop,numbers[i])
      elif i > s:
        SVGtext += '<text x="%0.8f" y="%0.8f" text-anchor="middle" font-family="Verdana" font-size="7" fill="rgb(0%% 0%%, 0%%)" >%d</text>\n' % (x,y+numberShiftBottom,numbers[i])

  SVGtext = '<svg width="%d" height="70" version="1.1" xmlns="http://www.w3.org/2000/svg">\n\n' % (width+leftX*2) + SVGtext
  SVGtext = '"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">\n' + SVGtext
  SVGtext = '<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" ' + SVGtext
  SVGtext = '<?xml version="1.0" encoding="UTF-8"?>\n' + SVGtext

  SVGtext += "</svg>"


  return SVGtext