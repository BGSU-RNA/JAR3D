
# about jar3dalign_2015-04-03.jar
# Main-Class: edu.bgsu.rna.jar3d.JAR3DCorrespondences
#       String fastaFileName = args[0];
#       String modelDirPath = args[1];
#       String modelName = args[2];
#       int rotation = Integer.parseInt(args[3]);
#       String outputFileName = args[4];

# http://varna.lri.fr/index.php?lang=en&css=varna&page=command
# java -cp VARNAvX-Y.jar fr.orsay.lri.varna.applications.VARNAcmd [-i inputFile|-sequenceDBN XXX -structureDBN YYY] -o outFile [opts]
# java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd -sequenceDBN ACGU -structureDBN .... -o test.png

import os
from correspondenceToSVG import correspondenceToInteractionList
from correspondenceToSVG import interactionListToSVG

JAR3DPath = "C:\\Users\\zirbel\\Documents\\JAR3D\\"
JAR3DAligner = "jar3dalign_2015-04-03.jar"
viroidPath = "C:\\Users\\zirbel\\Documents\\JAR3D\\viroid\\"
viroidDataFile = viroidPath + "WT sequences plus model IDs.csv"
fastaFile = viroidPath + "temp.fasta"

with open(viroidDataFile,'r') as fh:
    viroidData = fh.readlines()

#viroidData = viroidData[2:]

html = "<html><body>"

for data in viroidData:
    print(data)
    fields = data.split(",")
    loopNum = fields[0]
    sequence = fields[1]

    if (loopNum == "1" or loopNum == "27"):
        continue

    numbers = []
    if "*" in sequence:            # internal loop
        s = sequence.find("*")
        for i in range(0,len(sequence)):
            if i < s:
                numbers.append(int(fields[2])+i)
            elif i > s:
                numbers.append(int(fields[3])+i-s-1)
            else:
                numbers.append(0)    # * should not have a number
    else:                          # hairpin loop
        a = int(fields[2])
        b = int(fields[3])
        if a < b:                  # regular hairpin
            for i in range(0,len(sequence)):
                numbers.append(int(fields[2])+i)
        else:                      # unusual case of a hairpin at the end of a circular RNA
            for i in range(0,len(sequence)):
                if len(sequence) - i > b:
                    numbers.append(a+i)
                else:
                    numbers.append(b-len(sequence)+i)

    with open(fastaFile,'w') as fh:
        fh.write("> header\n")          # note that the header is necessary!
        fh.write(sequence+"\n")

    for m in range(1,len(fields)-3):    # counter through model numbers 1 to 10
        modelName = fields[m+3].replace("\n","").replace(" ","")

        html += "Loop %s Model %d " % (loopNum, m)
        html += '<a href="http://rna.bgsu.edu/rna3dhub/motif/view/%s">%s</a><br>' % (modelName,modelName)

        if "IL" in modelName:
            maxRotation = 1
            type = "IL"
        else:
            maxRotation = 0
            type = "HL"

        bestScore = -float('inf')
        bestRotation = 0
        for rotation in range(0,maxRotation+1):
            align = "java -jar " + JAR3DPath + JAR3DAligner + " "
            align += fastaFile + " "
            align += JAR3DPath + type + "\\1.18\\lib\\" + " "
            align += modelName + " "
            align += str(rotation) + " "
            correspondenceFile = viroidPath + "correspondences_" + str(rotation) + ".txt"

            align += correspondenceFile
            print("Loop " + loopNum + " " + str(m) + " " + fastaFile + " " + modelName + " " + str(rotation))
            os.system(align)

            if type == "IL":
                with open(correspondenceFile,'r') as fh:
                    correspondenceData = fh.readlines()

                for line in correspondenceData:
                    if "Sequence_1 has_score " in line:
                        lineFields = line.split(" ")
                        currentScore = float(lineFields[2])

                        if currentScore > bestScore:
                            bestRotation = rotation
                            bestScore = currentScore

        rotation = bestRotation
        print("Best rotation is "+str(rotation))
        correspondenceFile = viroidPath + "correspondences_" + str(rotation) + ".txt"

        imageFile = viroidPath + "Loop_" + loopNum + "_Model_" + "%02d" % m + "_" + modelName + "_" + str(rotation) + ".svg"

        html += '<img src="%s"><br>' % ("Loop_" + loopNum + "_Model_" + "%02d" % m + "_" + modelName + "_" + str(rotation) + ".svg")

        interactionList = correspondenceToInteractionList(JAR3DPath + type + "\\1.18\\lib",modelName,correspondenceFile)
        if rotation == 1:
            s = sequence.find("*")
            map = range(s+1,len(sequence)) + [s] + range(0,s)   # works for internal loops; junctions harder
            for i in range(0,len(interactionList)):
                interactionList[i]["firstSequencePosition"] = map[interactionList[i]["firstSequencePosition"]]
                interactionList[i]["secondSequencePosition"] = map[interactionList[i]["secondSequencePosition"]]

        SVGtext = interactionListToSVG(sequence,interactionList,numbers)

#        print(SVGtext)

        with open(imageFile,'w') as fh:
            fh.write(SVGtext)

        html += "</body></html>"

        with open(viroidPath + "alldiagrams.html",'w') as fh:
            fh.write(html)

#        VARNA = correspondenceToVARNA(JAR3DPath + type + "\\1.18\\lib",modelName,correspondenceFile,VARNAFile)
#        os.system(VARNA)

