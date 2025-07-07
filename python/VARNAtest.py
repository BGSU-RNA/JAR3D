
# VARNA options at: http://varna.lri.fr/index.php?lang=en&page=command&css=varna

FASTASequence = "ACGU"

outputFileName = "test.png"

option = {}
option["sequenceDBN"] = FASTASequence
option["structureDBN"] = "."*len(FASTASequence)

allBPs = '(1,4):edge5="h",edge3="wc",stericity="trans";'
allBPs += "(2,3)"

option["auxBPs"] = allBPs
option["rotation"] = "90.0"
option["basenum"] = "#FFFFFF"
option["eror"] = "True"
#option["drawBackbone"] = "false"
#option["annotations"] = "27:anchor=1,type=B"

# (%i,%i):edge5=%s,edge3=%s,stericity=%s;',i,j,edge5,edge3,stericity

VARNA = "java -cp VARNAv3-93.jar fr.orsay.lri.varna.applications.VARNAcmd "
for a in option:
	VARNA += " -" + a + " " + option[a]
VARNA += " -o " + outputFileName

print(VARNA)
