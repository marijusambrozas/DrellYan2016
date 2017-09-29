import sys, json

Lumi = 35867.060
print "Luminosity: ", Lumi

JSONFileName = sys.argv[1]
print "[JSONFileName]", JSONFileName
with open(JSONFileName) as SampleInfo_File:
	SampleInfo = json.load(SampleInfo_File)

for Sample in SampleInfo["Sample"]:
	IsMC = Sample["IsMC"]
	if IsMC:
		XSec = Sample["CrossSection"]
		SumWeights = Sample["SumWeights"]
		NormFactor = (Lumi * XSec) / SumWeights
		print "\n[%s (IsMC = %d)] (XSec, Sum(weights), Norm.factor) = (%e, %lf, %e)" % (Sample["Tag"], Sample["IsMC"], XSec, SumWeights, NormFactor)
	else:
		print "\n[%s]" % (Sample["Tag"])

	for Dir in Sample["List_Dir"]:
		print "\t%s" % Dir