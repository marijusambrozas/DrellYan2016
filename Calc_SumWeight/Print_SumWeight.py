import os, sys
from ROOT import TFile, TH1D

if len(sys.argv) < 2:
	print "Usage: python %s [root file name]" % sys.argv[0]
	sys.exit()

InputFileName = sys.argv[1]
f_input = TFile(InputFileName)

for key in f_input.GetListOfKeys():
	if key.IsFolder():
		DirName = key.GetName()

		for subkey in key.ReadObj().GetListOfKeys():
			if subkey.GetName() == "h_SumWeight":
				h = subkey.ReadObj()
				SumWeight = h.GetBinContent(1)
				print "[%s] Sum(weight) = %.1lf\n" % (DirName, SumWeight)
				break

print "Finished"


