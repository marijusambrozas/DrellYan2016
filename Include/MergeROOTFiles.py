import os
import sys
from ROOT import TFile, TH1D

if len(sys.argv) < 2:
	print "Usage: python MergyROOTFiles.py [output root file name]"
	sys.exit()

OutputFileName = sys.argv[1]

# -- pre-merge binned samples -- #
List_PreMergedSample = ["DYPowheg", "DYMuMu", "DYEE"]


for PreMergedSample in List_PreMergedSample:
	isFound = False
	for File in os.listdir("."):
		if ".root" in File and "ROOTFile_%s_M" % PreMergedSample in File:
			isFound = True
			break

	if isFound:
		cmd_hadd = "hadd ROOTFile_%s.root ROOTFile_%s_M*.root" % (PreMergedSample, PreMergedSample)
		print cmd_hadd
		os.system( cmd_hadd )

############################################

# -- collect the list of root files that will be merged -- #
List_ROOTFile = []
print "List of ROOT files that will be merged: "
for File in os.listdir("."):
	ext = os.path.splitext(File)[-1]
	if ext == ".root":
		print "\t"+File
		List_ROOTFile.append ( File )

if OutputFileName in os.listdir("."):
	print "%s is already exists! ... please check"
	sys.exit()

# -- output root file -- #
f_output = TFile(OutputFileName, "RECREATE")

print "Loop over list of ROOT files ..."
List_ROOTFile.sort()
for ROOTFile in List_ROOTFile:

	# -- make a directory -- #
	DirName = ROOTFile.split(".")[0].split("ROOTFile_")[-1]
	print "\tDirName = %s" % (DirName)
	f_output.cd() # -- return to top directory -- #
	f_output.mkdir( DirName )

	f_input = TFile( ROOTFile )
	f_input.cd()
	# -- loop over object in TFile -- #
	for key in f_input.GetListOfKeys():
		obj = key.ReadObj()

		f_output.cd( DirName ) # -- move to sub directory -- #
		obj.Write()
		
		f_input.cd()

	f_input.Close()
	# cmd_mv = "mv %s ./Local" % (ROOTFile)
	# print "Move: %s -> ./Local" % (ROOTFile)
	# os.system(cmd_mv)

print "Merging root files are finished: Output file = ", f_output.GetName()
f_output.Close()