from ROOT import gROOT, TColor, TROOT, TFile, TH1D, TCanvas, TPad, TLatex
from array import array
import math
import os, sys

class ReadFEWZOutputFor2D:
	def __init__(self):
		self.fileName = ""
		self.List_diRapBinEdge_M15to64 = [0, 0.7, 1.1, 1.9, 2.4, 1000] 
		self.List_diPtBinEdge_M15to64 = [0, 20, 30, 35, 40, 45, 50, 60, 90, 200, 1000]

		self.List_diRapBinEdge_M64to106 = [0, 0.7, 1.9, 1000.0] 
		self.List_diPtBinEdge_M64to106 = [0, 20, 30, 35, 40, 45, 50, 60, 90, 200, 1000]

		self.List_diRapBinEdge_M106to3000 = [0, 0.7, 1000.0] 
		self.List_diPtBinEdge_M106to3000 = [0, 20, 100, 1000]

		self.List_diRapBinEdge = []
		self.List_diPtBinEdge = []

	def SetFileName( self, value ):
		self.fileName = value

	def MakeHist_diPt( self, f_output ):
		self.FindBinEdge()
		histName = self.MakeHistName()
		print "histName: ", histName

		h_xSec = TH1D(histName, "", len(self.List_diPtBinEdge)-1, array("d", self.List_diPtBinEdge) )
		h_relUnc_integ = TH1D(histName+"_relUnc_integ", "", len(self.List_diPtBinEdge)-1, array("d", self.List_diPtBinEdge) )
		h_relUnc_pdf = TH1D(histName+"_relUnc_pdf", "", len(self.List_diPtBinEdge)-1, array("d", self.List_diPtBinEdge) )
		h_relUnc_tot = TH1D(histName+"_relUnc_tot", "", len(self.List_diPtBinEdge)-1, array("d", self.List_diPtBinEdge) )

		List_xSecInfo = self.ReadFEWZOutput_diPtBin()
		if len(List_xSecInfo) != len(self.List_diPtBinEdge)-1:
			print "# cross section values != # diPt bin edges-1"
			sys.exit()

		for i in range(0, len(List_xSecInfo) ):
			i_bin = i+1

			xSecInfo = List_xSecInfo[i]

			binCenter = xSecInfo[0]
			xSec = xSecInfo[1]
			IntegErr = xSecInfo[2]
			pdfErrP = xSecInfo[3]
			# pdfErrM = xSecInfo[4]

			# pdfErr = 0
			# if pdfErrP > pdfErrM:
			# 	pdfErr = pdfErrP
			# else:
			# 	pdfErr = pdfErrM

			pdfErr = pdfErrP

			totErr = math.sqrt( IntegErr*IntegErr + pdfErr*pdfErr )

			h_xSec.SetBinContent( i_bin, xSec )
			h_xSec.SetBinError( i_bin, totErr )
			print "[%02d bin] (xsec, error) = (%lf, %lf)" % (i_bin, xSec, totErr )

			relUnc_integ = IntegErr / xSec;
			relUnc_pdf = pdfErr / xSec;
			relUnc_tot = totErr / xSec;

			h_relUnc_integ.SetBinContent( i_bin, relUnc_integ )
			h_relUnc_integ.SetBinError( i_bin, 0 )

			h_relUnc_pdf.SetBinContent( i_bin, relUnc_pdf )
			h_relUnc_pdf.SetBinError( i_bin, 0 )

			h_relUnc_tot.SetBinContent( i_bin, relUnc_tot )
			h_relUnc_tot.SetBinError( i_bin, 0 )

		f_output.cd()
		h_xSec.Write()
		h_relUnc_integ.Write()
		h_relUnc_pdf.Write()
		h_relUnc_tot.Write()

	def ReadFEWZOutput_diPtBin( self ):
		f = open(self.fileName, "r")
		lines = f.readlines()

		i = 0;
		for line in lines:
			if "Z/W pT" in line:
				i_begin = i+2

			if "Z/W rapidity" in line:
				i_end = i-2
				break

			i += 1

		List_xSecInfo = []
		# -- Z pT part -- #
		for i in range(i_begin, i_end+1):
			print "lines[%d]" % i, lines[i]

			Numbers = []
			for item in lines[i].split():
				if self.IsNumber(item):
					Numbers += [float(item)]

			List_xSecInfo.append( Numbers )

		return List_xSecInfo


	def FindBinEdge( self ):
		if "M15to64" in self.fileName:
			self.List_diRapBinEdge = self.List_diRapBinEdge_M15to64
			self.List_diPtBinEdge = self.List_diPtBinEdge_M15to64
		elif "M64to106" in self.fileName:
			self.List_diRapBinEdge = self.List_diRapBinEdge_M64to106
			self.List_diPtBinEdge = self.List_diPtBinEdge_M64to106
		elif "M106to3000" in self.fileName:
			self.List_diRapBinEdge = self.List_diRapBinEdge_M106to3000
			self.List_diPtBinEdge = self.List_diPtBinEdge_M106to3000

	def MakeHistName( self ):
		items = self.fileName.split("_")
		histName = "%s_%s" % (items[1], items[2].split(".dat")[0])
		return histName

	def IsNumber(self, s):
		try:
			float(s)
			return True
		except ValueError:
			return False

if __name__ == '__main__':
	f_output = TFile("ROOTFile_FEWZ_diPt.root", "RECREATE")

	List_outputFile = []
	for fileName in os.listdir("."):
		if ".dat" in fileName:
			List_outputFile.append( fileName )

	List_outputFile.sort()

	for outputFile in List_outputFile:
		print "processing ", outputFile, "\n"
		reader = ReadFEWZOutputFor2D()
		reader.SetFileName(outputFile)
		reader.MakeHist_diPt( f_output )
