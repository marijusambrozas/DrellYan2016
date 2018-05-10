import time
import sys

class Generator:
	def __init__( self ):
		self.TIME = time.strftime('%Y%m%d', time.localtime(time.time()))
		self.alpha0 = 0 # -- for PI contribution -- #
		self.alphaEff = 0.007756146746
		self.pQCDOrder = 1 # -- 0=LO, 1=NLO, 2=NNLO -- #
		self.EWControl = 0 # -- 0: include all EWK correction -- #

		# -- kinematic cuts. default = no restriction -- #
		self.massMin = 0
		self.massMax = 9999
		self.diRapMin = -2000
		self.diRapMax = 2000
		self.diPtMin = 0
		self.diPtMax = 100000
		self.leadLepPtMin = 0
		self.leadLepPtMax = 100000
		self.leadLepAbsEtaMin = 0
		self.leadLepAbsEtaMax = 100
		self.subLepPtMin = 0
		self.subLepPtMax = 100000
		self.subLepAbsEtaMin = 0
		self.subLepAbsEtaMax = 100

		self.PDFSet = "NNPDF30_nnlo_as_0118"

		self.relUnc = 0.1

		self.nCore = 24

		self.Tag = "" # -- used for name -- #

	def TurnOnPI( self ):
		self.alpha0 = 0.007297352568

	def SetAcceptanceDY2016( self ):
		self.leadLepPtMin = 28
		self.subLepPtMin = 17
		self.leadLepAbsEtaMax = 2.4
		self.subLepAbsEtaMax = 2.4

	def SetpQCDOrder( self, value ):
		self.pQCDOrder = value

	def SetnCore( self, value ):
		self.nCore = value

	def SetEWControl( self, value ):
		self.EWControl = value

	def SetPDFSet( self, value ):
		self.PDFSet = value

	def SetMassRange( self, minValue, maxValue ):
		self.massMin = minValue
		self.massMax = maxValue

	def SetDiRapRange( self, minValue, maxValue ):
		self.diRapMin = minValue
		self.diRapMax = maxValue

	def SetDiPtRange( self, minValue, maxValue ):
		self.diPtMin = minValue
		self.diPtMax = maxValue

	def SetTag( self, value ):
		self.Tag = value

	def GenerateParameterInput( self ):
		fileName_param = self.MakeFileName( 'param' )
		scale = int( (self.massMin + self.massMax) / 2.0 )

		f = open(fileName_param, "w")
		f.write(
"""=============================================
'CMS collision energy (GeV)    = ' 13000d0
=============================================
'Factorization scale  (GeV)    = ' {_scale:0.0f}.0d0
'Renormalization scale  (GeV)  = ' {_scale:0.0f}.0d0
=============================================
'Z production (pp=1,ppbar=2)   = ' 1
=============================================
Set to Alpha QED for incoming photon to 0 to turn off photon-induced (photon PDF dependent) channels
'Alpha QED for incoming photon = ' {_alpha0}d0
'Alpha QED effective           = ' {_alphaEff}d0
'Fermi constant (1/Gev^2)      = ' 1.16637d-5
=============================================
'Lepton mass (GeV)             = ' 1.05d-1
'W mass (GeV)                  = ' 80.403d0
'W width (GeV)                 = ' 2.141d0
'Z mass (GeV)                  = ' 91.1876d0
'Z width (GeV)                 = ' 2.4952d0
'Top mass (GeV)                = ' 172.9d0
'Higgs mass (GeV)              = ' 125d0
=============================================
Only QED corrections is on if the input scheme is manual
Input scheme: 0. Manual input; 1. Gmu scheme; 2. AlphaMz scheme; 3. alpha0 scheme
'Which input scheme:           = ' 1
'sin^2(theta_w)                = ' 0.22255d0
'up quark charge               = ' 0.6666667d0
'down quark charge             = ' -0.3333333d0
'lepton charge                 = ' -1d0
'up quark vector coupling      = ' 0.4091d0
'down quark vector coupling    = ' -0.7045d0
'lepton vector coupling        = ' -0.11360d0
'up quark axial coupling       = ' -1d0
'down quark axial coupling     = ' 1d0
'lepton axial coupling         = ' 1d0
=============================================
Vegas Parameters
'Relative accuracy (in %)           = ' {_relUnc:.1f}d0
'Absolute accuracy                  = ' 0d0
'Number of calls per iteration      = ' 1000000
'Number of increase calls per iter. = ' 500000
'Maximum number of evaluations      = ' 200000000
'Random number seed for Vegas       = ' 11
=============================================
'QCD Perturb. Order (0=LO, 1=NLO, 2=NNLO) = ' {_pQCDOrder:.0f}
'EW Perturb. Order (0=LO, 1=NLO)    = ' 1
'Z pole focus (1=Yes, 0=No)	= ' 0
'EW control (leave 0 to keep all on) = ' {_EWControl:.0f} 
'Turn off photon (1=Yes, 0=No, disabled if weak corr. is on) = ' 0
=============================================
'Lepton-pair invariant mass minimum = ' {_massMin:.0f}.0d0
'Lepton-pair invariant mass maximum = ' {_massMax:.0f}.0d0
'Transverse mass minimum           = ' 0d0
'Transverse mass maximum           = ' 100000d0
'Z pT minimum                       = ' {_diPtMin:.0f}d0
'Z pT maximum                       = ' {_diPtMax:.0f}d0
'Z rapidity minimum                 = ' {_diRapMin}d0
'Z rapidity maximum                 = ' {_diRapMax}d0
'Lepton pT minimum                  = ' 0d0
'Lepton pT maximum                  = ' 100000d0
'Anti-lepton pT minimum             = ' 0.0d0
'Anti-lepton pT maximum             = ' 100000d0
'pT min for softer lepton           = ' {_subLepPtMin:.0f}.0d0
'pT max for softer lepton           = ' {_subLepPtMax:.0f}d0
'pT min for harder lepton           = ' {_leadLepPtMin:.0f}.0d0
'pT max for harder lepton           = ' {_leadLepPtMax:.0f}d0
Taking absolute value of lepton pseudorapidity?
'(yes = 1, no = 0) 		    = ' 1
'Lepton pseudorapidity minimum      = ' 0d0
'Lepton pseudorapidity maximum      = ' 100.0d0
Taking absolute value of anti-lepton pseudorapidity?
'(yes = 1, no = 0) 		    = ' 1
'Anti-lepton pseudorapidity minimum = ' 0d0
'Anti-lepton pseudorapidity maximum = ' 100.0d0
Taking absolute value of soft lepton pseudorapidity?
'(yes = 1, no = 0)                  = ' 1
'Softer lepton pseudorapidity min   = ' {_subLepAbsEtaMin}d0 
'Softer Lepton pseudorapidity max   = ' {_subLepAbsEtaMax}d0
Taking absolute value of hard lepton pseudorapidity?
'(yes = 1, no = 0)                  = ' 1
'Harder lepton pseudorapidity min   = ' {_leadLepAbsEtaMin}d0
'Harder Lepton pseudorapidity max   = ' {_leadLepAbsEtaMax}d0
PHOTON RECOMBINATION-----------------------------
'DeltaR sep. for photon recomb.     = ' 0.1d0
'Minimum pT for observable photon   = ' 10d0
'Maximum eta for observable photon  = ' 2.5d0
PHOTON CUTS--------------------------------------
'Minimum Number of Photon           = ' 0
'Maximum Number of Photon           = ' 1
JET DEFINITION-------------------------------
Jet Algorithm & Cone Size ('ktal'=kT algorithm, 'aktal'=anti-kT algorithm, 'cone'=cone)
'ktal, aktal or cone		    = ' ktal
'Jet algorithm cone size (deltaR)   = ' 0.4d0
'DeltaR separation for cone algo    = ' 1.3
'Minimum pT for observable jets     = ' 20d0
'Maximum eta for observable jets    = ' 4.5d0
JET CUTS--------------------------------------
'Minimum Number of Jets		    = ' 0
'Maximum Number of Jets		    = ' 2
'Min. leading jet pT                = ' 0d0
ISOLATION CUTS-------------------------------
'Lep-Anti-lep deltaR minimum        = ' 0.0d0
'Lep-Anti-lep deltaPhi min	    = ' 0.0d0
'Lep-Anti-lep deltaPhi max	    = ' 4.0d0
'Lep-Jet deltaR minimum             = ' 0.0d0
'Lep-Photon deltaR minimum          = ' 0.0d0
=============================================
Cut on Z rapidity for well-defined Collins-Soper Angles at pp Collider
'Z rapidity cutoff for CS frame     = ' 0.0d0
=============================================
(See manual for complete listing)
'PDF set =                        ' '{_PDFSet}'
'Turn off PDF error (1=Yes, 0=No)    = ' 0
(Active for MSTW2008 only, if PDF error is on:)
(Compute PDF+as errors: 1; just PDF errors: 0)
'Which alphaS                       = ' 0
(Active for MSTW2008 only; 0: 90 CL for PDFs+alphas, 1: 68 CL)
'PDF+alphas confidence level        = ' 1
=============================================""".format( 
		_scale=scale, _alpha0=self.alpha0, _alphaEff=self.alphaEff, 
		_relUnc=self.relUnc, _pQCDOrder=self.pQCDOrder, _EWControl=self.EWControl, 
		_massMin=self.massMin, _massMax=self.massMax, _diPtMin=self.diPtMin,_diPtMax=self.diPtMax,
		_diRapMin=self.diRapMin,_diRapMax=self.diRapMax,
		_subLepPtMin=self.subLepPtMin, _subLepPtMax=self.subLepPtMax,
		_leadLepPtMin=self.leadLepPtMin, _leadLepPtMax=self.leadLepPtMax,
		_subLepAbsEtaMin=self.subLepAbsEtaMin, _subLepAbsEtaMax=self.subLepAbsEtaMax,
		_leadLepAbsEtaMin=self.leadLepAbsEtaMin, _leadLepAbsEtaMax=self.leadLepAbsEtaMax,
		_PDFSet=self.PDFSet )
		)
		f.close()

	def GenerateBatchJobScript( self ):
		fileName_param = self.MakeFileName( 'param' )
		fileName_hist = self.MakeFileName( 'hist' )
		fileName_output = self.MakeFileName( 'output' )
		dirName = self.MakeFileName( 'dir' )

		orderName = "LO"
		if self.pQCDOrder == 1: orderName = "NLO"
		elif self.pQCDOrder == 2: orderName = "NNLO"

		fileName_script = self.MakeFileName( 'script' )
		f = open(fileName_script, "w")
		f.write(
"""#!/bin/bash

#########################################################
# -- qsub commands: #$ + command (details: man qsub) -- #
#########################################################
# -- shell used for executing the script -- #
#$ -S /bin/sh

# -- combine standard output & error file -- #
#$ -j y

# -- send the mail when the job is aborted or ended -- #
#$ -m ae -M kplee@cern.ch

# -- stay in the directory where qsub command is executed -- #
#$ -cwd

cwd=$(pwd)

export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# -- CMSSW enviornment -- #
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_9_2_0
cmsenv

cd ${{cwd}}

echo "run: ./local_run.sh z {_dirName} {_fileName_param} {_fileName_hist} {_fileName_output} ../ {_nCore}"
./local_run.sh z {_dirName} {_fileName_param} {_fileName_hist} {_fileName_output} ../ {_nCore}

echo "run: ./finish.sh {_dirName} {_orderName}.{_fileName_output}"
./finish.sh {_dirName} {_orderName}.{_fileName_output}

echo "job is completed"

# -- &>log: "Invalid null command" Error occurs. please use >&log. -- #

# -- PLEASE ENTER AFTER THE LAST LINE! ... IF YOU DON'T, LAST LINE WILL NOT BE EXECUTED BY BATCH JOB -- # 
""".format(
		_dirName=dirName, _fileName_param=fileName_param,
		_fileName_hist=fileName_hist, _fileName_output=fileName_output,
		_nCore=self.nCore, _orderName=orderName)
		)
		
		print "[%s is created]" % fileName_script
		f.close()

		return fileName_script


	# -- to have common format for the file name -- #
	def MakeFileName( self, fileType ):
		List_fileType = ["param", "hist", "script", "dir", "output", "qsub"]
		if fileType not in List_fileType:
			print "fileType = %s is not the correct one. check the details"
			sys.exit()

		fileName = ""
		fileName_base = "v%s_%s" % (self.TIME, self.Tag)

		if fileType == "dir":
			fileName = fileName_base
		elif fileType == "output": fileName = "%s.dat" % fileName_base
		elif fileType == "script": fileName = "%s.sh" % fileName_base
		elif fileType == "qsub": fileName = "%s_fullqSub.sh" % fileName_base
		elif fileType == "param" or fileType == "hist":
			fileName = "%s_%s.txt" % (fileName_base, fileType )

		return fileName

class GeneratorFor2D:
	def __init__( self ):
		self.generator = Generator()
		# -- common settings -- #
		self.generator.SetpQCDOrder( 1 )
		self.generator.SetEWControl( 7 ) # -- 7: QED off -- #
		self.generator.SetPDFSet( 'NNPDF31_nnlo_as_0118' )
		self.generator.SetnCore( 20 )

		# -- edges -- #
		self.List_diRapBinEdge_M15to64 = [0, 0.7, 1.1, 1.9, 2.4, 1000] 
		self.List_diPtBinEdge_M15to64 = [0, 20, 30, 35, 40, 45, 50, 60, 90, 200, 1000]

		self.List_diRapBinEdge_M64to106 = [0, 0.7, 1.9, 1000.0] 
		self.List_diPtBinEdge_M64to106 = [0, 20, 30, 35, 40, 45, 50, 60, 90, 200, 1000]

		self.List_diRapBinEdge_M106to3000 = [0, 0.7, 1000.0] 
		self.List_diPtBinEdge_M106to3000 = [0, 20, 100, 1000]

		self.List_fileName_script = []

	def GenerateInputs_All( self ):
		self.GenerateInputs_EachMassRange( 15, 64, self.List_diRapBinEdge_M15to64, self.List_diPtBinEdge_M15to64 )
		self.GenerateInputs_EachMassRange( 64, 106, self.List_diRapBinEdge_M64to106, self.List_diPtBinEdge_M64to106 )
		self.GenerateInputs_EachMassRange( 106, 3000, self.List_diRapBinEdge_M106to3000, self.List_diPtBinEdge_M106to3000 )


	def GenerateInputs_EachMassRange( self, massMin, massMax, List_diRapBinEdge, List_diPtBinEdge ):
		self.generator.SetMassRange( massMin, massMax )

		diPtMin = List_diPtBinEdge[0]
		diPtMax = List_diPtBinEdge[-1]
		self.generator.SetDiPtRange( diPtMin, diPtMax )

		for i_diRap in range( 0, len(List_diRapBinEdge)-1 ):
			diRapMin = List_diRapBinEdge[i_diRap]
			diRapMax = List_diRapBinEdge[i_diRap+1]
			self.generator.SetDiRapRange( diRapMin, diRapMax )

			tag = "M%.0lfto%.0lf_DiRapBin%02d" % (massMin, massMax, i_diRap)
			self.generator.SetTag( tag )

			self.generator.GenerateParameterInput()

			self.GenerateInput_HistogramInput( tag, massMin, massMax, diRapMin, diRapMax, List_diPtBinEdge )

			fileName_script = self.generator.GenerateBatchJobScript()
			self.List_fileName_script.append( fileName_script )

	def GenerateInput_HistogramInput( self, tag, massMin, massMax, diRapMin, diRapMax, List_diPtBinEdge ):

		fileName_diPtBin = self.generator.MakeFileName( 'hist' )
		fileName_diPtBin = fileName_diPtBin.replace( '_hist', '_diPtBin' )
		f_diPtBin = open(fileName_diPtBin, "w" )
		for binEdge in List_diPtBinEdge:
			f_diPtBin.write('%.1lf\n' % binEdge)
		f_diPtBin.close()

		fileName_hist = self.generator.MakeFileName( 'hist' );
		f = open(fileName_hist, "w")
		f.write(
"""HISTOGRAMS------(Order of Histograms Can Not Be Changed)----------------------
Name (DO NOT CHANGE)    Num Bins (<30)	Lower Bound	Upper Bound	Write Out (1=hist, 2=cuml, 3=both1&2, 4=rev-cuml, 5=both1&4, 0=none)
'1.  Z/W pT           ' {_fileName_diPtBin}		0d0		250d0		1
'2.  Z/W eta          ' 20		{_diRapMin:.0f}d0		{_diRapMax:.0f}d0		1
'3.  Q_ll inv mass    ' 20		{_massMin:.0f}d0		{_massMax:.0f}d0		1
'4.  l-/lep. pT       ' 25		0d0		100d0		1
'5.  l-/lep. eta      ' 24		-9.6d0		9.6d0		1
'6.  l+/neut. pT      ' 25		0d0		100d0		1
'7.  l+/neut. eta     ' 24		-9.6d0		9.6d0		1
'8.  jet 1 pT         ' 20		0d0		100d0		1
'9.  jet 1 eta        ' 20		-5d0		5d0		1
'10. jet 2 pT         ' 20		0d0		100d0		1
'11. jet 2 eta        ' 20		-5d0		5d0		1
'12. photon pT        ' 20              0d0             100d0           1
'13. photon eta       ' 20              -5d0            5d0             1
'14. beam thrust      ' 20              0d0             100d0           1
'15. dR_l-l+          ' 20		0d0		5d0		1
'16. dR_j1,l-	      '	20		0d0		5d0		1
'17. dR_j1,l+	      '	20		0d0		5d0		1
'18. dR_j2,l-	      '	20		0d0		5d0		1
'19. dR_j2,l+	      '	20		0d0		5d0		1
'20. dR_j1j2          ' 20		0d0		5d0		1
'21. dR_phot,l+-      ' 20              0d0             5d0             1
'22. H_T	      ' 20		0d0		200d0		1
'23. Mass_T           ' 20		0d0		1000d0		1
'24. A_FB vs Q_ll     ' 10              {_massMin:.0f}d0            {_massMax:.0f}d0           0
====================================================================================
Moments (A_0, A_1, A_2) related to Collins-Soper Angles
'25. A_i vs Z pT      ' 10		0d0		100d0		1
'26. phi (CS Frame)   '	10		-3.14159265d0	3.14159265d0	1
'27. cos(theta) (CS)  '	10		-1d0		1d0		1
'28. dPhi_l-l+        ' 10		0d0		3.1415927d0	1
====================================================================================
Smoothing parameters
'Level (0 = none)     '  2
'Bin fraction (< 0.5) ' 1d-1
Method of combining iterations (0 = more reliable err. estimation; 1 = more consistent with tot. x-section)
'Method choice = ' 0
Histogram bin display (0 = bin central value, -1 = bin low edge, 1 = bin upper edge)
'Display option =' 0
""".format(
		_fileName_diPtBin="'"+fileName_diPtBin+"'",
		_diRapMin=diRapMin, _diRapMax=diRapMax,
		_massMin=massMin, _massMax=massMax )
		)
		f.close()

	def GenerateQsubScript( self ):
		fileName = self.generator.MakeFileName("qsub")
		f = open( fileName, "w")
		f.write("#!bin/bash\n\n")
		for fileName_script in self.List_fileName_script:
			f.write("qsub -V %s\n" % fileName_script)

		f.write("\n")
		f.write('echo "finished"\n')

# -- main part -- #
generatorFor2D = GeneratorFor2D()
generatorFor2D.GenerateInputs_All()
generatorFor2D.GenerateQsubScript()

