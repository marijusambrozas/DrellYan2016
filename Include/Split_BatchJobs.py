import os, sys, json
from getopt import gnu_getopt as getopt

def usage():
    print sys.argv[0], " : split the jobs"
    print "  Mandatory options :"
    print "   --code CODE.cxx                  C++ Code file name: should be absolute path"
    print "   --sample SAMPLENAME    		   Sample Name"
    print "   --njob N                9         Total number of jobs"
    print "   --lumi LUMI                      Integrated lumionsity in pb"
    print "   --outdir Dir                     Directory where outputs are stored"
    print "  Optional options :"
    print "   --queue queueName                      queue name (default: fastq or bigq)"
    print "   --json JSONFILENAME    		   json file name containing all sample information"
    sys.exit()

# Parse arguments
if len(sys.argv) < 2: usage()
try:
    opts, args = getopt(sys.argv[1:], 'n', ["code=", "sample=", "njob=", "lumi=", "queue=", "outdir=", "--json"])
    opts = dict(opts)
except:
    print "!!! Error parsing arguments"
    usage()


class SplitJobs:
	def __init__(self, _opts):
		self.CodeFullPath = _opts['--code']
		self.nJob = int(_opts['--njob'])
		self.Lumi = float(_opts['--lumi'])
		self.OutDir = _opts['--outdir']
		self.CodeName = self.CodeFullPath.split('/')[-1]

		# -- get sample info from input json -- #
		self.IncludePath = os.getenv("KP_INCLUDE_PATH")
		self.JSONFileName = self.IncludePath + "/SampleInfo_Latest.json"
		if '--json' in _opts:
			self.JSONFileName = _opts['--json']
		self.SampleInfo = self.Get_SampleInfo( _opts['--sample'] )

		self.queue = "fastq"
		if os.getenv("HOSTNAME") == "tamsa2.snu.ac.kr":
			self.queue = "bigq"
		if '--queue' in _opts:
			self.queue = _opts['--queue']

		_TYPE = ""
		if self.SampleInfo["IsMC"] == 1: _TYPE = "MC"
		else: _TYPE = "Data"

		print "+" * 100
		print "Create %d jobs to run %s with lumi = %lf on %s -> queue name = %s" % (self.nJob, self.CodeName, self.Lumi, _TYPE, self.queue)
		print "Output directory: %s" % (self.OutDir)
		print "+" * 100
		
		# -- get all list of root files corresponding to this sample -- #
		self.List_ROOTFiles = self.GetListOfROOTFiles()

		# -- get normalization factor if it is MC sample (for data, norm.factor is set as 1) -- #
		if self.SampleInfo["IsMC"] == 1:
			XSec = self.SampleInfo["CrossSection"]
			SumWeights = self.SampleInfo["SumWeights"]

			if XSec < 0 or SumWeights < 0:
				print "Sample info is not complete: (Xsec, Sum(weight)) = (%e, %e). exit()" % (XSec, SumWeights)
				sys.exit()
			self.NormFactor = (self.Lumi * XSec) / SumWeights
		else:
			self.NormFactor = 1

		os.chdir( self.OutDir )


	def Get_SampleInfo( self, SampleName ):
		JSONFileName = self.JSONFileName
		print "[Used json for full sample info.]", JSONFileName
		with open(JSONFileName) as FullSampleInfo_File:
			FullSampleInfo = json.load(FullSampleInfo_File)

		Found = False
		FoundSampleInfo = set()
		for Sample in FullSampleInfo["Sample"]:
			if SampleName == Sample["Tag"]:
				Found = True
				FoundSampleInfo = Sample
				break

		if not Found:
			print "No correspond info. for %s in the %s!" % (SampleName, JSONFileName)
			sys.exit()

		return FoundSampleInfo

	def CreateWorkSpace( self ):
		DirName = "%s" % (self.SampleInfo["Tag"])
		
		List_File_cwd = os.listdir( "." )
		if DirName in List_File_cwd:
			print "%s is already exists!" % (DirName)
			sys.exit()

		os.mkdir( DirName )

		nROOTFile = len(self.List_ROOTFiles)
		if self.nJob > nROOTFile:
			print "nJob > nROOTFile -> nJob is set as same with nROOTFile"
			self.nJob = nROOTFile

		nROOTFilePerJob = int( float(nROOTFile) / float(self.nJob) )
		print "nJob = %d, nROOTFile = %d -> nROOTFilePerJob = %d\n" % (self.nJob, nROOTFile, nROOTFilePerJob)

		List_cmd_qsub = []
		List_cmd_hadd = []
		for i in range(0, self.nJob):
			List_ROOTFilesPerJob = []
			if i == self.nJob-1:
				List_ROOTFilesPerJob = self.List_ROOTFiles[int(i*nROOTFilePerJob):]
			else:
				List_ROOTFilesPerJob = self.List_ROOTFiles[int(i*nROOTFilePerJob):int((i+1)*nROOTFilePerJob)]
			
			# print "List_ROOTFilesPerJob"
			# for rootfile in List_ROOTFilesPerJob:
			# 	print "%s" % (rootfile)
			# print "\n"
			
			self.CreateWorkSpace_PerJob( i, List_ROOTFilesPerJob, List_cmd_qsub, List_cmd_hadd )

		self.CreateScript_qsub( List_cmd_qsub )
		self.CreateScript_hadd( List_cmd_hadd )

	def CreateWorkSpace_PerJob( self, _iter, _List_ROOTFilesPerJob, List_cmd_qsub, List_cmd_hadd ):
		DirName = "Job_%d" % (_iter)
		os.mkdir( "./"+self.SampleInfo["Tag"]+"/"+DirName )

		FileName = "ROOTFileList_%d.txt" % (_iter)
		f_ROOTFileList = open("./"+self.SampleInfo["Tag"]+"/"+DirName+"/"+FileName, "w")
		for rootfilepath in _List_ROOTFilesPerJob:
			f_ROOTFileList.write( rootfilepath )
			f_ROOTFileList.write( "\n" )
		f_ROOTFileList.close()

		cmd_cp = "cp %s ./%s/%s" % (self.CodeFullPath, self.SampleInfo["Tag"], DirName)
		# print cmd_cp
		os.system( cmd_cp )

		cmd_execute = "root -l -b -q '"+ self.CodeName + '++("%s", "%s", %d, %.10e)' % (FileName, self.SampleInfo["Tag"], self.SampleInfo["IsMC"], self.NormFactor) + "'";
		print cmd_execute

		BatchFileName = self.CreateBatchJobScript( _iter, DirName, cmd_execute )

		cmd_cd = "cd ${cwd}/%s" % (DirName)
		cmd_sub = "qsub -V -q %s %s" % (self.queue, BatchFileName)
		if self.queue == "NoBatch":
			cmd_sub = "source %s >&log.txt&" % BatchFileName
		List_cmd_qsub.append( [cmd_cd, cmd_sub] )

		cmd_hadd = "${cwd}/%s/*.root \\" % (DirName)
		List_cmd_hadd.append( cmd_hadd )

	def CreateScript_qsub( self, _List_cmd_qsub ):
		f = open("./"+self.SampleInfo["Tag"]+"/qsub_all.sh", "w")
		f.write( "#!bin/bash\n" )
		f.write( "cwd=$(pwd)\n\n" )
		for cmd_qsub in _List_cmd_qsub:
			cmd_cd = cmd_qsub[0]
			cmd_sub = cmd_qsub[1]

			f.write( cmd_cd )
			f.write( "\n" )
			f.write( cmd_sub )
			f.write( "\n\n" )

		f.write( 'cd ${cwd}\n' )
		f.write( 'echo "finished"\n' )
		f.close()

		print "cd %s; source qsub_all.sh" % (self.SampleInfo["Tag"])

	def CreateScript_hadd( self, _List_cmd_hadd ):
		f = open("./"+self.SampleInfo["Tag"]+"/hadd_all.sh", "w")
		f.write( "#!bin/bash\n" )
		f.write( "cwd=$(pwd)\n\n" )
		f.write( "hadd ROOTFile_%s.root \\\n" % (self.SampleInfo["Tag"]) )
		for cmd_hadd in _List_cmd_hadd:
			f.write( cmd_hadd )
			f.write( "\n" )
		
		f.write( "\n" )
		f.write( 'echo "finished"\n' )
		f.close()

	def CreateBatchJobScript( self, iter, _DirName, cmd_execute ):
		BatchFileName = "%s_%d.sh" % (self.SampleInfo["Tag"], iter)
		f = open("./"+self.SampleInfo["Tag"]+"/"+_DirName+"/"+BatchFileName, "w")
		f.write(
"""#!/bin/bash

#########################################################
# -- qsub commands: #$ + command (details: man qsub) -- #
#########################################################
# -- shell used for executing the script -- #
#$ -S /bin/sh

# -- combine standard output & error file -- #
#$ -j y

## -- send the mail when the job is aborted or ended -- #
##$ -m ae -M kplee@cern.ch

# -- stay in the directory where qsub command is executed -- #
#$ -cwd

cwd=$(pwd)

export SCRAM_ARCH=slc6_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# -- CMSSW enviornment -- #
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc530/cms/cmssw/CMSSW_8_0_26
cmsenv

cd ${{cwd}}

{_cmd_execute}

echo "job is completed"

# -- &>log: "Invalid null command" Error occurs. please use >&log. -- #

# -- PLEASE ENTER AFTER THE LAST LINE! ... IF YOU DON'T, LAST LINE WILL NOT BE EXECUTED BY BATCH JOB -- # 
""".format(_cmd_execute=cmd_execute)
		)
		print "[%s is created]" % BatchFileName
		f.close()

		return BatchFileName

	def GetListOfROOTFiles( self ):
		List_FullPath = []
		BasePath = os.getenv('KP_DATA_PATH')
		for DirName in self.SampleInfo["List_Dir"]:
			List_FullPath.append( "%s/%s" % (BasePath, DirName) )

		List_ROOTFiles = []
		for fullpath in List_FullPath:
			FileList = os.listdir( fullpath )

			for file in FileList:
				if ".root" in file:
					List_ROOTFiles.append( fullpath+"/"+file )

		if len(List_ROOTFiles) == 0:
			print "There is no available root files ... check the directory name"
			print "Path: ", List_FullPath
			sys.exit()

		return List_ROOTFiles

if __name__ == "__main__":
	split = SplitJobs( opts )
	split.CreateWorkSpace()