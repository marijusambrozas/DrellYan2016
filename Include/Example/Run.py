import os, sys, time, json
from os.path import expanduser

if len(sys.argv) < 2:
	print "Please put the json file name as first arguement\n\tex>python Run.py BatchJobInfo.json\n"
	sys.exit()

JSONFileName = sys.argv[1]
print "[JSONFileName]", JSONFileName
with open(JSONFileName) as BatchJobInfo_File:
	BatchJobInfo = json.load(BatchJobInfo_File)

def MakeScript_Sub(_OutDir, _List_JobType, _TIME):
	f = open("%s/script_qsub_ALL.sh" % _OutDir, "w")
	f.write( "#!bin/bash\n" )
	f.write( "cwd2=$(pwd)\n\n" )
	# f.write( "cd %s\n" % (_OutDir) )
	for JobType in _List_JobType:
		SampleName = JobType["Tag"]
		f.write( "cd %s; source qsub_all.sh\n" % (SampleName) )
		f.write( "cd ../\n")

	f.write( "\n\n" )

	f.write( 'cd ${cwd2}\n' )
	f.write( 'echo "full submission is finished"\n' )
	f.close()

	print "cd %s; source script_qsub_ALL.sh" % _OutDir

def MakeScript_GetOutput(_OutDir, _List_JobType, _TIME):
	f = open("%s/script_GetOutput.sh" % _OutDir, "w")
	f.write( "#!bin/bash\n" )
	f.write( "cwd2=$(pwd)\n\n" )
	# f.write( "cd %s\n" % (_OutDir) )

	for JobType in _List_JobType:
		SampleName = JobType["Tag"]
		f.write( "cd %s; source hadd_all.sh\n" % (SampleName) )
		f.write( "cp *.root ${cwd2}\n" )
		f.write( "cd ../\n")
		f.write( "\n" )

	f.write( "\n\n" )

	f.write( 'cd ${cwd2}\n' )
	f.write( 'echo "finished"\n' )
	f.close()

##############
# -- main -- #
##############

# -- output directory -- #
AnalyzerPATH = os.environ['KP_ANALYZER_PATH']
TIME = time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time()))
OutDir = "Local/v%s_%s" % (TIME, BatchJobInfo["CodeName"].split('.cxx')[0] )
# OutDir = os.path.expanduser( OutDir )
os.makedirs( OutDir ) # -- make recursively -- #

# -- convert path for output directory to absolute path -- #
OutDirAbsPath = os.path.abspath( OutDir )

# -- convert path for code to absolute path -- #
CodeAbsPath = os.path.abspath( BatchJobInfo["CodeName"] )

# -- make script -- #
MakeScript_Sub( OutDir, BatchJobInfo["List_JobType"], TIME )
MakeScript_GetOutput( OutDir, BatchJobInfo["List_JobType"], TIME )

sys.path.append( "%s/Include" % AnalyzerPATH )
from Split_BatchJobs import SplitJobs

for JobType in BatchJobInfo["List_JobType"]:
	SampleName = JobType["Tag"]
	nJob = (int)(JobType["nJob"])

	opts = dict()
	opts['--code'] = CodeAbsPath
	opts['--sample'] = SampleName
	opts['--njob'] = nJob
	opts['--lumi'] = (float)(BatchJobInfo["Lumi"])
	opts['--outdir'] = OutDirAbsPath
	opts['--queue'] = BatchJobInfo["Queue"]

	split = SplitJobs( opts )
	split.CreateWorkSpace()

print "+"*100
print "[Job submission] cd %s; source script_qsub_ALL.sh" % OutDir
print "+"*100