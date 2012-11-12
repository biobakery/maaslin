#!/usr/bin/env python

"""
Author: George Weingart
Description: Wrapper program for maaslin
"""

#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#####################################################################################

__author__ = "George Weingart"
__copyright__ = "Copyright 2012"
__credits__ = ["George Weingart"]
__license__ = "MIT"
__maintainer__ = "George Weingart"
__email__ = "george.weingart@gmail.com"
__status__ = "Development"

from cStringIO import StringIO
import sys,string
import os
import tempfile 
from pprint import pprint
import argparse

######################################################################################
#  Parse input parms                                                                 #
######################################################################################
def read_params(x):
	parser = argparse.ArgumentParser(description='MaAsLin Argparser')
	parser.add_argument('--lastmeta', action="store", dest='lastmeta',nargs='?')
	parser.add_argument('--input', action="store", dest='input',nargs='?')
	parser.add_argument('--output', action="store", dest='output',nargs='?')
	parser.add_argument('--zip_file', action="store", dest='zip_file',nargs='?')
	parser.add_argument('--alpha', action="store", type=float,default=0.05,dest='alpha',nargs='?')
	parser.add_argument('--min_abd', action="store", type=float,default=0.0001,dest='min_abd',nargs='?')
	parser.add_argument('--min_samp', action="store", type=float,default=0.01,dest='min_samp',nargs='?')
	parser.add_argument('--tool_option1', action="store", dest='tool_option1',nargs='?')
	parser.add_argument('--root_dir', action="store", dest='root_dir',nargs='?')
	return  parser



######################################################################################
#  Build read config file                                                            #
######################################################################################
def build_read_config_file(strTempDir,results, DSrc, DMaaslin):
	fname = results.input
	input_file = open(fname)
	input_lines = input_file.readlines()
	LenInput = len(input_lines)
	input_file.close()
	TopLimit = int(results.lastmeta)
	ReadConfigFileName = os.path.join(strTempDir,"Test.read.config")
	Q = "'"

	WorkingDir = os.getcwd()
	os.chdir(DMaaslin)

	Limit1 = Q + "2-" + str(TopLimit )  + Q 
	ReadConfigTb1 = [
 	os.path.join(DSrc,"CreateReadConfigFile.R"),
	"-c",
	Limit1,
        ReadConfigFileName,
	"Metadata"
	">/dev/null",\
	"2>&1"
	]

	cmd_config1 = " ".join(ReadConfigTb1)

	os.system(cmd_config1)

	Limit2 = Q + str(TopLimit +1 )  + '-' + Q 
	ReadConfigTb2 = [
 	os.path.join(DSrc,"CreateReadConfigFile.R"),
	"-a",
	"-c",
	Limit2,
        ReadConfigFileName,
	"Abundance"
	">/dev/null",\
	"2>&1"
	]

	cmd_config2 = " ".join(ReadConfigTb2)
	os.system(cmd_config2)

	os.chdir(WorkingDir)
	return  ReadConfigFileName


######################################################################################
#   Main  Program                                                                    #
######################################################################################

# Parse commandline in
parser = read_params( sys.argv )
results = parser.parse_args()

### If option 2 is selected inform user on 2 outputs
if results.tool_option1 == "2":
	print "***Please note: 2 output files are  generated: Complete zipped results + Summary  ***"

### Project name
strProjectName = os.path.splitext(os.path.basename(results.input))[0]

### Define directory locations
D = os.path.join(results.root_dir,"tools","maaslin","maaslin")
DSrc = os.path.join(results.root_dir,"tools","maaslin","src")
DInput = os.path.join(results.root_dir,"tools","maaslin","input")
DMaaslin = os.path.join(results.root_dir,"tools","maaslin")
DMaaslinGalaxy = os.path.join(results.root_dir,"tools","maaslin","galaxy")


### Make temporary folder to work in
### Change permissions to make useable 
strTempDir = tempfile.mkdtemp()
cmd_chmod = "chmod 755 /" + strTempDir
os.system(cmd_chmod)
cmd_mkdir1 = "mkdir -m 755 " +  os.path.join(strTempDir,strProjectName)
os.system(cmd_mkdir1)

### Transpose the pcl file to a tsv file
TbCmdTranspose = [\
	"python",
	DMaaslinGalaxy  + "/transpose.py<" + str(results.input) +  ">" +  os.path.join(strTempDir,"output.tsv")\
	]
cmd_transpose = " ".join(TbCmdTranspose)
os.system(cmd_transpose)

### Make path for target output file
OutputFile = os.path.join(strTempDir,strProjectName,strProjectName+".txt")

### Make read config file
ReadConfigFileName = build_read_config_file(strTempDir,results, DSrc, DMaaslin)

### Build MaAsLin comamnd
CmdsArray = [\
os.path.join(DSrc,"Maaslin.R"),  \
"-d", str(results.alpha),\
"-r", str(results.min_abd),\
"-p", str(results.min_samp), \
"-i", \
ReadConfigFileName, \
OutputFile, \
os.path.join(strTempDir,"output.tsv"), \
"-v",\
"ERROR",\
">/dev/null",\
"2>&1"
]
invoke_maaslin_cmd = " ".join(CmdsArray)


### Write to directory cmd line used for troubleshooting
#CmdFileName = os.path.join(strTempDir,"cmdfile.txt")
#OutFile = open(CmdFileName,"w")
#OutputString = invoke_maaslin_cmd + "\n"
#OutFile.write(OutputString)
#OutFile.close()

### Call MaAsLin
os.system(invoke_maaslin_cmd)


### Copy output file to make available to galaxy
cmd_copy = "cp " + os.path.join(strTempDir,strProjectName+"/output.txt") + " " + results.output
MsgFileName = os.path.join(strTempDir,strProjectName+"/output.txt") 

if  not os.path.isfile(MsgFileName):
	cmd_copy = "cp " + os.path.join(strTempDir,strProjectName+"/output.txt") + " " + results.output
	OutFile = open(MsgFileName,"w")
	OutputString = "A MaAsLin error has occurred\n"
	OutputString = OutputString + "It typically happens when incorrect 'Last metadata row' was selected\n"
	OutputString = OutputString + "For demo data please choose 'Weight'\n"
	OutFile.write(OutputString)
	OutFile.close()

os.system(cmd_copy)

### Zip up output folder
cmd_zip = "zip -jr " + os.path.join(strTempDir,strProjectName+".zip") + " " + os.path.join(strTempDir,strProjectName) + ">/dev/null 2>&1"

os.system(cmd_zip)

### Copy output folder to make available to galaxy
cmd_copy_zip = "cp " + os.path.join(strTempDir,strProjectName+".zip") + " " + results.zip_file
os.system(cmd_copy_zip)

### Delete temp directory
cmd_del_tempdir = "rm -r " + strTempDir
######os.system(cmd_del_tempdir)
