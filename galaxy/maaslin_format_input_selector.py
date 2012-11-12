#!/usr/bin/env python

"""
Author: George Weingart
Description: Dynamically read columns from input file for UI
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

import sys,string,time
from pprint import pprint

def red(st,l):
	if len(st) <= l: return st 
	l1,l2 = l/2,l/2
	return st[:l1]+".."+st[len(st)-l2:]

def get_cols(data,full_names):
	if data == "": return []
	max_len =32 
        fname = data.dataset.file_name
	input_file = open(fname)
	input_lines = input_file.readlines()
	input_file.close()
	table_lines = []
	for x in input_lines:
		first_column = x.split('\t')[0]
		table_lines.append(first_column)

	opt = []
	rc = ''
	lines = []
        try:
		lines = [(red((rc+v.split()[0]),max_len),'%d' % (i+1),False) for i,v in enumerate(table_lines) if v]

	except:
		l1 = '*ALL*'
		l2 = 1
		l3 = False
		MyList = [l1,l2,l3]
		lines.append(MyList)
	return opt+lines

def get_cols_add_line(data,full_names,lastmeta):
	if data == "": return []
	display_to = 1
	try:
		display_to = int(lastmeta)
	except:		
		pass

	max_len = 32 
        fname = data.dataset.file_name
	input_file = open(fname)
	input_lines = input_file.readlines()
	input_file.close()
	table_lines = []
	for x in input_lines:
		first_column = x.split('\t')[0]
		table_lines.append(first_column)
	table_lines.insert(0,'-')
	if  not display_to == 1:
		del  table_lines[display_to + 1:]


	opt = []
	rc = ''
	lines = []
        try:
		lines = [(red((rc+v.split()[0]),max_len),'%d' % (i+1),False) for i,v in enumerate(table_lines) if v]

	except:
		l1 = '*ALL*'
		l2 = 1
		l3 = False
		MyList = [l1,l2,l3]
		lines.append(MyList)
	return opt+lines

def get_cols_features(data,full_names,lastmeta):
	if data == "": return []
	display_from = 1
	try:
		display_from = int(lastmeta)
	except:		
		pass
	max_len = 32 
        fname = data.dataset.file_name
	input_file = open(fname)
	input_lines = input_file.readlines()
	input_file.close()
	table_lines = []
	for x in input_lines:
		first_column = x.split('\t')[0]
		table_lines.append(first_column)
	
	opt = []
	rc = ''
	del table_lines[:display_from]
	lines = []
        try:
		lines = [(red((rc+v.split()[0]),max_len),'%d' % (i+1),False) for i,v in enumerate(table_lines) if v]

	except:
		l1 = '*ALL*'
		l2 = 1
		l3 = False
		MyList = [l1,l2,l3]
		lines.append(MyList)
	return opt+lines
