#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
import pyparsing as pp
import re
from os import path

def main():    
	
	fn = "/home/schuetto/data/nightly/manual/cp2k/cp2k/src/input_cp2k_motion.F"
	instrument_file(fn)
##	for (path, dirs, files) in os.walk("./cp2k/cp2k/src/"):
#	for (path, dirs, files) in os.walk("./"):
#		for fn in files:
#			#if("." not in fn): continue
#			if(fn.lower().rsplit(".",1)[-1] not in ('f90','f')):
#				continue
#			instrument_file(path+fn)

#===============================================================================
def instrument_file(filename):
	print("Instrumenting: %s"%filename)
	routine_call_parser = create_routine_call_parser()
	routines_to_instrument = ("section_create", "keyword_create", "cp_print_key_section_create","add_format_keyword")
	                
	content = open(filename).read()
	#print content
	routine_num = 0 
	next_pos = 0
	
#	for i, c in enumerate(content):
#		print i, c
	
	while(True):
		pos = next_pos
		matches = list(routine_call_parser.scanString(content[pos:], maxMatches=1))
		if(len(matches)==0): break
		assert(len(matches)==1)
		routine,begin,end = list(matches)[0]
		routine_num +=1
		next_pos = pos + end
	
		#print routine.name
		
		if(routine.name not in routines_to_instrument): continue
		
		print("Instrumenting:"+routine.name)
		
		# read named args
		named_args = dict()
		for a in routine.args:
			if(len(a)==1): continue
			named_args[a.key] = a.value

		
		keyword_descr = named_args['description']
		
		if(keyword_descr.inner[0] != '"'):
			print "Found odd description:"+keyword_descr.inner
			continue
		
		a, b = pos+keyword_descr.begin, pos+keyword_descr.end
		old_descr = content[a: b]
		new_descr = '"FOOBAR file:%s routine:%d"'%(path.basename(filename),routine_num)
		content = content[:a] + new_descr + content[b:]
		next_pos = b  - len(old_descr) + len(new_descr)
		
	#print content
	f = open(filename, "w")
	f.write(content)
	f.close()


#===============================================================================
def create_routine_call_parser():
	name = pp.Word(pp.alphanums+"_.").setResultsName("name")
	constant = pp.Word(pp.alphanums+"_.") | pp.quotedString
	marker = pp.Empty().setParseAction(lambda locn,tokens: locn)
	
	line_break = pp.Suppress(pp.Optional(pp.Combine(pp.Optional("&") +  pp.lineEnd)))
	def CommaSeparated(e): return( pp.Optional(e + pp.ZeroOrMore(pp.Suppress(",")+ e)) )
	def markered(e): return(pp.Group(marker.setResultsName("begin") + pp.Combine(e).setResultsName("inner") + marker.setResultsName("end"))) 
	
	expr = pp.Forward()
	binary_op = pp.Literal(r"//")
	liste = "(/" + CommaSeparated(expr) + "/)"
	arg = pp.Group(pp.Optional(line_break + name.setResultsName("key") + pp.Suppress("=")) + markered(expr).setResultsName("value"))
	arg_list = ("(" + pp.Group(CommaSeparated(arg)).setResultsName("args") + ")")
	function_call = name + line_break + arg_list
	
	expr << line_break + (function_call ^ liste ^ name ^ constant ) + line_break + pp.Optional( binary_op + expr) 
	
	routine_call = pp.Keyword("CALL") + name + arg_list
	routine_call.setWhitespaceChars(" \t")
	routine_call.parseWithTabs()
	#routine_call.enablePackrat()
	
	return(routine_call)
			
#===============================================================================		
if(__name__=="__main__"): main()
#EOF
