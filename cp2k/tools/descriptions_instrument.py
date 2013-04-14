#!/usr/bin/python
# -*- coding: utf-8 -*-

import os, sys
import pyparsing as pp
import re
from os import path


def main():    
	
	#fn = "/home/schuetto/data/nightly/manual/cp2k/cp2k/src/input_cp2k_motion.F"
	#instrument_file(fn)
	#return
	for (path, dirs, files) in os.walk("/home/schuetto/data/nightly/manual/cp2k/cp2k/src"):
		for fn in files:
			if(fn.lower().rsplit(".",1)[-1] not in ('f90','f')):
				continue
			instrument_file(path+"/"+fn)

#===============================================================================
def instrument_file(filename):
	print("Scanning file: %s"%filename)
	routine_call_parser = create_routine_call_parser()
	routines_to_instrument = ("section_create", "keyword_create", "cp_print_key_section_create","add_format_keyword")
	
	content = open(filename).read()
	
	
	if(not any([r in content for r in routines_to_instrument])): return
	
		
	#print content
	call_num = 0 
	next_pos = 0
	
#	for i, c in enumerate(content):
#		print i, c
	
	pp_hits = set()
	
	while(True):
		pos = next_pos
		pp_matches = list(routine_call_parser.scanString(content[pos:], maxMatches=1))
		
		if(len(pp_matches)==0): break
		assert(len(pp_matches)==1)
		routine,match_begin,match_end = list(pp_matches)[0]
		
		next_pos = pos + match_end
				
		if(routine.name not in routines_to_instrument): continue
		
		#print("Instrumenting call to:%s"%routine.name)
		
		call_num +=1
		pp_hits.add(pos+match_begin)
		
		call_id = "CALLID file:%s num:%d"%(path.basename(filename), call_num)
		
		# read named args
		named_args = dict()
		for arg in reversed(routine.args):
			if(len(arg)==1): continue
			a, b = pos+arg.value.begin, pos+arg.value.end
			old_value = content[a:b]
			new_value = None
			if(arg.key == "description"):
				new_value = change_description(old_value, arg.value, call_id)
			elif(arg.key == "enum_desc"):
				new_value = change_enum_desc(old_value, arg.value, call_id)
			if(new_value):
				content = content[:a] + new_value + content[b:]
				next_pos +=  - len(old_value) + len(new_value)
			
		

	re_pattern = re.compile("\n\s+CALL\s+("+("|".join(routines_to_instrument))+")")
	re_hits = set([m.start() for m in re_pattern.finditer(content, re.IGNORECASE)])
	
	
	#print sorted(re_hits)
	#print sorted(pp_hits)
	#print len(re_hits-pp_hits)
	#sys.exit()
	
	if(len(re_hits-pp_hits) > 0):
		print("Misses:")
		for i in re_hits-pp_hits:
			print("="*20)
			print(content[i:i+500])
		sys.exit(1)
	
	
	#print content
	#sys.exit()
	f = open(filename, "w")
	f.write(content)
	f.close()

#===============================================================================
def change_description(old_desc,tokens,call_id):
	if(old_desc[0] not in ('"', "'")):
		print("Found odd description: ")
		print(old_desc)
		return None
				
	new_desc = '"<<<%s field:description>>>"&\n'%call_id
	return(new_desc)


#===============================================================================
def change_enum_desc(old_code,tokens,call_id):
	if(tokens.inner[0] != "s2a"):
		print("Found odd enum_desc")
		print old_code
		return None
		
	n_items = len(tokens.inner.args)
	new_descs = ['"<<<%s field:enum_%d>>>"'%(call_id,i) for i in range(n_items)]
	new_code = "s2a(" + (',&\n'.join(new_descs)) + ")"
	return(new_code)
	
#===============================================================================
def create_routine_call_parser():
	name = pp.Word(pp.alphanums+"_.").setResultsName("name")
	constant = pp.Word(pp.alphanums+"-_.") | pp.quotedString
	marker = pp.Empty().setParseAction(lambda locn,tokens: locn)
	
	line_break = pp.Suppress(pp.Optional("&" + (pp.lineEnd | ("!"+pp.SkipTo(pp.lineEnd))) + pp.ZeroOrMore("!"+pp.SkipTo(pp.lineEnd)) + pp.Optional("&") ))
	def CommaSeparated(e): return( pp.Optional(e + pp.ZeroOrMore(pp.Suppress(",")+ e)) )
	def markered(e): return(pp.Group(marker.setResultsName("begin") + pp.Group(e).setResultsName("inner") + marker.setResultsName("end"))) 
	
	expr = pp.Forward()
	binary_op = pp.Literal("//") | pp.Literal("+") | pp.Literal("-") | pp.Literal("*") | pp.Literal("/")
	liste = "(/" + CommaSeparated(expr) + "/)"
	arg = pp.Group(pp.Optional(line_break + name.setResultsName("key") + pp.Suppress("=")) + markered(expr).setResultsName("value"))
	arg_list = ("(" + pp.Group(CommaSeparated(arg)).setResultsName("args") + ")")
	function_call = name + line_break + arg_list
	
	expr << line_break + (function_call ^ liste ^ name ^ constant ) + line_break + pp.Optional( binary_op + expr) 
	
	routine_call = "\n" + pp.Keyword("CALL") + name + arg_list
	routine_call.setWhitespaceChars(" \t")
	routine_call.parseWithTabs()
	#routine_call.enablePackrat()
	
	return(routine_call)
			
#===============================================================================		
if(__name__=="__main__"): main()
#EOF
