#!/usr/bin/env python

import sys
import re
import os, os.path
import normalizeUse
import replacer
import addSynopsis
from sys import argv

operatorsStr=r"\.(?:and|eqv?|false|g[et]|l[et]|n(?:e(?:|qv)|ot)|or|true)\."

keywordsStr="(?:a(?:llocat(?:able|e)|ssign(?:|ment))|b(?:ackspace|lock)|c(?:a(?:ll|se)|haracter|lose|o(?:m(?:mon|plex)|nt(?:ains|inue))|ycle)|d(?:ata|eallocate|imension|o(?:|uble))|e(?:lse(?:|if|where)|n(?:d(?:|do|file|if)|try)|quivalence|x(?:it|ternal))|f(?:or(?:all|mat)|unction)|goto|i(?:f|mplicit|n(?:clude|quire|t(?:e(?:ger|nt|rface)|rinsic)))|logical|module|n(?:amelist|one|ullify)|o(?:nly|p(?:en|erator|tional))|p(?:a(?:rameter|use)|ointer|r(?:ecision|i(?:nt|vate)|o(?:cedure|gram))|ublic)|re(?:a[dl]|cursive|sult|turn|wind)|s(?:ave|e(?:lect|quence)|top|ubroutine)|t(?:arget|hen|ype)|use|w(?:h(?:ere|ile)|rite))"

intrinsic_procStr=r"(?:a(?:bs|c(?:har|os)|djust[lr]|i(?:mag|nt)|ll(?:|ocated)|n(?:int|y)|s(?:in|sociated)|tan2?)|b(?:it_size|test)|c(?:eiling|har|mplx|o(?:njg|sh?|unt)|shift)|d(?:ate_and_time|ble|i(?:gits|m)|ot_product|prod)|e(?:oshift|psilon|xp(?:|onent))|f(?:loor|raction)|huge|i(?:a(?:char|nd)|b(?:clr|its|set)|char|eor|n(?:dex|t)|or|shftc?)|kind|l(?:bound|en(?:|_trim)|g[et]|l[et]|og(?:|10|ical))|m(?:a(?:tmul|x(?:|exponent|loc|val))|erge|in(?:|exponent|loc|val)|od(?:|ulo)|vbits)|n(?:earest|int|ot)|p(?:ack|r(?:e(?:cision|sent)|oduct))|r(?:a(?:dix|n(?:dom_(?:number|seed)|ge))|e(?:peat|shape)|rspacing)|s(?:ca(?:le|n)|e(?:lected_(?:int_kind|real_kind)|t_exponent)|hape|i(?:gn|nh?|ze)|p(?:acing|read)|qrt|um|ystem_clock)|t(?:anh?|iny|r(?:ans(?:fer|pose)|im))|u(?:bound|npack)|verify)(?= *\()"

toUpcaseRe=re.compile("(?<![A-Za-z0-9_%#])(?<!% )(?P<toUpcase>"+operatorsStr+
                      "|"+ keywordsStr +"|"+ intrinsic_procStr +
                      ")(?![A-Za-z0-9_%])",flags=re.IGNORECASE)
linePartsRe=re.compile("(?P<commands>[^\"'!]*)(?P<comment>!.*)?"+
                       "(?P<string>(?P<qchar>[\"']).*?(?P=qchar))?")

def upcaseStringKeywords(line):
    """Upcases the fortran keywords, operators and intrinsic routines
    in line"""
    res=""
    start=0
    while start<len(line):
        m=linePartsRe.match(line[start:])
        if not m: raise SyntaxError("Syntax error, open string")
        res=res+toUpcaseRe.sub(lambda match: match.group("toUpcase").upper(),
                               m.group("commands"))
        if m.group("comment"):
            res=res+m.group("comment")
        if m.group("string"):
            res=res+m.group("string")
        start=start+m.end()
    return res

def upcaseKeywords(infile,outfile,logFile=sys.stdout):
    """Writes infile to outfile with all the fortran keywords upcased"""
    lineNr=0
    try:
        while 1:
            line=infile.readline()
            lineNr=lineNr+1
            if not line: break
            outfile.write(upcaseStringKeywords(line))
    except SyntaxError, e:
        e.lineno=lineNr
        e.text=line
        raise

def prettifyFile(infile,outfile,normalize_use=1, upcase_keywords=1,
             interfaces_dir=None,replace=None,logFile=sys.stdout):
    """prettifyes the fortran source in infile into outfile
    if normalize_use normalizes the use statements (defaults to true)
    if upcase_keywords upcases the keywords (defaults to true)
    if interfaces_dir is defined (and contains the directory with the
    interfaces) updates the synopsis
    if replace does the replacements contained in replacer.py (defaults
    to false)"""
    ifile=infile
    try:
        if normalize_use:
            tmpfile=os.tmpfile()
            normalizeUse.rewriteUse(ifile,tmpfile,logFile)
            tmpfile.seek(0)
            ifile.close()
            ifile=tmpfile
        if replace:
            tmpfile=os.tmpfile()
            replacer.replaceWords(ifile,tmpfile,logFile)
            tmpfile.seek(0)
            ifile.close()
            ifile=tmpfile
        if upcase_keywords:
            tmpfile=os.tmpfile()
            upcaseKeywords(ifile,tmpfile,logFile)
            tmpfile.seek(0)
            ifile.close()
            ifile=tmpfile
        if interfaces_dir:
            fileName=os.path.basename(infile.name)
            fileName=fileName[:fileName.rfind(".")]
            try:
                interfaceFile=open(os.path.join(interfaces_dir,
                                                fileName+".int"),"r")
            except:
                logFile.write("error opening file "+
                              os.path.join(interfaces_dir,
                                           fileName+".int")+"\n")
                logFile.write("skipping addSynopsis step for "+fileName+"\n")
                interfaceFile=None
            if interfaceFile:
                tmpfile=os.tmpfile()
                addSynopsis.addSynopsisToFile(interfaceFile,ifile,
                                              tmpfile,logFile=logFile)
                tmpfile.seek(0)
                ifile.close()
                ifile=tmpfile
        while 1:
            line=ifile.readline()
            if not line: break
            outfile.write(line)
        ifile.close()
    except:
        logFile.write("error processing file '"+infile.name+"'\n")
        outfile.close()
        os.rename(outfile.name,outfile.name+".err")
        raise

if __name__ == '__main__':
    usageDesc=("usage:\n"+sys.argv[0]+ """
    [--[no-]upcase] [--[no-]normalize-use] [--[no-]replace]
    [--interface-dir=~/cp2k/obj/platform/target] [--help]
    out_dir file1 [file2 ...]

    Writes file1,... to outdir after performing on them upcase of the
    fortran keywords, and normalizion the use statements.
    If the interface direcory is given updates also the synopsis.
    If requested the replacements performed by the replacer.py script
    are also preformed""")
    argStart=1
    normalize_use=1
    upcase_keywords=1,
    interfaces_dir=None
    replace=None
    if "--help" in sys.argv:
        print usageDesc
        sys.exit(0)
    for i in range(1,len(sys.argv)):
        if (sys.argv[i][:2]!="--"): break
        argStart=i+1
        m=re.match(r"--(no-)?(normalize-use|upcase)",sys.argv[i])
        if m:
            if (m.groups()[1]=="upcase"):
                upcase_keywords=not m.groups()[0]
            else:
                normalize_use=not m.groups()[0]
        else:
            m=re.match(r"--interface-dir=(.*)",sys.argv[i])
            if m:
                interfaces_dir=os.path.expanduser(m.groups()[0])
            else:
                print "ignoring unknown directive",sys.argv[i]
    if len(sys.argv)-argStart<2:
        print usageDesc
    else:
        outDir=sys.argv[argStart]
        if not os.path.isdir(outDir):
            print "out_dir must be a directory"
            print usageDesc
        else:
            for fileName in sys.argv[argStart+1:]:
                infile=open(fileName,'r')
                outfile=open(os.path.join(outDir,
                                          os.path.basename(fileName)),'w')
                prettifyFile(infile,outfile,normalize_use=normalize_use,
                             upcase_keywords=upcase_keywords,
                             interfaces_dir=interfaces_dir,
                             replace=replace)
