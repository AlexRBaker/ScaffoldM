#!/usr/bin/env python
###############################################################################
#                                                                             #
#    dataloader.py                                                            #
#                                                                             #
#    A class for loading bamm files and contigs (FASTA)                       #
#                                                                             #
#    Copyright (C) Alexander Baker                                            #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Alexander Baker"
__copyright__ = "Copyright 2015"
__credits__ = ["Alexander Baker"]
__license__ = "GPLv3"
__maintainer__ = "Alexander Baker"
__email__ = "Alexander.baker@uqconnect.edu.au"

###############################################################################
###############################################################################
###############################################################################
###############################################################################
from bamm.bamParser import BamParser
from bamm.bamFile import BM_coverageType
from bamm.cWrapper import *
import numpy as np
import os
import sys
cov_type = BM_coverageType(CT.P_MEAN_OUTLIER, 1, 1)
BP = BamParser(cov_type)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Few goals in this file:
# Load bamms - use BamM to give cov and link information
# be able to extract links and coverage
# Check at least some contigs in contig file are present in bam
# if not abort - indicate they don't match
# 
###############################################################################
class DataLoader(object):
    """Utilities wrapper"""
    def __init__(self,
                 cov=False,
                 links=True,
                 bammloc,
                 contigloc):
	####
	# Variables:
	# coverage - binary - store coverage info
	# links - binary - store links info
	# bammloc - location of bamm file/s
	# contig loc - location of contig files
	####
	# Attributes
	# One for all var (same name)
	# Storage of Bamm processing
	# Extraction of useful Bamm bits
	# Some functions for getting info
	####
        if cov:
            self.cov=1
        else:
            self.cov=0
        if links:
            self.links=1
        else:
            self.links=0
        self.bammloc=bammloc
        self.contigloc=contigloc
        self.bammparse=BP.parseBams(bammloc,doLinks=self.links, /
        doCovs=self.cov,threads=min(len(bammloc),CompThreads)
        if self.cov:
            self.coverages=self.bammparse.BFI.coverages

        self.contigNames=self.bammparse.BFI.contigNames
        

#### Not storing contig or bamm in memory - accessing on call
##Might remove get functions if not used or decide to use numpy
##Will move to separate utility/general functions file later


###All for parsing tsv file of links (BamM command line output)
    def getcolumn(matrix,index):
		try:
			return [[row[i] for i in index] if len(index)>=2 \
			else row[index] for row in matrix]
		except TypeError:
			print "The index must be a index"
        
    def parsetsv(textfile="links.tsv"):
        ###Need to extend to attempt to detect delimiter
        ##Should be existing package for this
        import csv
        with open(textfile) as tsv:
            return [line for line in csv.reader(tsv,delimiter="\t")]
###Will parselinks file as intermediate until I fix BamM install
    if (False):
		self.linksfile=parsetsv()
            
    def findIDind(linkmatrix,ID='cid'):
        try:
            varis=getrow(linkmatrix,0)
            Colno=[i for (i, j) in enumerate(varis) if j.find(ID)>=0]
            return Colno[0]
        except SyntaxError or TypeError:
            print "The Id does not appear to be a string or \n the linkmatrix is not iterable"
        
    def getlinks(self,contig1,linksfile,contig2=False):
		morecnames=getcolumn(linksfile,1)
        if contig2==False
            #Indices needed to check second column
            LinkIndices=[i for (i, j) in enumerate(getcolumn(linksfile,0))\
             if j==contig1 or getcolumn(linksfile,1)[i]==contig1] 
        ## Cover both cases of link being at cid 1 or cid2
             
            Links=[row for (i,row) in enumerate(linksfile)\
             if j==contig1 or  in LinkIndices]
             
        else:
            contig2names=getcolumn(linksfile,1) ##Gets 2nd col of contigs in links file
            #Indices uneeded
            LinkIndices=[i for i, j in enumerate(getcolumn(linksfile,0))#ContigsNamesfromLinkfile)\
             if j==contig1 and getcolumn(linksfile,1)[i]==contig2]
             ## Force to and from contigs to match specified
             
            Links=[row for (i,row) in enumerate(linksfile)\
             if i in LinkIndices]
             
        return Links
	##Does work - extracts linking reads and passes them out
    ##If just contig1 then gives all reads linking contig1
    ##If both then return reads linking contig 1 and contig 2
    
    def getcov(self,contig1):
        return [row for row in self.coverages if row[0]==contig1]
	##Does work -extracts and returns contig coverage info
	##maybe allow for separating on mutiply mapped libraries via BamM 
	
		

###Now functions for parsing BamM Python module output  
###Expect to be ~50-100 lines
### First pythonize C link files from BamM
    if self.links:
		self.links=self.bammparse.pythonizeLinks\
		(self.bammparse.BFI, os.path.basename(bamloc.strip(os.sep))
#In bamM pythonize links calls makekey to create a key from linkpairs class	
# returns return "%d,%d" % (self.cid1, self.cid2) 
# unsure how cid remains unique
#need to read more on how Mike constructed his keys and values
#need more info on the cids for that.
	def getdictlinks(self, contig1, contig2=False):
		try:
			linkkeys=self.links.keys()
			if contig2=False:
				###Current understanding of key is a string of cid1,2 sep by comma
				###Therefore search cid for desired contig name, extract if present
				return [[key,self.links[key]] for key in linkkeys if contig1 in key]
			else:
				return [[key,self.links[key]] for key in linkkeys if contig1 in key and if contig2 in key]
		except AttributeError:
			print "Likely that links is not a _dict"
		
		
		
	def 
    def runCommand(self, cmd):
        """Run a command and take care of stdout

        expects 'cmd' to be a string like "foo -b ar"

        returns (stdout, stderr)
        """
        from multiprocessing import Pool
        from subprocess import Popen, PIPE

        p = Popen(cmd.split(' '), stdout=PIPE)
        return p.communicate()

    def parseFile(self, filename):
        """parse a file"""
        import sys
        try:
            with open(filename, "r") as fh:
                for line in fh:
                    print line
        except:
            print "Error opening file:", filename, sys.exc_info()[0]
            raise 
