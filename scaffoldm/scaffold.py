#!/usr/bin/env python
###############################################################################
#                                                                             #
#    sacffold  .py                                                            #
#                                                                             #
#    Class for storing and printing scaffolds and assoicated summary info     #
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
#Overall scaffold structre
#Include how to print scaffold
#Will have dict detailing necessary info on scaffolds
#involved contig and scaffold names will be extracted
# Key - Tuple - contig and scaffold
# item - start point and end point in contig, orientation, downstream contig
###############################################################################

import dataloader
import dataparser
import fancyalgorithm
import gapestimator
import sys
import os

class Scaffold(dataloader):
    """Utilities wrapper"""
    def __init__(self,
				 bammloc,
                 contigloc,
                 links=True,
                 cov=False):

		dataloader.__init__(self,bammloc,contigloc,links,cov)
		
		##The scaffold will have the following structure
		##The key - (scaffold name, contig name)
		##Value - start point, end point, orientation
		##Cover duplicates with forced name addition
		self.Scaffold={}
		
		self.contigNames=[]
		
		for key in self.scaffold:
			self.contigNames.append(key[0])
			
	###Fancy algorithm will be replaced with the function which determines
	###assignment of contigs to scaffolds
	###Picks start,end and orientation
	fancyalgorithm(self)
	###Gapestimator - implementation of sahlini et al 2013 gap estimation
	###Appends gap size to all self.Scaffold list entries for
	### scaffold tig pairs
	gapestimator(self)
	
	
	def arrangetigs(self, contigname,contigspec):
		##Simple, takes slice and rearranges contig for scaffodl making
		##Orientation - 0 normal, 1 flipped
		tigseq=extractcontigs(self,contigname,self.contigloc)
		start=contigspec[0]
		end=contigspec[1]
		orientation=contigspec[2]
		if not orientation:
			return tigseq[start:end:1]
		elif orientation:
			if start==1:
				return tigseq[end::-1]
			else:
				return tigseq[end:start-1:-1]
		
	def makescaffold():
		
		Scaffolds={}
		ContigsSets=[]
		for key in self.scaffolds:
			if key not in Scaffolds.iterkeys():
				Scaffolds[key[0]]=[key[1]]
		
		
		
		Scaffolds=list(Scaffolds).sort()
		Contigs=self.contigNames.sort()
		
		for fold in Scaffolds:
			for tig in fold:
				self.scaffold[(

				
				
		
    def extractcontigs(self,contigname,contigloc):
    ###Just assigns contigs file via contigloc.
	###Temporary just for use when making scaffold
    ###Will extract the text for that contig
        import sys
        try:
            with open(contigfile,'r') as Contigs:
                head=Contigs.readline()
                if not head.startswith('>'):
                    raise TypeError("Not a FASTA file:")
                Contigs.seek(0)
                title=head[1:].rstrip() ##Strips whitespace and >
                record=0
                contigseq=[]
                for line in Contigs:
                    if line.startswith('>') and line.find(contigname)>=0:
                        record=1
                        contigseq.append(line)
                    elif line.startswith('>') and contigname not in line:
                        record=0
                    elif record==1:
                        contigseq.append(line)
                    else:
                        pass
                seq=''.join(contigseq)
                return seq
        except:
            print "Error opening file:", contigfile,sys.exc_info()[0]
            raise
	#all_keys = set().union(*(d.keys() for d in mylist)) gets all
	# from list of dicts

	def printscaffold(self):
		
		
		
		
		
		
		
		
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
