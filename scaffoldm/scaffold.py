#!/usr/bin/env python
###############################################################################
#                                                                             #
#    scaffold  .py                                                            #
#                                                                             #
#    Class for storing and printing scaffolds and associated summary info     #
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
		##The key - scaffold name
		##value - dictionary with contig names as key
		##Inner value -start,end, orientation, position in scaffold
		##Cover duplicates with forced name addition
		self.scaffold={}
		self.name=bamname+"scaffold"
		self.contigNames=[]
		self.linesize=70
		
		for key in self.scaffold[scaffoldname]:
			self.contigNames.append(key)
			
	###Fancy algorithm will be replaced with the function which determines
	###assignment of contigs to scaffolds
	###Picks start,end and orientation
	fancyalgorithm(self)
	###Gapestimator - implementation of sahlini et al 2013 gap estimation
	###Appends gap size to all self.Scaffold list entries for
	### scaffold tig pairs
	gapestimator(self)
	
	def extractcontigs(contigname,contigloc):
		'''Just assigns contigs file via contigloc.
		Temporary just for use when making scaffold
		Will extract the text for that contig'''
		import sys
		try:
			with open(contigloc,'r') as Contigs:
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
	def arrangetigs(self, contigname,contigspec):
		##Simple, takes slice and rearranges contig for scaffold making
		##Orientation - 0 normal, 1 flipped
		tigseq=extractcontigs(contigname,self.contigloc)
		start=contigspec[0]
		end=contigspec[1]
		orientation=contigspec[2]
		if start==end:
			return tigseq[::(1-2*orientation)]
		elif not orientation:
			return tigseq[start:end:1]
		elif orientation:
			if start==0:
				return tigseq[end::-1]
			else:
				return tigseq[end:start-1:-1]
			
	def printscaffold(self):
			'''contigspec gonna probably be 5 entries - start,end,orientaiton.
		gap,order
		'''
		scaffoldname=self.name
		##First check if file exists
		filename=scaffoldname+".fasta"
		gapsize=3
		linesize=70
		OrderInScaffold=4
		if not os.path.isfile(filename):
			scaffile=open(filename,'w')
			scaffile.close()
		##Sorts the items in the scaffold by the order in scaffold parameter
		##Then, loops over this sorted list of tuples and extracts key
		##Then, return list of contig name in sorted order in scaffold
		Contigs=[item[0] for item in sorted(self.scaffold[scaffoldname].items(),\
		key=lambda input:input[1][OrderInScaffold])]
		
		with open(filename,'a') as scaffile:
			for tig in Contigs:
				for segment in chunker(arrangetigs(self,tig,\
					self.scaffold[scaffoldname][tig])+\
					self.scaffold[scaffoldname][tig][gapsize]*"N",self.linesize,"\n"):
					scaffile.write(segment)
					#scaffile.write(arrangetigs(self,tig,self.scaffold[scaffoldname][tig]))
					#scaffile.write(self.scaffold[scaffoldname][tig][gapsize]*"N")
	
def chunker(string,chunksize,end):
	''' Creates chunks from string and appends an end term
	Might need to make a cleaner function to remove all end terms in string'''
	'''Note this return generator - should be iterated over'''
	'''I should implement this in my test for increased speed'''
	try:
		stringmod=string.translate(None,end)
		for i in xrange(0,len(stringmod),chunksize):
			if len(stringmod)>=i+chunksize:
				yield stringmod[i:i+chunksize]+end
			else:
				yield stringmod[i:i+chunksize]
	except TypeError:
		print "end must be concatenable to string, \
		intended use is type(str) for both"
'''
test from string import maketrans

translate goodness for cleaning string
def test1(string,rep,char):
	return (string*rep).translate(None,char)
	
	
	
map = maketrans('aeiou', '0' * 5)
def str_translate(s, map):
    return s.translate(map)
    as far faster method of replacing characters''' 
'''
Just line used to experiment with yield and chunker.
for i in chunker("".join([tig[random.randrange(0,4,1)] for i in xrange(20000)]),70,""):
	print i
	'''
