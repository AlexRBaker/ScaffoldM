#!/usr/bin/env python
###############################################################################
#                                                                             #
#    scaffold.py                                                              #
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
import sys
import os

class Scaffold(object):
    """Utilities wrapper"""
    def __init__(self,
                scaffold,
                scaffoldname,
                linesize=70,
                contigloc='Dupes.fna'):        
        ##The scaffold will have the following structure
        ##The key - scaffold name
        ##value - dictionary with contig names as key
        ##Inner value -start,end, orientation, position in scaffold
        ##Cover duplicates with forced name addition
        self.scaffold=scaffold
        self.name=scaffoldname
        OrderInScaffold=4 #Ind contig place in scaffold
        if scaffold!={}:
            self.contigNames=[item[0] for item in sorted(self.scaffold[scaffoldname].items(),\
            key=lambda input:input[1][OrderInScaffold])]
        else:
            self.contigNames=None
        self.linesize=linesize
        self.contigloc=contigloc
            
    ###Fancy algorithm will be replaced with the function which determines
    ###assignment of contigs to scaffolds
    ###Picks start,end and orientation
    ###Gapestimator - implementation of sahlini et al 2013 gap estimation
    ###Appends gap size to all self.Scaffold list entries for
    ### scaffold tig pairs
    
    def extractcontigs(self,contigname,contigloc,header=True):
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
                        if header:
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
            print "Error opening file:", contigloc,sys.exc_info()[0]
            raise
    #all_keys = set().union(*(d.keys() for d in mylist)) gets all
    # from list of dicts
    def arrangetigs(self, contigname,contigspec,header=False):
        ##Simple, takes slice and rearranges contig for scaffold making
        ##Orientation - 0 normal, 1 flipped
        tigseq=self.extractcontigs(contigname,self.contigloc,header)
        start=contigspec[0]
        end=contigspec[1]
        orientation=contigspec[2]
        if start==end: #Screens for default case - no cuts to contig
            return tigseq[::(1-2*orientation)]
        elif not orientation: #Orientation determinines whether or not to flip
            return tigseq[start:end:1]
        elif orientation:
            if start==0:
                return tigseq[end::-1]
            else:
                return tigseq[end:start:-1]
            
    def printscaffold(self,filename=None,header=False):
        '''contigspec gonna probably be 5 entries - start,end,orientaiton.
        gap,order
        '''
        ##Can then optionally choose 
        #to write all to one file or one for each scaffodl
        scaffoldname=self.name
        if filename==None:
            filename=scaffoldname+".fasta"
        else:
            pass
        ##First check if file exists
        contigloc=self.contigloc
        gapsize=3 #Index for gap between contigs
        linesize=70 #Length of printed line
        if not os.path.isfile(filename):
            scaffile=open(filename,'w')
            scaffile.close()
        ##Sorts the items in the scaffold by the order in scaffold parameter
        ##Then, loops over this sorted list of tuples and extracts key
        ##Then, return list of contig name in sorted order in scaffold
        Contigs=self.contigNames
        with open(filename,'a+') as scaffile:
            ##Need to add check to ensure header is written on its own line
            scafheader=">{0}|Contigs:{1}\n".format(scaffoldname,",".join(Contigs))
            scaffile.write(scafheader)
            for tig in Contigs:
                scaffile.seek(0,0)
                offset=linesize
                ContigSeq=(self.arrangetigs(tig,\
                        self.scaffold[scaffoldname][tig])+\
                        self.scaffold[scaffoldname][tig][gapsize]*"N")
                #Only seems to work if more than one line with \n
                #Therefore, check if more than one line
                scaffile.seek(0,0)
                scaffile.readline() #Moving past first line
                line2=scaffile.readline()
                scaffile.seek(0,0)
                if line2=='':
                    modseq=ContigSeq[:linesize+2].translate(None,'\n')
                    if offset==linesize:
                        pass
                    else:
                        scaffile.write(modseq[:(linesize-offset)]+"\n")
                    for segment in self.chunker(ContigSeq[(linesize-offset):],linesize,"\n"):
                        scaffile.write(segment)
                else:
                    #Now - extra line break in newfile after header
                    scaffile.seek(0,2)
                    while scaffile.read(1)!="\n":
                        scaffile.seek(-2,1)
                    fileend=scaffile.readline()
                    offset=len(fileend.strip('\n'))
                    modseq=ContigSeq[:linesize+2].translate(None,'\n')
                    scaffile.seek(0,2)
                    if fileend==scafheader:
                        offset=0
                    else:
                        pass
                    if offset==linesize:
                        scaffile.write(modseq[:(linesize-offset)])
                    else:
                        scaffile.write(modseq[:(linesize-offset)]+"\n")
                    for segment in self.chunker(ContigSeq[(linesize-offset):],linesize,"\n"):
                        scaffile.write(segment)
            scaffile.seek(0,2)
            while scaffile.read(1)!="\n":
                scaffile.seek(-2,1)
            if not scaffile.readline().endswith('\n'):
                scaffile.write('\n')

    def chunker(self,string,chunksize,end):
        ''' Creates chunks from string and appends an end term
        Note this return generator - should be iterated over'''
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
            
    def loadtigsequence(self, contig1):
        return
    def loadscafsequence(self):
        return 
