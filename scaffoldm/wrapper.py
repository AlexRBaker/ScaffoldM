#!/usr/bin/env python
###############################################################################
#                                                                             #
#    wrapper.py                                                               #
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
#This files is the overall wrapper for the project
#It will call all the classes, dataloader, dataparser and scaffold  
#(Will maybe make a contig class)
#Progess is expected as follows
#dataloader is called - passes data to data parsers
#data parser calls gapest and algo to create the scaffold dicitonary
#Scaffold class is called and created for each scaffold
#Handles commandline input and python package organisation
###############################################################################
import os
import argparse
from dataloader import DataLoader
from dataparser import DataParser
from scaffold import Scaffold
###############################################################################
if __name__ == "__main__": ###Check if arguments coming in from command line
    parser = argparse.ArgumentParser(description='CreateNewContigs.')
    parser.add_argument('-b','--bam', type=str, nargs='*', \
		help='The .bamfiles to be processed',default='Dupes.MG1655refS100E20000Complete-Empirical.bam')
    parser.add_argument('-f','--fasta', type=str, nargs='?', \
	default='Dupes.fna', \
		help='name of links file')
    parser.add_argument('-l','--links', type=str, nargs='?', \
	default='links.tsv', \
		help='name of links file')
    parser.add_argument('-c','--covs', type=str, nargs='?', \
	default='covs.tsv', \
		help='name of coverage file')
    parser.add_argument('-i','--inserts', type=str, nargs='?', \
	default='inserts.tsv', \
		help='name of insert file')
    parser.add_argument('-n','--librarynumber', type=int, nargs='*', \
	default=1, \
		help='name of insert file')


    args = parser.parse_args()
    if type(args.bam)==str:
        bamnames=[args.bam]
    else:
        bamnames=args.bam
    contigloc=args.fasta
    linksnames=args.links
    covnames=args.covs
    insertnames=args.inserts
    libnos=args.librarynumber
    libnos=[str(ele) for ele in libnos]
    print libnos
    print args
    #print "Looking for a bug", os.getcwd()
    ##Loads up the bams and passes them onto BamM to get links, coverage etc.
    data=DataLoader(bamnames,
                 contigloc,
                 libno=libnos,
                 linksname=linksnames,
                 covname=covnames,
                 insertname=insertnames)
    print data.bamnames, "THese are the bamnames"
    ###Note dataparser will have inside it all the scaffolds made from the data.
    parser=DataParser(data.links,
                            data.coverages,
                            data.inserts,
                            data.bamnames,
                            data.contigNames,
                            data.contigloc,
                            data.tiglen)
    #Tells the parser to process the data
    parser.parse()
    #parser.output() #Prints some summary information used in quality assessments.
    #namely gap sizes and a scaffold summary file for easy parsing.
    
    ###Visualise will include some coverage based graphs,
    ###A representation of the network of linked contigs- shall be done via cytoscape
    ###And, maybe some graphs of gap size(These are Done)
    ###Also, have an N(x) graph - can crudely represent distribution of scaffodl lengths
    #parser.visualise()
