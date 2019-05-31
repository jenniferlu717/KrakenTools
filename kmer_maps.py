#! /usr/bin/env python
################################################################
#kmer_maps.py takes a Bracken databaseXmers.kraken file
#and provides the overall classification of Xmers among
#the main taxonomic levels.
#Copyright (C) 2019 Jennifer Lu, jlu26@jhmi.edu
#
#This file is part of Kraken-Tools
#Kraken-Tools is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or 
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.
#################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 05/29/2019
#
#This program takes a Bracken databaseXmers.kraken file
#and provides the overall classification of Xmers among
#the main taxonomic levels.
#
#Parameters
#   -h, --help................show help message.
#   -i, --input...............Input databaseXmers.kraken file 
#   -n, --nodes...............nodes.dmp file with taxonomy IDs --> taxonomy levels
#
#Input Kraken file format (tab-delimited)
#   - sequence ID 
#   - taxid
#   - taxid
#   - number of X length reads
#   - mapping of taxid:reads  
#Input nodes.dmp file format (tab|tab-delimited)
#   - taxid
#   - parent taxid
#   - level type 
#
#No output file: Output will be printed to stdout
# 
#################################################################
import os, argparse, sys
#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', required=True, dest='i_file',
        help='Input Kraken file produced by Bracken - databaseXmers.kraken')
    parser.add_argument('-n','--nodes', required=True, dest='n_file',
        help='Taxonomy nodes.dmp file')
    args=parser.parse_args()
    
    #Get taxid2lvl map
    taxid2lvl = {}  
    lvl2perc = {} 
    n_file = open(args.n_file,'r')
    for line in n_file:
        l_vals = line.strip().split("\t|\t")
        taxid = l_vals[0]
        rank = l_vals[2] 
        if rank not in lvl2perc:
            lvl2perc[rank] = []
        taxid2lvl[taxid] = rank   
    n_file.close()    
    lvl2perc['u'] = []

    #Make counts
    i_file = open(args.i_file,'r')
    for line in i_file:
        l_vals = line.strip().split("\t") 
        kmer_map = l_vals[4] 
        total = 0
        #Read in the mapping of reads to taxids and change to reads to lvls 
        curr_lvl2counts = {} 
        for val in kmer_map.split(" "):
            [taxid, count] = val.split(":")
            count = float(count)
            if taxid == '0':
                lvl = "u"
            else:
                lvl = taxid2lvl[taxid]  
            total += count 
            curr_lvl2counts[lvl] = count 
        #Save percentage of total reads mapping to each lvl 
        for lvl in lvl2perc:
            if lvl in curr_lvl2counts:
                perc = curr_lvl2counts[lvl]/float(total)
            else:
                perc = 0.0
            lvl2perc[lvl].append(perc)
    i_file.close() 

    #Print averages
    for lvl in lvl2perc:
        avg = sum(lvl2perc[lvl])/float(len(lvl2perc[lvl]))
        sys.stdout.write("%s\t%0.6f\n" % (lvl,avg)) 
#################################################################
if __name__ == "__main__":
    main()
########################END OF PROGRAM###########################

