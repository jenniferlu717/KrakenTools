#!/usr/bin/env python
################################################################
#combine_kreports.py takes multiple kraken-style reports and combines
#them into a single report file
#Copyright (C) 2019-2020 Jennifer Lu, jennifer.lu717@gmail.com
#
#This file is part of KrakenTools
#KrakenTools is free software; you can redistribute it and/or modify
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
#Updated: 05/16/2019
#
#This program reads in multiple Kraken report files and generates
#a combined Kraken report with columns for read counts and summarized
#read counts for each sample, along with two columns for across-sample sums
#
#Parameters:
#   -h, --help................show help message.
#   -r X, --report-file X.....all input kraken reports (separated by spaces)
#   -o X, --output X..........output kraken report filename
#   --display-headers.........includes header lines mapping samples to abbreviated names
#                             [default:true]
#   --no-headers..............do not include header lines [default:false] 
#   --sample-names............sample names for each kraken report (separated by spaces)
#                             [if none are given, each sample is given names S1, S2, etc] 
#Each Input report file format (tab-delimited)
#   - percentage of total reads
#   - number of reads (including reads within subtree)
#   - number of reads (only at this level)
#   - taxonomic classification level (U, D, P, C, O, F, G, S,...etc)
#   - NCBI taxonomic ID
#   - name of level
#Output file format (tab-delimited)
#   - percentage of total reads (for summed reads)
#   - combined number of reads (including reads within subtree)
#   - combined number of reads (only at this level)
#   - S1_all_reads, S1_lvl_reads, S2_all_reads, S2_lvl_reads, ...etc.
#   - taxonomic classification level (U, D, P, C, O, F, G, S,...etc)
#   - NCBI taxonomic ID
#   - name of level
#Methods 
#   - main
#   - process_kraken_report
####################################################################
import os, sys, argparse
import operator
from time import gmtime 
from time import strftime 

#Tree Class 
#usage: tree node used in constructing a taxonomy tree
#   including only the taxonomy levels and genomes identified in the Kraken report
class Tree(object):
    'Tree node.'
    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads, children=None, parent=None):
        self.name = name
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.tot_all = all_reads
        self.tot_lvl = lvl_reads
        self.all_reads = {}
        self.lvl_reads = {}
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self,node):
        assert isinstance(node,Tree)
        self.children.append(node)
    def add_reads(self, sample, all_reads, lvl_reads):
        self.all_reads[sample] = all_reads
        self.lvl_reads[sample] = lvl_reads
        self.tot_all += all_reads
        self.tot_lvl += lvl_reads
    def __lt__(self,other):
        return self.tot_all < other.tot_all
         
####################################################################
#process_kraken_report
#usage: parses a single line in the kraken report and extracts relevant information
#input: kraken report file with the following tab delimited lines 
#   - percent of total reads   
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level 
#       (U, - (root), - (cellular org), D, P, C, O, F, G, S) 
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
#   - spaces + name 
#returns:
#   - classification/genome name
#   - taxonomy ID for this classification
#   - level for this classification (number)
#   - level name (U, -, D, P, C, O, F, G, S)
#   - all reads classified at this level and below in the tree
#   - reads classified only at this level
def process_kraken_report(curr_str):
    split_str = curr_str.strip().split('\t')
    if len(split_str) < 5:
        return []
    try:
        int(split_str[1])
    except ValueError:
        return []
    #Extract relevant information
    all_reads =  int(split_str[1])
    level_reads = int(split_str[2])
    level_type = split_str[3]
    taxid = split_str[4] 
    #Get name and spaces
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1 
        else:
            break 
    #Determine which level based on number of spaces
    level_num = int(spaces/2)
    return [name, taxid, level_num, level_type, all_reads, level_reads]
    
####################################################################
#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--report-file','--report-files',
        '--report','--reports', required=True,dest='r_files',nargs='+',
        help='Input kraken report files to combine (separate by spaces)') 
    parser.add_argument('-o','--output', required=True,dest='output',
        help='Output kraken report file with combined information')
    parser.add_argument('--display-headers',required=False,dest='headers',
        action='store_true', default=True,
        help='Include header lines')
    parser.add_argument('--no-headers',required=False,dest='headers',
        action='store_false',default=True,
        help='Do not include header lines')
    parser.add_argument('--sample-names',required=False,nargs='+',
        dest='s_names',default=[],help='Sample names to use as headers in the new report')
    parser.add_argument('--only-combined', required=False, dest='c_only',
        action='store_true', default=False, 
        help='Include only the total combined reads column, not the individual sample cols')
    args=parser.parse_args()
    

    #Initialize combined values 
    main_lvls = ['U','R','D','K','P','C','O','F','G','S']
    map_lvls = {'kingdom':'K', 'superkingdom':'D','phylum':'P','class':'C','order':'O','family':'F','genus':'G','species':'S'}
    count_samples = 0
    num_samples = len(args.r_files)
    sample_names = args.s_names
    root_node = -1 
    prev_node = -1
    curr_node = -1
    u_reads = {0:0} 
    total_reads = {0:0} 
    taxid2node = {}

    #Check input values 
    if len(sample_names) > 0 and len(sample_names) != num_samples: 
        sys.stderr.write("Number of sample names provided does not match number of reports\n")
        sys.exit(1)
    #Map names
    id2names = {} 
    id2files = {} 
    if len(sample_names) == 0:
        for i in range(num_samples):
            id2names[i+1] = "S" + str(i+1)
            id2files[i+1] = ""
    else:
        for i in range(num_samples):
            id2names[i+1] = sample_names[i] 
            id2files[i+1] = ""
    
    #################################################
    #STEP 1: READ IN REPORTS
    #Iterate through reports and make combined tree! 
    sys.stdout.write(">>STEP 1: READING REPORTS\n")
    sys.stdout.write("\t%i/%i samples processed" % (count_samples, num_samples))
    sys.stdout.flush()
    for r_file in args.r_files:
        count_samples += 1 
        sys.stdout.write("\r\t%i/%i samples processed" % (count_samples, num_samples))
        sys.stdout.flush()
        id2files[count_samples] = r_file
        #Open File 
        curr_file = open(r_file,'r')
        for line in curr_file: 
            report_vals = process_kraken_report(line)
            if len(report_vals) < 5:
                continue
            [name, taxid, level_num, level_id, all_reads, level_reads] = report_vals
            if level_id in map_lvls:
                level_id = map_lvls[level_id]
            #Total reads 
            total_reads[0] += level_reads
            total_reads[count_samples] = level_reads 
            #Unclassified 
            if level_id == 'U' or taxid == '0':
                u_reads[0] += level_reads
                u_reads[count_samples] = level_reads 
                continue
            #Tree Root 
            if taxid == '1': 
                if count_samples == 1:
                    root_node = Tree(name, taxid, level_num, 'R', 0,0)
                    taxid2node[taxid] = root_node 
                root_node.add_reads(count_samples, all_reads, level_reads) 
                prev_node = root_node
                continue 
            #Move to correct parent
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent
            #IF NODE EXISTS 
            if taxid in taxid2node: 
                taxid2node[taxid].add_reads(count_samples, all_reads, level_reads) 
                prev_node = taxid2node[taxid]
                continue 
            #OTHERWISE
            #Determine correct level ID
            if level_id == '-' or len(level_id)> 1:
                if prev_node.level_id in main_lvls:
                    level_id = prev_node.level_id + '1'
                else:
                    num = int(prev_node.level_id[-1]) + 1
                    level_id = prev_node.level_id[:-1] + str(num)
            #Add node to tree
            curr_node = Tree(name, taxid, level_num, level_id, 0, 0, None, prev_node)
            curr_node.add_reads(count_samples, all_reads, level_reads)
            taxid2node[taxid] = curr_node
            prev_node.add_child(curr_node)
            prev_node = curr_node 
        curr_file.close()

    sys.stdout.write("\r\t%i/%i samples processed\n" % (count_samples, num_samples))
    sys.stdout.flush()

    #################################################
    #STEP 2: SETUP OUTPUT FILE
    sys.stdout.write(">>STEP 2: WRITING NEW REPORT HEADERS\n")
    o_file = open(args.output,'w') 
    #Lines mapping sample ids to filenames
    if args.headers: 
        o_file.write("#Number of Samples: %i\n" % num_samples) 
        o_file.write("#Total Number of Reads: %i\n" % total_reads[0])
        for i in id2names:
            o_file.write("#")
            o_file.write("%s\t" % id2names[i])
            o_file.write("%s\n" % id2files[i])
        #Report columns
        o_file.write("#perc\ttot_all\ttot_lvl")
        if not args.c_only:
            for i in id2names:
                o_file.write("\t%s_all" % i)
                o_file.write("\t%s_lvl" % i)
        o_file.write("\tlvl_type\ttaxid\tname\n")
    #################################################
    #STEP 3: PRINT TREE
    sys.stdout.write(">>STEP 3: PRINTING REPORT\n")
    #Print line for unclassified reads
    o_file.write("%0.4f\t" % (float(u_reads[0])/float(total_reads[0])*100))
    for i in u_reads:
        if i == 0 or (i > 0 and not args.c_only):
            o_file.write("%i\t" % u_reads[i])
            o_file.write("%i\t" % u_reads[i])
    o_file.write("U\t0\tunclassified\n")
    #Print for all remaining reads 
    all_nodes = [root_node]
    curr_node = -1
    curr_lvl = 0
    prev_node = -1
    while len(all_nodes) > 0:
        #Remove node and insert children
        curr_node = all_nodes.pop()
        if len(curr_node.children) > 0:
            curr_node.children.sort()
            for node in curr_node.children:
                all_nodes.append(node)
        #Print information for this node 
        o_file.write("%0.4f\t" % (float(curr_node.tot_all)/float(total_reads[0])*100))
        o_file.write("%i\t" % curr_node.tot_all)
        o_file.write("%i\t" % curr_node.tot_lvl)
        if not args.c_only:
            for i in range(num_samples):
                if (i+1) not in curr_node.all_reads: 
                    o_file.write("0\t0\t")
                else:
                    o_file.write("%i\t" % curr_node.all_reads[i+1])
                    o_file.write("%i\t" % curr_node.lvl_reads[i+1])
        o_file.write("%s\t" % curr_node.level_id)
        o_file.write("%s\t" % curr_node.taxid)
        o_file.write(" "*curr_node.level_num*2)
        o_file.write("%s\n" % curr_node.name)
    o_file.close() 
####################################################################
if __name__ == "__main__":
    main()
