#! /usr/bin/env python
####################################################################
#kreport2mpa.py converts a Kraken-style report into mpa [MetaPhlAn) format
#Copyright (C) 2017-2020 Jennifer Lu, jennifer.lu717@gmail.com

#This file is part of KrakenTools.
#KrakenTools is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.

####################################################################
#Jennifer Lu, jlu26@jhmi.edu
#11/06/2017
#Updated: 07/12/2020
#
#This program reads in a Kraken report file and generates
#an mpa-format (MetaPhlAn) style report. Each line represents
#a possible taxon classification. The first column is lists the 
#domain, kingdom, phyla, etc, leading up to each taxon.
#The levels are separated by the | delimiter, with the type of 
#level specified before each name with a single letter and underscore
#(d_ for domain, k_ for kingdom, etc). 
#The second column is the number of reads classified within 
#that taxon's subtree.
#
#Input file:
#   - Kraken report file generates from the kraken raw output file
#Input Parameters to Specify [OPTIONAL]:
#   - header_line = prints a header line in mpa-report 
#       [Default: no header]
#   - intermediate-ranks = includes non-traditional taxon levels
#       (traditional levels: domain, kingdom, phylum, class, order, 
#       family, genus, species)
#       [Default: no intermediate ranks]
#Output file format (tab-delimited)
#   - Taxonomy tree levels |-delimited, with level type [d,k,p,c,o,f,g,s,x]
#   - Number of reads within subtree of the specified level
#
#Methods
#   - main
#   - process_kraken_report
#
import os, sys, argparse

#process_kraken_report
#usage: parses a single line in the kraken report and extracts relevant information
#input: kraken report file with the following tab delimited lines
#   - percent of total reads
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level
#       (U, D, P, C, O, F, G, S, -)
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria,...etc)
#   - spaces + name
#returns:
#   - classification/genome name
#   - level name (U, -, D, P, C, O, F, G, S)
#   - reads classified at this level and below in the tree
def process_kraken_report(curr_str):
    split_str = curr_str.strip().split('\t')
    try:
        int(split_str[1])
    except ValueError:
        return []
    percents = float(split_str[0])
    all_reads = int(split_str[1])
    level_type = split_str[3]
    #Get name and spaces 
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1
        else:
            break
    name = name.replace(' ','_')
    #Determine level based on number of spaces
    level_num = spaces/2
    return [name, level_num, level_type, all_reads, percents]

#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--report-file', '--report', required=True,
        dest='r_file', help='Input kraken report file for converting')
    parser.add_argument('-o', '--output', required=True,
        dest='o_file', help='Output mpa-report file name')
    parser.add_argument('--display-header', action='store_true', 
        dest='add_header', default=False, required=False,
        help='Include header [Kraken report filename] in mpa-report file [default: no header]') 
    parser.add_argument('--read_count', action='store_true',
        dest='use_reads', default=True, required=False,
        help='Use read count for output [default]')
    parser.add_argument('--percentages', action='store_false',
        dest='use_reads', default=True, required=False,
        help='Use percentages for output [instead of reads]')
    parser.add_argument('--intermediate-ranks', action='store_true',
        dest='x_include', default=False, required=False,
        help='Include non-traditional taxonomic ranks in output')
    parser.add_argument('--no-intermediate-ranks', action='store_false',
        dest='x_include', default=False, required=False,
        help='Do not include non-traditional taxonomic ranks in output [default]')
    args=parser.parse_args()

    #Process report file and output 
    curr_path = [] 
    prev_lvl_num = -1
    r_file = open(args.r_file, 'r')
    o_file = open(args.o_file, 'w')
    #Print header
    if args.add_header:
        o_file.write("#Classification\t" + os.path.basename(args.r_file) + "\n")
    
    #Read through report file 
    main_lvls = ['R','K','D','P','C','O','F','G','S']
    for line in r_file:
        report_vals = process_kraken_report(line)
        #If header line, skip
        if len(report_vals) < 5: 
            continue
        #Get relevant information from the line 
        [name, level_num, level_type, all_reads, percents] = report_vals
        if level_type == 'U':
            continue
        #Create level name 
        if level_type not in main_lvls:
            level_type = "x"
        elif level_type == "K":
            level_type = "k"
        elif level_type == "D":
            level_type = "k"
        level_str = level_type.lower() + "__" + name
        #Determine full string to add
        if prev_lvl_num == -1:
            #First level
            prev_lvl_num = level_num
            curr_path.append(level_str)
        else:
            #Move back if needed
            while level_num != (prev_lvl_num + 1):
                prev_lvl_num -= 1
                curr_path.pop()
            #Print if at non-traditional level and that is requested
            if (level_type == "x" and args.x_include) or level_type != "x":
                #Print all ancestors of current level followed by |
                for string in curr_path:
                    if (string[0] == "x" and args.x_include) or string[0] != "x":
                        if string[0] != "r": 
                            o_file.write(string + "|")
                #Print final level and then number of reads
                if args.use_reads:
                    o_file.write(level_str + "\t" + str(all_reads) + "\n")
                else:
                    o_file.write(level_str + "\t" + str(percents) + "\n")
            #Update
            curr_path.append(level_str)
            prev_lvl_num = level_num
    o_file.close()
    r_file.close()

if __name__ == "__main__":
    main()
