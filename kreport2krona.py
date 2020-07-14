#!/usr/bin/env python
####################################################################
#kreport2krona.py converts a Kraken-style report into Krona-compatible format
#Copyright (C) 2019-2020 Jennifer Lu, jennifer.lu717@gmail.com

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
#Updated: 05/30/2019
#
#This program reads in a Kraken report file and generates
#a Krona-compatible Text report that can be imported into Krona using ktImportText. 
#
#Parameters:
#   -h, --help...............show help message
#   -r X, --report-file X....input kraken report filename 
#   -o X, --output X.........output Krona text filename 
#   --intermediate-ranks.....include non-traditional taxonomic ranks 
#   --no-intermediate-ranks..do not include non-traditional ranks [default] 
#Input file format (tab-delimited)
#   - percentage of total reads
#   - number of reads (including reads within subtree)
#   - number of reads (only at this level)
#   - taxonomic classification of level (U, D, P, C, O, F, G, S,...etc) 
#   - NCBI taxonomic ID 
#   - name of level   
#Output file format (tab-delimited)
#   - Number of reads within subtree of the specified level
#   - Taxonomy tree levels tab-delimited, with level type [k,p,c,o,f,g,s,x]
#
#Methods
#   - main
#   - process_kraken_report
####################################################################
import os, sys, argparse

####################################################################
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
    all_reads = int(split_str[1])
    lvl_reads = int(split_str[2])
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
    return [name, level_num, level_type, lvl_reads]

###################################################################
#kreport2krona_all
#usage: prints all levels for a kraken report 
#input: kraken report file and output krona file names 
#returns: none 
def kreport2krona_all(report_file, out_file):
    #Process report file and output 
    curr_path = [] 
    prev_lvl_num = -1
    r_file = open(report_file, 'r')
    o_file = open(out_file, 'w')
    #Read through report file 
    main_lvls = ['D','P','C','O','F','G','S']
    for line in r_file:
        report_vals = process_kraken_report(line)
        #If header line, skip
        if len(report_vals) < 4: 
            continue
        #Get relevant information from the line 
        [name, level_num, level_type, lvl_reads] = report_vals
        if level_type == 'U':
            o_file.write(str(lvl_reads) + "\tUnclassified\n")
            continue
        #Create level name 
        if level_type not in main_lvls:
            level_type = "x"
        elif level_type == "D":
            level_type = "K"
        level_str = level_type.lower() + "__" + name
        #Determine full string to add
        if prev_lvl_num == -1:
            #First level
            prev_lvl_num = level_num
            curr_path.append(level_str)
            o_file.write(str(lvl_reads) + "\t" + level_str + "\n")
        else:
            o_file.write(str(lvl_reads))
            #Move back if needed
            while level_num != (prev_lvl_num + 1):
                prev_lvl_num -= 1
                curr_path.pop()
            #Print all ancestors of current level followed by |
            for string in curr_path:
                if string[0] != "r": 
                    o_file.write("\t" + string)
            #Print final level and then number of reads
            o_file.write("\t" + level_str + "\n")
            #Update
            curr_path.append(level_str)
            prev_lvl_num = level_num
    o_file.close()
    r_file.close()
    
###################################################################
#kreport2krona_main
#usage: prints only main taxonomy levels for a kraken report 
#input: kraken report file and output krona file names 
#returns: none 
def kreport2krona_main(report_file, out_file):
    #Process report file and output 
    main_lvls = ['D','P','C','O','F','G','S']
    curr_path = [] 
    prev_lvl_num = -1
    num2path = {} 
    path2reads = {} 
    line_num = -1
    #Read through report file 
    r_file = open(report_file, 'r')
    for line in r_file:
        line_num += 1
        #########################################
        report_vals = process_kraken_report(line)
        #If header line, skip
        if len(report_vals) < 4: 
            continue
        #Get relevant information from the line 
        [name, level_num, level_type, lvl_reads] = report_vals
        if level_type == 'U':
            num2path[line_num] = ["Unclassified"]
            path2reads["Unclassified"] = lvl_reads 
            continue
        #########################################
        #Create level name 
        if level_type not in main_lvls:
            level_type = "x"
        elif level_type == "D":
            level_type = "K"
        level_str = level_type.lower() + "__" + name
        #########################################
        #Determine full string to add
        if prev_lvl_num == -1:
            #First level
            prev_lvl_num = level_num
            curr_path.append(level_str)
            #Save
            if curr_path[-1][0] == "x":
                num2path[line_num] = ""
            else:
                path2reads[curr_path[-1]] = lvl_reads
                num2path[line_num] = []
                for i in curr_path:
                    num2path[line_num].append(i)
            continue
        else:
            #########################################
            #Move back if needed
            while level_num != (prev_lvl_num + 1):
                prev_lvl_num -= 1
                curr_path.pop()
            #Update the list 
            curr_path.append(level_str)
            prev_lvl_num = level_num
            #########################################
            #IF AT NON-TRADITIONAL LEVEL, ADD TO PARENT
            if level_type == "x":
                test_num = len(curr_path) - 1
                while(test_num >= 0):
                    if curr_path[test_num][0] != "x":
                        path2reads[curr_path[test_num]] += lvl_reads 
                        test_num = -1
                    test_num = test_num - 1 
                num2path[line_num] = ""
            #IF AT TRADITIONAL LEVEL, SAVE 
            if level_type != "x":
                path2reads[curr_path[-1]] = lvl_reads
                num2path[line_num] = []
                for i in curr_path:
                    num2path[line_num].append(i)
    r_file.close() 
    
    #WRITE OUTPUT FILE
    o_file = open(out_file, 'w')
    for i in range(0,line_num):
        #Get values
        curr_path = num2path[i] 
        if len(curr_path) > 0:
            curr_reads = path2reads[curr_path[-1]] 
            if curr_path[-1][0] != "x":
                o_file.write("%i" % curr_reads)
            for name in curr_path:
                if name[0] != "r" and name[0] != "x":
                    o_file.write("\t%s" % name)
            o_file.write("\n")
    o_file.close()

######################################################################
#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--report-file', '--report', required=True,
        dest='r_file', help='Input kraken report file for converting')
    parser.add_argument('-o', '--output', required=True,
        dest='o_file', help='Output krona-report file name')
    parser.add_argument('--intermediate-ranks', action='store_true',
        dest='x_include', default=False, required=False,
        help='Include non-traditional taxonomic ranks in output')
    parser.add_argument('--no-intermediate-ranks', action='store_false',
        dest='x_include', default=False, required=False,
        help='Do not include non-traditional taxonomic ranks in output [default: no intermediate ranks]')
    args=parser.parse_args()

    #Determine which krona report to make 
    if args.x_include:
        kreport2krona_all(args.r_file,args.o_file)
    else:
        kreport2krona_main(args.r_file,args.o_file) 

#################################################################
if __name__ == "__main__":
    main()
#########################END OF PROGRAM##########################
