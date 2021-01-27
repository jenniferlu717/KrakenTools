#!/usr/bin/env python
######################################################################
#make_taxonomy.py takes in three files: nodes.dmp, names.dmp, and 
#seqid2taxid.map to create a condensed taxonomy file 
#Copyright (C) 2020 Jennifer Lu, jennifer.lu717@gmail.com
#
#This file is part of KrakenTools
#KrakenTools is free software; oyu can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.
#
######################################################################
#Jennifer Lu, jlu26@jhmi.edu
#Updated: 04/14/2020
#
#This program creates a condensed taxonomy file for a given Kraken database
#To create a taxonomy, the nodes.dmp, names.dmp, and seqid2taxid.map files
#must be provided
#
#Required Parameters:
#   --nodes X...........................nodes.dmp file
#   --names X...........................names.dmp file 
#   --seqid2taxid X.....................seqid2taxid.map file
#   -o, --output X......................output file with taxonomy info
#Optional Parameters:
#   -h, --help..........................show help message.
#################################################################################
import os, sys, argparse
from time import gmtime
from time import strftime 
#################################################################################
#Tree Class
#usage: tree node used in constructing taxonomy tree
class Tree(object):
    'Tree node.'
    def __init__(self,  taxid, level_rank, parent=None,children=None):
        self.taxid = taxid
        self.level_rank= level_rank
        #Other attributes for later
        self.name = ''
        self.level_num = -1
        self.p_taxid = -1
        #Parent/children attributes
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self, node):
        assert isinstance(node,Tree)
        self.children.append(node)
#################################################################################
#Main method
def main():
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--nodes',dest='nodes_file', required=True,
        help='nodes.dmp file from taxonomy')
    parser.add_argument('--names',dest='names_file', required=True,
        help='names.dmp file from taxonomy')
    parser.add_argument('--seqid2taxid',dest='s2t_file', required=True,
        help='seqid2taxid.map file')
    parser.add_argument('-o','--output',dest='out_file', required=True,
        help='output taxonomy file')
    args = parser.parse_args()

    #Start Program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')

    map_ranks = {'superkingdom':'D',
        'phylum':'P',
        'class':'C',
        'order':'O',
        'family':'F',
        'genus':'G',
        'species':'S'} 
    #STEP 1/5: PARSE NODES.DMP FILE
    root_node = -1
    taxid2node = {} 
    p_notsaved = {} 
    nodes_f = open(args.nodes_file,'r')
    sys.stdout.write(">> STEP 1/5: Reading %s\n" % args.nodes_file)
    sys.stdout.write("\t%0 nodes read")
    count_nodes = 0 
    for line in nodes_f:
        count_nodes += 1
        if (count_nodes % 100) == 0:
            sys.stdout.write("\r\t%i nodes read" % count_nodes)
            sys.stdout.flush()
        [curr_taxid,parent_taxid,rank] = line.strip().split("\t|\t")[0:3]
        #Make/Save node
        newrank = '-'
        if rank in map_ranks:
            newrank = map_ranks[rank]
        curr_node = Tree(curr_taxid, newrank)
        taxid2node[curr_taxid] = curr_node
        #root
        if curr_taxid == "1":
            curr_node.level_rank = 'R'
            root_node = curr_node
        elif parent_taxid in taxid2node:
            #save parent
            curr_node.parent = taxid2node[parent_taxid]
            curr_node.p_taxid = parent_taxid
            taxid2node[parent_taxid].add_child(curr_node)
        else:
            #parent not linked
            p_notsaved[curr_taxid] = curr_node
            curr_node.p_taxid = parent_taxid 
    nodes_f.close() 
    sys.stdout.write("\r\t%i nodes read\n" % count_nodes)
    sys.stdout.flush()
    #Fix parents
    for taxid in p_notsaved:
        p_taxid = p_notsaved[taxid].p_taxid
        if p_taxid not in taxid2node:
            sys.stderr.write("ERROR: %s not found in nodes.dmp file\n" % p_taxid)
            continue 
        p_node = taxid2node[p_taxid]
        p_notsaved[taxid].parent = p_node
        p_node.add_child(p_notsaved[taxid])
    #STEP 2/5: PARSE SEQID2TAXID FILE TO GET LEAF NODES 
    sys.stdout.write(">> STEP 2/5: Reading %s\n" % args.s2t_file)
    sys.stdout.write("\t%0 taxids read")
    count_taxids = 0
    s2t_f = open(args.s2t_file,'r')
    leaves = {}
    for line in s2t_f:
        taxid = line.strip().split('\t')[-1]
        #taxid not yet saved
        if taxid not in leaves:
            #node exists for taxid
            if taxid in taxid2node:
                leaves[taxid] = taxid2node[taxid]
                count_taxids += 1
                if (count_taxids % 100 == 0):
                    sys.stdout.write("\r\t%i taxids read" % count_taxids)
                    sys.stdout.flush()
            else:
                sys.stderr.write("ERROR: node for %s not found\n" % taxid)
    s2t_f.close()
    sys.stdout.write("\r\t%i taxids read\n" % count_taxids)
    sys.stdout.flush()
    #STEP 3/5: CONDENSE TAXONOMY TO ONLY INCLUDE TAXIDS IN SEQID2TAXID.MAP
    sys.stdout.write(">> STEP 3/5: Condensing taxonomy\n")
    sys.stdout.write("\t%saving 0 taxids")
    count_final = 0
    save_taxids = {} 
    for leaf in leaves:
        if leaf in save_taxids:
            continue
        save_taxids[leaf] = leaves[leaf]
        p_node = taxid2node[leaf].parent
        count_final += 1
        sys.stdout.write("\r\tsaving %i taxids" % count_final)
        sys.stdout.flush()
        while(p_node != None):
            if p_node.taxid not in save_taxids:
                #travel up path to root 
                save_taxids[p_node.taxid] = p_node
                p_node = p_node.parent
                count_final += 1
                sys.stdout.write("\r\tsaving %i taxids" % count_final)
                sys.stdout.flush()
            else: 
                #Parent path already parsed 
                p_node = None
    sys.stdout.write("\r\tsaving %i taxids\n" % count_final)
    sys.stdout.flush()
    #STEP 4/5: PARSE NAMES.DMP TO GET NAMES FOR TAXIDS IN TREE 
    sys.stdout.write(">> STEP 4/5: Reading %s\n" % args.names_file)
    count_names = 0
    sys.stdout.write("\t%i/%i names found" % (count_names, count_final))
    names_f = open(args.names_file,'r')
    for line in names_f:
        [taxid,name] = line.strip().split('\t|\t')[0:2]
        if taxid in save_taxids:
            if save_taxids[taxid].name == '':
                save_taxids[taxid].name = name
                count_names += 1
                sys.stdout.write("\r\t%i/%i names found" % (count_names, count_final))
                sys.stdout.flush()
            elif "scientific name" in line:
                save_taxids[taxid].name = name
    names_f.close()
    sys.stdout.write("\r\t%i/%i names found\n" % (count_names, count_final))
    sys.stdout.flush()
    #STEP 5/5: PRINT NEW TAXONOMY 
    sys.stdout.write(">> STEP 5/5: Printing final taxonomy to %s\n" % args.out_file)
    print_count = 0
    sys.stdout.write("\t%i nodes printed" % print_count)
    sys.stdout.flush()
    o_file = open(args.out_file,'w')
    parse_nodes = [root_node]
    root_node.level_num = 0
    printed_taxids = {"1":root_node}
    while len(parse_nodes) > 0:
        #Get the first node in list
        curr_node = parse_nodes.pop(0)
        print_count += 1
        sys.stdout.write("\r\t%i nodes printed" % print_count)
        sys.stdout.flush()
        #Print current node
        o_file.write("%s\t|\t" % curr_node.taxid)
        if curr_node.taxid == "1":
            o_file.write("1\t|\t")
        else:
            o_file.write("%s\t|\t" % curr_node.p_taxid)
        o_file.write("%s\t|\t" % curr_node.level_rank)
        o_file.write("%s\t|\t" % curr_node.level_num)
        o_file.write("%s\n" % curr_node.name)
        #Parse through children 
        for child in curr_node.children: 
            #Only save taxids needed
            if child.taxid in save_taxids:
                child.level_num = curr_node.level_num + 1
                #Fix children ranks
                if child.level_rank == '-': 
                    if len(curr_node.level_rank) == 1:
                        child.level_rank = curr_node.level_rank + "1"
                    else:
                        new_num = int(curr_node.level_rank[1:]) + 1
                        child.level_rank = curr_node.level_rank[0] + str(new_num)    
                parse_nodes.append(child)
                printed_taxids[child.taxid] = child 
    o_file.close() 
    sys.stdout.write("\r\t%i nodes printed\n" % print_count)
    sys.stdout.flush()
    #Error check
    for taxid in save_taxids:
        if taxid not in printed_taxids:
            sys.stderr.write("ERROR: %s not linked to root\n" % taxid) 
    #End of program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM END TIME: " + time + '\n')
    sys.exit(0)

#################################################################################
if __name__ == "__main__":
    main()
#################################################################################
##################################END OF PROGRAM#################################
#################################################################################

