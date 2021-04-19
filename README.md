# Kraken Tools
For news and updates, refer to the github page: https://github.com/jenniferlu717/KrakenTools/

KrakenTools is a suite of scripts to be used for post-analysis of 
Kraken/KrakenUniq/Kraken2/Bracken results. Please cite the relevant paper
if using KrakenTools with any of the listed programs. 

## Links to Kraken github pages
1. [Kraken 1](https://github.com/DerrickWood/kraken)
2. [Kraken 2](https://github.com/DerrickWood/kraken2)
3. [KrakenUniq](https://github.com/fbreitwieser/krakenuniq)
4. [Bracken](https://github.com/jenniferlu717/Bracken) 

For issues with any of the above programs, 
please open a github issue on their respective github pages. 
This github repository is dedicated to only the scripts provided here. 

---------------------------------------------------------
## Scripts included in KrakenTools 
1. [extract\_kraken\_reads.py](#extract\_kraken\_readspy)
2. [combine\_kreports.py](#combine\_kreportspy)
3. [kreport2krona.py](#kreport2kronapy)
4. [kreport2mpa.py](#kreport2mpapy)
5. [combine\_mpa.py](#combine\_mpapy)
6. [filter\_bracken\_out.py](#filter\_bracken\_outy)
7. [fix\_unmapped.py](#fix\_unmappedpy)
8. [make\_ktaxonomy.py](#make\_ktaxonomypy)
9. [make\_kreport.py](#make\_kreportpy)

# Running Scripts:
No installation required. 
All scripts are run on the command line as described.

Users can make scripts executable by running

    chmod +x myscript.py
    ./myscript.py -h 

---------------------------------------------------------
# extract\_kraken\_reads.py

This program extract reads classified at any user-specified taxonomy IDs. User
must specify the Kraken output file, the sequence file(s), and at least one
taxonomy ID. Additional options are specified below. As of April 19, 2021, this script is compatible with KrakenUniq/Kraken2Uniq reports.

## 1. extract\_kraken\_reads.py usage/options 
    
`python extract_kraken_reads.py`
*   `-k, --kraken MYFILE.KRAKEN.............`Kraken output file
*   `-s, -s1, -1, -U SEQUENCE.FILE..........`FASTA/FASTQ sequence file (may be gzipped)
*   `-s2, -2 SEQUENCE2.FILE.................`FASTA/FASTQ sequence file (for paired reads, may be gzipped)
*   `-o, --output2 OUTPUT.FASTA.............`output FASTA/Q file with extracted seqs
*   `-t, --taxid TID TID2 etc...............`list of taxonomy IDs to extract (separated by spaces)        

Optional:
*   `-o2, --output2 OUTPUT.FASTA.............`second output FASTA/Q file with extracted seqs (for paired reads)
*   `--fastq-output..........................`Instead of producing FASTA files, print FASTQ files (requires FASTQ input)
*   `--exclude...............................`Instead of finding reads matching specified taxids, finds reads NOT matching specified taxids.
*   `-r, --report MYFILE.KREPORT.............`Kraken report file (required if specifying --include-children or --include-parents)
*   `--include-children......................`include reads classified at more specific levels than specified taxonomy ID levels. 
*   `--include-parents.......................`include reads classified at all taxonomy levels between root and the specified taxonomy ID levels.
*   `--max #.................................`maximum number of reads to save.
*   `--append................................`if output file exists, appends reads
*   `--noappend..............................`[default] rewrites existing output file
    
## 2. extract\_kraken\_reads.py input files

Input sequence files must be either FASTQ or FASTA files. Input files
can be gzipped or not. The program will automatically detect whether
the file is gzipped and whether it is FASTQ or FASTA formatted based on
the first character in the file (">" for FASTA, "@" for FASTQ)

## 3. extract\_kraken\_reads.py paired input/output 
    
Users that ran Kraken using paired reads should input both read files into
extract\_kraken\_reads.py as follows:
        
    extract_kraken_reads.py -k myfile.kraken -s1 read1.fq -s2 reads2.fq
    
Given paired reads, the script requires users to provide two output file 
names to contain extracted reads:

    extract_kraken_reads.py -k myfile.kraken -s1 read1.fq -s2 reads2.fq -o extracted1.fq -o2 extracted2.fq
    
The delimiter (`--delimiter` or `-d`) option has been removed.
    
    `extract_kraken_reads.py -k myfile.kraken ... -o reads_S1.fa -o2 reads_s2.fa

## 4. extract\_kraken\_reads.py --exclude flag
    
By default, reads classified at specified taxonomy IDs will be extracted (and any taxids selected using `--include-parents`/`--include-children`. However, specifying `--exclude` will cause the reads NOT classified at any specified taxonomy IDs.

For example: 
1. `extract_kraken_reads.py -k myfile.kraken ... --taxid 9606 --exclude` ==> extract all reads NOT classified as Human (taxid 9606). 
2. `extract_kraken_reads.py -k myfile.kraken ... --taxid 2 --exclude --include-children` ==> extract all reads NOT classified as Bacteria (taxid 2) or any classification in the Bacteria subtree. 
3. `extract_kraken_reads.py -k myfile.kraken ... --taxid 9606 --exclude --include-parents` ==> extract all reads NOT classified as Human or any classification in the direct ancestry of Human (e.g. will exclude reads classified at the Primate, Chordata, or Eukaryota levels). 
    
## 5. extract\_kraken\_reads.py --include-parents/--include-children flags

By default, only reads classified exactly at the specified taxonomy IDs
will be extracted. Options --include-children and --include parents can be
used to extract reads classified within the same lineage as a specified
taxonomy ID. For example, given a Kraken report containing the following:

        [%]     [reads] [lreads][lvl]   [tid]       [name]
        100     1000    0       R       1           root
        100     1000    0       R1      131567        cellular organisms
        100     1000    50      D       2               Bacteria
        0.95    950     0       P       1224              Proteobacteria
        0.95    950     0       C       1236                Gammaproteobacteria
        0.95    950     0       O       91347                 Enterobacterales
        0.95    950     0       F       543                     Enterobacteriaceae
        0.95    950     0       G       561                       Escherichia
        0.95    950     850     S       562                         Escherichia coli
        0.05    50      50      S1      498388                        Escherichia coli C
        0.05    50      50      S1      316401                        Escherichia coli ETEC

    
1.  `extract_kraken_reads.py  [options] -t 562` ==> 850 reads classified as _E. coli_ will be extracted
2.  `extract_kraken_reads.py  [options] -t 562 --include-parents` ==> 900 reads classified as _E. coli_ or Bacteria will be extracted
3.  `extract_kraken_reads.py  [options] -t 562 --include-children` ==> 950 reads classified as _E. coli_, _E. coli C_, or _E. coli ETEC_ will be extracted
4.  `extract_kraken_reads.py  [options] -t 498388` ==> 50 reads classified as _E. coli C_ will be extracted
5.  `extract_kraken_reads.py  [options] -t 498388 --include-parents` ==> 950 reads classified as _E. coli C_, _E. coli_, or Bacteria will be extracted
6.  `extract_kraken_reads.py  [options] -t 1 --include-children` ==> All classified reads will be extracted 

---------------------------------------------------------
# combine\_kreports.py 

This script combines multiple Kraken reports into a combined report file.

## 1. combine\_kreports.py usage/options
    
`python complete_kreports.py`
*    `-r 1.KREPORT 2.KREPORT........................`Kraken-style reports to combine 
*    `-o COMBINED.KREPORT...........................`Output file 

Optional:
*   `--display-headers..............................`include headers describing the samples and columns [all headers start with #]
*   `--no-headers...................................`do not include headers in output
*   `--sample-names.................................`give abbreviated names for each sample [default: S1, S2, ... etc]
*   `--only-combined................................`output uses exact same columns as a single Kraken-style report file. Only total numbers for read counts and percentages will be used. Reads from individual reports will not be included.

## 2. combine\_kreports.py output
Percentage is only reported for the summed read counts, not for each individual sample. 

The output file therefore contains the following tab-delimited columns:
*    `perc............`percentage of total reads rooted at this clade 
*    `tot_all ........`total reads rooted at this clade (including reads at more specific clades) 
*    `tot_lvl.........`total reads at this clade  (not including reads at more specific clades)
*    `1_all...........`reads from Sample 1 rooted at this clade 
*    `1_lvl...........`reads from Sample 1 at this clade 
*    `2_all...........`""
*    `2_lvl...........`""
*    etc..
*    `lvl_type........`Clade level type (R, D, P, C, O, F, G, S....) 
*    `taxid...........`taxonomy ID of this clade
*    `name............`name of this clade 

---------------------------------------------------------
# kreport2krona.py 

This program takes a Kraken report file and prints out a krona-compatible TEXT file

## 1. kreport2krona.py usage/options
    
`python kreport2krona.py`
*    `-r/--report MYFILE.KREPORT........`Kraken report file 
*    `-o/--output MYFILE.KRONA..........`Output Krona text file
    
Optional:
*    `--no-intermediate-ranks...........`[default]only output standard levels [D,P,C,O,F,G,S] 
*    `--intermediate-ranks..............`include non-standard levels

## 2. kreport2krona.py example usage
    
    kraken2 --db KRAKEN2DB --threads THREADNUM --report MYSAMPLE.KREPORT \
        --paired SAMPLE_1.FASTA SAMPLE_2.FASTA > MYSAMPLE.KRAKEN2
    python kreport2krona.py -r MYSAMPLE.KREPORT -o MYSAMPLE.krona 
    ktImportText MYSAMPLE.krona -o MYSAMPLE.krona.html
    
Krona information: see https://github.com/marbl/Krona. 

## 3. kreport2krona.py example output 

`--no-intermediate-ranks`
        
        6298        Unclassified
        8           k__Bacteria
        4           k__Bacteria     p_Proteobacteria
        6           k__Bacteria     p_Proteobacteria    c__Gammaproteobacteria
        ...


`--intermediate-ranks`

        6298        Unclassified
        79          x__root
        0           x__root     x__cellular_organisms
        8           x__root     x__cellular organisms   k__Bacteria
        4           x__root     x__cellular organisms   k__Bacteria     p__Proteobacteria
        6           x__root     x__cellular organisms   k__Bacteria     p__Proteobacteria   c__Gammaproteobacteria
        ....

---------------------------------------------------------
# kreport2mpa.py 

This program takes a Kraken report file and prints out a mpa (MetaPhlAn) -style TEXT file

## 1. kreport2mpa.py usage/options
    
`python kreport2mpa.py`
*    `-r/--report MYFILE.KREPORT........`Kraken report file 
*    `-o/--output MYFILE.MPA.TXT........`Output MPA-STYLE text file
    
Optional:
*    `--display-header..................`display header line (#Classification, MYFILE.KREPORT) [default: no header]
*    `--no-intermediate-ranks...........`[default] only output standard levels [D,P,C,O,F,G,S] 
*    `--intermediate-ranks..............`include non-standard levels
*    `--read-count......................`[default] use read count for output
*    `--percentages.....................`use percentage of total reads for output 

## 2. kreport2mpa.py example usage
    
    kraken2 --db KRAKEN2DB --threads THREADNUM --report MYSAMPLE.KREPORT \
        --paired SAMPLE_1.FASTA SAMPLE_2.FASTA > MYSAMPLE.KRAKEN2
    python kreport2mpa.py -r MYSAMPLE.KREPORT -o MYSAMPLE.MPA.TXT 
    
## 3. kreport2mpa.py example output 

The output will contain one tab character inbetween the classification and the read count.

`--no-intermediate-ranks/--read-count` 
        
        #Classification                                           SAMPLE.KREPORT
        k__Bacteria                                               36569
        k__Bacteria|p__Proteobacteria                             21001
        k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria      11648
        ... 

`--intermediate-ranks/--read-count`
        
        #Classification                                           SAMPLE.KREPORT
        x__cellular_organisms                                     38462
        x__cellular_organisms|k__Bacteria                         36569 
        x__cellular_organisms|k__Bacteria|p__Proteobacteria       21001
        ... 

---------------------------------------------------------
# combine\_mpa.py 

This program combines multiple outputs from [kreport2mpa.py](#kreport2mpapy).
Files to be combined must have been generated using the same kreport2mpa.py options.

**Important:** 
1. Input files to combine\_mpa.py cannot be a mix of intermediate/no intermediate rank outputs.
2. Input files should be generated using the same Kraken database. 
3. Input files cannot be a mix of read counts/percentage kreport2mpa.py outputs.
**combine**\_**mpa.py** will not test the input files prior to combining. 

If no header is in a given sample file, the program will number the files "Sample #1", "Sample #2", etc. 

## 1. combine\_mpa.py usage/options
    
`python combine_mpa.py`
*    `-i/--input MYFILE1.MPA MYFILE2.MPA.......`Multiple MPA-STYLE text files (separated by spaces) 
*    `-o/--output MYFILE.COMBINED.MPA..........`Output MPA-STYLE text file
    
## 2. combine\_mpa.py example output 

        #Classification                                           Sample #1    Sample #2
        k__Bacteria                                               36569         20034
        k__Bacteria|p__Proteobacteria                             21001         18023
        k__Bacteria|p__Proteobacteria|c__Gammaproteobacteria      11648         15000


---------------------------------------------------------
# filter\_bracken\_out.py

This program takes the output file of a Bracken report and filters the desired taxonomy IDs. 

## 1. filter\_bracken\_out.py usage/options
    
`python filter_bracken_out.py`
*   `-i/--input MYFILE.BRACKEN..........`Bracken output file
*   `-o/--output MYFILE.BRACKEN_NEW.....`Bracken-style output file with filtered taxids
*   `--include TID TID2.................`taxonomy IDs to include in output file [space-delimited]
*   `--exclude TID TID2.................`taxonomy IDs to exclude in output file [space-delimited]

User should specify either taxonomy IDs with `--include` or `--exclude`. If
both are specified, taxonomy IDs should not be in both lists and only
taxonomies to include will be evaluated. 
    
When specifying the --include flag, only lines for the included taxonomy
IDs will be extracted to the filtered output file. The percentages in the
filtered file will be re-calculated so the total percentage in the output 
file will sum to 100%. 

When specifying the --exclude flag alone, all lines in the Bracken file
will be preserved EXCEPT for the lines matching taxonomy IDs provided. 
    
## 2. filter\_bracken\_out.py example usage

This program can be useful for isolating a subset of species to 
better understand the distribution of those particular species in the sample. 
    
For example:

* `python filter_bracken_out.py [options] --include 1764 1769 1773 1781
  39689` will allow users to get the relative percentages of 
  _Mycobacterium avium, marinum, tuberculosis, leprae, and gallinarum_ in their samples.

In other cases, users may want to focus on the distribution of all species
that are NOT the host species in a given sample. This program can then
recalculate percentage distributions for species when excluding reads for
the host.
    
For example, given this output:
        
        name                     tax_id      tax_lvl     kraken....  added...   new.... fraction...
        Homo sapiens             9606        S           ...         ....       999000  0.999000
        Streptococcus pyogenes   1314        S           ...         ....       10      0.000001
        Streptococcus agalactiae 1311        S           ...         ....       5       0.000000
        Streptococcus pneumoniae 1313        S           ...         ....       3       0.000000
        Bordetella pertussis     520         S           ...         ....       20      0.000002
        ...
    
Users may not be interested in the 999,000 reads that are host DNA, but
would rather know the percentage of non-host reads for each of the non-host
species.  Using `python filter_bracken_out.py [options] --exclude 9606`
allows better resolution of the non-host species, allowing each of the
fraction of reads to be recalculated out of 1,000 instead of 1,000,000
reads in the above example. The output would then be:

        name                     tax_id      tax_lvl     kraken....  added...   new.... fraction...
        Streptococcus pyogenes   1314        S           ...         ....       10      0.01000
        Streptococcus agalactiae 1311        S           ...         ....       5       0.05000
        Streptococcus pneumoniae 1313        S           ...         ....       3       0.03000
        Bordetella pertussis     520         S           ...         ....       200     0.20000
        ...
    
---------------------------------------------------------
# fix\_unmapped.py 
When building a Kraken database, an "unmapped.txt" file may be generated if a 
taxonomy for a given sequence is not found. This script can search through 
any accession2taxid files provided and the unmapped.txt file and generate a 
seqid2taxid.map file to be appended to the one already generated. 

## 1. fix\_unmapped.py usage/options
    
`python fix_unmapped.py`
*    `-i/--input unmapped.txt...........`Any file containing accession IDs to map
*    `--accession2taxid REF_FILES.......`Any tab-delimited file with 4 columns, accessions = column 1, taxonomy IDs = column 3
*    `-o/--output OUT_FILE..............`Output tab-delimited file with 2 columns: accessions and taxids
    
Optional:
*    `-r/--remaining....................`file containing any unmapped accession IDs after search [default: `still_unmapped.txt`]

## 2. fix\_unmapped.py example usage
    
    rm *.k2d 
    mv seqid2taxid.map seqid2taxid_1.map
    python fix_unmapped -i unmapped.txt --accession2taxid taxonomy/*accession2taxid -o seqid2taxid_temp.map 
    cat seqid2taxid_1.map seqid2taxid_temp.map
    kraken2-build --build --db . --threads 4

---------------------------------------------------------
# make\_ktaxonomy.py
For future KrakenTools scripts, this program generates a single text file
that contains all of the taxonomy information required. This program is intended to 
generate a single text taxonomy file for any Kraken 1, Kraken 2, or KrakenUniq database.

*Important:* The output of this program does not replace any Kraken database file 
(do not replace your taxo.k2d or .db files). 

## 1. make\_ktaxonomy.py usage/options

`python make_ktaxonomy.py`
*   `--nodes taxonomy/nodes.dmp...........`nodes.dmp file in Kraken DB taxonomy/ folder 
*   `--names taxonomy/names.dmp...........`names.dmp file in Kraken DB taxonomy/ folder 
*   `--seqid2taxid seqid2taxid.map........`seqid2taxid.map file generated by kraken-build/kraken2-build/krakenuniq-build when building the database. This is a 2-column tab-delimited file containing sequence IDs and taxonomy IDs.
*   `-o/--output OUT_FILE.................`Output text file. More details below

The program will inform users if a taxonomy ID is listed in the `seqid2taxid.map` 
file but not in either the `nodes.dmp` or the `names.dmp` files. 

### 2. make\_ktaxonomy.py output file format
The output file is similar to the nodes.dmp/names.dmp file format, but not identical. 
Each of the following columns is separated by a tab-vertical line-tab (e.g. `\t|\t`).

1. taxonomy ID 
2. parent taxonomy ID
3. rank type (R = root, D = domain/superkingdom, P = phylum, etc.)
4. level number (distance from root)
5. name

For ranks outside of the traditional taxonomy ranks (R, D, P, C, O, F, G, S),
the rank type will be assigned based on the closest parent, with a number to specify
distance from that parent. For example, the strains will be labeled with `S1` while 
ranks inbetween Genus and Species will be labeled with `G1, G2, etc`. 
 
Currently, names for each node are selected based on the first name listed in 
the `names.dmp` file or the name designated as `scientific name`. 
`scientific names` will be preferred over all others.    

## 3. make\_ktaxonomy.py required input 
1. taxonomy/nodes.dmp
2. taxonomy/names.dmp
3. seqid2taxid.map 

## 4. KrakenTools scripts requiring make\_ktaxonomy.py output:
1. [make\_kreport.py](#make\_kreportpy)

---------------------------------------------------------
# make\_kreport.py 
This program will generate a kraken-style report file from the kraken output file. 
Currently, this only generates reports for Kraken 1 or Kraken 2. This program
does not currently work for KrakenUniq output files (to be completed in a future project). 

This program requires that users first generate the taxonomy file 
created by [make\_ktaxonomy.py](#make\ktaxonomypy). 

## 1. make\_kreport.py usage/options 
`python make_kreport.py`
*   `-i/-k/--input KRAKEN_FILE........`default Kraken output file (5 tab-delimited columns, taxid in third column)
*   `-t/--taxonomy TAXONOMY_FILE......`output from make\_ktaxonomy.py
*   `-o/--output REPORT_FILE..........`output Kraken report file (6 tab-delimited columns)

Optional
*   `--use-read-len...................`make report using summed read lengths instead of read counts

## 2. make\_kreport.py example
Given a Kraken 2 database `KRAKENDB/` and sample file `EXAMPLE_READS.fq`, 
the following commands can be used to generate a Kraken report file 
with this script.

    python make_ktaxonomy.py --nodes KRAKENDB/taxonomy/nodes.dmp --names KRAKENDB/taxonomy/names.dmp --seqid2taxid KRAKENDB/seqid2taxid.map -o KRAKENDB/mydb_taxonomy.txt 
    kraken2 --db KRAKENDB --threads 4 EXAMPLE_READS.fq > EXAMPLE.kraken2 
    python make_kreport.py -i EXAMPLE.kraken2 -t KRAKENDB/mydb_taxonomy.txt -o EXAMPLE.kreport2 

## 3. make\_kreport.py --use-read-len option
By default, the output Kraken report will list read counts for each taxonomy ID. However,
if all read lengths are not the same, users can add the `--use-read-len` option, which will
result in reporting summed read lengths for each taxon. 

## 4. make\_kreport.py output format
The output format for kreport.py is identical to the format generated by 
`kraken-report` or the `--report` switch with `kraken2`. The output
file contains 6 tab-delimited columns as follows:

1. Percentage of total reads
2. Reads classified within sub-tree
3. Reads classified at this specific node (reads cannot be more specifically classified)
4. Level type (R = root, K = kingdom, P = phylum, etc)
5. Taxonomy ID 
6. Name (preceeded by spaces to indicate distance from root) 


---------------------------------------------------------
# Author Information 
Jennifer Lu
jennifer.lu717@gmail.com
jlu26@jhmi.edu 
Page Updated: 2020/05/10
