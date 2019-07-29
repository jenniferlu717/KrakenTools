# Kraken Tools
For news and updates, refer to the github page: https://github.com/jenniferlu717/KrakenTools/

Kraken-Tools is a suite of scripts to be used for post-analysis of 
Kraken/KrakenUniq/Kraken2/Bracken results. Please cite the relevant paper
if using Kraken-Tools with any of the listed programs. 

# Running Scripts:
No installation required. 
All scripts are scripts that are run on the command line as follows.

Users can make scripts executable by running

    chmod +x myscript.py
    ./myscript.py -h 

## extract\_kraken\_reads.py

This program extract reads classified at any user-specified taxonomy IDs. User
must specify the Kraken output file, the sequence file(s), and at least one
taxonomy ID. Additional options are specified below.

1. USAGE/OPTIONS
    
    python extract\_kraken\_reads.py
    *   `-k, --kraken MYFILE.KRAKEN.............`Kraken output file
    *   `-s, -s1, -1, -U SEQUENCE.FILE..........`FASTA/FASTQ sequence file (may be gzipped)
    *   `-s2, -2 SEQUENCE2.FILE.................`FASTA/FASTQ sequence file (for paired reads, may be gzipped)
    *   `-o, --output OUTPUT.FASTA..............`output FASTA file with extracted seqs
    *   `-t, --taxid TID TID2 etc...............`list of taxonomy IDs to extract (separated by spaces)        

    Optional:
    *   `-r, --report MYFILE.KREPORT.............`Kraken report file (required if specifying --include-children or --include-parents)
    *   `--include-children......................`include reads classified at more specific levels than specified taxonomy ID levels. 
    *   `--include-parents.......................`include reads classified at all taxonomy levels between root and the specified taxonomy ID levels.
    *   `--max #.................................`maximum number of reads to save.
    *   `--append................................`if output file exists, appends reads
    *   `--noappend..............................`[default] rewrites existing output file
    
2. INPUT FILES: SEQUENCE FILES

    Input sequence files must be either FASTQ or FASTA files. Input files
    can be gzipped or not. The program will automatically detect whether
    the file is gzipped and whether it is FASTQ or FASTA formatted based on
    the first character in the file (">" for FASTA, "@" for FASTQ)

3. PAIRED INPUT/OUTPUT
    
    Users that ran Kraken using paired reads should input both read files into extract_kraken_reads.py as follows:
        ./extract\_kraken\_reads.py -k myfile.kraken -s1 read1.fq -s2 reads2.fq


4. --include-parents/--include-children flags
    
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

## combine\_kreports.py 

This script combines multiple Kraken reports into a combined report file.

1. USAGE/OPTIONS
    
    `python complete_kreports.py`
    *    `-r 1.KREPORT 2.KREPORT........................`Kraken-style reports to combine 
    *    `-o COMBINED.KREPORT...........................`Output file 

    Optional:
    *   `--display-headers..............................`include headers describing the samples and columns [all headers start with #]
    *   `--no-headers...................................`do not include headers in output
    *   `--sample-names.................................`give abbreviated names for each sample [default: S1, S2, ... etc]
    *   `--only-combined................................`output uses exact same columns as a single Kraken-style report file. Only total numbers for read counts and percentages will be used. Reads from individual reports will not be included.

2. OUTPUT 
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

## kreport2krona.py 

This program takes a Kraken report file and prints out a krona-compatible TEXT file

1. USAGE/OPTIONS
    
    `python kreport2krona.py`
    *    `-r/--report ................`Kraken report file 
    *    `-o/--output ................`Output Krona text file
    
    Optional:
    *    `--no-intermediate-ranks.....`only output standard levels [D,P,C,O,F,G,S] 
    *    `--intermediate-ranks........`[default] include non-standard levels

2. EXAMPLE USAGE 
    
        kraken2 --db KRAKEN2DB --threads THREADNUM --report MYSAMPLE.KREPORT --paired SAMPLE\_1.FASTA SAMPLE\_2.FASTA > MYSAMPLE.KRAKEN2
        python kreport2krona.py -r MYSAMPLE.KREPORT -o MYSAMPLE.krona 
        ktImportText MYSAMPLE.krona -o MYSAMPLE.krona.html
    
    Krona information: see https://github.com/marbl/Krona. 

