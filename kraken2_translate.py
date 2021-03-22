#!/usr/bin/python3
#################################################################
# kraken2_translate.py allows users to combine a kraken2 report and
# classification file to get a full taxonomic classification for each read.
# Copyright (C) 2019-2020 Zachary Munro, zacharymunro2@gmail.com
#
# This file is part of Kraken-Tools.
# Kraken-Tools is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the license, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
####################################################################
# Zachary Munro, zacharymunro2@gmail.com
# Updated: 09/22/2020
#
# This program reads in a kraken2 report and kraken2 classification
# output file and it writes out a CSV with columns:
#   read_id, read_length, taxonomy
#
# Parameters:
#   --classification, .....local path to classification file
#   --report, .............local path to report file
#   --output, ..................name of file to output result to
#   --mpa-format, ..............whether to only use mpa taxons
# Output file format (comma separated)
#   - read name/id
#   - read length in bases
#   - complete taxonomic classification produced by kraken2
#######################################################################
import os
import time
from pathlib import Path
from typing import Dict, List
import argparse


# Unofficial kraken report column names
KRAKEN2_REPORT_HEADERS = [
    "percent_under_taxon",
    "count_under_taxon",
    "count_at_taxon",
    "rank_id",
    "ncbi_id",
    "name",
]

# The lower the rank code, the lower in the taxonomic tree that taxon occurs
TAXON_ORDER = {"R": 0, "D": 1, "K": 2, "P": 3, "C": 4, "O": 5, "F": 6, "G": 7, "S": 8}


class ReportRow:
    """ReportRow encodes the rank data for a report row
        e.g. For a rank of R1:
            rank code: R
            subrank_code: 1
            rank_id: R1
    """

    rank_code: int
    subrank_code: int
    rank_id: str
    ncbi_id: str
    name: str

    def __init__(self, row):
        self.rank_code = TAXON_ORDER[row["rank_id"][0]]
        if row["rank_id"] in TAXON_ORDER.keys():
            self.subrank_code = 0
        else:
            self.subrank_code = int(row["rank_id"][1:])
        self.rank_id = row["rank_id"]
        self.ncbi_id = row["ncbi_id"]
        self.name = row["name"]

    def __lt__(self, other):
        # p1 < p2 calls p1.__lt__(p2)
        if self.rank_code == other.rank_code:
            return self.subrank_code < other.subrank_code
        return self.rank_code < other.rank_code

    def __gt__(self, other):
        # p1 > p2 calls p1.__gt__(p2)
        if self.rank_code == other.rank_code:
            return self.subrank_code > other.subrank_code
        return self.rank_code > other.rank_code

    def __eq__(self, other):
        # p1 == p2 calls p1.__eq__(p2)
        if self.rank_code == other.rank_code:
            return self.subrank_code == other.subrank_code
        return False

    def __ge__(self, other):
        # p1 >= p2 calls p1.__ge__(p2)
        return self > other or self == other


def report_to_data_rows(report_file: str) -> List[Dict[str, str]]:
    """report_to_data_rows parses the kraken2 report file into a list of the row data
    Args:
        report_file (str) - path to kraken2 report file
    Returns:
        (List[Dict[str, str]]) list of Dictionaries containing the column names mapped to the data
            in that column for each row
    """
    all_parsed_lines = []
    with open(report_file, "r") as report_f:
        for line in report_f.readlines():
            parsed_line = line.strip().split(maxsplit=5)
            all_parsed_lines.append(parsed_line)
    data = [dict(zip(KRAKEN2_REPORT_HEADERS, row)) for row in all_parsed_lines]
    return data


def build_taxonomies_from_kraken2_report(
    row_dicts: List[Dict[str, str]], mpa_format: bool
) -> Dict[str, List[str]]:
    """build_taxonomies_from_kraken2_report uses the lines of the kraken2 report to create a mapping
    of read_id to complete taxonomic classification string
    Args:
        row_dicts (List[Dict[str, str]]) - list of the rows from the kraken report represented as
            dictionaries, which map the column name to the data in that row,column
        mpa_format (bool) - whether only the mpa taxons should by used
    Returns:
        (Dict[str, List[str]]) mapping of Read ID to an ordered list of the the taxon
            classifications for that read
    """
    ncbi_id_to_taxonomy_map = {}
    taxonomic_hierarcy = []  # stack
    for row in row_dicts:
        # Ignore unclassified rows for now
        if row["rank_id"][0] not in TAXON_ORDER.keys():
            continue

        current_row = ReportRow(row)

        # Starting from root
        if len(taxonomic_hierarcy) == 0:
            taxonomic_hierarcy.append(current_row)
            ncbi_id_taxonomy = [
                f"{x.rank_id}__{x.name}"
                for x in taxonomic_hierarcy
                if x.rank_id in TAXON_ORDER.keys()
            ]
            ncbi_id_to_taxonomy_map[int(current_row.ncbi_id)] = ncbi_id_taxonomy
            continue

        # Able to use comparison operators to see if new row is lower or higher in taxonomic
        #   hierarchy because of overridden comparison operators in the ReportRow class
        if taxonomic_hierarcy[-1] < current_row:
            taxonomic_hierarcy.append(current_row)
        elif taxonomic_hierarcy[-1] == current_row:
            taxonomic_hierarcy.pop()
            taxonomic_hierarcy.append(current_row)
        else:
            # Need to pop from current taxonomy hierarchy list
            while len(taxonomic_hierarcy) != 0 and taxonomic_hierarcy[-1] >= current_row:
                taxonomic_hierarcy.pop()
            taxonomic_hierarcy.append(current_row)

        if mpa_format:
            ncbi_id_taxonomy = [
                f"{x.rank_id}__{x.name}"
                for x in taxonomic_hierarcy
                if x.rank_id in TAXON_ORDER.keys() and x.subrank_code == 0
            ]
        else:
            ncbi_id_taxonomy = [f"{x.rank_id}__{x.name}" for x in taxonomic_hierarcy]
        ncbi_id_to_taxonomy_map[int(current_row.ncbi_id)] = ncbi_id_taxonomy
    return ncbi_id_to_taxonomy_map


def parse_kraken2_classification(
    ncbi_id_to_taxonomy_map: Dict[int, List[str]], classification_file: str, output_file: str
) -> None:
    """parse_kraken2_classification writes out the final translated file
    Args:
        ncbi_id_to_taxonomy_map (Dict[int, List[str]]) - dictionary mapping ncbi_id integer to the
            taxonomy reported by kraken2 split into an ordered list by taxon
        classification_file (str) - path to classification file
        output_file (str) - path to output file
    """
    Path(output_file).touch()
    output_data_lines = []
    with open(classification_file, "r") as classification_f:
        for line in classification_f.readlines():
            line_parts = line.strip().split()
            read_id = line_parts[1]
            ncbi_id = line_parts[2]
            read_length = line_parts[3]
            if str(ncbi_id) == "0":
                taxonomy = "unclassified"
            else:
                taxonomy = "|".join(ncbi_id_to_taxonomy_map[int(ncbi_id)])
            output_data_lines.append(",".join([read_id, read_length, taxonomy]))
    with open(output_file, "w") as out_f:
        out_f.write("read_id,read_length,taxonomy")
        for line in output_data_lines:
            out_f.write("\n" + line)


def translate_kraken2_output(
    classification_file: str, report_file: str, output_file: str, mpa_format: bool
) -> None:
    """translate_kraken2_output uses kraken2 classification file and report file to create a csv
    similar to kraken1's kraken-translate output
    Args:
        classification_file (str) - path to kraken2 classification file
        report_file (str) - path to kraken2 report file
        output_file (str) - path to output file
        mpa_format (bool) - whether only the mpa taxons should by used
    """
    df = report_to_data_rows(report_file)
    taxonomy_map = build_taxonomies_from_kraken2_report(df, mpa_format)
    parse_kraken2_classification(taxonomy_map, classification_file, output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--classification",
        dest="classification_file",
        type=str,
        help="path to kraken2 classification file",
    )
    parser.add_argument(
        "--report", dest="report_file", type=str, help="path to kraken2 report file"
    )
    parser.add_argument(
        "--mpa-format",
        dest="mpa_format",
        action="store_true",
        help="whether to only output standard taxons",
    )
    parser.add_argument("--output", type=str, help="path to output destination")
    args = parser.parse_args()
    translate_kraken2_output(
        args.classification_file, args.report_file, args.output, args.mpa_format
    )

