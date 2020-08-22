#!/usr/bin/python3
"""
This script reads eCIS-screen output file, finds input GenBank file for each found CIS, 
finds contig containing CIS genes and reports all gene names in CIS locus and their neighbors, 
with coordinates and protein sequences. 
"""
import os
import sys
import argparse
import gzip
import Bio
from Bio import GenBank #, SeqIO, SeqFeature
#from Bio.SeqRecord import SeqRecord
from collections import defaultdict


def get_args():
    """Returns command-line arguments"""
    desc = '''This script reports genes encoding CIS components and their neighbors.'''
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', dest='in_file', type=str, default=None, help='Path to eCIS-screen output file')
    parser.add_argument('-g', dest='gbk_dir', type=str, default=None, help='Path to genomes directory (gzipped GenBank files for eCIS-screen)')
    parser.add_argument('-o', dest='out_file', type=str, default=None, help='Output file name')
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return args


def read_cis_summary(in_file):
    """Reads eCIS-screen output file """
    result = defaultdict(dict)
    with open(in_file, 'r') as infile:
        # Read first line and get list of afp families
        header = infile.readline()
        protein_families = header.rstrip('\n\r').split('\t')[1:]
        for line in infile:
            # Skip all comment lines
            if line.startswith('#'):
                continue
            row = line.rstrip('\n\r').split()
            # First element contains joined genome assembly ID and contig ID. Let's separate them.
            assembly_id = row[0].split('_')[0]
            contig_id = '_'.join(row[0].split('_')[1:])
            # For each contig, make a dictinary with afp family name as kay and list of proteins as value
            result[assembly_id][contig_id] = {}
            for fam_index, family in enumerate(row[1:]):
                if family == '-':
                    continue
                for protein in family.split(','):
                    result[assembly_id][contig_id][protein] = protein_families[fam_index]
    return result


def get_neighbors(assembly, seq_record, targets):
    """ Finds genes of interest in GenBank sequence record"""
    # This parameter defines how many genes we take upstream and downstream of the locus
    neighbors_count = 5

    all_genes = {}
    # In a GenBank sequnce record, find all CDS features
    for feature in seq_record.features:
        if feature.key == 'CDS':
            # Parse location string
            if feature.location.startswith('complement'):
                gene_strand = -1
                location = feature.location[11:-1]
            else:
                gene_strand = 1
                location = feature.location
            gene_start, gene_end = location.split('..')
            gene_start = int(gene_start)
            gene_end = int(gene_end)
            # Print warning if there are more than one genes with the same start. Just in case. Normally, it should not happen.
            if gene_start in all_genes:
                print('WARNING: gene with start', gene_start, 'already exists', str(all_genes[gene_start]))
            # Search for protein sequence, locus tag and function in qualifiers 
            product = ''
            locus_tag = ''
            translation = ''
            for qualifier in feature.qualifiers:
                if qualifier.key == '/locus_tag=':
                    locus_tag = qualifier.value[1:-1]
                elif qualifier.key == '/product=':
                    product = qualifier.value[1:-1]
                elif qualifier.key == '/translation=':
                    translation = qualifier.value[1:-1]
            # Put gene data into dictionary: assembly ID, contig ID, locus tag, start, end, strand, function, protein sequence
            all_genes[gene_start] = [assembly, 
                seq_record.accession[0],
                locus_tag,
                gene_start,
                gene_end,
                gene_strand,
                product,
                translation
                ]
    # Sort all genes by start coordinate
    all_genes_sorted = []
    for gene_start in sorted(all_genes.keys()):
        all_genes_sorted.append(all_genes[gene_start])
    # In the sorted list of genes, find CIS gene indices
    cis_indices = []
    for gene_index, gene in enumerate(all_genes_sorted):
        if gene[2] in targets:
            cis_indices.append(gene_index)
    # If CIS genes not found, return empty list
    if not cis_indices:
        print('No genes found in', seq_record.accession[0], assembly)
        return []
    # Find the first and the last gene of interest
    min_index = min(cis_indices) - 5
    if min_index < 0:
        min_index = 0
    max_index = max(cis_indices) + 5
    if max_index >= len(all_genes_sorted):
        max_index = len(all_genes_sorted) - 1
    result = []
    # Find all genes of interest, put into list and return
    for gene_index, gene in enumerate(all_genes_sorted):
        if gene_index >= min_index and gene_index <= max_index:
            # Add gene index (gene number in the contig)
            gene.append(gene_index)
            if gene_index in cis_indices:
                # For CIS genes, add predicted afp family 
                gene.append(targets[gene[2]])
            else:
                # For neighbor genes, add "other"
                gene.append('other')
            result.append('\t'.join([str(x) for x in gene]))
    return result


def main():
    """ Main function"""
    # Get command-line arguments
    args = get_args()
    # If any of the arguments is missing, fill in with default values
    if args.in_file is None:
        in_file = 'ecis-screen_summary.txt'
    else:
        in_file = args.in_file
    if args.out_file is None:
        out_file = 'cis_locus_genes.txt'
    else:
        out_file = args.out_file
    if args.gbk_dir is None:
        gbk_dir = 'genomes'
    else:
        gbk_dir = args.gbk_dir
    # read eCIS-screen output file
    cis_data = read_cis_summary(in_file)
    candidate_proteins = []
    with open(out_file, 'w') as outfile:
        for assembly in cis_data:
            # Make GenBank file path
            gbk_file = os.path.join(gbk_dir, assembly, assembly + '_genomic.gbff.gz')
            print ('Now reading', gbk_file)
            for seq_record in GenBank.parse(gzip.open(gbk_file, 'rt')):
                if seq_record.accession[0] not in cis_data[assembly]:
                    continue
                print('Sequence found', seq_record.accession[0])
                # Get list of genes
                neighbors = get_neighbors(assembly, seq_record, cis_data[assembly][seq_record.accession[0]])
                # Write list of genes to the output file
                outfile.write('\n'.join(neighbors) + '\n\n')


if __name__=='__main__':
    main()
