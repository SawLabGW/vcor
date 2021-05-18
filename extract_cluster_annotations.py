#!/usr/bin/env python

import re
import argparse
import pandas as pd
from Bio import SeqIO

def get_annotation(cluster, genbank, strain):
    df = pd.read_csv(cluster, sep="\t")
    recs = [i for i in SeqIO.parse(genbank, "genbank")]
    gene_coords = {}
    for rec in recs:
        for f in rec.features:
            if f.type == "CDS":
                if 'pseudo' not in f.qualifiers.keys():
                    ## Needed to truncate the first amino acid residue because Anvio somehow didn't preserve the same coordinates as the original Genbank file
                    ## needed to start from [1:] instead of all of it
                    x = f.qualifiers['translation'][0].replace("-", "")
                    modified = x[1:]
                    gene_coords[modified] = (f.qualifiers['locus_tag'][0], f.location, f.qualifiers['product'][0])
    
    hits = []
    for i, j, k in zip(df.aa_sequence, df.gene_cluster_id, df.genome_name):
        x = i.replace("-", "")
        mod_aa = x[1:]
        if mod_aa in gene_coords:
            hits.append((gene_coords[mod_aa][0], gene_coords[mod_aa][2], j, strain))

    sets = set(hits)

    for j in sorted(sets):
        print('{0}\t{1}\t{2}\t{3}'.format(j[2], j[3], j[0], j[1]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script parses for annotation of gene clusters identified by Anvi'o")
    parser.add_argument("-c", "--gene_cluster_info", required=True, help="gene cluster info file")
    parser.add_argument("-g", "--genbank_file", required=True, help="original Genbank file from NCBI")
    parser.add_argument("-s", "--strain", required=True, help="strain name")
    args = parser.parse_args()
    get_annotation(args.gene_cluster_info, args.genbank_file, args.strain)

