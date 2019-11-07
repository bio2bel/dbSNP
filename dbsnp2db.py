"""
@author: laurendelong
"""
import urllib.request
import os
import json
import bz2
import re
from dataclasses import dataclass
from typing import List
import itertools
import sqlite3

conn = sqlite3.connect('dbsnp.db')
c = conn.cursor()


@dataclass
class SnpInfo:
    dbsnp_id: str
    assembly_id: str
    gen_type: str
    gene_name: str
    symbol: str
    entrez_id: str
    dna_change: List[str]
    rnas: str
    rna_type: str
    rna_change: List[str]
    proteins: str
    prot_type: str
    aa_change: List[str]


def main():
    # Here we begin the downloading of JSON files from the dbSNP database:

    # Create tables
    c.execute('CREATE TABLE dbsnp_id_MT (dbsnp_id)')
    c.execute('CREATE TABLE dna_change_MT (dna_change, dbsnp_id)')
    c.execute('CREATE TABLE gene_MT (gene_name, symbol, entrez_id)')
    c.execute('CREATE TABLE snp2gene_MT (dbsnp_id, entrez_id)')

    url = 'ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chrMT.json.bz2'
    path = '/home/llong/Downloads/refsnp/refsnp-chrMT.json.bz2'
    if not os.path.exists(path):
        print('Beginning file download with urllib2...')
        urllib.request.urlretrieve(url, path)
        print('...Finished file download with urllib2.')

    id_list = []
    dna_change_list = []
    gene_list = []
    snp2gene_list = []

    # Here we parse through the files:
    print('Now decompressing and reading JSON.bz2 files with *bz2* and *json* ...')
    with bz2.BZ2File(path, 'rb') as f_in:
        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            dbsnp_id = rs_obj['refsnp_id']  # the dbsnp id
            # Make the dbsnp_id table for database
            id_list.append((dbsnp_id,))

            all_ann_list_raw = rs_obj['primary_snapshot_data'][
                'allele_annotations']  # these are the assembly annotations

            if len(all_ann_list_raw) >= 2:  # if it has sufficient info
                assembl_ann_list_raw = all_ann_list_raw[1]['assembly_annotation']  # against each assembly
                if len(assembl_ann_list_raw) != 0:  # if it contains gene info
                    gene_list_raw = assembl_ann_list_raw[0][
                        'genes']  # and each of the genes affected within each assembly
                    if len(gene_list_raw) > 0:
                        # Here I start extracting gene info:
                        for x, y, z in itertools.product(range(len(all_ann_list_raw)),
                                                         range(len(assembl_ann_list_raw)),
                                                         range(len(gene_list_raw))):
                            assembly_id = all_ann_list_raw[x]['assembly_annotation'][y]['seq_id']
                            if assembly_id[0:2] == 'AC':
                                gen_type = 'Complete genomic molecule, usually alternate assembly'
                            elif assembly_id[0:2] == 'NC':
                                gen_type = 'Complete genomic molecule, usually reference assembly'
                            elif assembly_id[0:2] == 'NG':
                                gen_type = 'Incomplete genomic region'
                            elif assembly_id[0:2] == 'NT':
                                gen_type = 'Contig or scaffold, clone-based or WGSa'
                            elif assembly_id[0:2] == 'NW':
                                gen_type = 'Contig or scaffold, primarily WGSa'
                            elif assembly_id[0:2] == 'NZ':
                                gen_type = 'Complete genomes and unfinished WGS data'
                            else:
                                gen_type = ''

                            gene_name = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['name']
                            symbol = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
                            entrez_id = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['id']

                            # Here I make the gene table
                            gene_list.append((gene_name, symbol, entrez_id))
                            # Here I make an intermediary table
                            snp2gene_list.append((dbsnp_id, entrez_id))

                            rna_list_raw = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['rnas']

                            for nuc in rna_list_raw:
                                if 'id' in nuc:
                                    rnas = nuc['id']  # the rna transcript affected by the mutation
                                    if rnas[0:2] == 'NM':
                                        rna_type = 'Protein-coding transcripts (usually curated)'
                                    elif rnas[0:2] == 'NR':
                                        rna_type = 'Non-protein-coding transcripts'
                                    elif rnas[0:2] == 'XM':
                                        rna_type = 'Predicted model protein-coding transcript'
                                    elif rnas[0:2] == 'XR':
                                        rna_type = 'Predicted model non-protein-coding transcript'
                                    else:
                                        rna_type = ''
                                else:
                                    rnas = ''
                                    rna_type = ''
                                if 'product_id' in nuc:
                                    proteins = nuc['product_id']  # the protein affected by the mutation
                                    if proteins[0:2] == 'AP':
                                        prot_type = 'Annotated on AC_ alternate assembly'
                                    elif proteins[0:2] == 'NP':
                                        prot_type = 'Associated with an NM_ or NC_ accession'
                                    elif proteins[0:2] == 'YP':
                                        prot_type = 'Annotated on genomic molecules without an instantiated ' \
                                                    'transcript record '
                                    elif proteins[0:2] == 'XP':
                                        prot_type = 'Predicted model, associated with an XM_ accession'
                                    elif proteins[0:2] == 'WP':
                                        prot_type = 'Non-redundant across multiple strains and species'
                                    else:
                                        prot_type = ''
                                else:
                                    proteins = ''
                                    prot_type = ''

                                # Here I parse through each hgvs entry and assign it to either a nuc. change or a.a. change
                                hgvs_entries = rs_obj['primary_snapshot_data']['placements_with_allele']
                                dna_change = []
                                aa_change = []
                                rna_change = []
                                for entry in hgvs_entries:
                                    for variant in entry['alleles']:
                                        hgvs = re.split(":[cgmnopr].", variant['hgvs'])
                                        if len(hgvs) > 1:
                                            if hgvs[0] == assembly_id:
                                                dna_change.append(hgvs[1])
                                                # Here I make the DNA change table
                                                dna_change_list.append((hgvs[1], dbsnp_id))
                                            elif hgvs[0] == proteins:
                                                aa_change.append(hgvs[1])
                                            elif hgvs[0] == rnas:
                                                rna_change.append(hgvs[1])
                                            else:
                                                continue

    c.executemany('INSERT INTO dbsnp_id_MT VALUES (?)', id_list)
    conn.commit()
    c.executemany('INSERT INTO dna_change_MT VALUES (?,?)', dna_change_list)
    conn.commit()
    c.executemany('INSERT INTO gene_MT VALUES (?,?,?)', gene_list)
    conn.commit()
    c.executemany('INSERT INTO snp2gene_MT VALUES (?,?)', snp2gene_list)
    conn.commit()

    print("Finished writing files to database")
    conn.close()


if __name__ == '__main__':
    main()
