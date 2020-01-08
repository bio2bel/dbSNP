"""
@author: laurendelong
"""
import urllib.request
import os
import json
import bz2
import re
import itertools
import sqlite3
import time

conn = sqlite3.connect('dbsnp.db')
c = conn.cursor()

chroms = list(range(1, 23)) + ['X', 'Y', 'MT']

""" dbsnp2db downloads data from dbSNP on NCBI, parses it, and organizes the information into a relational database
which is easier to query for large data with many SNPs """

# Here we create tables for our database
c.execute('CREATE TABLE dbsnp_id (dbsnp_id, chromosome)')
c.execute('CREATE TABLE snp2assembly (dbsnp_id, assembly_id)')
c.execute('CREATE TABLE assembly_id (assembly_id, genome_type)')
c.execute('CREATE TABLE snp2gene (dbsnp_id, entrez_id)')
c.execute('CREATE TABLE gene (gene_name, symbol, entrez_id)')
c.execute('CREATE TABLE dna_change (dna_change, dbsnp_id)')
c.execute('CREATE TABLE rna (rnas, rna_type, entrez_id)')
c.execute('CREATE TABLE rna_change (rna_change, dbsnp_id, rnas)')
c.execute('CREATE TABLE proteins (proteins, prot_type, rnas, entrez_id)')
c.execute('CREATE TABLE aminoacid_change (aa_change, dbsnp_id, proteins)')

def main(chrom):
    t0 = time.time()
    # Here we begin the downloading of JSON files from the dbSNP database:

    url = 'ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chr{}.json.bz2'.format(chrom)
    path = '/home/llong/Downloads/refsnp/refsnp-chr{}.json.bz2'.format(chrom)
    if not os.path.exists(path):
        print('Beginning file download with urllib2...')
        urllib.request.urlretrieve(url, path)
        print('...Finished file download with urllib2.')

    id_list = []
    snp2assembly_list = []
    assembly_list = []
    dna_change_list = []
    rna_change_list = []
    aa_change_list = []
    gene_list = []
    snp2gene_list = []
    rna_list = []
    prot_list = []

    # Here we parse through the files:
    print('Now decompressing and reading JSON.bz2 files from chromosome {} with *bz2* and *json* ...'.format(chrom))
    with bz2.BZ2File(path, 'rb') as f_in:
        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            dbsnp_id = rs_obj['refsnp_id']  # the dbsnp id
            # Make the dbsnp_id table for database
            id_list.append((dbsnp_id, chrom))

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
                            # Here I make the snp2assem table
                            snp2assem_tuple = (dbsnp_id, assembly_id)
                            if snp2assem_tuple not in snp2assembly_list:
                                snp2assembly_list.append(snp2assem_tuple)
                            else:
                                continue

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

                            # Here I make the assembly table
                            assembly_tuple = (assembly_id, gen_type)
                            if assembly_tuple not in assembly_list:
                                assembly_list.append(assembly_tuple)

                            gene_name = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['name']
                            symbol = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
                            entrez_id = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['id']

                            # Here I make the gene table
                            gene_tuple = (gene_name, symbol, entrez_id)
                            if gene_tuple not in gene_list:
                                gene_list.append(gene_tuple)
                            # Here I make an intermediary table
                            id_tuple = (dbsnp_id, entrez_id)
                            if id_tuple not in snp2gene_list:
                                snp2gene_list.append(id_tuple)

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
                                # Here I make the rna transcript table
                                rna_tuple = (rnas, rna_type, entrez_id)
                                if rna_tuple not in rna_list:
                                    rna_list.append(rna_tuple)

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
                                # Here I make the protein table
                                prot_tuple = (proteins, prot_type, rnas, entrez_id)
                                if prot_tuple not in prot_list:
                                    prot_list.append(prot_tuple)

                                # Here I parse through each hgvs entry and assign it to either a nuc. change or a.a. change
                                hgvs_entries = rs_obj['primary_snapshot_data']['placements_with_allele']
                                for entry in hgvs_entries:
                                    for variant in entry['alleles']:
                                        hgvs = re.split(":[cgmnopr].", variant['hgvs'])
                                        if len(hgvs) > 1:
                                            if hgvs[0] == assembly_id:
                                                # Here I make the DNA change table
                                                dna_tuple = (hgvs[1], dbsnp_id)
                                                if dna_tuple not in dna_change_list:
                                                    dna_change_list.append(dna_tuple)
                                            elif hgvs[0] == proteins:
                                                # Here I make the Amino Acid change table
                                                aa_tuple = (hgvs[1], dbsnp_id, proteins)
                                                if aa_tuple not in aa_change_list:
                                                    aa_change_list.append(aa_tuple)
                                            elif hgvs[0] == rnas:
                                                # Here I make the RNA change table
                                                rna_tuple = (hgvs[1], dbsnp_id, rnas)
                                                if rna_tuple not in rna_change_list:
                                                    rna_change_list.append(rna_tuple)
                                            else:
                                                continue

    c.executemany('INSERT INTO dbsnp_id VALUES (?, ?)', id_list)
    conn.commit()
    c.executemany('INSERT INTO snp2assembly VALUES (?,?)', snp2assembly_list)
    conn.commit()
    c.executemany('INSERT INTO assembly_id VALUES (?,?)', assembly_list)
    conn.commit()
    c.executemany('INSERT INTO snp2gene VALUES (?,?)', snp2gene_list)
    conn.commit()
    c.executemany('INSERT INTO gene VALUES (?,?,?)', gene_list)
    conn.commit()
    c.executemany('INSERT INTO dna_change VALUES (?,?)', dna_change_list)
    conn.commit()
    c.executemany('INSERT INTO rna VALUES (?,?,?)', rna_list)
    conn.commit()
    c.executemany('INSERT INTO rna_change VALUES (?,?,?)', rna_change_list)
    conn.commit()
    c.executemany('INSERT INTO proteins VALUES (?,?,?,?)', prot_list)
    conn.commit()
    c.executemany('INSERT INTO aminoacid_change VALUES (?,?,?)', aa_change_list)
    conn.commit()
    t1 = time.time()
    totaltime = t1 - t0
    print("Finished writing files from chromosome {} to the database in {} seconds.".format(chrom, totaltime))


if __name__ == '__main__':
    for chrom in chroms:
        main(chrom)

