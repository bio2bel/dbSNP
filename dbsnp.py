import urllib.request
import os
import json
import bz2
import re
import numpy as np
from dataclasses import dataclass
from typing import List

@dataclass
class SnpInfo:
    dbsnp_id: str
    assembly_id: str
    dna_change: str
    aa_change: str
    rna_change: str
    gene_name: str
    gene_abbr: str
    gene_id: str
    protein_access: List[str]

def main():

    # Here we begin the downloading of JSON files from the dbSNP database:

    url = 'ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chrMT.json.bz2'
    path = '/home/llong/Downloads/refsnp-chrMT.json.bz2'
    if not os.path.exists(path):
        print('Beginning file download with urllib2...')
        urllib.request.urlretrieve(url, path)
        print('...Finished file download with urllib2.')

    # Here we parse through the files:
    print('Now decompressing and reading JSON.bz2 files with *bz2* and *json* ...')
    with bz2.BZ2File(path, 'rb') as f_in, open('/home/llong/Downloads/refsnp-chrMT.csv', 'w') as output:
        print('dbsnp_id', 'assembly_id', 'dna_change', 'aa_change', 'gene_name',
              'gene_abbr', 'gene_id', 'protein_access', sep='\t', file=output)

        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            dbsnp_id = rs_obj['refsnp_id'] # the dbsnp id

            all_ann_list_raw = rs_obj['primary_snapshot_data']['allele_annotations'] # these are the assembly annotations
            gene_list = []

            # Here I parse through each hgvs entry and assign it to either a nuc. change or a.a. change
            hgvs_entries = rs_obj['primary_snapshot_data']['placements_with_allele']
            for entry in hgvs_entries:
                seq_type = entry['placement_annot']['seq_type']
                if len(entry['alleles']) > 1:
                    hgvs = entry['alleles'][1]['hgvs']
                    hgvs = re.split(":[cgmnopr].", hgvs)
                    if len(hgvs) > 1:
                        if (seq_type == 'refseq_chromosome') or (seq_type == 'refseq_genomic'):
                            dna_change = hgvs[1]
                        elif seq_type == 'refseq_prot':
                            aa_change = hgvs[1]
                        elif seq_type == 'refseq_mrna':
                            rna_change = hgvs[1]
                    else:
                        dna_change = ''
                        aa_change = ''
                        rna_change = ''

            snp_infos = []
            if len(all_ann_list_raw) < 2: # for insufficient information

                snp_infos.append(SnpInfo(
                    dbsnp_id=dbsnp_id,
                    assembly_id='',
                    dna_change=dna_change,
                    aa_change=aa_change,
                    rna_change=rna_change,
                    gene_name='',
                    gene_abbr='',
                    gene_id='',
                    protein_access=[],
                ))

            else:
                assembl_ann_list_raw = all_ann_list_raw[1]['assembly_annotation']  # against each assembly
                if len(assembl_ann_list_raw) == 0:
                    gene_list.append(np.NaN)
                    snp_infos.append(SnpInfo(
                        dbsnp_id=dbsnp_id,
                        assembly_id='',
                        dna_change=dna_change,
                        aa_change=aa_change,
                        rna_change=rna_change,
                        gene_name='',
                        gene_abbr='',
                        gene_id='',
                        protein_access=[],
                    ))
                else:

                    gene_list_raw = assembl_ann_list_raw[0]['genes'] # and each of the genes affected within each assembly
                    if len(gene_list_raw) > 0:
                        for x in range(len(all_ann_list_raw)):
                            for y in range(len(assembl_ann_list_raw)):
                                for z in range(len(gene_list_raw)):
                                    gene_name = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['name']
                                    gene_abbr = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
                                    gene_id = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['id']
                                    assembly_id = all_ann_list_raw[x]['assembly_annotation'][y]['seq_id']
                                    rna_list_raw = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['rnas']

                                    proteins = []
                                    for t in range(len(rna_list_raw)):
                                        if 'product_id' in rna_list_raw[t]:
                                            protein = rna_list_raw[t]['product_id'] # the protein affected by the mutation
                                            proteins.append(protein)

                                    snp_infos.append(SnpInfo(
                                        dbsnp_id=dbsnp_id,
                                        assembly_id=assembly_id,
                                        dna_change=dna_change,
                                        aa_change=aa_change,
                                        rna_change='',
                                        gene_name=gene_name,
                                        gene_abbr=gene_abbr,
                                        gene_id=gene_id,
                                        protein_access=proteins,
                                    ))

            for i in zip(snp_infos):
                print(dbsnp_id, assembly_id, dna_change, aa_change, gene_name, gene_abbr, gene_id,
                      sep='\t', file=output)
    print("Finished writing files to CSV.")

if __name__ == '__main__':
    main()