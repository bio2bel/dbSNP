import urllib.request
import os
import json
import bz2
import re
from dataclasses import dataclass
from typing import List
import itertools


@dataclass
class SnpInfo:
    dbsnp_id: str
    assembly_id: str
    gene_name: str
    gene_abbr: str
    gene_id: str
    dna_change: List[str]
    rnas: str
    rna_change: List[str]
    proteins: str
    aa_change: List[str]

def snp2csv(cl, out):
    print(cl.dbsnp_id,
          cl.assembly_id,
          cl.gene_name,
          cl.gene_abbr,
          cl.gene_id,
          cl.dna_change,
          cl.rnas,
          cl.rna_change,
          cl.proteins,
          cl.aa_change,
          sep='\t',
          file=out
          )


def main():
    # Here we begin the downloading of JSON files from the dbSNP database:

    url = 'ftp://ftp.ncbi.nih.gov/snp/latest_release/JSON/refsnp-chrY.json.bz2'
    path = '/home/llong/Downloads/refsnp-chrY.json.bz2'
    if not os.path.exists(path):
        print('Beginning file download with urllib2...')
        urllib.request.urlretrieve(url, path)
        print('...Finished file download with urllib2.')

    # Here we parse through the files:
    print('Now decompressing and reading JSON.bz2 files with *bz2* and *json* ...')
    with bz2.BZ2File(path, 'rb') as f_in, open('/home/llong/Downloads/refsnp-chrY.csv', 'w') as output:
        print('dbsnp_id', 'assembly_id', 'gene_name', 'gene_abbr', 'gene_id',
              'dna_change', 'rna', 'rna_change', 'proteins', 'aa_change', sep='\t', file=output)

        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            dbsnp_id = rs_obj['refsnp_id']  # the dbsnp id

            all_ann_list_raw = rs_obj['primary_snapshot_data']['allele_annotations']  # these are the assembly annotations

            if len(all_ann_list_raw) >= 2:  # if it has sufficient info
                assembl_ann_list_raw = all_ann_list_raw[1]['assembly_annotation']  # against each assembly
                if len(assembl_ann_list_raw) != 0: # if it contains gene info
                    gene_list_raw = assembl_ann_list_raw[0]['genes']  # and each of the genes affected within each assembly
                    if len(gene_list_raw) > 0:
                        # Here I start extracting gene info:
                        for x, y, z in itertools.product(range(len(all_ann_list_raw)),
                                                         range(len(assembl_ann_list_raw)),
                                                         range(len(gene_list_raw))):
                            assembly_id = all_ann_list_raw[x]['assembly_annotation'][y]['seq_id']
                            gene_name = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['name']
                            gene_abbr = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
                            gene_id = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['id']
                            rna_list_raw = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['rnas']

                            for nuc in rna_list_raw:
                                if 'id' in nuc:
                                    rnas = nuc['id'] # the rna transcript affected by the mutation
                                else:
                                    rnas = ''
                                if 'product_id' in nuc:
                                    proteins = nuc['product_id']  # the protein affected by the mutation
                                else:
                                    proteins = ''

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
                                            elif hgvs[0] == proteins:
                                                aa_change.append(hgvs[1])
                                            elif hgvs[0] == rnas:
                                                rna_change.append(hgvs[1])
                                            else:
                                                continue

                                max_list = max(len(dna_change), len(rna_change), len(aa_change))
                                for i in [dna_change, rna_change, aa_change]:
                                    diff = abs(max_list - len(i))
                                    i.extend(list(itertools.repeat('', diff)))

                                for n in range(max_list):
                                    snp_infos = SnpInfo(dbsnp_id,
                                                        assembly_id,
                                                        gene_name,
                                                        gene_abbr,
                                                        gene_id,
                                                        dna_change[n],
                                                        rnas,
                                                        rna_change[n],
                                                        proteins,
                                                        aa_change[n]
                                                        )

                                    snp2csv(snp_infos, output)

    print("Finished writing files to CSV.")


if __name__ == '__main__':
    main()