import urllib.request
import os
import json
import bz2
import re
import numpy as np
import itertools

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
        print('dbsnp_id', 'assembly_id', 'dna_change', 'aa_change', 'rna_change', 'gene_name',
              'gene_abbr', 'gene_id', 'protein_access', sep='\t', file=output)
        for line in f_in:
            rs_obj = json.loads(line.decode('utf-8'))
            dbsnp_id = rs_obj['refsnp_id'] # the dbsnp id

            # Here I parse through each hgvs entry and assign it to either a nuc. change or a.a. change
            hgvs_list = []
            aa_list = []
            rna_list = []
            seq_id_list = []
            hgvs_entries = rs_obj['primary_snapshot_data']['placements_with_allele']
            for entry in range(len(hgvs_entries)):
                seq_type = hgvs_entries[entry]['placement_annot']['seq_type']
                if len(hgvs_entries[entry]['alleles']) < 2:
                    hgvs_list.append(np.NaN)
                elif (seq_type == 'refseq_chromosome') or (seq_type == 'refseq_genomic'):
                    seq_id_list.append(hgvs_entries[entry]['seq_id'])
                    # the sequence IDs of the known sequences compared to the SNP
                    hgvs = hgvs_entries[entry]['alleles'][1]['hgvs']
                    hgvs = re.split(":[cgmnopr].", hgvs) # the actual nucleotide change
                    if len(hgvs) > 1:
                        hgvs_list.append(hgvs[1])
                    else:
                        hgvs_list.append(np.NaN)
                elif seq_type == 'refseq_prot':
                    hgvs_prot = hgvs_entries[entry]['alleles'][1]['hgvs']
                    hgvs_prot = re.split(":[cgmnopr].", hgvs_prot) # amino acid change
                    if len(hgvs_prot) > 1:
                        aa_list.append(hgvs_prot[1])
                    else:
                        aa_list.append(np.NaN)
                elif seq_type == 'refseq_mrna':
                    hgvs_mrna = hgvs_entries[entry]['alleles'][1]['hgvs']
                    hgvs_mrna = re.split(":[cgmnopr].", hgvs_mrna) # RNA nucleotide change
                    if len(hgvs_mrna) > 1:
                        rna_list.append(hgvs_mrna[1])
                    else:
                        rna_list.append(np.NaN)
                else:
                    break

            # Here I parse through for the genes affected
            all_ann_list_raw = rs_obj['primary_snapshot_data']['allele_annotations'] # these are the assembly annotations
            gene_list = []
            loci_list = []
            gene_id_list = []
            prot_list = []

            if len(all_ann_list_raw) < 2: # for insufficient information
                gene_list.append(np.NaN)
            else:
                assembl_ann_list_raw = all_ann_list_raw[1]['assembly_annotation']  # against each assembly
                if len(assembl_ann_list_raw) == 0:
                    gene_list.append(np.NaN)
                else:
                    gene_list_raw = assembl_ann_list_raw[0]['genes'] # and each of the genes affected within each assembly
                    if len(gene_list_raw) > 0:
                        for x in range(len(all_ann_list_raw)):
                            for y in range(len(assembl_ann_list_raw)):
                                for z in range(len(gene_list_raw)):
                                    gene_name = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['name']
                                    gene_abbr = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['locus']
                                    gene_id = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['id']
                                    gene = (gene_name + " (" + gene_abbr + " - " + str(gene_id) + ")")  # the genes it falls in
                                    rna_list_raw = all_ann_list_raw[x]['assembly_annotation'][y]['genes'][z]['rnas']
                                    for t in range(len(rna_list_raw)):
                                        if 'product_id' in rna_list_raw[t]:
                                            protein = rna_list_raw[t]['product_id'] # the protein affected by the mutation
                                            prot_list.append(protein)
                                        else:
                                            prot_list.append(np.NaN)
                                    gene_list.append(gene_name)
                                    loci_list.append(gene_abbr)
                                    gene_id_list.append(gene_id)

            # Now I have one string containing dnsnp id, and nine lists containing the other info:

            # Here I make the dbsnp_id list the same length as the max list
            max_list = max(len(seq_id_list), len(hgvs_list), len(aa_list), len(rna_list),
                           len(gene_list), len(loci_list), len(gene_id_list), len(prot_list))
            dbsnp_id = list(itertools.repeat(dbsnp_id, max_list))
            # And now the other lists
            for i in [seq_id_list, hgvs_list, aa_list, rna_list, gene_list, loci_list, gene_id_list, prot_list]:
                diff = abs(max_list - len(i))
                i.extend(list(itertools.repeat(np.NaN, diff)))

            for i in zip(dbsnp_id, seq_id_list, hgvs_list, aa_list, rna_list,
                         gene_list, loci_list, gene_id_list, prot_list):
                print(list(i)[0], list(i)[1], list(i)[2], list(i)[3], list(i)[4], list(i)[5],
                      list(i)[6], list(i)[7], list(i)[8], sep='\t', file=output)
    print("Finished writing files to CSV.")

if __name__ == '__main__':
    main()