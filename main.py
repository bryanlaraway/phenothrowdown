#AUTHOR: Bryan Laraway
#PROJECT: Comparison of OWLSim and Phenologs for the identification of models of human disease and gene candidates for human disease.
#PURPOSE: This script will call of the functions/methods/scripts for performing the processing required for this analysis.

import json
import urllib.request
import codecs
import time
from socket import *
import os
import re
import csv
import pickle
from decimal import Decimal, getcontext
from numpy import *
from scipy.stats import hypergeom
import math
import matplotlib.pyplot as plt

start_time = time.time()

hu_disease_to_phenotype_hash = {'disease_id': {}}
mouse_genotype_to_phenotype_hash = {'genotype_id': {}}
zfin_genotype_to_phenotype_hash = {'genotype_id': {}}
getcontext().prec = 500
#print(getcontext())
#Selected distinct PANTHER IDs from the NIF/DISCO tables.
#TODO: See about getting these numbers from the Panther table to allow for dynamic updating with file updates.
#total_human_mouse_orthologs = 5625
#total_human_zebrafish_orthologs = 5212
#total_mouse_zebrafish_orthologs = 5210
#TODO: Need to pass phenotype/gene labels for identification, or look them up upon final output.


class main():

    # NOTE: Could either include the fetch code to retrieve the data from the resources,
    # or retrieve them and have the code just open local files, already retrieved.

    # Required table from NIF/DISCO
    tables = [
        'dvp.pr_nlx_151835_1',  # HPO: Annoations:DiseasePhenotypes view
        'dvp.pr_nlx_151835_2',  # HPO: Annoations:Phenotype to gene view
        'dvp.pr_nlx_151835_3',  # HPO: Annoations:Disease to gene
        'dvp.pr_nif_0000_00096_5',  # MGI:MouseGenotypes
        'dvp.pr_nif_0000_00096_6',  # MGI:MousePhenotypes
        'dvp.pr_nif_0000_21427_10',  # ZFIN:Genotype-Phenotype
        'dvp.pr_nif_0000_21427_11',  # ZFIN:OrganismGenotypes
        'dvp.pr_nlx_84521_1'  # PANTHER:Orthologs,
        'dvp.pr_nif_0000_02550_3' # ANIMAL QTL DB: Traits
    ]

    files = {
        'aqtl' : {'file' : 'dvp.pr_nif_0000_02550_3'},
        'mgi' : {'file' : 'dvp.pr_nif_0000_00096_6'}
    }

    ####### SCIGRAPH DATA ASSEMBLY #######

    def _assemble_human_disease_to_phenotype(self, limit=None):
        print('INFO: Assembling human disease to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/diseases.csv'
        inter = 'inter/hpo/human_disease_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human diseases to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                #print(row)
                disease_id = row[0]
                #print(disease_id)
                disease_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/diseases/'+disease_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(disease_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))
                    #print(data)
                    pheno_ids = data['nodes']
                    #print(pheno_ids)
                    for rs in pheno_ids:
                        if disease_id not in hu_disease_to_phenotype_hash:
                            hu_disease_to_phenotype_hash[disease_id] = [rs['id']]
                            #print(hu_disease_to_phenotype_hash[disease_id])
                        else:
                            hu_disease_to_phenotype_hash[disease_id].append(rs['id'])
                            #print(hu_disease_to_phenotype_hash[disease_id])
                except Exception:
                    print('Retrieval of '+disease_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break
                    #
                #print(len(hu_disease_to_phenotype_hash.keys()))
        #print(hu_disease_to_phenotype_hash)
                #if disease_id not in hu_disease_to_phenotype_hash:
                    #hu_disease_to_phenotype_hash['disease_id'] =
        with open(inter, 'wb') as handle:
            pickle.dump(hu_disease_to_phenotype_hash, handle)
        print('INFO: Done assembling human disease to phenotype data.')
        print('INFO: '+str(len(hu_disease_to_phenotype_hash.keys()))+' human diseases processed.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    def _assemble_mouse_genotype_to_phenotype(self, limit=None):
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/genotypes.csv'
        inter = 'inter/mgi/mouse_geno_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' mouse genotypes to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                genotype_id = row[0]
                genotype_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/features/'+genotype_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(genotype_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))
                    #print(data)
                    print(genotype_id)
                    pheno_ids = data['nodes']
                    #print(pheno_ids)
                    for rs in pheno_ids:
                        if genotype_id not in mouse_genotype_to_phenotype_hash:
                            mouse_genotype_to_phenotype_hash[genotype_id] = [rs['id']]
                            #print(mouse_genotype_to_phenotype_hash[genotype_id])
                        else:
                            mouse_genotype_to_phenotype_hash[genotype_id].append(rs['id'])
                            #print(mouse_genotype_to_phenotype_hash[genotype_id])
        #print(hu_disease_to_phenotype_hash)
                #if disease_id not in hu_disease_to_phenotype_hash:
                    #hu_disease_to_phenotype_hash['disease_id'] =
                except Exception:
                    print('Retrieval of '+genotype_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break

        with open(inter, 'wb') as handle:
            pickle.dump(mouse_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling mouse genotype to phenotype data.')
        print('INFO: '+str(len(mouse_genotype_to_phenotype_hash.keys()))+' mouse genotypes present.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    def assemble_zebrafish_genotype_to_phenotype(self, limit=None):
        print('INFO:Assembling zebrafish genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/genotypes.csv'
        inter = 'inter/zfin/zebrafish_geno_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish genotypes to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                genotype_id = row[0]
                genotype_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/features/'+genotype_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(genotype_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))
                    #print(data)
                    pheno_ids = data['nodes']
                    #print(pheno_ids)
                    for rs in pheno_ids:
                        if genotype_id not in zfin_genotype_to_phenotype_hash:
                            zfin_genotype_to_phenotype_hash[genotype_id] = [rs['id']]
                            #print(zfin_genotype_to_phenotype_hash[genotype_id])
                        else:
                            zfin_genotype_to_phenotype_hash[genotype_id].append(rs['id'])
                            #print(zfin_genotype_to_phenotype_hash[genotype_id])
                    #print(len(zfin_genotype_to_phenotype_hash.keys()))
                except Exception:
                    print('Retrieval of '+genotype_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish genotype to phenotype data.')
        print('INFO: '+str(len(zfin_genotype_to_phenotype_hash.keys()))+' zebrafish genotypes present.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    # Will require adjustment once the SciGraph REST call is operational again.
    def assemble_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling zebrafish genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/genotypes.csv'
        inter = 'inter/zfin/zebrafish_geno_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish genotypes to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                genotype_id = row[0]
                genotype_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/features/'+genotype_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(genotype_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))
                    #print(data)
                    pheno_ids = data['nodes']
                    #print(pheno_ids)
                    for rs in pheno_ids:
                        if genotype_id not in zfin_genotype_to_phenotype_hash:
                            zfin_genotype_to_phenotype_hash[genotype_id] = [rs['id']]
                            #print(zfin_genotype_to_phenotype_hash[genotype_id])
                        else:
                            zfin_genotype_to_phenotype_hash[genotype_id].append(rs['id'])
                            #print(zfin_genotype_to_phenotype_hash[genotype_id])
                    #print(len(zfin_genotype_to_phenotype_hash.keys()))
                except Exception:
                    print('Retrieval of '+genotype_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish genotype to phenotype data.')
        print('INFO: '+str(len(zfin_genotype_to_phenotype_hash.keys()))+' zebrafish genotypes present.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    ####### NIF DATA ASSEMBLY #######

    ####### PHENOLOG PHENOTYPE TO GENE #######

    #FIXME: Some of the ZFIN phenotypes are coming through as blanks.

    def assemble_nif_zfin_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling zebrafish phenotype to ortholog data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        inter1 = 'inter/zfin/zebrafish_pheno_gene_hash.txt'
        inter2 = 'inter/zfin/zebrafish_pheno_ortholog_hash.txt'
        zfin_phenotype_to_gene_hash = {}
        zfin_phenotype_to_ortholog_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
            row_count = limit
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, effective_genotype_id, effective_genotype_label, effective_genotype_label_html,
                 intrinsic_genotype_id, intrinsic_genotype_num, intrinsic_genotype_label, intrinsic_genotype_label_html,
                 extrinsic_genotype_id, extrinsic_genotype_label, extrinsic_genotype_label_html, phenotype_id,
                 phenotype_label, phenotype_modifier, implicated_gene_ids, implicated_gene_labels, start_stage_id,
                 start_stage_zfin_id, start_stage_label, end_stage_id ,end_stage_zfin_id, end_stage_label, stages,
                 genomic_background_id, genomic_background_num, genomic_background_label,
                 affected_structure_or_process_1_superterm_id, affected_structure_or_process_1_superterm_name,
                 affected_structure_or_process_1_subterm_id, affected_structure_or_process_1_subterm_name,
                 quality_id, quality_label, affected_structure_or_process_2_superterm_id,
                 affected_structure_or_process_2_superterm_name, affected_structure_or_process_2_subterm_id,
                 affected_structure_or_process_2_subterm_name, environment_id, environment_label, evidence_code_id,
                 evidence_code_symbol, evidence_code_label, publication_id, publication_label, publication_url,
                 taxon_id, taxon_label, v_uid, v_uuid, v_lastmodified) = row
                print('INFO: Processing phenotype '+str(line_counter)+' out of '+str(row_count)+'.')

                if phenotype_id == '' or phenotype_id is None:
                    continue
                elif implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue
                else:
                    print(phenotype_id)
                    #FIXME: Going to need to convert the ZFIN Gene IDs to NCBIGene IDs.
                    #TODO: Need to handle phenotypes with no associated genes.
                    genes = implicated_gene_ids.split(',')
                    print(genes)
                    if phenotype_id not in zfin_phenotype_to_gene_hash:
                        zfin_phenotype_to_gene_hash[phenotype_id] = genes
                        #print(zfin_genotype_to_phenotype_hash[genotype_id])
                        for gene in genes:
                            panther_id = self.get_ortholog(gene, 'inter/panther/panther_zebrafish.txt')
                            if panther_id != 'fail':
                                #print('found ortholog')
                                if phenotype_id not in zfin_phenotype_to_ortholog_hash:
                                    zfin_phenotype_to_ortholog_hash[phenotype_id]= [panther_id]
                                elif panther_id not in zfin_phenotype_to_ortholog_hash[phenotype_id]:
                                    zfin_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)
                            #elif panther_id == 'fail':
                                #print('No ortholog found.')
                    else:
                        for gene in genes:
                            zfin_phenotype_to_gene_hash[phenotype_id].append(gene)
                            #print(zfin_genotype_to_phenotype_hash[genotype_id])
                            #print(len(zfin_genotype_to_phenotype_hash.keys()))
                            print('Repeat phenotype: '+phenotype_id)
                            panther_id = self.get_ortholog(gene, 'inter/panther/panther_zebrafish.txt')
                            if panther_id != 'fail':
                                #print('found ortholog')
                                if phenotype_id not in zfin_phenotype_to_ortholog_hash:
                                    zfin_phenotype_to_ortholog_hash[phenotype_id]= [panther_id]
                                elif panther_id not in zfin_phenotype_to_ortholog_hash[phenotype_id]:
                                    zfin_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break
        with open(inter1, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_ortholog_hash, handle)
        print('INFO: Done assembling zebrafish phenotype to gene data.')
        print('INFO: '+str(len(zfin_phenotype_to_gene_hash.keys()))+' zebrafish phenotypes present.')
        return

    def assemble_nif_mgi_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling mouse phenotype to gene ortholog data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        inter1 = 'inter/mgi/mouse_pheno_gene_hash.txt'
        inter2 = 'inter/mgi/mouse_pheno_ortholog_hash.txt'
        mgi_phenotype_to_gene_hash = {}
        mgi_phenotype_to_ortholog_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' mouse phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.')
            row_count = limit
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (annotation_id, effective_genotype_id, effective_genotype_label, effective_genotype_label_html,
                 intrinsic_genotype_id, intrinsic_genotype_label, intrinsic_genotype_label_html,
                 genomic_variation_complement_id, genomic_variation_complement_label,
                 genomic_variation_complement_label_html, implicated_gene_ids, implicated_gene_labels,
                 implicated_sequence_alteration_ids, implicated_sequence_alteration_labels, genomic_background_id,
                 genomic_background_label, phenotype_id, phenotype_label, phenotype_description_free_text,
                 phenotype_modifier, evidence_code_id, evidence_code_symbol, evidence_code_label, environment_id,
                 environment_label, publication_id, publication_label, publication_url, taxon_id, taxon_label,
                 e_uid, v_uid, v_uuid, v_lastmodified) = row
                print('INFO: Processing phenotype '+str(line_counter)+' out of '+str(row_count)+'.')

                if phenotype_id == '' or phenotype_id is None:
                    continue
                elif implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue

                #print(phenotype_id)
                genes = implicated_gene_ids.split(',')
                #print(genes)
                if phenotype_id not in mgi_phenotype_to_gene_hash:
                    mgi_phenotype_to_gene_hash[phenotype_id] = genes
                    #print(mgi_phenotype_to_gene_hash[phenotype_id])
                    for gene in genes:
                        panther_id = self.get_ortholog(gene,'inter/panther/panther_mouse.txt')
                        if panther_id != 'fail':
                            #print('found ortholog')
                            if phenotype_id not in mgi_phenotype_to_ortholog_hash:
                                mgi_phenotype_to_ortholog_hash[phenotype_id]= [panther_id]
                            elif panther_id not in mgi_phenotype_to_ortholog_hash[phenotype_id]:
                                mgi_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)
                else:
                    for gene in genes:
                        mgi_phenotype_to_gene_hash[phenotype_id].append(gene)
                        #print(mgi_phenotype_to_gene_hash[phenotype_id])
                        #print(len(mgi_phenotype_to_gene_hash.keys()))
                        #print('Repeat phenotype: '+phenotype_id)
                        panther_id = self.get_ortholog(gene, 'inter/panther/panther_mouse.txt')
                        if panther_id != 'fail':
                            #print('found ortholog')
                            if phenotype_id not in mgi_phenotype_to_ortholog_hash:
                                mgi_phenotype_to_ortholog_hash[phenotype_id]= [panther_id]
                            elif panther_id not in mgi_phenotype_to_ortholog_hash[phenotype_id]:
                                mgi_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter1, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling mouse phenotype to gene data.')
        print('INFO: '+str(len(mgi_phenotype_to_gene_hash.keys()))+' mouse phenotypes present.')
        with open(inter2, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_ortholog_hash, handle)
        print('INFO: Done assembling mouse phenotype to ortholog data.')
        return

    def assemble_nif_hpo_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling human phenotype to gene data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/dvp.pr_nlx_151835_2'
        inter1 = 'inter/hpo/human_pheno_gene_hash.txt'
        inter2 = 'inter/hpo/human_pheno_ortholog_hash.txt'
        hpo_phenotype_to_gene_hash = {}
        hpo_phenotype_to_ortholog_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.')
            row_count = limit
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, phenotype_id, phenotype_label, gene_id, gene_num,
                 gene_label, v_uid, v_uuid, v_lastmodified) = row

                print('INFO: Processing human phenotype row '+str(line_counter)+' out of '+str(row_count)+'.')

                if phenotype_id == '' or phenotype_id is None:
                    continue
                elif gene_id == '' or gene_id is None:
                    continue

                gene_id = re.sub('NCBI_gene:', 'NCBIGene:', gene_id)
                print(phenotype_id)

                #print(genes)
                if phenotype_id not in hpo_phenotype_to_gene_hash:
                    hpo_phenotype_to_gene_hash[phenotype_id] = [gene_id]
                    #print(hpo_phenotype_to_gene_hash[genotype_id])
                    panther_id = self.get_ortholog(gene_id,'inter/panther/panther_human.txt')
                    if panther_id != 'fail':
                        #print('found ortholog')
                        hpo_phenotype_to_ortholog_hash[phenotype_id] = [panther_id]
                    #elif panther_id == 'fail':
                        #print('No ortholog found.')
                else:
                    hpo_phenotype_to_gene_hash[phenotype_id].append(gene_id)
                    #print(hpo_phenotype_to_gene_hash[genotype_id])
                    #print(len(hpo_phenotype_to_gene_hash.keys()))
                    print('Repeat phenotype: '+phenotype_id)
                    panther_id = self.get_ortholog(gene_id, 'inter/panther/panther_mouse.txt')
                    if panther_id != 'fail':
                        #print('found ortholog')
                        if phenotype_id not in hpo_phenotype_to_ortholog_hash:
                            hpo_phenotype_to_ortholog_hash[phenotype_id] = [panther_id]
                        elif panther_id not in hpo_phenotype_to_ortholog_hash[phenotype_id]:
                            hpo_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)
                    #elif panther_id == 'fail':
                        #print('No ortholog found.')
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter1, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling human phenotype to gene data.')
        print('INFO: '+str(len(hpo_phenotype_to_gene_hash.keys()))+' human phenotypes present.')
        with open(inter2, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_ortholog_hash, handle)
        print('INFO: Done assembling human phenotype to ortholog data.')
        return

    def assemble_nif_animalqtl_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling animalQTLdb phenotype to gene data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/animalqtldb/dvp.pr_nif_0000_02550_3'
        inter = 'inter/aqtl/aqtl_pheno_gene_hash.txt'
        aqtl_phenotype_to_gene_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' AQTL phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.')
            row_count = limit
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, qtl_id, qtl_num, qtl_label, qtl_url, qtl_source, qtl_type_label, association_type_label,
                 trait_category_label, trait_id, trait_num, trait_label, trait_url, vto_id, vto_label, cmo_id,
                 cmo_label, pto_id, pto_label, gene_id, gene_label, gene_url, gene_num, gene_id_src, taxon_label,
                 taxon_id, breed, organism, map_type, model, test_base, significance, genomic_location, chromosome,
                 start_bp, stop_bp, frame, strand, score, genetic_location, peak_cm, range_cm, combo_flankmarkers,
                 flankmarkers, flankmark_a1, flankmark_a2, flankmark_b1, flankmark_b2, peak_mark, publication_num,
                 publication_url, publication_id, p_value, variance, bayes_value, f_statistics, lod_score,
                 additive_effect, dominance_effect, likelihood_ratio, ls_means, v_uid, v_uuid, v_lastmodified) = row

                print('INFO: Processing phenotype '+str(line_counter)+' out of '+str(row_count)+'.')

                if trait_id == '' or trait_id is None:
                    continue
                elif gene_id == '' or gene_id is None:
                    continue

                print(trait_id)
                #TODO: Separate by species.

                #print(genes)
                if trait_id not in aqtl_phenotype_to_gene_hash:
                    aqtl_phenotype_to_gene_hash[trait_id] = [gene_id]
                    #print(aqtl_phenotype_to_gene_hash[genotype_id])
                else:
                    aqtl_phenotype_to_gene_hash[trait_id].append(gene_id)
                    #print(aqtl_phenotype_to_gene_hash[genotype_id])
                    #print(len(aqtl_phenotype_to_gene_hash.keys()))
                    print('Repeat phenotype: '+trait_id)
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(aqtl_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling human phenotype to gene data.')
        print('INFO: '+str(len(aqtl_phenotype_to_gene_hash.keys()))+' human phenotypes present.')
        return

    ####### OWLSIM GENOTYPE TO PHENOTYPE #######

    def assemble_nif_zfin_genotype_to_phenotype(self, limit=None):
        #TODO: Assuming want to filter out to intrinsic genotypes only?
        # Can filter on extrinsic genotype = ''
        print('INFO:Assembling zebrafish genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        inter = 'inter/zfin/zebrafish_genotype_phenotype_hash.txt'
        zfin_genotype_to_phenotype_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish genotype-phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, effective_genotype_id, effective_genotype_label, effective_genotype_label_html,
                 intrinsic_genotype_id, intrinsic_genotype_num, intrinsic_genotype_label, intrinsic_genotype_label_html,
                 extrinsic_genotype_id, extrinsic_genotype_label, extrinsic_genotype_label_html, phenotype_id,
                 phenotype_label, phenotype_modifier, implicated_gene_ids, implicated_gene_labels, start_stage_id,
                 start_stage_zfin_id, start_stage_label, end_stage_id ,end_stage_zfin_id, end_stage_label, stages,
                 genomic_background_id, genomic_background_num, genomic_background_label,
                 affected_structure_or_process_1_superterm_id, affected_structure_or_process_1_superterm_name,
                 affected_structure_or_process_1_subterm_id, affected_structure_or_process_1_subterm_name,
                 quality_id, quality_label, affected_structure_or_process_2_superterm_id,
                 affected_structure_or_process_2_superterm_name, affected_structure_or_process_2_subterm_id,
                 affected_structure_or_process_2_subterm_name, environment_id, environment_label, evidence_code_id,
                 evidence_code_symbol, evidence_code_label, publication_id, publication_label, publication_url,
                 taxon_id, taxon_label, v_uid, v_uuid, v_lastmodified) = row

                if extrinsic_genotype_id != '' or extrinsic_genotype_id is not None:
                    print('Skipping genotype with extrinsic modifiers: '+effective_genotype_id)
                    continue
                elif extrinsic_genotype_id == '' or extrinsic_genotype_id is None:
                    #print(phenotype_id)
                    #genes = implicated_gene_ids.split()
                    #print(genes)
                    if effective_genotype_id not in zfin_genotype_to_phenotype_hash:
                        zfin_genotype_to_phenotype_hash[effective_genotype_id] = [phenotype_id]
                        #print(zfin_genotype_to_phenotype_hash[genotype_id])
                    else:
                        zfin_genotype_to_phenotype_hash[effective_genotype_id].append(phenotype_id)
                        #print(zfin_genotype_to_phenotype_hash[genotype_id])
                        #print(len(zfin_genotype_to_phenotype_hash.keys()))
                        print('Repeat genotype: '+effective_genotype_id)
                    if limit is not None and line_counter > limit:
                        break
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish phenotype to gene data.')
        print('INFO: '+str(len(zfin_genotype_to_phenotype_hash.keys()))+' zebrafish phenotypes present.')
        return

    def assemble_nif_mgi_genotype_to_phenotype(self, limit=None):
        #TODO: Assuming want to filter out to intrinsic genotypes only?
        # Can filter on extrinsic genotype = ''
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        inter = 'inter/mgi/mouse_genotype_phenotype_hash.txt'
        mgi_genotype_to_phenotype_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' mouse genotype-phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (annotation_id, effective_genotype_id, effective_genotype_label, effective_genotype_label_html,
                 intrinsic_genotype_id, intrinsic_genotype_label, intrinsic_genotype_label_html,
                 genomic_variation_complement_id, genomic_variation_complement_label,
                 genomic_variation_complement_label_html, implicated_gene_ids, implicated_gene_labels,
                 implicated_sequence_alteration_ids, implicated_sequence_alteration_labels, genomic_background_id,
                 genomic_background_label, phenotype_id, phenotype_label, phenotype_description_free_text,
                 phenotype_modifier, evidence_code_id, evidence_code_symbol, evidence_code_label, environment_id,
                 environment_label, publication_id, publication_label, publication_url, taxon_id,
                 taxon_label, e_uid, v_uid, v_uuid, v_lastmodified) = row

                #print(phenotype_id)
                #genes = implicated_gene_ids.split()
                #print(genes)
                if effective_genotype_id not in mgi_genotype_to_phenotype_hash:
                    mgi_genotype_to_phenotype_hash[effective_genotype_id] = [phenotype_id]
                    #print(mgi_genotype_to_phenotype_hash[genotype_id])
                else:
                    mgi_genotype_to_phenotype_hash[effective_genotype_id].append(phenotype_id)
                    #print(mgi_genotype_to_phenotype_hash[genotype_id])
                    #print(len(mgi_genotype_to_phenotype_hash.keys()))
                    print('Repeat genotype: '+effective_genotype_id)
                if limit is not None and line_counter > limit:
                    break
        with open(inter, 'wb') as handle:
            pickle.dump(mgi_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling mouse genotype to phenotype data.')
        print('INFO: '+str(len(mgi_genotype_to_phenotype_hash.keys()))+' mouse phenotypes present.')
        return

    def assemble_nif_hpo_disease_to_gene(self, limit=None):
        print('INFO:Assembling human disease to gene data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/dvp.pr_nlx_151835_3'
        inter = 'inter/hpo/human_disease_gene_hash.txt'
        hpo_disease_to_gene_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human disease to gene rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, disease_id, disorder_name, disorder_database_link, gene_id,
                 gene_num, gene_label, v_uid, v_uuid, v_lastmodified) = row
                print(disease_id)

                #print(genes)
                if disease_id not in hpo_disease_to_gene_hash:
                    hpo_disease_to_gene_hash[disease_id] = [gene_id]
                    #print(hpo_phenotype_to_gene_hash[genotype_id])
                else:
                    hpo_disease_to_gene_hash[disease_id].append(gene_id)
                    #print(hpo_disease_to_gene_hash[disease_id])
                    #print(len(hpo_disease_to_gene_hash.keys()))
                    print('Repeat disease: '+disease_id)
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(hpo_disease_to_gene_hash, handle)
        print('INFO: Done assembling human phenotype to gene data.')
        print('INFO: '+str(len(hpo_disease_to_gene_hash.keys()))+' human phenotypes present.')
        return

    def assemble_nif_hpo_disease_to_phenotype(self, limit=None):
        print('INFO:Assembling human disease to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/dvp.pr_nlx_151835_1'
        inter = 'inter/hpo/nif_human_disease_phenotype_hash.txt'
        hpo_disease_to_phenotype_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human disease to phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, disorder_id, disorder_database_prefix, disorder_id_num, disorder_database_link, disorder_name,
                 disorder_qualifier, phenotype_id, phenotype_label, publication_id, evidence_code_id,
                 evidence_code_symbol, evidence_code_label, onset_id, onset_label, frequency, aspect, aspect_text,
                 synonyms, v_uid, v_uuid, v_lastmodified) = row

                print(disorder_id)

                #print(genes)
                if disorder_id not in hpo_disease_to_phenotype_hash:
                    hpo_disease_to_phenotype_hash[disorder_id] = [phenotype_id]
                    #print(hpo_disease_to_phenotype_hash[disease)id])
                else:
                    hpo_disease_to_phenotype_hash[disorder_id].append(phenotype_id)
                    #print(hpo_disease_to_phenotype_hash[disorder_id])
                    #print(len(hpo_disease_to_phenotype_hash.keys()))
                    print('Repeat disease: '+disorder_id)
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(hpo_disease_to_phenotype_hash, handle)
        print('INFO: Done assembling human disease to phenotype data.')
        print('INFO: '+str(len(hpo_disease_to_phenotype_hash.keys()))+' human diseases present.')
        return

    def assemble_nif_mgi_gene_to_phenotype(self, limit=None):
        print('INFO:Assembling mouse gene to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        inter = 'inter/mgi/mouse_gene_phenotype_hash.txt'
        mgi_gene_to_phenotype_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' mouse gene to phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (annotation_id, effective_genotype_id, effective_genotype_label, effective_genotype_label_html,
                 intrinsic_genotype_id, intrinsic_genotype_label, intrinsic_genotype_label_html,
                 genomic_variation_complement_id, genomic_variation_complement_label,
                 genomic_variation_complement_label_html, implicated_gene_ids, implicated_gene_labels,
                 implicated_sequence_alteration_ids, implicated_sequence_alteration_labels, genomic_background_id,
                 genomic_background_label, phenotype_id, phenotype_label, phenotype_description_free_text,
                 phenotype_modifier, evidence_code_id, evidence_code_symbol, evidence_code_label, environment_id,
                 environment_label, publication_id, publication_label, publication_url, taxon_id,
                 taxon_label, e_uid, v_uid, v_uuid, v_lastmodified) = row

                #print(implicated_gene_ids)

                if not re.match('.*,.*',implicated_gene_ids):
                    print(implicated_gene_labels)
                    #print(genes)
                    if implicated_gene_ids not in mgi_gene_to_phenotype_hash:
                        mgi_gene_to_phenotype_hash[implicated_gene_ids] = [phenotype_id]
                    #print(mgi_gene_to_phenotype_hash[gene_id])
                    else:
                        mgi_gene_to_phenotype_hash[implicated_gene_ids].append(phenotype_id)
                        #print(mgi_gene_to_phenotype_hash[gene_id])
                        #print(len(mgi_gene_to_phenotype_hash.keys()))
                        print('Repeat gene: '+implicated_gene_ids)
                else:
                    print('Skipping multi-gene genotype: '+effective_genotype_label)
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(mgi_gene_to_phenotype_hash, handle)
        print('INFO: Done assembling mouse gene to phenotype data.')
        print('INFO: '+str(len(mgi_gene_to_phenotype_hash.keys()))+' human genes present.')
        return

    def assemble_nif_zfin_gene_to_phenotype(self, limit=None):
        print('INFO:Assembling zebrafish gene to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        inter = 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt'
        zfin_gene_to_phenotype_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish gene to phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, effective_genotype_id, effective_genotype_label, effective_genotype_label_html,
                 intrinsic_genotype_id, intrinsic_genotype_num, intrinsic_genotype_label, intrinsic_genotype_label_html,
                 extrinsic_genotype_id, extrinsic_genotype_label, extrinsic_genotype_label_html, phenotype_id,
                 phenotype_label, phenotype_modifier, implicated_gene_ids, implicated_gene_labels, start_stage_id,
                 start_stage_zfin_id, start_stage_label, end_stage_id ,end_stage_zfin_id, end_stage_label, stages,
                 genomic_background_id, genomic_background_num, genomic_background_label,
                 affected_structure_or_process_1_superterm_id, affected_structure_or_process_1_superterm_name,
                 affected_structure_or_process_1_subterm_id, affected_structure_or_process_1_subterm_name,
                 quality_id, quality_label, affected_structure_or_process_2_superterm_id,
                 affected_structure_or_process_2_superterm_name, affected_structure_or_process_2_subterm_id,
                 affected_structure_or_process_2_subterm_name, environment_id, environment_label, evidence_code_id,
                 evidence_code_symbol, evidence_code_label, publication_id, publication_label, publication_url,
                 taxon_id, taxon_label, v_uid, v_uuid, v_lastmodified) = row

                if phenotype_id == '' or phenotype_id == None:
                    continue
                else:

                    #print(implicated_gene_ids)
                    #FIXME: Going to need to convert the ZFIN Gene IDs to NCBIGene IDs.

                    if not re.match('.*,.*',implicated_gene_ids):
                        print(implicated_gene_labels)
                        #print(genes)
                        if implicated_gene_ids not in zfin_gene_to_phenotype_hash:
                            zfin_gene_to_phenotype_hash[implicated_gene_ids] = [phenotype_id]
                        #print(zfin_gene_to_phenotype_hash[gene_id])
                        else:
                            zfin_gene_to_phenotype_hash[implicated_gene_ids].append(phenotype_id)
                            #print(zfin_gene_to_phenotype_hash[gene_id])
                            #print(len(zfin_gene_to_phenotype_hash.keys()))
                            print('Repeat gene: '+implicated_gene_ids)
                    else:
                        print('Skipping multi-gene genotype: '+effective_genotype_label)
                    if limit is not None and line_counter > limit:
                        break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_gene_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish gene to phenotype data.')
        print('INFO: '+str(len(zfin_gene_to_phenotype_hash.keys()))+' human phenotypes present.')
        return

    ####### OWLSIM DATA PROCESSING #######

    def perform_owlsim_queries(self, raw1, raw2, out, limit=None):
        print('INFO: Performing OWLSim queries.')
        line_counter = 0
        comparison_count = 0
        failure_counter = 0
        if limit is not None:
            print('Only querying first '+str(limit)+' phenotypic profile pairs.')
            comparison_count = limit
        #raw1 = 'inter/hpo/nif_human_disease_phenotype_hash.txt'
        #raw2 = 'inter/mgi/mouse_genotype_phenotype_hash.txt'
        data1 = open(raw1, 'rb')
        organism_a_hash = pickle.load(data1)
        data1.close()
        data2 = open(raw2, 'rb')
        organism_b_hash = pickle.load(data2)
        data2.close()
        if limit is None:
            comparison_count = len(organism_a_hash) * len(organism_b_hash)
            print('INFO: '+str(comparison_count)+' phenotypic profile comparisons to process.')
        #data2 = open(raw2,'r', encoding="iso-8859-1")
        #with open(raw1, 'r', encoding="iso-8859-1") as handle1:
            #organism_a_hash = pickle.loads(handle1.read())
        #with open(raw2, 'r', encoding="iso-8859-1") as handle2:
            #organism_b_hash = pickle.loads(handle2.read())
        #print(organism_a_hash)
        base_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?'
        #print(organism_a_hash)
        #print(organism_b_hash)
        with open(out, 'w', newline='') as outfile:
            #wlsimwriter = csv.writer(csvfile, delimiter='\t', quotechar="'")
            for i in organism_a_hash:
                entity_a = i
                entity_a_attributes = organism_a_hash[i]
                #print(attributes)
                #print(entity_a_attributes)
                phenotypic_profile_a = 'a='+('&a=').join(entity_a_attributes)
                for j in organism_b_hash:
                    entity_b = j
                    entity_b_attributes = organism_b_hash[j]
                    #print(entity_b_attributes)
                    phenotypic_profile_b = '&b='+('&b=').join(entity_b_attributes)
                    query_url = base_url+phenotypic_profile_a+phenotypic_profile_b
                    #print(query_url)
                    line_counter += 1
                    print('INFO: Processing phenotypic profile comparison '+str(line_counter)+' out of '+str(comparison_count)+'.')
                    try:
                        response = urllib.request.urlopen(query_url, timeout=5)
                        reader = codecs.getreader("utf-8")
                        data = json.load(reader(response))
                        #print(data)
                        #print('#####')
                        #print('query successful')
                        results = data['results']
                        maxIC = data['results'][0]['maxIC']
                        simJ = data['results'][0]['simJ']
                        #FIXME: Is this the correct variable to grab for the ICCS?
                        ICCS = data['results'][0]['bmaSymIC']
                        #FIXME: Is this the correct variable to grab for the simIC?
                        simIC = data['results'][0]['simGIC']
                        #print(results)
                        #FIXME: Queries are working, need to adjust writing output to file.

                        query_flag = 'success'
                        sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                        json.dump(sequence, outfile)
                        outfile.write('\n')

                        #print(sequence)
                        #print('failed here')
                        #row = str.join(sequence)

                        #print(row)
                        #owlsimwriter.writerow(row)
                        #print('query processing completed')

                    except Exception:
                        #print('Processing of OWLSim query failed.')
                        #Creating an empty set of metrics for failed queries (queries with unresolved IRIs).
                        #FIXME: May want to run a set with this and without this, as the 0s will effect averages.
                        maxIC = 0
                        simJ = 0
                        ICCS = 0
                        simIC = 0
                        query_flag = 'fail'

                        sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                        json.dump(sequence, outfile)
                        outfile.write('\n')
                        continue

        #entity_a = 'entity_1'
        #entity_a_attributes = ['attr1','attr2','attr3']
        #entity_b = 'genotype_id'
        #entity_b_attributes = ['attrb1','attrb2','attrb3']
        #print(str(row_count)+' human diseases to process.')


        #FIXME: Need to adjust attribute handling for first attribute and all following attributes (a= vs &a=)

        #query_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878'
        #phenotypic_profile_a = 'a='+('&a=').join(entity_a_attributes)
        #phenotypic_profile_b = '&b='+('&b=').join(entity_b_attributes)
        #combined_url = base_url+phenotypic_profile_a+phenotypic_profile_b
        #print(combined_url)
        #query_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0003731&b=HP:0000580'
        #query_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0003731&a=MP:0000599&a=MP:0005331&b=HP:0000580&b=HP:0002240&b=HP:0000831'


        return

    ####### PHENOLOG DATA PROCESSING #######

    def trim_panther_data(self, inter, taxons):
        print('INFO: Trimming PANTHER data.')
        line_counter = 0
        output_line_counter = 0
        raw = 'raw/panther/dvp.pr_nlx_84521_1'
        #inter = 'inter/panther/panther.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' PANTHER rows to process.')
        with open(inter, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            with open(raw, 'r', encoding="iso-8859-1") as csvfile:
                filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
                next(filereader,None)
                for row in filereader:
                    line_counter += 1
                    (panther_speciesa, tax_id_a, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,
                    proteina, panther_speciesb, tax_id_b, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,
                    gene_label_b, proteinb, orthology_class, orthology_class_label, ancestor_taxon, panther_id,
                    e_uid, v_uid, v_uuid, v_lastmodified) = row

                    #Currently filtering on the big three taxons, and ortholog relations only.
                    if (taxon_id_a in taxons or taxon_id_b in taxons) and orthology_class_label == 'Ortholog':
                        output_row = (panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,
                        proteina, panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,
                        gene_label_b, proteinb, orthology_class, orthology_class_label, panther_id)
                        #print('found one')
                        output_line_counter += 1
                        csvwriter.writerow(output_row)


        print('PANTHER file trimmed to '+str(output_line_counter)+' rows.')

        return

    def get_common_orthologs(self, inter, taxons):
        print('INFO: Getting common orthologs between species.')
        line_counter = 0
        ortholog_counter = 0
        raw = 'raw/panther/dvp.pr_nlx_84521_1'
        common_orthologs = []
        #inter = 'inter/panther/common_orthologs_'+species_a+'_vs_'+species_b+'.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' PANTHER rows to process.')

        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (panther_speciesa, tax_id_a, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,
                proteina, panther_speciesb, tax_id_b, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,
                gene_label_b, proteinb, orthology_class, orthology_class_label, ancestor_taxon, panther_id,
                e_uid, v_uid, v_uuid, v_lastmodified) = row

                if (taxon_id_a in taxons or taxon_id_b in taxons) and orthology_class_label == 'Ortholog':
                    if panther_id not in common_orthologs:
                        common_orthologs.append(panther_id)
                        ortholog_counter += 1

        with open(inter, 'wb') as handle:
            pickle.dump(common_orthologs, handle)

        print(str(ortholog_counter)+' common orthologs found.')

        return

    def get_ortholog(self, query_gene_id, panther):

        with open(panther, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            #print(str(row_count)+' zebrafish gene to phenotype rows to process.')
            print('query= '+query_gene_id)
            for row in filereader:
                #print(row)
                (panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,proteina,
                 panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,gene_label_b, proteinb,
                 orthology_class, orthology_class_label, panther_id) = row
                #print(panther_speciesa)
                #print(panther_id)

                if query_gene_id in [genea, gene_id_a, geneb, gene_id_b]:
                    result_panther_id = panther_id
                    print('found ortholog for '+query_gene_id+'.')
                    return(result_panther_id)

        print('no ortholog found for '+query_gene_id+'.')
        return('fail')

    def perform_phenolog_calculations(self, inter1, inter2, out, shared_orthologs, fdr_cutoff):
        print('INFO: Performing phenolog calculations.')
        line_counter = 0
        failure_counter = 0
        total_ortholog_matches = 0
        total_ortholog_nonmatches = 0
        ortholog_matches = 0
        ortholog_non_matches = 0
        phenotype_a_ortholog_count = 0
        phenotype_b_ortholog_count = 0

        with open(inter1, 'rb') as handle:
            species_a_pheno_gene_hash = pickle.load(handle)
        #print(species_a_pheno_gene_hash)
        with open(inter2, 'rb') as handle:
            species_b_pheno_gene_hash = pickle.load(handle)
        #print(species_b_pheno_gene_hash)
        with open(shared_orthologs, 'rb') as handle:
            num_shared_orthologs = len(pickle.load(handle))
        with open(out, 'w', newline='') as outfile:
            for i in species_a_pheno_gene_hash:
                species_a_phenotype_id = i
                species_a_orthologs = species_a_pheno_gene_hash[i]
                #print(species_a_orthologs)
                phenotype_a_ortholog_count = len(species_a_orthologs)

                for j in species_b_pheno_gene_hash:
                    species_b_phenotype_id = j
                    species_b_orthologs = species_b_pheno_gene_hash[j]
                    #print(species_b_orthologs)
                    ortholog_matches = 0
                    ortholog_non_matches = 0

                    phenotype_b_ortholog_count = len(species_b_orthologs)
                    for k in species_a_orthologs:

                        species_a_ortholog = k
                        for l in species_b_orthologs:
                            species_b_ortholog = l

                            if species_a_ortholog == species_b_ortholog:
                                #print('species a ortholog:'+species_a_ortholog+' matches species b ortholog:'+species_b_ortholog)
                                ortholog_matches += 1
                                total_ortholog_matches += 1
                            else:
                                #print('species a ortholog:'+species_a_ortholog+' does not match species b ortholog:'+species_b_ortholog)
                                ortholog_non_matches += 1
                                total_ortholog_nonmatches += 1

                        if ortholog_matches > 0:
                            #print('Matches: '+str(ortholog_matches))
                            #print('Non-matches: '+str(ortholog_non_matches))
                            m = float(phenotype_b_ortholog_count)
                            n = float(phenotype_a_ortholog_count)
                            N = float(num_shared_orthologs)
                            c = float(ortholog_matches)
                            prb = float(hypergeom.pmf(c, N, m, n))
                            if prb <= fdr_cutoff:
                                significance = 'Significant'
                            else:
                                significance = 'Not Significant'
                            #print(prb)
                            # Required output : phenotype a/b, species a/b, gene list a/b, probability, fdr adjusted probability?
                            sequence = (species_a_phenotype_id, species_a_orthologs, phenotype_a_ortholog_count, species_b_phenotype_id, species_b_orthologs, phenotype_b_ortholog_count, num_shared_orthologs, ortholog_matches, prb, fdr_cutoff, significance)
                            json.dump(sequence, outfile)
                            outfile.write('\n')

            print('Total Matches: '+str(total_ortholog_matches))
            print('Total non-matches: '+str(total_ortholog_nonmatches))

            #prb = "{:.2E}".format(Decimal(hypergeom.cdf(24, 5000, 47, 174)))
            #print(prb)
            #prb = 1 - hypergeom.cdf(0, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(1, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(2, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(3, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(4, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(5, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(6, 5000, 7, 12)
            #print(prb)
            #prb = 1 - hypergeom.cdf(7, 5000, 7, 12)
            #print(prb)
                # After the number of matching orthologs has been tallied, perform the
                # hypergeometric probability calculation for the phenotype-gene data objects,
                # then write the results to the output file.

                # Essentially, in comparing the ortholog lists we are seeing how many matches we get for a given set of
                # orthologs from phenotype B, with a certain number of draws. So, does this calculation need to be run in both directions,
                # given that the ortholog lists for species a and species b may be of different sizes?
                #TODO: Consult the phenolog paper on the question above.
                # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom

                # N = total number of orthologs shared between species
                # n = nummber of orthologs in species A phenotype
                # m = nummber of orthologs in species B phenotype
                # c = number of common orthologs between phenotypes (ortholog matches)

                # c =
                # M = Total number of objects. (Global number of orthologs for species A? All possible phenologs that could be drawn?)
                # n = Total number of type I objects. (Total number of orthologs in the list for phenotype B?)
                # N = Number of type I objects drawn. (Number of matching orhtologs.)
                #prb = hypergeom.cdf(x, M, n, N)
        return

    def generate_random_data(self, pheno_gene_hash, common_orthologs, out_dir, limit=1000):
        print('INFO: Creating random data sets.')
        #Take the phenotype-ortholog hashes I have created.
        #Remove all orthologs and place in a list, replace with 0s or 1s or something.
        #Randomly shuffle the ortholog list
        #Iterate through the hash and replace the 0s with an ortholog ID IF that ortholog ID is not present in the hash.

        #Question: How to be sure that the data set is random, and that I’m not creating 1000 identical data sets?
        # Should the random data set have the same presence of orthologs as the test set?
        # Meaning, if ortholog X is present in the test data set Y times, should it also be present in the test data set Y times?
        #Need to have a way to make the data set creation fail if we get to the end and can only put an ortholog with a phenotype that already has that ortholog with it.

        # Question: Is it appropriate to use the orthologs that are in the data set/associated with phenotypes, or
        # would it make more sense to populate from the total orthologs shared between two species?
        # What about getting random.sample from the total list of orthologs, with replacement?
        counter = 1
        while counter <= limit:
            #print('INFO: Creating random data set '+str(counter)+' out of '+str(limit)+'.')
            test_pheno_ortholog_hash = {}

            orthologs = []
            list_length = 0
            with open(pheno_gene_hash, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            with open(common_orthologs, 'rb') as handle:
                orthologs = pickle.load(handle)
            # Adjusting this to instead call from the full pool of common orthologs.
            #for i in pheno_ortholog_hash:
                #for j in pheno_ortholog_hash[i]:
                    #orthologs.append(j)
            #FIXME: How random is this? Is it sufficiently random?
            random.shuffle(orthologs)

            for i in pheno_ortholog_hash:
                test_pheno_ortholog_hash[i] = []
                ortholog_list_length = len(pheno_ortholog_hash[i])
                #print(ortholog_list_length)
                ortholog_draw = orthologs
                random.shuffle(ortholog_draw)
                #random_orthologs = random.sample(orthologs)
                #test_pheno_ortholog_hash[i].append(random_orthologs)

                for j in pheno_ortholog_hash[i]:
                    test_pheno_ortholog_hash[i].append(ortholog_draw.pop())
                    random.shuffle(ortholog_draw)
                    #list_length = len(orthologs)
                    #print(list_length)

                    #if orthologs[(1 - list_length)] not in test_pheno_ortholog_hash[i]:
                        #test_pheno_ortholog_hash[i].append(orthologs.pop())
                    #else:
                        #while orthologs[(1 - list_length)] in test_pheno_ortholog_hash[i]:
                            #random.shuffle(orthologs)
                        #test_pheno_ortholog_hash[i].append(orthologs.pop())
                    #if orthologs.pop
            #print('Completed randomization successfully!')
            with open((out_dir+'random_'+str(counter)+'.txt'), 'wb') as handle:
                pickle.dump(test_pheno_ortholog_hash, handle)
            counter += 1
        return

    def set_stage_for_FDR_calculation(self, limit=1000):
        print('INFO: Setting stage for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        main.generate_random_data('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_mouse.txt', 'inter/random/human_vs_mouse/mouse/')
        main.generate_random_data('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_mouse.txt', 'inter/random/human_vs_mouse/human/')

        main.generate_random_data('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', 'inter/random/human_vs_zebrafish/zebrafish/')
        main.generate_random_data('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', 'inter/random/human_vs_zebrafish/human/')

        main.generate_random_data('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', 'inter/random/mouse_vs_zebrafish/zebrafish/')
        main.generate_random_data('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', 'inter/random/mouse_vs_zebrafish/mouse/')
        print('INFO: Done with random data generation.')
        fdr_global_p_value_list = []

        hvm_human_dir = 'inter/random/human_vs_mouse/human/'
        hvm_mouse_dir = 'inter/random/human_vs_mouse/mouse/'
        hvz_human_dir = 'inter/random/human_vs_zebrafish/human/'
        hvz_zebrafish_dir = 'inter/random/human_vs_zebrafish/zebrafish/'
        mvz_mouse_dir = 'inter/random/mouse_vs_zebrafish/mouse/'
        mvz_zebrafish_dir = 'inter/random/mouse_vs_zebrafish/zebrafish/'

        with open('inter/panther/common_orthologs_human_mouse.txt', 'rb') as handle:
            hvm_common_orthologs = len(pickle.load(handle))
        with open('inter/panther/common_orthologs_human_zebrafish.txt', 'rb') as handle:
            hvz_common_orthologs = len(pickle.load(handle))
        with open('inter/panther/common_orthologs_mouse_zebrafish.txt', 'rb') as handle:
            mvz_common_orthologs = len(pickle.load(handle))
        #print(hvm_common_orthologs)
        #print(hvz_common_orthologs)
        #print(mvz_common_orthologs)

        random_counter = 1
        # Switch to 'while' when ready for full set testing.
        while random_counter < limit:
            random_counter += 1
            print('INFO: Performing phenolog calculation on random data set '+str(random_counter)+' out of '+str(limit)+'.')
            fdr_p_value_list = []
            hvm_human_file = hvm_human_dir+'random_'+str(random_counter)+'.txt'
            hvm_mouse_file = hvm_mouse_dir+'random_'+str(random_counter)+'.txt'
            hvz_human_file = hvz_human_dir+'random_'+str(random_counter)+'.txt'
            hvz_zebrafish_file = hvz_zebrafish_dir+'random_'+str(random_counter)+'.txt'
            mvz_mouse_file = mvz_mouse_dir+'random_'+str(random_counter)+'.txt'
            mvz_zebrafish_file = mvz_zebrafish_dir+'random_'+str(random_counter)+'.txt'

            with open(hvm_human_file, 'rb') as handle:
                hvm_human_pheno_ortholog_hash = pickle.load(handle)
            with open(hvm_mouse_file, 'rb') as handle:
                hvm_mouse_pheno_ortholog_hash = pickle.load(handle)
            with open(hvz_human_file, 'rb') as handle:
                hvz_human_pheno_ortholog_hash = pickle.load(handle)
            with open(hvz_zebrafish_file, 'rb') as handle:
                hvz_zebrafish_pheno_ortholog_hash = pickle.load(handle)
            with open(mvz_mouse_file, 'rb') as handle:
                mvz_mouse_pheno_ortholog_hash = pickle.load(handle)
            with open(mvz_zebrafish_file, 'rb') as handle:
                mvz_zebrafish_pheno_ortholog_hash = pickle.load(handle)

            #TODO: Need to return all of the p-values from the hypergeometric probability calculation for sorting and 5% cutoff.
            hvm_phenolog_p_values = main.perform_phenolog_calculations_for_FDR(hvm_human_pheno_ortholog_hash, hvm_mouse_pheno_ortholog_hash, hvm_common_orthologs)
            #print(hvm_phenolog_p_values)
            fdr_p_value_list.extend(hvm_phenolog_p_values)
            hvz_phenolog_p_values = main.perform_phenolog_calculations_for_FDR(hvz_human_pheno_ortholog_hash, hvz_zebrafish_pheno_ortholog_hash, hvz_common_orthologs)
            #print(hvz_phenolog_p_values)
            fdr_p_value_list.extend(hvz_phenolog_p_values)
            mvz_phenolog_p_values = main.perform_phenolog_calculations_for_FDR(mvz_mouse_pheno_ortholog_hash, mvz_zebrafish_pheno_ortholog_hash, mvz_common_orthologs)
            #print(mvz_phenolog_p_values)
            fdr_p_value_list.extend(mvz_phenolog_p_values)

            # After grabbing the p-values from each function, assemble and sort.
            # Select the p-value that resides at the 0.05 percentile and add it to a list.


            #print('fdr p value list: '+str(len(fdr_p_value_list)))
            fdr_p_value_list.sort()
            #print(fdr_p_value_list)
            cutoff_position = math.ceil((len(fdr_p_value_list))*0.05) - 1
            #print(fdr_p_value_list[cutoff_position])
            fdr_global_p_value_list.append(fdr_p_value_list[cutoff_position])

        fdr_global_p_value_list.sort()
        global_cutoff_position = math.ceil((len(fdr_global_p_value_list))*0.05) - 1
        print('The empirical FDR adjustment cutoff is '+str(fdr_global_p_value_list[global_cutoff_position])+'.')

        return fdr_global_p_value_list[global_cutoff_position]


    def perform_phenolog_calculations_for_FDR(self, species_a_po_hash, species_b_po_hash, shared_orthologs):
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        line_counter = 0
        failure_counter = 0
        total_ortholog_matches = 0
        total_ortholog_nonmatches = 0
        ortholog_matches = 0
        ortholog_non_matches = 0
        phenotype_a_ortholog_count = 0
        phenotype_b_ortholog_count = 0
        total_hyp_calcs = 0
        phenolog_p_value_list = []

        species_a_pheno_gene_hash = species_a_po_hash
        species_b_pheno_gene_hash = species_b_po_hash

        for i in species_a_pheno_gene_hash:
            # Phenotype for species A
            species_a_phenotype_id = i
            species_a_orthologs = species_a_pheno_gene_hash[i]
            #print(species_a_orthologs)
            phenotype_a_ortholog_count = len(species_a_orthologs)

            for j in species_b_pheno_gene_hash:
                # Phenotype for species B
                species_b_phenotype_id = j
                species_b_orthologs = species_b_pheno_gene_hash[j]
                #print(species_b_orthologs)
                #ortholog_matches = 0
                #ortholog_non_matches = 0

                ortholog_matches = 0
                ortholog_non_matches = 0

                phenotype_b_ortholog_count = len(species_b_orthologs)
                for k in species_a_orthologs:
                    # Orthologs for species A
                    #ortholog_matches = 0
                    #ortholog_non_matches = 0

                    species_a_ortholog = k
                    for l in species_b_orthologs:
                        # Orthologs for species B
                        species_b_ortholog = l
                        if species_a_ortholog == species_b_ortholog:
                            #print('species a ortholog:'+species_a_ortholog+' matches species b ortholog:'+species_b_ortholog)
                            ortholog_matches += 1
                            #print(species_a_ortholog+' == '+species_b_ortholog)
                            total_ortholog_matches += 1
                        else:
                            #print('species a ortholog:'+species_a_ortholog+' does not match species b ortholog:'+species_b_ortholog)
                            ortholog_non_matches += 1
                            total_ortholog_nonmatches += 1

                if ortholog_matches > 0:
                    #print('Matches: '+str(ortholog_matches))
                    #print('Non-matches: '+str(ortholog_non_matches))
                    m = float(phenotype_b_ortholog_count)
                    n = float(phenotype_a_ortholog_count)
                    N = float(shared_orthologs)
                    c = float(ortholog_matches)
                    prb = float(hypergeom.pmf(c, N, m, n))
                    #print(str(c)+', '+str(N)+', '+str(m)+', '+str(n))
                    #print(prb)
                    phenolog_p_value_list.append(prb)
                    total_hyp_calcs += 1
        print('Total Matches: '+str(total_ortholog_matches))
        print('Total non-matches: '+str(total_ortholog_nonmatches))
        print('Total phenolog calculations: '+str(total_hyp_calcs))

            # After the number of matching orthologs has been tallied, perform the
            # hypergeometric probability calculation for the phenotype-gene data objects,
            # then write the results to the output file.

            # Essentially, in comparing the ortholog lists we are seeing how many matches we get for a given set of
            # orthologs from phenotype B, with a certain number of draws. So, does this calculation need to be run in both directions,
            # given that the ortholog lists for species a and species b may be of different sizes?
            #TODO: Consult the phenolog paper on the question above.
            # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom

            # N = total number of orthologs shared between species
            # n = nummber of orthologs in species A phenotype
            # m = nummber of orthologs in species B phenotype
            # c = number of common orthologs between phenotypes (ortholog matches)

            # c =
            # M = Total number of objects. (Global number of orthologs for species A? All possible phenologs that could be drawn?)
            # n = Total number of type I objects. (Total number of orthologs in the list for phenotype B?)
            # N = Number of type I objects drawn. (Number of matching orhtologs.)
            #prb = hypergeom.cdf(x, M, n, N)

        return phenolog_p_value_list

###MAIN####

limit = 500
main = main()

#Trim the PANTHER data set for each taxon.
#main.trim_panther_data('inter/panther/panther_human.txt', ['NCBITaxon:9606'])
#main.trim_panther_data('inter/panther/panther_mouse.txt', ['NCBITaxon:10090'])
#main.trim_panther_data('inter/panther/panther_zebrafish.txt', ['NCBITaxon:7955'])
#main.trim_panther_data('inter/panther/panther_hmz_trio.txt', ['NCBITaxon:9606', 'NCBITaxon:7955', 'NCBITaxon:10090'])

#main.get_common_orthologs('inter/panther/common_orthologs_human_zebrafish.txt', ['NCBITaxon:9606', 'NCBITaxon:7955'])
#main.get_common_orthologs('inter/panther/common_orthologs_human_mouse.txt', ['NCBITaxon:9606', 'NCBITaxon:10090'])
#main.get_common_orthologs('inter/panther/common_orthologs_mouse_zebrafish.txt', ['NCBITaxon:10090', 'NCBITaxon:7955'])

### Data assembly via SciGraph ###
#main._assemble_human_disease_to_phenotype(limit)
#main._assemble_mouse_genotype_to_phenotype(limit)
#FIXME: Note that the zebrafish data is not currently available through REST services.
#main.assemble_zebrafish_genotype_to_phenotype(500)

### Data assembly via NIF/DISCO ###
#main.assemble_nif_zfin_phenotype_to_gene(limit)
#main.assemble_nif_mgi_phenotype_to_gene(limit)
#main.assemble_nif_hpo_phenotype_to_gene(limit)
#main.assemble_nif_animalqtl_phenotype_to_gene(limit)

#main.assemble_nif_hpo_disease_to_gene(limit)
#main.assemble_nif_zfin_genotype_to_phenotype(limit)
#main.assemble_nif_mgi_genotype_to_phenotype(limit)
#main.assemble_nif_mgi_gene_to_phenotype(limit)
#main.assemble_nif_zfin_gene_to_phenotype(limit)
#main.assemble_nif_hpo_disease_to_phenotype(limit)


####### OWLSIM COMPARISONS #######

# Compare human disease phenotypic profiles & mouse genotype phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human diseases vs mouse genotypes')
#main.perform_owlsim_queries('inter/hpo/nif_human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt','out/owlsim/human_disease_mouse_genotype.txt')
#print('INFO: Done processing human diseases vs mouse genotypes')

# Compare human disease phenotypic profiles & zebrafish genotype phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human disease vs zebrafish genotype')
#main.perform_owlsim_queries('inter/hpo/nif_human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt','out/owlsim/human_disease_zebrafish_genotype.txt')
#print('INFO: Done processing human disease vs zebrafish genotype')

# Compare mouse genotype phenotypic profiles & zebrafish genotype phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing mouse genotype vs zebrafish genotypes')
#main.perform_owlsim_queries('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt','out/owlsim/mouse_genotype_zebrafish_genotype.txt')
#print('INFO: Done processing mouse genotype vs zebrafish genotypes')

# Compare human disease phenotypic profiles & mouse gene phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human disease vs mouse genes')
#main.perform_owlsim_queries('inter/hpo/nif_human_disease_phenotype_hash.txt', 'inter/mgi/mouse_gene_phenotype_hash.txt','out/owlsim/human_disease_mouse_gene.txt')
#print('INFO: Done processing human disease vs mouse genes')

# Compare human disease phenotypic profiles & zebrafish gene phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human disease vs zebrafish genes')
#main.perform_owlsim_queries('inter/hpo/nif_human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt','out/owlsim/human_disease_zebrafish_gene.txt')
#print('INFO: Done processing human disease vs zebrafish genes')

# Compare mouse gene phenotypic profiles & zebrafish gene phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing mouse genes vs zebrafish genes')
#main.perform_owlsim_queries('inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt','out/owlsim/mouse_gene_zebrafish_gene.txt')
#print('INFO: Done processing mouse genes vs zebrafish genes')



#TODO: Data assembly
# Would it be easier/more efficient to add the phenolog data to the OWLSim output table,
# or do separate output tables and then match them together into a combined table? Essentially need to do an SQL join.


####### FDR CALCULATION #######
fdr_cutoff = main.set_stage_for_FDR_calculation()

####### PHENOLOG COMPARISONS #######
# NOTE: Either run the FDR calculations or set an FDR cutoff before running the phenolog calculations.
#fdr_cutoff = 0.000229
main.perform_phenolog_calculations('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt', 'out/phenolog/human_vs_mouse.txt', 'inter/panther/common_orthologs_human_mouse.txt', fdr_cutoff)
main.perform_phenolog_calculations('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'out/phenolog/mouse_vs_zebrafish.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', fdr_cutoff)
main.perform_phenolog_calculations('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'out/phenolog/human_vs_zebrafish.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', fdr_cutoff)



elapsed_time = time.time() - start_time
print('Processing completed in '+str(elapsed_time)+' seconds.')

#TODO: Make sure and have the ability to filter between single-gene genotypes and multi-gene genotypes.

#http://owlsim.monarchinitiative.org/compareAttributeSets?a=HP:0001263&b=MP:0010864
#Format for mutliple:
#http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878




