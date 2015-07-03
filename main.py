#!/usr/bin/env python3

#AUTHOR: Bryan Laraway
#PROJECT: Comparison of OWLSim and Phenologs for the identification of models of human disease and gene candidates for human disease.
#PURPOSE: This script will call of the functions/methods/scripts for performing the processing required for this analysis.

import json
import urllib.request
import codecs
import time
import gc
from socket import *
import os
import re
import csv
import pickle
from decimal import Decimal, getcontext
import cProfile
import numpy
from numpy import random
from scipy.stats import hypergeom, pearsonr
import math
import heapq
import multiprocessing
#from multiprocessing import Pool, Process, Manager, Lock, Value
import itertools
from threading import Thread
#import threading
from ctypes import c_int
from queue import Queue
#import matplotlib.pyplot as plt

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
    #NOTE: While this code works, acquiring data via SciGraph is rate-limiting due to the URL calls.
    #I have instead downloaded the flat files from NIF/DISCO

    def _assemble_human_disease_to_phenotype(self, limit=None):
        """This function takes a list of diseases from HPO and retrieves the associated phenotypes from SciGraph."""

        print('INFO: Assembling human disease to phenotype data.')

        # Set counters and open files for processing.
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/diseases.csv'
        inter = 'inter/hpo/human_disease_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = (sum(1 for row in filereader)) - 1
            print(str(row_count)+' human diseases to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:

                # Query SciGraph using the disease_id and retrieve the phenotypes in a JSON object.
                line_counter += 1
                disease_id = row[0]
                disease_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/diseases/'+disease_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(disease_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))

                    # Parse the JSON object and add the disease-phenotype associations to the hash.
                    pheno_ids = data['nodes']
                    for rs in pheno_ids:
                        if disease_id not in hu_disease_to_phenotype_hash:
                            hu_disease_to_phenotype_hash[disease_id] = [rs['id']]
                        else:
                            hu_disease_to_phenotype_hash[disease_id].append(rs['id'])
                except Exception:
                    print('Retrieval of '+disease_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break

        # Dump data to file.
        with open(inter, 'wb') as handle:
            pickle.dump(hu_disease_to_phenotype_hash, handle)
        print('INFO: Done assembling human disease to phenotype data.')
        print('INFO: '+str(len(hu_disease_to_phenotype_hash.keys()))+' human diseases processed.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    def _assemble_mouse_genotype_to_phenotype(self, limit=None):
        """This function takes a list of genotypes from MGI and retrieves the associated phenotypes from SciGraph."""

        print('INFO:Assembling mouse genotype to phenotype data.')

        # Set counters and open files for processing.
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/genotypes.csv'
        inter = 'inter/mgi/mouse_geno_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = (sum(1 for row in filereader)) - 1
            print(str(row_count)+' mouse genotypes to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:

                # Query SciGraph using the genotype_id and retrieve the phenotypes in a JSON object.
                line_counter += 1
                genotype_id = row[0]
                genotype_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/features/'+genotype_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(genotype_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))

                    # Parse the JSON object and add the genotype-phenotype associations to the hash.
                    pheno_ids = data['nodes']
                    for rs in pheno_ids:
                        if genotype_id not in mouse_genotype_to_phenotype_hash:
                            mouse_genotype_to_phenotype_hash[genotype_id] = [rs['id']]
                        else:
                            mouse_genotype_to_phenotype_hash[genotype_id].append(rs['id'])
                except Exception:
                    print('Retrieval of '+genotype_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break

        # Dump data to file.
        with open(inter, 'wb') as handle:
            pickle.dump(mouse_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling mouse genotype to phenotype data.')
        print('INFO: '+str(len(mouse_genotype_to_phenotype_hash.keys()))+' mouse genotypes present.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    def assemble_zebrafish_genotype_to_phenotype(self, limit=None):
        """This function takes a list of genotypes from ZFIN and retrieves the associated phenotypes from SciGraph."""

        print('INFO:Assembling zebrafish genotype to phenotype data.')

        # Set counters and open files for processing.
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/genotypes.csv'
        inter = 'inter/zfin/zebrafish_geno_pheno_hash.txt'
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = (sum(1 for row in filereader)) - 1
            print(str(row_count)+' zebrafish genotypes to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:

                # Query SciGraph using the genotype_id and retrieve the phenotypes in a JSON object.
                line_counter += 1
                genotype_id = row[0]
                genotype_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/features/'+genotype_id+'/phenotypes/targets'
                try:
                    response = urllib.request.urlopen(genotype_url, timeout=5)
                    reader = codecs.getreader("utf-8")
                    data = json.load(reader(response))

                    # Parse the JSON object and add the genotype-phenotype associations to the hash.
                    pheno_ids = data['nodes']
                    #print(pheno_ids)
                    for rs in pheno_ids:
                        if genotype_id not in zfin_genotype_to_phenotype_hash:
                            zfin_genotype_to_phenotype_hash[genotype_id] = [rs['id']]
                        else:
                            zfin_genotype_to_phenotype_hash[genotype_id].append(rs['id'])
                except Exception:
                    print('Retrieval of '+genotype_id+' failed.')
                    failure_counter += 1
                    continue
                if limit is not None and line_counter > limit:
                    break

        # Dump data to file.
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish genotype to phenotype data.')
        print('INFO: '+str(len(zfin_genotype_to_phenotype_hash))+' zebrafish genotypes present.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return

    ####### NIF DATA ASSEMBLY #######

    ####### PHENOLOG PHENOTYPE TO GENE #######

    def assemble_nif_zfin_phenotype_to_gene(self, limit=None):
        """This function assembles zebrafish phenotype to gene associations from the NIF/DISCO flat data file"""

        print('INFO:Assembling zebrafish phenotype to ortholog data.')

        # Set up counters and open required files.
        line_counter = 0
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

                # Read in a row and split into individual variables.
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

                # Skip phenotypes without IDs and phenotypes with no associated genes.
                if phenotype_id == '' or phenotype_id is None:
                    continue
                if implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue

                # Split the implicated genes list.
                genes = implicated_gene_ids.split(',')

                # If phenotype is not in the phenotype to gene hash, add phenotype to hash.
                if phenotype_id not in zfin_phenotype_to_gene_hash:
                    zfin_phenotype_to_gene_hash[phenotype_id] = []

                # If phenotype is not in the phenotype to ortholog hash, add phenotype to hash.
                if phenotype_id not in zfin_phenotype_to_ortholog_hash:
                    zfin_phenotype_to_ortholog_hash[phenotype_id] = []

                # If gene is not in the phenotype to gene hash, add gene to hash.
                for gene in genes:
                    if gene not in zfin_phenotype_to_gene_hash[phenotype_id]:
                        zfin_phenotype_to_gene_hash[phenotype_id].append(gene)

                        # Convert genes to orthologs using zebrafish-trimmed PANTHER table as lookup.
                        panther_id = self.get_ortholog(gene, 'inter/panther/panther_zebrafish.txt')

                        # If ortholog is not in the phenotype to ortholog hash, add ortholog to hash.
                        if panther_id != 'fail' and panther_id not in zfin_phenotype_to_ortholog_hash[phenotype_id]:
                            zfin_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break

        # Dump data to files.
        with open(inter1, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_ortholog_hash, handle)

        print('INFO: Done assembling zebrafish phenotype to gene/ortholog data.')
        print('INFO: '+str(len(zfin_phenotype_to_gene_hash.keys()))+' zebrafish phenotypes present.')
        return

    def assemble_nif_mgi_phenotype_to_gene(self, limit=None):
        """This function assembles mouse phenotype to gene associations from the NIF/DISCO flat data file"""

        print('INFO:Assembling mouse phenotype to gene/ortholog data.')

        # Set up counters and open required files.
        line_counter = 0
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

                # Read in a row and split into individual variables
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

                # Skip phenotypes without IDs and phenotypes with no associated genes.
                if phenotype_id == '' or phenotype_id is None:
                    continue
                if implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue

                # Split the implicated genes list.
                genes = implicated_gene_ids.split(',')

                # If phenotype is not in the phenotype to gene hash, add phenotype to hash.
                if phenotype_id not in mgi_phenotype_to_gene_hash:
                    mgi_phenotype_to_gene_hash[phenotype_id] = []

                # If phenotype is not in the phenotype to ortholog hash, add phenotype to hash.
                if phenotype_id not in mgi_phenotype_to_ortholog_hash:
                    mgi_phenotype_to_ortholog_hash[phenotype_id] = []

                # If gene is not in the phenotype to gene hash, add gene to hash.
                for gene in genes:
                    if gene not in mgi_phenotype_to_gene_hash[phenotype_id]:
                        mgi_phenotype_to_gene_hash[phenotype_id].append(gene)

                        # Convert genes to orthologs using mouse-trimmed PANTHER table as lookup.
                        panther_id = self.get_ortholog(gene,'inter/panther/panther_mouse.txt')

                        # If ortholog is not in the phenotype to ortholog hash, add ortholog to hash.
                        if panther_id != 'fail' and panther_id not in mgi_phenotype_to_ortholog_hash[phenotype_id]:
                            mgi_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break

        # Dump data to files.
        with open(inter1, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_ortholog_hash, handle)

        print('INFO: Done assembling mouse phenotype to gene/ortholog data.')
        print('INFO: '+str(len(mgi_phenotype_to_gene_hash.keys()))+' mouse phenotypes present.')

        return


    def assemble_nif_hpo_phenotype_to_gene(self, limit=None):
        """This function assembles human phenotype to gene associations from the NIF/DISCO flat data file"""

        print('INFO:Assembling human phenotype to gene data.')

        # Set up counters and open required files.
        line_counter = 0
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
            next(filereader, None)
            for row in filereader:

                # Read in a row and split into individual variables
                line_counter += 1
                (e_uid, phenotype_id, phenotype_label, gene_id, gene_num,
                 gene_label, v_uid, v_uuid, v_lastmodified) = row
                print('INFO: Processing human phenotype row '+str(line_counter)+' out of '+str(row_count)+'.')

                if phenotype_id == '' or phenotype_id is None:
                    continue
                if gene_id == '' or gene_id is None:
                    continue

                # Convert NCBIGene ID prefix.
                gene_id = re.sub('NCBI_gene:', 'NCBIGene:', gene_id)

                # If phenotype is not in the phenotype to gene hash, add phenotype to hash.
                if phenotype_id not in hpo_phenotype_to_gene_hash:
                    hpo_phenotype_to_gene_hash[phenotype_id] = []

                # If phenotype is not in the phenotype to ortholog hash, add phenotype to hash.
                if phenotype_id not in hpo_phenotype_to_ortholog_hash:
                    hpo_phenotype_to_ortholog_hash[phenotype_id] = []

                # If gene is not in the phenotype to gene hash, add gene to hash.
                if gene_id not in hpo_phenotype_to_gene_hash[phenotype_id]:
                    hpo_phenotype_to_gene_hash[phenotype_id].append(gene_id)

                    # Convert genes to orthologs using human-trimmed PANTHER table as lookup.
                    panther_id = self.get_ortholog(gene_id,'inter/panther/panther_human.txt')

                    # If ortholog is not in the phenotype to ortholog hash, add ortholog to hash.
                    if panther_id != 'fail' and panther_id not in hpo_phenotype_to_ortholog_hash[phenotype_id]:
                            hpo_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break

        # Dump files to disk.
        with open(inter1, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_ortholog_hash, handle)

        print('INFO: Done assembling human phenotype to gene/ortholog data.')
        print('INFO: '+str(len(hpo_phenotype_to_gene_hash.keys()))+' human phenotypes present.')
        return

    def assemble_nif_animalqtl_phenotype_to_gene(self, limit=None):
        """ Currently not used """
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
                    if effective_genotype_id not in zfin_genotype_to_phenotype_hash:
                        zfin_genotype_to_phenotype_hash[effective_genotype_id] = [phenotype_id]
                    else:
                        zfin_genotype_to_phenotype_hash[effective_genotype_id].append(phenotype_id)
                    if limit is not None and line_counter > limit:
                        break
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish genotype to phenotype data.')
        print('INFO: '+str(len(zfin_genotype_to_phenotype_hash))+' zebrafish genotypes present.')
        return

    def assemble_nif_mgi_genotype_to_phenotype(self, limit=None):
        #TODO: Assuming want to filter out to intrinsic genotypes only?
        # Can filter on extrinsic genotype = ''
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
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

                if effective_genotype_id not in mgi_genotype_to_phenotype_hash:
                    mgi_genotype_to_phenotype_hash[effective_genotype_id] = [phenotype_id]
                else:
                    mgi_genotype_to_phenotype_hash[effective_genotype_id].append(phenotype_id)
                if limit is not None and line_counter > limit:
                    break
        with open(inter, 'wb') as handle:
            pickle.dump(mgi_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling mouse genotype to phenotype data.')
        print('INFO: '+str(len(mgi_genotype_to_phenotype_hash))+' mouse genotypes present.')
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
        print('INFO: '+str(len(hpo_disease_to_gene_hash))+' human phenotypes present.')
        return

    def assemble_nif_hpo_disease_to_phenotype(self, limit=None):
        print('INFO:Assembling human disease to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/dvp.pr_nlx_151835_1'
        inter = 'inter/hpo/human_disease_phenotype_hash.txt'
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
                    if phenotype_id not in hpo_disease_to_phenotype_hash[disorder_id]:
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
        print('INFO: '+str(len(hpo_disease_to_phenotype_hash))+' human diseases present.')
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
                if implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue
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
                    if implicated_gene_ids == '' or implicated_gene_ids is None:
                        continue

                    elif not re.match('.*,.*',implicated_gene_ids):
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

    def perform_owlsim_queries_threaded(self, raw1, raw2, out, limit=None):
        print('INFO: Performing OWLSim queries.')
        line_counter = 0
        comparison_count = 0
        failure_counter = 0
        comparison_hash = {}
        comparison_list = []
        if limit is not None:
            print('Only querying first '+str(limit)+' phenotypic profile pairs.')
            comparison_count = limit
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
        #base_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?'
        base_url = 'http://0.0.0.0:9031/compareAttributeSets?'
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
                    #print('INFO: Assembling phenotypic profile comparison query '+str(line_counter)+' out of '+str(comparison_count)+'.')
                    comparison_id = entity_a+'_'+entity_b
                    comparison_list.append((comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes))
            print('INFO: Done assembling phenotypic profile comparison queries.')

            ###### THREADING INSERT ######
            queue = Queue(maxsize=50)
            num_threads = len(comparison_list)
            for sequence in comparison_list:
                queue.put(tuple)

                worker = Thread(target=multithread_owlsim_queries, args=(sequence))


            if __name__ == '__main__':

                print('INFO: Multithreading started')


                #multiprocessing.Semaphore(cores)
                #jobs = []
                #phenotype_iterable = []
                #phenotype_counter = 0

                #(comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes) = tuple
                #phenotype_counter += 1
                #print('Working on phenotype '+str(phenotype_counter)+' out of '+str(len(phenotype_list))+'.')
                print('Creating pool.')
                #results = [pool.apply_async(multiprocess_owlsim_queries, args=(tuple)) for tuple in comparison_list]
                print('Processing results.')
                comparison_list = []
                for p in results:
                    sequence  = p.get()
                    #sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                    json.dump(sequence, outfile)
                    outfile.write('\n')
                print('Done processing results.')

                print('INFO: Multiprocessing completed')
                ###### END THREADING INSERT ######

        return

    def assemble_owlsim_queries(self, raw1, raw2, interfile_directory, interfile_prefix, limit=None):

        line_counter = 0
        line_counter_max = 5000000
        file_counter = 1
        comparison_count = 0
        failure_counter = 0
        comparison_hash = {}
        comparison_list = []
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
        #base_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?'
        base_url = 'http://0.0.0.0:9031/compareAttributeSets?'
        #print(organism_a_hash)
        #print(organism_b_hash)
        print('INFO: Assembling phenotypic profile comparison queries.')
        file_name = interfile_directory+'/'+interfile_prefix+'_'+str(file_counter)+'.txt'
        csvfile = open(file_name, 'w', newline='')
        csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
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

                #print('INFO: Assembling phenotypic profile comparison query '+str(line_counter)+' out of '+str(comparison_count)+'.')
                comparison_id = entity_a+'_'+entity_b
                #comparison_tuple = ([comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes])
                output_row = (comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes)


                csvwriter.writerow(output_row)
                line_counter += 1
                if line_counter == line_counter_max:
                    line_counter = 0
                    file_counter += 1
                    csvfile.close()
                    file_name = interfile_directory+'/'+interfile_prefix+'_'+str(file_counter)+'.txt'
                    csvfile = open(file_name, 'w', newline='')
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                if limit is not None and line_counter > limit:
                    break
            if limit is not None and line_counter > limit:
                break

        print('INFO: Done assembling phenotypic profile comparison queries.')
        return

    def perform_owlsim_queries(self, raw1, raw2, interfile_directory, interfile_prefix, outfile_directory, outfile_prefix, num_files, limit=None):
        print('INFO: Performing OWLSim queries.')
        line_counter = 0

        if limit is not None:
            print('Only querying first '+str(limit)+' phenotypic profile pairs.')
            comparison_count = limit

        data1 = open(raw1, 'rb')
        organism_a_hash = pickle.load(data1)
        data1.close()
        data2 = open(raw2, 'rb')
        organism_b_hash = pickle.load(data2)
        data2.close()
        if limit is None:
            comparison_count = len(organism_a_hash) * len(organism_b_hash)
            print('INFO: '+str(comparison_count)+' phenotypic profile comparisons to process.')
        # Dump the hash files from memory.
        organism_a_hash = {}
        organism_b_hash = {}
        #with open(inter, 'r', encoding="iso-8859-1") as csvfile:
            #filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            #row_count = sum(1 for row in filereader)
            #print(str(row_count)+' rows to process.')
        for x in range(11, num_files+1):
            interfile = interfile_directory+'/'+interfile_prefix+'_'+str(x)+'.txt'
            outfile = outfile_directory+'/'+outfile_prefix+'_'+str(x)+'.txt'

            with open(outfile, 'w', newline='') as outfile:
                with open(interfile, 'r', encoding="iso-8859-1") as csvfile:
                    filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')

                    ###### MULTIPROCESSING INSERT ######
                    if __name__ == '__main__':
                        #lock = multiprocessing.Lock

                        print('INFO: Multiprocessing started')
                        cores = (multiprocessing.cpu_count()-1)
                        #cores = 100
                        pool = multiprocessing.Pool(processes=cores)
                        num_chunks = 10000
                        chunks = itertools.groupby(filereader, keyfunc)
                        while True:
                            queries = [list(chunk) for key, chunk in itertools.islice(chunks, num_chunks)]

                            if queries:
                                for result in pool.imap(multiprocess_owlsim_queries, queries):
                                    json.dump(result, outfile)
                                    outfile.write('\n')
                            else:
                                break
                        #multiprocessing.Semaphore(cores)
                        #jobs = []
                        #phenotype_iterable = []
                        #phenotype_counter = 0

                        #(comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes) = tuple
                        #phenotype_counter += 1
                        #print('Working on phenotype '+str(phenotype_counter)+' out of '+str(len(phenotype_list))+'.')

                        #pool.apply(multiprocess_owlsim_queries, args=(lock, row, outfile)) for row in filereader
                        #results = [pool.apply(multiprocess_owlsim_queries, args=(lock, row, outfile)) for row in filereader]


                        '''print('Processing results.')
                        comparison_list = []
                        for x in results:
                            (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)  = x
                            sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                            json.dump(sequence, outfile)
                            outfile.write('\n')
                        print('Done processing results.')'''

                    print('INFO: Multiprocessing completed')
                    ###### END MULTIPROCESSING INSERT ######

        return

    ####### PHENOLOG DATA PROCESSING #######

    def trim_panther_data(self, inter, taxons):
        """ This function trims the PANTHER flat file from NIF/DISCO for a given taxon. """

        print('INFO: Trimming PANTHER data.')

        # Set counters and open required files.
        line_counter = 0
        output_line_counter = 0
        raw = 'raw/panther/dvp.pr_nlx_84521_1'
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

                    #
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
                        #FIXME: less than or equal or just less than?
                        if prb <= fdr_cutoff:
                            significance = 'Significant'
                        else:
                            significance = 'Not Significant'
                        #print(prb)
                        # Required output : phenotype a/b, species a/b, gene list a/b, probability, fdr adjusted probability?
                        #sequence = (species_a_phenotype_id, phenotype_a_ortholog_count, species_b_phenotype_id, phenotype_b_ortholog_count, num_shared_orthologs, ortholog_matches, prb, fdr_cutoff, significance)
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
        # Take the phenotype-ortholog hashes I have created.
        # Remove all orthologs and place in a list, replace with 0s or 1s or something.
        # Randomly shuffle the ortholog list
        # Iterate through the hash and replace the 0s with an ortholog ID IF that ortholog ID is not present in the hash.

        # Question: How to be sure that the data set is random, and that I’m not creating 1000 identical data sets?
        # Should the random data set have the same presence of orthologs as the test set?
        # Meaning, if ortholog X is present in the test data set Y times,
        # should it also be present in the test data set Y times?
        # Need to have a way to make the data set creation fail if we get to the end and
        # can only put an ortholog with a phenotype that already has that ortholog with it.

        # Question: Is it appropriate to use the orthologs that are in the data set/associated with phenotypes, or
        # would it make more sense to populate from the total orthologs shared between two species?
        # What about getting random.sample from the total list of orthologs, with replacement?
        counter = 1
        while counter <= limit:
            print('INFO: Creating random data set '+str(counter)+' out of '+str(limit)+'.')
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
                    random.shuffle(ortholog_draw)
                    test_pheno_ortholog_hash[i].append(ortholog_draw[0])
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

    def set_stage_for_fdr_calculation(self):
        print('INFO: Setting stage for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.


        fdr_global_p_value_list = []

        #print(hvm_common_orthologs)
        #print(hvz_common_orthologs)
        #print(mvz_common_orthologs)

        ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':


            cores = (multiprocessing.cpu_count()-1)
            #cores = 1
            pool = multiprocessing.Pool(processes=cores)

            #multiprocessing.Semaphore(cores)
            #jobs = []
            #phenotype_iterable = []
            #phenotype_counter = 0

            #(comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes) = tuple
            #phenotype_counter += 1
            #print('Working on phenotype '+str(phenotype_counter)+' out of '+str(len(phenotype_list))+'.')
            print('INFO: Multiprocessing started')
            #random_data_set_counter = []
            #for i in range(1,5):
                #random_data_set_counter.append(i)

            fdr_global_p_value_list = [pool.map(multiprocess_fdr_calculation, range(1,1001))]


            #print('Processing results.')
            comparison_list = []
            #for p in results:
                #fdr_p_value = p.get()
                #fdr_global_p_value_list.append(fdr_p_value)

            #print('Done processing results.')

            print('INFO: Multiprocessing completed')
            ###### END MULTIPROCESSING INSERT ######
            fdr_global_p_value_list.sort()
            with open('inter/random/fdr/fdr_global_p_value_list.txt', 'wb') as handle:
                pickle.dump(fdr_global_p_value_list, handle)

            global_cutoff_position = math.ceil((len(fdr_global_p_value_list))*0.05) - 1
            print('The empirical FDR adjustment cutoff is '+str(fdr_global_p_value_list[global_cutoff_position])+'.')

        return fdr_global_p_value_list[global_cutoff_position]

    def assemble_partial_fdr(self):

        with open('inter/random/fdr/fdr_999_global_p_value_list.txt', 'rb') as handle:
            fdr_global_p_value_list = pickle.load(handle)
        with open('inter/random/fdr/fdr_1_global_p_value_list.txt', 'rb') as handle:
            fdr_extra_p_value = pickle.load(handle)

        print(fdr_global_p_value_list[0])
        print(fdr_extra_p_value[0])
        fdr_global_p_value_list = fdr_global_p_value_list[0]
        fdr_extra_p_value = fdr_extra_p_value[0]
        fdr_global_p_value_list.append(fdr_extra_p_value[0])
        print(fdr_global_p_value_list)
        print(len(fdr_global_p_value_list))
        #fdr_global_p_value_list.sort()
        p_value_sum = float(0)
        for x in range(0,1000):
            print(fdr_global_p_value_list[x])
            p_value_sum += fdr_global_p_value_list[x]
        fdr_cutoff = float(p_value_sum)/float(len(fdr_global_p_value_list))
        #global_cutoff_position = math.ceil((len(fdr_global_p_value_list))*0.05) - 1
        print('The empirical FDR adjustment cutoff is '+str(fdr_cutoff)+'.')

        with open('inter/random/fdr/fdr_global_p_value_list.txt', 'wb') as handle:
            pickle.dump(fdr_global_p_value_list, handle)

        return


    def perform_phenolog_calculations_for_fdr(self, species_a_po_hash, species_b_po_hash, shared_orthologs):
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

    def assemble_significant_phenologs(self):
        #json = open('out/phenolog/human_vs_mouse.txt')
        #data = json.loads(json)
        #json_lines = []
        with open('inter/phenolog/all_significant_phenologs.txt', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            with open('inter/phenolog/hvm_significant_phenologs.txt', 'w', newline='') as csvfile2:
                csvwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='\"')
                with open('out/phenolog/human_vs_mouse.txt', 'r') as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        phenolog_data = json.loads(json_line)
                        if phenolog_data[10] == 'Significant':
                            phenotype_a = phenolog_data[0]
                            phenotype_b = phenolog_data[3]
                            combo_ab = phenotype_a+'_'+phenotype_b
                            combo_ba = phenotype_b+'_'+phenotype_a
                            output_row = (phenotype_a, phenotype_b, combo_ab, combo_ba)
                            #print('found one')
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)
            with open('inter/phenolog/hvz_significant_phenologs.txt', 'w', newline='') as csvfile2:
                csvwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='\"')
                with open('out/phenolog/human_vs_zebrafish.txt', 'r') as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        phenolog_data = json.loads(json_line)
                        if phenolog_data[10] == 'Significant':
                            phenotype_a = phenolog_data[0]
                            phenotype_b = phenolog_data[3]
                            combo_ab = phenotype_a+'_'+phenotype_b
                            combo_ba = phenotype_b+'_'+phenotype_a
                            output_row = (phenotype_a, phenotype_b, combo_ab, combo_ba)
                            #print('found one')
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)
            with open('inter/phenolog/mvz_significant_phenologs.txt', 'w', newline='') as csvfile2:
                csvwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='\"')
                with open('out/phenolog/mouse_vs_zebrafish.txt', 'r') as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        phenolog_data = json.loads(json_line)
                        if phenolog_data[10] == 'Significant':
                            phenotype_a = phenolog_data[0]
                            phenotype_b = phenolog_data[3]
                            combo_ab = phenotype_a+'_'+phenotype_b
                            combo_ba = phenotype_b+'_'+phenotype_a
                            output_row = (phenotype_a, phenotype_b, combo_ab, combo_ba)
                            #print('found one')
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)

        return


    def assemble_significant_phenologs_with_scores(self):
        #json = open('out/phenolog/human_vs_mouse.txt')
        #data = json.loads(json)
        #json_lines = []
        with open('out/phenolog/all_significant_phenologs_scored.txt', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            with open('out/phenolog/hvm_significant_phenologs_scored.txt', 'w', newline='') as csvfile2:
                csvwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='\"')
                with open('out/phenolog/human_vs_mouse.txt', 'r') as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        phenolog_data = json.loads(json_line)
                        if phenolog_data[10] == 'Significant':
                            phenotype_a = phenolog_data[0]
                            phenotype_a_ortholog_count = phenolog_data[2]
                            phenotype_b = phenolog_data[3]
                            phenotype_b_ortholog_count = phenolog_data[5]
                            ortholog_matches = phenolog_data[7]
                            probability = phenolog_data[8]
                            output_row = (phenotype_a, phenotype_a_ortholog_count, phenotype_b, phenotype_b_ortholog_count, ortholog_matches, probability)
                            #print('found one')
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)
            print('INFO: Done assembling human-mouse phenologs.')
            with open('out/phenolog/hvz_significant_phenologs_scored.txt', 'w', newline='') as csvfile2:
                csvwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='\"')
                with open('out/phenolog/human_vs_zebrafish.txt', 'r') as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        phenolog_data = json.loads(json_line)
                        if phenolog_data[10] == 'Significant':
                            phenotype_a = phenolog_data[0]
                            phenotype_a_ortholog_count = phenolog_data[2]
                            phenotype_b = phenolog_data[3]
                            phenotype_b_ortholog_count = phenolog_data[5]
                            ortholog_matches = phenolog_data[7]
                            probability = phenolog_data[8]
                            output_row = (phenotype_a, phenotype_a_ortholog_count, phenotype_b, phenotype_b_ortholog_count, ortholog_matches, probability)
                            #print('found one')
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)
            print('INFO: Done assembling human-zebrafish phenologs.')
            with open('out/phenolog/mvz_significant_phenologs_scored.txt', 'w', newline='') as csvfile2:
                csvwriter2 = csv.writer(csvfile2, delimiter='\t', quotechar='\"')
                with open('out/phenolog/mouse_vs_zebrafish.txt', 'r') as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        phenolog_data = json.loads(json_line)
                        if phenolog_data[10] == 'Significant':
                            phenotype_a = phenolog_data[0]
                            phenotype_a_ortholog_count = phenolog_data[2]
                            phenotype_b = phenolog_data[3]
                            phenotype_b_ortholog_count = phenolog_data[5]
                            ortholog_matches = phenolog_data[7]
                            probability = phenolog_data[8]
                            output_row = (phenotype_a, phenotype_a_ortholog_count, phenotype_b, phenotype_b_ortholog_count, ortholog_matches, probability)
                            #print('found one')
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)
            print('INFO: Done assembling mouse-zebrafish phenologs.')

        return

    def assemble_hvm_phenologs(self):
        print('INFO: Assembling human-mouse phenologs.')
        hvm_phenologs = []
        with open('inter/phenolog/hvm_significant_phenologs.txt', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (phenotype_a, phenotype_b, combo_ab, combo_ba) = row
                hvm_phenologs.append(combo_ab)
        print(len(hvm_phenologs))
        unique_hvm_phenologs = set(hvm_phenologs)
        print(len(unique_hvm_phenologs))
        with open('inter/phenolog/hvm_phenolog_combo.txt', 'wb') as handle:
            pickle.dump(unique_hvm_phenologs, handle)
        print('INFO: Assembled human-mouse phenologs.')
        return

    def assemble_hvz_phenologs(self):
        print('INFO: Assembling human-zebrafish phenologs.')
        hvz_phenologs = []
        with open('inter/phenolog/hvz_significant_phenologs.txt', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (phenotype_a, phenotype_b, combo_ab, combo_ba) = row
                hvz_phenologs.append(combo_ab)
        print(len(hvz_phenologs))
        unique_hvz_phenologs = set(hvz_phenologs)
        print(len(unique_hvz_phenologs))
        with open('inter/phenolog/hvz_phenolog_combo.txt', 'wb') as handle:
            pickle.dump(unique_hvz_phenologs, handle)
        print('INFO: Assembled human-zebrafish phenologs.')
        return

    def assemble_mvz_phenologs(self):
        print('INFO: Assembling mouse-zebrafish phenologs.')
        mvz_phenologs = []
        with open('inter/phenolog/mvz_significant_phenologs.txt', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (phenotype_a, phenotype_b, combo_ab, combo_ba) = row
                mvz_phenologs.append(combo_ab)
        print(len(mvz_phenologs))
        unique_mvz_phenologs = set(mvz_phenologs)
        print(len(unique_mvz_phenologs))
        with open('inter/phenolog/mvz_phenolog_combo.txt', 'wb') as handle:
            pickle.dump(unique_mvz_phenologs, handle)
        print('INFO: Assembled mouse-zebrafish phenologs.')
        return

    def set_stage_for_extension_fdr_calculation(self, limit=1000):
        print('INFO: Setting stage for second FDR estimation.')
        # Need to calculate phenolog extension for each pairwise species and combine in order to get a full
        # set of 'genologs' (?) for proper estimation of FDR.

        ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':


            #cores = (multiprocessing.cpu_count()-1)
            cores = 2
            pool = multiprocessing.Pool(processes=cores)

            print('INFO: Multiprocessing started')

            fdr_global_p_value_list = [pool.map(multiprocess_ext_fdr_calculation, range(1, 11))]
            with open('inter/random/fdr_ext/fdr_ext_global_p_value_list.txt', 'wb') as handle:
                pickle.dump(fdr_global_p_value_list, handle)
            print('INFO: Multiprocessing completed')
            ###### END MULTIPROCESSING INSERT ######
            '''
            p_value_sum = float(0)
            for x in range(0,1000):
                #print(fdr_global_p_value_list[x])
                p_value_sum += fdr_global_p_value_list[x]
            fdr_cutoff = float(p_value_sum)/float(len(fdr_global_p_value_list))
            #global_cutoff_position = math.ceil((len(fdr_global_p_value_list))*0.05) - 1
            print('The empirical FDR adjustment cutoff is '+str(fdr_cutoff)+'.')
            '''
            #fdr_global_p_value_list.sort()
            #global_cutoff_position = math.ceil((len(fdr_global_p_value_list))*0.05) - 1
            #if global_cutoff_position < 0:
                #global_cutoff_position = 0
            #print('The empirical FDR adjustment cutoff is '+str(fdr_global_p_value_list[global_cutoff_position])+'.')

        return #fdr_global_p_value_list[global_cutoff_position]



    def generate_human_random_ext_data(self):
        print('INFO: Creating random data sets.')


            ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)
            print('INFO: Multiprocessing started')

            pool.map(multiprocess_generate_random_human_ext_data, range(1,1001))

            print('INFO: Multiprocessing completed')
            ###### END MULTIPROCESSING INSERT ######

        return

    def generate_mouse_random_ext_data(self):
        print('INFO: Creating random data sets.')


            ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)
            print('INFO: Multiprocessing started')

            pool.map(multiprocess_generate_random_mouse_ext_data, range(1,1001))

            print('INFO: Multiprocessing completed')
            ###### END MULTIPROCESSING INSERT ######

        return

    def generate_zebrafish_random_ext_data(self):
        print('INFO: Creating random data sets.')


            ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)
            print('INFO: Multiprocessing started')

            pool.map(multiprocess_generate_random_zebrafish_ext_data, range(1,1001))

            print('INFO: Multiprocessing completed')
            ###### END MULTIPROCESSING INSERT ######

        return


    def perform_phenolog_calculations_for_ext_fdr(self, species_a_gp_hash, species_b_gp_hash, shared_phenologs):
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        line_counter = 0
        failure_counter = 0
        comparison_counter = 0
        total_phenotype_matches = 0
        total_phenotype_nonmatches = 0
        phenotype_matches = 0
        phenotype_non_matches = 0
        genotype_a_phenotype_count = 0
        genotype_b_phenotype_count = 0
        total_hyp_calcs = 0
        phenolog_ext_p_value_list = []

        species_a_geno_pheno_hash = species_a_gp_hash
        species_b_geno_pheno_hash = species_b_gp_hash

        total_comparisons = (len(species_a_geno_pheno_hash))*(len(species_b_geno_pheno_hash))
        print('INFO: '+str(total_comparisons)+' total comparisons to perform.')

        for i in species_a_geno_pheno_hash:
            # Genotype for species A
            species_a_genotype_id = i
            species_a_phenotypes = species_a_geno_pheno_hash[i]
            #print(species_a_phenotypes)
            genotype_a_phenotype_count = len(species_a_phenotypes)

            for j in species_b_geno_pheno_hash:
                # Genotype for species B
                species_b_genotype_id = j
                species_b_phenotypes = species_b_geno_pheno_hash[j]
                #print(species_b_phenotypes)
                #phenotype_matches = 0
                #phenotype_non_matches = 0

                phenotype_matches = 0
                phenotype_non_matches = 0

                genotype_b_phenotype_count = len(species_b_phenotypes)

                comparison_counter += 1
                if re.match('.*000000', str(comparison_counter)):
                    print('INFO: Processing comparison '+str(comparison_counter)+' out of '+str(total_comparisons)+' total comparisons to perform.')
                for k in species_a_phenotypes:
                    # Orthologs for species A
                    #ortholog_matches = 0
                    #ortholog_non_matches = 0

                    species_a_phenotype = k
                    for l in species_b_phenotypes:
                        # Orthologs for species B
                        species_b_phenotype = l

                        ab_combo = species_a_phenotype+'_'+species_b_phenotype
                        ba_combo = species_b_phenotype+'_'+species_a_phenotype
                        if ab_combo in shared_phenologs or ba_combo in shared_phenologs:
                            #print('species a ortholog:'+species_a_ortholog+' matches species b ortholog:'+species_b_ortholog)
                            phenotype_matches += 1
                            #print(species_a_ortholog+' == '+species_b_ortholog)
                            total_phenotype_matches += 1
                        else:
                            #print('species a ortholog:'+species_a_ortholog+' does not match species b ortholog:'+species_b_ortholog)
                            phenotype_non_matches += 1
                            total_phenotype_nonmatches += 1

                if phenotype_matches > 0:
                    #print('Matches: '+str(ortholog_matches))
                    #print('Non-matches: '+str(ortholog_non_matches))
                    m = float(genotype_b_phenotype_count)
                    n = float(genotype_a_phenotype_count)
                    N = float(len(shared_phenologs))
                    c = float(phenotype_matches)
                    prb = float(hypergeom.pmf(c, N, m, n))
                    #print(str(c)+', '+str(N)+', '+str(m)+', '+str(n))
                    #print(prb)
                    phenolog_ext_p_value_list.append(prb)
                    total_hyp_calcs += 1
        print('Total Matches: '+str(total_phenotype_matches))
        print('Total non-matches: '+str(total_phenotype_nonmatches))
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

        return phenolog_ext_p_value_list

    def perform_phenolog_calculations_for_ext_fdr_alternate(self, species_a_gp_hash, species_b_gp_hash, out_file):
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        #line_counter = 0
        #failure_counter = 0
        #comparison_counter = 0
        #total_phenotype_matches = 0
        #total_phenotype_nonmatches = 0
        #phenotype_matches = 0
        #phenotype_non_matches = 0
        #genotype_a_phenotype_count = 0
        #genotype_b_phenotype_count = 0
        #total_hyp_calcs = 0

        phenolog_ext_p_value_list = []

        species_a_geno_pheno_hash = species_a_gp_hash
        species_b_geno_pheno_hash = species_b_gp_hash

        total_comparisons = (len(species_a_geno_pheno_hash))*(len(species_b_geno_pheno_hash))
        print('INFO: '+str(total_comparisons)+' total comparisons to perform.')

        print('INFO: Assembling disease/genotype comparison list.')

        comparison_list = []

        for element in itertools.product(species_a_geno_pheno_hash,species_b_geno_pheno_hash):
            comparison_list.append(element)

        print('INFO: Done assembling disease/genotype comparison list.')

        if __name__ == '__main__':

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)

            print('INFO: Starting multiprocessing.')

            for results in pool.imap_unordered(multiprocess_ext_fdr_calculation_alternate, comparison_list, chunksize=100):
                if results is not None:
                    phenolog_ext_p_value_list.append(results)
                    #print(results)


            #print('INFO: Done processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')

        with open(out_file, 'wb') as handle:
            pickle.dump(phenolog_ext_p_value_list, handle)
        #print('Total Matches: '+str(total_phenotype_matches))
        #print('Total non-matches: '+str(total_phenotype_nonmatches))
        #print('Total phenolog calculations: '+str(total_hyp_calcs))

        phenolog_ext_p_value_list = []
        comparison_list = []
        gc.collect()
        return #phenolog_ext_p_value_list


    def perform_phenolog_ext_calculations(self, inter1, inter2, out, shared_phenologs, ext_fdr_cutoff):
        print('INFO: Performing phenolog extension calculations.')
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        line_counter = 0
        failure_counter = 0
        total_phenotype_matches = 0
        total_phenotype_nonmatches = 0
        phenotype_matches = 0
        phenotype_non_matches = 0
        genotype_a_phenotype_count = 0
        genotype_b_phenotype_count = 0
        total_hyp_calcs = 0
        phenolog_ext_p_value_list = []

        num_shared_phenologs = len(shared_phenologs)
        with open(inter1, 'rb') as handle:
            species_a_geno_pheno_hash = pickle.load(handle)
        #print(species_a_pheno_gene_hash)
        with open(inter2, 'rb') as handle:
            species_b_geno_pheno_hash = pickle.load(handle)
        with open(out, 'w', newline='') as outfile:

            for i in species_a_geno_pheno_hash:
                # Genotype for species A
                species_a_genotype_id = i
                species_a_phenotypes = species_a_geno_pheno_hash[i]
                #print(species_a_phenotypes)
                genotype_a_phenotype_count = len(species_a_phenotypes)

                for j in species_b_geno_pheno_hash:
                    # Genotype for species B
                    species_b_genotype_id = j
                    species_b_phenotypes = species_b_geno_pheno_hash[j]
                    #print(species_b_phenotypes)
                    #phenotype_matches = 0
                    #phenotype_non_matches = 0

                    phenotype_matches = 0
                    phenotype_non_matches = 0

                    genotype_b_phenotype_count = len(species_b_phenotypes)
                    for k in species_a_phenotypes:
                        # Orthologs for species A
                        #ortholog_matches = 0
                        #ortholog_non_matches = 0

                        species_a_phenotype = k
                        for l in species_b_phenotypes:
                            # Orthologs for species B
                            species_b_phenotype = l

                            ab_combo = species_a_phenotype+'_'+species_b_phenotype
                            ba_combo = species_b_phenotype+'_'+species_a_phenotype
                            if ab_combo in shared_phenologs or ba_combo in shared_phenologs:
                                #print('species a ortholog:'+species_a_ortholog+' matches species b ortholog:'+species_b_ortholog)
                                phenotype_matches += 1
                                #print(species_a_ortholog+' == '+species_b_ortholog)
                                total_phenotype_matches += 1
                            else:
                                #print('species a ortholog:'+species_a_ortholog+' does not match species b ortholog:'+species_b_ortholog)
                                phenotype_non_matches += 1
                                total_phenotype_nonmatches += 1

                    if phenotype_matches > 0:
                        #print('Matches: '+str(ortholog_matches))
                        #print('Non-matches: '+str(ortholog_non_matches))
                        m = float(genotype_b_phenotype_count)
                        n = float(genotype_a_phenotype_count)
                        N = float(num_shared_phenologs)
                        c = float(phenotype_matches)
                        prb = float(hypergeom.pmf(c, N, m, n))
                        #print(str(c)+', '+str(N)+', '+str(m)+', '+str(n))
                        #print(prb)
                        if prb <= ext_fdr_cutoff:
                            significance = 'Significant'
                        else:
                            significance = 'Not Significant'

                        sequence = (species_a_genotype_id, species_a_phenotypes, genotype_a_phenotype_count, species_b_genotype_id, species_b_phenotypes, genotype_b_phenotype_count, num_shared_phenologs, phenotype_matches, prb, ext_fdr_cutoff, significance)
                        json.dump(sequence, outfile)
                        outfile.write('\n')
                        total_hyp_calcs += 1
        print('Total matches: '+str(total_phenotype_matches))
        print('Total non-matches: '+str(total_phenotype_nonmatches))
        print('Total phenolog extension calculations: '+str(total_hyp_calcs))
        return

    def parse_zp(self, raw, inter):
        zp_hash = {}
        line_counter = 0
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                line_counter += 1
                (zp_id, zp_label, extra_1, extra_2, extra_3, extra_4, extra_5, extra_6) = row

                if zp_id not in zp_hash:
                    zp_hash[zp_id] = zp_label

        with open(inter, 'wb') as handle:
            pickle.dump(zp_hash, handle)



        return

    def parse_mp(self, raw, inter):
        mp_hash = {}
        with open(raw, 'r') as handle:
            for line in handle:
                #print(line)
                if re.match('id:.*', line.rstrip()):
                    mp_id = line.rstrip()
                    mp_id = re.sub('id: ', '', mp_id)
                    #print(mp_id)
                    next_line = next(handle)
                    if re.match('name:.*', next_line.rstrip()):
                        mp_label = next_line.rstrip()
                        mp_label = re.sub('name: ', '', mp_label)
                        #print(mp_label)

                        if mp_id not in mp_hash:
                            mp_hash[mp_id] = mp_label

        with open(inter, 'wb') as handle:
            pickle.dump(mp_hash, handle)

        return

    def parse_hp(self, raw, inter):
        hp_hash = {}
        with open(raw, 'r') as handle:
            for line in handle:
                #print(line)
                if re.match('id:.*', line.rstrip()):
                    hp_id = line.rstrip()
                    hp_id = re.sub('id: ', '', hp_id)
                    #print(hp_id)
                    next_line = next(handle)
                    if re.match('name:.*', next_line.rstrip()):
                        hp_label = next_line.rstrip()
                        hp_label = re.sub('name: ', '', hp_label)
                        #print(hp_label)

                        if hp_id not in hp_hash:
                            hp_hash[hp_id] = hp_label

        with open(inter, 'wb') as handle:
            pickle.dump(hp_hash, handle)

        return

    ####### PHENOLOG ORTHOLOG PHENOTYPE MATRICES #######

    def assemble_ortholog_phenotype_matrix(self, limit=None):
        # To speed things up, going to assemble from already assembled phenotype-ortholog hash

        print('INFO: Assembling phenotype-ortholog matrices.')
        #How to do this?



        # Create gene list
        # Convert to ortholog list (unique)
        # Create phenotype-ortholog matrix filled with zeroes
        # Cycle through table, set the phenotype-ortholog matrix

        # Create gene-ortholog hash (Create list of genes, update with ortholog as value)
        # Create phenotype list
        # Cycle through
        phenotype_list = []
        ortholog_list = []
        phenotype_limit = limit
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            #print(species_a_pheno_gene_hash)

            phenotype_counter = 0
            for i in pheno_ortholog_hash:
                #if phenotype_counter < phenotype_limit:
                    #phenotype_counter += 1
                if i not in phenotype_list:
                    phenotype_list.append(i)
                for j in pheno_ortholog_hash[i]:
                    if j not in ortholog_list:
                        ortholog_list.append(j)
                #else:
                    #continue
        phenotype_list.sort()
        ortholog_list.sort()
        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))
        #print(phenotype_list[5])
        #print(phenotype_list.index('ZP:0000007'))
        #print(phenotype_list[5])
        ortholog_phenotype_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            phenotype_counter = 0
            for i in pheno_ortholog_hash:
                #if phenotype_counter < phenotype_limit:
                    #phenotype_counter += 1
                phenotype_index = phenotype_list.index(i)
                for j in pheno_ortholog_hash[i]:
                    ortholog_index = ortholog_list.index(j)
                    ortholog_phenotype_matrix[phenotype_index][ortholog_index] = 1
                #else:
                    #continue

        #print(phenotype_list[0])
        #print(ortholog_phenotype_matrix)
        #print(len(phenotype_ortholog_matrix))
        #a = numpy.matrix(phenotype_list, ortholog_list)
        #print(str(numpy.sum(phenotype_ortholog_matrix)))

        print('INFO: Done assembling phenotype-ortholog matrices.')


        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
            pickle.dump(phenotype_list, handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
            pickle.dump(ortholog_list, handle)
        numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
        numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)

        return (ortholog_phenotype_matrix, phenotype_list, ortholog_list)

    def assemble_ortholog_phenotype_matrix_alternate(self, limit=None):
        # To speed things up, going to assemble from already assembled phenotype-ortholog hash

        print('INFO: Assembling phenotype-ortholog matrices.')
        #How to do this?



        # Create gene list
        # Convert to ortholog list (unique)
        # Create phenotype-ortholog matrix filled with zeroes
        # Cycle through table, set the phenotype-ortholog matrix

        # Create gene-ortholog hash (Create list of genes, update with ortholog as value)
        # Create phenotype list
        # Cycle through
        phenotype_list = []
        ortholog_list = []
        phenotype_limit = limit
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            #print(species_a_pheno_gene_hash)

            phenotype_counter = 0
            for i in pheno_ortholog_hash:
                if phenotype_counter < phenotype_limit:
                    phenotype_counter += 1
                    if i not in phenotype_list:
                        phenotype_list.append(i)
                    for j in pheno_ortholog_hash[i]:
                        if j not in ortholog_list:
                            ortholog_list.append(j)
                else:
                    continue
        phenotype_list.sort()
        ortholog_list.sort()
        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))
        #print(phenotype_list[5])
        #print(phenotype_list.index('ZP:0000007'))
        #print(phenotype_list[5])
        ortholog_phenotype_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            phenotype_counter = 0
            for i in pheno_ortholog_hash:
                if phenotype_counter < phenotype_limit:
                    phenotype_counter += 1
                    phenotype_index = phenotype_list.index(i)
                    for j in pheno_ortholog_hash[i]:
                        ortholog_index = ortholog_list.index(j)
                        ortholog_phenotype_matrix[phenotype_index][ortholog_index] = 1
                else:
                    continue

        #print(phenotype_list[0])
        #print(ortholog_phenotype_matrix)
        #print(len(phenotype_ortholog_matrix))
        #a = numpy.matrix(phenotype_list, ortholog_list)
        #print(str(numpy.sum(phenotype_ortholog_matrix)))

        print('INFO: Done assembling phenotype-ortholog matrices.')


        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
            pickle.dump(phenotype_list, handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
            pickle.dump(ortholog_list, handle)
        numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
        numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)

        return (ortholog_phenotype_matrix, phenotype_list, ortholog_list)


    ####### PHENOLOG GENE CANDIDATE PREDICTION ALGORITHM #######
    # Multiple things to do here.
    #

    def create_phenolog_gene_candidate_matrices(self):

        #Testing - read in the matrix, compare matrix columns

        # Set the number of nearest neighbor phenotypes to consider for predictions.
        k = 11
        # Added file dumps to main.assemble_ortholog_phenotype_matrix, so this function call is not necessary if already run.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()
        # If testing, set a limit to something reasonable like 100.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix(500)

        ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')

        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)

        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))
        #Have the matrix, need to get the sum of the k (10) nearest neighbors, weight by the pearson correlation coefficient.
        # Pearson correlation to determine the k nearest neighbors. So I need to calculate the similarity of phenotypes
        # for all pair-wise phenotype combinations. So need a similarity score matrix in addition to the weight matrix.
        # Use the hypergeometric CDF to provide scores for the weight matrix.

        #Creating a small test matrix for testing.
        test_matrix = numpy.random.randint(2, size=(10,10))
        #ortholog_phenotype_matrix = test_matrix
        test_phenotype_list = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']
        #phenotype_list = test_phenotype_list
        test_ortholog_list = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O10']
        #ortholog_list = test_ortholog_list
        total_phenotypes = len(phenotype_list)


        #print(phenotype_list[0])
        distance_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        distance_matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        distance_matrix_counter = 0
        weight_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        weight_matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        weight_matrix_counter = 0
        total_orthologs = len(ortholog_list)

        print('INFO: '+str(distance_matrix_comparisons)+' matrix comparisons to process.')

        #Need to multi-process this segment.


        if __name__ == '__main__':
            #with Manager() as manager:

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)

            #multiprocessing.Semaphore(cores)
            #jobs = []
            #phenotype_iterable = []
            phenotype_counter = 0

            #Takes ~65 seconds to reach this point.
            print('INFO: Assembling phenotype matrix coordinates.')

            for i in phenotype_list:
                #i = phenotype_list[0]
                input_phenotype_index_i = phenotype_list.index(i)

                print('INFO: Processing phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                matrix_coordinates = []
                for j in phenotype_list:
                    input_phenotype_index_j = phenotype_list.index(j)
                    matrix_coordinates.append([input_phenotype_index_i,input_phenotype_index_j])
                print('INFO: Done assembling phenotype matrix coordinates.')
                print('INFO: Starting multiprocessing.')
                results = pool.map(multiprocess_matrix_comparisons, matrix_coordinates)
                #print(results)
                print('INFO: Processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                for x in results:
                    #print(x)
                    (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient) = x
                    distance_matrix[phenotype_index_i][phenotype_index_j] = coefficient
                    weight_matrix[phenotype_index_i][phenotype_index_j] = hyp_prob
                    #print(coefficient)
                print('INFO: Done processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')

            # Dump all of the files to disk.
            #numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)
            print('INFO: Dumping distance matrix to disk.')
            #numpy.save('inter/phenolog_gene_cand/distance_matrix.npy', distance_matrix)
            #FYI: The human readable matrix file is 3X the size of the .npy file.
            #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
            print('INFO: Dumping weight matrix to disk.')
            #numpy.save('inter/phenolog_gene_cand/weight_matrix.npy', weight_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/weight_matrix.txt', weight_matrix)
            #with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
                #pickle.dump(phenotype_list, handle)
            #print(phenotype_list[0])
            #with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
                #pickle.dump(ortholog_list, handle)
            #print(ortholog_list[0])
            print('DONE!')
        return


    def create_empty_phenolog_gene_candidate_matrices(self):

        #Testing - read in the matrix, compare matrix columns

        # Set the number of nearest neighbor phenotypes to consider for predictions.
        k = 11
        # Added file dumps to main.assemble_ortholog_phenotype_matrix, so this function call is not necessary if already run.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()
        # If testing, set a limit to something reasonable like 100.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()

        #ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')

        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)

        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))
        #Have the matrix, need to get the sum of the k (10) nearest neighbors, weight by the pearson correlation coefficient.
        # Pearson correlation to determine the k nearest neighbors. So I need to calculate the similarity of phenotypes
        # for all pair-wise phenotype combinations. So need a similarity score matrix in addition to the weight matrix.
        # Use the hypergeometric CDF to provide scores for the weight matrix.

        distance_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        distance_matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        weight_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))

        print('INFO: '+str(distance_matrix_comparisons)+' matrix comparisons to process.')

        # Dump all of the files to disk.
        print('INFO: Dumping empty distance matrix to disk.')
        numpy.save('inter/phenolog_gene_cand/distance_matrix.npy', distance_matrix)
        #FYI: The human readable matrix file is 3X the size of the .npy file.
        #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
        print('INFO: Dumping empty weight matrix to disk.')
        numpy.save('inter/phenolog_gene_cand/weight_matrix.npy', weight_matrix)

        return

    def populate_phenolog_gene_candidate_matrices(self):

        #Testing - read in the matrix, compare matrix columns

        # Set the number of nearest neighbor phenotypes to consider for predictions.
        k = 11
        # Added file dumps to main.assemble_ortholog_phenotype_matrix, so this function call is not necessary if already run.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()
        # If testing, set a limit to something reasonable like 100.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()

        ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')
        #print(ortholog_phenotype_matrix[1][1])

        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)

        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))
        #Have the matrix, need to get the sum of the k (10) nearest neighbors, weight by the pearson correlation coefficient.
        # Pearson correlation to determine the k nearest neighbors. So I need to calculate the similarity of phenotypes
        # for all pair-wise phenotype combinations. So need a similarity score matrix in addition to the weight matrix.
        # Use the hypergeometric CDF to provide scores for the weight matrix.

        #Creating a small test matrix for testing.
        #test_matrix = numpy.random.randint(2, size=(10,10))
        #ortholog_phenotype_matrix = test_matrix
        #test_phenotype_list = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']
        #phenotype_list = test_phenotype_list
        #test_ortholog_list = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O10']
        #ortholog_list = test_ortholog_list
        total_phenotypes = len(phenotype_list)


        #print(phenotype_list[0])
        #distance_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        distance_matrix_counter = 0
        #weight_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        weight_matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        weight_matrix_counter = 0
        total_orthologs = len(ortholog_list)

        print('INFO: '+str(matrix_comparisons)+' matrix comparisons to process.')

        #Need to multi-process this segment.


        if __name__ == '__main__':
            #with Manager() as manager:

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)

            #multiprocessing.Semaphore(cores)
            #jobs = []
            #phenotype_iterable = []
            phenotype_counter = 0

            #Takes ~65 seconds to reach this point.
            print('INFO: Assembling phenotype matrix coordinates.')

            for x in range(0, len(phenotype_list)): #len(phenotype_list)
                i = phenotype_list[x]
                input_phenotype_index_i = phenotype_list.index(i)

                print('INFO: Processing phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                print('INFO: Processing phenotype '+i+'.')
                matrix_coordinates = []
                for j in phenotype_list:
                    input_phenotype_index_j = phenotype_list.index(j)
                    matrix_coordinates.append([input_phenotype_index_i,input_phenotype_index_j])
                print('INFO: Done assembling phenotype matrix coordinates.')
                print('INFO: Starting multiprocessing.')
                results = pool.map(multiprocess_matrix_comparisons, matrix_coordinates)
                #print(results)
                print('INFO: Processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
                weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')


                for x in results:
                    #print(x)
                    (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient) = x
                    distance_matrix[phenotype_index_i][phenotype_index_j] = coefficient
                    weight_matrix[phenotype_index_i][phenotype_index_j] = hyp_prob
                    #print(coefficient)
                print('INFO: Done processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')

                # Dump all of the files to disk.
                numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
                #numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)
                print('INFO: Dumping distance matrix to disk.')
                numpy.save('inter/phenolog_gene_cand/distance_matrix.npy', distance_matrix)
                #FYI: The human readable matrix file is 3X the size of the .npy file.
                #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
                print('INFO: Dumping weight matrix to disk.')
                numpy.save('inter/phenolog_gene_cand/weight_matrix.npy', weight_matrix)
                #numpy.savetxt('inter/phenolog_gene_cand/weight_matrix.txt', weight_matrix)

            #numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/weight_matrix.txt', weight_matrix)
            #with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
                #pickle.dump(phenotype_list, handle)
            #print(phenotype_list[0])
            #with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
                #pickle.dump(ortholog_list, handle)
            #print(ortholog_list[0])
            print('DONE!')
        return

    def populate_phenolog_gene_candidate_matrices_alternate(self):

        #Testing - read in the matrix, compare matrix columns

        # Set the number of nearest neighbor phenotypes to consider for predictions.
        k = 11
        # Added file dumps to main.assemble_ortholog_phenotype_matrix, so this function call is not necessary if already run.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()
        # If testing, set a limit to something reasonable like 100.
        #(ortholog_phenotype_matrix, phenotype_list, ortholog_list) = main.assemble_ortholog_phenotype_matrix()

        ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')
        #print(ortholog_phenotype_matrix[1][1])

        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)

        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))
        #Have the matrix, need to get the sum of the k (10) nearest neighbors, weight by the pearson correlation coefficient.
        # Pearson correlation to determine the k nearest neighbors. So I need to calculate the similarity of phenotypes
        # for all pair-wise phenotype combinations. So need a similarity score matrix in addition to the weight matrix.
        # Use the hypergeometric CDF to provide scores for the weight matrix.

        #Creating a small test matrix for testing.
        #test_matrix = numpy.random.randint(2, size=(10,10))
        #ortholog_phenotype_matrix = test_matrix
        #test_phenotype_list = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']
        #phenotype_list = test_phenotype_list
        #test_ortholog_list = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O10']
        #ortholog_list = test_ortholog_list
        total_phenotypes = len(phenotype_list)


        #print(phenotype_list[0])
        #distance_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        distance_matrix_counter = 0
        #weight_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        weight_matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        weight_matrix_counter = 0
        total_orthologs = len(ortholog_list)

        print('INFO: '+str(matrix_comparisons)+' matrix comparisons to process.')

        #Need to multi-process this segment.


        if __name__ == '__main__':
            #with Manager() as manager:

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)

            #multiprocessing.Semaphore(cores)
            #jobs = []
            #phenotype_iterable = []
            phenotype_counter = 0

            print('INFO: Loading distance matrix.')
            distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
            print('INFO: Loading weight matrix.')
            weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')

            #Takes ~65 seconds to reach this point.
            print('INFO: Assembling phenotype matrix coordinates.')


            xcounter = 0
            for x in range(22000, len(phenotype_list)): #len(phenotype_list)
                xcounter += 1
                i = phenotype_list[x]
                input_phenotype_index_i = phenotype_list.index(i)

                print('INFO: Processing phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                print('INFO: Processing phenotype '+i+'.')
                matrix_coordinates = []
                for j in phenotype_list:
                    input_phenotype_index_j = phenotype_list.index(j)
                    matrix_coordinates.append([input_phenotype_index_i,input_phenotype_index_j])
                print('INFO: Done assembling phenotype matrix coordinates.')
                print('INFO: Starting multiprocessing.')

                #ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')
                #for results in pool.imap_unordered(multiprocess_matrix_comparisons, matrix_coordinates, chunksize=100):
                #distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
                #weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')
                for results in pool.imap_unordered(multiprocess_matrix_comparisons, matrix_coordinates, chunksize=1000):
                    (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient) = results
                    distance_matrix[phenotype_index_i][phenotype_index_j] = coefficient
                    weight_matrix[phenotype_index_i][phenotype_index_j] = hyp_prob

                print('INFO: Done processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                if xcounter == 1000:
                    print('INFO: Dumping distance matrix to disk.')
                    numpy.save('inter/phenolog_gene_cand/distance_matrix.npy', distance_matrix)
                    #FYI: The human readable matrix file is 3X the size of the .npy file.
                    #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
                    print('INFO: Dumping weight matrix to disk.')
                    numpy.save('inter/phenolog_gene_cand/weight_matrix.npy', weight_matrix)
                    xcounter = 0

            # Dump all of the files to disk.
            numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)
            print('INFO: Dumping distance matrix to disk.')
            numpy.save('inter/phenolog_gene_cand/distance_matrix.npy', distance_matrix)
            #FYI: The human readable matrix file is 3X the size of the .npy file.
            #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
            print('INFO: Dumping weight matrix to disk.')
            numpy.save('inter/phenolog_gene_cand/weight_matrix.npy', weight_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/weight_matrix.txt', weight_matrix)

            #numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/distance_matrix.txt', distance_matrix)
            #numpy.savetxt('inter/phenolog_gene_cand/weight_matrix.txt', weight_matrix)
            #with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
                #pickle.dump(phenotype_list, handle)
            #print(phenotype_list[0])
            #with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
                #pickle.dump(ortholog_list, handle)
            #print(ortholog_list[0])
            print('DONE!')
        return



    def merge_matrices(self): #, input_matrix_1, input_matrix_2, output_matrix
        input_matrix_1 = numpy.random.randint(2, size=(10,10))
        print(input_matrix_1)
        input_matrix_2 = numpy.random.randint(2, size=(10,10))
        print(input_matrix_2)
        output_matrix = numpy.zeros((len(input_matrix_1), len(input_matrix_1)))
        #print(merged_matrix)
        #print(len(input_matrix_1))
        #print(input_matrix_1[0][1])
        for x in range(0, len(input_matrix_1)):
            for y in range(0, len(input_matrix_1)):
                output_matrix[x][y] = max(input_matrix_1[x][y], input_matrix_2[x][y])

        print(output_matrix)
        return


    def create_phenolog_gene_candidate_prediction_matrix(self):


        #FIXME: This needs to be adjusted to ensure that the nearest neighbor phenotypes are not in the same species!!!
        # What would work better? Including a prefix check for the nearest neighbor,
        # or splitting the main matrix into a few submatrices for the comparison phenotypes?
        distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
        weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')
        ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')


        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)
        #print(phenotype_list[0])
        #print(ortholog_list[0])

        phenotype_ortholog_prediction_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))


        # Want to get the 10 nearest neighbors for a given phenotype.
        # Pass in phenotype indice
        # Get the slice for the phenotype indice.
        # Create a clone of the slice, sort, and get the value of the top k entries.
        #

        #Test phenotype index choose 0
        #y = 7
        for y in range(0, len(phenotype_list)):

            test_phenotype = ortholog_phenotype_matrix[y]
            test_distance_slice = distance_matrix[y]
            intermediate_nearest_neighbors = heapq.nlargest(11, range(len(test_distance_slice)), test_distance_slice.take)

            #test_distance_slice[0] = 'X'
            test_weight_slice = weight_matrix[0]
            #test_weight_slice[0] = 'X'

            #print(test_phenotype)
            #print(distance_matrix[y])

            #print(weight_matrix[y])
            nearest_neighbors = []
            # This will print out the phenotype IDs for the k nearest neighbors
            for z in intermediate_nearest_neighbors:
                if z != y:
                     nearest_neighbors.append(z)
            #print(intermediate_nearest_neighbors)
            #print(nearest_neighbors)
            print('Input phenotype: '+phenotype_list[y])
            #print(ortholog_phenotype_matrix[y])
            #for m in nearest_neighbors:
                #print('Nearest neighbor phenotype: '+phenotype_list[m])
                #print(ortholog_phenotype_matrix[m])
            #print(neighbors)


            # Next I need to take those k nearest neighbor phenotypes, and calculate the probability that
            # ortholog i is associated with phenotype j based on those k phenotypes,
            # using the phenotype matrix and weight matrix.
            # Take the sum

            # i in nearest neighbors is the phenotype index
            for i in nearest_neighbors:
                nearby_phenotype = ortholog_phenotype_matrix[i]
                for j in range(0, len(nearby_phenotype)):
                    if ortholog_phenotype_matrix[i][j] != 0:
                        phenotype_ortholog_prediction_matrix[y][j] += weight_matrix[i][j]*ortholog_phenotype_matrix[i][j]

            #print(phenotype_ortholog_prediction_matrix[y])


        numpy.save('inter/phenolog_gene_cand/phenotype_ortholog_prediction_matrix.npy', phenotype_ortholog_prediction_matrix)
        numpy.savetxt('inter/phenolog_gene_cand/phenotype_ortholog_prediction_matrix.txt', phenotype_ortholog_prediction_matrix)



        return

    def assemble_phenolog_gene_candidate_predictions(self):
        #Have ortholog predictions in matrix for each phenotype, now need to take the matrix data and output to a table
        # that is human readable for each phenotype for comparison/assembly with OWLSim output. Maybe write to JSON as well?

        phenotype_ortholog_prediction_matrix = numpy.load('inter/phenolog_gene_cand/phenotype_ortholog_prediction_matrix.npy')
        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)
        with open('out/phenolog_gene_cand/phenolog_ortholog_candidate_predictions.txt', 'w', newline='\n') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            phenotype_ortholog_candidate_hash = {}
            # For each phenotype, grab the corresponding slice from the prediction matrix. If an ortholog entry is not == 0,
            # then assemble the PANTHER IDs in an output table/file.
            #test_phenotype_index = 5
            phenotype_index_counter = 0
            #test_phenotype_array = phenotype_ortholog_prediction_matrix[test_phenotype_index]
            for phenotype_array in phenotype_ortholog_prediction_matrix:
                ortholog_predictions = []
                ortholog_index_counter = 0
                #print(phenotype_array)
                phenotype_id = phenotype_list[phenotype_index_counter]
                phenotype_ortholog_candidate_hash[phenotype_id] = {}
                for additive_probability in phenotype_array:
                    #if phenotype_ortholog_candidate_hash[phenotype_id] not in phenotype_ortholog_candidate_hash:

                    if additive_probability != 0:
                        #ortholog_index = phenotype_list.index(additive_probability)
                        ortholog = ortholog_list[ortholog_index_counter]
                        ortholog_predictions.append(ortholog)
                        phenotype_ortholog_candidate_hash[phenotype_id][ortholog] = additive_probability
                    ortholog_index_counter += 1
                print(phenotype_list[phenotype_index_counter])
                print(ortholog_predictions)
                phenotype_index_counter += 1
                #print(phenotype_ortholog_candidate_hash)



                output_row = (phenotype_id, ortholog_predictions, phenotype_ortholog_candidate_hash[phenotype_id])
                csvwriter.writerow(output_row)

        with open('out/phenolog_gene_cand/phenolog_ortholog_candidate_prediction_hash.txt', 'wb') as handle:
            pickle.dump(phenotype_ortholog_candidate_hash, handle)

        #with open('out/phenolog_gene_cand/mouse_genotype_zebrafish_genotype.txt', 'w', newline='') as outfile:

            #line_counter += 1
            #print('INFO: Processing phenotypic profile comparison '+str(line_counter)+' out of '+str(comparison_count)+'.')

            #sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
            #json.dump(sequence, outfile)
            #outfile.write('\n')


        return

counter = multiprocessing.Value(c_int)
counter_lock = multiprocessing.Lock()


def increment():
    """ This function provides a shared counter for multiprocessing funtctions for tracking progress."""
    with counter_lock:
        counter.value += 1
        print('INFO: Processing comparison '+str(counter.value))


def keyfunc(row):
    return row[0]


def multiprocess_matrix_comparisons(matrix_coordinates):
    #increment()
    #print(matrix_coordinates)
    #with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
        #ortholog_list = pickle.load(handle)

    #Total number of orthologs is 2905. Hard coding to remove multiple openings of the ortholog_list.txt file.
    len_ortholog_list = 2905
    # For testing with smaller phenotype sets, set number of orthologs from saved file
    #with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            #ortholog_list = pickle.load(handle)
    #len_ortholog_list = len(ortholog_list)
    #Comment out when done testing.
    #test_phenotype_list = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'P9', 'P10']
    #phenotype_list = test_phenotype_list
    #test_ortholog_list = ['O1', 'O2', 'O3', 'O4', 'O5', 'O6', 'O7', 'O8', 'O9', 'O10']
    #ortholog_list = test_ortholog_list

    phenotype_index_i = matrix_coordinates[0]
    phenotype_index_j = matrix_coordinates[1]
    #distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
    #weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')
    #print('Loading ortholog phenotype matrix.')
    #ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')
    #print(len(ortholog_phenotype_matrix))
    #print('Done loading ortholog phenotype matrix.')

    ortholog_counter = 0
    ortholog_match = 0
    #print(len(ortholog_list))

    (coefficient, p_value) = pearsonr(read_only_ortholog_phenotype_matrix[phenotype_index_i], read_only_ortholog_phenotype_matrix[phenotype_index_j])
    #print(str(coeffecient)+'_'+str(p_value))
    #distance_matrix[phenotype_index_i][phenotype_index_j] = coefficient
    #distance_matrix_counter += 1
    #print('INFO: Completed matrix comparison '+str(distance_matrix_counter)+' out of '+str(distance_matrix_comparisons)+'.')
    #phenotype_i_draws = numpy.sum(ortholog_phenotype_matrix[phenotype_index_i])
    #print('Phenotype I draws = '+ str(phenotype_i_draws))
    #phenotype_j_draws = numpy.sum(ortholog_phenotype_matrix[phenotype_index_j])
    #print('Phenotype J draws = '+ str(phenotype_j_draws))
    for x in range(0, (len_ortholog_list)):
    #while ortholog_counter < len(ortholog_list):
        if read_only_ortholog_phenotype_matrix[phenotype_index_i][x] == 1 and read_only_ortholog_phenotype_matrix[phenotype_index_j][x] == 1:
            ortholog_match += 1
        ortholog_counter += 1
    #(ortholog_match, ortholog_counter) = map(ortholog_matches(ortholog_phenotype_matrix, phenotype_index_i, phenotype_index_j, x), x in range(0, len_ortholog_list))
    #print('Ortholog Matches: '+str(ortholog_match))
    #print('Ortholog Counter: '+str(ortholog_counter))

    # N = total number of orthologs shared between species
    # n = nummber of orthologs in species A phenotype
    # m = nummber of orthologs in species B phenotype
    # c = number of common orthologs between phenotypes (ortholog matches)

    m = float(numpy.sum(read_only_ortholog_phenotype_matrix[phenotype_index_i]))
    n = float(numpy.sum(read_only_ortholog_phenotype_matrix[phenotype_index_j]))
    N = float(len_ortholog_list) # Should this be the length of the ortholog list, or total orthologs shared between the three species?
    c = float(ortholog_match)
    hyp_prob = (hypergeom.cdf(c, N, m, n))
    #weight_matrix[phenotype_index_i][phenotype_index_j] = hyp_prob
    return (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient)

def ortholog_matches(ortholog_phenotype_matrix, phenotype_index_i, phenotype_index_j, x):
    ortholog_counter = 0
    ortholog_match = 0
    if ortholog_phenotype_matrix[phenotype_index_i][x] == 1 and ortholog_phenotype_matrix[phenotype_index_j][x] == 1:
        ortholog_match += 1
    ortholog_counter += 1
    return ortholog_match, ortholog_counter

def multiprocess_owlsim_queries(row):

    increment()
    (comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes) = row[0]
    try:
        response = urllib.request.urlopen(query_url, timeout=1000)
        reader = codecs.getreader("utf-8")
        data = json.load(reader(response))
        #print(data)
        #print('#####')
        #print('query successful')
        results = data['results']
        maxIC = data['results'][0]['maxIC']
        simJ = data['results'][0]['simJ']
        ICCS = data['results'][0]['bmaSymIC']
        simIC = data['results'][0]['simGIC']
        #print(results)
        query_flag = 'success'
        sequence = (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
        #lock.acquire()
        #json.dump(sequence, outfile)
        #outfile.write('\n')
        #lock.release()

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

        sequence = (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
        #lock.acquire()
        #json.dump(sequence, outfile)
        #outfile.write('\n')
        #lock.release()

    return (sequence)


class multithread_owlsim_queries(Thread):
    def __init__(self, url, name):
        Thread.__init__(self)
        self.name = name
        self.url = url
        return


    def run(self):

        print('INFO: Performing an OWLSim query.')
        (comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes) = tuple
        try:
            response = urllib.request.urlopen(query_url, timeout=1000)
            reader = codecs.getreader("utf-8")
            data = json.load(reader(response))
            #print(data)
            #print('#####')
            #print('query successful')
            results = data['results']
            maxIC = data['results'][0]['maxIC']
            simJ = data['results'][0]['simJ']
            ICCS = data['results'][0]['bmaSymIC']
            simIC = data['results'][0]['simGIC']
            #print(results)
            query_flag = 'success'
            sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
            #json.dump(sequence, outfile)
            #outfile.write('\n')

            #print(sequence)
            #print('failed here')
            #row = str.join(sequence)

            #print(row)
            #owlsimwriter.writerow(row)
            #print('query processing completed')

        except Exception:
            #print('Processing of OWLSim query failed.')
            #Creating an empty set of metrics for failed queries (queries with unresolved IRIs).
            maxIC = 0
            simJ = 0
            ICCS = 0
            simIC = 0
            query_flag = 'fail'

            sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
            #json.dump(sequence, outfile)
            #outfile.write('\n')



        return (sequence)

def multiprocess_fdr_calculation(i):

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

    print('INFO: Performing phenolog calculation on random data set '+str(i)+' out of 1000.')
    fdr_p_value_list = []
    hvm_human_file = hvm_human_dir+'random_'+str(i)+'.txt'
    hvm_mouse_file = hvm_mouse_dir+'random_'+str(i)+'.txt'
    hvz_human_file = hvz_human_dir+'random_'+str(i)+'.txt'
    hvz_zebrafish_file = hvz_zebrafish_dir+'random_'+str(i)+'.txt'
    mvz_mouse_file = mvz_mouse_dir+'random_'+str(i)+'.txt'
    mvz_zebrafish_file = mvz_zebrafish_dir+'random_'+str(i)+'.txt'

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
    hvm_phenolog_p_values = main.perform_phenolog_calculations_for_fdr(hvm_human_pheno_ortholog_hash, hvm_mouse_pheno_ortholog_hash, hvm_common_orthologs)
    #print(hvm_phenolog_p_values)
    print('INFO: Completed human vs mouse calculations for random data set '+str(i)+'.')
    fdr_p_value_list.extend(hvm_phenolog_p_values)
    hvz_phenolog_p_values = main.perform_phenolog_calculations_for_fdr(hvz_human_pheno_ortholog_hash, hvz_zebrafish_pheno_ortholog_hash, hvz_common_orthologs)
    #print(hvz_phenolog_p_values)
    print('INFO: Completed human vs zebrafish calculations for random data set '+str(i)+'.')
    fdr_p_value_list.extend(hvz_phenolog_p_values)
    mvz_phenolog_p_values = main.perform_phenolog_calculations_for_fdr(mvz_mouse_pheno_ortholog_hash, mvz_zebrafish_pheno_ortholog_hash, mvz_common_orthologs)
    #print(mvz_phenolog_p_values)
    fdr_p_value_list.extend(mvz_phenolog_p_values)
    print('INFO: Completed mouse vs zebrafish calculations for random data set '+str(i)+'.')

    # After grabbing the p-values from each function, assemble and sort.
    # Select the p-value that resides at the 0.05 percentile and add it to a list.

    #print('fdr p value list: '+str(len(fdr_p_value_list)))
    print('INFO: Sorting p-values for random data set '+str(i)+'.')
    fdr_p_value_list.sort()
    with open('inter/random/fdr/fdr_p_value_list_random_set_'+str(i)+'.txt', 'wb') as handle:
        pickle.dump(fdr_p_value_list, handle)
    #print(fdr_p_value_list)
    cutoff_position = math.ceil((len(fdr_p_value_list))*0.05) - 1
    #print(fdr_p_value_list[cutoff_position])
    fdr_cutoff_value = fdr_p_value_list[cutoff_position]

    print('INFO: Processing for random data set '+str(i)+' completed.')
    return fdr_cutoff_value

def multiprocess_generate_random_human_ext_data(x):

    random_geno_pheno_hash = {}

    phenotypes = []

    with open('inter/ontologies/hp_hash.txt', 'rb') as handle:
        phenotype_hash = pickle.load(handle)
    for i in phenotype_hash:
        if i not in phenotypes:
            phenotypes.append(i)
    with open('inter/hpo/human_disease_phenotype_hash.txt', 'rb') as handle:
        geno_pheno_hash = pickle.load(handle)

    random.shuffle(phenotypes)

    for i in geno_pheno_hash:
        random_geno_pheno_hash[i] = []
        ortholog_list_length = len(geno_pheno_hash[i])
        #print(ortholog_list_length)
        phenotype_draw = phenotypes
        random.shuffle(phenotype_draw)
        #random_orthologs = random.sample(orthologs)
        #test_pheno_ortholog_hash[i].append(random_orthologs)

        for j in geno_pheno_hash[i]:
            random.shuffle(phenotype_draw)
            random_geno_pheno_hash[i].append(phenotype_draw[0])


    with open(('inter/random/human/random_ext_'+str(x)+'.txt'), 'wb') as handle:
        pickle.dump(random_geno_pheno_hash, handle)
    print('Completed human random data set '+str(x)+' out of 1000.')
    return

def multiprocess_generate_random_mouse_ext_data(x):

    random_geno_pheno_hash = {}

    phenotypes = []

    with open('inter/ontologies/mp_hash.txt', 'rb') as handle:
        phenotype_hash = pickle.load(handle)
    for i in phenotype_hash:
        if i not in phenotypes:
            phenotypes.append(i)
    with open('inter/mgi/mouse_genotype_phenotype_hash.txt', 'rb') as handle:
        geno_pheno_hash = pickle.load(handle)

    random.shuffle(phenotypes)

    for i in geno_pheno_hash:
        random_geno_pheno_hash[i] = []
        ortholog_list_length = len(geno_pheno_hash[i])
        #print(ortholog_list_length)
        phenotype_draw = phenotypes
        random.shuffle(phenotype_draw)
        #random_orthologs = random.sample(orthologs)
        #test_pheno_ortholog_hash[i].append(random_orthologs)

        for j in geno_pheno_hash[i]:
            random.shuffle(phenotype_draw)
            random_geno_pheno_hash[i].append(phenotype_draw[0])


    with open(('inter/random/mouse/random_ext_'+str(x)+'.txt'), 'wb') as handle:
        pickle.dump(random_geno_pheno_hash, handle)
    print('Completed mouse random data set '+str(x)+' out of 1000.')
    return

def multiprocess_generate_random_zebrafish_ext_data(x):

    random_geno_pheno_hash = {}

    phenotypes = []

    with open('inter/ontologies/zp_hash.txt', 'rb') as handle:
        phenotype_hash = pickle.load(handle)
    for i in phenotype_hash:
        if i not in phenotypes:
            phenotypes.append(i)
    with open('inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'rb') as handle:
        geno_pheno_hash = pickle.load(handle)

    random.shuffle(phenotypes)

    for i in geno_pheno_hash:
        random_geno_pheno_hash[i] = []
        ortholog_list_length = len(geno_pheno_hash[i])
        #print(ortholog_list_length)
        phenotype_draw = phenotypes
        random.shuffle(phenotype_draw)
        #random_orthologs = random.sample(orthologs)
        #test_pheno_ortholog_hash[i].append(random_orthologs)

        for j in geno_pheno_hash[i]:
            random.shuffle(phenotype_draw)
            random_geno_pheno_hash[i].append(phenotype_draw[0])


    with open(('inter/random/zebrafish/random_ext_'+str(x)+'.txt'), 'wb') as handle:
        pickle.dump(random_geno_pheno_hash, handle)
    print('Completed zebrafish random data set '+str(x)+' out of 1000.')
    return

def multiprocess_ext_fdr_calculation(i):

    processing_start_time = time.time()
    #print('INFO: Setting stage for second FDR estimation.')
    # Need to calculate phenolog extension for each pairwise species and combine in order to get a full
    # set of 'genologs' (?) for proper estimation of FDR.

    human_dir = 'inter/random/human/'
    mouse_dir = 'inter/random/mouse/'
    zebrafish_dir = 'inter/random/zebrafish/'



    fdr_p_value_list = []
    #FIXME: Remove this on full data set running.
    #fdr_p_value_list.append(0.002222)
    human_file = human_dir+'random_ext_'+str(i)+'.txt'
    mouse_file = mouse_dir+'random_ext_'+str(i)+'.txt'
    zebrafish_file = zebrafish_dir+'random_ext_'+str(i)+'.txt'

    with open(human_file, 'rb') as handle:
        human_geno_pheno_hash = pickle.load(handle)
    with open(mouse_file, 'rb') as handle:
        mouse_geno_pheno_hash = pickle.load(handle)
    with open(zebrafish_file, 'rb') as handle:
        zebrafish_geno_pheno_hash = pickle.load(handle)
    with open('inter/phenolog/hvm_phenolog_combo.txt', 'rb') as handle:
        hvm_phenologs = pickle.load(handle)
        hvm_phenologs = set(hvm_phenologs)
        #print('INFO: There are '+str(len(hvm_phenologs))+' human-mouse phenologs present.')
    with open('inter/phenolog/hvz_phenolog_combo.txt', 'rb') as handle:
        hvz_phenologs = pickle.load(handle)
        hvz_phenologs = set(hvz_phenologs)
        #print('INFO: There are '+str(len(hvz_phenologs))+' human-zebrafish phenologs present.')
    with open('inter/phenolog/mvz_phenolog_combo.txt', 'rb') as handle:
        mvz_phenologs = pickle.load(handle)
        mvz_phenologs = set(mvz_phenologs)
        #print('INFO: There are '+str(len(mvz_phenologs))+' mouse-zebrafish phenologs present.')
    print('INFO: Performing phenolog extension calculation on random data set '+str(i)+' out of 1000.')

    #TODO: Need to return all of the p-values from the hypergeometric probability calculation for sorting and 5% cutoff.
    hvm_phenolog_p_values = main.perform_phenolog_calculations_for_ext_fdr(human_geno_pheno_hash, mouse_geno_pheno_hash, hvm_phenologs)
    #print(hvm_phenolog_p_values)
    fdr_p_value_list.extend(hvm_phenolog_p_values)
    print('INFO: Completed human vs mouse phenolog extension calculations on random data set '+str(i)+' out of 1000.')
    hvz_phenolog_p_values = main.perform_phenolog_calculations_for_ext_fdr(human_geno_pheno_hash, zebrafish_geno_pheno_hash, hvz_phenologs)
    #print(hvz_phenolog_p_values)
    fdr_p_value_list.extend(hvz_phenolog_p_values)
    print('INFO: Completed human vs zebrafish phenolog extension calculations on random data set '+str(i)+' out of 1000.')
    mvz_phenolog_p_values = main.perform_phenolog_calculations_for_ext_fdr(mouse_geno_pheno_hash, zebrafish_geno_pheno_hash, mvz_phenologs)
    #print(mvz_phenolog_p_values)
    fdr_p_value_list.extend(mvz_phenolog_p_values)
    print('INFO: Completed mouse vs zebrafish phenolog extension calculations on random data set '+str(i)+' out of 1000.')

    # After grabbing the p-values from each function, assemble and sort.
    # Select the p-value that resides at the 0.05 percentile and add it to a list.
    print('INFO: Sorting p-values for random data set '+str(i)+' out of 1000.')
    #print('fdr p value list: '+str(len(fdr_p_value_list)))
    fdr_p_value_list.sort()
    with open('inter/random/fdr_ext/fdr_ext_p_value_list_random_set_'+str(i)+'.txt', 'wb') as handle:
        pickle.dump(fdr_p_value_list, handle)
    #print(fdr_p_value_list)
    cutoff_position = math.ceil((len(fdr_p_value_list))*0.05) - 1
    #print(fdr_p_value_list[cutoff_position])
    if cutoff_position < 0:
        cutoff_position = 0
    fdr_cutoff_value = fdr_p_value_list[cutoff_position]

    print('INFO: Phenolog extension calculation on random data set '+str(i)+' completed.')
    processing_time = time.time() - processing_start_time
    print('Processing completed in '+str(processing_time)+' seconds.')

    return fdr_cutoff_value

def multiprocess_ext_fdr_calculation_alternate(comparison_list):
    #increment()

    total_phenotype_matches = 0
    total_phenotype_nonmatches = 0

    species_a_genotype_id = comparison_list[0]
    species_a_phenotypes = read_only_human_geno_pheno_hash[comparison_list[0]]
    #print(species_a_phenotypes)
    genotype_a_phenotype_count = len(species_a_phenotypes)

    # Genotype for species B
    species_b_genotype_id = comparison_list[1]
    species_b_phenotypes = read_only_zebrafish_geno_pheno_hash[comparison_list[1]]
    #print(species_b_phenotypes)
    phenotype_matches = 0
    phenotype_non_matches = 0


    genotype_b_phenotype_count = len(species_b_phenotypes)

    for k in species_a_phenotypes:
        # Orthologs for species A
        #ortholog_matches = 0
        #ortholog_non_matches = 0

        species_a_phenotype = k
        for l in species_b_phenotypes:
            # Orthologs for species B
            species_b_phenotype = l

            ab_combo = species_a_phenotype+'_'+species_b_phenotype
            ba_combo = species_b_phenotype+'_'+species_a_phenotype
            if ab_combo in read_only_hvz_phenologs or ba_combo in read_only_hvz_phenologs:
                #print('species a ortholog:'+species_a_ortholog+' matches species b ortholog:'+species_b_ortholog)
                phenotype_matches += 1
                #print(species_a_ortholog+' == '+species_b_ortholog)
                total_phenotype_matches += 1
            else:
                #print('species a ortholog:'+species_a_ortholog+' does not match species b ortholog:'+species_b_ortholog)
                phenotype_non_matches += 1
                total_phenotype_nonmatches += 1

    if phenotype_matches > 0:
        #print('Matches: '+str(ortholog_matches))
        #print('Non-matches: '+str(ortholog_non_matches))
        m = float(genotype_b_phenotype_count)
        n = float(genotype_a_phenotype_count)
        N = float(len(read_only_hvz_phenologs))
        c = float(phenotype_matches)
        prb = float(hypergeom.pmf(c, N, m, n))
        #print(str(c)+', '+str(N)+', '+str(m)+', '+str(n))
        #print(prb)
        #phenolog_ext_p_value_list.append(prb)
        #total_hyp_calcs += 1

        return prb
    else:
        return


####### MAIN #######

limit = None

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
#main.assemble_nif_zfin_phenotype_to_gene(limit)  # Completed in 3.22 days, 85118 rows processed.
#main.assemble_nif_mgi_phenotype_to_gene(limit)  # # Completed on full data set in 175.3 hours (7.3 days)
#main.assemble_nif_hpo_phenotype_to_gene(limit)  # Completed on full data set in 75.5 hours.
#main.assemble_nif_animalqtl_phenotype_to_gene(limit)

#main.assemble_nif_hpo_disease_to_gene(limit)
#main.assemble_nif_zfin_genotype_to_phenotype(limit)
#main.assemble_nif_mgi_genotype_to_phenotype(limit)
#main.assemble_nif_mgi_gene_to_phenotype(limit)
#main.assemble_nif_zfin_gene_to_phenotype(limit)
#main.assemble_nif_hpo_disease_to_phenotype(limit)


####### ASSEMBLE OWLSIM COMPARISON QUERIES #######
### CAUTION: These assemblies will take a large amount of time and file space!
#main.assemble_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/owlsim/human_disease_mouse_genotype', 'human_disease_mouse_genotype_queries', limit)
#main.assemble_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'inter/owlsim/human_disease_zebrafish_genotype', 'human_disease_zebrafish_genotype_queries', limit)
#main.assemble_owlsim_queries('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'inter/owlsim/mouse_genotype_zebrafish_genotype', 'mouse_genotype_zebrafish_genotype_queries', limit)
#main.assemble_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/owlsim/human_disease_mouse_gene', 'human_disease_mouse_gene_queries', limit)

#main.assemble_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'inter/owlsim/human_disease_zebrafish_gene', 'human_disease_zebrafish_gene_queries', limit)

#main.assemble_owlsim_queries('inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'inter/owlsim/mouse_gene_zebrafish_gene', 'mouse_gene_zebrafish_gene_queries', limit)



####### OWLSIM COMPARISONS #######

#OWLSim url calls take about 3 hours for 100,000 comparisons.
#Current implementation runs at 185/second, or 666,000/hour?
# 42 million would take 63 hours, or 2.5 days

#latest test: 320/second

#Processing completed in  hours,  comparisons.
#Human Diseases = 9214
#Mouse Genotypes = 56427
#Total comparisons = 519,918,378
# Compare human disease phenotypic profiles & mouse genotype phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human diseases vs mouse genotypes')
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/owlsim/human_disease_mouse_genotype', 'human_disease_mouse_genotype_queries', 'out/owlsim/human_disease_mouse_genotype', 'human_disease_mouse_genotype_results', 104)
#print('INFO: Done processing human diseases vs mouse genotypes')

#Processing completed in  hours,  comparisons.
#Human Diseases = 9214
#zebrafish Genotype = 8535
#Total comparisons = 78,641,490
# Compare human disease phenotypic profiles & zebrafish genotype phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human disease vs zebrafish genotype')
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'inter/owlsim/human_disease_zebrafish_genotype', 'human_disease_zebrafish_genotype_queries', 'out/owlsim/human_disease_zebrafish_genotype', 'human_disease_zebrafish_genotype_results', 16)
#print('INFO: Done processing human disease vs zebrafish genotype')

#Processing completed in  hours,  comparisons. Estimated to take 669 days?
#Mouse genotype = 56427
#zebrafish genotype = 8535
#Total comparisons = 481,604,445
# Compare mouse genotype phenotypic profiles & zebrafish genotype phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing mouse genotype vs zebrafish genotypes')
#main.perform_owlsim_queries('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'inter/owlsim/mouse_genotype_zebrafish_genotype', 'mouse_genotype_zebrafish_genotype_queries', 'out/owlsim/mouse_genotype_zebrafish_genotype', 'mouse_genotype_zebrafish_genotype_results', 97)
#print('INFO: Done processing mouse genotype vs zebrafish genotypes')

#Processing completed!
#Human Diseases = 9214
#Mouse genes = 13102
#Total comparisons = 120,712,614
# Compare human disease phenotypic profiles & mouse gene phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human disease vs mouse genes')
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/owlsim/human_disease_mouse_gene','human_disease_mouse_gene_queries', 'out/owlsim/human_disease_mouse_gene', 'human_disease_mouse_gene_results', 25)
#print('INFO: Done processing human disease vs mouse genes')

#Processing completed!
#Human Diseases = 9214
#zebrafish Genes = 4580
#Total comparisons = 42,190,906
# Compare human disease phenotypic profiles & zebrafish gene phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing human disease vs zebrafish genes')
#(self, raw1, raw2, interfile_directory, interfile_prefix, outfile_directory, , outfile_prefix, num_files, limit=None)
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'inter/owlsim/human_disease_zebrafish_gene', 'human_disease_zebrafish_gene_queries', 'out/owlsim/human_disease_zebrafish_gene', 'human_disease_zebrafish_gene_results', 9)
#main.perform_owlsim_queries_threaded('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt','out/owlsim/human_disease_zebrafish_gene.txt')
#print('INFO: Done processing human disease vs zebrafish genes')

#Processing completed!
#Mouse Genes = 13102
#zebrafish Genes = 4580
#Total comparisons = 59,989,479
# Compare mouse gene phenotypic profiles & zebrafish gene phenotypic profiles via OWLSim.
#print('INFO: OWLSim processing mouse genes vs zebrafish genes')
#main.perform_owlsim_queries('inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'inter/owlsim/mouse_gene_zebrafish_gene','mouse_gene_zebrafish_gene_queries', 'out/owlsim/mouse_gene_zebrafish_gene', 'mouse_gene_zebrafish_gene_results', 12)
#print('INFO: Done processing mouse genes vs zebrafish genes')

#with open('inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'rb') as handle:
    #species_b_pheno_gene_hash = pickle.load(handle)
#print(len(species_b_pheno_gene_hash.keys()))

#TODO: Data assembly
# Would it be easier/more efficient to add the phenolog data to the OWLSim output table,
# or do separate output tables and then match them together into a combined table? Essentially need to do an SQL join.


####### PHENOLOG FDR CALCULATION #######

# Each random data generation takes ~8-10 hours.
#main.generate_random_data('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_mouse.txt', 'inter/random/human_vs_mouse/mouse/')
#main.generate_random_data('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_mouse.txt', 'inter/random/human_vs_mouse/human/')

#main.generate_random_data('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', 'inter/random/human_vs_zebrafish/zebrafish/')
#main.generate_random_data('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', 'inter/random/human_vs_zebrafish/human/')

#main.generate_random_data('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', 'inter/random/mouse_vs_zebrafish/zebrafish/')
#main.generate_random_data('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', 'inter/random/mouse_vs_zebrafish/mouse/')
#print('INFO: Done with random data generation.')

#main.set_stage_for_fdr_calculation()
#fdr_cutoff = main.set_stage_for_fdr_calculation()
#main.assemble_partial_fdr()
#print(fdr_cutoff)

####### PHENOLOG COMPARISONS #######
# NOTE: Either run the FDR calculations or set an FDR cutoff before running the phenolog calculations.
# Cutoff below is the average from the 1000 random data sets.
fdr_cutoff = 0.004426898733810069
#main.perform_phenolog_calculations('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt', 'out/phenolog/human_vs_mouse.txt', 'inter/panther/common_orthologs_human_mouse.txt', fdr_cutoff)
#main.perform_phenolog_calculations('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'out/phenolog/mouse_vs_zebrafish.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', fdr_cutoff)
#main.perform_phenolog_calculations('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'out/phenolog/human_vs_zebrafish.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', fdr_cutoff)

#main.assemble_significant_phenologs()
#main.assemble_significant_phenologs_with_scores()

####### PHENOLOG EXTENSION FDR CALCULATION #######

#main.parse_hp('raw/ontologies/hp.obo', 'inter/ontologies/hp_hash.txt')
#main.parse_mp('raw/ontologies/MPheno_OBO.ontology', 'inter/ontologies/mp_hash.txt')
#main.parse_zp('raw/ontologies/zp_mapping.txt', 'inter/ontologies/zp_hash.txt')

#Creation of random data sets for phenlog extension FDR.
#main.generate_human_random_ext_data()
#main.generate_mouse_random_ext_data()
#main.generate_zebrafish_random_ext_data()

#Assemble the phenlog lookup files.
#main.assemble_hvm_phenologs()
#main.assemble_hvz_phenologs()
#main.assemble_mvz_phenologs()


#main.set_stage_for_extension_fdr_calculation()

for i in range(212, 220):

    with open('inter/phenolog/hvz_phenolog_combo.txt', 'rb') as handle:
        read_only_hvz_phenologs = set(pickle.load(handle))
    with open('inter/random/human/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_human_geno_pheno_hash = pickle.load(handle)
    with open('inter/random/zebrafish/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_zebrafish_geno_pheno_hash = pickle.load(handle)
    print('INFO: Processing human vs zebrafish random data set '+str(i)+'.')
    p_value_out_file = 'inter/phenolog_ext/hvz_p_values/hvz_p_values_'+str(i)+'.txt'
    main.perform_phenolog_calculations_for_ext_fdr_alternate(read_only_human_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash, p_value_out_file)

    read_only_hvz_phenologs = []
    read_only_human_geno_pheno_hash = {}
    read_only_zebrafish_geno_pheno_hash = {}
    gc.collect()

    print('INFO: Done processing human vs zebrafish random data set '+str(i)+'.')

'''
for i in range(1, 1001):

    with open('inter/phenolog/hvm_phenolog_combo.txt', 'rb') as handle:
        read_only_hvm_phenologs = set(pickle.load(handle))
    with open('inter/random/human/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_human_geno_pheno_hash = pickle.load(handle)
    with open('inter/random/mouse/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_mouse_geno_pheno_hash = pickle.load(handle)
    print('INFO: Processing human vs mouse random data set '+str(i)+'.')
    p_value_out_file = 'inter/phenolog_ext/hvm_p_values/hvm_p_values_'+str(i)+'.txt'
    main.perform_phenolog_calculations_for_ext_fdr_alternate(read_only_human_geno_pheno_hash, read_only_mouse_geno_pheno_hash, p_value_out_file)

    read_only_phenologs = []
    read_only_human_geno_pheno_hash = {}
    read_only_mouse_geno_pheno_hash = {}

    print('INFO: Done processing human vs zebrafish random data set '+str(i)+'.')

for i in range(1, 1001):

    with open('inter/phenolog/mvz_phenolog_combo.txt', 'rb') as handle:
        read_only_mvz_phenologs = set(pickle.load(handle))
    with open('inter/random/mouse/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_mouse_geno_pheno_hash = pickle.load(handle)
    with open('inter/random/zebrafish/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_zebrafish_geno_pheno_hash = pickle.load(handle)
    print('INFO: Processing mouse vs zebrafish random data set '+str(i)+'.')
    p_value_out_file = 'inter/phenolog_ext/mvz_p_values/hvz_p_values_'+str(i)+'.txt'
    main.perform_phenolog_calculations_for_ext_fdr_alternate(read_only_mouse_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash, p_value_out_file)

    read_only_mvz_phenologs = []
    read_only_mouse_geno_pheno_hash = {}
    read_only_zebrafish_geno_pheno_hash = {}

    print('INFO: Done processing human vs zebrafish random data set '+str(i)+'.')
'''
#main.perform_hvm_phenolog_calculations_for_ext_fdr_alternate(read_only_human_geno_pheno_hash, read_only_mouse_geno_pheno_hash)
#main.perform_mvz_phenolog_calculations_for_ext_fdr_alternate(read_only_mouse_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash)



#ext_fdr_cutoff = 0.00022089684117479534
#main.perform_phenolog_ext_calculations('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt', 'out/phenolog_ext/human_vs_mouse.txt', 'inter/phenolog/hvm_significant_phenologs.txt', ext_fdr_cutoff)
#main.perform_phenolog_ext_calculations('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt', 'out/phenolog_ext/human_vs_mouse.txt', 'inter/phenolog/all_significant_phenologs.txt', ext_fdr_cutoff)
#main.perform_phenolog_ext_calculations('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'out/phenolog_ext/human_vs_zebrafish.txt', 'inter/phenolog/hvz_significant_phenologs.txt', ext_fdr_cutoff)
#main.perform_phenolog_ext_calculations('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'out/phenolog_ext/mouse_vs_zebrafish.txt', 'inter/phenolog/mvz_significant_phenologs.txt', ext_fdr_cutoff)


####### PHENOLOG GENE CANDIDATE PREDICTIONS #######

#NOTE: Must be assembled after the nif phenotype to gene assembly has been performed.
#NOTE: These three functions are now combined and run in the create_phenolog_gene_candidate_matrices function.
#(human_ortholog_phenotype_matrix, human_phenotype_list, human_ortholog_list) = main.assemble_ortholog_phenotype_matrix('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/hpo/human_pheno_ortholog_matrix.txt')
#(mouse_ortholog_phenotype_matrix, mouse_phenotype_list, mouse_ortholog_list) = main.assemble_ortholog_phenotype_matrix('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_matrix.txt')
#(zebrafish_ortholog_phenotype_matrix, zebrafish_phenotype_list, zebrafish_ortholog_list) = main.assemble_ortholog_phenotype_matrix('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_matrix.txt')

#This process requires multi-processing due to the large number of comparisons that need to be performed.
#cProfile.run()
#main.assemble_ortholog_phenotype_matrix()
#main.assemble_ortholog_phenotype_matrix_alternate(100)
#main.create_phenolog_gene_candidate_matrices()
#main.create_empty_phenolog_gene_candidate_matrices()
#main.populate_phenolog_gene_candidate_matrices()
#main.populate_phenolog_gene_candidate_matrices_alternate()

#read_only_ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')
#main.populate_phenolog_gene_candidate_matrices_alternate()

#main.populate_phenolog_gene_candidate_matrices_alternate()
#read_only_ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')
#main.populate_phenolog_gene_candidate_matrices_alternate()
#main.merge_matrices()
#main.create_phenolog_gene_candidate_prediction_matrix()
#main.assemble_phenolog_gene_candidate_predictions()
#test_matrix = numpy.zeros((5, 2))
#print(test_matrix)
#test_matrix[0][0] = 1
#test_matrix[3][1] = 5
#print(test_matrix)

#with open('inter/owlsim/mouse_gene_zebrafish_gene_queries.txt', 'r', encoding="iso-8859-1") as csvfile:
    #filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
    #row_count = sum(1 for row in filereader)
    #print(str(row_count)+' rows to process.')

'''
a = ['A', 'B', 'C']
b = ['X', 'Y', 'Z']
for element in itertools.product(a,b):
    print(element)
    print(element[0])
    print(element[1])
'''


elapsed_time = time.time() - start_time
print('Processing completed in '+str(elapsed_time)+' seconds.')

#TODO: Make sure and have the ability to filter between single-gene genotypes and multi-gene genotypes.

#http://owlsim.monarchinitiative.org/compareAttributeSets?a=HP:0001263&b=MP:0010864
#Format for mutliple:
#http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878


