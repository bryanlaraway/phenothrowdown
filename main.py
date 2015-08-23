#!/usr/bin/env python3

#AUTHOR: Bryan Laraway
#PROJECT: Comparison of OWLSim and Phenologs for the identification of models of human disease and gene candidates for human disease.
#PURPOSE: This script will call of the functions/methods/scripts for performing the processing required for this analysis.

import json
import urllib.request
import codecs
import time
import gc
#from socket import *
#import os
import sys
import re
import csv
import pickle
from decimal import Decimal, getcontext
#import cProfile
import numpy
from numpy import random
from scipy.stats import hypergeom, pearsonr
import math
import heapq
import multiprocessing
#from memory_profiler import profile
#from multiprocessing import Pool, Process, Manager, Lock, Value
import itertools
from threading import Thread
#import threading
from ctypes import c_int
from queue import Queue
from collections import *
from functools import reduce
import matplotlib.pyplot as plt


start_time = time.time()

hu_disease_to_phenotype_hash = {'disease_id': {}}
mouse_genotype_to_phenotype_hash = {'genotype_id': {}}
zfin_genotype_to_phenotype_hash = {'genotype_id': {}}
getcontext().prec = 500
#print(getcontext())
#Selected distinct PANTHER IDs from the NIF/DISCO tables.
#TODO: See about getting these numbers from the Panther table to allow for dynamic updating with file updates.
#total_human_mouse_orthologs = 5625, with LDO = 5729
#total_human_zebrafish_orthologs = 5212, with LDO = 5750
#total_mouse_zebrafish_orthologs = 5210, with LDO = 5748
#TODO: Need to pass phenotype/gene labels for identification, or look them up upon final output.


class main():

    # NOTE: Could either include the fetch code to retrieve the data from the resources,
    # or retrieve them and have the code just open local files, already retrieved.

    # Required table from NIF/DISCO
    tables = [
        'dvp.pr_nlx_151835_1',  # HPO: Annotations:DiseasePhenotypes view
        'dvp.pr_nlx_151835_2',  # HPO: Annotations:Phenotype to gene view
        'dvp.pr_nlx_151835_3',  # HPO: Annotations:Disease to gene
        'dvp.pr_nif_0000_00096_5',  # MGI:MouseGenotypes
        'dvp.pr_nif_0000_00096_6',  # MGI:MousePhenotypes
        'dvp.pr_nif_0000_21427_10',  # ZFIN:Genotype-Phenotype
        'dvp.pr_nif_0000_21427_11',  # ZFIN:OrganismGenotypes
        'dvp.pr_nlx_84521_1'  # PANTHER:Orthologs,
        'dvp.pr_nif_0000_02550_3' # ANIMAL QTL DB: Traits
    ]

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

    ####### PHENOTYPE ID TO LABEL ASSEMBLY #######

    def assemble_nif_hpo_phenotype_id_to_label(self, limit=None):
        """This function assembles a hash for human phenotype IDs and their labels from the NIF/DISCO flat data file"""

        print('INFO:Assembling human phenotype ID to label hash.')

        # Set up counters and open required files.
        line_counter = 0
        raw = 'raw/hpo/dvp.pr_nlx_151835_2'
        inter = 'inter/hpo/human_phenotype_id_to_label_hash.txt'
        hpo_phenotype_id_to_label_hash = {}
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
                if phenotype_id not in hpo_phenotype_id_to_label_hash:
                    hpo_phenotype_id_to_label_hash[phenotype_id] = phenotype_label

                if limit is not None and line_counter > limit:
                    break

        # Dump files to disk.
        with open(inter, 'wb') as handle:
            pickle.dump(hpo_phenotype_id_to_label_hash, handle)


        print('INFO: Done assembling human phenotype ID to label hash.')
        print('INFO: '+str(len(hpo_phenotype_id_to_label_hash.keys()))+' human phenotypes present.')
        return

    def assemble_nif_zfin_phenotype_id_to_label(self, limit=None):
        """This function assembles a hash for zebrafish phenotype IDs and their labels from the NIF/DISCO flat data file"""

        print('INFO:Assembling zebrafish phenotype to ortholog data.')

        # Set up counters and open required files.
        line_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        inter = 'inter/zfin/zebrafish_phenotype_id_to_label_hash.txt'
        zfin_phenotype_id_to_label_hash = {}
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

                # If phenotype is not in the phenotype to gene hash, add phenotype to hash.
                if phenotype_id not in zfin_phenotype_id_to_label_hash:
                    zfin_phenotype_id_to_label_hash[phenotype_id] = phenotype_label

                if limit is not None and line_counter > limit:
                    break

        # Dump data to files.
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_phenotype_id_to_label_hash, handle)

        print('INFO: Done assembling zebrafish phenotype to gene/ortholog data.')
        print('INFO: '+str(len(zfin_phenotype_id_to_label_hash.keys()))+' zebrafish phenotypes present.')
        return

    def assemble_nif_mgi_phenotype_id_to_label(self, limit=None):
        """This function assembles a hash for mouse phenotype IDs and their labels from the NIF/DISCO flat data file"""

        print('INFO:Assembling mouse phenotype to gene/ortholog data.')

        # Set up counters and open required files.
        line_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        inter = 'inter/mgi/mouse_phenotype_id_to_label_hash.txt'
        mgi_phenotype_id_to_label_hash = {}
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
                if phenotype_id not in mgi_phenotype_id_to_label_hash:
                    mgi_phenotype_id_to_label_hash[phenotype_id] = phenotype_label

                if limit is not None and line_counter > limit:
                    break

        # Dump data to files.
        with open(inter, 'wb') as handle:
            pickle.dump(mgi_phenotype_id_to_label_hash, handle)

        print('INFO: Done assembling mouse phenotype to gene/ortholog data.')
        print('INFO: '+str(len(mgi_phenotype_id_to_label_hash.keys()))+' mouse phenotypes present.')

        return


    ####### GENE ID TO LABEL ASSEMBLY #######

    def assemble_nif_mgi_gene_id_to_label(self):
        """This function assembles mouse gene id to gene labels from the NIF/DISCO flat data file"""

        print('INFO:Assembling mouse gene ID to label hash.')

        # Set up counters and open required files.
        line_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        inter = 'inter/mgi/mouse_gene_id_to_label_hash.txt'
        gene_id_to_label_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' mouse rows to process.')
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
                #print('INFO: Processing phenotype '+str(line_counter)+' out of '+str(row_count)+'.')

                if implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue

                # Split the implicated genes list.
                genes = implicated_gene_ids.split(',')
                gene_labels = implicated_gene_labels.split(',')
                #print(genes)
                # If gene is not in the gene ID to label hash, add gene to hash.
                for gene in genes:
                    if gene not in gene_id_to_label_hash:
                        gene_index = genes.index(gene)
                        #print(gene_index)
                        gene_id_to_label_hash[gene] = gene_labels[gene_index]


        # Dump data to files.
        with open(inter, 'wb') as handle:
            pickle.dump(gene_id_to_label_hash, handle)

        print('INFO: Done assembling mouse gene ID to label hash.')

        return

    def assemble_nif_zfin_gene_id_to_label(self):
        print('INFO:Assembling zebrafish gene ID to label hash.')
        line_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        inter = 'inter/zfin/zebrafish_gene_id_to_label_hash.txt'
        gene_id_to_label_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish gene to phenotype rows to process.')
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

                # Skip phenotypes without IDs and phenotypes with no associated genes.
                if phenotype_id == '' or phenotype_id is None:
                    continue
                if implicated_gene_ids == '' or implicated_gene_ids is None:
                    continue

                # Split the implicated genes list.
                genes = implicated_gene_ids.split(',')

                gene_labels = implicated_gene_labels.split(',')
                #print(genes)
                # If gene is not in the gene ID to label hash, add gene to hash.
                for gene in genes:
                    if gene not in gene_id_to_label_hash:
                        gene_index = genes.index(gene)
                        #print(gene_index)
                        gene_id_to_label_hash[gene] = gene_labels[gene_index]

        # Dump data to files.
        with open(inter, 'wb') as handle:
            pickle.dump(gene_id_to_label_hash, handle)

        print('INFO: Done assembling zebrafish gene to label hash.')

        return

    def assemble_nif_hpo_gene_id_to_label(self):
        print('INFO:Assembling human gene ID to label hash.')
        line_counter = 0
        failure_counter = 0
        raw1 = 'raw/hpo/dvp.pr_nlx_151835_3'
        raw2 = 'raw/hpo/dvp.pr_nlx_151835_2'
        inter = 'inter/hpo/human_gene_id_to_label_hash.txt'
        gene_id_to_label_hash = {}
        with open(raw1, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human disease to gene rows to process.')
        with open(raw1, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, disease_id, disorder_name, disorder_database_link, gene_id,
                 gene_num, gene_label, v_uid, v_uuid, v_lastmodified) = row
                #print(disease_id)

                # Convert NCBIGene ID prefix.
                gene_id = re.sub('NCBI_gene:', 'NCBIGene:', gene_id)
                #print(genes)
                if gene_id not in gene_id_to_label_hash:
                    gene_id_to_label_hash[gene_id] = gene_label
                    #print(hpo_phenotype_to_gene_hash[genotype_id])

        # Set up counters and open required files.
        line_counter = 0
        with open(raw2, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human phenotype rows to process.')
        with open(raw2, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:

                # Read in a row and split into individual variables
                line_counter += 1
                (e_uid, phenotype_id, phenotype_label, gene_id, gene_num,
                 gene_label, v_uid, v_uuid, v_lastmodified) = row

                if phenotype_id == '' or phenotype_id is None:
                    continue
                if gene_id == '' or gene_id is None:
                    continue

                # Convert NCBIGene ID prefix.
                gene_id = re.sub('NCBI_gene:', 'NCBIGene:', gene_id)
                if gene_id not in gene_id_to_label_hash:
                    gene_id_to_label_hash[gene_id] = gene_label

        # Dump files to disk.
        with open(inter, 'wb') as handle:
            pickle.dump(gene_id_to_label_hash, handle)

        print('INFO: Done assembling human gene ID to label hash.')
        return


    ####### PHENOLOG PHENOTYPE TO GENE/ORTHOLOG #######

    def trim_panther_data(self, inter, taxons):
        """
        This function trims the PANTHER flat file from NIF/DISCO for a given taxon,
        which speeds up data assembly when converting from genes to orthologs.
        :param inter: Directory and file name for saving the trimmed PANTHER file.
        :param taxons: taxon IDs for filtering.
        :return:
        """

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
                    line_counter += 1
                    (panther_speciesa, tax_id_a, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,
                    proteina, panther_speciesb, tax_id_b, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,
                    gene_label_b, proteinb, orthology_class, orthology_class_label, ancestor_taxon, panther_id,
                    e_uid, v_uid, v_uuid, v_lastmodified) = row

                    #Currently filtering on the big three taxons, and ortholog relations only.
                    if (taxon_id_a in taxons or taxon_id_b in taxons) and (orthology_class_label == 'Least Diverged Ortholog' or orthology_class_label == 'Ortholog'):
                        output_row = (panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,
                        proteina, panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,
                        gene_label_b, proteinb, orthology_class, orthology_class_label, panther_id)
                        output_line_counter += 1
                        csvwriter.writerow(output_row)

        print('PANTHER file trimmed to '+str(output_line_counter)+' rows.')
        return

    def get_common_orthologs(self, inter, taxons):
        """
        This function takes as input a list of taxons and a the PANTHER flat file,
        filters the PANTHER flat file using the provided taxon IDs to obtain
        the common orthologs between the taxons, and writes the list to a new output file.
        :param inter: directory & file name for the output file.
        :param taxons: taxon IDs for filtering the PANTHER flat file.
        :return:
        """
        print('INFO: Getting common orthologs between species.')
        line_counter = 0
        ortholog_counter = 0
        raw = 'raw/panther/dvp.pr_nlx_84521_1'
        common_orthologs = []

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

                if (taxon_id_a in taxons or taxon_id_b in taxons) and (orthology_class_label == 'Least Diverged Ortholog' or orthology_class_label == 'Ortholog'):
                    if panther_id not in common_orthologs:
                        common_orthologs.append(panther_id)
                        ortholog_counter += 1

        with open(inter, 'wb') as handle:
            pickle.dump(common_orthologs, handle)

        print(str(ortholog_counter)+' common orthologs found.')

        return

    def get_ortholog(self, query_gene_id, panther):
        """
        This function is used when creating the phenotype-ortholog hashes.
        It takes a gene ID and, if there is an ortholog match, return the PANTHER ID of the ortholog.
        :param query_gene_id: Gene ID used to query for an ortholog.
        :param panther: PANTHER file trimmed for the specific taxon.
        :return: PANTHER ID if successfail, fail flag if unsuccessful.
        """
        with open(panther, 'r') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            print('query= '+query_gene_id)
            for row in filereader:
                (panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,proteina,
                 panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,gene_label_b, proteinb,
                 orthology_class, orthology_class_label, panther_id) = row
                if query_gene_id in [genea, gene_id_a, geneb, gene_id_b]:
                    result_panther_id = panther_id
                    print('found ortholog for '+query_gene_id+'.')
                    return(result_panther_id)
        print('no ortholog found for '+query_gene_id+'.')
        return('fail')

    def assemble_nif_zfin_phenotype_to_gene(self, limit=None):
        """This function assembles zebrafish phenotype to gene associations from the NIF/DISCO flat data file"""

        print('INFO:Assembling zebrafish phenotype to ortholog data.')
        gene_to_ortholog_hash = {}
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
                        if gene not in gene_to_ortholog_hash:

                            # Convert genes to orthologs using zebrafish-trimmed PANTHER table as lookup.
                            panther_id = self.get_ortholog(gene, 'inter/panther/panther_zebrafish.txt')
                            gene_to_ortholog_hash[gene] = panther_id

                            # If ortholog is not in the phenotype to ortholog hash, add ortholog to hash.
                            if panther_id != 'fail' and panther_id not in zfin_phenotype_to_ortholog_hash[phenotype_id]:
                                zfin_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)
                        else:
                            panther_id = gene_to_ortholog_hash[gene]
                            if panther_id != 'fail' and panther_id not in zfin_phenotype_to_ortholog_hash[phenotype_id]:
                                zfin_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)


                if limit is not None and line_counter > limit:
                    break

        # Dump data to files.
        with open(inter1, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_ortholog_hash, handle)
        with open('inter/zfin/zebrafish_gene_to_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(gene_to_ortholog_hash, handle)

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
        gene_to_ortholog_hash = {}
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

                        if gene not in gene_to_ortholog_hash:

                            # Convert genes to orthologs using mouse-trimmed PANTHER table as lookup.
                            panther_id = self.get_ortholog(gene,'inter/panther/panther_mouse.txt')
                            gene_to_ortholog_hash[gene] = panther_id

                            # If ortholog is not in the phenotype to ortholog hash, add ortholog to hash.
                            if panther_id != 'fail' and panther_id not in mgi_phenotype_to_ortholog_hash[phenotype_id]:
                                mgi_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)
                        else:
                            panther_id = gene_to_ortholog_hash[gene]
                            if panther_id != 'fail' and panther_id not in mgi_phenotype_to_ortholog_hash[phenotype_id]:
                                mgi_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break

        # Dump data to files.
        with open(inter1, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_ortholog_hash, handle)
        with open('inter/mgi/mouse_gene_to_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(gene_to_ortholog_hash, handle)
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
        gene_to_ortholog_hash = {}
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

                    if gene_id not in gene_to_ortholog_hash:

                        # Convert genes to orthologs using human-trimmed PANTHER table as lookup.
                        panther_id = self.get_ortholog(gene_id,'inter/panther/panther_human.txt')
                        gene_to_ortholog_hash[gene_id] = panther_id

                        # If ortholog is not in the phenotype to ortholog hash, add ortholog to hash.
                        if panther_id != 'fail' and panther_id not in hpo_phenotype_to_ortholog_hash[phenotype_id]:
                            hpo_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)
                    else:
                        panther_id = gene_to_ortholog_hash[gene_id]
                        if panther_id != 'fail' and panther_id not in hpo_phenotype_to_ortholog_hash[phenotype_id]:
                            hpo_phenotype_to_ortholog_hash[phenotype_id].append(panther_id)

                if limit is not None and line_counter > limit:
                    break

        # Dump files to disk.
        with open(inter1, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_gene_hash, handle)
        with open(inter2, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_ortholog_hash, handle)
        with open('inter/hpo/human_gene_to_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(gene_to_ortholog_hash, handle)

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


    ####### OWLSIM GENOTYPE/GENE TO PHENOTYPE #######

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


                if extrinsic_genotype_id != '':
                    print('Skipping genotype with extrinsic modifiers: '+effective_genotype_id)
                    continue
                elif extrinsic_genotype_id == '' or extrinsic_genotype_id is None:
                    if phenotype_id != '' and phenotype_id is not None:
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

                if effective_genotype_id not in mgi_genotype_to_phenotype_hash and phenotype_id != '':
                    mgi_genotype_to_phenotype_hash[effective_genotype_id] = [phenotype_id]
                elif phenotype_id != '':
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

                # Convert NCBIGene ID prefix.
                gene_id = re.sub('NCBI_gene:', 'NCBIGene:', gene_id)

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
                if disorder_id not in hpo_disease_to_phenotype_hash and phenotype_id != '':
                    hpo_disease_to_phenotype_hash[disorder_id] = [phenotype_id]
                    #print(hpo_disease_to_phenotype_hash[disease)id])
                elif phenotype_id != '':
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

                if phenotype_id == '' or phenotype_id is None:
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


    ####### OWLSIM DATA PROCESSING ####### DOCUMENTATION COMPLETED

    def perform_owlsim_queries_threaded(self, raw1, raw2, out, limit=None):
        """
        This is an unused testing function for multi-threading of OWLSim queries. Multiprocessing is sufficiently fast.
        :param raw1:
        :param raw2:
        :param out:
        :param limit:
        :return:
        """

        print('INFO: Performing OWLSim queries.')
        line_counter = 0
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
        base_url = 'http://0.0.0.0:9031/compareAttributeSets?'
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
                #for p in results:
                    #sequence  = p.get()
                    #sequence = (entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                    #json.dump(sequence, outfile)
                    #outfile.write('\n')
                print('Done processing results.')
                print('INFO: Multiprocessing completed')
                ###### END THREADING INSERT ######

        return

    def assemble_owlsim_queries(self, raw1, raw2, interfile_directory, interfile_prefix, limit=None):
        """
        This function assembles the comparison ID, disease/genotype/gene IDs, and URL query for the OWLSim server,
        and saves each set as a row in a text file.
        :param raw1:
        :param raw2:
        :param interfile_directory:
        :param interfile_prefix:
        :param limit:
        :return:
        """

        # Set up line counters and max line count for splitting to separate files.
        line_counter = 0
        line_counter_max = 5000000
        file_counter = 1

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

        # Example URL format: http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878
        # Local server URL format: http://0.0.0.0:9031/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878
        base_url = 'http://0.0.0.0:9031/compareAttributeSets?'
        print('INFO: Assembling phenotypic profile comparison queries.')
        file_name = interfile_directory+'/'+interfile_prefix+'_'+str(file_counter)+'.txt'
        csvfile = open(file_name, 'w', newline='')
        csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')

        # Cycle through each organism's disease/genotype to phenotype hash and
        # assemble the comparison information and OWLSim query for
        # each phenotypic profile comparison, then write to file.
        for i in organism_a_hash:
            entity_a = i
            entity_a_attributes = organism_a_hash[i]
            phenotypic_profile_a = 'a='+('&a=').join(entity_a_attributes)
            for j in organism_b_hash:
                entity_b = j
                entity_b_attributes = organism_b_hash[j]
                phenotypic_profile_b = '&b='+('&b=').join(entity_b_attributes)
                query_url = base_url+phenotypic_profile_a+phenotypic_profile_b
                comparison_id = entity_a+'_'+entity_b
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

    def perform_owlsim_queries(self, raw1, raw2, interfile_directory, interfile_prefix, outfile_directory, outfile_prefix, num_files):
        """
        This function takes as input the assembled OWLSim query files, one row at a time, and performs the query.
        It then parses the results from the OWLSim server, adding a success/fail flag, and writes to a file.
        :param raw1: Used to obtain the number of diseases/genotypes/genes to compare from organism 1.
        :param raw2: Used to obtain the number of diseases/genotypes/genes to compare from organism .
        :param interfile_directory: Directory of the OWLSim query text files.
        :param interfile_prefix: Prefix of the OWLSim query text files.
        :param outfile_directory: Directory of the OWLSim output text files.
        :param outfile_prefix: Prefix of the OWLSim output text files.
        :param num_files: Total number of OWLSim query text files.
        :param limit: Optional limit for testing, not currently used
        :return:
        """

        print('INFO: Performing OWLSim queries.')
        with open(raw1, 'rb') as data1:
            organism_a_hash = pickle.load(data1)
        with open(raw2, 'rb') as data2:
            organism_b_hash = pickle.load(data2)
        comparison_count = len(organism_a_hash) * len(organism_b_hash)
        print('INFO: '+str(comparison_count)+' phenotypic profile comparisons to process.')

        # Cycle through the OWLSim query files to feed queries to the OWLSim server.
        for x in range(10, num_files+1): #1, num_files+1
            interfile = interfile_directory+'/'+interfile_prefix+'_'+str(x)+'.txt'
            outfile = outfile_directory+'/'+outfile_prefix+'_'+str(x)+'.txt'

            with open(outfile, 'w', newline='') as outfile:
                with open(interfile, 'r', encoding="iso-8859-1") as csvfile:
                    filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')

                    ###### MULTIPROCESSING INSERT ######
                    # This section will chunk the queries to allow more efficient splitting across processors,
                    # as the results are written to a file one at a time, but the queries can be
                    # sent/receieved without waiting for the results to be written.
                    # This approach avoids having to use the multiprocessing lock.
                    if __name__ == '__main__':
                        #lock = multiprocessing.Lock

                        print('INFO: Multiprocessing started')
                        cores = (multiprocessing.cpu_count()-1)
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

                    print('INFO: Multiprocessing completed')
                    ###### END MULTIPROCESSING INSERT ######

        return

    def trim_owlsim_output(self, input_dir, output_dir):
        """
        This function will take the results from the OWLSim queries, trim any unsuccessful queries,
        and return a file or files of successful queries.
        :return:
        """

        with open(output_file, 'w', newline='') as outfile:
            for i in range(1, 2):
                input_file = input_dir+str(i)+'.txt'
                output_file = output_dir+str(i)+'.txt'

                with open(input_file) as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        owlsim_data = json.loads(json_line)
                        (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag) = owlsim_data
                        if query_flag == 'success':
                            write_data = (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                            json.dump(write_data, outfile)
                            outfile.write('\n')
                            #print(query_flag)


        return

    def assemble_owlsim_gene_candidates(self, input_dir, output_file):
        """
        This function will take a trimmed set of OWLSim queries,
        assemble gene candidate scores for a disease/genotype,
        sort the results, and write the results to an output file.
        :return:
        """
        #output_file = output_dir+str(i)+'.txt'
        with open(output_file, 'w', newline='') as outfile:
            human_disease_gene_prediction_hash = {}
            for i in range(1, 26):
                print('Processing file number '+str(i))
                input_file = input_dir+str(i)+'.txt'


                with open(input_file) as handle:
                    for line in handle:
                        json_line = line.rstrip()
                        owlsim_data = json.loads(json_line)
                        (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag) = owlsim_data
                        if query_flag == 'success':

                            # Need a hash for the gene and

                            if entity_a not in human_disease_gene_prediction_hash:
                                entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                human_disease_gene_prediction_hash[entity_a] = {entity_b : entity_b_hash}
                                #human_disease_gene_prediction_hash[entity_a].update(entity_b_hash)
                            else:
                                entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                human_disease_gene_prediction_hash[entity_a].update({entity_b : entity_b_hash})



                #for i in human_disease_gene_prediction_hash:
                    #print('Disease: '+str(i)+' Match: '+str(human_disease_gene_prediction_hash[i]))
                            #write_data = (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                            #json.dump(write_data, outfile)
                            #outfile.write('\n')
                            #print(query_flag)


        return


    def assemble_owlsim_gene_candidate_alternate(self):
        """
        This function will take a trimmed set of OWLSim queries,
        assemble gene candidate scores for a disease/genotype,
        sort the results, and write the results to an output file.
        :return:
        """
        #output_file = 'out/owlsim/human_disease_gene_candidate_predictions/'+str()+'.txt'
        #with open(output_file, 'w', newline='') as outfile:
        human_disease_list_file = 'out/owlsim/human_disease_gene_candidate_predictions/human_disease_list.txt'
        human_disease_list = []
        current_disease = ''
        hvm_directory =  'out/owlsim/human_disease_mouse_gene/human_disease_mouse_gene_results_'
        hvz_directory = 'out/owlsim/human_disease_zebrafish_gene/human_disease_zebrafish_gene_results_'
        for i in range(1, 26):
            print('Processing file number '+str(i))
            input_file = hvm_directory+str(i)+'.txt'

            with open(input_file) as handle:
                for line in handle:
                    json_line = line.rstrip()
                    owlsim_data = json.loads(json_line)
                    (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag) = owlsim_data
                    entity_a = re.sub(':', '_', entity_a)

                    if query_flag == 'success':
                        if current_disease == '':
                            current_disease = entity_a
                            output_file = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(entity_a)+'.txt'
                            human_disease_list.append(entity_a)
                            human_disease_gene_prediction_hash = {}
                            entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                            human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                        else:
                            if current_disease == entity_a:
                                entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                            else:
                                print('Switching to disease '+str(entity_a))
                                with open(output_file, 'wb') as outfile:
                                    pickle.dump(human_disease_gene_prediction_hash, outfile)
                                current_disease = entity_a
                                output_file = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(entity_a)+'.txt'
                                if entity_a not in human_disease_list:
                                    human_disease_list.append(entity_a)
                                    human_disease_gene_prediction_hash = {}
                                    entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                    human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                                    #human_disease_gene_prediction_hash[entity_a] = {entity_b : entity_b_hash}
                                    #human_disease_gene_prediction_hash[entity_a].update(entity_b_hash)
                                else:
                                    with open(output_file, 'rb') as outfile:
                                        human_disease_gene_prediction_hash = pickle.load(outfile)
                                    entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                    human_disease_gene_prediction_hash[entity_b] = entity_b_hash


        for i in range(1, 10):
            print('Processing file number '+str(i))
            input_file = hvz_directory+str(i)+'.txt'

            with open(input_file) as handle:
                for line in handle:
                    json_line = line.rstrip()
                    owlsim_data = json.loads(json_line)
                    (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag) = owlsim_data
                    entity_a = re.sub(':', '_', entity_a)

                    if query_flag == 'success':
                        if current_disease == '':
                            current_disease = entity_a
                            output_file = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(entity_a)+'.txt'
                            human_disease_list.append(entity_a)
                            human_disease_gene_prediction_hash = {}
                            entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                            human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                        else:
                            if current_disease == entity_a:
                                entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                            else:
                                print('Switching to disease '+str(entity_a))
                                with open(output_file, 'wb') as outfile:
                                    pickle.dump(human_disease_gene_prediction_hash, outfile)
                                current_disease = entity_a
                                output_file = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(entity_a)+'.txt'
                                if entity_a not in human_disease_list:
                                    human_disease_list.append(entity_a)
                                    human_disease_gene_prediction_hash = {}
                                    entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                    human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                                    #human_disease_gene_prediction_hash[entity_a] = {entity_b : entity_b_hash}
                                    #human_disease_gene_prediction_hash[entity_a].update(entity_b_hash)
                                else:
                                    with open(output_file, 'rb') as outfile:
                                        human_disease_gene_prediction_hash = pickle.load(outfile)
                                    entity_b_hash = {'maxIC' : maxIC, 'simJ' : simJ, 'ICCS' : ICCS, 'simIC' : simIC}
                                    human_disease_gene_prediction_hash[entity_b] = entity_b_hash
                                    #with open(output_file, 'wb') as outfile:
                                        #pickle.dump(human_disease_gene_prediction_hash, outfile)


                #for i in human_disease_gene_prediction_hash:
                    #print('Disease: '+str(i)+' Match: '+str(human_disease_gene_prediction_hash[i]))
                            #write_data = (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)
                            #json.dump(write_data, outfile)
                            #outfile.write('\n')
                            #print(query_flag)
            with open(human_disease_list_file, 'wb') as outfile:
                pickle.dump(human_disease_list, outfile)

        return

    def assemble_owlsim_top_20_gene_candidates(self):

        owlsim_disease_list = []

        disease_subset = read_only_disease_subset # ['OMIM_157900', 'ORPHANET_904', 'ORPHANET_84', 'ORPHANET_46348', 'OMIM_272120', 'ORPHANET_2812', 'ORPHANET_791', 'ORPHANET_478', 'ORPHANET_110', 'OMIM_209900', 'OMIM_614592', 'ORPHANET_1873', 'OMIM_305400']

        human_disease_list_file = 'out/owlsim/human_disease_gene_candidate_predictions/human_disease_list.txt'
        with open(human_disease_list_file, 'rb') as handle:
            human_disease_list = pickle.load(handle)


        with open('inter/hpo/human_gene_id_to_label_hash.txt', 'rb') as handle:
            human_gene_id_to_label_hash = pickle.load(handle)
        with open('inter/mgi/mouse_gene_id_to_label_hash.txt', 'rb') as handle:
            mouse_gene_id_to_label_hash = pickle.load(handle)
        with open('inter/zfin/zebrafish_gene_id_to_label_hash.txt', 'rb') as handle:
            zebrafish_gene_id_to_label_hash = pickle.load(handle)
        gene_id_to_label_hash = {}
        gene_id_to_label_hash.update(human_gene_id_to_label_hash)
        gene_id_to_label_hash.update(mouse_gene_id_to_label_hash)
        gene_id_to_label_hash.update(zebrafish_gene_id_to_label_hash)

        #print(human_disease_gene_prediction_hash)
        #disease_id = 'ORPHANET:2824'
        for disease_id in disease_subset:
            try:
                disease_id = re.sub(':', '_', disease_id)
                print('Processing disease ID '+str(disease_id)+'.')
                with open('out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(disease_id)+'.txt', 'rb') as handle:
                    human_disease_gene_prediction_hash = pickle.load(handle)
                #disease_id = re.sub('_', ':', x)

                top_20_max_ic = heapq.nlargest(20, human_disease_gene_prediction_hash, key=lambda k:human_disease_gene_prediction_hash[k]['maxIC'])
                top_20_iccs = heapq.nlargest(20, human_disease_gene_prediction_hash, key=lambda k:human_disease_gene_prediction_hash[k]['ICCS'])
                top_20_sim_ic = heapq.nlargest(20, human_disease_gene_prediction_hash, key=lambda k:human_disease_gene_prediction_hash[k]['simIC'])
                top_20_sim_j = heapq.nlargest(20, human_disease_gene_prediction_hash, key=lambda k:human_disease_gene_prediction_hash[k]['simJ'])
                #print(top_20_max_ic)
                with open('out/owlsim/human_disease_gene_candidate_predictions/top_twenty_genes/maxic/'+str(disease_id)+'.txt', 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in top_20_max_ic:
                        gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        max_ic = '('+str(round(human_disease_gene_prediction_hash[gene_candidate_id]['maxIC'], 2))+')'
                        #output_row = (gene_candidate_id, gene_candidate_label, max_ic)
                        if re.match('MGI:.*', gene_candidate_id):
                            species = '[Mouse]'
                        elif re.match('ZFIN:.*', gene_candidate_id):
                            species = '[Zebrafish]'
                        else:
                            print('Species not found.')
                            species = 'Unknown'
                        output_row = (gene_candidate_label, max_ic, species)
                        csvwriter.writerow(output_row)
                        #print('Gene candidate: '+str(gene_candidate)+', MaxIC:'+str(human_disease_gene_prediction_hash[gene_candidate]['maxIC']))
                with open('out/owlsim/human_disease_gene_candidate_predictions/top_twenty_genes/iccs/'+str(disease_id)+'.txt', 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in top_20_iccs:
                        if re.match('MGI:.*', gene_candidate_id):
                            species = '[Mouse]'
                        elif re.match('ZFIN:.*', gene_candidate_id):
                            species = '[Zebrafish]'
                        else:
                            print('Species not found.')
                            species = 'Unknown'
                        gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        iccs = '('+str(round(human_disease_gene_prediction_hash[gene_candidate_id]['ICCS'], 2))+')'
                        #output_row = (gene_candidate_id, gene_candidate_label, iccs)
                        output_row = (gene_candidate_label, iccs, species)
                        csvwriter.writerow(output_row)
                with open('out/owlsim/human_disease_gene_candidate_predictions/top_twenty_genes/simic/'+str(disease_id)+'.txt', 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in top_20_sim_ic:
                        if re.match('MGI:.*', gene_candidate_id):
                            species = '[Mouse]'
                        elif re.match('ZFIN:.*', gene_candidate_id):
                            species = '[Zebrafish]'
                        else:
                            print('Species not found.')
                            species = 'Unknown'
                        gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        sim_ic = '('+str(round(human_disease_gene_prediction_hash[gene_candidate_id]['simIC'], 2))+')'
                        #output_row = (gene_candidate_id, gene_candidate_label, sim_ic)
                        output_row = (gene_candidate_label, sim_ic, species)
                        csvwriter.writerow(output_row)
                with open('out/owlsim/human_disease_gene_candidate_predictions/top_twenty_genes/simj/'+str(disease_id)+'.txt', 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in top_20_sim_j:
                        if re.match('MGI:.*', gene_candidate_id):
                            species = '[Mouse]'
                        elif re.match('ZFIN:.*', gene_candidate_id):
                            species = '[Zebrafish]'
                        else:
                            print('Species not found.')
                            species = 'Unknown'
                        gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        sim_j = '('+str(round(human_disease_gene_prediction_hash[gene_candidate_id]['simJ'], 2))+')'
                        #output_row = (gene_candidate_id, gene_candidate_label, sim_j)
                        output_row = (gene_candidate_label, sim_j, species)
                        csvwriter.writerow(output_row)
                owlsim_disease_list.append(disease_id)
            except:
                print('OWLSim candidate processing failed for disease ID '+str(disease_id)+'.')


        with open('inter/omim/owlsim_disease_list.txt', 'wb') as handle:
            pickle.dump(owlsim_disease_list, handle)


        return


    ####### PHENOLOG DATA PROCESSING ####### DOCUMENTATION COMPLETED

    def generate_random_data(self, pheno_gene_hash, common_orthologs, out_dir, limit=1000):
        """
        This function creates random data sets for the determination of the FDR cutoff for the significant phenologs.
        It takes as input the real phenotype-ortholog hash and the common orthologs between two species.
        It creates a test phenotype-ortholog hash using the real phenotype-ortholog hash as a guide.
        This allows for the creation of a random data set of similar size and shape
        (same number of phenotypes with the same number of associated orthologs for each phenotype).

        :param pheno_gene_hash: Phenotype-ortholog hash file for a given species.
        :param common_orthologs: Table of common orthologs between two species for the random data set.
        :param out_dir: Directory for saving the random data set files.
        :param limit: Total number of random data sets to create.
        :return:
        """
        print('INFO: Creating random data sets.')
        counter = 1
        while counter <= limit:
            print('INFO: Creating random data set '+str(counter)+' out of '+str(limit)+'.')
            test_pheno_ortholog_hash = {}

            with open(pheno_gene_hash, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            with open(common_orthologs, 'rb') as handle:
                orthologs = pickle.load(handle)
            random.shuffle(orthologs)
            for i in pheno_ortholog_hash:
                test_pheno_ortholog_hash[i] = []
                ortholog_draw = orthologs
                random.shuffle(ortholog_draw)

                for j in pheno_ortholog_hash[i]:
                    random.shuffle(ortholog_draw)
                    test_pheno_ortholog_hash[i].append(ortholog_draw[0])

            with open((out_dir+'random_'+str(counter)+'.txt'), 'wb') as handle:
                pickle.dump(test_pheno_ortholog_hash, handle)
            counter += 1
        return

    def set_stage_for_fdr_calculation(self):
        """
        This function manages the calculation of the FDR for the phenolog calculations.
        It sets up the multiprocessing of the 1000 random data sets,
        calling out to the multiprocess_fdr_calculation function,
        assembles the p-value for each random data set where 5% of the identified phenologs are significant,
        and then returns the average p-value from this set, which will be used as the FDR cutoff for the
        phenolog calculations on the real data set.
        :return:
        """
        print('INFO: Setting stage for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        fdr_global_p_value_list = []

        ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)
            print('INFO: Multiprocessing started')
            fdr_global_p_value_list = [pool.map(multiprocess_fdr_calculation, range(1,1001))]

            print('INFO: Multiprocessing completed')
            ###### END MULTIPROCESSING INSERT ######

            with open('inter/random/fdr/fdr_global_p_value_list.txt', 'wb') as handle:
                pickle.dump(fdr_global_p_value_list, handle)
            p_value_sum = float(0)
            for x in fdr_global_p_value_list:
                p_value_sum += x
            fdr_cutoff = float(p_value_sum)/float(len(fdr_global_p_value_list))
            print('The empirical FDR adjustment cutoff is '+str(fdr_cutoff)+'.')

        return fdr_cutoff

    def assemble_partial_fdr(self):
        """
        This function is used to assemble the p-values from saved files,
        should the set_stage_for_fdr_calculation fail during processing or
        the FDR calculations be divided across machines.
        This function is used to assemble the p-values for each random data set
        where 5% of the identified phenologs are significant,
        and then returns the average p-value from this set, which will be used as the FDR cutoff for the
        phenolog calculations on the real data set.
        :return:
        """
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
        p_value_sum = float(0)
        for x in range(0,1000):
            print(fdr_global_p_value_list[x])
            p_value_sum += fdr_global_p_value_list[x]
        fdr_cutoff = float(p_value_sum)/float(len(fdr_global_p_value_list))
        print('The empirical FDR adjustment cutoff is '+str(fdr_cutoff)+'.')

        with open('inter/random/fdr/fdr_global_p_value_list.txt', 'wb') as handle:
            pickle.dump(fdr_global_p_value_list, handle)

        return

    def perform_phenolog_calculations_for_fdr(self, species_a_po_hash, species_b_po_hash, shared_orthologs):
        """
        This function performs the phenolog calculations between
        the phenotypes of two species from the random data sets.
        For each cross-species pair of phenotypes, the associated orthologs are compared for matches.
        If two phenotypes have one or more matching orthologs, the hypergeometric probability is calculated.
        P-values are added to a list, which is returned after all phenotype comparisons are complete.
        :param species_a_po_hash: The phenotype-ortholog hash for species A.
        :param species_b_po_hash: The phenotype-ortholog hash for species B.
        :param shared_orthologs: The file containing the orthologs shared between the two compared species.
        :return: List of p-values from the hypergeometric probability calculation.
        """

        total_ortholog_matches = 0
        total_ortholog_nonmatches = 0

        total_hyp_calcs = 0
        phenolog_p_value_list = []

        species_a_pheno_gene_hash = species_a_po_hash
        species_b_pheno_gene_hash = species_b_po_hash
        # Iterate through the phenotypes for each species,
        # determining the number of ortholog matches between the orthologs associated with each phenotype.
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
                ortholog_matches = 0
                ortholog_non_matches = 0
                phenotype_b_ortholog_count = len(species_b_orthologs)
                for k in species_a_orthologs:
                    # Orthologs for species A
                    species_a_ortholog = k
                    for l in species_b_orthologs:
                        # Orthologs for species B
                        species_b_ortholog = l
                        if species_a_ortholog == species_b_ortholog:
                            ortholog_matches += 1
                            total_ortholog_matches += 1
                        else:
                            ortholog_non_matches += 1
                            total_ortholog_nonmatches += 1

                if ortholog_matches > 0:
                    # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
                    # N = total number of orthologs shared between species
                    # n = nummber of orthologs in species A phenotype
                    # m = nummber of orthologs in species B phenotype
                    # c = number of common orthologs between phenotypes (ortholog matches)
                    m = float(phenotype_b_ortholog_count)
                    n = float(phenotype_a_ortholog_count)
                    N = float(shared_orthologs)
                    c = float(ortholog_matches)
                    prb = float(hypergeom.pmf(c, N, m, n))
                    phenolog_p_value_list.append(prb)
                    total_hyp_calcs += 1
        print('Total Matches: '+str(total_ortholog_matches))
        print('Total non-matches: '+str(total_ortholog_nonmatches))
        print('Total phenolog calculations: '+str(total_hyp_calcs))

        return phenolog_p_value_list

    def perform_phenolog_calculations(self, inter1, inter2, out, shared_orthologs, fdr_cutoff):
        """
        This function performs the calculations to determine whether two phenotypes are phenologs
        using the previously determined FDR cutoff.
        :param inter1:
        :param inter2:
        :param out:
        :param shared_orthologs:
        :param fdr_cutoff:
        :return:
        """
        print('INFO: Performing phenolog calculations.')
        total_ortholog_matches = 0
        total_ortholog_nonmatches = 0

        with open(inter1, 'rb') as handle:
            species_a_pheno_gene_hash = pickle.load(handle)
        with open(inter2, 'rb') as handle:
            species_b_pheno_gene_hash = pickle.load(handle)
        with open(shared_orthologs, 'rb') as handle:
            num_shared_orthologs = len(pickle.load(handle))
        with open(out, 'w', newline='') as outfile:
            # Iterate through the phenotypes for each species,
            # determining the number of ortholog matches between the orthologs associated with each phenotype.
            for i in species_a_pheno_gene_hash:
                species_a_phenotype_id = i
                species_a_orthologs = species_a_pheno_gene_hash[i]
                phenotype_a_ortholog_count = len(species_a_orthologs)
                for j in species_b_pheno_gene_hash:
                    species_b_phenotype_id = j
                    species_b_orthologs = species_b_pheno_gene_hash[j]
                    ortholog_matches = 0
                    ortholog_non_matches = 0
                    phenotype_b_ortholog_count = len(species_b_orthologs)
                    for k in species_a_orthologs:
                        species_a_ortholog = k
                        for l in species_b_orthologs:
                            species_b_ortholog = l
                            if species_a_ortholog == species_b_ortholog:
                                ortholog_matches += 1
                                total_ortholog_matches += 1
                            else:
                                ortholog_non_matches += 1
                                total_ortholog_nonmatches += 1

                    if ortholog_matches > 0:
                        # Perform the hypergeometric probability calculation.
                        # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
                        # N = total number of orthologs shared between species
                        # n = number of orthologs in species A phenotype
                        # m = number of orthologs in species B phenotype
                        # c = number of common orthologs between phenotypes (ortholog matches)

                        m = float(phenotype_b_ortholog_count)
                        n = float(phenotype_a_ortholog_count)
                        N = float(num_shared_orthologs)
                        c = float(ortholog_matches)
                        prb = float(hypergeom.pmf(c, N, m, n))
                        if prb <= fdr_cutoff:
                            significance = 'Significant'
                        else:
                            significance = 'Not Significant'

                        # For each calculation, write the output to a file in JSON format.
                        sequence = (species_a_phenotype_id, species_a_orthologs, phenotype_a_ortholog_count, species_b_phenotype_id, species_b_orthologs, phenotype_b_ortholog_count, num_shared_orthologs, ortholog_matches, prb, fdr_cutoff, significance)
                        json.dump(sequence, outfile)
                        outfile.write('\n')

            print('Total Matches: '+str(total_ortholog_matches))
            print('Total non-matches: '+str(total_ortholog_nonmatches))

        return

    def assemble_significant_phenologs(self):
        """
        This function assembles the significant phenologs from the files of phenologs for species pairs.
        :return:
        """
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
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)

        return

    def assemble_significant_phenologs_with_scores(self):
        """
        This function assembles the significant phenologs from the files of phenologs for species pairs,
        with the associated hypergeometric probabilty p-values.
        :return:
        """
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
                            csvwriter.writerow(output_row)
                            csvwriter2.writerow(output_row)
            print('INFO: Done assembling mouse-zebrafish phenologs.')

        return

    def assemble_hvm_phenologs(self):
        """
        This function iterates through the significant human-mouse phenologs,
        creates a combo field using the two phenotype IDs that will be used
        as a lookup for the phenolog extension calculations.
        :return:
        """
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
        """
        This function iterates through the significant human-zebrafish phenologs,
        creates a combo field using the two phenotype IDs that will be used
        as a lookup for the phenolog extension calculations.
        :return:
        """
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
        """
        This function iterates through the significant mouse-zebrafish phenologs,
        creates a combo field using the two phenotype IDs that will be used
        as a lookup for the phenolog extension calculations.
        :return:
        """
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

    def annotate_significant_phenologs_with_scores_and_labels(self):
        with open('inter/ontologies/hp_hash.txt', 'rb') as handle:
            hp_hash = pickle.load(handle)
        with open('raw/hpo/dvp.pr_nlx_151835_2', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:
                # Read in a row and split into individual variables
                (e_uid, phenotype_id, phenotype_label, gene_id, gene_num,
                 gene_label, v_uid, v_uuid, v_lastmodified) = row
                if phenotype_id not in hp_hash:
                    hp_hash[phenotype_id] = phenotype_label
        with open('inter/ontologies/mp_hash.txt', 'rb') as handle:
            mp_hash = pickle.load(handle)
        with open('inter/ontologies/zp_hash.txt', 'rb') as handle:
            zp_hash = pickle.load(handle)

        with open('out/phenolog/all_significant_phenologs_labelled_scored.txt', 'w', newline='') as csvoutfile:
            csvwriter = csv.writer(csvoutfile, delimiter='\t', quotechar='\"')
            with open('out/phenolog/all_significant_phenologs_scored.txt', 'r') as csvreadfile:
                filereader = csv.reader(csvreadfile, delimiter='\t', quotechar='\"')
                for row in filereader:

                    (phenotype_a, phenotype_a_ortholog_count, phenotype_b, phenotype_b_ortholog_count,
                     ortholog_matches, probability) = row
                    if re.match('HP:.*', phenotype_a):
                        phenotype_a_label = hp_hash[phenotype_a]
                    elif re.match('MP:.*', phenotype_a):
                        phenotype_a_label = mp_hash[phenotype_a]
                    elif re.match('ZP:.*', phenotype_a):
                        phenotype_a_label = zp_hash[phenotype_a]
                    else:
                        phenotype_a_label = ''
                        print('INFO: No label found for '+phenotype_a+'.')

                    if re.match('HP:.*', phenotype_b):
                        phenotype_b_label = hp_hash[phenotype_b]
                    elif re.match('MP:.*', phenotype_b):
                        phenotype_b_label = mp_hash[phenotype_b]
                    elif re.match('ZP:.*', phenotype_b):
                        phenotype_b_label = zp_hash[phenotype_b]
                    else:
                        phenotype_b_label = ''
                        print('INFO: No label found for '+phenotype_b+'.')

                    output_row = (phenotype_a, phenotype_a_label, phenotype_a_ortholog_count, phenotype_b, phenotype_b_label, phenotype_b_ortholog_count, ortholog_matches, probability)
                    csvwriter.writerow(output_row)

        return


    ####### PHENOLOG EXTENSION DATA PROCESSING #######

    def parse_zp(self, raw, inter):
        """
        This function parses the zebrafish phenotype ontology to
        obtain a hash of zebrafish phenotypes and their associated labels.
        This phenotype hash will be used to generate random zebrafish genotype-phenotype data sets.
        :param raw: The zebrafish phenotype ontology.
        :param inter: The output file for the zebrafish phenotype hash.
        :return:
        """
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
        """
        This function parses the mouse phenotype ontology to
        obtain a hash of mouse phenotypes and their associated labels.
        This phenotype hash will be used to generate random mouse genotype-phenotype data sets.
        :param raw: The mouse phenotype ontology.
        :param inter: The output file for the mouse phenotype hash.
        :return:
        """
        mp_hash = {}
        with open(raw, 'r') as handle:
            for line in handle:
                if re.match('id:.*', line.rstrip()):
                    mp_id = line.rstrip()
                    mp_id = re.sub('id: ', '', mp_id)
                    next_line = next(handle)
                    if re.match('name:.*', next_line.rstrip()):
                        mp_label = next_line.rstrip()
                        mp_label = re.sub('name: ', '', mp_label)
                        if mp_id not in mp_hash:
                            mp_hash[mp_id] = mp_label

        with open(inter, 'wb') as handle:
            pickle.dump(mp_hash, handle)

        return

    def parse_hp(self, raw, inter):
        """
        This function parses the human phenotype ontology to
        obtain a hash of human phenotypes and their associated labels.
        This phenotype hash will be used to generate random human genotype-phenotype data sets.
        :param raw: The human phenotype ontology.
        :param inter: The output file for the human phenotype hash.
        :return:
        """

        hp_hash = {}
        with open(raw, 'r') as handle:
            for line in handle:
                if re.match('id:.*', line.rstrip()):
                    hp_id = line.rstrip()
                    hp_id = re.sub('id: ', '', hp_id)
                    next_line = next(handle)
                    if re.match('name:.*', next_line.rstrip()):
                        hp_label = next_line.rstrip()
                        hp_label = re.sub('name: ', '', hp_label)
                        if hp_id not in hp_hash:
                            hp_hash[hp_id] = hp_label

        with open(inter, 'wb') as handle:
            pickle.dump(hp_hash, handle)

        return

    def set_stage_for_extension_fdr_calculation(self, limit=1000):
        """

        :param limit:
        :return:
        """
        print('INFO: Setting stage for second FDR estimation.')
        # Need to calculate phenolog extension for each pairwise species and combine in order to get a full
        # set of 'genologs' (?) for proper estimation of FDR.

        ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':


            cores = (multiprocessing.cpu_count()-1)
            #cores = 2
            pool = multiprocessing.Pool(processes=cores)

            print('INFO: Multiprocessing started')

            fdr_global_p_value_list = [pool.map(multiprocess_ext_fdr_calculation, range(1, 1001))]
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
        """
        This function handles the multiprocessing setup of the creation of
        random human data sets for the phenolog extension FDR calculations.
        It calls out to the multiprocess_generate_random_human_ext_data function.
        :return:
        """
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
        """
        This function handles the multiprocessing setup of the creation of
        random mouse data sets for the phenolog extension FDR calculations.
        It calls out to the multiprocess_generate_random_mouse_ext_data function.
        :return:
        """
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
        """
        This function handles the multiprocessing setup of the creation of
        random zebrafish data sets for the phenolog extension FDR calculations.
        It calls out to the multiprocess_generate_random_zebrafish_ext_data function.
        :return:
        """
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
        """

        :param species_a_gp_hash:
        :param species_b_gp_hash:
        :param shared_phenologs:
        :return:
        """
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

        comparison_counter = 0
        total_phenotype_matches = 0
        total_phenotype_nonmatches = 0
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
            genotype_a_phenotype_count = len(species_a_phenotypes)

            for j in species_b_geno_pheno_hash:
                # Genotype for species B
                species_b_genotype_id = j
                species_b_phenotypes = species_b_geno_pheno_hash[j]
                phenotype_matches = 0
                phenotype_non_matches = 0
                genotype_b_phenotype_count = len(species_b_phenotypes)
                comparison_counter += 1
                if re.match('.*000000', str(comparison_counter)):
                    print('INFO: Processing comparison '+str(comparison_counter)+' out of '+str(total_comparisons)+' total comparisons to perform.')
                for k in species_a_phenotypes:
                    # Orthologs for species A
                    species_a_phenotype = k
                    for l in species_b_phenotypes:
                        # Orthologs for species B
                        species_b_phenotype = l
                        ab_combo = species_a_phenotype+'_'+species_b_phenotype
                        ba_combo = species_b_phenotype+'_'+species_a_phenotype
                        if ab_combo in shared_phenologs or ba_combo in shared_phenologs:
                            phenotype_matches += 1
                            total_phenotype_matches += 1
                        else:
                            phenotype_non_matches += 1
                            total_phenotype_nonmatches += 1

                if phenotype_matches > 0:
                    # Relevent SciPy documentation: http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html#scipy.stats.hypergeom
                    # N = total number of phenologs shared between species
                    # n = number of phenotypes in species A disease/genotype
                    # m = number of phenotypes in species B disease/genotype
                    # c = number of common phenotypes between disease/genotype (phenotype matches)
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



        return phenolog_ext_p_value_list

    def perform_phenolog_calculations_for_ext_fdr_hvz(self, species_a_gp_hash, species_b_gp_hash, out_file):
        """
        This function
        and sends the
        :param species_a_gp_hash:
        :param species_b_gp_hash:
        :param out_file:
        :return:
        """

        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

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

        ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':
            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)
            print('INFO: Starting multiprocessing.')
            for results in pool.imap_unordered(multiprocess_ext_fdr_calculation_hvz, comparison_list, chunksize=100):
                if results is not None:
                    phenolog_ext_p_value_list.append(results)
                del results
        ###### END MULTIPROCESSING INSERT ######

        print('INFO: Done with multiprocessing.')
        with open(out_file, 'wb') as handle:
            pickle.dump(phenolog_ext_p_value_list, handle)
        '''
        print('INFO: Sorting p-values for random data set '+str(i)+'.')

        fdr_p_value_list.sort()
        with open('inter/random/fdr/fdr_p_value_list_random_set_'+str(i)+'.txt', 'wb') as handle:
            pickle.dump(fdr_p_value_list, handle)
        #print(fdr_p_value_list)
        cutoff_position = math.ceil((len(fdr_p_value_list))*0.05) - 1
        #print(fdr_p_value_list[cutoff_position])
        fdr_cutoff_value = fdr_p_value_list[cutoff_position]
        '''

        del phenolog_ext_p_value_list
        del comparison_list
        gc.collect()
        return #phenolog_ext_p_value_list

    def perform_phenolog_calculations_for_ext_fdr_hvz_alternate(self, species_a_gp_hash, species_b_gp_hash, out_file):
        """
        This function
        and sends the
        :param species_a_gp_hash:
        :param species_b_gp_hash:
        :param out_file:
        :return:
        """

        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

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

        ###### MULTIPROCESSING INSERT ######
        if __name__ == '__main__':
            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)
            print('INFO: Starting multiprocessing.')
            with open(out_file, 'w', newline='') as csvfile:
                csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')


                for results in pool.imap_unordered(multiprocess_ext_fdr_calculation_hvz, comparison_list, chunksize=100):
                    if results is not None:
                        output_row = (results)
                        csvwriter.writerow(output_row)
                    del results
        ###### END MULTIPROCESSING INSERT ######

        print('INFO: Done with multiprocessing.')
        #with open(out_file, 'wb') as handle:
            #pickle.dump(phenolog_ext_p_value_list, handle)
        '''
        print('INFO: Sorting p-values for random data set '+str(i)+'.')

        fdr_p_value_list.sort()
        with open('inter/random/fdr/fdr_p_value_list_random_set_'+str(i)+'.txt', 'wb') as handle:
            pickle.dump(fdr_p_value_list, handle)
        #print(fdr_p_value_list)
        cutoff_position = math.ceil((len(fdr_p_value_list))*0.05) - 1
        #print(fdr_p_value_list[cutoff_position])
        fdr_cutoff_value = fdr_p_value_list[cutoff_position]
        '''

        del phenolog_ext_p_value_list
        del comparison_list
        gc.collect()
        return #phenolog_ext_p_value_list

    def perform_phenolog_calculations_for_ext_fdr_hvm(self, species_a_gp_hash, species_b_gp_hash, out_file):
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

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

            for results in pool.imap_unordered(multiprocess_ext_fdr_calculation_hvm, comparison_list, chunksize=100):
                if results is not None:
                    phenolog_ext_p_value_list.append(results)
                    #print(results)


            #print('INFO: Done processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
        print('INFO: Done with multiprocessing.')
        with open(out_file, 'wb') as handle:
            pickle.dump(phenolog_ext_p_value_list, handle)
        #print('Total Matches: '+str(total_phenotype_matches))
        #print('Total non-matches: '+str(total_phenotype_nonmatches))
        #print('Total phenolog calculations: '+str(total_hyp_calcs))

        phenolog_ext_p_value_list = None
        comparison_list = None
        gc.collect()
        return #phenolog_ext_p_value_list

    def perform_phenolog_calculations_for_ext_fdr_mvz(self, species_a_gp_hash, species_b_gp_hash, out_file):
        #print('INFO: Performing phenolog calculations for FDR estimation.')
        # Need to calculate phenologs for each pairwise species and combine in order to get a full
        # set of phenologs for proper estimation of FDR.

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

            for results in pool.imap_unordered(multiprocess_ext_fdr_calculation_mvz, comparison_list, chunksize=100):
                if results is not None:
                    phenolog_ext_p_value_list.append(results)
                    #print(results)


            #print('INFO: Done processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
        print('INFO: Done with multiprocessing.')
        with open(out_file, 'wb') as handle:
            pickle.dump(phenolog_ext_p_value_list, handle)
        #print('Total Matches: '+str(total_phenotype_matches))
        #print('Total non-matches: '+str(total_phenotype_nonmatches))
        #print('Total phenolog calculations: '+str(total_hyp_calcs))

        phenolog_ext_p_value_list = None
        comparison_list = None
        gc.collect()
        return #phenolog_ext_p_value_list

    def identify_significance_threshold_for_random_data_sets(self):
        """
        This function opens saved p-value lists for each species pair from
        a processed random phenolog extension data set, sorts the list,
        identifies the p-value where 5% of the matches would be significant,
        and saves this value to a new file.
        :return:
        """

        for i in range(2, 3):
            print('INFO: Starting cutoff determination for data set '+str(i)+'.')
            phenolog_ext_p_value_list = []
            hvm_file = 'inter/phenolog_ext/hvm_p_values/hvm_p_values_'+str(i)+'.txt'
            hvz_file = 'inter/phenolog_ext/hvz_p_values/hvz_p_values_'+str(i)+'.txt'
            mvz_file = 'inter/phenolog_ext/mvz_p_values/mvz_p_values_'+str(i)+'.txt'
            out_file = 'inter/phenolog_ext/p_value_cutoffs/threshold_p_value_'+str(i)+'.txt'
            with open(hvm_file, 'rb') as handle:
               hvm_p_value_list = pickle.load(handle)
            with open(hvz_file, 'rb') as handle:
               hvz_p_value_list = pickle.load(handle)
            with open(mvz_file, 'rb') as handle:
               mvz_p_value_list = pickle.load(handle)
            for x in hvm_p_value_list:
                phenolog_ext_p_value_list.append(float(x))
                #print(x)
            for y in hvz_p_value_list:
                phenolog_ext_p_value_list.append(float(y))
                #print(y)
            for z in mvz_p_value_list:
                phenolog_ext_p_value_list.append(float(z))
                #print(z)
            print('INFO: Sorting p values.')
            phenolog_ext_p_value_list.sort()
            #print(fdr_p_value_list)
            cutoff_position = math.ceil((len(phenolog_ext_p_value_list))*0.05) - 1
            print(phenolog_ext_p_value_list[cutoff_position])
            print('Initial cutoff posiiton: '+str(cutoff_position)+'.')
            fdr_cutoff_value = float(phenolog_ext_p_value_list[cutoff_position])
            while fdr_cutoff_value == 0.0:
                cutoff_position += 1
                fdr_cutoff_value = float(phenolog_ext_p_value_list[cutoff_position])
            print('Total number of p-values: '+str(len(phenolog_ext_p_value_list))+'.')

            print(phenolog_ext_p_value_list[cutoff_position])
            with open(out_file, 'wb') as handle:
                pickle.dump(fdr_cutoff_value, handle)
            print('INFO: Finished cutoff determination for data set '+str(i)+'.')
            with open(out_file, 'rb') as handle:
                final_cutoff = pickle.load(handle)
            print('Adjusted cutoff positon: '+str(cutoff_position)+'.')
            print(final_cutoff)
        return

    def count_zeroes(self):

        hvz_file = 'inter/phenolog_ext/hvz_p_values/hvz_p_values_2.txt'
        with open(hvz_file, 'rb') as handle:
           hvz_p_value_list = pickle.load(handle)
        count = 0
        for y in hvz_p_value_list:
            if y == 0.0:
                count += 1
        print('Number of zeroes: '+str(count)+'.')
        print('5 percent cutoff: '+str(len(hvz_p_value_list) * .05)+'.')
        print('Total p-values:'+str(len(hvz_p_value_list))+'.')
        return

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


    ####### PHENOLOG ORTHOLOG PHENOTYPE MATRICES ####### DOCUMENTATION COMPLETED

    def assemble_ortholog_phenotype_matrix(self):
        """
        This function assembles the ortholog-phenotype matrix to be used in gene candidate predictions.
        It utilizes the already created phenotype-ortholog hashes for each species.
        It then saves the phenotype list, the ortholog list, and the ortholog-phenotype matrix to disk.
        """

        print('INFO: Assembling phenotype-ortholog matrices.')
        # Create the phenotype and ortholog lists from the existing phenotype-ortholog hashes for each species.
        phenotype_list = []
        ortholog_list = []
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            for i in pheno_ortholog_hash:
                if i not in phenotype_list:
                    phenotype_list.append(i)
                for j in pheno_ortholog_hash[i]:
                    if j not in ortholog_list:
                        ortholog_list.append(j)

        phenotype_list.sort()
        ortholog_list.sort()
        total_phenotypes = len(phenotype_list)
        print('INFO: Total number of phenotypes: '+str(total_phenotypes))
        total_orthologs = len(ortholog_list)
        print('INFO: Total number of orthologs: '+str(total_orthologs))

        # Create the ortholog-phenotype matrix as a matrix of zeroes.
        ortholog_phenotype_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))
        # Now, for each species, enter a 1 at each matrix coordinate corresponding to a phenotype and an associated ortholog.
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)
            for i in pheno_ortholog_hash:
                phenotype_index = phenotype_list.index(i)
                for j in pheno_ortholog_hash[i]:
                    ortholog_index = ortholog_list.index(j)
                    ortholog_phenotype_matrix[phenotype_index][ortholog_index] = 1
        print('INFO: Done assembling phenotype-ortholog matrices.')

        # Save all files to disk.
        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
            pickle.dump(phenotype_list, handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
            pickle.dump(ortholog_list, handle)
        numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
        numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)

        return (ortholog_phenotype_matrix, phenotype_list, ortholog_list)

    def assemble_ortholog_phenotype_matrix_alternate(self, limit=None):
        """
        This function assembles the ortholog-phenotype matrix to be used in gene candidate predictions.
        It utilizes the already created phenotype-ortholog hashes for each species.
        It then saves the phenotype list, the ortholog list, and the ortholog-phenotype matrix to disk.
        This version includes an optional limit for testing purposes.
        """

        print('INFO: Assembling phenotype-ortholog matrices.')
        # Create the phenotype and ortholog lists from the existing phenotype-ortholog hashes for each species.
        phenotype_list = []
        ortholog_list = []
        phenotype_limit = limit
        for x in ['inter/hpo/human_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_hash.txt']:
            with open(x, 'rb') as handle:
                pheno_ortholog_hash = pickle.load(handle)

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

        # Create the ortholog-phenotype matrix as a matrix of zeroes.
        ortholog_phenotype_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))

        # Now, for each species, enter a 1 at each matrix coordinate corresponding to a phenotype and an associated ortholog.
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

        print('INFO: Done assembling phenotype-ortholog matrices.')

        # Save all files to disk.
        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
            pickle.dump(phenotype_list, handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
            pickle.dump(ortholog_list, handle)
        numpy.save('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy', ortholog_phenotype_matrix)
        numpy.savetxt('inter/phenolog_gene_cand/ortholog_phenotype_matrix.txt', ortholog_phenotype_matrix)

        return (ortholog_phenotype_matrix, phenotype_list, ortholog_list)


    ####### PHENOLOG GENE CANDIDATE PREDICTION ALGORITHM #######

    def create_phenolog_gene_candidate_matrices(self):

        # Set the number of nearest neighbor phenotypes to consider for predictions.
        k = 11

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

        distance_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))
        distance_matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        weight_matrix = numpy.zeros((len(phenotype_list), len(phenotype_list)))

        print('INFO: '+str(distance_matrix_comparisons)+' matrix comparisons to process.')

        if __name__ == '__main__':

            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)

            #Takes ~65 seconds to reach this point.
            print('INFO: Assembling phenotype matrix coordinates.')

            for i in phenotype_list:
                input_phenotype_index_i = phenotype_list.index(i)

                print('INFO: Processing phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                matrix_coordinates = []
                for j in phenotype_list:
                    input_phenotype_index_j = phenotype_list.index(j)
                    matrix_coordinates.append([input_phenotype_index_i,input_phenotype_index_j])
                print('INFO: Done assembling phenotype matrix coordinates.')
                print('INFO: Starting multiprocessing.')
                results = pool.map(multiprocess_matrix_comparisons, matrix_coordinates)
                print('INFO: Processing results for phenotype '+str(input_phenotype_index_i+1)+' out of '+str(len(phenotype_list))+'.')
                for x in results:
                    (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient) = x
                    distance_matrix[phenotype_index_i][phenotype_index_j] = coefficient
                    weight_matrix[phenotype_index_i][phenotype_index_j] = hyp_prob
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
            with open('inter/phenolog_gene_cand/phenotype_list.txt', 'wb') as handle:
                pickle.dump(phenotype_list, handle)
            with open('inter/phenolog_gene_cand/ortholog_list.txt', 'wb') as handle:
                pickle.dump(ortholog_list, handle)

            print('DONE!')
        return

    def create_empty_phenolog_gene_candidate_matrices(self):


        # Set the number of nearest neighbor phenotypes to consider for predictions.
        k = 11

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
                    (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient) = x
                    distance_matrix[phenotype_index_i][phenotype_index_j] = coefficient
                    weight_matrix[phenotype_index_i][phenotype_index_j] = hyp_prob
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

        matrix_comparisons = (len(phenotype_list)*len(phenotype_list))
        print('INFO: '+str(matrix_comparisons)+' matrix comparisons to process.')

        if __name__ == '__main__':
            cores = (multiprocessing.cpu_count()-1)
            pool = multiprocessing.Pool(processes=cores)

            print('INFO: Loading distance matrix.')
            distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
            print('INFO: Loading weight matrix.')
            weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')

            print('INFO: Assembling phenotype matrix coordinates.')

            xcounter = 0
            for x in range(0, len(phenotype_list)): #len(phenotype_list)
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

            print('DONE!')
        return

    def merge_matrices(self): #, input_matrix_1, input_matrix_2, output_matrix
        """
        """
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
        """
        """

        distance_matrix = numpy.load('inter/phenolog_gene_cand/distance_matrix.npy')
        weight_matrix = numpy.load('inter/phenolog_gene_cand/weight_matrix.npy')
        ortholog_phenotype_matrix = numpy.load('inter/phenolog_gene_cand/ortholog_phenotype_matrix.npy')


        with open('inter/phenolog_gene_cand/phenotype_list.txt', 'rb') as handle:
            phenotype_list = pickle.load(handle)
        with open('inter/phenolog_gene_cand/ortholog_list.txt', 'rb') as handle:
            ortholog_list = pickle.load(handle)
        with open('inter/hpo/human_pheno_gene_hash.txt', 'rb') as handle:
            human_phenotype_gene_hash = pickle.load(handle)
        with open('inter/mgi/mouse_pheno_gene_hash.txt', 'rb') as handle:
            mouse_phenotype_gene_hash = pickle.load(handle)
        with open('inter/zfin/zebrafish_pheno_gene_hash.txt', 'rb') as handle:
            zebrafish_phenotype_gene_hash = pickle.load(handle)
        phenotype_gene_hash = {}
        phenotype_gene_hash.update(human_phenotype_gene_hash)
        phenotype_gene_hash.update(mouse_phenotype_gene_hash)
        phenotype_gene_hash.update(zebrafish_phenotype_gene_hash)


        with open('inter/hpo/human_phenotype_id_to_label_hash.txt', 'rb') as handle:
            human_phenotype_id_to_label_hash = pickle.load(handle)
        with open('inter/mgi/mouse_phenotype_id_to_label_hash.txt', 'rb') as handle:
            mouse_phenotype_id_to_label_hash = pickle.load(handle)
        with open('inter/zfin/zebrafish_phenotype_id_to_label_hash.txt', 'rb') as handle:
            zebrafish_phenotype_id_to_label_hash = pickle.load(handle)
        phenotype_id_to_label_hash = {}
        phenotype_id_to_label_hash.update(human_phenotype_id_to_label_hash)
        phenotype_id_to_label_hash.update(mouse_phenotype_id_to_label_hash)
        phenotype_id_to_label_hash.update(zebrafish_phenotype_id_to_label_hash)
        nearest_neighbor_hash = {}


        phenotype_ortholog_prediction_matrix = numpy.zeros((len(phenotype_list), len(ortholog_list)))

        # Want to get the 10 nearest neighbors for a given phenotype.
        # Pass in phenotype indice
        # Get the slice for the phenotype indice.
        # Create a clone of the slice, sort, and get the value of the top k entries.
        with open('out/phenolog_gene_cand/nearest_neighbor_phenotypes.txt', 'w', newline='\n') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            for y in range(0, len(phenotype_list)):

                test_phenotype_id = phenotype_list[y]
                print('Test phenotype: '+test_phenotype_id)
                if re.match('HP:.*',test_phenotype_id):
                    phenotype_filter = 'HP'
                elif re.match('MP:.*',test_phenotype_id):
                    phenotype_filter = 'MP'
                elif re.match('ZP:.*',test_phenotype_id):
                    phenotype_filter = 'ZP'
                else:
                    print('ERROR: Unknown phenotype prefix for '+str(test_phenotype_id)+'.')
                    break

                test_distance_slice = distance_matrix[y]

                # The following code will set distance values to zero in the test distance slice
                # if the matching phenotype is from the same species as the test phenotype.
                for x in range(0, len(phenotype_list)):
                    if x != y:
                        match_phenotype_id = phenotype_list[x]
                        match_phenotype_prefix = re.sub(':.*', '', match_phenotype_id)
                        if phenotype_filter == match_phenotype_prefix:
                            test_distance_slice[x] = -1

                intermediate_nearest_neighbors = heapq.nlargest(11, range(len(test_distance_slice)), test_distance_slice.take)
                phenotype_id = phenotype_list[y]
                phenotype_label = phenotype_id_to_label_hash[phenotype_id]
                nearest_neighbors = []
                nearest_neighbor_ids = []
                nearest_neighbor_labels = []
                # This will print out the phenotype IDs for the k nearest neighbors
                for z in intermediate_nearest_neighbors:
                    if z != y:
                        nearest_neighbors.append(z)
                        nearest_neighbor_ids.append(phenotype_list[z])
                for i in nearest_neighbor_ids:
                    #nearest_neighbor_label = 0
                    nearest_neighbor_labels.append(phenotype_id_to_label_hash[i])
                print('Input phenotype: '+phenotype_list[y])
                print('Nearest neighbor phenotypes: '+str(nearest_neighbor_ids)+'.')
                print(nearest_neighbors)
                #for i in nearest_neighbor_ids:
                    #nearest_neighbor_labels.append(phenotype_id_to_label_hash[i])


                # For nearest neighbor output file: phenotype_id, phenotype_label, nn-phenotype_ids, nn-phenotype_labels

                if phenotype_id not in nearest_neighbor_hash:
                    nearest_neighbor_hash[phenotype_id] = nearest_neighbor_ids
                output_row = (phenotype_id, phenotype_label, nearest_neighbor_ids, nearest_neighbor_labels)
                csvwriter.writerow(output_row)

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

        with open('inter/phenolog_gene_cand/nearest_neighbor_hash.txt', 'wb') as handle:
            pickle.dump(nearest_neighbor_hash, handle)
        numpy.save('inter/phenolog_gene_cand/phenotype_ortholog_prediction_matrix.npy', phenotype_ortholog_prediction_matrix)
        numpy.savetxt('inter/phenolog_gene_cand/phenotype_ortholog_prediction_matrix.txt', phenotype_ortholog_prediction_matrix)

        return

    def assemble_phenolog_gene_candidate_predictions_for_phenotypes(self):
        """

        """
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
            phenotype_index_counter = 0
            for phenotype_array in phenotype_ortholog_prediction_matrix:
                ortholog_predictions = []
                ortholog_index_counter = 0
                phenotype_id = phenotype_list[phenotype_index_counter]
                phenotype_ortholog_candidate_hash[phenotype_id] = {}
                for additive_probability in phenotype_array:

                    if additive_probability != 0:
                        ortholog = ortholog_list[ortholog_index_counter]
                        ortholog_predictions.append(ortholog)
                        phenotype_ortholog_candidate_hash[phenotype_id][ortholog] = additive_probability
                    ortholog_index_counter += 1
                print(phenotype_list[phenotype_index_counter])
                print(ortholog_predictions)
                phenotype_index_counter += 1

                output_row = (phenotype_id, ortholog_predictions, phenotype_ortholog_candidate_hash[phenotype_id])
                csvwriter.writerow(output_row)

        with open('out/phenolog_gene_cand/phenolog_ortholog_candidate_prediction_hash.txt', 'wb') as handle:
            pickle.dump(phenotype_ortholog_candidate_hash, handle)

        return
    #TODO: Write function for ranking gene candidate predictions using phenotype_ortholog_candidate_hash.
    #TODO: Write function for assembling multiple gene candidate predictions for disease/genotype phenotypic profiles.

    def assemble_model_level_phenolog_gene_candidate_predictions(self):
        """
        This function assembles phenolog gene candidate predictions for single phenotypes into a collection of
        gene candidate predictions for a disease/genotype by assembling all gene candidates from the phenotypes
        present in the disease/genotype phenotypic profile.
        """
        raw = 'out/phenolog_gene_cand/phenolog_ortholog_candidate_prediction_hash.txt'
        with open(raw, 'rb') as handle:
            phenolog_ortholog_prediction_hash = pickle.load(handle)

        human_file = 'inter/hpo/human_disease_phenotype_hash.txt'
        with open(human_file, 'rb') as handle:
            human_disease_phenotype_hash = pickle.load(handle)
        mouse_file = 'inter/mgi/mouse_genotype_phenotype_hash.txt'
        with open(mouse_file, 'rb') as handle:
            mouse_genotype_phenotype_hash = pickle.load(handle)
        zebrafish_file = 'inter/zfin/zebrafish_genotype_phenotype_hash.txt'
        with open(zebrafish_file, 'rb') as handle:
            zebrafish_genotype_phenotype_hash = pickle.load(handle)

        human_outfile = 'out/phenolog_gene_cand/human_disease_ortholog_candidates.txt'
        mouse_outfile = 'out/phenolog_gene_cand/mouse_genotype_ortholog_candidates.txt'
        zebrafish_outfile = 'out/phenolog_gene_cand/zebrafish_genotype_ortholog_candidates.txt'

        with open(human_outfile, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            for i in human_disease_phenotype_hash:
                disease_id = i
                phenotype_ids = human_disease_phenotype_hash[disease_id]
                phenotype_gene_candidate_hash = {}
                phenotype_counter = 0

                for j in phenotype_ids:
                    try:
                        if phenotype_counter == 0:
                            phenotype_gene_candidate_hash = phenolog_ortholog_prediction_hash[j]
                            phenotype_counter += 1
                        else:
                            phenotype_gene_candidate_hash.update(phenolog_ortholog_prediction_hash[j])
                    except:
                        print('No gene candidate predictions for phenotype '+str(j)+'.')

                phenotype_gene_candidate_hash = OrderedDict(sorted(phenotype_gene_candidate_hash.items(), reverse=True, key=lambda t: t[1]))
                output_row = (disease_id, phenotype_gene_candidate_hash)
                csvwriter.writerow(output_row)

        with open(mouse_outfile, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            for i in mouse_genotype_phenotype_hash:
                genotype_id = i
                phenotype_ids = mouse_genotype_phenotype_hash[genotype_id]
                phenotype_gene_candidate_hash = {}
                phenotype_counter = 0

                for j in phenotype_ids:
                    try:
                        if phenotype_counter == 0:
                            phenotype_gene_candidate_hash = phenolog_ortholog_prediction_hash[j]
                            phenotype_counter += 1
                        else:
                            phenotype_gene_candidate_hash.update(phenolog_ortholog_prediction_hash[j])
                    except:
                        print('No gene candidate predictions for phenotype '+str(j)+'.')

                phenotype_gene_candidate_hash = OrderedDict(sorted(phenotype_gene_candidate_hash.items(), reverse=True, key=lambda t: t[1]))
                output_row = (genotype_id, phenotype_gene_candidate_hash)
                csvwriter.writerow(output_row)

        with open(zebrafish_outfile, 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            for i in zebrafish_genotype_phenotype_hash:
                genotype_id = i
                phenotype_ids = zebrafish_genotype_phenotype_hash[genotype_id]
                phenotype_gene_candidate_hash = {}
                phenotype_counter = 0

                for j in phenotype_ids:
                    try:
                        if phenotype_counter == 0:
                            phenotype_gene_candidate_hash = phenolog_ortholog_prediction_hash[j]
                            phenotype_counter += 1
                        else:
                            phenotype_gene_candidate_hash.update(phenolog_ortholog_prediction_hash[j])
                    except:
                        print('No gene candidate predictions for phenotype '+str(j)+'.')

                phenotype_gene_candidate_hash = OrderedDict(sorted(phenotype_gene_candidate_hash.items(), reverse=True, key=lambda t: t[1]))
                output_row = (genotype_id, phenotype_gene_candidate_hash)
                csvwriter.writerow(output_row)

        #print(phenolog_ortholog_prediction_hash['HP:0000002'])
        return

    def assemble_nearest_neighbor_phenotypes_hash(self):
        # NOTE: This results in a string-based list, and not the python list expected.
        nearest_neighbor_hash = {}

        with open('out/phenolog_gene_cand/nearest_neighbor_phenotypes.txt', 'r', newline='\n') as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (phenotype_id, phenotype_label, nearest_neighbor_ids, nearest_neighbor_labels) = row
                if re.match('HP:.*', phenotype_id):
                    if phenotype_id not in nearest_neighbor_hash:
                        nearest_neighbor_hash[phenotype_id] = nearest_neighbor_ids

        with open('inter/phenolog_gene_cand/nearest_neighbor_hash.txt', 'wb') as handle:
            pickle.dump(nearest_neighbor_hash, handle)
        return

    def assemble_phenolog_gene_candidates_for_diseases(self):

        # Need to take a human disease, grab the associated phenotypes from the disease-phenotype hash,
        # then grab the gene candidate predictions for those phenotypes,
        # combine the gene candidate predictions (selectively update the score if it is greater than the prior score),
        # and then return the top 20 scoring gene candidates.

        # Does it make sense to drop out the genes that are currently associated with the disease of interest?

        phenolog_disease_list = []

        with open('inter/phenolog_gene_cand/nearest_neighbor_hash.txt', 'rb') as handle:
            nearest_neighbor_hash = pickle.load(handle)
        with open('inter/hpo/human_pheno_gene_hash.txt', 'rb') as handle:
            human_phenotype_gene_hash = pickle.load(handle)
        with open('inter/mgi/mouse_pheno_gene_hash.txt', 'rb') as handle:
            mouse_phenotype_to_gene_hash = pickle.load(handle)
        with open('inter/zfin/zebrafish_pheno_gene_hash.txt', 'rb') as handle:
            zebrafish_phenotype_to_gene_hash = pickle.load(handle)
        phenotype_to_gene_hash = {}
        phenotype_to_gene_hash.update(human_phenotype_gene_hash)
        phenotype_to_gene_hash.update(mouse_phenotype_to_gene_hash)
        phenotype_to_gene_hash.update(zebrafish_phenotype_to_gene_hash)


        with open('inter/zfin/zebrafish_gene_to_ortholog_hash.txt', 'rb') as handle:
            zebrafish_gene_to_ortholog_hash = pickle.load(handle)
        with open('inter/mgi/mouse_gene_to_ortholog_hash.txt', 'rb') as handle:
            mouse_gene_to_ortholog_hash = pickle.load(handle)
        with open('inter/hpo/human_gene_to_ortholog_hash.txt', 'rb') as handle:
            human_gene_to_ortholog_hash = pickle.load(handle)

        gene_to_ortholog_hash = {}
        gene_to_ortholog_hash.update(zebrafish_gene_to_ortholog_hash)
        gene_to_ortholog_hash.update(mouse_gene_to_ortholog_hash)
        gene_to_ortholog_hash.update(human_gene_to_ortholog_hash)

        with open('inter/hpo/human_gene_id_to_label_hash.txt', 'rb') as handle:
            human_gene_id_to_label_hash = pickle.load(handle)
        with open('inter/mgi/mouse_gene_id_to_label_hash.txt', 'rb') as handle:
            mouse_gene_id_to_label_hash = pickle.load(handle)
        with open('inter/zfin/zebrafish_gene_id_to_label_hash.txt', 'rb') as handle:
            zebrafish_gene_id_to_label_hash = pickle.load(handle)
        gene_id_to_label_hash = {}
        gene_id_to_label_hash.update(human_gene_id_to_label_hash)
        gene_id_to_label_hash.update(mouse_gene_id_to_label_hash)
        gene_id_to_label_hash.update(zebrafish_gene_id_to_label_hash)

        raw = 'out/phenolog_gene_cand/phenolog_ortholog_candidate_prediction_hash.txt'
        with open(raw, 'rb') as handle:
            phenolog_ortholog_prediction_hash = pickle.load(handle)

        human_file = 'inter/hpo/human_disease_phenotype_hash.txt'
        with open(human_file, 'rb') as handle:
            human_disease_phenotype_hash = pickle.load(handle)

        disease_subset = read_only_disease_subset #['ORPHANET_904', 'ORPHANET_84', 'ORPHANET_46348', 'OMIM_272120', 'ORPHANET_2812', 'ORPHANET_791', 'ORPHANET_478', 'ORPHANET_110', 'OMIM_614592', 'ORPHANET_1873', 'OMIM_305400'] #['ORPHANET_904', 'ORPHANET_84', 'ORPHANET_46348', 'OMIM_272120', 'ORPHANET_2812', 'ORPHANET_791', 'ORPHANET_478', 'ORPHANET_110', 'OMIM_614592', 'ORPHANET_1873', 'OMIM_305400']

        for i in disease_subset:
            try:
                i = re.sub(':', '_', i)
                print('Processing disease ID '+str(i)+'.')
                human_max_outfile = 'out/phenolog_gene_cand/human_disease_gene_candidate_predictions/top_twenty_genes_max/'+str(i)+'.txt'
                human_total_max_outfile = 'out/phenolog_gene_cand/human_disease_gene_candidate_predictions/all_genes/max_scores/'+str(i)+'.txt'

                human_additive_outfile = 'out/phenolog_gene_cand/human_disease_gene_candidate_predictions/top_twenty_genes_additive/'+str(i)+'.txt'
                human_total_additive_outfile = 'out/phenolog_gene_cand/human_disease_gene_candidate_predictions/all_genes/additive_scores/'+str(i)+'.txt'
                disease_id = re.sub('_', ':', i)

                phenotype_ids = human_disease_phenotype_hash[disease_id]
                #print(phenotype_ids)
                max_phenotype_gene_candidate_hash = {}
                additive_phenotype_gene_candidate_hash = {}
                phenotype_counter = 0
                gene_list = []
                cleared_phenotype_ids = []
                nearest_neighbor_phenotypes = []
                for x in phenotype_ids:
                    try:
                        nearest_neighbor_phenotypes = nearest_neighbor_hash[x]
                        cleared_phenotype_ids.append(x)
                        for y in nearest_neighbor_phenotypes:
                            #print(y)
                            associated_genes = phenotype_to_gene_hash[y]
                            #print(associated_genes)
                            for z in associated_genes:
                                #print(z)
                                if z not in gene_list:
                                    gene_list.append(z)
                    except:

                        print('No nearest neighbor phenotypes for phenotype '+str(x)+'.')
                orthogroup_to_gene_hash = {}
                #print(gene_list)
                for x in gene_list:
                    try:
                        panther_id = gene_to_ortholog_hash[x]
                    except:
                        if re.match('MGI:.*', x):
                            panther_id =  self.get_ortholog(x, 'inter/panther/panther_mouse.txt')
                        elif re.match('ZFIN:.*', x):
                            panther_id =  self.get_ortholog(x, 'inter/panther/panther_zebrafish.txt')
                        elif re.match('NCBIGene:.*', x):
                            panther_id =  self.get_ortholog(x, 'inter/panther/panther_human.txt')
                    if panther_id != 'fail' and panther_id not in orthogroup_to_gene_hash:
                        orthogroup_to_gene_hash[panther_id] = [x]
                    elif panther_id != 'fail' and x not in orthogroup_to_gene_hash[panther_id]:
                        orthogroup_to_gene_hash[panther_id].append(x)
                    #print(panther_id)


                '''
                nearest_neighbor_gene_hash = {}
                nearest_neighbor_gene_to_ortholog_hash = {}
                for x in nearest_neighbor_phenotypes:
                    nearest_neighbor_gene_hash[x] = human_phenotype_gene_hash[x]
                    for y in nearest_neighbor_gene_hash[x]:
                        nearest_neighbor_gene_to_ortholog_hash[y] = self.get_ortholog(y, 'inter/panther/panther_hmz_trio.txt')
                '''

                for j in cleared_phenotype_ids:
                    for k in phenolog_ortholog_prediction_hash[j]:

                        if k not in max_phenotype_gene_candidate_hash:
                            max_phenotype_gene_candidate_hash[k] = phenolog_ortholog_prediction_hash[j][k]
                            additive_phenotype_gene_candidate_hash[k] = phenolog_ortholog_prediction_hash[j][k]
                        else:
                            max_phenotype_gene_candidate_hash[k] = max(max_phenotype_gene_candidate_hash[k], phenolog_ortholog_prediction_hash[j][k])
                            additive_phenotype_gene_candidate_hash[k] = (additive_phenotype_gene_candidate_hash[k]+phenolog_ortholog_prediction_hash[j][k])
                        #try
                        #except:
                            #print('No gene candidate predictions for phenotype '+str(j)+'.')


                with open(human_total_max_outfile, 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for orthogroup_candidate_id in max_phenotype_gene_candidate_hash:
                        gene_candidate_ids = orthogroup_to_gene_hash[orthogroup_candidate_id]
                        gene_candidate_labels = []
                        for x in gene_candidate_ids:
                            gene_candidate_labels.append(gene_id_to_label_hash[x])
                        score = '('+str(round(max_phenotype_gene_candidate_hash[orthogroup_candidate_id], 2))+')'
                        output_row =  (gene_candidate_ids, gene_candidate_labels, score)
                        #output_row =  (gene_candidate_labels, score) #(gene_candidate_ids, gene_candidate_labels, score)
                        csvwriter.writerow(output_row)
                with open(human_total_additive_outfile, 'w', newline='') as additive_csvfile:
                    csvwriter = csv.writer(additive_csvfile, delimiter='\t', quotechar='\"')
                    for orthogroup_candidate_id in additive_phenotype_gene_candidate_hash:
                        gene_candidate_ids = orthogroup_to_gene_hash[orthogroup_candidate_id]
                        gene_candidate_labels = []
                        for x in gene_candidate_ids:
                            gene_candidate_labels.append(gene_id_to_label_hash[x])
                        score = '('+str(round(additive_phenotype_gene_candidate_hash[orthogroup_candidate_id], 2))+')'
                        output_row =  (gene_candidate_ids, gene_candidate_labels, score)
                        #output_row = (gene_candidate_labels, score) #(gene_candidate_ids, gene_candidate_labels, score)
                        csvwriter.writerow(output_row)


                top_20_max = heapq.nlargest(20, max_phenotype_gene_candidate_hash, key=lambda k:max_phenotype_gene_candidate_hash[k])
                top_20_additive = heapq.nlargest(20, additive_phenotype_gene_candidate_hash, key=lambda k:additive_phenotype_gene_candidate_hash[k])

                with open(human_max_outfile, 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for orthogroup_candidate_id in top_20_max:
                        gene_candidate_ids = orthogroup_to_gene_hash[orthogroup_candidate_id]
                        gene_candidate_labels = []
                        for x in gene_candidate_ids:
                            gene_candidate_labels.append(gene_id_to_label_hash[x])
                        score = '('+str(round(max_phenotype_gene_candidate_hash[orthogroup_candidate_id], 2))+')'
                        output_row =  (gene_candidate_ids, gene_candidate_labels, score)
                        #output_row =  (gene_candidate_labels, score) #(gene_candidate_ids, gene_candidate_labels, score)
                        csvwriter.writerow(output_row)
                with open(human_additive_outfile, 'w', newline='') as additive_csvfile:
                    csvwriter = csv.writer(additive_csvfile, delimiter='\t', quotechar='\"')
                    for orthogroup_candidate_id in top_20_additive:
                        gene_candidate_ids = orthogroup_to_gene_hash[orthogroup_candidate_id]
                        gene_candidate_labels = []
                        for x in gene_candidate_ids:
                            gene_candidate_labels.append(gene_id_to_label_hash[x])
                        score = '('+str(round(additive_phenotype_gene_candidate_hash[orthogroup_candidate_id], 2))+')'
                        output_row =  (gene_candidate_ids, gene_candidate_labels, score)
                        #output_row = (gene_candidate_labels, score) #(gene_candidate_ids, gene_candidate_labels, score)
                        csvwriter.writerow(output_row)
                if i not in phenolog_disease_list:
                    phenolog_disease_list.append(i)
            except:
                print('Phenolog gene candidate processing failed for disease ID '+str(i)+'.')

        with open('inter/omim/phenolog_disease_list.txt', 'wb') as handle:
            pickle.dump(phenolog_disease_list, handle)
        return

    def assemble_phenolog_orthogroup_candidates_for_diseases(self):

        # Need to take a human disease, grab the associated phenotypes from the disease-phenotype hash,
        # then grab the gene candidate predictions for those phenotypes,
        # combine the gene candidate predictions (selectively update the score if it is greater than the prior score),
        # and then return the top 20 scoring gene candidates.

        with open('inter/hpo/human_gene_id_to_label_hash.txt', 'rb') as handle:
            human_gene_id_to_label_hash = pickle.load(handle)
        with open('inter/mgi/mouse_gene_id_to_label_hash.txt', 'rb') as handle:
            mouse_gene_id_to_label_hash = pickle.load(handle)
        with open('inter/zfin/zebrafish_gene_id_to_label_hash.txt', 'rb') as handle:
            zebrafish_gene_id_to_label_hash = pickle.load(handle)
        gene_id_to_label_hash = {}
        gene_id_to_label_hash.update(human_gene_id_to_label_hash)
        gene_id_to_label_hash.update(mouse_gene_id_to_label_hash)
        gene_id_to_label_hash.update(zebrafish_gene_id_to_label_hash)

        raw = 'out/phenolog_gene_cand/phenolog_ortholog_candidate_prediction_hash.txt'
        with open(raw, 'rb') as handle:
            phenolog_ortholog_prediction_hash = pickle.load(handle)

        human_file = 'inter/hpo/human_disease_phenotype_hash.txt'
        with open(human_file, 'rb') as handle:
            human_disease_phenotype_hash = pickle.load(handle)

        disease_subset = read_only_disease_subset #['ORPHANET_904', 'ORPHANET_84', 'ORPHANET_46348', 'OMIM_272120', 'ORPHANET_2812', 'ORPHANET_791', 'ORPHANET_478', 'ORPHANET_110', 'OMIM_614592', 'ORPHANET_1873', 'OMIM_305400']

        for i in disease_subset:
            try:
                i = re.sub(':', '_', i)
                print('Processing disease ID '+str(i)+'.')
                human_max_outfile = 'out/phenolog_gene_cand/human_disease_orthogroup_candidate_predictions/top_twenty_orthogroups_max/'+str(i)+'.txt'
                human_additive_outfile = 'out/phenolog_gene_cand/human_disease_orthogroup_candidate_predictions/top_twenty_orthogroups_additive/'+str(i)+'.txt'
                human_total_max_outfile = 'out/phenolog_gene_cand/human_disease_orthogroup_candidate_predictions/all_genes/max_scores/'+str(i)+'.txt'
                human_total_additive_outfile = 'out/phenolog_gene_cand/human_disease_orthogroup_candidate_predictions/all_genes/additive_scores/'+str(i)+'.txt'

                disease_id = re.sub('_', ':', i)

                phenotype_ids = human_disease_phenotype_hash[disease_id]
                max_phenotype_gene_candidate_hash = {}
                additive_phenotype_gene_candidate_hash = {}
                phenotype_counter = 0
                cleared_phenotype_ids = []
                for x in phenotype_ids:
                    try:
                        test_hash = phenolog_ortholog_prediction_hash[x]
                        cleared_phenotype_ids.append(x)
                    except:
                        print('No nearest neighbor phenotype available for phenotype '+str(x)+'.')


                for j in cleared_phenotype_ids:

                    for k in phenolog_ortholog_prediction_hash[j]:

                        if k not in max_phenotype_gene_candidate_hash:
                            max_phenotype_gene_candidate_hash[k] = phenolog_ortholog_prediction_hash[j][k]
                            additive_phenotype_gene_candidate_hash[k] = phenolog_ortholog_prediction_hash[j][k]
                        else:
                            max_phenotype_gene_candidate_hash[k] = max(max_phenotype_gene_candidate_hash[k], phenolog_ortholog_prediction_hash[j][k])
                            additive_phenotype_gene_candidate_hash[k] = (additive_phenotype_gene_candidate_hash[k]+phenolog_ortholog_prediction_hash[j][k])

                with open(human_total_max_outfile, 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in max_phenotype_gene_candidate_hash:
                        #gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        score = '('+str(round(max_phenotype_gene_candidate_hash[gene_candidate_id], 2))+')'
                        output_row = (gene_candidate_id, score)
                        csvwriter.writerow(output_row)
                with open(human_total_additive_outfile, 'w', newline='') as additive_csvfile:
                    csvwriter = csv.writer(additive_csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in additive_phenotype_gene_candidate_hash:
                        #gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        score = '('+str(round(additive_phenotype_gene_candidate_hash[gene_candidate_id], 2))+')'
                        output_row = (gene_candidate_id, score)
                        csvwriter.writerow(output_row)

                top_20_max = heapq.nlargest(20, max_phenotype_gene_candidate_hash, key=lambda k:max_phenotype_gene_candidate_hash[k])
                top_20_additive = heapq.nlargest(20, additive_phenotype_gene_candidate_hash, key=lambda k:additive_phenotype_gene_candidate_hash[k])

                with open(human_max_outfile, 'w', newline='') as csvfile:
                    csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in top_20_max:
                        #gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        score = '('+str(round(max_phenotype_gene_candidate_hash[gene_candidate_id], 2))+')'
                        output_row = (gene_candidate_id, score)
                        csvwriter.writerow(output_row)
                with open(human_additive_outfile, 'w', newline='') as additive_csvfile:
                    csvwriter = csv.writer(additive_csvfile, delimiter='\t', quotechar='\"')
                    for gene_candidate_id in top_20_additive:
                        #gene_candidate_label = gene_id_to_label_hash[gene_candidate_id]
                        score = '('+str(round(additive_phenotype_gene_candidate_hash[gene_candidate_id], 2))+')'
                        output_row = (gene_candidate_id, score)
                        csvwriter.writerow(output_row)
            except:
                print('Phenolog ortholog candidate processing failed for disease ID '+str(i)+'.')

        return

    ####### OMIM ASSERTED GENE PROCESSING #######

    def create_zfin_gene_to_ncbi_hash(self):

        zfin_gene_to_ncbi_gene_hash = {}
        ncbi_gene_to_zfin_gene_hash = {}
        with open('raw/panther/zfin_to_ncbi.csv', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:
                (zfin_gene_id, ncbi_gene_id) = row
                if zfin_gene_id not in zfin_gene_to_ncbi_gene_hash:
                    zfin_gene_to_ncbi_gene_hash[zfin_gene_id] = ncbi_gene_id
                if ncbi_gene_id not in ncbi_gene_to_zfin_gene_hash:
                    ncbi_gene_to_zfin_gene_hash[ncbi_gene_id] = zfin_gene_id

        with open('inter/zfin/zfin_gene_to_ncbi_gene_hash.txt', 'wb') as handle:
            pickle.dump(zfin_gene_to_ncbi_gene_hash, handle)
        with open('inter/zfin/ncbi_gene_to_zfin_gene_hash.txt', 'wb') as handle:
            pickle.dump(ncbi_gene_to_zfin_gene_hash, handle)

        with open('inter/zfin/zebrafish_gene_id_to_label_hash.txt', 'rb') as handle:
            zebrafish_gene_id_to_label_hash = pickle.load(handle)
        for gene_id in zebrafish_gene_id_to_label_hash:
            try:
                converted_id = zfin_gene_to_ncbi_gene_hash[gene_id]
            except:
                print('No NCBIGene ID for '+str(gene_id)+'.')
        return

    def create_mgi_gene_to_ncbi_hash(self):

        mgi_gene_to_ncbi_gene_hash = {}
        ncbi_gene_to_mgi_gene_hash = {}
        with open('raw/panther/mgi_to_ncbi.csv', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:
                (mgi_gene_id, ncbi_gene_id) = row
                if mgi_gene_id not in mgi_gene_to_ncbi_gene_hash:
                    mgi_gene_to_ncbi_gene_hash[mgi_gene_id] = ncbi_gene_id
                if ncbi_gene_id not in ncbi_gene_to_mgi_gene_hash:
                    ncbi_gene_to_mgi_gene_hash[ncbi_gene_id] = mgi_gene_id

        with open('inter/mgi/mgi_gene_to_ncbi_gene_hash.txt', 'wb') as handle:
            pickle.dump(mgi_gene_to_ncbi_gene_hash, handle)
        with open('inter/mgi/ncbi_gene_to_mgi_gene_hash.txt', 'wb') as handle:
            pickle.dump(ncbi_gene_to_mgi_gene_hash, handle)

        return

    def create_mim_to_gene_hash(self):

        mim_to_ncbi_gene_hash = {}
        ncbi_gene_to_mim_hash = {}
        with open('raw/omim/mim2gene.txt', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader, None)
            for row in filereader:
                (mim_number, mim_entry_type, ncbi_gene_id, gene_symbol, ensembl_gene_id) = row
                if mim_number not in mim_to_ncbi_gene_hash:
                    mim_to_ncbi_gene_hash[mim_number] = ncbi_gene_id
                if ncbi_gene_id not in ncbi_gene_to_mim_hash:
                    ncbi_gene_to_mim_hash[ncbi_gene_id] = mim_number

        with open('inter/omim/mim_to_ncbi_gene_hash.txt', 'wb') as handle:
            pickle.dump(mim_to_ncbi_gene_hash, handle)
        with open('inter/omim/ncbi_gene_to_mim_hash.txt', 'wb') as handle:
            pickle.dump(ncbi_gene_to_mim_hash, handle)

        return

    def convert_morbid_map_to_ncbi_gene(self):

        with open('inter/omim/mim_to_ncbi_gene_hash.txt', 'rb') as handle:
            mim_to_ncbi_gene_hash = pickle.load(handle)
        disorder_list = []
        disorder_gene_association_id_list = []
        morbid_disease_to_gene_hash = {}
        with open('inter/omim/morbid_disease_to_gene.csv', 'w', encoding="iso-8859-1") as csvfile1:
            csvwriter = csv.writer(csvfile1, delimiter='\t', quotechar='\"')
            with open('raw/omim/morbidmap', 'r', encoding="iso-8859-1") as csvfile2:
                filereader = csv.reader(csvfile2, delimiter='|', quotechar='\"')
                next(filereader, None)
                for row in filereader:
                    #print(row)
                    (disorder_number, gene_symbol, gene_mim_number, cytogenetic_location) = row
                    #print(disorder_id)
                    #Need to properly parse the disorder ID field.
                    p = re.findall(r'(\d{6})', disorder_number)
                    if len(p) == 1:
                        disorder_id = 'OMIM:'+str(p[0])
                        #print(disorder_id)
                        try:
                            gene_id = 'NCBIGene:'+str(mim_to_ncbi_gene_hash[gene_mim_number])
                            disease_gene_association_id = disorder_id+'_'+gene_id
                            output_row = (disease_gene_association_id, disorder_id, gene_id)
                            csvwriter.writerow(output_row)
                            if disease_gene_association_id not in disorder_gene_association_id_list:
                                disorder_gene_association_id_list.append(disease_gene_association_id)
                            if disorder_id not in disorder_list:
                                disorder_list.append(disorder_id)
                        except:
                            print('No NCBIGene ID for '+str(gene_mim_number)+'.')


                    elif len(p) == 0:
                        continue
                    else:
                        print('More than one ID number.'+str(p))
                        for disorder in p:
                            disorder_id = 'OMIM:'+str(disorder)
                            try:
                                gene_id = 'NCBIGene:'+str(mim_to_ncbi_gene_hash[gene_mim_number])
                                disease_gene_association_id = disorder_id+'_'+gene_id
                                output_row = (disease_gene_association_id, disorder_id, gene_id)
                                csvwriter.writerow(output_row)
                                if disease_gene_association_id not in disorder_gene_association_id_list:
                                    disorder_gene_association_id_list.append(disease_gene_association_id)
                                if disorder_id not in disorder_list:
                                    disorder_list.append(disorder_id)
                            except:
                                print('No NCBIGene ID for '+str(gene_mim_number)+'.')
        with open('inter/omim/disorder_list.txt', 'wb') as handle:
            pickle.dump(disorder_list, handle)
        with open('inter/omim/disorder_gene_association_id_list.txt', 'wb') as handle:
            pickle.dump(disorder_gene_association_id_list, handle)

        return

    def convert_morbid_map_genes_to_orthologs(self):

        with open('inter/hpo/human_gene_to_ortholog_hash.txt', 'rb') as handle:
            human_gene_to_ortholog_hash = pickle.load(handle)

        disease_ortholog_association_id_list = []
        with open('inter/omim/morbid_disease_to_ortholog.csv', 'w', newline='') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t', quotechar='\"')
            with open('inter/omim/morbid_disease_to_gene.csv', 'r', encoding="iso-8859-1") as csvfile:
                filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
                for row in filereader:
                    (disease_gene_association_id, disorder_id, gene_id) = row
                    try:
                        ortholog_id = human_gene_to_ortholog_hash[gene_id]
                        if ortholog_id == 'fail':
                            continue
                        else:
                            disease_ortholog_association_id = disorder_id+'_'+ortholog_id
                            output_row = (disease_ortholog_association_id, disorder_id, ortholog_id)
                            csvwriter.writerow(output_row)
                            if disease_ortholog_association_id not in disease_ortholog_association_id_list:
                                disease_ortholog_association_id_list.append(disease_ortholog_association_id)
                    except:
                        try:
                            ortholog_id = self.get_ortholog(gene_id,'inter/panther/panther_human.txt')
                            if ortholog_id == 'fail':
                                continue
                            else:
                                if gene_id not in human_gene_to_ortholog_hash:
                                    human_gene_to_ortholog_hash[gene_id] = ortholog_id
                                disease_ortholog_association_id = disorder_id+'_'+ortholog_id
                                output_row = (disease_ortholog_association_id, disorder_id, ortholog_id)
                                csvwriter.writerow(output_row)
                                if disease_ortholog_association_id not in disease_ortholog_association_id_list:
                                    disease_ortholog_association_id_list.append(disease_ortholog_association_id)
                        except:
                            print('No ortholog found for gene ID '+gene_id+'.')

        with open('inter/omim/disorder_ortholog_association_id_list.txt', 'wb') as handle:
            pickle.dump(disease_ortholog_association_id_list, handle)
        with open('inter/hpo/human_gene_to_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(human_gene_to_ortholog_hash, handle)

        return

    def assemble_complete_disease_list(self):
        with open('inter/omim/owlsim_disease_list.txt', 'rb') as handle:
            owlsim_disease_list = pickle.load(handle)
        print(str(len(owlsim_disease_list)))
        with open('inter/omim/phenolog_disease_list.txt', 'rb') as handle:
            phenolog_disease_list = pickle.load(handle)
        print(str(len(phenolog_disease_list)))
        #common_diseases = set(owlsim_disease_list) ^ set(phenolog_disease_list)
        common_diseases = []
        for disease_id in owlsim_disease_list:
            if disease_id in phenolog_disease_list and disease_id not in common_diseases:
                common_diseases.append(disease_id)
        for disease_id in phenolog_disease_list:
            if disease_id in owlsim_disease_list and disease_id not in common_diseases:
                common_diseases.append(disease_id)
        print(str(len(common_diseases)))
        print(common_diseases)
        with open('inter/omim/common_disease_list.txt', 'wb') as handle:
            pickle.dump(common_diseases, handle)
        return

    def assemble_ortholog_disease_list(self):
        with open('inter/omim/owlsim_disease_list.txt', 'rb') as handle:
            owlsim_disease_list = pickle.load(handle)
        print(str(len(owlsim_disease_list)))
        with open('inter/omim/phenolog_disease_list.txt', 'rb') as handle:
            phenolog_disease_list = pickle.load(handle)
        print(str(len(phenolog_disease_list)))
        #common_diseases = set(owlsim_disease_list) ^ set(phenolog_disease_list)
        common_diseases = []
        for disease_id in owlsim_disease_list:
            if disease_id in phenolog_disease_list and disease_id not in common_diseases:
                common_diseases.append(disease_id)
        for disease_id in phenolog_disease_list:
            if disease_id in owlsim_disease_list and disease_id not in common_diseases:
                common_diseases.append(disease_id)
        print(str(len(common_diseases)))
        print(common_diseases)
        with open('inter/omim/common__ortholog_disease_list.txt', 'wb') as handle:
            pickle.dump(common_diseases, handle)
        return

    def create_ortholog_lookup_hashes(self):

        human_to_mouse_ldo_hash = {}
        human_to_mouse_ortholog_hash = {}
        human_to_zebrafish_ldo_hash = {}
        human_to_zebrafish_ortholog_hash = {}
        mouse_to_human_ldo_hash = {}
        mouse_to_human_ortholog_hash = {}
        mouse_to_zebrafish_ldo_hash = {}
        mouse_to_zebrafish_ortholog_hash = {}
        zebrafish_to_human_ldo_hash = {}
        zebrafish_to_human_ortholog_hash = {}
        zebrafish_to_mouse_ldo_hash = {}
        zebrafish_to_mouse_ortholog_hash = {}
        row_count = 0
        total_row_count = 0
        with open('inter/panther/panther_hmz_trio.txt', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                total_row_count += 1
        print(str(total_row_count)+' rows to process.')
        with open('inter/panther/panther_hmz_trio.txt', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                row_count += 1
                print('Processing row '+str(row_count)+' out of '+str(total_row_count)+'.')
                (panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,proteina,
                 panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,gene_label_b, proteinb,
                 orthology_class, orthology_class_label, panther_id) = row


                if taxon_id_a == 'NCBITaxon:9606':
                    if taxon_id_b == 'NCBITaxon:10090':
                        if orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                human_to_mouse_ldo_hash[gene_id_a] = gene_id_b
                                mouse_to_human_ldo_hash[gene_id_b] = gene_id_a
                        if orthology_class == 'O' or orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                if gene_id_a not in human_to_mouse_ortholog_hash:
                                    human_to_mouse_ortholog_hash[gene_id_a] = [gene_id_b]
                                else:
                                    human_to_mouse_ortholog_hash[gene_id_a].append(gene_id_a)
                                if gene_id_b not in mouse_to_human_ortholog_hash:
                                    mouse_to_human_ortholog_hash[gene_id_b] = [gene_id_a]
                                else:
                                    mouse_to_human_ortholog_hash[gene_id_b].append(gene_id_a)

                    elif taxon_id_b == 'NCBITaxon:7955':
                        if orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                human_to_zebrafish_ldo_hash[gene_id_a] = gene_id_b
                                zebrafish_to_human_ldo_hash[gene_id_b] = gene_id_a
                        if orthology_class == 'O' or orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                if gene_id_a not in human_to_zebrafish_ortholog_hash:
                                    human_to_zebrafish_ortholog_hash[gene_id_a] = [gene_id_b]
                                else:
                                    human_to_zebrafish_ortholog_hash[gene_id_a].append(gene_id_a)
                                if gene_id_b not in zebrafish_to_human_ortholog_hash:
                                    zebrafish_to_human_ortholog_hash[gene_id_b] = [gene_id_a]
                                else:
                                    zebrafish_to_human_ortholog_hash[gene_id_b].append(gene_id_a)
                    else:
                        continue

                elif taxon_id_a == 'NCBITaxon:10090':
                    if taxon_id_b == 'NCBITaxon:7955':
                        if orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                mouse_to_zebrafish_ldo_hash[gene_id_a] = gene_id_b
                                zebrafish_to_mouse_ldo_hash[gene_id_b] = gene_id_a
                        if orthology_class == 'O' or orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                if gene_id_a not in mouse_to_zebrafish_ortholog_hash:
                                    mouse_to_zebrafish_ortholog_hash[gene_id_a] = [gene_id_b]
                                else:
                                    mouse_to_zebrafish_ortholog_hash[gene_id_a].append(gene_id_a)
                                if gene_id_b not in zebrafish_to_mouse_ortholog_hash:
                                    zebrafish_to_mouse_ortholog_hash[gene_id_b] = [gene_id_a]
                                else:
                                    zebrafish_to_mouse_ortholog_hash[gene_id_b].append(gene_id_a)

                    elif taxon_id_b == 'NCBITaxon:9606':
                        if orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                mouse_to_human_ldo_hash[gene_id_a] = gene_id_b
                                human_to_mouse_ldo_hash[gene_id_b] = gene_id_a
                        if orthology_class == 'O' or orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                if gene_id_a not in mouse_to_human_ortholog_hash:
                                    mouse_to_human_ortholog_hash[gene_id_a] = [gene_id_b]
                                else:
                                    mouse_to_human_ortholog_hash[gene_id_a].append(gene_id_a)
                                if gene_id_b not in human_to_mouse_ortholog_hash:
                                    human_to_mouse_ortholog_hash[gene_id_b] = [gene_id_a]
                                else:
                                    human_to_mouse_ortholog_hash[gene_id_b].append(gene_id_a)
                    else:
                        continue

                elif taxon_id_a == 'NCBITaxon:7955':
                    if taxon_id_b == 'NCBITaxon:10090':
                        if orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                zebrafish_to_mouse_ldo_hash[gene_id_a] = gene_id_b
                                mouse_to_zebrafish_ldo_hash[gene_id_b] = gene_id_a
                        if orthology_class == 'O' or orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                if gene_id_a not in zebrafish_to_mouse_ortholog_hash:
                                    zebrafish_to_mouse_ortholog_hash[gene_id_a] = [gene_id_b]
                                else:
                                    zebrafish_to_mouse_ortholog_hash[gene_id_a].append(gene_id_a)
                                if gene_id_b not in mouse_to_zebrafish_ortholog_hash:
                                    mouse_to_zebrafish_ortholog_hash[gene_id_b] = [gene_id_a]
                                else:
                                    mouse_to_zebrafish_ortholog_hash[gene_id_b].append(gene_id_a)

                    elif taxon_id_b == 'NCBITaxon:9606':
                        if orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                zebrafish_to_human_ldo_hash[gene_id_a] = gene_id_b
                                human_to_zebrafish_ldo_hash[gene_id_b] = gene_id_a
                        if orthology_class == 'O' or orthology_class == 'LDO':
                            if re.match('NCBIGene:.*', gene_id_a) and re.match('NCBIGene:.*', gene_id_b):
                                if gene_id_a not in zebrafish_to_human_ortholog_hash:
                                    zebrafish_to_human_ortholog_hash[gene_id_a] = [gene_id_b]
                                else:
                                    zebrafish_to_human_ortholog_hash[gene_id_a].append(gene_id_a)
                                if gene_id_b not in human_to_zebrafish_ortholog_hash:
                                    human_to_zebrafish_ortholog_hash[gene_id_b] = [gene_id_a]
                                else:
                                    human_to_zebrafish_ortholog_hash[gene_id_b].append(gene_id_a)
                    else:
                        continue

                else:
                    continue

        with open('inter/panther/human_to_mouse_ldo_hash.txt', 'wb') as handle:
            pickle.dump(human_to_mouse_ldo_hash, handle)
        with open('inter/panther/human_to_mouse_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(human_to_mouse_ortholog_hash, handle)
        with open('inter/panther/human_to_zebrafish_ldo_hash.txt', 'wb') as handle:
            pickle.dump(human_to_zebrafish_ldo_hash, handle)
        with open('inter/panther/human_to_zebrafish_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(human_to_zebrafish_ortholog_hash, handle)
        with open('inter/panther/mouse_to_human_ldo_hash.txt', 'wb') as handle:
            pickle.dump(mouse_to_human_ldo_hash, handle)
        with open('inter/panther/mouse_to_human_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(mouse_to_human_ortholog_hash, handle)
        with open('inter/panther/zebrafish_to_human_ldo_hash.txt', 'wb') as handle:
            pickle.dump(zebrafish_to_human_ldo_hash, handle)
        with open('inter/panther/zebrafish_to_human_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(zebrafish_to_human_ortholog_hash, handle)
        with open('inter/panther/mouse_to_zebrafish_ldo_hash.txt', 'wb') as handle:
            pickle.dump(mouse_to_zebrafish_ldo_hash, handle)
        with open('inter/panther/mouse_to_zebrafish_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(mouse_to_zebrafish_ortholog_hash, handle)
        with open('inter/panther/zebrafish_to_mouse_ldo_hash.txt', 'wb') as handle:
            pickle.dump(zebrafish_to_mouse_ldo_hash, handle)
        with open('inter/panther/zebrafish_to_mouse_ortholog_hash.txt', 'wb') as handle:
            pickle.dump(zebrafish_to_mouse_ortholog_hash, handle)

        return

    def convert_omim_ncbi_gene_to_mgi_and_zfin(self):

        with open('inter/omim/morbid_disease_to_gene.csv', 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            for row in filereader:
                (panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,proteina,
                 panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,gene_label_b, proteinb,
                 orthology_class, orthology_class_label, panther_id) = row


        return

    def assemble_owlsim_data_for_ROC(self):

        # Have a common disease list for OMIM, OWLSim, and Phenologs.
        # OWLSim and Phenolog gene predictions are in mouse/zebrafish gene IDs.
        # Will need to convert those gene IDs to human NCBI Gene IDs through PANTHER.
        # However, it would be faster to simply convert the OMIM Gene IDs to the corresponding MGI/ZFIN Gene IDs, right?

        with open('inter/omim/common_disease_list.txt', 'rb') as handle:
            common_diseases_list = pickle.load(handle)
        with open('inter/mgi/mgi_gene_to_ncbi_gene_hash.txt', 'rb') as handle:
            mgi_gene_to_ncbi_gene_hash = pickle.load(handle)
        with open('inter/zfin/zfin_gene_to_ncbi_gene_hash.txt', 'rb') as handle:
            zfin_gene_to_ncbi_gene_hash = pickle.load(handle)
        with open('inter/omim/ncbi_gene_to_mim_hash.txt', 'rb') as handle:
            ncbi_gene_to_mim_hash = pickle.load(handle)

        for disease in common_diseases_list:
            try:
                file_prefix = re.sub(':', '_', disease)
                filename = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(file_prefix)+'.txt'
                #print(filename)
                with open(filename, 'rb') as handle:
                    human_disease_gene_prediction_hash = pickle.load(handle)




                    #print(human_disease_gene_prediction_hash)
            except:
                print('No data available for disorder ID '+str(disease)+'.')

        return


    def trim_morbid_disease_to_gene(self):
        with open('inter/omim/morbid_disease_to_gene_trimmed.csv', 'w', newline='') as csvfile1:
            csvwriter = csv.writer(csvfile1, delimiter='\t', quotechar='\"')
            with open('inter/omim/morbid_disease_to_gene.csv', 'r', encoding="iso-8859-1") as csvfile2:
                filereader = csv.reader(csvfile2, delimiter='\t', quotechar='\"')
                for row in filereader:
                    (disease_gene_association_id, disorder_id, gene_id) = row
                    if gene_id == 'NCBIGene:-':
                        continue
                    else:
                        output_row = (disease_gene_association_id, disorder_id, gene_id)
                        csvwriter.writerow(output_row)
        return

    '''

    def assemble_ROC_score_lists(self):

        max_ic_list = []
        iccs_list = []
        sim_ic_list = []
        sim_j_list = []
        phenolog_gene_score_list = []
        phenolog_ortholog_score_list = []

        with open('inter/omim/common_disease_list.txt', 'rb') as handle:
            common_diseases_list = pickle.load(handle)
        with open('inter/mgi/mgi_gene_to_ncbi_gene_hash.txt', 'rb') as handle:
            mgi_gene_to_ncbi_gene_hash = pickle.load(handle)
        with open('inter/mgi/ncbi_gene_to_mgi_gene_hash.txt', 'rb') as handle:
            ncbi_gene_to_mgi_gene_hash = pickle.load(handle)
        with open('inter/zfin/zfin_gene_to_ncbi_gene_hash.txt', 'rb') as handle:
            zfin_gene_to_ncbi_gene_hash = pickle.load(handle)
        with open('inter/zfin/ncbi_gene_to_zfin_gene_hash.txt', 'rb') as handle:
            ncbi_gene_to_zfin_gene_hash = pickle.load(handle)
        with open('inter/omim/ncbi_gene_to_mim_hash.txt', 'rb') as handle:
            ncbi_gene_to_mim_hash = pickle.load(handle)
        with open('inter/omim/mim_to_ncbi_gene_hash.txt', 'rb') as handle:
            mim_to_ncbi_gene_hash = pickle.load(handle)


        with open('inter/panther/human_to_mouse_ldo_hash.txt', 'rb') as handle:
            human_to_mouse_ldo_hash = pickle.load(handle)
        with open('inter/panther/human_to_mouse_ortholog_hash.txt', 'rb') as handle:
            human_to_mouse_ortholog_hash = pickle.load(handle)
        with open('inter/panther/human_to_zebrafish_ldo_hash.txt', 'rb') as handle:
            human_to_zebrafish_ldo_hash = pickle.load(handle)
        with open('inter/panther/human_to_zebrafish_ortholog_hash.txt', 'rb') as handle:
            human_to_zebrafish_ortholog_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_human_ldo_hash.txt', 'rb') as handle:
            mouse_to_human_ldo_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_human_ortholog_hash.txt', 'rb') as handle:
            mouse_to_human_ortholog_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_human_ldo_hash.txt', 'rb') as handle:
            zebrafish_to_human_ldo_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_human_ortholog_hash.txt', 'rb') as handle:
            zebrafish_to_human_ortholog_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_zebrafish_ldo_hash.txt', 'rb') as handle:
            mouse_to_zebrafish_ldo_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_zebrafish_ortholog_hash.txt', 'rb') as handle:
            mouse_to_zebrafish_ortholog_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_mouse_ldo_hash.txt', 'rb') as handle:
            zebrafish_to_mouse_ldo_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_mouse_ortholog_hash.txt', 'rb') as handle:
            zebrafish_to_mouse_ortholog_hash = pickle.load(handle)

        with open('inter/omim/morbid_disease_predictions.csv', 'w') as csvfile1:
            csvwriter = csv.writer(csvfile1, delimiter='\t', quotechar='\"')
            with open('inter/omim/morbid_disease_to_gene_trimmed.csv', 'r', encoding="iso-8859-1") as csvfile2:
                filereader = csv.reader(csvfile2, delimiter='\t', quotechar='\"')
                for row in filereader:
                    (disease_gene_association_id, disease_id, gene_id) = row
                    print(disease_id)
                    disease_id_underscored = re.sub(':', '_', disease_id)
                    if disease_id_underscored not in common_diseases_list:
                        print('Disease ID '+str(disease_id)+' not in common disease list.')
                        continue
                    else:
                        zebrafish_max_ic = ''
                        zebrafish_iccs = ''
                        zebrafish_sim_ic = ''
                        zebrafish_sim_j = ''

                        mouse_max_ic = ''
                        mouse_iccs = ''
                        mouse_sim_ic = ''
                        mouse_sim_j = ''


                        try:
                            filename = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(disease_id_underscored)+'.txt'
                            #print(filename)
                            with open(filename, 'rb') as handle:
                                human_disease_gene_prediction_hash = pickle.load(handle)
                            print('OWLSim file open.')
                            try:
                                zebrafish_ncbi_gene_id = human_to_zebrafish_ldo_hash[gene_id]
                                print('Zebrafish LDO found.')
                                zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                print('ZFIN gene ID: '+str(zfin_gene_id))
                                print(human_disease_gene_prediction_hash[zfin_gene_id])
                                zebrafish_max_ic = human_disease_gene_prediction_hash[zfin_gene_id]['maxIC']
                                zebrafish_iccs = human_disease_gene_prediction_hash[zfin_gene_id]['ICCS']
                                zebrafish_sim_ic = human_disease_gene_prediction_hash[zfin_gene_id]['simIC']
                                zebrafish_sim_j = human_disease_gene_prediction_hash[zfin_gene_id]['simJ']
                            except:
                                print('No zebrafish OWLSim scores found for LDO.')
                                try:
                                    zebrafish_gene_ids = human_to_zebrafish_ortholog_hash[gene_id]
                                    print('Zebrafish ortholog found.')
                                    if len(zebrafish_gene_ids) == 1:
                                        zebrafish_ncbi_gene_id = zebrafish_gene_ids[0]
                                        zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                        print('ZFIN gene ID: '+str(zfin_gene_id))
                                        print(human_disease_gene_prediction_hash[zfin_gene_id])
                                        zebrafish_max_ic = human_disease_gene_prediction_hash[zfin_gene_id]['maxIC']
                                        zebrafish_iccs = human_disease_gene_prediction_hash[zfin_gene_id]['ICCS']
                                        zebrafish_sim_ic = human_disease_gene_prediction_hash[zfin_gene_id]['simIC']
                                        zebrafish_sim_j = human_disease_gene_prediction_hash[zfin_gene_id]['simJ']
                                    elif len(zebrafish_gene_ids) > 1:
                                        print('More than one ortholog gene ID found!')
                                        # TODO: Need to handle multiple genes.
                                        for x in zebrafish_gene_ids:
                                            print(x)
                                            zfin_gene_id = ncbi_gene_to_zfin_gene_hash[x]
                                            print('ZFIN gene ID: '+str(zfin_gene_id))
                                            print(human_disease_gene_prediction_hash[zfin_gene_id])

                                except:
                                    print('No zebrafish OWLSim scores found.')

                            try:
                                mouse_ncbi_gene_id = human_to_mouse_ldo_hash[gene_id]
                                print('Mouse LDO found.')
                                mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                                print('MGI gene ID: '+str(mgi_gene_id))
                                print(human_disease_gene_prediction_hash[mgi_gene_id])
                                mouse_max_ic = human_disease_gene_prediction_hash[mgi_gene_id]['maxIC']
                                mouse_iccs = human_disease_gene_prediction_hash[mgi_gene_id]['ICCS']
                                mouse_sim_ic = human_disease_gene_prediction_hash[mgi_gene_id]['simIC']
                                mouse_sim_j = human_disease_gene_prediction_hash[mgi_gene_id]['simJ']
                            except:
                                print('No mouse OWLSim scores found for LDO.')
                                try:
                                    mouse_gene_ids = human_to_mouse_ortholog_hash[gene_id]
                                    print('Mouse ortholog found.')
                                    if len(mouse_gene_ids) == 1:
                                        mouse_ncbi_gene_id = mouse_gene_ids[0]
                                        mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                                        print('MGI gene ID: '+str(mgi_gene_id))
                                        print(human_disease_gene_prediction_hash[mgi_gene_id])
                                        mouse_max_ic = human_disease_gene_prediction_hash[mgi_gene_id]['maxIC']
                                        mouse_iccs = human_disease_gene_prediction_hash[mgi_gene_id]['ICCS']
                                        mouse_sim_ic = human_disease_gene_prediction_hash[mgi_gene_id]['simIC']
                                        mouse_sim_j = human_disease_gene_prediction_hash[mgi_gene_id]['simJ']
                                    elif len(mouse_gene_ids) > 1:
                                        print('More than one ortholog gene ID found!')
                                        # TODO: Need to handle multiple genes.
                                        for x in mouse_gene_ids:
                                            print(x)
                                            mgi_gene_id = ncbi_gene_to_mgi_gene_hash[x]
                                            print('MGI gene ID: '+str(mgi_gene_id))
                                            print(human_disease_gene_prediction_hash[mgi_gene_id])
                                except:
                                    print('No mouse OWLSim scores found.')


                        except:
                            print('Trouble with OWLSim data for disease '+str(disease_id)+'.')


                    # TODO: Need to add phenolog data lookup here.

                    phenolog_score = ''

                    try:
                        zebrafish_ncbi_gene_id = human_to_zebrafish_ldo_hash[gene_id]
                        print('Zebrafish LDO found.')
                        zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                        print('ZFIN gene ID: '+str(zfin_gene_id))

                    except:
                        print('No zebrafish LDO found.')
                        zebrafish_gene_ids = human_to_zebrafish_ortholog_hash[gene_id]
                        #print('Zebrafish ortholog found.')
                        if len(zebrafish_gene_ids) == 1:
                            zebrafish_ncbi_gene_id = zebrafish_gene_ids[0]
                    try:
                        zebrafish_gene_ids = human_to_zebrafish_ortholog_hash[gene_id]
                        print('Zebrafish orthologs found.')
                        if len(zebrafish_gene_ids) == 1:
                            zebrafish_ncbi_gene_id = zebrafish_gene_ids[0]
                            zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                            print('ZFIN gene ID: '+str(zfin_gene_id))
                    except:
                        print('No zebrafish orthologs found.')


                    try:
                        mouse_ncbi_gene_id = human_to_mouse_ldo_hash[gene_id]
                        print('Mouse LDO found.')
                        mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                        print('MGI gene ID: '+str(mgi_gene_id))
                    except:
                        print('No mouse LDO found.')
                    try:
                        mouse_gene_ids = human_to_mouse_ortholog_hash[gene_id]
                        print('Mouse orthologs found.')
                        if len(mouse_gene_ids) == 1:
                            mouse_ncbi_gene_id = mouse_gene_ids[0]
                            mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                            print('MGI gene ID: '+str(mgi_gene_id))
                    except:
                        print('No mouse orthologs found.')





                        filename = 'out/phenolog_gene_cand/human_disease_gene_candidate_predictions/all_genes/max_scores/'+str(disease_id_underscored)+'.txt'
                        #print(filename)

                        with open(filename, 'r', encoding="iso-8859-1") as csvfile2:
                            filereader = csv.reader(csvfile2, delimiter='\t', quotechar='\"')
                            for row in filereader:
                                (gene_candidate_ids, gene_candidate_labels, score) = row
                                if len(gene_candidate_ids) == 1:
                                    gene_id = gene_candidate_ids[0]
                                    try:
                                        zebrafish_ncbi_gene_id = human_to_zebrafish_ldo_hash[gene_id]
                                        print('Zebrafish LDO found.')
                                        zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                        print('ZFIN gene ID: '+str(zfin_gene_id))
                                        print(human_disease_gene_prediction_hash[zfin_gene_id])
                                        zebrafish_phenolog_score = ''
                                        zebrafish_max_ic = human_disease_gene_prediction_hash[zfin_gene_id]['maxIC']
                                        zebrafish_iccs = human_disease_gene_prediction_hash[zfin_gene_id]['ICCS']
                                        zebrafish_sim_ic = human_disease_gene_prediction_hash[zfin_gene_id]['simIC']
                                        zebrafish_sim_j = human_disease_gene_prediction_hash[zfin_gene_id]['simJ']



                                if len(gene_candidate_ids) > 1:
                                    for gene_id in gene_candidate_ids:
                                        try:

                                        except:
                                elif len(gene_candidate_ids) == 0:
                                    continue



                        print('Phenolog file open.')
                    except:
                        print('Trouble with Phenolog data for disease '+str(disease_id)+'.')


                    output_row = (disease_gene_association_id, disease_id, gene_id,
                                  zebrafish_max_ic, zebrafish_iccs, zebrafish_sim_ic, zebrafish_sim_j,
                                  mouse_max_ic, mouse_iccs, mouse_sim_ic, mouse_sim_j)
                    csvwriter.writerow(output_row)



        with open('out/roc/max_ic_list.txt', 'wb') as handle:
            pickle.dump(max_ic_list, handle)
        with open('out/roc/iccs_list.txt', 'wb') as handle:
            pickle.dump(iccs_list, handle)
        with open('out/roc/sim_ic_list.txt', 'wb') as handle:
            pickle.dump(sim_ic_list, handle)
        with open('out/roc/sim_j_list.txt', 'wb') as handle:
            pickle.dump(sim_j_list, handle)
        with open('out/roc/phenolog_gene_score_list.txt', 'wb') as handle:
            pickle.dump(phenolog_gene_score_list, handle)
        with open('out/roc/phenolog_ortholog_score_list.txt', 'wb') as handle:
            pickle.dump(phenolog_ortholog_score_list, handle)

        return

    '''

    def assemble_ROC_score_lists_alternate(self):

        max_ic_list = []
        iccs_list = []
        sim_ic_list = []
        sim_j_list = []
        phenolog_gene_score_list = []
        phenolog_ortholog_score_list = []

        with open('inter/omim/common_disease_list.txt', 'rb') as handle:
            common_diseases_list = pickle.load(handle)
        with open('inter/mgi/mgi_gene_to_ncbi_gene_hash.txt', 'rb') as handle:
            mgi_gene_to_ncbi_gene_hash = pickle.load(handle)
        with open('inter/mgi/ncbi_gene_to_mgi_gene_hash.txt', 'rb') as handle:
            ncbi_gene_to_mgi_gene_hash = pickle.load(handle)
        with open('inter/zfin/zfin_gene_to_ncbi_gene_hash.txt', 'rb') as handle:
            zfin_gene_to_ncbi_gene_hash = pickle.load(handle)
        with open('inter/zfin/ncbi_gene_to_zfin_gene_hash.txt', 'rb') as handle:
            ncbi_gene_to_zfin_gene_hash = pickle.load(handle)
        with open('inter/omim/ncbi_gene_to_mim_hash.txt', 'rb') as handle:
            ncbi_gene_to_mim_hash = pickle.load(handle)
        with open('inter/omim/mim_to_ncbi_gene_hash.txt', 'rb') as handle:
            mim_to_ncbi_gene_hash = pickle.load(handle)


        with open('inter/panther/human_to_mouse_ldo_hash.txt', 'rb') as handle:
            human_to_mouse_ldo_hash = pickle.load(handle)
        with open('inter/panther/human_to_mouse_ortholog_hash.txt', 'rb') as handle:
            human_to_mouse_ortholog_hash = pickle.load(handle)
        with open('inter/panther/human_to_zebrafish_ldo_hash.txt', 'rb') as handle:
            human_to_zebrafish_ldo_hash = pickle.load(handle)
        with open('inter/panther/human_to_zebrafish_ortholog_hash.txt', 'rb') as handle:
            human_to_zebrafish_ortholog_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_human_ldo_hash.txt', 'rb') as handle:
            mouse_to_human_ldo_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_human_ortholog_hash.txt', 'rb') as handle:
            mouse_to_human_ortholog_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_human_ldo_hash.txt', 'rb') as handle:
            zebrafish_to_human_ldo_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_human_ortholog_hash.txt', 'rb') as handle:
            zebrafish_to_human_ortholog_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_zebrafish_ldo_hash.txt', 'rb') as handle:
            mouse_to_zebrafish_ldo_hash = pickle.load(handle)
        with open('inter/panther/mouse_to_zebrafish_ortholog_hash.txt', 'rb') as handle:
            mouse_to_zebrafish_ortholog_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_mouse_ldo_hash.txt', 'rb') as handle:
            zebrafish_to_mouse_ldo_hash = pickle.load(handle)
        with open('inter/panther/zebrafish_to_mouse_ortholog_hash.txt', 'rb') as handle:
            zebrafish_to_mouse_ortholog_hash = pickle.load(handle)

        with open('inter/omim/morbid_disease_predictions.csv', 'w') as csvfile1:
            csvwriter = csv.writer(csvfile1, delimiter='\t', quotechar='\"')
            with open('inter/omim/morbid_disease_to_gene_trimmed.csv', 'r', encoding="iso-8859-1") as csvfile2:
                filereader = csv.reader(csvfile2, delimiter='\t', quotechar='\"')
                for row in filereader:
                    (disease_gene_association_id, disease_id, gene_id) = row
                    print(disease_id)
                    disease_id_underscored = re.sub(':', '_', disease_id)
                    if disease_id_underscored not in common_diseases_list:
                        print('Disease ID '+str(disease_id)+' not in common disease list.')
                        continue
                    else:

                        query_zfin_ldo_gene_id = ''
                        query_zfin_ortholog_gene_ids = ''

                        query_mgi_ldo_gene_id = ''
                        query_mgi_ortholog_gene_ids = ''

                        zebrafish_max_ic = ''
                        zebrafish_iccs = ''
                        zebrafish_sim_ic = ''
                        zebrafish_sim_j = ''

                        mouse_max_ic = ''
                        mouse_iccs = ''
                        mouse_sim_ic = ''
                        mouse_sim_j = ''


                        try:
                            zebrafish_ncbi_gene_id = human_to_zebrafish_ldo_hash[gene_id]
                            print('Zebrafish LDO found.')
                            query_zfin_ldo_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                            print('ZFIN gene ID: '+str(query_zfin_ldo_gene_id))

                        except:
                            print('No zebrafish LDO found.')
                            zebrafish_gene_ids = human_to_zebrafish_ortholog_hash[gene_id]
                            #print('Zebrafish ortholog found.')
                            if len(zebrafish_gene_ids) == 1:
                                query_zfin_ortholog_gene_ids = zebrafish_gene_ids[0]
                        try:
                            zebrafish_gene_ids = human_to_zebrafish_ortholog_hash[gene_id]
                            print('Zebrafish orthologs found.')
                            if len(zebrafish_gene_ids) == 1:
                                zebrafish_ncbi_gene_id = zebrafish_gene_ids[0]
                                zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                print('ZFIN gene ID: '+str(zfin_gene_id))
                        except:
                            print('No zebrafish orthologs found.')


                        try:
                            mouse_ncbi_gene_id = human_to_mouse_ldo_hash[gene_id]
                            print('Mouse LDO found.')
                            query_mgi_ldo_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                            print('MGI gene ID: '+str(query_mgi_ldo_gene_id))
                        except:
                            print('No mouse LDO found.')
                        try:
                            mouse_gene_ids = human_to_mouse_ortholog_hash[gene_id]
                            print('Mouse orthologs found.')
                            if len(mouse_gene_ids) == 1:
                                mouse_ncbi_gene_id = mouse_gene_ids[0]
                                mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                                print('MGI gene ID: '+str(mgi_gene_id))
                        except:
                            print('No mouse orthologs found.')





                        try:
                            filename = 'out/owlsim/human_disease_gene_candidate_predictions/all_genes/'+str(disease_id_underscored)+'.txt'
                            #print(filename)
                            with open(filename, 'rb') as handle:
                                human_disease_gene_prediction_hash = pickle.load(handle)
                            print('OWLSim file open.')
                            try:
                                zebrafish_ncbi_gene_id = human_to_zebrafish_ldo_hash[gene_id]
                                print('Zebrafish LDO found.')
                                zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                print('ZFIN gene ID: '+str(zfin_gene_id))
                                print(human_disease_gene_prediction_hash[zfin_gene_id])
                                zebrafish_max_ic = human_disease_gene_prediction_hash[zfin_gene_id]['maxIC']
                                zebrafish_iccs = human_disease_gene_prediction_hash[zfin_gene_id]['ICCS']
                                zebrafish_sim_ic = human_disease_gene_prediction_hash[zfin_gene_id]['simIC']
                                zebrafish_sim_j = human_disease_gene_prediction_hash[zfin_gene_id]['simJ']
                            except:
                                print('No zebrafish OWLSim scores found for LDO.')
                                try:
                                    zebrafish_gene_ids = human_to_zebrafish_ortholog_hash[gene_id]
                                    print('Zebrafish ortholog found.')
                                    if len(zebrafish_gene_ids) == 1:
                                        zebrafish_ncbi_gene_id = zebrafish_gene_ids[0]
                                        zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                        print('ZFIN gene ID: '+str(zfin_gene_id))
                                        print(human_disease_gene_prediction_hash[zfin_gene_id])
                                        zebrafish_max_ic = human_disease_gene_prediction_hash[zfin_gene_id]['maxIC']
                                        zebrafish_iccs = human_disease_gene_prediction_hash[zfin_gene_id]['ICCS']
                                        zebrafish_sim_ic = human_disease_gene_prediction_hash[zfin_gene_id]['simIC']
                                        zebrafish_sim_j = human_disease_gene_prediction_hash[zfin_gene_id]['simJ']
                                    elif len(zebrafish_gene_ids) > 1:
                                        print('More than one ortholog gene ID found!')
                                        # TODO: Need to handle multiple genes.
                                        for x in zebrafish_gene_ids:
                                            print(x)
                                            zfin_gene_id = ncbi_gene_to_zfin_gene_hash[x]
                                            print('ZFIN gene ID: '+str(zfin_gene_id))
                                            print(human_disease_gene_prediction_hash[zfin_gene_id])

                                except:
                                    print('No zebrafish OWLSim scores found.')

                            try:
                                mouse_ncbi_gene_id = human_to_mouse_ldo_hash[gene_id]
                                print('Mouse LDO found.')
                                mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                                print('MGI gene ID: '+str(mgi_gene_id))
                                print(human_disease_gene_prediction_hash[mgi_gene_id])
                                mouse_max_ic = human_disease_gene_prediction_hash[mgi_gene_id]['maxIC']
                                mouse_iccs = human_disease_gene_prediction_hash[mgi_gene_id]['ICCS']
                                mouse_sim_ic = human_disease_gene_prediction_hash[mgi_gene_id]['simIC']
                                mouse_sim_j = human_disease_gene_prediction_hash[mgi_gene_id]['simJ']
                            except:
                                print('No mouse OWLSim scores found for LDO.')
                                try:
                                    mouse_gene_ids = human_to_mouse_ortholog_hash[gene_id]
                                    print('Mouse ortholog found.')
                                    if len(mouse_gene_ids) == 1:
                                        mouse_ncbi_gene_id = mouse_gene_ids[0]
                                        mgi_gene_id = ncbi_gene_to_mgi_gene_hash[mouse_ncbi_gene_id]
                                        print('MGI gene ID: '+str(mgi_gene_id))
                                        print(human_disease_gene_prediction_hash[mgi_gene_id])
                                        mouse_max_ic = human_disease_gene_prediction_hash[mgi_gene_id]['maxIC']
                                        mouse_iccs = human_disease_gene_prediction_hash[mgi_gene_id]['ICCS']
                                        mouse_sim_ic = human_disease_gene_prediction_hash[mgi_gene_id]['simIC']
                                        mouse_sim_j = human_disease_gene_prediction_hash[mgi_gene_id]['simJ']
                                    elif len(mouse_gene_ids) > 1:
                                        print('More than one ortholog gene ID found!')
                                        # TODO: Need to handle multiple genes.
                                        for x in mouse_gene_ids:
                                            print(x)
                                            mgi_gene_id = ncbi_gene_to_mgi_gene_hash[x]
                                            print('MGI gene ID: '+str(mgi_gene_id))
                                            print(human_disease_gene_prediction_hash[mgi_gene_id])
                                except:
                                    print('No mouse OWLSim scores found.')


                        except:
                            print('Trouble with OWLSim data for disease '+str(disease_id)+'.')


                    # TODO: Need to add phenolog data lookup here.

                    phenolog_score = ''







                        filename = 'out/phenolog_gene_cand/human_disease_gene_candidate_predictions/all_genes/max_scores/'+str(disease_id_underscored)+'.txt'
                        #print(filename)

                        with open(filename, 'r', encoding="iso-8859-1") as csvfile2:
                            filereader = csv.reader(csvfile2, delimiter='\t', quotechar='\"')
                            for row in filereader:
                                (gene_candidate_ids, gene_candidate_labels, score) = row
                                if len(gene_candidate_ids) == 1:
                                    gene_id = gene_candidate_ids[0]
                                    try:
                                        zebrafish_ncbi_gene_id = human_to_zebrafish_ldo_hash[gene_id]
                                        print('Zebrafish LDO found.')
                                        zfin_gene_id = ncbi_gene_to_zfin_gene_hash[zebrafish_ncbi_gene_id]
                                        print('ZFIN gene ID: '+str(zfin_gene_id))
                                        print(human_disease_gene_prediction_hash[zfin_gene_id])
                                        zebrafish_phenolog_score = ''
                                        zebrafish_max_ic = human_disease_gene_prediction_hash[zfin_gene_id]['maxIC']
                                        zebrafish_iccs = human_disease_gene_prediction_hash[zfin_gene_id]['ICCS']
                                        zebrafish_sim_ic = human_disease_gene_prediction_hash[zfin_gene_id]['simIC']
                                        zebrafish_sim_j = human_disease_gene_prediction_hash[zfin_gene_id]['simJ']



                                if len(gene_candidate_ids) > 1:
                                    for gene_id in gene_candidate_ids:
                                        try:

                                        except:
                                elif len(gene_candidate_ids) == 0:
                                    continue



                        print('Phenolog file open.')
                    except:
                        print('Trouble with Phenolog data for disease '+str(disease_id)+'.')


                    output_row = (disease_gene_association_id, disease_id, gene_id,
                                  zebrafish_max_ic, zebrafish_iccs, zebrafish_sim_ic, zebrafish_sim_j,
                                  mouse_max_ic, mouse_iccs, mouse_sim_ic, mouse_sim_j)
                    csvwriter.writerow(output_row)



        with open('out/roc/max_ic_list.txt', 'wb') as handle:
            pickle.dump(max_ic_list, handle)
        with open('out/roc/iccs_list.txt', 'wb') as handle:
            pickle.dump(iccs_list, handle)
        with open('out/roc/sim_ic_list.txt', 'wb') as handle:
            pickle.dump(sim_ic_list, handle)
        with open('out/roc/sim_j_list.txt', 'wb') as handle:
            pickle.dump(sim_j_list, handle)
        with open('out/roc/phenolog_gene_score_list.txt', 'wb') as handle:
            pickle.dump(phenolog_gene_score_list, handle)
        with open('out/roc/phenolog_ortholog_score_list.txt', 'wb') as handle:
            pickle.dump(phenolog_ortholog_score_list, handle)

        return


#disease_subset = ['ORPHANET_904', 'ORPHANET_84', 'ORPHANET_46348', 'OMIM_272120', 'ORPHANET_2812', 'ORPHANET_791', 'ORPHANET_478', 'ORPHANET_110', 'OMIM_614592', 'ORPHANET_1873', 'OMIM_305400']


counter = multiprocessing.Value(c_int)
counter_lock = multiprocessing.Lock()

def increment():
    """
    This function provides a shared counter for multiprocessing functions for tracking progress.
    :return:
    """
    with counter_lock:
        counter.value += 1
        print('INFO: Processing comparison '+str(counter.value))

def keyfunc(row):
    """
    This function is used by itertools in the chunking of OWLSim queries.
    :param row:
    :return:
    """
    return row[0]



####### OWLSIM QUERY MULTIPROCESSING #######

def multiprocess_owlsim_queries(row):
    """
    This function is used by the process_owlsim_queries function.
    It takes the passed in row and divides into individual variables,
    feeds the URL request to the local instance of the OWLSim server
    :param row:
    :return:
    """

    increment()
    (comparison_id, query_url, entity_a, entity_a_attributes, entity_b, entity_b_attributes) = row[0]
    try:
        response = urllib.request.urlopen(query_url, timeout=1000)
        reader = codecs.getreader("utf-8")
        data = json.load(reader(response))
        results = data['results']
        maxIC = data['results'][0]['maxIC']
        simJ = data['results'][0]['simJ']
        ICCS = data['results'][0]['bmaSymIC']
        simIC = data['results'][0]['simGIC']
        query_flag = 'success'
        sequence = (comparison_id, entity_a, entity_a_attributes, entity_b, entity_b_attributes, maxIC, simJ, ICCS, simIC, query_flag)

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

    return (sequence)

class multithread_owlsim_queries(Thread):
    """
    This class is for the multithreading of OWLSim queries, currently unused.
    """

    def __init__(self, url, name):
        Thread.__init__(self)
        self.name = name
        self.url = url
        return

    def run(self):
        """
        This is a currently unused function used in the multithreading approach for OWLSim queries.
        """
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

####### PHENOLOG FDR MULTIPROCESSING #######

def multiprocess_fdr_calculation(i):
    """

    :param i:
    :return:
    """

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


####### PHENOLOG EXTENSION MULTIPROCESSING #######

def multiprocess_generate_random_human_ext_data(x):
    """
    This function creates random data sets for the determination of the FDR cutoff for
    the significant genotype/disease comparisons of the phenolog extension.
    It takes as input the real species-specific pheontype hash and the disease-phenotype hash for human.
    It creates a random disease-phenotype hash using the real disease-phenotype hash as a guide.
    This allows for the creation of a random data set of similar size and shape
    (same number of diseases with the same number of associated phenotypes for each disease).
    :param x: Number for the random data set.
    :return:
    """
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
        phenotype_draw = phenotypes
        random.shuffle(phenotype_draw)
        for j in geno_pheno_hash[i]:
            random.shuffle(phenotype_draw)
            random_geno_pheno_hash[i].append(phenotype_draw[0])
    with open(('inter/random/human/random_ext_'+str(x)+'.txt'), 'wb') as handle:
        pickle.dump(random_geno_pheno_hash, handle)
    print('Completed human random data set '+str(x)+' out of 1000.')
    return

def multiprocess_generate_random_mouse_ext_data(x):
    """
    This function creates random data sets for the determination of the FDR cutoff for
    the significant genotype/disease comparisons of the phenolog extension.
    It takes as input the real species-specific pheontype hash and the genotype-phenotype hash for mouse.
    It creates a random genotype-phenotype hash using the real genotype-phenotype hash as a guide.
    This allows for the creation of a random data set of similar size and shape
    (same number of genotypes with the same number of associated phenotypes for each genotype).
    :return:
    """
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
        phenotype_draw = phenotypes
        random.shuffle(phenotype_draw)
        for j in geno_pheno_hash[i]:
            random.shuffle(phenotype_draw)
            random_geno_pheno_hash[i].append(phenotype_draw[0])
    with open(('inter/random/mouse/random_ext_'+str(x)+'.txt'), 'wb') as handle:
        pickle.dump(random_geno_pheno_hash, handle)
    print('Completed mouse random data set '+str(x)+' out of 1000.')
    return

def multiprocess_generate_random_zebrafish_ext_data(x):
    """
    This function creates random data sets for the determination of the FDR cutoff for
    the significant genotype/disease comparisons of the phenolog extension.
    It takes as input the real species-specific pheontype hash and the genotype-phenotype hash for zebrafish.
    It creates a random genotype-phenotype hash using the real genotype-phenotype hash as a guide.
    This allows for the creation of a random data set of similar size and shape
    (same number of genotypes with the same number of associated phenotypes for each genotype).
    :param x: Number for the random data set.
    :return:
    """
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
        phenotype_draw = phenotypes
        random.shuffle(phenotype_draw)
        for j in geno_pheno_hash[i]:
            random.shuffle(phenotype_draw)
            random_geno_pheno_hash[i].append(phenotype_draw[0])
    with open(('inter/random/zebrafish/random_ext_'+str(x)+'.txt'), 'wb') as handle:
        pickle.dump(random_geno_pheno_hash, handle)
    print('Completed zebrafish random data set '+str(x)+' out of 1000.')
    return

def multiprocess_ext_fdr_calculation(i):
    """

    :param i:
    :return:
    """

    processing_start_time = time.time()
    #print('INFO: Setting stage for second FDR estimation.')
    # Need to calculate phenolog extension for each pairwise species and combine in order to get a full
    # set of 'genologs' (?) for proper estimation of FDR.

    human_dir = 'inter/random/human/'
    mouse_dir = 'inter/random/mouse/'
    zebrafish_dir = 'inter/random/zebrafish/'

    fdr_p_value_list = []
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

def multiprocess_ext_fdr_calculation_hvz(comparison_list):
    """
    This function is the multi-processed portion of the phenolog extension FDR calculation.
    If there is at least one phenolog match, will perform the hypergeometric probability calculation and return the probability.
    :param comparison_list: One pair-wise combination of human and zebrafish genotypes from the comparison list.
    :return:
    """
    # Activate to track comparisons.
    #increment()

    total_phenotype_matches = 0
    total_phenotype_nonmatches = 0

    species_a_genotype_id = comparison_list[0]
    species_a_phenotypes = read_only_human_geno_pheno_hash[comparison_list[0]]
    genotype_a_phenotype_count = len(species_a_phenotypes)

    # Genotype for species B
    species_b_genotype_id = comparison_list[1]
    species_b_phenotypes = read_only_zebrafish_geno_pheno_hash[comparison_list[1]]
    phenotype_matches = 0
    phenotype_non_matches = 0


    genotype_b_phenotype_count = len(species_b_phenotypes)

    for k in species_a_phenotypes:
        # Orthologs for species A
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
        m = float(genotype_b_phenotype_count)
        n = float(genotype_a_phenotype_count)
        N = float(len(read_only_hvz_phenologs))
        c = float(phenotype_matches)
        prb = float(hypergeom.pmf(c, N, m, n))

        return prb
    else:
        return

def multiprocess_ext_fdr_calculation_hvm(comparison_list):
    """
    This function is the multi-processed portion of the phenolog extension FDR calculation.
    If there is at least one phenolog match, will perform the hypergeometric probability calculation and return the probability.
    :param comparison_list: One pair-wise combination of human and mouse genotypes from the comparison list.
    :return:
    """
    # Activate to track comparisons.
    #increment()

    total_phenotype_matches = 0
    total_phenotype_nonmatches = 0

    species_a_genotype_id = comparison_list[0]
    species_a_phenotypes = read_only_human_geno_pheno_hash[comparison_list[0]]
    #print(species_a_phenotypes)
    genotype_a_phenotype_count = len(species_a_phenotypes)

    # Genotype for species B
    species_b_genotype_id = comparison_list[1]
    species_b_phenotypes = read_only_mouse_geno_pheno_hash[comparison_list[1]]
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
            if ab_combo in read_only_hvm_phenologs or ba_combo in read_only_hvm_phenologs:
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
        N = float(len(read_only_hvm_phenologs))
        c = float(phenotype_matches)
        prb = float(hypergeom.pmf(c, N, m, n))
        #print(str(c)+', '+str(N)+', '+str(m)+', '+str(n))
        #print(prb)
        #phenolog_ext_p_value_list.append(prb)
        #total_hyp_calcs += 1

        return prb
    else:
        return

def multiprocess_ext_fdr_calculation_mvz(comparison_list):
    """
    This function is the multi-processed portion of the phenolog extension FDR calculation.
    If there is at least one phenolog match, will perform the hypergeometric probability calculation and return the probability.
    :param comparison_list: One pair-wise combination of mouse and zebrafish genotypesfrom the comparison list.
    :return:
    """
    # Activate to track comparisons.
    #increment()

    total_phenotype_matches = 0
    total_phenotype_nonmatches = 0

    species_a_genotype_id = comparison_list[0]
    species_a_phenotypes = read_only_mouse_geno_pheno_hash[comparison_list[0]]
    genotype_a_phenotype_count = len(species_a_phenotypes)

    # Genotype for species B
    species_b_genotype_id = comparison_list[1]
    species_b_phenotypes = read_only_zebrafish_geno_pheno_hash[comparison_list[1]]
    phenotype_matches = 0
    phenotype_non_matches = 0
    genotype_b_phenotype_count = len(species_b_phenotypes)

    for k in species_a_phenotypes:
        # Orthologs for species A
        species_a_phenotype = k
        for l in species_b_phenotypes:
            # Orthologs for species B
            species_b_phenotype = l

            ab_combo = species_a_phenotype+'_'+species_b_phenotype
            ba_combo = species_b_phenotype+'_'+species_a_phenotype
            if ab_combo in read_only_mvz_phenologs or ba_combo in read_only_mvz_phenologs:
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
        N = float(len(read_only_mvz_phenologs))
        c = float(phenotype_matches)
        prb = float(hypergeom.pmf(c, N, m, n))
        #print(str(c)+', '+str(N)+', '+str(m)+', '+str(n))
        #print(prb)
        #phenolog_ext_p_value_list.append(prb)
        #total_hyp_calcs += 1

        return prb
    else:
        return


####### PHENOLOG GENE CANDIDATE MULTIPROCESSING #######

def multiprocess_matrix_comparisons(matrix_coordinates):
    """
    This function processes
    :param matrix_coordinates:
    :return:
    """

    #Total number of orthologs is 3572. Hard coding to remove multiple openings of the ortholog_list.txt file.
    len_ortholog_list = 3572

    phenotype_index_i = matrix_coordinates[0]
    phenotype_index_j = matrix_coordinates[1]

    ortholog_counter = 0
    ortholog_match = 0

    (coefficient, p_value) = pearsonr(read_only_ortholog_phenotype_matrix[phenotype_index_i], read_only_ortholog_phenotype_matrix[phenotype_index_j])
    for x in range(0, (len_ortholog_list)):
        if read_only_ortholog_phenotype_matrix[phenotype_index_i][x] == 1 and read_only_ortholog_phenotype_matrix[phenotype_index_j][x] == 1:
            ortholog_match += 1
        ortholog_counter += 1

    # N = total number of orthologs shared between species
    # n = nummber of orthologs in species A phenotype
    # m = nummber of orthologs in species B phenotype
    # c = number of common orthologs between phenotypes (ortholog matches)

    m = float(numpy.sum(read_only_ortholog_phenotype_matrix[phenotype_index_i]))
    n = float(numpy.sum(read_only_ortholog_phenotype_matrix[phenotype_index_j]))
    N = float(len_ortholog_list)
    c = float(ortholog_match)
    hyp_prob = (hypergeom.cdf(c, N, m, n))
    return (phenotype_index_i, phenotype_index_j, hyp_prob, coefficient)

def ortholog_matches(ortholog_phenotype_matrix, phenotype_index_i, phenotype_index_j, x):
    """

    :param ortholog_phenotype_matrix:
    :param phenotype_index_i:
    :param phenotype_index_j:
    :param x:
    :return:
    """
    ortholog_counter = 0
    ortholog_match = 0
    if ortholog_phenotype_matrix[phenotype_index_i][x] == 1 and ortholog_phenotype_matrix[phenotype_index_j][x] == 1:
        ortholog_match += 1
    ortholog_counter += 1
    return ortholog_match, ortholog_counter



####### MAIN #######

# Here is where I'm calling each of the functions individually, as needed.

limit = None

main = main()

####### ASSEMBLY OF CROSS-SPECIES ORTHOLOG LISTS #######

# Trim the PANTHER data set for each taxon.
#main.trim_panther_data('inter/panther/panther_human.txt', ['NCBITaxon:9606'])
#main.trim_panther_data('inter/panther/panther_mouse.txt', ['NCBITaxon:10090'])
#main.trim_panther_data('inter/panther/panther_zebrafish.txt', ['NCBITaxon:7955'])
#main.trim_panther_data('inter/panther/panther_hmz_trio.txt', ['NCBITaxon:9606', 'NCBITaxon:7955', 'NCBITaxon:10090'])

# Assemble ortholog lookup-tables for each pair of species.
#main.get_common_orthologs('inter/panther/common_orthologs_human_zebrafish.txt', ['NCBITaxon:9606', 'NCBITaxon:7955'])
#main.get_common_orthologs('inter/panther/common_orthologs_human_mouse.txt', ['NCBITaxon:9606', 'NCBITaxon:10090'])
#main.get_common_orthologs('inter/panther/common_orthologs_mouse_zebrafish.txt', ['NCBITaxon:10090', 'NCBITaxon:7955'])



####### DATA ASSEMBLY VIA SCIGRAPH #######
# I have abandoned using SciGraph in favor of the flat tables from NIF/DISCO, as the URL queries were too slow.
#main._assemble_human_disease_to_phenotype(limit)
#main._assemble_mouse_genotype_to_phenotype(limit)
#main.assemble_zebrafish_genotype_to_phenotype(limit)


####### DATA ASSEMBLY VIA  NIF/DISCO #######

# Assemble the phenotype ID to label files.
#main.assemble_nif_hpo_phenotype_id_to_label()
#main.assemble_nif_mgi_phenotype_id_to_label()
#main.assemble_nif_zfin_phenotype_id_to_label()

# Assemble the gene ID to label files.
#main.assemble_nif_mgi_gene_id_to_label()
#main.assemble_nif_zfin_gene_id_to_label()
#main.assemble_nif_hpo_gene_id_to_label()

# Assemble the phenotype to gene files for phenologs.
#main.assemble_nif_zfin_phenotype_to_gene(limit)  # Completed in 3.22 days, 85118 rows processed. With hash speed improvement, 8 hours.
#main.assemble_nif_mgi_phenotype_to_gene(limit)  # # Completed on full data set in 175.3 hours (7.3 days) With hash speed improvement, 11.5 hours.
#main.assemble_nif_hpo_phenotype_to_gene(limit)  # Completed on full data set in 75.5 hours. With hash speed improvement, 128 minutes.
#main.assemble_nif_animalqtl_phenotype_to_gene(limit)

# Assemble the files for OWLSim queries and phenolog extension.
#main.assemble_nif_zfin_genotype_to_phenotype(limit)
#main.assemble_nif_mgi_genotype_to_phenotype(limit)e
#main.assemble_nif_mgi_gene_to_phenotype(limit)
#main.assemble_nif_zfin_gene_to_phenotype(limit)
#main.assemble_nif_hpo_disease_to_phenotype(limit)
#main.assemble_nif_hpo_disease_to_gene(limit)



####### ASSEMBLE OWLSIM COMPARISON QUERIES #######
# These scripts assemble the OWLSim queries and save them to files, 5 million queries/file.
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
#latest test: 320/second, or 100K/5 minutes

#Processing...
#Human Diseases = 9214
#Mouse Genotypes = 56427
#Total comparisons = 519,918,378
# Compare human disease phenotypic profiles & mouse genotype phenotypic profiles via OWLSim.
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/owlsim/human_disease_mouse_genotype', 'human_disease_mouse_genotype_queries', 'out/owlsim/human_disease_mouse_genotype', 'human_disease_mouse_genotype_results', 104)

#Process completed!
#Human Diseases = 9214
#zebrafish Genotype = 8535
#Total comparisons = 78,641,490
# Compare human disease phenotypic profiles & zebrafish genotype phenotypic profiles via OWLSim.
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'inter/owlsim/human_disease_zebrafish_genotype', 'human_disease_zebrafish_genotype_queries', 'out/owlsim/human_disease_zebrafish_genotype', 'human_disease_zebrafish_genotype_results', 16)

#Processing completed!
#Mouse genotype = 56427
#zebrafish genotype = 8535
#Total comparisons = 481,604,445
# Compare mouse genotype phenotypic profiles & zebrafish genotype phenotypic profiles via OWLSim.
#main.perform_owlsim_queries('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'inter/owlsim/mouse_genotype_zebrafish_genotype', 'mouse_genotype_zebrafish_genotype_queries', 'out/owlsim/mouse_genotype_zebrafish_genotype', 'mouse_genotype_zebrafish_genotype_results', 97)

#Processing completed!
#Human Diseases = 9214
#Mouse genes = 13102
#Total comparisons = 120,712,614
# Compare human disease phenotypic profiles & mouse gene phenotypic profiles via OWLSim.
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/owlsim/human_disease_mouse_gene','human_disease_mouse_gene_queries', 'out/owlsim/human_disease_mouse_gene', 'human_disease_mouse_gene_results', 25)

#Processing completed!
#Human Diseases = 9214
#zebrafish Genes = 4580
#Total comparisons = 42,190,906
# Compare human disease phenotypic profiles & zebrafish gene phenotypic profiles via OWLSim.
#(self, raw1, raw2, interfile_directory, interfile_prefix, outfile_directory, , outfile_prefix, num_files, limit=None)
#main.perform_owlsim_queries('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'inter/owlsim/human_disease_zebrafish_gene', 'human_disease_zebrafish_gene_queries', 'out/owlsim/human_disease_zebrafish_gene', 'human_disease_zebrafish_gene_results', 9)

#Processing completed!
#Mouse Genes = 13102
#zebrafish Genes = 4580
#Total comparisons = 59,989,479
# Compare mouse gene phenotypic profiles & zebrafish gene phenotypic profiles via OWLSim.
#main.perform_owlsim_queries('inter/mgi/mouse_gene_phenotype_hash.txt', 'inter/zfin/zebrafish_gene_to_phenotype_hash.txt', 'inter/owlsim/mouse_gene_zebrafish_gene','mouse_gene_zebrafish_gene_queries', 'out/owlsim/mouse_gene_zebrafish_gene', 'mouse_gene_zebrafish_gene_results', 12)

#main.trim_owlsim_output('out/owlsim/human_disease_mouse_gene/human_disease_mouse_gene_results_', 'out/owlsim/human_disease_mouse_gene/human_disease_mouse_gene_trimmed_results_')
#main.assemble_owlsim_gene_candidates('out/owlsim/human_disease_mouse_gene/human_disease_mouse_gene_results_', 'out/owlsim/human_disease_mouse_gene/human_disease_mouse_gene_predictions.txt')

#main.assemble_owlsim_gene_candidate_alternate()
#main.assemble_owlsim_top_20_gene_candidates()


####### PHENOLOG FDR CALCULATION #######

# Generate random data sets for each organism using common orthologs between the other organisms.
#main.generate_random_data('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_mouse.txt', 'inter/random/human_vs_mouse/mouse/')
#main.generate_random_data('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_mouse.txt', 'inter/random/human_vs_mouse/human/')

#main.generate_random_data('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', 'inter/random/human_vs_zebrafish/zebrafish/')
#main.generate_random_data('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_human_zebrafish.txt', 'inter/random/human_vs_zebrafish/human/')

#main.generate_random_data('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', 'inter/random/mouse_vs_zebrafish/zebrafish/')
#main.generate_random_data('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/panther/common_orthologs_mouse_zebrafish.txt', 'inter/random/mouse_vs_zebrafish/mouse/')

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
#main.annotate_significant_phenologs_with_scores_and_labels()


####### PHENOLOG EXTENSION FDR CALCULATION #######

# Parse each of the phenotype ontology files to obtain a set of all possible phenotypes for each species.
#main.parse_hp('raw/ontologies/hp.obo', 'inter/ontologies/hp_hash.txt')
#main.parse_mp('raw/ontologies/MPheno_OBO.ontology', 'inter/ontologies/mp_hash.txt')
#main.parse_zp('raw/ontologies/zp_mapping.txt', 'inter/ontologies/zp_hash.txt')

# Create random data sets for phenlog extension FDR from all possible phenotypes and
# the significant phenologs between species pair.
#main.generate_human_random_ext_data()
#main.generate_mouse_random_ext_data()
#main.generate_zebrafish_random_ext_data()

#Assemble the phenolog lookup files.
#main.assemble_hvm_phenologs()
#main.assemble_hvz_phenologs()
#main.assemble_mvz_phenologs()


#main.set_stage_for_extension_fdr_calculation()
#gc.set_debug(gc.DEBUG_LEAK)


# The next three code snippets process the phenolog extension calculations to determine the FDR.
# However, I have switched to using the three code snippets after this section in combination with a bash script,
# as there is an issue with Python not deleting lists of floats even when specifically removed.
'''
with open('inter/phenolog/hvz_phenolog_combo.txt', 'rb') as handle:
    read_only_hvz_phenologs = set(pickle.load(handle))
for i in range(1, 1001):
    with open('inter/random/human/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_human_geno_pheno_hash = pickle.load(handle)
    with open('inter/random/zebrafish/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_zebrafish_geno_pheno_hash = pickle.load(handle)
    print('INFO: Processing human vs zebrafish random data set '+str(i)+'.')
    p_value_out_file = 'inter/phenolog_ext/hvz_p_values/hvz_p_values_'+str(i)+'.txt'
    main.perform_phenolog_calculations_for_ext_fdr_hvz(read_only_human_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash, p_value_out_file)
    del read_only_human_geno_pheno_hash
    del read_only_zebrafish_geno_pheno_hash
    gc.collect()
    print('INFO: Done processing human vs zebrafish random data set '+str(i)+'.')


with open('inter/phenolog/hvm_phenolog_combo.txt', 'rb') as handle:
    read_only_hvm_phenologs = set(pickle.load(handle))
for i in range(1, 1001):
    with open('inter/random/human/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_human_geno_pheno_hash = pickle.load(handle)
    with open('inter/random/mouse/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_mouse_geno_pheno_hash = pickle.load(handle)
    print('INFO: Processing human vs mouse random data set '+str(i)+'.')
    p_value_out_file = 'inter/phenolog_ext/hvm_p_values/hvm_p_values_'+str(i)+'.txt'
    main.perform_phenolog_calculations_for_ext_fdr_hvm(read_only_human_geno_pheno_hash, read_only_mouse_geno_pheno_hash, p_value_out_file)
    read_only_human_geno_pheno_hash = None
    read_only_mouse_geno_pheno_hash = None
    gc.collect()
    print('INFO: Done processing human vs mouse random data set '+str(i)+'.')


with open('inter/phenolog/mvz_phenolog_combo.txt', 'rb') as handle:
    read_only_mvz_phenologs = set(pickle.load(handle))
for i in range(1, 1001):
    with open('inter/random/mouse/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_mouse_geno_pheno_hash = pickle.load(handle)
    with open('inter/random/zebrafish/random_ext_'+str(i)+'.txt', 'rb') as handle:
        read_only_zebrafish_geno_pheno_hash = pickle.load(handle)
    print('INFO: Processing mouse vs zebrafish random data set '+str(i)+'.')
    p_value_out_file = 'inter/phenolog_ext/mvz_p_values/mvz_p_values_'+str(i)+'.txt'
    main.perform_phenolog_calculations_for_ext_fdr_mvz(read_only_mouse_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash, p_value_out_file)
    read_only_mouse_geno_pheno_hash = None
    read_only_zebrafish_geno_pheno_hash = None
    gc.collect()
    print('INFO: Done processing mouse vs zebrafish random data set '+str(i)+'.')
'''

# The next three code snippets process the phenolog extension calculations to determine the FDR using an external bash script.
# This gets around the memory management issue.
'''
with open('inter/phenolog/hvz_phenolog_combo.txt', 'rb') as handle:
    read_only_hvz_phenologs = set(pickle.load(handle))
with open('inter/random/human/random_ext_'+str(sys.argv[1])+'.txt', 'rb') as handle:
    read_only_human_geno_pheno_hash = pickle.load(handle)
with open('inter/random/zebrafish/random_ext_'+str(sys.argv[1])+'.txt', 'rb') as handle:
    read_only_zebrafish_geno_pheno_hash = pickle.load(handle)
print('INFO: Processing human vs zebrafish random data set '+str(sys.argv[1])+'.')
p_value_out_file = 'inter/phenolog_ext/hvz_p_values/hvz_p_values_'+str(sys.argv[1])+'.txt'
main.perform_phenolog_calculations_for_ext_fdr_hvz(read_only_human_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash, p_value_out_file)
del read_only_human_geno_pheno_hash
del read_only_zebrafish_geno_pheno_hash
gc.collect()
print('INFO: Done processing human vs zebrafish random data set '+str(sys.argv[1])+'.')
'''
'''
with open('inter/phenolog/hvm_phenolog_combo.txt', 'rb') as handle:
    read_only_hvm_phenologs = set(pickle.load(handle))
with open('inter/random/human/random_ext_'+str(sys.argv[1])+'.txt', 'rb') as handle:
    read_only_human_geno_pheno_hash = pickle.load(handle)
with open('inter/random/mouse/random_ext_'+str(sys.argv[1])+'.txt', 'rb') as handle:
    read_only_mouse_geno_pheno_hash = pickle.load(handle)
print('INFO: Processing human vs mouse random data set '+str(sys.argv[1])+'.')
p_value_out_file = 'inter/phenolog_ext/hvm_p_values/hvm_p_values_'+str(sys.argv[1])+'.txt'
main.perform_phenolog_calculations_for_ext_fdr_hvm(read_only_human_geno_pheno_hash, read_only_mouse_geno_pheno_hash, p_value_out_file)
read_only_human_geno_pheno_hash = None
read_only_mouse_geno_pheno_hash = None
gc.collect()
print('INFO: Done processing human vs mouse random data set '+str(sys.argv[1])+'.')
'''
'''
with open('inter/phenolog/mvz_phenolog_combo.txt', 'rb') as handle:
    read_only_mvz_phenologs = set(pickle.load(handle))
with open('inter/random/mouse/random_ext_'+str(sys.argv[1])+'.txt', 'rb') as handle:
    read_only_mouse_geno_pheno_hash = pickle.load(handle)
with open('inter/random/zebrafish/random_ext_'+str(sys.argv[1])+'.txt', 'rb') as handle:
    read_only_zebrafish_geno_pheno_hash = pickle.load(handle)
print('INFO: Processing mouse vs zebrafish random data set '+str(sys.argv[1])+'.')
p_value_out_file = 'inter/phenolog_ext/mvz_p_values/mvz_p_values_'+str(sys.argv[1])+'.txt'
main.perform_phenolog_calculations_for_ext_fdr_mvz(read_only_mouse_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash, p_value_out_file)
read_only_mouse_geno_pheno_hash = None
read_only_zebrafish_geno_pheno_hash = None
gc.collect()
print('INFO: Done processing mouse vs zebrafish random data set '+str(sys.argv[1])+'.')
'''

#main.identify_significance_threshold_for_random_data_sets()
#main.count_zeroes()

#main.perform_phenolog_calculations_for_ext_fdr(read_only_human_geno_pheno_hash, read_only_mouse_geno_pheno_hash)
#main.perform_hvm_phenolog_calculations_for_ext_fdr_alternate(read_only_human_geno_pheno_hash, read_only_mouse_geno_pheno_hash)
#main.perform_mvz_phenolog_calculations_for_ext_fdr_alternate(read_only_mouse_geno_pheno_hash, read_only_zebrafish_geno_pheno_hash)

#with open('inter/phenolog/hvm_phenolog_combo.txt', 'rb') as handle:
#    read_only_hvz_phenologs = set(pickle.load(handle))
#print(len(read_only_hvz_phenologs))
# Made up FDR for testing purposes.
#ext_fdr_cutoff = 0.00022089684117479534

# Once the FDR has been determined, can run the actual phenolog extension calculations using the FDR cutoff.
#main.perform_phenolog_ext_calculations('inter/hpo/human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt', 'out/phenolog_ext/human_vs_mouse.txt', 'inter/phenolog/hvm_significant_phenologs.txt', ext_fdr_cutoff)
#main.perform_phenolog_ext_calculations('inter/hpo/human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'out/phenolog_ext/human_vs_zebrafish.txt', 'inter/phenolog/hvz_significant_phenologs.txt', ext_fdr_cutoff)
#main.perform_phenolog_ext_calculations('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt', 'out/phenolog_ext/mouse_vs_zebrafish.txt', 'inter/phenolog/mvz_significant_phenologs.txt', ext_fdr_cutoff)


####### PHENOLOG GENE CANDIDATE PREDICTIONS #######

#NOTE: Must be assembled after the nif phenotype to gene assembly has been performed.
#NOTE: These three functions are now combined and run in the create_phenolog_gene_candidate_matrices function.
#(human_ortholog_phenotype_matrix, human_phenotype_list, human_ortholog_list) = main.assemble_ortholog_phenotype_matrix('inter/hpo/human_pheno_ortholog_hash.txt', 'inter/hpo/human_pheno_ortholog_matrix.txt')
#(mouse_ortholog_phenotype_matrix, mouse_phenotype_list, mouse_ortholog_list) = main.assemble_ortholog_phenotype_matrix('inter/mgi/mouse_pheno_ortholog_hash.txt', 'inter/mgi/mouse_pheno_ortholog_matrix.txt')
#(zebrafish_ortholog_phenotype_matrix, zebrafish_phenotype_list, zebrafish_ortholog_list) = main.assemble_ortholog_phenotype_matrix('inter/zfin/zebrafish_pheno_ortholog_hash.txt', 'inter/zfin/zebrafish_pheno_ortholog_matrix.txt')

#This process requires multi-processing due to the large number of comparisons that need to be performed.

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
#main.assemble_phenolog_gene_candidate_predictions_for_phenotypes()




#main.assemble_model_level_phenolog_gene_candidate_predictions()

#DON'T USE!
#main.assemble_nearest_neighbor_phenotypes_hash()

#main.assemble_phenolog_gene_candidates_for_diseases()
#main.assemble_phenolog_orthogroup_candidates_for_diseases()

#with open('inter/hpo/human_pheno_gene_hash.txt', 'rb') as handle:
    #human_pheno_ortholog_hash = pickle.load(handle)
#print(human_pheno_ortholog_hash['HP:0005532'])
#main.assemble_gene_candidates_for_diseases_max_score()


#read_only_disease_subset = ['OMIM_157900', 'OMIM_167400', 'OMIM_260530', 'ORPHANET_904', 'ORPHANET_84', 'ORPHANET_46348', 'OMIM_272120', 'ORPHANET_2812', 'ORPHANET_791', 'ORPHANET_478', 'ORPHANET_110', 'OMIM_614592', 'ORPHANET_1873', 'OMIM_305400', 'OMIM_157900']

with open('inter/omim/disorder_list.txt', 'rb') as handle:
    read_only_disease_subset = pickle.load(handle)

#print(str(len(read_only_disease_subset)))
#main.assemble_owlsim_top_20_gene_candidates()
#main.assemble_phenolog_gene_candidates_for_diseases()
#main.assemble_phenolog_orthogroup_candidates_for_diseases()

#main.assemble_complete_disease_list()

####### OMIM ASSERTED MODELS #######


#main.create_zfin_gene_to_ncbi_hash()
#main.create_mgi_gene_to_ncbi_hash()
#main.create_mim_to_gene_hash()

#main.convert_morbid_map_to_ncbi_gene()
#main.trim_morbid_disease_to_gene()
#main.convert_morbid_map_genes_to_orthologs()
#main.create_ortholog_lookup_hashes()
#main.assemble_owlsim_data_for_ROC()
#main.assemble_ROC_score_lists()
main.assemble_ROC_score_lists_alternate()

elapsed_time = time.time() - start_time
print('Processing completed in '+str(elapsed_time)+' seconds.')

#TODO: Make sure and have the ability to filter between single-gene genotypes and multi-gene genotypes.
#
# URL Format for OWLSim queries:
#http://owlsim.monarchinitiative.org/compareAttributeSets?a=HP:0001263&b=MP:0010864
# URL format for mutliple phenotypes:
#http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878


#http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0002169&b=ZP:0000411&b=ZP:0002713&b=ZP:&b=ZP:&b=ZP:0000692&b=ZP:0000054&b=ZP:0000043&b=ZP:0000038&b=ZP:0000737&b=ZP:&b=ZP:0001192&b=ZP:&b=3300
# Local server URL format: http://0.0.0.0:9031/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878