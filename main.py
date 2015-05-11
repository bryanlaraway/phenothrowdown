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

start_time = time.time()

hu_disease_to_phenotype_hash = {'disease_id': {}}
mouse_genotype_to_phenotype_hash = {'genotype_id': {}}
zfin_genotype_to_phenotype_hash = {'genotype_id': {}}


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

    # Need list of phenotypes with all associated genes for human, mouse, zebrafish

    #FIXME: Creating a general phenotype-gene processing script.
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

    def assemble_nif_zfin_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling zebrafish genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        inter = 'inter/zfin/zebrafish_pheno_gene_hash.txt'
        zfin_phenotype_to_gene_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish phenotype rows to process.')
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

                print(phenotype_id)
                #FIXME: Going to need to convert the ZFIN Gene IDs to NCBIGene IDs.
                genes = implicated_gene_ids.split()
                print(genes)
                if phenotype_id not in zfin_phenotype_to_gene_hash:
                    zfin_phenotype_to_gene_hash[phenotype_id] = genes
                    #print(zfin_genotype_to_phenotype_hash[genotype_id])
                else:
                    zfin_phenotype_to_gene_hash[phenotype_id].append(genes)
                    #print(zfin_genotype_to_phenotype_hash[genotype_id])
                    #print(len(zfin_genotype_to_phenotype_hash.keys()))
                    print('Repeat phenotype: '+phenotype_id)
                if limit is not None and line_counter > limit:
                    break
        with open(inter, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling zebrafish phenotype to gene data.')
        print('INFO: '+str(len(zfin_phenotype_to_gene_hash.keys()))+' zebrafish phenotypes present.')
        return

    def assemble_nif_mgi_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        inter = 'inter/mgi/mouse_pheno_gene_hash.txt'
        mgi_phenotype_to_gene_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' zebrafish phenotype rows to process.')
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
                 environment_label, publication_id, publication_label, publication_url, taxon_id, taxon_label,
                 e_uid, v_uid, v_uuid, v_lastmodified) = row

                print(phenotype_id)
                #FIXME: Going to need to convert the MGI Gene IDs to NCBIGene IDs.
                genes = implicated_gene_ids.split()
                print(genes)
                if phenotype_id not in mgi_phenotype_to_gene_hash:
                    mgi_phenotype_to_gene_hash[phenotype_id] = genes
                    #print(mgi_phenotype_to_gene_hash[genotype_id])
                else:
                    mgi_phenotype_to_gene_hash[phenotype_id].append(genes)
                    #print(mgi_phenotype_to_gene_hash[genotype_id])
                    #print(len(mgi_phenotype_to_gene_hash.keys()))
                    print('Repeat phenotype: '+phenotype_id)
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling mouse phenotype to gene data.')
        print('INFO: '+str(len(mgi_phenotype_to_gene_hash.keys()))+' mouse phenotypes present.')
        return


    def assemble_nif_hpo_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/dvp.pr_nlx_151835_2'
        inter = 'inter/hpo/human_pheno_gene_hash.txt'
        hpo_phenotype_to_gene_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            next(filereader,None)
            for row in filereader:
                line_counter += 1
                (e_uid, phenotype_id, phenotype_label, gene_id, gene_num,
                 gene_label, v_uid, v_uuid, v_lastmodified) = row

                print(phenotype_id)
                #FIXME: Going to need to convert the MGI Gene IDs to NCBIGene IDs.

                #print(genes)
                if phenotype_id not in hpo_phenotype_to_gene_hash:
                    hpo_phenotype_to_gene_hash[phenotype_id] = [gene_id]
                    #print(hpo_phenotype_to_gene_hash[genotype_id])
                else:
                    hpo_phenotype_to_gene_hash[phenotype_id].append(gene_id)
                    #print(hpo_phenotype_to_gene_hash[genotype_id])
                    #print(len(hpo_phenotype_to_gene_hash.keys()))
                    print('Repeat phenotype: '+phenotype_id)
                if limit is not None and line_counter > limit:
                    break
        #TODO: Need to filter out phenotypes that don't have any associated genes.
        with open(inter, 'wb') as handle:
            pickle.dump(hpo_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling human phenotype to gene data.')
        print('INFO: '+str(len(hpo_phenotype_to_gene_hash.keys()))+' human phenotypes present.')
        return

    def assemble_nif_animalqtl_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/animalqtldb/dvp.pr_nif_0000_02550_3'
        inter = 'inter/aqtl/aqtl_pheno_gene_hash.txt'
        aqtl_phenotype_to_gene_hash = {}
        with open(raw, 'r', encoding="iso-8859-1") as csvfile:
            filereader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            row_count = sum(1 for row in filereader)
            row_count = row_count - 1
            print(str(row_count)+' human phenotype rows to process.')
        if limit is not None:
            print('Only parsing first '+str(limit)+' rows.' )
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
                if extrinsic_genotype_id == '' or extrinsic_genotype_id is None:
                    print('Skipping genotype with extrinsic modifiers: '+effective_genotype_id)
                    #print(phenotype_id)
                    #FIXME: Going to need to convert the ZFIN Gene IDs to NCBIGene IDs.
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
                #FIXME: Going to need to convert the ZFIN Gene IDs to NCBIGene IDs.
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
                #FIXME: Going to need to convert the MGI Gene IDs to NCBIGene IDs.

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
                #FIXME: Going to need to convert the MGI Gene IDs to NCBIGene IDs.

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
                #FIXME: Going to need to convert the MGI Gene IDs to NCBIGene IDs.

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
        failure_counter = 0
        if limit is not None:
            print('Only querying first '+str(limit)+' phenotypic profile pairs.')
        #raw1 = 'inter/hpo/nif_human_disease_phenotype_hash.txt'
        #raw2 = 'inter/mgi/mouse_genotype_phenotype_hash.txt'
        data1 = open(raw1, 'rb')
        organism_a_hash = pickle.load(data1)
        data1.close()
        data2 = open(raw2, 'rb')
        organism_b_hash = pickle.load(data2)
        data2.close()
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
            csvwriter = csv.writer(csvfile, delimiter=' ')
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
                        output_row = [panther_speciesa, taxon_id_a, speciesa, taxon_label_a, genea, gene_id_a, gene_label_a,
                        proteina, panther_speciesb, taxon_id_b, speciesb, taxon_label_b, geneb, gene_id_b,
                        gene_label_b, proteinb, orthology_class, orthology_class_label, panther_id]
                        #print('found one')
                        output_line_counter += 1
                        csvwriter.writerow(output_row)

        print('PANTHER file trimmed to '+str(output_line_counter)+' rows.')

        return


    def get_gene_orthologs(self):
        # TODO: Would it make sense to do a pre-processing step where we filter out based on the taxon? Might speed up calculations.

        return





    def perform_phenolog_calculations(self, raw1, raw2, out, limit=None):
        print('INFO: Performing phenolog calculations.')
        line_counter = 0
        failure_counter = 0
        if limit is not None:
            print('Only querying first '+str(limit)+' phenotypic profile pairs.')
        #raw1 = 'inter/hpo/nif_human_disease_phenotype_hash.txt'
        #raw2 = 'inter/mgi/mouse_genotype_phenotype_hash.txt'
        data1 = open(raw1, 'rb')
        organism_a_hash = pickle.load(data1)
        data1.close()
        data2 = open(raw2, 'rb')
        organism_b_hash = pickle.load(data2)
        data2.close()
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

        #query_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878'
        #phenotypic_profile_a = 'a='+('&a=').join(entity_a_attributes)
        #phenotypic_profile_b = '&b='+('&b=').join(entity_b_attributes)
        #combined_url = base_url+phenotypic_profile_a+phenotypic_profile_b
        #print(combined_url)
        #query_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0003731&b=HP:0000580'
        #query_url = 'http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0003731&a=MP:0000599&a=MP:0005331&b=HP:0000580&b=HP:0002240&b=HP:0000831'


        return






###MAIN####

limit = 100
main = main()


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

# Compare human diseases & mouse genotypes via OWLSim.
#print('OWLSim processing human vs mouse')
#main.perform_owlsim_queries('inter/hpo/nif_human_disease_phenotype_hash.txt', 'inter/mgi/mouse_genotype_phenotype_hash.txt','out/owlsim_human_disease_mouse_genotype.txt')
#print('Done processing human vs mouse')
# Compare human diseases & zebrafish genotypes via OWLSim.
#print('OWLSim processing human vs zebrafish')
#main.perform_owlsim_queries('inter/hpo/nif_human_disease_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt','out/owlsim_human_disease_zebrafish_genotype.txt')
#print('Done processing human vs zebrafish')
# Compare mouse genotypes & zebrafish genotypes via OWLSim.
#print('OWLSim processing mouse vs zebrafish')
#main.perform_owlsim_queries('inter/mgi/mouse_genotype_phenotype_hash.txt', 'inter/zfin/zebrafish_genotype_phenotype_hash.txt','out/owlsim_mouse_genotype_zebrafish_genotype.txt')
#print('Done processing mouse vs zebrafish')

#print()

#Trim the PANTHER data set for
main.trim_panther_data('inter/panther/panther_human.txt', ['NCBITaxon:9606'])
main.trim_panther_data('inter/panther/panther_mouse.txt', ['NCBITaxon:10090'])
main.trim_panther_data('inter/panther/panther_zebrafish.txt', ['NCBITaxon:7955'])
main.trim_panther_data('inter/panther/panther_hmz_trio.txt', ['NCBITaxon:9606', 'NCBITaxon:7955', 'NCBITaxon:10090'])
#main.perform_phenolog_calculations()




elapsed_time = time.time() - start_time
print('Processing completed in '+str(elapsed_time)+' seconds.')

#TODO: Make sure and have the ability to filter between single-gene genotypes and multi-gene genotypes.

#http://owlsim.monarchinitiative.org/compareAttributeSets?a=HP:0001263&b=MP:0010864
#Format for mutliple:
#http://owlsim.crbs.ucsd.edu/compareAttributeSets?a=MP:0010864&b=HP:0001263&b=HP:0000878

##############    OBTAIN DATA FROM DATA SOURCES    #############
# PURPOSE: obtain data from the various resources that I will use
# in this study. Will need human, mouse, and zebrafish data.
# NOTE: May just compile the data sources off-line, removing requirements for database logins, etc.
# Then again, perhaps could just extract the NIF/DISCO data through the download services?

### GET HUMAN DATA ###
# Get HPO:Disease to Phenotype data from NIF (nlx_151835-1)
# Get HPO:Disease to Gene data from NIF (nlx_151835-2)
# Get HPO:Phenotype to Gene data from NIF (nlx_151835-3)

### GET MOUSE DATA ###
# Get MGI:MousePhenotypes data from NIF (nif-0000-00096-6) for genotype-phenotype associations.
# Get MGI:MouseGenotypes data from NIF (nif-0000-00096-5)
## NOTE:Genotype view likely not needed, as I can extract genotype, phenotype, and gene data from the Phenotypes view.

### GET ZEBRAFISH DATA ###
# Get Genotype to Phenotype data from NIF (nif-0000-21427-10) for genotype-phenotype associations.
# Get OrganismGenotype data from NIF (nif-0000-21427-11). Like with MGI, this likely won't be needed,
# as I can extract the implicated genes from the phenotype view for the gene-phenotype associations.

### GET GENE ORTHOLOG DATA ###
# Get Panther:Orthologs data from NIF (nlx_84521-1)
# NOTE: For optimization, will want to only grab data for human, mouse, and zebrafish,
# so select only rows where Species A/Species B in {Homo sapiens, Mus musculus, Danio rerio}.
# Should greatly reduce the size of the table.

### GET ADDITIONAL DATA ###
# Depending on time available/data available, will also need to get gene-phenotype data from other model organisms for the additional Phenologs work.





##############    SCRUB DATA    #############
# PURPOSE: Perform any pre-scrubbing that is necessary for each resource/data file.
# Data from NIF sources should fortunately be decently scrubbed (confirm this before proceeding),
# but other resources will likely require additional scrubbing before further processing.


#def _scrub_animal_qtl(self, limit=None):
#raw = ('/').join((self.rawdir,self.files['aqtl']['file']))
#out =



##############    PHENOLOG DATA PRE-PROCESSING    #############
#PURPOSE: Need to perform all of the preparatory steps to make the data ready for the phenologs calculations.
#Steps include





# HPO Annotations: Phenotype to Gene view ((nlx_151835-2)
# Columns available: e_uid, phenotype_id, phenotype_label, gene_id, gene_num, gene_label, v_uid, v_uuid, v_lastmodified
# phenotype_id = HP:1234567
# gene_id = NCBI_gene:1234567




# HPO Annotations: Disease Phenotypes view ((nlx_151835-1)
# Columns available: e_uid, disorder_id, disorder_database_prefix,disorder_database_link, disorder_name,
# disorder_qualifier,phenotype_id,phenotype_label, publication_id, evidence_code_id, evidence_code_label,
# onset_id,onset_label,frequency,aspect, aspect_text, synonyms, v_uid, v_uuid, v_lastmodified
# Sources: OMIM and, ORPHANET only. Can filter on just one source if necessary.
# disease_id = OMIM:1234567 or ORPHANET:1234567
# phenotype_id = HP:1234567


# HPO Annotations: Disease to Gene view ((nlx_151835-3)
# Columns available: e_uid, disease_id, disorder_name, disorder_database_link, gene_id, gene_num, gene_label, v_uid, v_uuid, v_lastmodified
# disease_id = OMIM:1234567 or ORPHANET:1234567 or DECIPHER:1234567
# gene_id = NCBI_gene:1234567
# Sources: OMIM, DECIPHER, ORPHANET. Can filter on just one source if necessary.




##############    ASSEMBLE DATA FOR PHENOLOGS    #############
# PURPOSE: To assemble the phenotype-gene lists necessary for performing the calculation of phenologs.
# NOTE: Create a method for this, as it should be reusable for all gene-phenotype tables for each species.
# CAVEATS: If a phenotype has less than 3 associated genes, then the phenotype must be removed from the analysis.
# Data input format: table with rows of gene_id - phenotype_id columns.
# Data output format: table with phenotype_id, [gene_id array] columns.

# Process for assembling lists:
# Create output table
# Open input table of scrubbed phenotype-gene associations.
# Iterate through each row of a scrubbed phenotype-gene resource table.
# If phenotype ID is not in output table
# Add phenotype ID to output table (column 0)
# Add gene ID to new array in output table (column 0, same row)
# ElseIf phenotype ID is in output table
# GoTo row with the phenotype ID
# If gene ID is not in array
# Append gene ID to array

# Close input table
# Close output table


##############    IDENTIFY ORTHOLOGS USING PANTHER    #############



##############    IDENTIFY PHENOLOGS    #############
# PURPOSE: To calculate the phenologs by comparing the phenotype-gene lists between species.


# Data input format: table with phenotype_id, [gene_id array] columns.
# Data output format: table with species_a, species_a_phenotype_id, species_b, species_b_phenotype_id,
# matched ortholog IDs, species A unmatched ortholog IDs, species B unmatched ortholog IDs


# For each pair of species, iterate through each possible phenotype pair between the two input tables.

# Open table from species A
# Set taxon_id for species A
# Open table from species B
# Set taxon_id for species B
# Create output table


# For row in species A table
# For row in species B table
# phenotype_pair_id = phenotype_a_id+'_'+phenotype_b_id
# (Necessary to make a full entry ID? Both phenotype IDs, both taxon IDs, ortholog IDs, etc?)
#TODO: How best to prepare variables for the hypergeometric probability calculation?
#test
# If there are matching orthologs
# Call hypergeometric probability calculation.




