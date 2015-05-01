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
        'dvp.pr_nlx_84521_1'  # PANTHER:Orthologs
    ]

    files = {
        'aqtl' : {'file' : 'dvp.pr_nif_0000_02550_3'},
        'mgi' : {'file' : 'dvp.pr_nif_0000_00096_6'}
    }


    def _assemble_human_disease_to_phenotype(self, limit=None):
        print('INFO: Assembling human disease to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/hpo/diseases.csv'
        out = 'out/hpo/human_disease_pheno_hash.txt'
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

                #print(len(hu_disease_to_phenotype_hash.keys()))
        #print(hu_disease_to_phenotype_hash)
                #if disease_id not in hu_disease_to_phenotype_hash:
                    #hu_disease_to_phenotype_hash['disease_id'] =
        with open(out, 'wb') as handle:
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
        out = 'out/mgi/mouse_geno_pheno_hash.txt'
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


        with open(out, 'wb') as handle:
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
        out = 'out/zfin/zebrafish_geno_pheno_hash.txt'
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
        with open(out, 'wb') as handle:
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
        out = 'out/zfin/zebrafish_geno_pheno_hash.txt'
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
        with open(out, 'wb') as handle:
            pickle.dump(zfin_genotype_to_phenotype_hash, handle)
        print('INFO: Done assembling zebrafish genotype to phenotype data.')
        print('INFO: '+str(len(zfin_genotype_to_phenotype_hash.keys()))+' zebrafish genotypes present.')
        print('INFO: '+str(failure_counter)+' failed to retrieve through SciGraph services.')
        return



    def assemble_nif_zfin_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling zebrafish genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/zfin/dvp.pr_nif_0000_21427_10'
        out = 'out/zfin/zebrafish_pheno_gene_hash.txt'
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
        with open(out, 'wb') as handle:
            pickle.dump(zfin_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling zebrafish phenotype to gene data.')
        print('INFO: '+str(len(zfin_phenotype_to_gene_hash.keys()))+' zebrafish phenotypes present.')
        return

    def assemble_nif_mgi_phenotype_to_gene(self, limit=None):
        print('INFO:Assembling mouse genotype to phenotype data.')
        line_counter = 0
        failure_counter = 0
        raw = 'raw/mgi/dvp.pr_nif_0000_00096_6'
        out = 'out/mgi/mouse_pheno_gene_hash.txt'
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
        with open(out, 'wb') as handle:
            pickle.dump(mgi_phenotype_to_gene_hash, handle)
        print('INFO: Done assembling mouse phenotype to gene data.')
        print('INFO: '+str(len(mgi_phenotype_to_gene_hash.keys()))+' mouse phenotypes present.')
        return



###MAIN####
limit = 100
#limit = 100
main = main()
#main._assemble_human_disease_to_phenotype(limit)
#main._assemble_mouse_genotype_to_phenotype(limit)
main.assemble_nif_zfin_phenotype_to_gene(limit)
main.assemble_nif_mgi_phenotype_to_gene(limit)
#FIXME: Note that the zebrafish data is not currently available through REST services.
#main.assemble_zebrafish_genotype_to_phenotype(500)


elapsed_time = time.time() - start_time
print('Processing completed in '+str(elapsed_time)+' seconds.')

#FIXME:OWLSim server call. Server is currently down.
#http://owlsim.monarchinitiative.org/compareAttributeSets?a=HP:0001263&b=MP:0010864

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

# If there are matching orthologs
# Call hypergeometric probability calculation.




