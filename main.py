#AUTHOR: Bryan Laraway
#PROJECT: Comparison of OWLSim and Phenologs for the identification of models of human disease and gene candidates for human disease.
#PURPOSE: This script will call of the functions/methods/scripts for performing the processing required for this analysis.

import json
import urllib.request
import codecs

import os
import re
import csv

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


### Testing using SciGraph Rest Services

mgi_gene_id = 'MGI:2182454'
mgi_gene_to_phenotype_hash = {}
human_disease_to_phenotype_hash = {}

disease_id = 'OMIM:267450'
disease_url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/diseases/'+disease_id+'/phenotypes/targets'
url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/features/MGI:2182454/phenotypes/targets'
with urllib.request.urlopen(url) as response:

#url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/diseases/OMIM:106240/phenotypes/targets'
#response = urllib.request(url)
    #print(response)
    reader = codecs.getreader("utf-8")
    data = json.load(reader(response))
    #print(data)
    pheno_ids = data['nodes']
    #print(pheno_ids)
    for rs in pheno_ids:
        if mgi_gene_id not in mgi_gene_to_phenotype_hash:
            mgi_gene_to_phenotype_hash[mgi_gene_id] = [rs['id']]
        else:
            mgi_gene_to_phenotype_hash[mgi_gene_id].append(rs['id'])
    print(mgi_gene_to_phenotype_hash)

with urllib.request.urlopen(disease_url) as response:
#url = 'http://rosie.crbs.ucsd.edu:9000/scigraph/dynamic/diseases/OMIM:106240/phenotypes/targets'
#response = urllib.request(url)
    #print(response)
    reader = codecs.getreader("utf-8")
    data = json.load(reader(response))
    #print(data)
    pheno_ids = data['nodes']
    #print(pheno_ids)
    for rs in pheno_ids:
        if disease_id not in human_disease_to_phenotype_hash:
            human_disease_to_phenotype_hash[disease_id] = [rs['id']]
        else:
            human_disease_to_phenotype_hash[disease_id].append(rs['id'])
    print(human_disease_to_phenotype_hash)



###MAIN####

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



