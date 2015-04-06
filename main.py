#AUTHOR: Bryan Laraway
#PROJECT: Comparison of OWLSim and Phenologs for the identification of models of human disease and gene candidates for human disease.
#PURPOSE: This script will call of the functions/methods/scripts for performing the processing required for this analysis.




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




##############    OBTAIN DATA FROM DATA SOURCES    #############
# PURPOSE: obtain data from the various resources that I will use
# in this study. Will need human, mouse, and zebrafish data.

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

# Data format: phenotype_id, [gene_id array]
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










