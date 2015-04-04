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
#PURPOSE: obtain data from the various resources that I will use
#in this study. Will need human, mouse, and zebrafish data.

### GET HUMAN DATA ###
#Get HPO:Disease to Phenotype data from NIF (nlx_151835-1)
#Get HPO:Disease to Gene data from NIF (nlx_151835-2)
#Get HPO:Phenotype to Gene data from NIF (nlx_151835-3)

### GET MOUSE DATA ###
#Get MGI:MousePhenotypes data from NIF (nif-0000-00096-6) for genotype-phenotype associations.
#Get MGI:MouseGenotypes data from NIF (nif-0000-00096-5)
##NOTE:Genotype view likely not needed, as I can extract genotype, phenotype, and gene data from the Phenotypes view.

### GET ZEBRAFISH DATA ###
#Get Genotype to Phenotype data from NIF (nif-0000-21427-10) for genotype-phenotype associations.
#Get OrganismGenotype data from NIF (nif-0000-21427-11). Like with MGI, this likely won't be needed, as I can extract the implicated genes from the phenotype view for the gene-phenotype associations.

### GET GENE ORTHOLOG DATA ###
#Get Panther:Orthologs data from NIF (nlx_84521-1)
#NOTE: For optimization, will want to only grab data for human, mouse, and zebrafish, so select only rows where Species A/Species B in {Homo sapiens, Mus musculus, Danio rerio}. Should greatly reduce the size of the table.

### GET ADDITIONAL DATA ###
#Depending on time available/data available, will also need to get gene-phenotype data from other model organisms for the additional Phenologs work.





##############    SCRUB DATA    #############
#PURPOSE: Perform any pre-scrubbing that is necessary for each resource/data file.
#Data from NIF sources should fortunately be decently scrubbed (confirm this before proceeding),
# but other resources will likely require additional scrubbing before further processing.


##############    PHENOLOG DATA PRE-PROCESSING    #############
#PURPOSE


