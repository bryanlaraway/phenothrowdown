/* OMIM_266300  */
insheet using "/Volumes/Time Machine/PycharmProjects/phenothrowdown/out/scatterplot_data/OMIM_266300.csv"

gen recall_max_ic_score = max_ic_score if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace max_ic_score = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_iccs_score = iccs_score if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace iccs_score = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_sim_ic_score = sim_ic_score if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace sim_ic_score = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_sim_j_score = sim_j_score if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace sim_j_score = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_phenolog_max_score = phenolog_max_score if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace phenolog_max_score = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_phenolog_additive_score = phenolog_additive_score if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace phenolog_additive_score = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_max_ic_rank = max_ic_rank if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace max_ic_rank = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_iccs_rank = iccs_rank if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace iccs_rank = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_sim_ic_rank = sim_ic_rank if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace sim_ic_rank = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_sim_j_rank = sim_j_rank if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace sim_j_rank = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_phenolog_max_rank = phenolog_max_rank if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace phenolog_max_rank = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")

gen recall_phenolog_additive_rank = phenolog_additive_rank if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")
replace phenolog_additive_rank = . if inlist(gene_candidate_id , "MGI:99456", "ZFIN:ZDB-GENE-030502-1")





twoway (scatter max_ic_score phenolog_max_score) (scatter recall_max_ic_score recall_phenolog_max_score), name(maxICvsPhenologMax) title("OWLSim MaxIC vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC) xtitle(Phenolog Max Score) 
twoway (scatter iccs_score phenolog_max_score) (scatter recall_iccs_score recall_phenolog_max_score), name(ICCSvsPhenologMax) title("OWLSim ICCS vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS) xtitle(Phenolog Max Score) 
twoway (scatter sim_ic_score phenolog_max_score) (scatter recall_sim_ic_score recall_phenolog_max_score), name(SimICvsPhenologMax) title("OWLSim SimIC vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC) xtitle(Phenolog Max Score) 
twoway (scatter sim_j_score phenolog_max_score) (scatter recall_sim_j_score recall_phenolog_max_score), name(SimJvsPhenologMax) title("OWLSim SimJ vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ) xtitle(Phenolog Max Score) 

twoway (scatter max_ic_score phenolog_additive_score) (scatter recall_max_ic_score recall_phenolog_additive_score), name(maxICvsPhenologAdd) title("OWLSim MaxIC vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC) xtitle(Phenolog Additive Score) 
twoway (scatter iccs_score phenolog_additive_score) (scatter recall_iccs_score recall_phenolog_additive_score), name(ICCSvsPhenologAdd) title("OWLSim ICCS vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS) xtitle(Phenolog Additive Score) 
twoway (scatter sim_ic_score phenolog_additive_score) (scatter recall_sim_ic_score recall_phenolog_additive_score), name(SimICvsPhenologAdd) title("OWLSim SimIC vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC) xtitle(Phenolog Additive Score) 
twoway (scatter sim_j_score phenolog_additive_score) (scatter recall_sim_j_score recall_phenolog_additive_score), name(SimJvsPhenologAdd) title("OWLSim SimJ vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ) xtitle(Phenolog Additive Score) 



twoway (scatter max_ic_rank phenolog_max_rank) (scatter recall_max_ic_rank recall_phenolog_max_rank), name(maxICRankvsPhenologMaxRank) title("OWLSim MaxIC Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC Rank) xtitle(Phenolog Max Rank) 
twoway (scatter iccs_rank phenolog_max_rank) (scatter recall_iccs_rank recall_phenolog_max_rank), name(ICCSRankvsPhenologMaxRank) title("OWLSim ICCS Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS Rank) xtitle(Phenolog Max Rank) 
twoway (scatter sim_ic_rank phenolog_max_rank) (scatter recall_sim_ic_rank recall_phenolog_max_rank), name(SimICRankvsPhenologMaxRank) title("OWLSim SimIC Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC Rank) xtitle(Phenolog Max Rank) 
twoway (scatter sim_j_rank phenolog_max_rank) (scatter recall_sim_j_rank recall_phenolog_max_rank), name(SimJRankvsPhenologMaxRank) title("OWLSim SimJ Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ Rank) xtitle(Phenolog Max Rank) 

twoway (scatter max_ic_rank phenolog_additive_rank) (scatter recall_max_ic_rank recall_phenolog_additive_rank), name(maxICRankvsPhenologAddRank) title("OWLSim MaxIC Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter iccs_rank phenolog_additive_rank) (scatter recall_iccs_rank recall_phenolog_additive_rank), name(ICCSRankvsPhenologAddRank) title("OWLSim ICCS Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter sim_ic_rank phenolog_additive_rank) (scatter recall_sim_ic_rank recall_phenolog_additive_rank), name(SimICRankvsPhenologAddRank) title("OWLSim SimIC Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter sim_j_rank phenolog_additive_rank) (scatter recall_sim_j_rank recall_phenolog_additive_rank), name(SimJRankvsPhenologAddRank) title("OWLSim SimJ Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ Rank) xtitle(Phenolog Additive Rank) 
