/* OMIM_166710  */
insheet using "/Volumes/Time Machine/PycharmProjects/phenothrowdown/out/scatterplot_data/OMIM_166710.csv"

gen recall_max_ic_score = max_ic_score if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_max_ic_score = max_ic_score if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace max_ic_score = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace max_ic_score = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")
gen zfish_max_ic_score = max_ic_score if strpos(gene_candidate_id, "ZFIN")
replace max_ic_score = . if strpos(gene_candidate_id, "ZFIN")

gen recall_iccs_score = iccs_score if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_iccs_score = iccs_score if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace iccs_score = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace iccs_score = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_sim_ic_score = sim_ic_score if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_sim_ic_score = sim_ic_score if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace sim_ic_score = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace sim_ic_score = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_sim_j_score = sim_j_score if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_sim_j_score = sim_j_score if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace sim_j_score = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace sim_j_score = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_phenolog_max_score = phenolog_max_score if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_phenolog_max_score = phenolog_max_score if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace phenolog_max_score = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace phenolog_max_score = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_phenolog_additive_score = phenolog_additive_score if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_phenolog_additive_score = phenolog_additive_score if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace phenolog_additive_score = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace phenolog_additive_score = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

gen zfish_phenolog_additive_score = phenolog_additive_score if strpos(gene_candidate_id, "ZFIN")
replace phenolog_additive_score = . if strpos(gene_candidate_id, "ZFIN")





gen recall_max_ic_rank = max_ic_rank if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_max_ic_rank = max_ic_rank if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace max_ic_rank = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace max_ic_rank = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_iccs_rank = iccs_rank if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_iccs_rank = iccs_rank if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace iccs_rank = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace iccs_rank = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_sim_ic_rank = sim_ic_rank if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_sim_ic_rank = sim_ic_rank if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace sim_ic_rank = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace sim_ic_rank = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_sim_j_rank = sim_j_rank if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_sim_j_rank = sim_j_rank if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace sim_j_rank = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace sim_j_rank = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_phenolog_max_rank = phenolog_max_rank if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_phenolog_max_rank = phenolog_max_rank if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace phenolog_max_rank = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace phenolog_max_rank = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")


gen recall_phenolog_additive_rank = phenolog_additive_rank if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace recall_phenolog_additive_rank = phenolog_additive_rank if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")

replace phenolog_additive_rank = . if inlist(gene_candidate_id , " MGI:88467", "ZFIN:ZDB-GENE-030131-9102", "ZFIN:ZDB-GENE-030131-4400", "MGI:88468", "ZFIN:ZDB-GENE-030131-8415", "MGI:1278315", "ZFIN:ZDB-GENE-050518-2")
replace phenolog_additive_rank = . if inlist(gene_candidate_id , "MGI:103076", "ZFIN:ZDB-GENE-080403-10", "MGI:101950", "ZFIN:ZDB-GENE-060503-420", "MGI:1353470", "ZFIN:ZDB-GENE-070308-5")



twoway (scatter max_ic_score phenolog_max_score, jitter(4) msize(.75)) (scatter recall_max_ic_score recall_phenolog_max_score, jitter(4) msize(.75)), name(maxICvsPhenologMax) title("OWLSim MaxIC vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC) xtitle(Phenolog Max Score) 
twoway (scatter iccs_score phenolog_max_score, jitter(4) msize(.75)) (scatter recall_iccs_score recall_phenolog_max_score, jitter(4) msize(.75)), name(ICCSvsPhenologMax) title("OWLSim ICCS vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS) xtitle(Phenolog Max Score) 
twoway (scatter sim_ic_score phenolog_max_score, jitter(4) msize(.75)) (scatter recall_sim_ic_score recall_phenolog_max_score, jitter(4) msize(.75)), name(SimICvsPhenologMax) title("OWLSim SimIC vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC) xtitle(Phenolog Max Score) 
twoway (scatter sim_j_score phenolog_max_score, jitter(4) msize(.75)) (scatter recall_sim_j_score recall_phenolog_max_score, jitter(4) msize(.75)), name(SimJvsPhenologMax) title("OWLSim SimJ vs Phenolog Max Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ) xtitle(Phenolog Max Score) 

twoway (scatter max_ic_score phenolog_additive_score, jitter(10) msize(1) mcolor(brown) msymbol(circle_hollow)) (scatter zfish_max_ic_score zfish_phenolog_additive_score, jitter(10) msize(1) mcolor(edkblue) msymbol(circle_hollow)) (scatter recall_max_ic_score recall_phenolog_additive_score, jitter(10) msize(1) mcolor(red) msymbol(circle_hollow)), name(MaxICvsPhenologAdd) title("OWLSim MaxIC vs Phenolog Additive Score") legend(label(1 "Mouse Gene Candidate Scores") label(2 "Zebrafish Gene Candidate Scores") label(3 "Known Disease Causative Genes")) ytitle(Max IC) xtitle(Phenolog Additive Score)
twoway (scatter iccs_score phenolog_additive_score, jitter(4) msize(.75)) (scatter recall_iccs_score recall_phenolog_additive_score, jitter(4) msize(.75)), name(ICCSvsPhenologAdd) title("OWLSim ICCS vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS) xtitle(Phenolog Additive Score) 
twoway (scatter sim_ic_score phenolog_additive_score, jitter(4) msize(.75)) (scatter recall_sim_ic_score recall_phenolog_additive_score, jitter(4) msize(.75)), name(SimICvsPhenologAdd) title("OWLSim SimIC vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC) xtitle(Phenolog Additive Score) 
twoway (scatter sim_j_score phenolog_additive_score, jitter(4) msize(.75)) (scatter recall_sim_j_score recall_phenolog_additive_score, jitter(4) msize(.75)), name(SimJvsPhenologAdd) title("OWLSim SimJ vs Phenolog Additive Score") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ) xtitle(Phenolog Additive Score) 



twoway (scatter max_ic_rank phenolog_max_rank, jitter(4) msize(.75)) (scatter recall_max_ic_rank recall_phenolog_max_rank, jitter(4) msize(.75)), name(maxICRankvsPhenologMaxRank) title("OWLSim MaxIC Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC Rank) xtitle(Phenolog Max Rank) 
twoway (scatter iccs_rank phenolog_max_rank, jitter(4) msize(.75)) (scatter recall_iccs_rank recall_phenolog_max_rank, jitter(4) msize(.75)), name(ICCSRankvsPhenologMaxRank) title("OWLSim ICCS Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS Rank) xtitle(Phenolog Max Rank) 
twoway (scatter sim_ic_rank phenolog_max_rank, jitter(4) msize(.75)) (scatter recall_sim_ic_rank recall_phenolog_max_rank, jitter(4) msize(.75)), name(SimICRankvsPhenologMaxRank) title("OWLSim SimIC Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC Rank) xtitle(Phenolog Max Rank) 
twoway (scatter sim_j_rank phenolog_max_rank, jitter(4) msize(.75)) (scatter recall_sim_j_rank recall_phenolog_max_rank, jitter(4) msize(.75)), name(SimJRankvsPhenologMaxRank) title("OWLSim SimJ Rank vs Phenolog Max Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ Rank) xtitle(Phenolog Max Rank) 

twoway (scatter max_ic_rank phenolog_additive_rank, jitter(4) msize(.75)) (scatter recall_max_ic_rank recall_phenolog_additive_rank, jitter(4) msize(.75)), name(maxICRankvsPhenologAddRank) title("OWLSim MaxIC Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Max IC Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter iccs_rank phenolog_additive_rank, jitter(4) msize(.75)) (scatter recall_iccs_rank recall_phenolog_additive_rank, jitter(4) msize(.75)), name(ICCSRankvsPhenologAddRank) title("OWLSim ICCS Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(ICCS Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter sim_ic_rank phenolog_additive_rank, jitter(4) msize(.75)) (scatter recall_sim_ic_rank recall_phenolog_additive_rank, jitter(4) msize(.75)), name(SimICRankvsPhenologAddRank) title("OWLSim SimIC Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(Sim IC Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter sim_j_rank phenolog_additive_rank, jitter(4) msize(.75)) (scatter recall_sim_j_rank recall_phenolog_additive_rank, jitter(4) msize(.75)), name(SimJRankvsPhenologAddRank) title("OWLSim SimJ Rank vs Phenolog Additive Rank") legend(label(1 "Gene Candidate Scores") label(2 "Known Disease Causative Genes")) ytitle(SimJ Rank) xtitle(Phenolog Additive Rank) 
