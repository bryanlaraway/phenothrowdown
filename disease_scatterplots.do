/*
Author: Bryan Laraway
*/

/* Laptop file location */
insheet using "/Users/larawayb/PycharmProjects/phenothrowdown/out/scatterplot_data/OMIM_260530.csv"

/* Desktop file location */
insheet using "/Volumes/Time Machine/PycharmProjects/phenothrowdown/out/scatterplot_data/OMIM_260530.csv"


twoway (scatter max_ic_score phenolog_max_score), name(maxICvsPhenologMax) title("maxIC vs Phenolog Max Score") ytitle(Max IC) xtitle(Phenolog Max Score) 
twoway (scatter iccs_score phenolog_max_score), name(ICCSvsPhenologMax) title("ICCS vs Phenolog Max Score") ytitle(ICCS) xtitle(Phenolog Max Score) 
twoway (scatter sim_ic_score phenolog_max_score), name(SimICvsPhenologMax) title("SimIC vs Phenolog Max Score") ytitle(Sim IC) xtitle(Phenolog Max Score) 
twoway (scatter sim_j_score phenolog_max_score), name(SimJvsPhenologMax) title("SimJ vs Phenolog Max Score") ytitle(SimJ) xtitle(Phenolog Max Score) 

twoway (scatter max_ic_score phenolog_additive_score), name(maxICvsPhenologAdd) title("maxIC vs Phenolog Additive Score") ytitle(Max IC) xtitle(Phenolog Additive Score) 
twoway (scatter iccs_score phenolog_additive_score), name(ICCSvsPhenologAdd) title("ICCS vs Phenolog Additive Score") ytitle(ICCS) xtitle(Phenolog Additive Score) 
twoway (scatter sim_ic_score phenolog_additive_score), name(SimICvsPhenologAdd) title("SimIC vs Phenolog Additive Score") ytitle(Sim IC) xtitle(Phenolog Additive Score) 
twoway (scatter sim_j_score phenolog_additive_score), name(SimJvsPhenologAdd) title("SimJ vs Phenolog Additive Score") ytitle(SimJ) xtitle(Phenolog Additive Score) 



twoway (scatter max_ic_rank phenolog_max_rank), name(maxICRankvsPhenologMaxRank) title("maxIC Rank vs Phenolog Max Rank") ytitle(Max IC Rank) xtitle(Phenolog Max Rank) 
twoway (scatter iccs_rank phenolog_max_rank), name(ICCSRankvsPhenologMaxRank) title("ICCS Rank vs Phenolog Max Rank") ytitle(ICCS Rank) xtitle(Phenolog Max Rank) 
twoway (scatter sim_ic_rank phenolog_max_rank), name(SimICRankvsPhenologMaxRank) title("SimIC Rank vs Phenolog Max Rank") ytitle(Sim IC Rank) xtitle(Phenolog Max Rank) 
twoway (scatter sim_j_rank phenolog_max_rank), name(SimJRankvsPhenologMaxRank) title("SimJ Rank vs Phenolog Max Rank") ytitle(SimJ Rank) xtitle(Phenolog Max Rank) 

twoway (scatter max_ic_rank phenolog_additive_rank), name(maxICRankvsPhenologAddRank) title("maxIC Rank vs Phenolog Additive Rank") ytitle(Max IC Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter iccs_rank phenolog_additive_rank), name(ICCSRankvsPhenologAddRank) title("ICCS Rank vs Phenolog Additive Rank") ytitle(ICCS Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter sim_ic_rank phenolog_additive_rank), name(SimICRankvsPhenologAddRank) title("SimIC Rank vs Phenolog Additive Rank") ytitle(Sim IC Rank) xtitle(Phenolog Additive Rank) 
twoway (scatter sim_j_rank phenolog_additive_rank), name(SimJRankvsPhenologAddRank) title("SimJ Rank vs Phenolog Additive Rank") ytitle(SimJ Rank) xtitle(Phenolog Additive Rank) 






twoway (scatter max_ic_score iccs_score), name(maxICvsICCS) title("maxIC vs ICCS") ytitle(ICCS) xtitle(Max IC) 
twoway (scatter max_ic_score iccs_score), name(maxICvsICCS) title("maxIC vs ICCS") ytitle(ICCS) xtitle(Max IC) 
twoway (scatter max_ic_score iccs_score), name(maxICvsICCS) title("maxIC vs ICCS") ytitle(ICCS) xtitle(Max IC) 
twoway (scatter max_ic_score iccs_score), name(maxICvsICCS) title("maxIC vs ICCS") ytitle(ICCS) xtitle(Max IC) 
twoway (scatter max_ic_score iccs_score), name(maxICvsICCS) title("maxIC vs ICCS") ytitle(ICCS) xtitle(Max IC) 
twoway (scatter max_ic_score iccs_score), name(maxICvsICCS) title("maxIC vs ICCS") ytitle(ICCS) xtitle(Max IC) 
