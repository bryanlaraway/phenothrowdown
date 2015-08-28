/*
Author: Bryan Laraway
*/


insheet using "/Users/larawayb/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions_with_rankings.csv"

/*
Assembly of top OWLSim score graphs
*/
by top_owlsim_max_ic_rank, sort: gen freq_1 = _N
by top_owlsim_max_ic_rank: gen cumfreq_1 = _N if _n == 1
gsort -top_owlsim_max_ic_rank
replace cumfreq_1 = sum(cumfreq_1)
gen cumperc_1 = abs(1-(cumfreq_1/n))
label variable cumfreq_1 "Top OWLSim Max IC Score"
twoway (line cumperc_1 top_owlsim_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_owlsim_iccs_rank, sort: gen freq_2 = _N
by top_owlsim_iccs_rank: gen cumfreq_2 = _N if _n == 1
gsort -top_owlsim_iccs_rank
replace cumfreq_2 = sum(cumfreq_2)
gen cumperc_2 = abs(1-(cumfreq_2/n))
label variable cumfreq_2 "Top OWLSim ICCS Score"
twoway (line cumperc_2 top_owlsim_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_owlsim_sim_ic_rank, sort: gen freq_3 = _N
by top_owlsim_sim_ic_rank: gen cumfreq_3 = _N if _n == 1
gsort -top_owlsim_sim_ic_rank
replace cumfreq_3 = sum(cumfreq_3)
gen cumperc_3 = abs(1-(cumfreq_3/n))
label variable cumfreq_3 "Top OWLSim SimIC Score"
twoway (line cumperc_3 top_owlsim_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_owlsim_sim_j_rank, sort: gen freq_4 = _N
by top_owlsim_sim_j_rank: gen cumfreq_4 = _N if _n == 1
gsort -top_owlsim_sim_j_rank
replace cumfreq_4 = sum(cumfreq_4)
gen cumperc_4 = abs(1-(cumfreq_4/n))
label variable cumfreq_4 "Top OWLSim SimJ Score"
twoway (line cumperc_4 top_owlsim_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined top OWLSim graphs */
twoway (scatter cumperc_1 top_owlsim_max_ic_rank) (scatter cumperc_2 top_owlsim_iccs_rank) (scatter cumperc_3 top_owlsim_sim_ic_rank) (scatter cumperc_4 top_owlsim_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_owlsim_max_ic_rank") label(2 "top_owlsim_iccs_rank") label(3 "top_owlsim_sim_ic_rank") label(4 "top_owlsim_sim_j_rank"))


/*
Assembly of top Phenolog score graphs
*/
by top_phenolog_additive_rank, sort: gen freq_5 = _N
by top_phenolog_additive_rank: gen cumfreq_5 = _N if _n == 1
replace cumfreq_5 = sum(cumfreq_5)
gen cumperc_5 = cumfreq_5/n
label variable cumfreq_5 "Top Phenolog Additive Score"
twoway (line cumperc_5 top_phenolog_additive_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_phenolog_max_rank, sort: gen freq_6 = _N
by top_phenolog_max_rank: gen cumfreq_6 = _N if _n == 1
replace cumfreq_6 = sum(cumfreq_6)
gen cumperc_6 = cumfreq_6/n
label variable cumfreq_6 "Top Phenolog Max Score"
twoway (line cumperc_6 top_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

/* Combined top Phenolog graphs */
twoway (scatter cumperc_6 top_phenolog_max_rank) (scatter cumperc_5 top_phenolog_additive_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_phenolog_max_rank") label(2 "top_phenolog_additive_rank"))

/* Combined OWLSim/Phenolog graphs */
twoway (scatter cumperc_1 top_owlsim_max_ic_rank) (scatter cumperc_2 top_owlsim_iccs_rank) (scatter cumperc_3 top_owlsim_sim_ic_rank) (scatter cumperc_4 top_owlsim_sim_j_rank) (scatter cumperc_6 top_phenolog_max_rank) (scatter cumperc_5 top_phenolog_additive_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_owlsim_max_ic_rank") label(2 "top_owlsim_iccs_rank") label(3 "top_owlsim_sim_ic_rank") label(4 "top_owlsim_sim_j_rank") label(5 "top_phenolog_max_rank") label(6 "top_phenolog_additive_rank"))




/*
Assembly of zebrafish LDO OWLSim score graphs
*/
by zebrafish_ldo_max_ic_rank, sort: gen freq_7 = _N
by zebrafish_ldo_max_ic_rank: gen cumfreq_7 = _N if _n == 1
gsort -zebrafish_ldo_max_ic_rank
replace cumfreq_7 = sum(cumfreq_7)
gen cumperc_7 = abs(1-(cumfreq_7/n))
label variable cumfreq_7 "Zebrafish LDO MaxIC Score"
twoway (line cumperc_7 zebrafish_ldo_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by zebrafish_ldo_iccs_rank, sort: gen freq_8 = _N
by zebrafish_ldo_iccs_rank: gen cumfreq_8 = _N if _n == 1
gsort -zebrafish_ldo_iccs_rank
replace cumfreq_8 = sum(cumfreq_8)
gen cumperc_8 = abs(1-(cumfreq_8/n))
label variable cumfreq_8 "Zebrafish LDO ICCS Score"
twoway (line cumperc_8 zebrafish_ldo_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by zebrafish_ldo_sim_ic_rank, sort: gen freq_9 = _N
by zebrafish_ldo_sim_ic_rank: gen cumfreq_9 = _N if _n == 1
gsort -zebrafish_ldo_sim_ic_rank
replace cumfreq_9 = sum(cumfreq_9)
gen cumperc_9 = abs(1-(cumfreq_9/n))
label variable cumfreq_9 "Zebrafish LDO SimIC Score"
twoway (line cumperc_9 zebrafish_ldo_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by zebrafish_ldo_sim_j_rank, sort: gen freq_10 = _N
by zebrafish_ldo_sim_j_rank: gen cumfreq_10 = _N if _n == 1
gsort -zebrafish_ldo_sim_j_rank
replace cumfreq_10 = sum(cumfreq_10)
gen cumperc_10 = abs(1-(cumfreq_10/n))
label variable cumfreq_10 "Zebrafish LDO SimJ Score"
twoway (line cumperc_10 zebrafish_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined zebrafish LDO OWLSim graphs */
twoway (scatter cumperc_7 zebrafish_ldo_max_ic_rank) (scatter cumperc_8 zebrafish_ldo_iccs_rank) (scatter cumperc_9 zebrafish_ldo_sim_ic_rank) (scatter cumperc_10 zebrafish_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "zebrafish_ldo_max_ic_rank") label(2 "zebrafish_ldo_iccs_rank") label(3 "zebrafish_ldo_sim_ic_rank") label(4 "zebrafish_ldo_sim_j_rank"))



/*
Assembly of mouse LDO OWLSim score graphs
*/
by mouse_ldo_max_ic_rank, sort: gen freq_11 = _N
by mouse_ldo_max_ic_rank: gen cumfreq_11 = _N if _n == 1
gsort -mouse_ldo_max_ic_rank
replace cumfreq_11 = sum(cumfreq_11)
gen cumperc_11 = abs(1-(cumfreq_11/n))
label variable cumfreq_11 "Mouse LDO MaxIC Score"
twoway (line cumperc_11 mouse_ldo_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ldo_iccs_rank, sort: gen freq_12 = _N
by mouse_ldo_iccs_rank: gen cumfreq_12 = _N if _n == 1
gsort -mouse_ldo_iccs_rank
replace cumfreq_12 = sum(cumfreq_12)
gen cumperc_12 = abs(1-(cumfreq_12/n))
label variable cumfreq_12 "Mouse LDO ICCS Score"
twoway (line cumperc_12 mouse_ldo_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ldo_sim_ic_rank, sort: gen freq_13 = _N
by mouse_ldo_sim_ic_rank: gen cumfreq_13 = _N if _n == 1
gsort -mouse_ldo_sim_ic_rank
replace cumfreq_13 = sum(cumfreq_13)
gen cumperc_13 = abs(1-(cumfreq_13/n))
label variable cumfreq_13 "Mouse LDO SimIC Score"
twoway (line cumperc_13 mouse_ldo_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ldo_sim_j_rank, sort: gen freq_14 = _N
by mouse_ldo_sim_j_rank: gen cumfreq_14 = _N if _n == 1
gsort -mouse_ldo_sim_j_rank
replace cumfreq_14 = sum(cumfreq_14)
gen cumperc_14 = abs(1-(cumfreq_14/n))
label variable cumfreq_14 "Mouse LDO SimJ Score"
twoway (line cumperc_14 mouse_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined mouse LDO OWLSim graphs */
twoway (scatter cumperc_11 mouse_ldo_max_ic_rank) (scatter cumperc_12 mouse_ldo_iccs_rank) (scatter cumperc_13 mouse_ldo_sim_ic_rank) (scatter cumperc_14 mouse_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "mouse_ldo_max_ic_rank") label(2 "mouse_ldo_iccs_rank") label(3 "mouse_ldo_sim_ic_rank") label(4 "mouse_ldo_sim_j_rank"))





/*
Assembly of zebrafish ortholog OWLSim score graphs
*/
by zebrafish_ortholog_max_ic_rank, sort: gen freq_15 = _N
by zebrafish_ortholog_max_ic_rank: gen cumfreq_15 = _N if _n == 1
gsort -zebrafish_ortholog_max_ic_rank
replace cumfreq_15 = sum(cumfreq_15)
gen cumperc_15 = abs(1-(cumfreq_15/n))
label variable cumfreq_15 "Zebrafish ortholog MaxIC Score"
twoway (line cumperc_15 zebrafish_ortholog_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by zebrafish_ortholog_iccs_rank, sort: gen freq_16 = _N
by zebrafish_ortholog_iccs_rank: gen cumfreq_16 = _N if _n == 1
gsort -zebrafish_ortholog_iccs_rank
replace cumfreq_16 = sum(cumfreq_16)
gen cumperc_16 = abs(1-(cumfreq_16/n))
label variable cumfreq_16 "Zebrafish ortholog ICCS Score"
twoway (line cumperc_16 zebrafish_ortholog_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by zebrafish_ortholog_sim_ic_rank, sort: gen freq_17 = _N
by zebrafish_ortholog_sim_ic_rank: gen cumfreq_17 = _N if _n == 1
gsort -zebrafish_ortholog_sim_ic_rank
replace cumfreq_17 = sum(cumfreq_17)
gen cumperc_17 = abs(1-(cumfreq_17/n))
label variable cumfreq_17 "Zebrafish ortholog SimIC Score"
twoway (line cumperc_17 zebrafish_ortholog_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by zebrafish_ortholog_sim_j_rank, sort: gen freq_18 = _N
by zebrafish_ortholog_sim_j_rank: gen cumfreq_18 = _N if _n == 1
gsort -zebrafish_ortholog_sim_j_rank
replace cumfreq_18 = sum(cumfreq_18)
gen cumperc_18 = abs(1-(cumfreq_18/n))
label variable cumfreq_18 "Zebrafish ortholog SimJ Score"
twoway (line cumperc_18 zebrafish_ortholog_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined zebrafish ortholog OWLSim graphs */
twoway (scatter cumperc_15 zebrafish_ortholog_max_ic_rank) (scatter cumperc_16 zebrafish_ortholog_iccs_rank) (scatter cumperc_17 zebrafish_ortholog_sim_ic_rank) (scatter cumperc_18 zebrafish_ortholog_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "zebrafish_ortholog_max_ic_rank") label(2 "zebrafish_ortholog_iccs_rank") label(3 "zebrafish_ortholog_sim_ic_rank") label(4 "zebrafish_ortholog_sim_j_rank"))



/*
Assembly of mouse ortholog OWLSim score graphs
*/
by mouse_ortholog_max_ic_rank, sort: gen freq_19 = _N
by mouse_ortholog_max_ic_rank: gen cumfreq_19 = _N if _n == 1
gsort -mouse_ortholog_max_ic_rank
replace cumfreq_19 = sum(cumfreq_19)
gen cumperc_19 = abs(1-(cumfreq_19/n))
label variable cumfreq_19 "Mouse ortholog MaxIC Score"
twoway (line cumperc_19 mouse_ortholog_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ortholog_iccs_rank, sort: gen freq_20 = _N
by mouse_ortholog_iccs_rank: gen cumfreq_20 = _N if _n == 1
gsort -mouse_ortholog_iccs_rank
replace cumfreq_20 = sum(cumfreq_20)
gen cumperc_20 = abs(1-(cumfreq_20/n))
label variable cumfreq_20 "Mouse ortholog ICCS Score"
twoway (line cumperc_20 mouse_ortholog_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ortholog_sim_ic_rank, sort: gen freq_21 = _N
by mouse_ortholog_sim_ic_rank: gen cumfreq_21 = _N if _n == 1
gsort -mouse_ortholog_sim_ic_rank
replace cumfreq_21 = sum(cumfreq_21)
gen cumperc_21 = abs(1-(cumfreq_21/n))
label variable cumfreq_21 "Mouse ortholog SimIC Score"
twoway (line cumperc_21 mouse_ortholog_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ortholog_sim_j_rank, sort: gen freq_22 = _N
by mouse_ortholog_sim_j_rank: gen cumfreq_22 = _N if _n == 1
gsort -mouse_ortholog_sim_j_rank
replace cumfreq_22 = sum(cumfreq_22)
gen cumperc_22 = abs(1-(cumfreq_22/n))
label variable cumfreq_22 "Mouse ortholog SimJ Score"
twoway (line cumperc_22 mouse_ortholog_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined mouse ortholog OWLSim graphs */
twoway (scatter cumperc_19 mouse_ortholog_max_ic_rank) (scatter cumperc_20 mouse_ortholog_iccs_rank) (scatter cumperc_21 mouse_ortholog_sim_ic_rank) (scatter cumperc_22 mouse_ortholog_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "mouse_ortholog_max_ic_rank") label(2 "mouse_ortholog_iccs_rank") label(3 "mouse_ortholog_sim_ic_rank") label(4 "mouse_ortholog_sim_j_rank"))





/*
Assembly of top LDO OWLSim score graphs
*/
by top_ldo_max_ic_rank, sort: gen freq_23 = _N
by top_ldo_max_ic_rank: gen cumfreq_23 = _N if _n == 1
gsort -top_ldo_max_ic_rank
replace cumfreq_23 = sum(cumfreq_23)
gen cumperc_23 = abs(1-(cumfreq_23/n))
label variable cumfreq_23 "Top LDO MaxIC Score"
twoway (line cumperc_23 top_ldo_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_ldo_iccs_rank, sort: gen freq_24 = _N
by top_ldo_iccs_rank: gen cumfreq_24 = _N if _n == 1
gsort -top_ldo_iccs_rank
replace cumfreq_24 = sum(cumfreq_24)
gen cumperc_24 = abs(1-(cumfreq_24/n))
label variable cumfreq_24 "Top LDO ICCS Score"
twoway (line cumperc_24 top_ldo_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_ldo_sim_ic_rank, sort: gen freq_25 = _N
by top_ldo_sim_ic_rank: gen cumfreq_25 = _N if _n == 1
gsort -top_ldo_sim_ic_rank
replace cumfreq_25 = sum(cumfreq_25)
gen cumperc_25 = abs(1-(cumfreq_25/n))
label variable cumfreq_25 "Top LDO SimIC Score"
twoway (line cumperc_25 top_ldo_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_ldo_sim_j_rank, sort: gen freq_26 = _N
by top_ldo_sim_j_rank: gen cumfreq_26 = _N if _n == 1
gsort -top_ldo_sim_j_rank
replace cumfreq_26 = sum(cumfreq_26)
gen cumperc_26 = abs(1-(cumfreq_26/n))
label variable cumfreq_26 "Top LDO SimJ Score"
twoway (line cumperc_26 top_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined top LDO OWLSim graphs */
twoway (scatter cumperc_23 top_ldo_max_ic_rank) (scatter cumperc_24 top_ldo_iccs_rank) (scatter cumperc_25 top_ldo_sim_ic_rank) (scatter cumperc_26 top_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_ldo_max_ic_rank") label(2 "top_ldo_iccs_rank") label(3 "top_ldo_sim_ic_rank") label(4 "top_ldo_sim_j_rank"))



/*
Assembly of top ortholog OWLSim score graphs
*/
by top_ortholog_max_ic_rank, sort: gen freq_27 = _N
by top_ortholog_max_ic_rank: gen cumfreq_27 = _N if _n == 1
gsort -top_ortholog_max_ic_rank
replace cumfreq_27 = sum(cumfreq_27)
gen cumperc_27 = abs(1-(cumfreq_27/n))
label variable cumfreq_27 "Top ortholog MaxIC Score"
twoway (line cumperc_27 top_ortholog_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_ortholog_iccs_rank, sort: gen freq_28 = _N
by top_ortholog_iccs_rank: gen cumfreq_28 = _N if _n == 1
gsort -top_ortholog_iccs_rank
replace cumfreq_28 = sum(cumfreq_28)
gen cumperc_28 = abs(1-(cumfreq_28/n))
label variable cumfreq_28 "Top ortholog ICCS Score"
twoway (line cumperc_28 top_ortholog_iccs_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_ldo_sim_ic_rank, sort: gen freq_29 = _N
by top_ldo_sim_ic_rank: gen cumfreq_29 = _N if _n == 1
gsort -top_ldo_sim_ic_rank
replace cumfreq_29 = sum(cumfreq_29)
gen cumperc_29 = abs(1-(cumfreq_29/n))
label variable cumfreq_29 "Top ortholog SimIC Score"
twoway (line cumperc_29 top_ldo_sim_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_ldo_sim_j_rank, sort: gen freq_30 = _N
by top_ldo_sim_j_rank: gen cumfreq_30 = _N if _n == 1
gsort -top_ldo_sim_j_rank
replace cumfreq_30 = sum(cumfreq_30)
gen cumperc_30 = abs(1-(cumfreq_30/n))
label variable cumfreq_30 "Top ortholog SimJ Score"
twoway (line cumperc_30 top_ldo_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/* Combined top ortholog OWLSim graphs */
twoway (scatter cumperc_27 top_ortholog_max_ic_rank) (scatter cumperc_28 top_ortholog_iccs_rank) (scatter cumperc_29 top_ortholog_sim_ic_rank) (scatter cumperc_30 top_ortholog_sim_j_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_ldo_max_ic_rank") label(2 "top_ldo_iccs_rank") label(3 "top_ldo_sim_ic_rank") label(4 "top_ldo_sim_j_rank"))





/*
Assembly of zebrafish max Phenolog score graphs
*/
by zebrafish_ldo_phenolog_max_rank, sort: gen freq_31 = _N
by zebrafish_ldo_phenolog_max_rank: gen cumfreq_31 = _N if _n == 1
replace cumfreq_31 = sum(cumfreq_31)
gen cumperc_31 = cumfreq_31/n
label variable cumfreq_31 "Zebrafish LDO Phenolog Max Score"
twoway (line cumperc_31 zebrafish_ldo_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by v63, sort: gen freq_32 = _N
by v63: gen cumfreq_32 = _N if _n == 1
replace cumfreq_32 = sum(cumfreq_32)
gen cumperc_32 = cumfreq_32/n
label variable cumfreq_32 "Zebrafish ortholog Phenolog Max Score"
twoway (line cumperc_32 v63), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

/*
Assembly of zebrafish additive Phenolog score graphs
*/

by v69, sort: gen freq_33 = _N
by v69: gen cumfreq_33 = _N if _n == 1
replace cumfreq_33 = sum(cumfreq_33)
gen cumperc_33 = cumfreq_33/n
label variable cumfreq_33 "Zebrafish LDO Phenolog Additive Score"
twoway (line cumperc_33 v69), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by v71, sort: gen freq_34 = _N
by v71: gen cumfreq_34 = _N if _n == 1
replace cumfreq_34 = sum(cumfreq_34)
gen cumperc_34 = cumfreq_34/n
label variable cumfreq_34 "Zebrafish ortholog Phenolog Additive Score"
twoway (line cumperc_34 v71), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


/*
Assembly of mouse max Phenolog score graphs
*/
by mouse_ldo_phenolog_max_rank, sort: gen freq_35 = _N
by mouse_ldo_phenolog_max_rank: gen cumfreq_35 = _N if _n == 1
replace cumfreq_35 = sum(cumfreq_35)
gen cumperc_35 = cumfreq_35/n
label variable cumfreq_35 "Mouse LDO Phenolog Max Score"
twoway (line cumperc_35 mouse_ldo_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by mouse_ortholog_phenolog_max_rank, sort: gen freq_36 = _N
by mouse_ortholog_phenolog_max_rank: gen cumfreq_36 = _N if _n == 1
replace cumfreq_36 = sum(cumfreq_36)
gen cumperc_36 = cumfreq_36/n
label variable cumfreq_36 "Mouse ortholog Phenolog Max Score"
twoway (line cumperc_36 mouse_ortholog_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

/*
Assembly of mouse additive Phenolog score graphs
*/

by mouse_ldo_phenolog_additive_rank, sort: gen freq_37 = _N
by mouse_ldo_phenolog_additive_rank: gen cumfreq_37 = _N if _n == 1
replace cumfreq_37 = sum(cumfreq_37)
gen cumperc_37 = cumfreq_37/n
label variable cumfreq_37 "Mouse LDO Phenolog Additive Score"
twoway (line cumperc_37 mouse_ldo_phenolog_additive_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by v75, sort: gen freq_38 = _N
by v75: gen cumfreq_38 = _N if _n == 1
replace cumfreq_38 = sum(cumfreq_38)
gen cumperc_38 = cumfreq_38/n
label variable cumfreq_38 "Mouse ortholog Phenolog Additive Score"
twoway (line cumperc_38 v75), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))












/*
sort top_phenolog_additive_rank
egen n = count(top_phenolog_additive_rank)
by top_phenolog_additive_rank, sort : egen float top_pheno_add_rank_group_size = count(top_phenolog_additive_rank)



by top_phenolog_additive_rank, sort: gen freq = _N
by top_phenolog_additive_rank: gen cumfreq = _N if _n == 1
replace cumfreq = sum(cumfreq)
tabdisp top_phenolog_additive_rank, cell(freq cumfreq)
gen cumperc = cumfreq/n
label variable cumfreq "Top Phenolog Additive Score"
twoway (line cumperc top_phenolog_additive_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits)


by top_phenolog_max_rank, sort: gen max_freq = _N
by top_phenolog_max_rank: gen max_cumfreq = _N if _n == 1
replace max_cumfreq = sum(max_cumfreq)
tabdisp top_phenolog_additive_rank, cell(max_freq max_cumfreq)
gen max_cumperc = max_cumfreq/n
label variable max_cumfreq "Top Phenolog Max Score"
twoway (line max_cumperc top_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits)


by top_owlsim_max_ic_rank, sort: gen maxic_freq = _N
by top_owlsim_max_ic_rank: gen maxic_cumfreq = _N if _n == 1
gsort -top_owlsim_max_ic_rank
replace max_cumfreq = sum(max_cumfreq)
gen max_cumperc = abs(1-(max_cumfreq/n))
label variable max_cumfreq "Top Phenolog Max Score"
twoway (line max_cumperc top_owlsim_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_owlsim_max_ic_rank, sort: gen max_freq = _N
by top_owlsim_max_ic_rank: gen max_cumfreq = _N if _n == 1
gsort -top_owlsim_max_ic_rank
replace max_cumfreq = sum(max_cumfreq)
gen max_cumperc = abs(1-(max_cumfreq/n))
label variable max_cumfreq "Top Phenolog Max Score"
twoway (line max_cumperc top_owlsim_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))


by top_owlsim_max_ic_rank, sort: gen max_freq = _N
by top_owlsim_max_ic_rank: gen max_cumfreq = _N if _n == 1
gsort -top_owlsim_max_ic_rank
replace max_cumfreq = sum(max_cumfreq)
gen max_cumperc = abs(1-(max_cumfreq/n))
label variable max_cumfreq "Top Phenolog Max Score"
twoway (line max_cumperc top_owlsim_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

by top_owlsim_max_ic_rank, sort: gen max_freq = _N
by top_owlsim_max_ic_rank: gen max_cumfreq = _N if _n == 1
gsort -top_owlsim_max_ic_rank
replace max_cumfreq = sum(max_cumfreq)
gen max_cumperc = abs(1-(max_cumfreq/n))
label variable max_cumfreq "Top Phenolog Max Score"
twoway (line max_cumperc top_owlsim_max_ic_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 1000))

gsort -top_owlsim_max_ic_score
by top_owlsim_max_ic_score: gen max_freq = _N
by top_owlsim_max_ic_score: gen top_owlsim_max_ic_cumfreq = _N if _n == 1
replace top_owlsim_max_ic_cumfreq = sum(top_owlsim_max_ic_cumfreq)
replace top_owlsim_max_ic_cumperc = top_owlsim_max_ic_cumfreq/n
label variable top_owlsim_max_ic_cumfreq "Top OWLSim MaxIC Score"
twoway (line top_owlsim_max_ic_cumperc top_owlsim_max_ic_score), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1))


by top_phenolog_max_rank, sort: gen max_freq = _N
by top_phenolog_max_rank: gen max_cumfreq = _N if _n == 1
replace max_cumfreq = sum(max_cumfreq)
tabdisp top_phenolog_additive_rank, cell(max_freq max_cumfreq)
gen max_cumperc = max_cumfreq/n
label variable max_cumfreq "Top Phenolog Max Score"
twoway (line max_cumperc top_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 Top Phenolog Max Rank))


*/
twoway (line max_cumperc top_phenolog_max_rank cumperc top_phenolog_additive_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits) 
twoway (line cumperc top_phenolog_additive_rank) (line max_cumperc top_phenolog_max_rank), ytitle(Percentage recall of human orthologs) xtitle(In top n hits)

/*
max variable length = 32
12345678901234567890123456789012
*/
