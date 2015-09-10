/*
Author: Bryan Laraway
*/

/* Laptop file location */
insheet using "/Users/larawayb/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions_with_rankings.csv"

/* Desktop file location */
insheet using "/Volumes/Time Machine/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions_with_rankings.csv"


/* Desktop file location */
insheet using "/Volumes/Time Machine/PycharmProjects/phenothrowdown/inter/omim/morbid_disease_predictions_with_rankings_k10.csv"

gen ones = 1
gen mouse_count = 3834
gen mouse_ldo_count = 3834
gen mouse_ortholog_count = 3699
gen zebrafish_count = 1694
gen zebrafish_ldo_count = 1662
gen zebrafish_ortholog_count = 1026
gen total_ortholog_count = 3843
gen total_ldo_count = 4014
gen total_record_match_count = 4015

egen n = count(ones)

/* Assembly of top OWLSim score graphs */

by top_owlsim_max_ic_rank, sort: gen freq_1 = _N
by top_owlsim_max_ic_rank: gen cumfreq_1 = _N if _n == 1
replace cumfreq_1 = sum(cumfreq_1)
gen cumperc_1 = abs((cumfreq_1/total_record_match_count))
label variable cumfreq_1 "Top OWLSim Max IC Score"

by top_owlsim_iccs_rank, sort: gen freq_2 = _N
by top_owlsim_iccs_rank: gen cumfreq_2 = _N if _n == 1
replace cumfreq_2 = sum(cumfreq_2)
gen cumperc_2 = abs((cumfreq_2/total_record_match_count))
label variable cumfreq_2 "Top OWLSim ICCS Score"

by top_owlsim_sim_ic_rank, sort: gen freq_3 = _N
by top_owlsim_sim_ic_rank: gen cumfreq_3 = _N if _n == 1
replace cumfreq_3 = sum(cumfreq_3)
gen cumperc_3 = abs((cumfreq_3/total_record_match_count))
label variable cumfreq_3 "Top OWLSim SimIC Score"

by top_owlsim_sim_j_rank, sort: gen freq_4 = _N
by top_owlsim_sim_j_rank: gen cumfreq_4 = _N if _n == 1
replace cumfreq_4 = sum(cumfreq_4)
gen cumperc_4 = abs((cumfreq_4/total_record_match_count))
label variable cumfreq_4 "Top OWLSim SimJ Score"


/* Assembly of top Phenolog score graphs */
by top_phenolog_additive_rank, sort: gen freq_5 = _N
by top_phenolog_additive_rank: gen cumfreq_5 = _N if _n == 1
replace cumfreq_5 = sum(cumfreq_5)
gen cumperc_5 = cumfreq_5/total_record_match_count
label variable cumfreq_5 "Top Phenolog Additive Score"

by top_phenolog_max_rank, sort: gen freq_6 = _N
by top_phenolog_max_rank: gen cumfreq_6 = _N if _n == 1
replace cumfreq_6 = sum(cumfreq_6)
gen cumperc_6 = cumfreq_6/total_record_match_count
label variable cumfreq_6 "Top Phenolog Max Score"


/* Assembly of zebrafish LDO OWLSim score graphs */
by zebrafish_ldo_max_ic_rank, sort: gen freq_7 = _N
by zebrafish_ldo_max_ic_rank: gen cumfreq_7 = _N if _n == 1
replace cumfreq_7 = sum(cumfreq_7)
gen cumperc_7 = abs((cumfreq_7/zebrafish_ldo_count))
label variable cumfreq_7 "Zebrafish LDO MaxIC Score"

by zebrafish_ldo_iccs_rank, sort: gen freq_8 = _N
by zebrafish_ldo_iccs_rank: gen cumfreq_8 = _N if _n == 1
replace cumfreq_8 = sum(cumfreq_8)
gen cumperc_8 = abs((cumfreq_8/zebrafish_ldo_count))
label variable cumfreq_8 "Zebrafish LDO ICCS Score"

by zebrafish_ldo_sim_ic_rank, sort: gen freq_9 = _N
by zebrafish_ldo_sim_ic_rank: gen cumfreq_9 = _N if _n == 1
replace cumfreq_9 = sum(cumfreq_9)
gen cumperc_9 = abs((cumfreq_9/zebrafish_ldo_count))
label variable cumfreq_9 "Zebrafish LDO SimIC Score"

by zebrafish_ldo_sim_j_rank, sort: gen freq_10 = _N
by zebrafish_ldo_sim_j_rank: gen cumfreq_10 = _N if _n == 1
replace cumfreq_10 = sum(cumfreq_10)
gen cumperc_10 = abs((cumfreq_10/zebrafish_ldo_count))
label variable cumfreq_10 "Zebrafish LDO SimJ Score"


/* Assembly of mouse LDO OWLSim score graphs */
by mouse_ldo_max_ic_rank, sort: gen freq_11 = _N
by mouse_ldo_max_ic_rank: gen cumfreq_11 = _N if _n == 1
replace cumfreq_11 = sum(cumfreq_11)
gen cumperc_11 = abs((cumfreq_11/mouse_ldo_count))
label variable cumfreq_11 "Mouse LDO MaxIC Score"

by mouse_ldo_iccs_rank, sort: gen freq_12 = _N
by mouse_ldo_iccs_rank: gen cumfreq_12 = _N if _n == 1
replace cumfreq_12 = sum(cumfreq_12)
gen cumperc_12 = abs((cumfreq_12/mouse_ldo_count))
label variable cumfreq_12 "Mouse LDO ICCS Score"

by mouse_ldo_sim_ic_rank, sort: gen freq_13 = _N
by mouse_ldo_sim_ic_rank: gen cumfreq_13 = _N if _n == 1
replace cumfreq_13 = sum(cumfreq_13)
gen cumperc_13 = abs((cumfreq_13/mouse_ldo_count))
label variable cumfreq_13 "Mouse LDO SimIC Score"

by mouse_ldo_sim_j_rank, sort: gen freq_14 = _N
by mouse_ldo_sim_j_rank: gen cumfreq_14 = _N if _n == 1
replace cumfreq_14 = sum(cumfreq_14)
gen cumperc_14 = abs((cumfreq_14/mouse_ldo_count))
label variable cumfreq_14 "Mouse LDO SimJ Score"


/* Assembly of zebrafish ortholog OWLSim score graphs */
by zebrafish_ortholog_max_ic_rank, sort: gen freq_15 = _N
by zebrafish_ortholog_max_ic_rank: gen cumfreq_15 = _N if _n == 1
replace cumfreq_15 = sum(cumfreq_15)
gen cumperc_15 = abs((cumfreq_15/zebrafish_ortholog_count))
label variable cumfreq_15 "Zebrafish ortholog MaxIC Score"

by zebrafish_ortholog_iccs_rank, sort: gen freq_16 = _N
by zebrafish_ortholog_iccs_rank: gen cumfreq_16 = _N if _n == 1
replace cumfreq_16 = sum(cumfreq_16)
gen cumperc_16 = abs((cumfreq_16/zebrafish_ortholog_count))
label variable cumfreq_16 "Zebrafish ortholog ICCS Score"

by zebrafish_ortholog_sim_ic_rank, sort: gen freq_17 = _N
by zebrafish_ortholog_sim_ic_rank: gen cumfreq_17 = _N if _n == 1
replace cumfreq_17 = sum(cumfreq_17)
gen cumperc_17 = abs((cumfreq_17/zebrafish_ortholog_count))
label variable cumfreq_17 "Zebrafish ortholog SimIC Score"

by zebrafish_ortholog_sim_j_rank, sort: gen freq_18 = _N
by zebrafish_ortholog_sim_j_rank: gen cumfreq_18 = _N if _n == 1
replace cumfreq_18 = sum(cumfreq_18)
gen cumperc_18 = abs((cumfreq_18/zebrafish_ortholog_count))
label variable cumfreq_18 "Zebrafish ortholog SimJ Score"


/* Assembly of mouse ortholog OWLSim score graphs */
by mouse_ortholog_max_ic_rank, sort: gen freq_19 = _N
by mouse_ortholog_max_ic_rank: gen cumfreq_19 = _N if _n == 1
replace cumfreq_19 = sum(cumfreq_19)
gen cumperc_19 = abs((cumfreq_19/mouse_ortholog_count))
label variable cumfreq_19 "Mouse ortholog MaxIC Score"

by mouse_ortholog_iccs_rank, sort: gen freq_20 = _N
by mouse_ortholog_iccs_rank: gen cumfreq_20 = _N if _n == 1
replace cumfreq_20 = sum(cumfreq_20)
gen cumperc_20 = abs((cumfreq_20/mouse_ortholog_count))
label variable cumfreq_20 "Mouse ortholog ICCS Score"

by mouse_ortholog_sim_ic_rank, sort: gen freq_21 = _N
by mouse_ortholog_sim_ic_rank: gen cumfreq_21 = _N if _n == 1
replace cumfreq_21 = sum(cumfreq_21)
gen cumperc_21 = abs((cumfreq_21/mouse_ortholog_count))
label variable cumfreq_21 "Mouse ortholog SimIC Score"

by mouse_ortholog_sim_j_rank, sort: gen freq_22 = _N
by mouse_ortholog_sim_j_rank: gen cumfreq_22 = _N if _n == 1
replace cumfreq_22 = sum(cumfreq_22)
gen cumperc_22 = abs((cumfreq_22/mouse_ortholog_count))
label variable cumfreq_22 "Mouse ortholog SimJ Score"


/* Assembly of top LDO OWLSim score graphs */
by top_ldo_max_ic_rank, sort: gen freq_23 = _N
by top_ldo_max_ic_rank: gen cumfreq_23 = _N if _n == 1
replace cumfreq_23 = sum(cumfreq_23)
gen cumperc_23 = abs((cumfreq_23/total_ldo_count))
label variable cumfreq_23 "Top LDO MaxIC Score"

by top_ldo_iccs_rank, sort: gen freq_24 = _N
by top_ldo_iccs_rank: gen cumfreq_24 = _N if _n == 1
replace cumfreq_24 = sum(cumfreq_24)
gen cumperc_24 = abs((cumfreq_24/total_ldo_count))
label variable cumfreq_24 "Top LDO ICCS Score"

by top_ldo_sim_ic_rank, sort: gen freq_25 = _N
by top_ldo_sim_ic_rank: gen cumfreq_25 = _N if _n == 1
replace cumfreq_25 = sum(cumfreq_25)
gen cumperc_25 = abs((cumfreq_25/total_ldo_count))
label variable cumfreq_25 "Top LDO SimIC Score"

by top_ldo_sim_j_rank, sort: gen freq_26 = _N
by top_ldo_sim_j_rank: gen cumfreq_26 = _N if _n == 1
replace cumfreq_26 = sum(cumfreq_26)
gen cumperc_26 = abs((cumfreq_26/total_ldo_count))
label variable cumfreq_26 "Top LDO SimJ Score"


/* Assembly of top ortholog OWLSim score graphs */
by top_ortholog_max_ic_rank, sort: gen freq_27 = _N
by top_ortholog_max_ic_rank: gen cumfreq_27 = _N if _n == 1
replace cumfreq_27 = sum(cumfreq_27)
gen cumperc_27 = abs((cumfreq_27/total_ortholog_count))
label variable cumfreq_27 "Top ortholog MaxIC Score"

by top_ortholog_iccs_rank, sort: gen freq_28 = _N
by top_ortholog_iccs_rank: gen cumfreq_28 = _N if _n == 1
replace cumfreq_28 = sum(cumfreq_28)
gen cumperc_28 = abs((cumfreq_28/total_ortholog_count))
label variable cumfreq_28 "Top ortholog ICCS Score"

by top_ortholog_sim_ic_rank, sort: gen freq_29 = _N
by top_ortholog_sim_ic_rank: gen cumfreq_29 = _N if _n == 1
replace cumfreq_29 = sum(cumfreq_29)
gen cumperc_29 = abs((cumfreq_29/total_ortholog_count))
label variable cumfreq_29 "Top ortholog SimIC Score"

by top_ortholog_sim_j_rank, sort: gen freq_30 = _N
by top_ortholog_sim_j_rank: gen cumfreq_30 = _N if _n == 1
replace cumfreq_30 = sum(cumfreq_30)
gen cumperc_30 = abs((cumfreq_30/total_ortholog_count))
label variable cumfreq_30 "Top ortholog SimJ Score"


/* Assembly of top zebrafish OWLSim score graphs */
by top_zebrafish_max_ic_rank, sort: gen freq_43 = _N
by top_zebrafish_max_ic_rank: gen cumfreq_43 = _N if _n == 1
replace cumfreq_43 = sum(cumfreq_43)
gen cumperc_43 = abs(cumfreq_43/zebrafish_count)
label variable cumfreq_43 "Top zebrafish MaxIC Score"

by top_zebrafish_iccs_rank, sort: gen freq_44 = _N
by top_zebrafish_iccs_rank: gen cumfreq_44 = _N if _n == 1
replace cumfreq_44 = sum(cumfreq_44)
gen cumperc_44 = abs((cumfreq_44/zebrafish_count))
label variable cumfreq_44 "Top zebrafish ICCS Score"

by top_zebrafish_sim_ic_rank, sort: gen freq_45 = _N
by top_zebrafish_sim_ic_rank: gen cumfreq_45 = _N if _n == 1
replace cumfreq_45 = sum(cumfreq_45)
gen cumperc_45 = abs((cumfreq_45/zebrafish_count))
label variable cumfreq_45 "Top zebrafish SimIC Score"

by top_zebrafish_sim_j_rank, sort: gen freq_46 = _N
by top_zebrafish_sim_j_rank: gen cumfreq_46 = _N if _n == 1
replace cumfreq_46 = sum(cumfreq_46)
gen cumperc_46 = abs((cumfreq_46/zebrafish_count))
label variable cumfreq_46 "Top zebrafish SimJ Score"


/* Assembly of top mouse OWLSim score graphs */
by top_mouse_max_ic_rank, sort: gen freq_47 = _N
by top_mouse_max_ic_rank: gen cumfreq_47 = _N if _n == 1
replace cumfreq_47 = sum(cumfreq_47)
gen cumperc_47 = abs((cumfreq_47/mouse_count))
label variable cumfreq_47 "Top mouse MaxIC Score"

by top_mouse_iccs_rank, sort: gen freq_48 = _N
by top_mouse_iccs_rank: gen cumfreq_48 = _N if _n == 1
replace cumfreq_48 = sum(cumfreq_48)
gen cumperc_48 = abs((cumfreq_48/mouse_count))
label variable cumfreq_48 "Top mouse ICCS Score"

by top_mouse_sim_ic_rank, sort: gen freq_49 = _N
by top_mouse_sim_ic_rank: gen cumfreq_49 = _N if _n == 1
replace cumfreq_49 = sum(cumfreq_49)
gen cumperc_49 = abs((cumfreq_49/mouse_count))
label variable cumfreq_49 "Top mouse SimIC Score"

by top_mouse_sim_j_rank, sort: gen freq_50 = _N
by top_mouse_sim_j_rank: gen cumfreq_50 = _N if _n == 1
replace cumfreq_50 = sum(cumfreq_50)
gen cumperc_50 = abs((cumfreq_50/mouse_count))
label variable cumfreq_50 "Top mouse SimJ Score"


/* Assembly of zebrafish max Phenolog score graphs */
by zebrafish_ldo_phenolog_max_rank, sort: gen freq_31 = _N
by zebrafish_ldo_phenolog_max_rank: gen cumfreq_31 = _N if _n == 1
replace cumfreq_31 = sum(cumfreq_31)
gen cumperc_31 = cumfreq_31/zebrafish_ldo_count
label variable cumfreq_31 "Zebrafish LDO Phenolog Max Score"

by v79, sort: gen freq_32 = _N
by v79: gen cumfreq_32 = _N if _n == 1
replace cumfreq_32 = sum(cumfreq_32)
gen cumperc_32 = cumfreq_32/zebrafish_ortholog_count
label variable cumfreq_32 "Zebrafish ortholog Phenolog Max Score"

/* Assembly of zebrafish additive Phenolog score graphs */
by v85, sort: gen freq_33 = _N
by v85: gen cumfreq_33 = _N if _n == 1
replace cumfreq_33 = sum(cumfreq_33)
gen cumperc_33 = cumfreq_33/zebrafish_ldo_count
label variable cumfreq_33 "Zebrafish LDO Phenolog Additive Score"

by v87, sort: gen freq_34 = _N
by v87: gen cumfreq_34 = _N if _n == 1
replace cumfreq_34 = sum(cumfreq_34)
gen cumperc_34 = cumfreq_34/zebrafish_ortholog_count
label variable cumfreq_34 "Zebrafish ortholog Phenolog Additive Score"


/* Assembly of mouse max Phenolog score graphs */
by mouse_ldo_phenolog_max_rank, sort: gen freq_35 = _N
by mouse_ldo_phenolog_max_rank: gen cumfreq_35 = _N if _n == 1
replace cumfreq_35 = sum(cumfreq_35)
gen cumperc_35 = cumfreq_35/mouse_ldo_count
label variable cumfreq_35 "Mouse LDO Phenolog Max Score"

by mouse_ortholog_phenolog_max_rank, sort: gen freq_36 = _N
by mouse_ortholog_phenolog_max_rank: gen cumfreq_36 = _N if _n == 1
replace cumfreq_36 = sum(cumfreq_36)
gen cumperc_36 = cumfreq_36/mouse_ortholog_count
label variable cumfreq_36 "Mouse ortholog Phenolog Max Score"


/* Assembly of mouse additive Phenolog score graphs */
by mouse_ldo_phenolog_additive_rank, sort: gen freq_37 = _N
by mouse_ldo_phenolog_additive_rank: gen cumfreq_37 = _N if _n == 1
replace cumfreq_37 = sum(cumfreq_37)
gen cumperc_37 = cumfreq_37/mouse_ldo_count
label variable cumfreq_37 "Mouse LDO Phenolog Additive Score"

by v91, sort: gen freq_38 = _N
by v91: gen cumfreq_38 = _N if _n == 1
replace cumfreq_38 = sum(cumfreq_38)
gen cumperc_38 = cumfreq_38/mouse_ortholog_count
label variable cumfreq_38 "Mouse ortholog Phenolog Additive Score"


/* Assembly of top zebrafish Phenolog score graphs */
by top_zebrafish_phenolog_max_rank, sort: gen freq_39 = _N
by top_zebrafish_phenolog_max_rank: gen cumfreq_39 = _N if _n == 1
replace cumfreq_39 = sum(cumfreq_39)
gen cumperc_39 = cumfreq_39/zebrafish_count
label variable cumfreq_39 "Top Zebrafish Phenolog Max Score"

by v95, sort: gen freq_40 = _N
by v95: gen cumfreq_40 = _N if _n == 1
replace cumfreq_40 = sum(cumfreq_40)
gen cumperc_40 = cumfreq_40/zebrafish_count
label variable cumfreq_40 "Top Zebrafish Phenolog Additive Score"


/* Assembly of top mouse Phenolog score graphs */
by top_mouse_phenolog_max_rank, sort: gen freq_41 = _N
by top_mouse_phenolog_max_rank: gen cumfreq_41 = _N if _n == 1
replace cumfreq_41 = sum(cumfreq_41)
gen cumperc_41 = cumfreq_41/mouse_count
label variable cumfreq_41 "Top Mouse Phenolog Max Score"

by top_mouse_phenolog_additive_rank, sort: gen freq_42 = _N
by top_mouse_phenolog_additive_rank: gen cumfreq_42 = _N if _n == 1
replace cumfreq_42 = sum(cumfreq_42)
gen cumperc_42 = cumfreq_42/mouse_count
label variable cumfreq_42 "Top Mouse Phenolog Additive Score"


/* Assembly of top Phenolog ldo score graphs */
by top_ldo_phenolog_max_rank, sort: gen freq_51 = _N
by top_ldo_phenolog_max_rank: gen cumfreq_51 = _N if _n == 1
replace cumfreq_51 = sum(cumfreq_51)
gen cumperc_51 = cumfreq_51/total_ldo_count
label variable cumfreq_51 "Top Phenolog LDO Max Score"

by top_ldo_phenolog_additive_rank, sort: gen freq_52 = _N
by top_ldo_phenolog_additive_rank: gen cumfreq_52 = _N if _n == 1
replace cumfreq_52 = sum(cumfreq_52)
gen cumperc_52 = cumfreq_52/total_ldo_count
label variable cumfreq_52 "Top Phenolog LDO Additive Score"


/* Assembly of top Phenolog ortholog score graphs */
by top_ortholog_phenolog_max_rank, sort: gen freq_53 = _N
by top_ortholog_phenolog_max_rank: gen cumfreq_53 = _N if _n == 1
replace cumfreq_53 = sum(cumfreq_53)
gen cumperc_53 = cumfreq_53/total_ortholog_count

label variable cumfreq_53 "Top Phenolog Ortholog Max Score"

by top_ortholog_phenolog_additive_r, sort: gen freq_54 = _N
by top_ortholog_phenolog_additive_r: gen cumfreq_54 = _N if _n == 1
replace cumfreq_54 = sum(cumfreq_54)
gen cumperc_54 = cumfreq_54/total_ortholog_count
label variable cumfreq_54 "Top Phenolog Ortholog Additive Score"


/* Execute to limit to rank >= 500 */
replace cumperc_1 = . if top_owlsim_max_ic_rank > 500
replace top_owlsim_max_ic_rank = . if top_owlsim_max_ic_rank > 500
replace cumperc_2 = . if top_owlsim_iccs_rank > 500
replace top_owlsim_iccs_rank = . if top_owlsim_iccs_rank > 500
replace cumperc_3 = . if top_owlsim_sim_ic_rank > 500
replace top_owlsim_sim_ic_rank = . if top_owlsim_sim_ic_rank > 500
replace cumperc_4 = . if top_owlsim_sim_j_rank > 500
replace top_owlsim_sim_j_rank = . if top_owlsim_sim_j_rank > 500
replace cumperc_5 = . if top_phenolog_additive_rank > 500
replace top_phenolog_additive_rank = . if top_phenolog_additive_rank > 500
replace cumperc_6 = . if top_phenolog_max_rank > 500
replace top_phenolog_max_rank = . if top_phenolog_max_rank > 500
replace cumperc_7 = . if zebrafish_ldo_max_ic_rank > 500
replace zebrafish_ldo_max_ic_rank = . if zebrafish_ldo_max_ic_rank > 500
replace cumperc_8 = . if zebrafish_ldo_iccs_rank > 500
replace zebrafish_ldo_iccs_rank = . if zebrafish_ldo_iccs_rank > 500
replace cumperc_9 = . if zebrafish_ldo_sim_ic_rank > 500
replace zebrafish_ldo_sim_ic_rank = . if zebrafish_ldo_sim_ic_rank > 500
replace cumperc_10 = . if zebrafish_ldo_sim_j_rank > 500
replace zebrafish_ldo_sim_j_rank = . if zebrafish_ldo_sim_j_rank > 500
replace cumperc_11 = . if mouse_ldo_max_ic_rank > 500
replace mouse_ldo_max_ic_rank = . if mouse_ldo_max_ic_rank > 500
replace cumperc_12 = . if mouse_ldo_iccs_rank > 500
replace mouse_ldo_iccs_rank = . if mouse_ldo_iccs_rank > 500
replace cumperc_13 = . if mouse_ldo_sim_ic_rank > 500
replace mouse_ldo_sim_ic_rank = . if mouse_ldo_sim_ic_rank > 500
replace cumperc_14 = . if mouse_ldo_sim_j_rank > 500
replace mouse_ldo_sim_j_rank = . if mouse_ldo_sim_j_rank > 500
replace cumperc_15 = . if zebrafish_ortholog_max_ic_rank > 500
replace zebrafish_ortholog_max_ic_rank = . if zebrafish_ortholog_max_ic_rank > 500
replace cumperc_16 = . if zebrafish_ortholog_iccs_rank > 500
replace zebrafish_ortholog_iccs_rank = . if zebrafish_ortholog_iccs_rank > 500
replace cumperc_17 = . if zebrafish_ortholog_sim_ic_rank > 500
replace zebrafish_ortholog_sim_ic_rank = . if zebrafish_ortholog_sim_ic_rank > 500
replace cumperc_18 = . if zebrafish_ortholog_sim_j_rank > 500
replace zebrafish_ortholog_sim_j_rank = . if zebrafish_ortholog_sim_j_rank > 500
replace cumperc_19 = . if mouse_ortholog_max_ic_rank > 500
replace mouse_ortholog_max_ic_rank = . if mouse_ortholog_max_ic_rank > 500
replace cumperc_20 = . if mouse_ortholog_iccs_rank > 500
replace mouse_ortholog_iccs_rank = . if mouse_ortholog_iccs_rank > 500
replace cumperc_21 = . if mouse_ortholog_sim_ic_rank > 500
replace mouse_ortholog_sim_ic_rank = . if mouse_ortholog_sim_ic_rank > 500
replace cumperc_22 = . if mouse_ortholog_sim_j_rank > 500
replace mouse_ortholog_sim_j_rank = . if mouse_ortholog_sim_j_rank > 500
replace cumperc_23 = . if top_ldo_max_ic_rank > 500
replace top_ldo_max_ic_rank = . if top_ldo_max_ic_rank > 500
replace cumperc_24 = . if top_ldo_iccs_rank > 500
replace top_ldo_iccs_rank = . if top_ldo_iccs_rank > 500
replace cumperc_25 = . if top_ldo_sim_ic_rank > 500
replace top_ldo_sim_ic_rank = . if top_ldo_sim_ic_rank > 500
replace cumperc_26 = . if top_ldo_sim_j_rank > 500
replace top_ldo_sim_j_rank = . if top_ldo_sim_j_rank > 500
replace cumperc_27 = . if top_ortholog_max_ic_rank > 500
replace top_ortholog_max_ic_rank = . if top_ortholog_max_ic_rank > 500
replace cumperc_28 = . if top_ortholog_iccs_rank > 500
replace top_ortholog_iccs_rank = . if top_ortholog_iccs_rank > 500
replace cumperc_29 = . if top_ortholog_sim_ic_rank > 500
replace top_ortholog_sim_ic_rank = . if top_ortholog_sim_ic_rank > 500
replace cumperc_30 = . if top_ortholog_sim_j_rank > 500
replace top_ortholog_sim_j_rank = . if top_ortholog_sim_j_rank > 500
replace cumperc_43 = . if top_zebrafish_max_ic_rank > 500
replace top_zebrafish_max_ic_rank = . if top_zebrafish_max_ic_rank > 500
replace cumperc_44 = . if top_zebrafish_iccs_rank > 500
replace top_zebrafish_iccs_rank = . if top_zebrafish_iccs_rank > 500
replace cumperc_45 = . if top_zebrafish_sim_ic_rank > 500
replace top_zebrafish_sim_ic_rank = . if top_zebrafish_sim_ic_rank > 500
replace cumperc_46 = . if top_zebrafish_sim_j_rank > 500
replace top_zebrafish_sim_j_rank = . if top_zebrafish_sim_j_rank > 500
replace cumperc_47 = . if top_mouse_max_ic_rank > 500
replace top_mouse_max_ic_rank = . if top_mouse_max_ic_rank > 500
replace cumperc_48 = . if top_mouse_iccs_rank > 500
replace top_mouse_iccs_rank = . if top_mouse_iccs_rank > 500
replace cumperc_49 = . if top_mouse_sim_ic_rank > 500
replace top_mouse_sim_ic_rank = . if top_mouse_sim_ic_rank > 500
replace cumperc_50 = . if top_mouse_sim_j_rank > 500
replace top_mouse_sim_j_rank = . if top_mouse_sim_j_rank > 500
replace cumperc_31 = . if zebrafish_ldo_phenolog_max_rank > 500
replace zebrafish_ldo_phenolog_max_rank = . if zebrafish_ldo_phenolog_max_rank > 500
replace cumperc_32 = . if v79 > 500
replace v79 = . if v79 > 500
replace cumperc_33 = . if v85 > 500
replace v85 = . if v85 > 500
replace cumperc_34 = . if v87 > 500
replace v87 = . if v87 > 500
replace cumperc_35 = . if mouse_ldo_phenolog_max_rank > 500
replace mouse_ldo_phenolog_max_rank = . if mouse_ldo_phenolog_max_rank > 500
replace cumperc_36 = . if mouse_ortholog_phenolog_max_rank > 500
replace mouse_ortholog_phenolog_max_rank = . if mouse_ortholog_phenolog_max_rank > 500
replace cumperc_37 = . if mouse_ldo_phenolog_additive_rank > 500
replace mouse_ldo_phenolog_additive_rank = . if mouse_ldo_phenolog_additive_rank > 500
replace cumperc_38 = . if v91 > 500
replace v91 = . if v91 > 500
replace cumperc_39 = . if top_zebrafish_phenolog_max_rank > 500
replace top_zebrafish_phenolog_max_rank = . if top_zebrafish_phenolog_max_rank > 500
replace cumperc_40 = . if v95 > 500
replace v95 = . if v95 > 500
replace cumperc_41 = . if top_mouse_phenolog_max_rank > 500
replace top_mouse_phenolog_max_rank = . if top_mouse_phenolog_max_rank > 500
replace cumperc_42 = . if top_mouse_phenolog_additive_rank > 500
replace top_mouse_phenolog_additive_rank = . if top_mouse_phenolog_additive_rank > 500
replace cumperc_51 = . if top_ldo_phenolog_max_rank > 500
replace top_ldo_phenolog_max_rank = . if top_ldo_phenolog_max_rank > 500
replace cumperc_52 = . if top_ldo_phenolog_additive_rank > 500
replace top_ldo_phenolog_additive_rank = . if top_ldo_phenolog_additive_rank > 500
replace cumperc_53 = . if top_ortholog_phenolog_max_rank > 500
replace top_ortholog_phenolog_max_rank = . if top_ortholog_phenolog_max_rank > 500
replace cumperc_54 = . if top_ortholog_phenolog_additive_r > 500
replace top_ortholog_phenolog_additive_r = . if top_ortholog_phenolog_additive_r > 500





set obs 4236
replace cumperc_1 = 0 in 4236
replace top_owlsim_max_ic_rank = 0 in 4236
replace cumperc_2 = 0 in 4236
replace top_owlsim_iccs_rank = 0 in 4236
replace cumperc_3 = 0 in 4236
replace top_owlsim_sim_ic_rank = 0 in 4236
replace cumperc_4 = 0 in 4236
replace top_owlsim_sim_j_rank = 0 in 4236
replace cumperc_5 = 0 in 4236
replace top_phenolog_additive_rank = 0 in 4236
replace cumperc_6 = 0 in 4236
replace top_phenolog_max_rank = 0 in 4236
replace cumperc_7 = 0 in 4236
replace zebrafish_ldo_max_ic_rank = 0 in 4236
replace cumperc_8 = 0 in 4236
replace zebrafish_ldo_iccs_rank = 0 in 4236
replace cumperc_9 = 0 in 4236
replace zebrafish_ldo_sim_ic_rank = 0 in 4236
replace cumperc_10 = 0 in 4236
replace zebrafish_ldo_sim_j_rank = 0 in 4236
replace cumperc_11 = 0 in 4236
replace mouse_ldo_max_ic_rank = 0 in 4236
replace cumperc_12 = 0 in 4236
replace mouse_ldo_iccs_rank = 0 in 4236
replace cumperc_13 = 0 in 4236
replace mouse_ldo_sim_ic_rank = 0 in 4236
replace cumperc_14 = 0 in 4236
replace mouse_ldo_sim_j_rank = 0 in 4236
replace cumperc_15 = 0 in 4236
replace zebrafish_ortholog_max_ic_rank = 0 in 4236
replace cumperc_16 = 0 in 4236
replace zebrafish_ortholog_iccs_rank = 0 in 4236
replace cumperc_17 = 0 in 4236
replace zebrafish_ortholog_sim_ic_rank = 0 in 4236
replace cumperc_18 = 0 in 4236
replace zebrafish_ortholog_sim_j_rank = 0 in 4236
replace cumperc_19 = 0 in 4236
replace mouse_ortholog_max_ic_rank = 0 in 4236
replace cumperc_20 = 0 in 4236
replace mouse_ortholog_iccs_rank = 0 in 4236
replace cumperc_21 = 0 in 4236
replace mouse_ortholog_sim_ic_rank = 0 in 4236
replace cumperc_22 = 0 in 4236
replace mouse_ortholog_sim_j_rank = 0 in 4236
replace cumperc_23 = 0 in 4236
replace top_ldo_max_ic_rank = 0 in 4236
replace cumperc_24 = 0 in 4236
replace top_ldo_iccs_rank = 0 in 4236
replace cumperc_25 = 0 in 4236
replace top_ldo_sim_ic_rank = 0 in 4236
replace cumperc_26 = 0 in 4236
replace top_ldo_sim_j_rank = 0 in 4236
replace cumperc_27 = 0 in 4236
replace top_ortholog_max_ic_rank = 0 in 4236
replace cumperc_28 = 0 in 4236
replace top_ortholog_iccs_rank = 0 in 4236
replace cumperc_29 = 0 in 4236
replace top_ldo_sim_ic_rank = 0 in 4236
replace cumperc_30 = 0 in 4236
replace top_ldo_sim_j_rank = 0 in 4236
replace cumperc_43 = 0 in 4236
replace top_zebrafish_max_ic_rank = 0 in 4236
replace cumperc_44 = 0 in 4236
replace top_zebrafish_iccs_rank = 0 in 4236
replace cumperc_45 = 0 in 4236
replace top_zebrafish_sim_ic_rank = 0 in 4236
replace cumperc_46 = 0 in 4236
replace top_zebrafish_sim_j_rank = 0 in 4236
replace cumperc_47 = 0 in 4236
replace top_mouse_max_ic_rank = 0 in 4236
replace cumperc_48 = 0 in 4236
replace top_mouse_iccs_rank = 0 in 4236
replace cumperc_49 = 0 in 4236
replace top_mouse_sim_ic_rank = 0 in 4236
replace cumperc_50 = 0 in 4236
replace top_mouse_sim_j_rank = 0 in 4236
replace cumperc_31 = 0 in 4236
replace zebrafish_ldo_phenolog_max_rank = 0 in 4236
replace cumperc_32 = 0 in 4236
replace v79 = 0 in 4236
replace cumperc_33 = 0 in 4236
replace v85 = 0 in 4236
replace cumperc_34 = 0 in 4236
replace v87 = 0 in 4236
replace cumperc_35 = 0 in 4236
replace mouse_ldo_phenolog_max_rank = 0 in 4236
replace cumperc_36 = 0 in 4236
replace mouse_ortholog_phenolog_max_rank = 0 in 4236
replace cumperc_37 = 0 in 4236
replace mouse_ldo_phenolog_additive_rank = 0 in 4236
replace cumperc_38 = 0 in 4236
replace v91 = 0 in 4236
replace cumperc_39 = 0 in 4236
replace top_zebrafish_phenolog_max_rank = 0 in 4236
replace cumperc_40 = 0 in 4236
replace v95 = 0 in 4236
replace cumperc_41 = 0 in 4236
replace top_mouse_phenolog_max_rank = 0 in 4236
replace cumperc_42 = 0 in 4236
replace top_mouse_phenolog_additive_rank = 0 in 4236
replace cumperc_51 = 0 in 4236
replace top_ldo_phenolog_max_rank = 0 in 4236
replace cumperc_52 = 0 in 4236
replace top_ldo_phenolog_additive_rank = 0 in 4236
replace cumperc_53 = 0 in 4236
replace top_ortholog_phenolog_max_rank = 0 in 4236
replace cumperc_54 = 0 in 4236
replace top_ortholog_phenolog_additive_r = 0 in 4236




/* Important Paper Graphs */


/* Combined OWLSim/Phenolog graphs */
twoway (line cumperc_1 top_owlsim_max_ic_rank, sort) (line cumperc_2 top_owlsim_iccs_rank, sort) (line cumperc_3 top_owlsim_sim_ic_rank, sort) (line cumperc_4 top_owlsim_sim_j_rank, sort) (line cumperc_6 top_phenolog_max_rank, sort) (line cumperc_5 top_phenolog_additive_rank, sort), name(Combo_OWLSim_Phenolog) title("Recall-Rank of Human Disease Genes" "by Prediction Method") ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ") label(5 "Phenolog - Max") label(6 "Phenolog - Additive")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined OWLSim/Phenolog graphs - Mouse Only */
twoway (line cumperc_47 top_mouse_max_ic_rank, sort) (line cumperc_48 top_mouse_iccs_rank, sort) (line cumperc_49 top_mouse_sim_ic_rank, sort) (line cumperc_50 top_mouse_sim_j_rank, sort) (line cumperc_41 top_mouse_phenolog_max_rank, sort) (line cumperc_42 top_mouse_phenolog_additive_rank, sort), name(Combo_Mouse_OWLSim_Phenolog) title("Recall-Rank of Human Disease Genes" "by Prediction Method from Mouse") ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ") label(5 "Phenolog - Max") label(6 "Phenolog - Additive")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined OWLSim/Phenolog graphs - Zebrafish Only */
twoway (line cumperc_43 top_zebrafish_max_ic_rank, sort) (line cumperc_44 top_zebrafish_iccs_rank, sort) (line cumperc_45 top_zebrafish_sim_ic_rank, sort) (line cumperc_46 top_zebrafish_sim_j_rank, sort) (line cumperc_39 top_zebrafish_phenolog_max_rank, sort) (line cumperc_40 v95, sort), name(Combo_Zebrafish_OWLSim_Phenolog) title("Recall-Rank of Human Disease Genes" "by Prediction Method from Zebrafish") ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ") label(5 "Phenolog - Max") label(6 "Phenolog - Additive")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined OWLSim/Phenolog graphs - LDO Only */
twoway (line cumperc_23 top_ldo_max_ic_rank, sort) (line cumperc_24 top_ldo_iccs_rank, sort) (line cumperc_25 top_ldo_sim_ic_rank, sort) (line cumperc_26 top_ldo_sim_j_rank, sort) (line cumperc_51 top_ldo_phenolog_max_rank, sort) (line cumperc_52 top_ldo_phenolog_additive_rank, sort), name(Combo_Top_LDO_OWLSim_Phenolog) title("Recall-Rank of Human Disease Genes" "by Prediction Method from LDOs") ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ") label(5 "Phenolog - Max") label(6 "Phenolog - Additive")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined OWLSim/Phenolog graphs - Ortholog Only */
twoway (line cumperc_27 top_ortholog_max_ic_rank, sort) (line cumperc_28 top_ortholog_iccs_rank, sort) (line cumperc_29 top_ortholog_sim_ic_rank, sort) (line cumperc_30 top_ortholog_sim_j_rank, sort) (line cumperc_53 top_ortholog_phenolog_max_rank, sort) (line cumperc_54 top_ortholog_phenolog_additive_r, sort), name(combo_top_ortho_OWLSim_Phenolog) title("Recall-Rank of Human Disease Genes" "by Prediction Method from Orthologs") ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ") label(5 "Phenolog - Max") label(6 "Phenolog - Additive")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)











/* Combined top OWLSim graphs */
twoway (line cumperc_1 top_owlsim_max_ic_rank, sort) (line cumperc_2 top_owlsim_iccs_rank, sort) (line cumperc_3 top_owlsim_sim_ic_rank, sort) (line cumperc_4 top_owlsim_sim_j_rank, sort), name(Top_OWLSim_Graphs) title(Recall-Rank of Human Disease Genes by OWLSim Metric) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined zebrafish LDO OWLSim graphs */
twoway (line cumperc_7 zebrafish_ldo_max_ic_rank, sort) (line cumperc_8 zebrafish_ldo_iccs_rank, sort) (line cumperc_9 zebrafish_ldo_sim_ic_rank, sort) (line cumperc_10 zebrafish_ldo_sim_j_rank, sort), name(Combo_zebrafish_LDO_OWLSim) title(Recall-Rank of Human Disease Genes from Zebrafish Least Diverged Orthologs by OWLSim Metrics) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "zebrafish_ldo_max_ic_rank") label(2 "zebrafish_ldo_iccs_rank") label(3 "zebrafish_ldo_sim_ic_rank") label(4 "zebrafish_ldo_sim_j_rank"))yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined mouse LDO OWLSim graphs */
twoway (line cumperc_11 mouse_ldo_max_ic_rank, sort) (line cumperc_12 mouse_ldo_iccs_rank, sort) (line cumperc_13 mouse_ldo_sim_ic_rank, sort) (line cumperc_14 mouse_ldo_sim_j_rank, sort), name(Combo_mouse_LDO_OWLSim) title(Recall-Rank of Human Disease Genes from Mouse Least Diverged Orthologs by OWLSim Metrics) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "mouse_ldo_max_ic_rank") label(2 "mouse_ldo_iccs_rank") label(3 "mouse_ldo_sim_ic_rank") label(4 "mouse_ldo_sim_j_rank"))yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined top Phenolog graphs */
twoway (line cumperc_6 top_phenolog_max_rank, sort) (line cumperc_5 top_phenolog_additive_rank, sort), name(Combo_top_Phenolog) title(Recall-Rank of Human Disease Genes from Zebrafish Least Diverged Orthologs by OWLSim Metrics) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_phenolog_max_rank") label(2 "top_phenolog_additive_rank")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined OWLSim/Phenolog graphs */
twoway (line cumperc_1 top_owlsim_max_ic_rank, sort) (line cumperc_2 top_owlsim_iccs_rank, sort) (line cumperc_3 top_owlsim_sim_ic_rank, sort) (line cumperc_4 top_owlsim_sim_j_rank, sort) (line cumperc_6 top_phenolog_max_rank, sort) (line cumperc_5 top_phenolog_additive_rank, sort), name(Combo_OWLSim_Phenolog) title(Recall-Rank of Human Disease Genes by Prediction Method) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "OWLSim - MaxIC") label(2 "OWLSim - ICCS") label(3 "OWLSim - SimIC") label(4 "OWLSim - SimJ") label(5 "Phenolog - Max") label(6 "Phenolog - Additive")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined zebrafish ortholog OWLSim graphs */
twoway (line cumperc_15 zebrafish_ortholog_max_ic_rank, sort) (line cumperc_16 zebrafish_ortholog_iccs_rank, sort) (line cumperc_17 zebrafish_ortholog_sim_ic_rank, sort) (line cumperc_18 zebrafish_ortholog_sim_j_rank, sort), name(Combo_zebrafish_ortho_OWLSim) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "zebrafish_ortholog_max_ic_rank") label(2 "zebrafish_ortholog_iccs_rank") label(3 "zebrafish_ortholog_sim_ic_rank") label(4 "zebrafish_ortholog_sim_j_rank")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined mouse ortholog OWLSim graphs */
twoway (line cumperc_19 mouse_ortholog_max_ic_rank, sort) (line cumperc_20 mouse_ortholog_iccs_rank, sort) (line cumperc_21 mouse_ortholog_sim_ic_rank, sort) (line cumperc_22 mouse_ortholog_sim_j_rank, sort), name(Combo_mouse_ortho_OWLSim) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "mouse_ortholog_max_ic_rank") label(2 "mouse_ortholog_iccs_rank") label(3 "mouse_ortholog_sim_ic_rank") label(4 "mouse_ortholog_sim_j_rank")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined top LDO OWLSim graphs */
twoway (line cumperc_23 top_ldo_max_ic_rank, sort) (line cumperc_24 top_ldo_iccs_rank, sort) (line cumperc_25 top_ldo_sim_ic_rank, sort) (line cumperc_26 top_ldo_sim_j_rank, sort), name(Combo_top_LDO_OWLSim) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_ldo_max_ic_rank") label(2 "top_ldo_iccs_rank") label(3 "top_ldo_sim_ic_rank") label(4 "top_ldo_sim_j_rank")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

/* Combined top ortholog OWLSim graphs */
twoway (line cumperc_27 top_ortholog_max_ic_rank, sort) (line cumperc_28 top_ortholog_iccs_rank, sort) (line cumperc_29 top_ortholog_sim_ic_rank, sort) (line cumperc_30 top_ortholog_sim_j_rank, sort), name(combo_top_ortho_OWLSim) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) legend(label(1 "top_ldo_max_ic_rank") label(2 "top_ldo_iccs_rank") label(3 "top_ldo_sim_ic_rank") label(4 "top_ldo_sim_j_rank")) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)





/* Individual graphs */
twoway (line cumperc_1 top_owlsim_max_ic_rank, sort), name(top_owlsim_max_ic_rank) title(Top OWLSim MaxIC Rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_2 top_owlsim_iccs_rank, sort), name(top_owlsim_iccs_rank) title(Top OWLSim ICCS Rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_3 top_owlsim_sim_ic_rank, sort), name(top_owlsim_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_4 top_owlsim_sim_j_rank, sort), name(top_owlsim_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_5 top_phenolog_additive_rank, sort), name(top_phenolog_additive_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_6 top_phenolog_max_rank, sort), name(top_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_7 zebrafish_ldo_max_ic_rank, sort), name(zebrafish_ldo_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_8 zebrafish_ldo_iccs_rank, sort), name(zebrafish_ldo_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_9 zebrafish_ldo_sim_ic_rank, sort), name(zebrafish_ldo_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_10 zebrafish_ldo_sim_j_rank, sort), name(zebrafish_ldo_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_11 mouse_ldo_max_ic_rank, sort), name(mouse_ldo_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_12 mouse_ldo_iccs_rank, sort), name(mouse_ldo_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_13 mouse_ldo_sim_ic_rank, sort), name(mouse_ldo_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_14 mouse_ldo_sim_j_rank, sort), name(mouse_ldo_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_15 zebrafish_ortholog_max_ic_rank, sort), name(zebrafish_ortholog_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_16 zebrafish_ortholog_iccs_rank, sort), name(zebrafish_ortholog_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_17 zebrafish_ortholog_sim_ic_rank, sort), name(zebrafish_ortholog_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_18 zebrafish_ortholog_sim_j_rank, sort), name(zebrafish_ortholog_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_19 mouse_ortholog_max_ic_rank, sort), name(mouse_ortholog_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_20 mouse_ortholog_iccs_rank, sort), name(mouse_ortholog_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_21 mouse_ortholog_sim_ic_rank, sort), name(mouse_ortholog_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_22 mouse_ortholog_sim_j_rank, sort), name(mouse_ortholog_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_23 top_ldo_max_ic_rank, sort), name(top_ldo_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_24 top_ldo_iccs_rank, sort), name(top_ldo_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_25 top_ldo_sim_ic_rank, sort), name(top_ldo_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_26 top_ldo_sim_j_rank, sort), name(top_ldo_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_27 top_ortholog_max_ic_rank, sort), name(top_ortholog_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_28 top_ortholog_iccs_rank, sort), name(top_ortholog_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_29 top_ortholog_sim_ic_rank, sort), name(top_ortholog_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_30 top_ortholog_sim_j_rank, sort), name(top_ortholog_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_43 top_zebrafish_max_ic_rank, sort), name(top_zebrafish_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_44 top_zebrafish_iccs_rank, sort), name(top_zebrafish_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_45 top_zebrafish_sim_ic_rank, sort), name(top_zebrafish_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_46 top_zebrafish_sim_j_rank, sort), name(top_zebrafish_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_47 top_mouse_max_ic_rank, sort), name(top_mouse_max_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_48 top_mouse_iccs_rank, sort), name(top_mouse_iccs_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_49 top_mouse_sim_ic_rank, sort), name(top_mouse_sim_ic_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_50 top_mouse_sim_j_rank, sort), name(top_mouse_sim_j_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_31 zebrafish_ldo_phenolog_max_rank, sort), name(zfish_ldo_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_32 v79, sort), name(zfish_ortho_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_33 v85, sort), name(zfish_ldo_phenolog_add_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_34 v87, sort), name(zfish_ortho_phenolog_add_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_35 mouse_ldo_phenolog_max_rank, sort), name(mouse_ldo_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_36 mouse_ortholog_phenolog_max_rank, sort), name(mouse_ortholog_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_37 mouse_ldo_phenolog_additive_rank, sort), name(mouse_ldo_phenolog_additive_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_38 v91, sort), name(mouse_ldo_ortholog_additive_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_39 top_zebrafish_phenolog_max_rank, sort), name(top_zebrafish_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_40 v95, sort), name(top_zfish_phenolog_add_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)

twoway (line cumperc_41 top_mouse_phenolog_max_rank, sort), name(top_mouse_phenolog_max_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)
twoway (line cumperc_42 top_mouse_phenolog_additive_rank, sort), name(top_mouse_phenolog_add_rank) ytitle(Percentage recall of human orthologs) xtitle(In top n hits) yscale(range(0 1)) xscale(range(0 500)) ylabel(#6) xlabel(#10)





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
