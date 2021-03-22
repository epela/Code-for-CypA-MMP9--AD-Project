use "master_dataset.dta", clear

*removing exclusions from UK Biobank - see https://data.bris.ac.uk/data/dataset/1ovaau5sxunp2cv8rcy88688v
merge 1:1 id_ieu using "exclusions_highly_related.dta", keep(1) nogen
merge 1:1 id_ieu using "exclusions_non_white_british.dta", keep(1) nogen
merge 1:1 id_ieu using "exclusions_recommended.dta", keep(1) nogen
merge 1:1 id_ieu using "exclusions_relateds.dta", keep(1) nogen

**removing withdrawals (from latest file provided by UKB)
merge 1:1 id_phe using "exclusions_withdrawals.dta"
drop if withdraw==1

*merging in the genetic risk score for cypa and MMP9 eQTLs and pQTLs -scale is SDs of pQTLs
drop if id_ieu==""
*eQTLs
merge 1:1 id_ieu using "cypa_mmp9_grs.dta", keep(3) nogen
drop _m
rename ensg00000196262 cypa_eqtls
rename ensg00000100985 mmp9_eqtls
*pQTLs
drop if id_ieu==""
merge 1:1 id_ieu using "pqtl_grs.dta", keep(3) nogen
drop _m
rename cypa cypa_pqtls
rename mmp9 mmp9_pqtls

*merging in the APOE4 SNP
drop if id_ieu==""
merge 1:1 id_ieu using "apoe_snp.dta", keep(3) nogen
gen APOE4=. 
replace APOE4=0 if apoe==0
replace APOE=1 if apoe>0 & apoe <0.3
replace APOE=2 if apoe>0.3


*generating AD variables
gen mother_AD=0
replace mother_AD=1 if n_20110_0_0==10
replace mother_AD=1 if n_20110_0_0==10
replace mother_AD=1 if n_20110_0_1 ==10
replace mother_AD=1 if n_20110_0_2 ==10
replace mother_AD=1 if n_20110_1_0 ==10
replace mother_AD=1 if n_20110_1_1 ==10
replace mother_AD=1 if n_20110_1_2 ==10
replace mother_AD=1 if n_20110_2_0 ==10
replace mother_AD=1 if n_20110_2_1 ==10
replace mother_AD=1 if n_20110_2_2 ==10
tab mother_AD
replace mother_AD=. if n_20110_0_0==. &  n_20110_0_1 ==. & n_20110_0_2==. &  n_20110_1_0==. & n_20110_1_1==. & n_20110_1_2==. & n_20110_2_0==. & n_20110_2_1==. & n_20110_2_2==.

gen father_AD=0
replace father_AD=1 if n_20107_0_0==10
replace father_AD=1 if n_20107_0_1==10
replace father_AD=1 if n_20107_0_2==10
replace father_AD=1 if n_20107_0_3==10
replace father_AD=1 if n_20107_1_0==10
replace father_AD=1 if n_20107_1_1==10
replace father_AD=1 if n_20107_1_2==10
replace father_AD=1 if n_20107_2_0==10
replace father_AD=1 if n_20107_2_1==10
replace father_AD=1 if n_20107_2_2==10
replace father_AD=1 if n_20107_2_3==10
tab father_AD

replace father_AD=. if  n_20107_0_0==. &  n_20107_0_1==. &  n_20107_0_2==. &  n_20107_0_3==. &  n_20107_1_0==. &  n_20107_1_1==. &  n_20107_1_2==. &  n_20107_2_0==. &  n_20107_2_1==. &  n_20107_2_2==. &  n_20107_2_3==. 

*generating marker of 0, 1, or 2 parents with AD
gen parent_ad=0 if mother_AD==0 & father_AD ==0
replace parent_ad=1 if mother_AD ==1 | father_AD==1
replace parent_ad=2 if mother_AD ==1 & father_AD ==1



******************************
********Cleaning data*********
******************************
*fluid intelligence*
summ n_20016*
replace n_20016_0_0= n_20016_1_0 if n_20016_0_0==.
replace n_20016_0_0= n_20016_1_0 if n_20016_0_0==.
summ n_20016_0_0
*standing FI scores
egen zFI_score = std(n_20016_0_0)
summ zFI_score 

*reaction time*
summ n_20023_0_0 n_20023_1_0 n_20023_2_0
replace n_20023_0_0= n_20023_1_0 if n_20023_0_0==.
replace n_20023_0_0= n_20023_2_0 if n_20023_0_0==.
summ n_20023_0_0
*standardizing reaction time scores
egen zreaction_time = std(n_20023_0_0)
summ zreaction_time

*visual memory*
replace n_399_0_1 = n_399_0_2 if n_399_0_1==.
replace n_399_0_1 = n_399_1_1 if n_399_0_1==.
replace n_399_0_1 = n_399_1_2 if n_399_0_1==.
replace n_399_0_1 = n_399_2_1 if n_399_0_1==.
replace n_399_0_1 = n_399_2_2 if n_399_0_1==.
replace n_399_0_1 = n_399_2_3 if n_399_0_1==.
*standardizing visual memory scores
summ n_399_0_1
egen zvisual_memory = std(n_399_0_1)
summ zvisual_memory

*generating tertiles for age stratified analyis
xtile age_tert = n_21022_0_0, nq(3)

*keeping eligible sample (participants with genetic data plus at least one outcome of interest - note no missing for age and sex)
keep if cypa!=. & mmp9!=.
keep if zFI_score!=. | zreaction_time!=. | zvisual_memory!=. | parent_ad!=.

save "analysis_dataset.dta", replace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

*************************************************************************************************************************************
**********Main analysis in whole sample for eQTLs on AD by proxy, FI, reaction time and visual memory, all adj for age and sex*******
*************************************************************************************************************************************

**analysis for AD-by-proxy as outcome (done separately to other othercomes as logistc not linear regression, and easier for storing results)
ologit parent_ad apoe n_21022_0_0 n_22001_0_0, or
ologit parent_ad cypa n_21022_0_0 n_22001_0_0, or
ologit parent_ad mmp9 n_21022_0_0 n_22001_0_0, or

*Effect of ApoE on outcomes
use "analysis_dataset.dta", clear

regress zFI_score apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0",detail(all) pval replace
regress zreaction_time apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_1",detail(all) pval replace
regress zvisual_memory apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_main_results", replace

*Effect of cypa eqtls on outcomes
use "analysis_dataset.dta", clear

regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_main_results", replace

*Effect of mmp9 eQTLs on outcomes
use "analysis_dataset.dta", clear

regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace
regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace
regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace


use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_main_results", replace




*****************************************************************************************
******Age stratified analysis for eQTLs on FI, reaction time and visual memory***********
*****************************************************************************************

//////////////////////////////////////////TERTILE ONE///////////////////////////////////////////////////////////////////
use "analysis_dataset.dta", clear
keep if age_tert==1

*ApoE analyis
regress zFI_score apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0",detail(all) pval replace
regress zreaction_time apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_1",detail(all) pval replace
regress zvisual_memory apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_youngest_tert", replace


*cypa analysis (adj age and sex)
use "analysis_dataset.dta", clear
keep if age_tert==1

regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_youngest_tert", replace


*mmp9 analysis
use "analysis_dataset.dta", clear
keep if age_tert==1

regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace
regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace
regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace


use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_youngest_tert", replace


//////////////////////////////////////////////////////////////////////////TERTILE TWO///////////////////////////////////////////////////////////////////////////////////////////////////
use "analysis_dataset.dta", clear
keep if age_tert==2

*ApoE analyis
regress zFI_score apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0",detail(all) pval replace
regress zreaction_time apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_1",detail(all) pval replace
regress zvisual_memory apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_middle_tert", replace


*cypa analysis (adj age and sex)
use "analysis_dataset.dta", clear
keep if age_tert==2

regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_middle_tert", replace

*mmp9 analysis
use "analysis_dataset.dta", clear
keep if age_tert==2

regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace
regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace
regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace


use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_middle_tert", replace


////////////////////////////////////////////////////////////////////////////TERTILE THREE//////////////////////////////////////////////////////////////////////////////////////////////////
use "analysis_dataset.dta", clear
keep if age_tert==3

*ApoE analyis
regress zFI_score apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0",detail(all) pval replace
regress zreaction_time apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_1",detail(all) pval replace
regress zvisual_memory apoe n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\apoe_oldest_tert", replace


*cypa analysis (adj age and sex)
use "analysis_dataset.dta", clear
keep if age_tert==3

regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_oldest_tert", replace


*mmp9 analysis
use "analysis_dataset.dta", clear
keep if age_tert==3

regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace
regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace
regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace


use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_oldest_tert", replace


************************************************************************************************************************************
***EXAMINING EFFECT OF CYPA AND MMP9 EQTLS ON OUTCOMES, STRATIFIED BY APOE4 CARRIER STATUS, ADJ AGE AND SEX*************************
************************************************************************************************************************************

///////////////////////////////////////////////ZERO APOE4 ALLELES/////////////////////////////////////////////
use "analysis_dataset.dta", clear
keep if APOE4==0

*cypa
regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_no_apoe", replace

use "analysis_dataset.dta", clear
keep if APOE4==0


*mmp9
regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace

regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace

regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_no_apoe", replace


////////////////////////////////////////////////////////////////////ONE APOE4 ALLELE//////////////////////////////////////////////////////////////////////////////////////////////////////////

use "analysis_dataset.dta", clear
keep if APOE4==1

*cypa
regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_hetero_apoe", replace


*mmp9

use "analysis_dataset.dta", clear
keep if APOE4==1

regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace

regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace

regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_hetero_apoe", replace


///////////////////////////////////////////////////////////////////TWO APOE4 ALLELES///////////////////////////////////////////////////////////////////////////////////////////////////////////

use "analysis_dataset.dta", clear
keep if APOE4==2

*cypa
regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace
regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace
regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_homo_apoe", replace


use "analysis_dataset.dta", clear
keep if APOE4==2


*mmp9
regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace

regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace

regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_homo_apoe", replace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

*outputting interaction P values
use "analysis_dataset.dta", clear

gen apoe_cypa=apoe*cypa
gen apoe_mmp9=apoe*mmp9

*cypa
regress zFI_score apoe_cypa apoe cypa n_21022_0_0 n_22001_0_0
regress zreaction_time apoe_cypa apoe cypa n_21022_0_0 n_22001_0_0
regress zvisual_memory apoe_cypa apoe cypa n_21022_0_0 n_22001_0_0

*mmp9
regress zFI_score apoe_mmp9 apoe mmp9 n_21022_0_0 n_22001_0_0
regress zreaction_time apoe_mmp9 apoe mmp9 n_21022_0_0 n_22001_0_0
regress zvisual_memory apoe_mmp9 apoe mmp9 n_21022_0_0 n_22001_0_0


*******************************************************************
*****Interactions between APOE and CYPA/MMP9 in oldest tertile*****
*******************************************************************
//////////////////////////////////////////////////////////////////////////////ZERO E4 ALLELES\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
*APOE stratified analysis
use "analysis_dataset.dta", clear
keep if age_tert==3
keep if APOE4==0

*cypa
regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace

regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace

regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_no_apoe_oldest_tert", replace

use "analysis_dataset.dta", clear
keep if age_tert==3
keep if APOE4==0


*mmp9
regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace

regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace

regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_no_apoe_oldest_tert", replace



//////////////////////////////////////////////////////////////////////////////ONE E4 ALLELE////////////////////////////////////////////////////////////////////////////////////////////////////

use "analysis_dataset.dta", clear
keep if age_tert==3
keep if APOE4==1

*cypa
regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace

regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace

regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_hetero_apoe_oldest_tert", replace


use "analysis_dataset.dta", clear
keep if age_tert==3
keep if APOE4==1


*mmp9
regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace

regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace

regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_hetero_apoe_oldest_tert", replace


/////////////////////////////////////////////////////////////////////////////TWO E4 ALLELES ////////////////////////////////////////////////////////////////////////////////////////////////

use "analysis_dataset.dta", clear
keep if age_tert==3
keep if APOE4==2

*cypa
regress zFI_score cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0",detail(all) pval replace

regress zreaction_time cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_1",detail(all) pval replace

regress zvisual_memory cypa n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\cypa_homo_apoe_oldest_tert", replace



use "analysis_dataset.dta", clear
keep if age_tert==3
keep if APOE4==2


*mmp9
regress zFI_score mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0",detail(all) pval replace

regress zreaction_time mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_1",detail(all) pval replace

regress zvisual_memory mmp9 n_21022_0_0 n_22001_0_0
regsave using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_2",detail(all) pval replace

use "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_0", clear

gen exposure="exp0"

forvalues i=1(1)2{
append using "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_`i'" 
replace exposure="exp`i'" if exposure==""
}

order exposure var
sort exposure 
bysort exposure: keep if _n==1

gen p25=(coef-stderr*1.96)
gen p75=(coef+stderr*1.96)

gen lci=p25
gen uci=p75

gen beta=coef
gen CI= string(beta ,"%9.4f") + " " + "(" +string( lci ,"%9.2f")+ "," +string(uci ,"%9.2f")+")" 
format pval %9.2f
order exposure var beta CI pval

save "C:\Users\epela\OneDrive - University of Bristol\MIGRATED O DRIVE\MRC Skills Development Fellowship\Projects\CYPA_MMP9\Data\pQTL analysis\regsave results\mmp9_homo_apoe_oldest_tert", replace

