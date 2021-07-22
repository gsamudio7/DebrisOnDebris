clear all
set more off
set scheme s1mono

cd "C:\Users\timothy.justicz\OneDrive - West Point\Documents\Q\spacecom"

import delimited "DebrisOnDebris.csv", asdouble

keep v48 v49 v9 v6 v12 v14 v15

label var v48 "Identifies all the records associated with a single event"
rename v48 eventid

label var v49 "The time to TCA (Time of Closest Approach)"
rename v49 daystoTCA

label var v9 "The best estimate of the probability of collision"
destring v9, replace force
rename v9 PcBest

label var v6 "If the size estimation approach described above in the treatment of field #9 is not possible, usually because RCS data are not available, a NAN will be present in field #9; in such cases, the Pc value here is field #6 should be used (this is the Pc based on a nominal value of 1.5m for both primary and secondary, 3m total."
rename v6 PcNom

label var v12 "If the probability of a catastrophic collision is above 0.5, then I would probably assume that it will be catastrophic."
destring v12, replace force
rename v12 ProbCatIfColl

label var v14 "This field gives the number of fragments 5cm or greater in size that would be expected to be produced in this particular collision were catastrophic.  You can use this field in conjunction with field #12 to determine the number of fragments expected:  if #12 is greater than 0.5, then #14 will give the number of fragments expected."
destring v14, replace force
rename v14 NumFragIfCatColl

label var v15 "If the collision is expected to be non-catastrophic, then this field gives the number of debris fragments expected."
destring v15, replace force
rename v15 NumFragIfNonCatColl

gen catastrophic=ProbCatIfColl>=0.5
replace catastrophic=0 if missing(ProbCatIfColl)

gen fragments=NumFragIfCatColl
replace fragments=NumFragIfNonCatColl if catastrophic==0
hist fragments, frac
graph export "fragments_hist.png", replace
replace fragments=0 if missing(fragments) // is this a valid assumption? RFI to POC

gen pc=PcBest
replace pc=PcNom if missing(pc)
hist pc, frac
graph export "pc_hist.png", replace

save "debrisondebris.dta", replace

export delimited "DebrisOnDebrisCleaned.csv", replace

gen pcnonzero=pc if pc!=0
hist pcnonzero, frac
graph export "pc_nonzero_hist.png", replace

gen logpc=log(pc)
hist logpc, frac
graph export "logpc_hist.png", replace

gen logfragments=log(fragments)
hist logfragments, frac
graph export "logfragments_hist.png", replace

twoway (scatter logfragments logpc) (lfit logfragments logpc)
graph export "scatter_logfragments_logpc.png", replace
