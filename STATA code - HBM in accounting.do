**Hierarchical Bayesian models in accounting: A tutorial 
**Bullipe R. Chintha and Sanjay Kallapur
**This do-file takes Compustat and risk-free rates data as input and creates Frequentist and Bayesian estimates of ERCs.

*Set working directory and create a "Temp" folder within Data folder for saving and using temp files
cd "Data/"

***Pre-processing (Risk-free rates for HBM3)

quietly{

**Import data of Long-Term Government Bond Yields: 10-year (from https://fred.stlouisfed.org/series/IRLTLT01USM156N)
import delimited "RF rates.csv", encoding(ISO-8859-2) clear 
gen date1=date(date,"YMD")
gen month=month(date1)
gen year=year(date1)
*Use December's rate
keep if month==12
rename year fyear
gen rf = irltlt01usm156n/100
keep fyear rf
save "Temp\RF rate.dta", replace

}

***Prepare data

quietly{

*Download Compustat data in .dta format[Standard filters: indfmt=="INDL" & datafmt=="STD" & popsrc=="D" & consol=="C"] 
use "Compustat.dta", replace

**Use fiscal year data [1980-2020]
keep if fyear>=1980
keep if fyear<=2020

**Filters
*Sanity Filters
keep if ceq>=0.5 
keep if at>=1.5
keep if prcc_f>=3
*Use December fiscalend observations
keep if fyr==12

**Fama French-12 classification
ffind sic, newvar(ffind_no) type(12)
*Drop financial and utility firms
drop if ffind_no==8 | ffind_no==11 

**Variable estimation
xtset gvkey fyear
gen mktcap=csho*prcc_f
gen Adj_P = prcc_f/ajex
gen Adj_D = dvpsp_f /ajex
gen return = (Adj_P-L1.Adj_P+Adj_D )/L1.Adj_P
gen earn = ib/L1.mktcap
gen delta_E = (ib-L1.ib)/L1.mktcap
gen lag_earn = L.earn

**Drop missing values
keep gvkey fyear return earn delta_E lag_earn
foreach v of var * { 
drop if missing(`v') 
}

**Check for gaps in time-series
egen maxfyear=max(fyear), by(gvkey)
egen minfyear=min(fyear), by(gvkey)
gen sample=maxfyear-minfyear+1
egen count=count(fyear), by (gvkey)
drop if sample!=count
drop maxfyear minfyear sample

**Restrict the sample to firms with more than 30 years data 
drop if count<30
distinct gvkey
*301 firms

**Split the sample into pre and post
by gvkey: gen year_no = _n
gen test_sample=1 if year_no <=15
replace test_sample =0 if test_sample ==.
tab test_sample

**Estimation of OLS ERCs
statsby _b _se, by(gvkey) saving(Temp/my_estimates_ERCdx, replace) : regress return delta_E if test_sample ==1
statsby _b _se, by(gvkey) saving(Temp/my_estimates_ERCdx_post, replace) : regress return delta_E if test_sample ==0

**Merging OLS coefficients 
merge m:1 gvkey using "Temp\my_estimates_ERCdx.dta"
drop _merge
rename _b_delta_E ERC_ols
rename _se_delta_E ERC_ols_sd
merge m:1 gvkey using "Temp\my_estimates_ERCdx_post.dta"
drop _merge
rename _b_delta_E ERC_ols_post
rename _se_delta_E ERC_ols_sd_post

**Merge risk-free rates
merge m:1 fyear using "Temp\RF rate.dta"
keep if _merge==3
drop _merge

**Define other variables for model-building
keep if test_sample==1
egen avg_RFrate = mean(rf), by(gvkey)
gen df_resid =13
gen ERC_var = ERC_ols_sd^2
sort gvkey fyear

**Save processed data files for model-building
export delimited "Temp/Input_for_HBM2_HBM3.csv", replace
duplicates drop gvkey, force
export delimited "Temp/Input_for_HBM1.csv", replace

}

***HBM 1: Input OLS ERCs

import delimited "Temp/Input_for_HBM1.csv", clear case(preserve)

*Input ERC_ols
bayesmh ERC_ols, noconstant reffects(gvkey) likelihood(t(ERC_var,13)) prior({ERC_ols:i.gvkey}, normal({mu},{sig2})) prior({mu}, normal(0, 1e6)) prior({sig2}, igamma(0.001, 0.001)) block({mu}) block({sig2}) adapt(tolerance(0.002)) showreffects({ERC_ols:i.gvkey}) nchains(3) initall({ERC_ols:} rnormal( 0, 100) {mu} rnormal(0, 100) {sig2} runiform(0, 100)) mcmcsize(5000) burnin(5000) rseed(149)

*For saving simulated data
bayesmh, saving(Temp/simdata_model1.dta, replace)

*MCMC convergence
bayesstats grubin

*MCMC efficiency
bayesstats ess

*Diagnostic plot - firmlevel ERC
bayesgraph diagnostic {ERC_ols: 1078.gvkey}

*Trace plot
bayesgraph trace {ERC_ols: 1078.gvkey}

*Density plot
bayesgraph kdensity {ERC_ols: 1078.gvkey}

*Diagnostic plot - hyperparameter mu
bayesgraph diagnostic {mu}


***HBM2: Varying-intercepts and varying-slopes model: inverse-Wishart prior 

*Input data and suppress the base level of “gvkey” in order to use factor notation
import delimited "Temp/Input_for_HBM2_HBM3.csv", clear case(preserve)
fvset base none gvkey

*Input raw data of return and delta_E
bayesmh return i.gvkey i.gvkey#c.delta_E, noconstant likelihood(normal({var_0})) prior ({return:i.gvkey i.gvkey#c.delta_E}, mvnormal(2, {return:_cons}, {return:delta_E}, {covar,m})) prior({var_0}, igamma(0.01, 0.01)) prior({covar,m}, iwishart(2, 3, I(2))) prior({return:delta_E}, normal(0, 100)) prior({return:_cons}, normal(0, 100)) block ({return: i.gvkey}, reffects) block ({return: i.gvkey#c.delta_E}, reffects) block({var_0}, gibbs) block({covar,m}, gibbs) block({return:_cons}) block({return:delta_E}) adapt(tolerance(0.002)) nchains(3) initall({return: i.gvkey} rnormal( 0, 100) {return: i.gvkey#c.delta_E} rnormal(0, 100) {var_0} runiform(0, 100) {return:delta_E} rnormal(0, 100) {return:_cons} rnormal(0, 100) {covar,m} I(2)) mcmcsize(5000) burnin(5000) rseed(149)

*Save the simulations
bayesmh, saving(Temp/simdata_model2.dta, replace)


***HBM3: Varying-intercepts and varying-slopes model: using information of risk-free rates 

*Input data and suppress the base level of “gvkey” in order to use factor notation
import delimited "Temp/Input_for_HBM2_HBM3.csv", clear case(preserve)
gen ratedelta_E = rf*delta_E
fvset base none gvkey

*Input raw data of return, delta_E and ratedelta_E
bayesmh return i.gvkey i.gvkey#c.delta_E ratedelta_E, noconstant likelihood(normal({var_0}))  prior ({return:i.gvkey i.gvkey#c.delta_E}, mvnormal(2, {return:_cons}, {return:delta_E}, {covar,m})) prior({var_0}, igamma(0.01, 0.01)) prior({covar,m}, iwishart(2, 3, I(2))) prior({return: delta_E}, normal(0, 100)) prior({return: _cons}, normal(0, 100)) prior({return: ratedelta_E}, normal(0, 100)) block ({return: i.gvkey}, reffects) block ({return: i.gvkey#c.delta_E}, reffects) block({var_0}, gibbs) block({covar,m}, gibbs) block({return:_cons}) block({return:delta_E}) block({return:ratedelta_E}) adapt(tolerance(0.002)) nchains(3) initall({return: i.gvkey} rnormal( 0, 100) {return: i.gvkey#c.delta_E} rnormal(0, 100) {var_0} runiform(0, 100) {return: delta_E} rnormal(0, 100) {return: _cons} rnormal(0, 100) {return: ratedelta_E} rnormal(0, 100) {covar,m} I(2)) mcmcsize(5000) burnin(5000) rseed(149)

*Save the simulations
bayesmh, saving(Temp/simdata_model3.dta, replace)


***Using simulated data: HBM1-3

quietly{

use "Temp/simdata_model1.dta", replace
collapse (mean) *
drop _chain _index _loglikelihood _logposterior _frequency
save "Temp/model1estimates.dta", replace

use "Temp/simdata_model2.dta", replace
collapse (mean) *
foreach var of varlist _chain - eq1_p301{
drop `var'
}
drop eq0_p1 eq0_p2 eq0_p3 eq0_p4 _frequency
save "Temp/model2estimates.dta", replace

use "Temp/simdata_model3.dta", replace
collapse (mean) *
foreach var of varlist _chain - eq1_p301{
drop `var'
}
drop eq0_p1 eq0_p2 eq0_p3 eq0_p4 _frequency
save "Temp/model3estimates.dta", replace

}

***


