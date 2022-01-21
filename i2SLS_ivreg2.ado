* 16/12 : change constant calculation to avoid a log of 0 & change eps.
* 19/12 change covariance matrix calculation for large data set
* 19/12 : add correction when no covariate is included.
* 21/12 : Manual iteration of 2SLS GMM + options to control nb iterations /convergence..
* 04/01 : retour de la constante + check de convergence 
* 21/01 : symmetric S.E. + correction de syntax + check singleton de PPML

cap program drop i2SLS_ivreg2
program define i2SLS_ivreg2, eclass
//syntax anything(fv ts numeric) [if] [in] [aweight pweight fweight iweight]  [, DELta(real 1) LIMit(real 0.00001) MAXimum(real 1000) Robust CLuster(string)  ]

syntax varlist [if] [in] [aweight pweight fweight iweight] [, DELta(real 1) LIMit(real 0.00001) ENDog(string) INSTR(string) MAXimum(real 1000) Robust CLuster(string)]           

marksample touse   
markout `touse'  `cluster', s  
	preserve 
quietly keep if `touse'
	
	if "`gmm2s'" !="" {
		local opt0 = "`gmm2s' "
	}
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "cluster(`cluster') "
	}
	local option = "`opt0'`opt1'`opt2'"
	*** Obtain lists of variables 
	local list_var `varlist'
	gettoken depvar list_var : list_var
	gettoken _rhs list_var : list_var, p("(")
*** check seperation : code from "ppml"
 tempvar logy                            																						// Creates regressand for first step
   qui gen `logy'=.  if (`touse')                                                                                // Creates regressand for first step
 quietly: replace `logy'=log(`depvar') if (`touse')&(`depvar'>0)
 quietly: reg `logy' `_rhs' `endog' if (`touse')&(`depvar'>0)	
 tempvar zeros                            																						// Creates regressand for first step
 quietly: gen `zeros'=1                                                                                                 // Initialize observations selector
 local _drop ""                                                                                                // List of regressors to exclude
 local indepvar ""     
    foreach x of varlist `_rhs' `endog' {  
      if (_se[`x']==0) {                                                                                       // Try to include regressors dropped
          qui summarize `x' if (`depvar'>0)&(`touse'), meanonly
          local _mean=r(mean)
          qui summarize `x' if (`depvar'==0)&(`touse')
          if (r(min)<`_mean')&(r(max)>`_mean'){                                            // Include regressor if conditions met and
		  if (`x'!=`endog'){
              local indepvar "`indepvar' `x'"     
		  } 
          }
          else{
              qui su `x' if `touse', d                                                                         // Otherwise, drop regressor
              local _mad=r(p50)
              qui inspect  `x'  if `touse'                                                                         
              qui replace `zeros'=0 if (`x'!=`_mad')&(r(N_unique)==2)&(`touse')                     // Mark observations to drop
              local _drop "`_drop' `x'"
          }
      }
      if (_se[`x']>0)&(`x'!=`endog') {                                                                  // Include safe regressors: LOOP 1.2
      local indepvar "`indepvar' `x'" 
      }
// End LOOP 1.2
    }   
 qui su `touse' if `touse', mean                                                                               // Summarize touse to obtain N
 local _enne=r(sum)                                                                                            // Save N
 qui replace `touse'=0 if (`zeros'==0)&("`keep'"=="")&(`depvar'==0)&(`touse')                                  // Drop observations with perfect fit
 di                                                                                                              // if keep is off
 local k_excluded : word count `_drop'                                                                         // Number of variables causing perfect fit
 di in green "Number of regressors excluded to ensure that the estimates exist: `k_excluded'" 
 if ("`_drop'" != "") di "Excluded regressors: `_drop'"                                                        // List dropped variables if any
 qui su `touse' if `touse', mean
 local _enne = `_enne' - r(sum)                                                                                // Number of observations dropped
 di in green "Number of observations excluded: `_enne'" 
 local _enne =  r(sum)
quietly keep if `touse'	
** drop collinear variables
	tempvar cste
	gen `cste' = 1
    _rmcoll `indepvar' `cste', forcedrop 
	local var_list `endog' `r(varlist)' `cste'  
	local instr_list `instr' `r(varlist)' `cste' 
	*** Initialisation de la boucle
	tempvar y_tild 
	quietly gen `y_tild' = log(`depvar' + 1)
	** prepare 2SLS
	*local var_list  `endog' `indepvar' `cste'
	*local instr_list `instr' `indepvar' `cste'
	mata : X=.
	mata : Z=.
	mata : y_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'")
	mata : st_view(Z,.,"`instr_list'")
	mata : st_view(y_tilde,.,"`y_tild'")
	mata : st_view(y,.,"`depvar'")
	mata : invPzX = invsym(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X))*cross(X,Z)*invsym(cross(Z,Z))
	mata : beta_initial = invPzX*cross(Z,y_tilde)
	mata : beta_initial = invPzX*cross(Z,y_tilde)
	mata : beta_t_1 = beta_initial // needed to initialize
	mata : beta_t_2 = beta_initial // needed to initialize
	mata : q_hat_m0 = 0
	local k = 1
	local eps = 1000	
	mata: q_hat = J(`maximum', 1, .)
	*** Iterations iOLS
	_dots 0
	while ((`k' < `maximum') & (`eps' > `limit' )) {
		* Nouveaux beta
	mata: alpha = log(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1]) ))
	mata: beta_initial[(cols(X)),1] = alpha
	mata: xb_hat = X*beta_initial
		* Update d'un nouveau y_tild et regression avec le nouvel y_tild
	mata: y_tilde = log(y + `delta'*exp(xb_hat)) :-mean(log(y + `delta'*exp(xb_hat))- xb_hat)
		* 2SLS 
	mata: beta_new = invPzX*cross(Z,y_tilde)
		* DiffÃ©rence entre les anciens betas et les nouveaux betas
	mata: criteria = mean(abs(beta_initial - beta_new):^(2))
mata: st_numscalar("eps", criteria)
mata: st_local("eps", strofreal(criteria))
		* safeguard for convergence.
	if `k'==`maximum'{
		  di "There has been no convergence so far: increase the number of iterations."  
	}
	if `k'>4{
	mata: q_hat[`k',1] = mean(log( abs(beta_new-beta_initial):/abs(beta_initial-beta_t_2)):/log(abs(beta_initial-beta_t_2):/abs(beta_t_2-beta_t_3)))	
	mata: check_3 = abs(mean(q_hat)-1)
		if mod(`k'-4,50)==0{
    mata: q_hat_m =  mm_median(q_hat[((`k'-49)..`k'),.] ,1)
	mata: check_1 = abs(q_hat_m - q_hat_m0)
	mata: check_2 = abs(q_hat_m-1)
	mata: st_numscalar("check_1", check_1)
	mata: st_local("check_1", strofreal(check_1))
	mata: st_numscalar("check_2", check_2)
	mata: st_local("check_2", strofreal(check_2))
	mata: st_numscalar("check_3", check_3)
	mata: st_local("check_3", strofreal(check_3))
	mata: q_hat_m0 = q_hat_m
		if ((`check_1'<1e-4)&(`check_2'>1e-2)) {
di "delta is too small to achieve convergence -- update to larger value"
	local k = `maximum'
		}
		if ((`check_3'>0.5) & (`k'>500)) {
	local k = `maximum'
di "q_hat too far from 1"
		}
					  }
	}
	if `k'>2 { // keep in memory the previous beta_hat for q_hat 
	mata:   beta_t_3 = beta_t_2
	mata:   beta_t_2 = beta_initial
	}
mata: beta_initial = beta_new
	local k = `k'+1
	_dots `k' 0
	}
	*** Calcul de la bonne matrice de variance-covariance
	* Calcul du "bon" rÃ©sidu
	mata: xb_hat = X*beta_new
	mata : y_tilde = log(y + `delta'*exp(xb_hat)) :-mean(log(y + `delta'*exp(xb_hat)) - xb_hat)
	* Retour en Stata 
	cap drop y_tild 
	quietly mata: st_addvar("double", "y_tild")
	mata: st_store(.,"y_tild",y_tilde)
	quietly ivreg2 y_tild `r(varlist)' (`endog' = `instr') [`weight'`exp'] if `touse', `option' 
	cap drop xb_hat
	quietly predict xb_hat, xb
	cap drop ui
	quietly gen ui = `depvar'*exp(-xb_hat)
	quietly replace ui = ui/(`delta' + ui)
	mata : ui= st_data(.,"ui")
	* Calcul de Sigma_0, de I-W, et de Sigma_tild
	matrix beta_final = e(b) // 	mata: st_matrix("beta_final", beta_new)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X):/rows(X))*Sigma_hat*(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X):/rows(X)) // recover original HAC 
	mata : invXpPzIWX = invsym(0.5:/rows(X)*cross(X,Z)*invsym(cross(Z,Z))*cross(Z,ui,X)+ 0.5:/rows(X)*cross(X,ui,Z)*invsym(cross(Z,Z))*cross(Z,X))
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
	mata : Sigma_tild = (Sigma_tild+Sigma_tild'):/2 
    	mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
	*** Stocker les resultats dans une matrice
	local names : colnames e(b)
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
    ereturn post beta_final Sigma_tild , obs(`e(N)') depname(`depvar') esample(`touse')  dof(`=e(df r)') 
	restore 
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar widstat = e(widstat)
ereturn scalar arf = e(arf)
ereturn local cmd "i2SLS_ivreg2"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
ereturn display
end

