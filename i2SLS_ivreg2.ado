* 16/12 : change constant calculation to avoid a log of 0 & change eps.
* 19/12 change covariance matrix calculation for large data set
* 19/12 : add correction when no covariate is included.
* 21/12 : Manual iteration of 2SLS GMM + options to control nb iterations /convergence..

program define i2SLS_ivreg2, eclass
	syntax [anything] [if]  [in] [aweight pweight fweight iweight]  [, DELta(real 1) LIMit(real 0.00001) MAXimum(real 1000) Robust CLuster(varlist numeric)  ]
	marksample touse
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
	local list_var `anything'
	gettoken depvar list_var : list_var
	if (strpos("`list_var'","(")==2){  // anormal case : no X, only Z 
	local list_var: subinstr local list_var "(" "", all
	local indepvar 
	gettoken endog list_var : list_var, bind
	gettoken endog instr_temp : endog , p("=")
	local list_var: subinstr local list_var "=" "", all
	gettoken instr list_var : list_var, bind
	}
	else{  // normal case : X exists 
	gettoken indepvar list_var : list_var, p("(")
	gettoken endog list_var : list_var, bind
	gettoken endog endog : endog, p("(")
	gettoken endog instr_temp : endog , p("=")
	gettoken equalsign instr_temp : instr_temp , p("=")
	gettoken instr instr_temp : instr_temp, p(")")
	}
			
	*** Initialisation de la boucle
	tempvar y_tild 
	quietly gen `y_tild' = log(`depvar' + 1)
	tempvar cste
	gen `cste' = 1
	** drop collinear variables
    _rmcoll `indepvar' `cste', forcedrop 
	local var_list `endog' `r(varlist)' `cste'  
	local instr_list `instr' `r(varlist)' `cste' 
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
	mata : beta_initial = invPzX*Z'*y_tilde
	local k = 0
	local eps = 1000	
	*** Iterations iOLS
	_dots 0
	while ((`k' < `maximum') & (`eps' > `limit' )) {
		* Nouveaux beta
	mata: xb_hat = X*beta_initial
		* Update d'un nouveau y_tild et regression avec le nouvel y_tild
	mata: y_tilde = log(y + `delta'*exp(xb_hat)) :-mean(log(y + `delta'*exp(xb_hat))- xb_hat)
		* 2SLS 
	mata: beta_new = invPzX*Z'*y_tilde
		* DiffÃ©rence entre les anciens betas et les nouveaux betas
	mata: criteria = mean(abs(beta_initial - beta_new))
mata: st_numscalar("eps", criteria)
mata: st_local("eps", strofreal(criteria))
mata: beta_initial = beta_new
	local k = `k'+1
	_dots `k' 0
	}
	*** Calcul de la bonne matrice de variance-covariance
	* Calcul du "bon" rÃ©sidu
	mata: xb_hat = X*beta_new
	mata : y_tilde = log(y + `delta'*exp(xb_hat)) :-mean(log(y + `delta'*exp(xb_hat)) - xb_hat)
	* 	gen `ui' = exp(`y_tild' + `c_hat' - `xb_hat') - `delta'
	mata: ui = y:*exp(-xb_hat)
	mata: ui = ui:/(`delta' :+ ui)
	* Retour en Stata 
	cap drop y_tild 
	quietly mata: st_addvar("double", "y_tild")
	mata: st_store(.,"y_tild",y_tilde)
	quietly ivreg2 y_tild `r(varlist)' (`endog' = `instr') [`weight'`exp'] if `touse', `option'  gmm2s
	* Calcul de Sigma_0, de I-W, et de Sigma_tild
	matrix beta_final = e(b) // 	mata: st_matrix("beta_final", beta_new)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X))*Sigma_hat*(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X)) // recover original HAC 
	mata : invXpPzIWX = invsym(0.5*cross(X,Z)*invsym(cross(Z,Z))*cross(Z,ui,X)+ 0.5*cross(X,ui,Z)*invsym(cross(Z,Z))*cross(Z,X))
*	mata : invXpPzIWX = invsym(0.5*X'*(Z*invsym(cross(Z,Z))*(Z':*ui')+ (ui:*Z)*invsym(cross(Z,Z))*Z')*X)
	*mata : invXpPzIWX = invsym(0.5*X'*(Pz*IW+IW*Pz)*X)
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
    mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
	*** Stocker les resultats dans une matrice
	local names : colnames e(b)
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`=e(df r)') 
	restore 
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar widstat = e(widstat)
ereturn scalar arf = e(arf)
ereturn local cmd "i2SLS_ivreg2"
ereturn local vcetype `option'
ereturn display
end

