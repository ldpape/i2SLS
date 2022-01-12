cap program drop i2SLS_ivreg2_test
program define i2SLS_ivreg2_test, eclass 
	syntax [anything] [if] [in] [aweight pweight fweight iweight] [, DELta(real 1) LIMit(real 0.00001) MAXimum(real 1000)]
	*marksample touse
	local list_var `anything'
	* get depvar and indepvar
	gettoken depvar list_var : list_var
	gettoken indepvar list_var : list_var, p("(")
    * get endogenous variables and instruments
	gettoken endog list_var : list_var, bind
	gettoken endog endog : endog, p("(")
    gettoken endog instr_temp : endog , p("=")
    gettoken equalsign instr_temp : instr_temp , p("=")
	gettoken instr instr_temp : instr_temp, p(")")
	* gen binary zero variable 
	cap drop dep_pos 
	cap drop *temp
	cap drop xb_temp
	gen dep_pos = `depvar'>0
       xi:  i2SLS_ivreg2 `depvar' `indepvar' (`endog' = `instr') , delta(`delta') maximum(`maximum') limit(`limit') robust
         *lhs of test
         predict xb_temp, xb
         gen u_hat_temp = `depvar'*exp(-xb_temp)
         gen lhs_temp = log(`delta'+u_hat_temp) - log(`delta') 
         * rhs of test
         gen temp = log(`depvar' + `delta'*exp(xb_temp)) - xb_temp
		 cap drop xb_temp
         egen c_hat_temp = mean(temp)
       xi: logit dep_pos `indepvar' `instr'
         predict p_hat_temp, pr
         gen lambda = (c_hat_temp-log(`delta'))/p_hat_temp
         * run the test
         reg lhs_temp lambda if dep_pos, nocons
         matrix b = e(b)
         ereturn post b
         * drop created variables
         cap drop *temp
		 cap drop xb_temp
		 cap drop dep_pos
		 cap drop lambda
******************************************************************************
*                   Return the information to STATA output		     		 *
******************************************************************************
end
 
