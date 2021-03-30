{smcl}
{* *! version 1.0 22march2021}{...}
{vieweralsosee "[R] poisson" "help poisson"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ppml" "help ppml"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{viewerjumpto "Syntax" "iOLS_ivreg2##syntax"}{...}
{viewerjumpto "Description" "iOLS_ivreg2##description"}{...}
{viewerjumpto "Citation" "iOLS_ivreg2##citation"}{...}
{viewerjumpto "Authors" "iOLS_ivreg2##contact"}{...}
{viewerjumpto "Examples" "iOLS_ivreg2##examples"}{...}
{viewerjumpto "Description" "iOLS_ivreg2##Testing"}{...}
{viewerjumpto "Stored results" "iOLS_ivreg2##results"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:iOLS_ivreg2} {hline 2}} Iterated Two Stage Least Squares i2SLS using ivreg2 {p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:iOLS_OLS}
{depvar} [{indepvars}]
{ifin} {it:{weight}} {cmd:,} [{help iOLS_ivreg2##options:options}] {p_end}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Model}
{syntab:SE:Robust/Cluster}
{synopt:{opt vce}{cmd:(}{help iOLS_OLS##opt_vce:vcetype}{cmd:)}}{it:vcetype}
may be {opt r:obust} (default) or {opt cl:uster} {help fvvarlist} (allowing two- and multi-way clustering){p_end}
{syntab:Delta}
{synopt:}{it:delta} is any strictly positive constant. {p_end}


{marker description}{...}
{title:Description}

{pstd} {cmd: i2SLS_ivreg2} estimates iOLS_delta, a solution to the problem of the log of zero.  This method relies on running the "regress" function iteratively.
This provides the reader with the final OLS estimates and allows the use the post-estimation commands available under regress (using Y_tilde = log(Y + delta*exp(xb))) as a 
dependent variable.  The benefit of using "regress" comes at the cost of limited capacity to deal with many fixed effects. In such a case, iOLS_hdfe may be more appropriate.

{pstd}This package:

{pmore} 1. relies on Stata's ivreg2 procedure for estimation.{p_end}

{pmore} 2. assumes the iOLS exogeneity condition with instrument Z E(Z'log(delta+U))=constant  {p_end}


{title:Background}

{pstd} i2SLS_delta is a solution to the problem of the log of zero with endogenous variables.

{marker absvar}{...}
{title:Syntax for absorbed variables}

{synoptset 22}{...}
{synopthdr: variables}
{synoptline}
{synopt:{it:depvar}} Dependent variable{p_end}
{synopt:{it:indepvars}} List of explanatory variables {p_end}
{synoptline}
{p2colreset}{...}


{marker caveats}{...}
{title:Caveats}

{pstd} Convergence is decided based on coefficients (sum of squared coefficients < 1e-6) and not on the modulus of the contraction mapping.


{pstd}The {help reg postestimation##predict:predict}, {help test}, and {help margins} postestimation commands are available after {cmd:iOLS_OLS}.


{marker contact}{...}
{title:Authors}

{pstd} Author {break}
CREST {break}
Email: {browse "mailto:author@author.edu":author@email.edu}
{p_end}




{marker citation}{...}
{title:Citation}

{pstd}
Citation to be defined. 


{marker examples}{...}
{title:Examples}

{pstd}First, we compare iOLS_ivreg2 with IV-Poisson (using ivpois)
{browse "http://www.haghish.com/statistics/stata-blog/stata-programming/download/ivpois.html":ivpois help file}.
{p_end}
{hline}
{phang2}{cmd:. use "http://www.stata-press.com/data/r14/airline"}{p_end}
{phang2}{cmd:. iOLS_ivreg2 injuries (XYZowned=n) , delta(1) vce(robust)}{p_end}
{phang2}{cmd:. ivpois injuries (XYZowned=n)}{p_end}
{hline}


{pstd} Second, we show how to test for the pattern of zeros with i2SLS. We use data on womens' participation into the labor force and look at returns to education.
{p_end}
{hline}
{phang2}{cmd:. clear all}{p_end}
{phang2}{cmd:. webuse womenwk }{p_end}
{phang2}{cmd:. replace wage = 0 if missing(wage) }{p_end}
{phang2}{cmd:. gen log_wage = log(wage) }{p_end}
{phang2}{cmd:. gen employment = wage!=0 }{p_end}

i2SLS_ivreg2 wage age children (education = married) , delta(1) robust
i2SLS_ivreg2 wage age children (education = married) , delta(1000) robust 
ivpois wage age children, endog(education) exog(married)


{phang2}{cmd:. i2SLS_ivreg2 wage age children (education = married) , delta(1) robust }{p_end}
{phang2}{cmd:. cap program drop i2SLS_bootstrap  }{p_end}
{phang2}{cmd:. program i2SLS_bootstrap, rclass  }{p_end}
{phang2}{cmd:. i2SLS_ivreg2 wage age children (education = married) , delta(1) robust  }{p_end}
{phang2}{cmd:. scalar delta = 1  }{p_end}
{phang2}{cmd:. *lhs of test  }{p_end}
{phang2}{cmd:. predict xb_temp, xb  }{p_end}
{phang2}{cmd:. gen u_hat_temp = wage*exp(-xb_temp)  }{p_end}
{phang2}{cmd:. gen lhs_temp = log(delta+u_hat_temp) - log(delta)  }{p_end}
{phang2}{cmd:. * rhs of test  }{p_end}
{phang2}{cmd:. gen temp = log(wage + delta*exp(xb_temp)) - xb_temp  }{p_end}
{phang2}{cmd:. egen c_hat_temp = mean(temp)   }{p_end}
{phang2}{cmd:. logit employment education age  }{p_end}
{phang2}{cmd:. predict p_hat_temp, pr  }{p_end}
{phang2}{cmd:. gen rhs_temp = (c_hat_temp-log(delta))/p_hat_temp  }{p_end}
{phang2}{cmd:. * run the test  }{p_end}
{phang2}{cmd:. reg lhs_temp rhs_temp if employment, nocons   }{p_end}
{phang2}{cmd:. matrix b = e(b)  }{p_end}
{phang2}{cmd:. ereturn post b  }{p_end}
{phang2}{cmd:. * drop created variables  }{p_end}
{phang2}{cmd:. cap drop *temp  }{p_end}
{phang2}{cmd:. end  }{p_end}

{phang2}{cmd:. bootstrap lambda = _b[rhs_temp] , reps(50): i2SLS_bootstrap  }{p_end}
{phang2}{cmd:. test lambda==1  }{p_end}
{hline}

{pstd} Third, you can convert your results into latex using esttab where "eps" provides the convergence criteria:
{p_end}
{hline}
{phang2}{cmd:. eststo clear}{p_end}
{phang2}{cmd:. eststo: i2SLS_ivreg2 wage age children (education = married) , delta(1) robust}{p_end}
{phang2}{cmd:. eststo: i2SLS_ivreg2 wage age children (education = married) , delta(10) robust}{p_end}
{phang2}{cmd:. eststo: esttab * using table_iv.tex,  scalars(delta eps) }{p_end}
{hline}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:iOLS_ivreg2} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(sample)}} marks the sample used for estimation {p_end}
{synopt:{cmd:e(eps)}} sum of the absolute differences between the parameters from the last two iterations of iOLS {p_end}
{synopt:{cmd:e(k)}} number of iterations of iOLS{p_end}

