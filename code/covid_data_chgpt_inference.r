
## inference with real data ## 
##############################

set.seed(1)

require(Rcpp)
require(RcppGSL)
require(RcppArmadillo)
require(RcppEigen)
require(RcppNumerical)
require(dplyr)
require(parallel)
require(geosphere)
require(ggplot2)
require(tidyr)
require(grid)
require(gridExtra)
require(boot)
require(rlang)
require(rbenchmark)


ste="GA"
area= "ATL"
n_iter = 2500000 # number of iterations for MCMC


if(ste=="GA"){

	dir_code = "x"
	dir_data = "x"
	dir_raw_data = "x"
	dir_results = "x"
	if(!dir.exists(dir_results)) dir.create(dir_results)

}


Rcpp::sourceCpp(paste(dir_code,"cpp_funcs_covid_data_inference.cpp",sep="")) 



## read data ##
###############

if(ste=='GA') observed_new_cases = read.csv(paste(dir_raw_data, "ga_new_cases_",area, ".csv",sep=""))

n_day_init_EI = 15

# n_age=2

n_age = ncol(observed_new_cases) -1

# init_EI = observed_new_cases[1:15,1:n_age]

if(n_age==1) init_EI = as.matrix(observed_new_cases[1:n_day_init_EI,1:n_age], ncol=n_age)
if(n_age>1) init_EI = observed_new_cases[1:n_day_init_EI,1:n_age]

for (i in 1:n_age){

	mean_i = mean(init_EI[,i][which(init_EI[,i]!=0)])
	init_EI[which(init_EI[,i]==0),i] = trunc(mean_i)


}




para_scalar = as.list(read.csv(file=paste(dir_data,"para_scalar.csv",sep="")))
para_vec_age = as.list(read.csv(file=paste(dir_data,"para_vec_age.csv",sep="")))
para_vec_sero = as.list(read.csv(file=paste(dir_data,"para_vec_sero.csv",sep="")))

epi_list = list()
for( i in 1:para_scalar$n_age){
	epi_list[[i]] = read.csv(file=paste(dir_data,"epi_age_gp_",i, ".csv", sep=""))

	epi_list[[i]]$n_SE =  -c(0,diff(epi_list[[i]]$S))


	epi_list[[i]]$n_EU = c( rep(0,para_scalar$EU_p1), head(epi_list[[i]]$EU1, nrow(epi_list[[i]])-para_scalar$EU_p1) )# n_EU(t) =EU1(t-D_EU)
	epi_list[[i]]$n_EI = c( rep(0,para_scalar$EI_p1), head(epi_list[[i]]$EI1, nrow(epi_list[[i]])-para_scalar$EI_p1) )# n_EI(t) =EI1(t-D_EI)


	n_tmp = 10

	epi_list[[i]]$n_IR = c( 0,(tail(init_EI[,i],para_scalar$IR_p1)), head(epi_list[[i]]$I1[-1], nrow(epi_list[[i]])-para_scalar$IR_p1-1) )# n_IR(t) =I1(t-D_IR)
	epi_list[[i]]$n_UR = c( 0,(trunc(tail(n_tmp*init_EI[,i],para_scalar$UR_p1))), head(epi_list[[i]]$U1[-1], nrow(epi_list[[i]])-para_scalar$UR_p1-1) )# n_UR(t) =U1(t-D_UR)




	##


	epi_list[[i]] = as.matrix(epi_list[[i]] %>% dplyr::select(age_gp, time, S, E, I, U, R, n_SE, n_EU, n_EI, n_UR, n_IR))

}

## check ##


iage = 1

n_SE = epi_list[[iage]][,8]
n_EU = epi_list[[iage]][,9]
n_EI = epi_list[[iage]][,10]
n_UR = epi_list[[iage]][,11]
n_IR = epi_list[[iage]][,12]
S0 = epi_list[[iage]][,3][1]
I0 = epi_list[[iage]][,5][1]
U0 = epi_list[[iage]][,6][1]
R0 = epi_list[[iage]][,7][1]

S = {}
for (t in 1:para_scalar$t_max){

	# if(t==1) S_prev = para_vec_age$N_age[iage]
	if(t==1) S_prev = S0
	if(t>1) S_prev = S[t-1];
	S[t] = S_prev - n_SE[t];

}
 

E_start = 0
E = {}
for (t in 1:para_scalar$t_max){

	if(t==1) E_prev =E_start
	if(t>1) E_prev = E[t-1]

	E[t] = E_prev + n_SE[t] - n_EI[t] - n_EU[t] 
}
 


I = {}
for (t in 1:para_scalar$t_max){

	if(t==1) I_prev = I0
	if(t>1) I_prev = I[t-1];
	I[t] = I_prev + n_EI[t] - n_IR[t] ;

}

U = {}
for (t in 1:para_scalar$t_max){

	if(t==1) U_prev = U0
	if(t>1) U_prev = U[t-1];
	U[t] = U_prev + n_EU[t] - n_UR[t] ;

}


R = {}
for (t in 1:para_scalar$t_max){

	if(t==1) R_prev = R0
	if(t>1) R_prev = R[t-1];
	R[t] = R_prev + n_IR[t] + n_UR[t] ;

}


table(epi_list[[iage]][,3]==S)
table(epi_list[[iage]][,4]==E)
table(epi_list[[iage]][,5]==I)
table(epi_list[[iage]][,6]==U)
table(epi_list[[iage]][,7]==R)



a = b = {}
for( t in 2:(para_scalar$t_max-para_scalar$EI_p1)){
	a[t] = n_SE[t]
	b[t] = n_EI[t+para_scalar$EI_p1] + n_EU[t+para_scalar$EI_p1]
}

table(a==b)


##

sero_list = list()
for( i in 1:length(para_vec_sero$t_sero)){
	sero_list[[i]] = as.matrix(read.csv(file=paste(dir_data,"sero_",i, ".csv", sep="")))
}


contact_mat = as.matrix(read.csv(file=paste(dir_data, "contact_mat.csv", sep="")))
contact_mat_A = as.matrix(read.csv(file=paste(dir_data, "contact_mat_A.csv", sep="")))
contact_mat_B = as.matrix(read.csv(file=paste(dir_data, "contact_mat_B.csv", sep="")))

mt_vec = as.numeric(read.csv(file=paste(dir_data, "mt_vec.csv", sep=""))$x)

## inference ##
###############


log_lh_current = llh_func_tmp( para_scalar[["prop_asym"]], epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results) # important: this will actually change permanently the par value!!


 mcmc_HH(n_iter, epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A,contact_mat_B, dir_results, log_lh_current)

