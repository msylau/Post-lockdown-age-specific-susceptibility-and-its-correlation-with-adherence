//#include <Rcpp.h> 
#include <gsl/gsl_rng.h> 
#include <gsl/gsl_randist.h>
#include <RcppArmadillo.h>
#include <omp.h>
#include <fstream>
#include <RcppNumerical.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>


using namespace Numer;
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]


//-- testing for uisng mop --//

// [[Rcpp::export]]
arma::mat inner(arma::vec x, arma::vec y){
  arma::mat ip=x.t()*y;
  return(ip);
}

double sq(double x){
  return(x*x);
}


// [[Rcpp::export]]
double sumsq_serial(arma::vec x)
{
  double sum = 0.0;
  for (int i=0; i<x.size(); i++){
    sum += sq(x(i));
  }  
  return sum;
}



// [[Rcpp::export]]
double sumsq_parallel(arma::vec x, int ncores)
{
  double sum = 0.0;
  omp_set_num_threads(ncores);
  #pragma omp parallel for shared(x) reduction(+:sum)
  for (int i=0; i<x.size(); i++){
	// Rcout << i << ", " << omp_get_thread_num() <<  " of " <<  omp_get_num_threads() << endl;
    sum += sq(x(i));
  }
  return sum;
}

//---

double func_time_beta (double t, double t_intervention, double t_intervention_2 , double t_intervention_3 ,double sum_beta, double alpha, double omega_1, double omega_2,double omega_3){
// model time-varying beta//

	double beta_t= 0.0;

	// if (t>=t_intervention) beta_t = (sum_beta)*exp(-1*(omega_1)*(t- t_intervention));
	// if(t<t_intervention) beta_t = sum_beta;

	// if(t<t_intervention) beta_t = sum_beta;
	// if (t>=t_intervention) beta_t = (sum_beta)*omega_1;


	if(t<t_intervention) beta_t = sum_beta;
	if (t>=t_intervention & t<t_intervention_2) beta_t = (sum_beta)*omega_1;
	if (t>=t_intervention_2) beta_t = (sum_beta)*omega_2;

	// if(t<t_intervention) beta_t = sum_beta;
	// if (t>=t_intervention & t<t_intervention_2) beta_t = (sum_beta)*omega_1;
	// if (t>=t_intervention_2 & t<t_intervention_3) beta_t = (sum_beta)*omega_2;
	// if (t>=t_intervention_3) beta_t = (sum_beta)*omega_3;


return(beta_t);

}


double compute_ij_beta(int& i_age, int& j_age, Rcpp::List& epi_list, NumericMatrix& contact_mat, NumericVector& beta_sym_vec,NumericVector& gamma_vec, int t, Rcpp::List& para_list){


	// for age gp i (infectee) //
	Rcpp::NumericMatrix epi_i = Rcpp::as<Rcpp::NumericMatrix>(epi_list[i_age]);


	NumericVector S_i = epi_i(_,2);
	NumericVector E_i = epi_i(_,3);
	NumericVector I_i = epi_i(_,4);
	NumericVector U_i = epi_i(_,5);
	NumericVector R_i = epi_i(_,6);


	// for age gp j (infector)  //
	Rcpp::NumericMatrix epi_j = Rcpp::as<Rcpp::NumericMatrix>(epi_list[j_age]);

	NumericVector S_j = epi_j(_,2);
	NumericVector E_j = epi_j(_,3);
	NumericVector I_j = epi_j(_,4);
	NumericVector U_j = epi_j(_,5);
	NumericVector R_j = epi_j(_,6);

	//

	double sum_beta = 0.0;

	// sum_beta = Rcpp::as<double>(para_list["beta_sym"])*S[t]*I[t] + Rcpp::as<double>(para_list["beta_sym"])*Rcpp::as<double>(para_list["scale_u"])*S[t]*U[t];

	double contact = contact_mat(i_age, j_age);

	// double tmp = (beta_sym_vec[j_age])*(gamma_vec[i_age])*contact*S_i[t-1];
	// sum_beta = tmp*I_j[t-1] + tmp*Rcpp::as<double>(para_list["scale_u"])*U_j[t-1];

	double tmp = (beta_sym_vec[j_age])*(gamma_vec[i_age])*contact*S_i[t];
	sum_beta = tmp*I_j[t] + tmp*Rcpp::as<double>(para_list["scale_u"])*U_j[t];



return(sum_beta);

}



//--

double llh_func( Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A, NumericMatrix& contact_mat_B, std::string dir_results){


	int t_max =  Rcpp::as<int>(para_scalar["t_max"]);
	int t_mask =  Rcpp::as<int>(para_scalar["t_mask"]);

	int n_age =  Rcpp::as<int>(para_scalar["n_age"]); // number of age gps

	int EU_p1 = Rcpp::as<int>(para_scalar["EU_p1"]);
	int EI_p1 = Rcpp::as<int>(para_scalar["EI_p1"]);
	int UR_p1 = Rcpp::as<int>(para_scalar["UR_p1"]);
	int IR_p1 = Rcpp::as<int>(para_scalar["IR_p1"]);


	NumericVector N_age =  para_vec_age["N_age"]; // pop in each age gp

	NumericVector beta_sym_vec = para_vec_age["beta_sym_vec"];
	NumericVector gamma_vec = para_vec_age["gamma_vec"];
	NumericVector gamma_vec_2 = para_vec_age["gamma_vec_2"];

	NumericVector t_sero = para_vec_sero["t_sero"];
	NumericVector n_sero = para_vec_sero["n_sero"];

	int m_sero = t_sero.size(); // number of serology testing time points

 	double llh = 0.0;

// omp_set_num_threads(2);
// #pragma omp parallel for shared(epi_list, para_scalar, para_vec_age, para_vec_sero,contact_mat, beta_sym_vec, gamma_vec, N_age, t_max, n_age) schedule(dynamic)

	for (int i_age=0; i_age<=(n_age-1); i_age++){

// double thread_llh =0.0;

		// epi data of a particular age gp: columns age_gp, time, S, E, I, U, R, n_SE, n_EU, n_EI, n_UR, n_IR //
		Rcpp::NumericMatrix epi = Rcpp::as<Rcpp::NumericMatrix>(epi_list[i_age]);

		// double gamma = gamma_vec[i_age];

		NumericVector S = epi(_,2);
		NumericVector E = epi(_,3);
		NumericVector I = epi(_,4);
		NumericVector U = epi(_,5);
		NumericVector R = epi(_,6);

		NumericVector n_SE = epi(_,7);
		NumericVector n_EU = epi(_,8);
		NumericVector n_EI = epi(_,9);
		NumericVector n_UR = epi(_,10);
		NumericVector n_IR = epi(_,11);

		// calculate likeihood //

		for (int t=1;t<=(t_max-1);t++){

			NumericVector gammas;
			NumericMatrix contact;

			switch(t<=t_mask){
				case 1:{
					gammas = gamma_vec;

					contact = contact_mat_A;
					// contact = contact_mat;

					break;
				}

				case 0:{
					gammas = gamma_vec_2;
					// gammas = gamma_vec;

					contact = contact_mat_B;
					// contact = contact_mat;

					break;
				}
			}

			double sum_beta_internal = 0.0;
			for(int j_age=0; j_age<=(n_age-1); j_age++){
				sum_beta_internal = sum_beta_internal  +  compute_ij_beta( i_age, j_age, epi_list, contact, beta_sym_vec, gammas, t-1, para_scalar)/N_age[j_age];


			}

			double sum_beta =  sum_beta_internal;


			sum_beta = func_time_beta(t-1, Rcpp::as<double>(para_scalar["t_intervention"]), Rcpp::as<double>(para_scalar["t_intervention_2"]),Rcpp::as<double>(para_scalar["t_intervention_3"]), sum_beta, Rcpp::as<double>(para_scalar["alpha"]), Rcpp::as<double>(para_scalar["omega_1"]), Rcpp::as<double>(para_scalar["omega_2"]),Rcpp::as<double>(para_scalar["omega_3"]));


			llh = llh + R::dpois( n_SE[t], sum_beta, TRUE); // R::dpois( x, lambda, log )
// thread_llh += R::dpois( n_SE[t], sum_beta, TRUE); // R::dpois( x, lambda, log )


			if( ((t+EU_p1)<=(t_max-1))  & (t>EU_p1)){


				double f_mt =0.0; 


				if(i_age==0) f_mt = Rcpp::as<double>(para_scalar["a"]) + Rcpp::as<double>(para_scalar["b"])*t; 
				if(i_age==1) f_mt = Rcpp::as<double>(para_scalar["a_2"]) + Rcpp::as<double>(para_scalar["b_2"])*t; 
				if(i_age==2) f_mt = Rcpp::as<double>(para_scalar["a_3"]) + Rcpp::as<double>(para_scalar["b_3"])*t; 
				if(i_age==3) f_mt = Rcpp::as<double>(para_scalar["a_4"]) + Rcpp::as<double>(para_scalar["b_4"])*t; 


				double p_i =  exp(f_mt)/(1+exp(f_mt)); 

				double p_u = 1.0 - p_i; 

// Rcout << n_EU[t] << "," << n_SE[t-EU_p1] << "," << p_u << "," <<  R::dbinom(n_EU[t] , n_SE[t-EU_p1], p_u, TRUE)<< endl;

				// llh = llh + R::dbinom(n_EU[t] , n_SE[t-EU_p1], p_u, TRUE);

				llh = llh + R::dbinom(n_EU[t+EU_p1] , n_SE[t], p_u, TRUE);
// thread_llh += R::dbinom(n_EU[t+EU_p1] , n_SE[t], p_u, TRUE);

			}



		} // end loop of  (int t=1;t<(t_max-1);t++)

// llh += thread_llh;

	} // end loop of (int i_age=0; i_age<=(n_age-1); i_age++)



return(llh);

}



// [[Rcpp::export]]
double llh_func_tmp( double par, Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A, NumericMatrix& contact_mat_B, std::string dir_results){


	double llh = 0.0;

	int t_max =  Rcpp::as<int>(para_scalar["t_max"]);

	int t_mask =  Rcpp::as<int>(para_scalar["t_mask"]);

	int n_age =  Rcpp::as<int>(para_scalar["n_age"]); // number of age gps

	int EU_p1 = Rcpp::as<int>(para_scalar["EU_p1"]);
	int EI_p1 = Rcpp::as<int>(para_scalar["EI_p1"]);
	int UR_p1 = Rcpp::as<int>(para_scalar["UR_p1"]);
	int IR_p1 = Rcpp::as<int>(para_scalar["IR_p1"]);


	NumericVector N_age =  para_vec_age["N_age"]; // pop in each age gp

	NumericVector beta_sym_vec = para_vec_age["beta_sym_vec"];
	NumericVector gamma_vec = para_vec_age["gamma_vec"];
	NumericVector gamma_vec_2 = para_vec_age["gamma_vec_2"];

	NumericVector t_sero = para_vec_sero["t_sero"];
	NumericVector n_sero = para_vec_sero["n_sero"];

	int m_sero = t_sero.size(); // number of serology testing time points


 para_scalar["prop_asym"]=  par; // the parameter to optimize


// NumericVector thread_llh_1(n_age,0.0);
// NumericVector thread_llh_2(n_age,0.0);
// omp_set_num_threads(2);
// #pragma omp parallel for 

		for (int i_age=0; i_age<=(n_age-1); i_age++){

			// epi data of a particular age gp: columns age_gp, time, S, E, I, U, R, n_SE, n_EU, n_EI, n_UR, n_IR //
			Rcpp::NumericMatrix epi = Rcpp::as<Rcpp::NumericMatrix>(epi_list[i_age]);

			// double gamma = gamma_vec[i_age];

			NumericVector S = epi(_,2);
			NumericVector E = epi(_,3);
			NumericVector I = epi(_,4);
			NumericVector U = epi(_,5);
			NumericVector R = epi(_,6);

			NumericVector n_SE = epi(_,7);
			NumericVector n_EU = epi(_,8);
			NumericVector n_EI = epi(_,9);
			NumericVector n_UR = epi(_,10);
			NumericVector n_IR = epi(_,11);

			//--calculate likeihood--//

			for (int t=1;t<=(t_max-1);t++){


				NumericVector gammas;
				NumericMatrix contact;

				switch(t<=t_mask){
					case 1:{
						gammas = gamma_vec;

						contact = contact_mat_A;
						// contact = contact_mat;

						break;
					}

					case 0:{
						gammas = gamma_vec_2;
						// gammas = gamma_vec;

						contact = contact_mat_B;
						// contact = contact_mat;

						break;
					}
				}

				double sum_beta_internal = 0.0;
				for(int j_age=0; j_age<=(n_age-1); j_age++){
					sum_beta_internal = sum_beta_internal  +  compute_ij_beta( i_age, j_age, epi_list, contact, beta_sym_vec, gammas, t-1, para_scalar)/N_age[j_age];
				}

				double sum_beta =  sum_beta_internal;
				// double sum_beta =  sum_beta_internal/sum(N_age);
				// double sum_beta =  sum_beta_internal/N_age[i_age];

				// sum_beta = func_time_beta(t-1, Rcpp::as<double>(para_scalar["t_intervention"]), sum_beta, Rcpp::as<double>(para_scalar["alpha"]), Rcpp::as<double>(para_scalar["omega_1"]), Rcpp::as<double>(para_scalar["omega_2"]));
				// sum_beta = func_time_beta(t-1, Rcpp::as<double>(para_scalar["t_intervention"]), Rcpp::as<double>(para_scalar["t_intervention_2"]), sum_beta, Rcpp::as<double>(para_scalar["alpha"]), Rcpp::as<double>(para_scalar["omega_1"]), Rcpp::as<double>(para_scalar["omega_2"]));
				sum_beta = func_time_beta(t-1, Rcpp::as<double>(para_scalar["t_intervention"]), Rcpp::as<double>(para_scalar["t_intervention_2"]),Rcpp::as<double>(para_scalar["t_intervention_3"]), sum_beta, Rcpp::as<double>(para_scalar["alpha"]), Rcpp::as<double>(para_scalar["omega_1"]), Rcpp::as<double>(para_scalar["omega_2"]),Rcpp::as<double>(para_scalar["omega_3"]));

				llh = llh + R::dpois( n_SE[t], sum_beta, TRUE); // R::dpois( x, lambda, log )
				// thread_llh_1[i_age] =  R::dpois( n_SE[t], sum_beta, TRUE); // R::dpois( x, lambda, log )

// Rcout <<  t << "," << sum_beta << "," << n_SE[t] << "," << R::dpois( n_SE[t], sum_beta, TRUE) << endl;

				// if(t>EU_p1){
				// if((t+EU_p1)<=(t_max-1)){
				if( ((t+EU_p1)<=(t_max-1))  & (t>EU_p1)){


				double f_mt =0.0; 
				if(i_age==0) f_mt = Rcpp::as<double>(para_scalar["a"]) + Rcpp::as<double>(para_scalar["b"])*t; 
				if(i_age==1) f_mt = Rcpp::as<double>(para_scalar["a_2"]) + Rcpp::as<double>(para_scalar["b_2"])*t; 
				if(i_age==2) f_mt = Rcpp::as<double>(para_scalar["a_3"]) + Rcpp::as<double>(para_scalar["b_3"])*t; 
				if(i_age==3) f_mt = Rcpp::as<double>(para_scalar["a_4"]) + Rcpp::as<double>(para_scalar["b_4"])*t; 

					// double p_i =  (1.0-Rcpp::as<double>(para_scalar["prop_asym"])) + Rcpp::as<double>(para_scalar["prop_asym"])*(exp(f_mt)/(1+exp(f_mt))); // the prob of detected: either has symptoms or no symptom but tested (exp(f_mt)/(1+exp(f_mt)) is the conidtional prob of tested given infection)
					// double p_i =  1.0-Rcpp::as<double>(para_scalar["prop_asym"]); 
					double p_i =  exp(f_mt)/(1+exp(f_mt)); 

					double p_u = 1.0 - p_i; 

	// Rcout << n_EU[t] << "," << n_SE[t-EU_p1] << "," << p_u << "," <<  R::dbinom(n_EU[t] , n_SE[t-EU_p1], p_u, TRUE)<< endl;

					// llh = llh + R::dbinom(n_EU[t] , n_SE[t-EU_p1], p_u, TRUE);

					llh = llh + R::dbinom(n_EU[t+EU_p1] , n_SE[t], p_u, TRUE);
					// thread_llh_2[i_age] = R::dbinom(n_EU[t+EU_p1] , n_SE[t], p_u, TRUE);

				}


			} // end loop of  (int t=1;t<(t_max-1);t++)

		} // end loop of (int i_age=0; i_age<=(n_age-1); i_age++)

// llh = sum(thread_llh_1) + sum(thread_llh_2);



return(llh);

}




//--- MCMC -- //

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}



void mcmc_gammas(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat,NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){



	NumericVector gamma_vec = para_vec_age["gamma_vec"];

	// arma::vec mu = {0.0,0.0,0.0.0.0};
	arma::vec mu = {gamma_vec[0], gamma_vec[1],gamma_vec[2]};

	// double cor_12 = 0.1;
	// double cor_13 =  0.1;
	// double cor_14 = 0.1;
	double cor_23 = 0.1;
	double cor_24 = 0.1;
	double cor_34 = 0.1;

	// arma::mat cor = { {1,cor_12, cor_13, cor_14},
	//  {cor_12, 1, cor_23, cor_24},
	//  {cor_13,cor_23,1,cor_34},
	//  {cor_14,cor_24,cor_34, 1} };

	arma::mat cor = { {1,cor_23, cor_24},
	 {cor_23, 1, cor_34},
	 {cor_24,cor_34, 1} };

	// double sd_1 = 0.01;
	double sd_2 = 0.1;
	double sd_3 = 0.1;
	double sd_4 = 0.1;

	// arma::mat sd = { {sd_1*sd_1,sd_1*sd_2, sd_1*sd_3, sd_1*sd_4},
	//  {sd_2*sd_1, sd_2*sd_2, sd_2*sd_3, sd_2*sd_4},
	//  {sd_3*sd_1,sd_3*sd_2,sd_3*sd_3,sd_3*sd_4},
	//  {sd_4*sd_1,sd_4*sd_2,sd_4*sd_3, sd_4*sd_4} };


	arma::mat sd = { {sd_2*sd_2,sd_2*sd_3, sd_2*sd_4},
	 				{sd_3*sd_2,sd_3*sd_3, sd_3*sd_4},
	 				{sd_4*sd_2,sd_4*sd_3, sd_4*sd_4} };


	arma::mat sigma = (0.01*sd) % cor; // element wise multiplication

	NumericVector gamma_vec_new = clone(gamma_vec);
	gamma_vec_new[3] = gamma_vec[3]; // this wont change

	gamma_vec_new[1] = -9999;
	while(gamma_vec_new[1]<0 | gamma_vec_new[2]<0 |gamma_vec_new[3]<0  ){

		arma::mat z_arma_mat = mvrnormArma(1, mu, sigma); 

		NumericMatrix z_mat = Rcpp::as<Rcpp::NumericMatrix>(wrap(z_arma_mat));

		NumericVector z = z_mat(0,_);

		gamma_vec_new[0] = z[0];
		gamma_vec_new[1] = z[1];
		gamma_vec_new[2] = z[2];


	}


	Rcpp::List para_vec_age_new =clone(para_vec_age);
	para_vec_age_new["gamma_vec"] = gamma_vec_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar, para_vec_age_new, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


// Rcout <<  log_lh_current << "," << log_lh_new << endl;


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

// Rcout <<  acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1)){

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				para_vec_age["gamma_vec"] = gamma_vec_new ;
			break;
			}
			case 0: {
			break;
			}
		}
	}
}





void mcmc_gammas_mask_full(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){



	NumericVector gamma_vec = para_vec_age["gamma_vec_2"];

	// arma::vec mu = {0.0,0.0,0.0.0.0};
	arma::vec mu = {gamma_vec[0], gamma_vec[1], gamma_vec[2],gamma_vec[3]};

	double cor_12 = 0.1;
	double cor_13 =  0.1;
	double cor_14 = 0.1;
	double cor_23 = 0.1;
	double cor_24 = 0.1;
	double cor_34 = 0.1;

	arma::mat cor = { {1,cor_12, cor_13, cor_14},
	 {cor_12, 1, cor_23, cor_24},
	 {cor_13,cor_23,1,cor_34},
	 {cor_14,cor_24,cor_34, 1} };

	// arma::mat cor = { {1,cor_23, cor_24},
	//  {cor_23, 1, cor_34},
	//  {cor_24,cor_34, 1} };

	double sd_1 = 0.1;
	double sd_2 = 0.1;
	double sd_3 = 0.1;
	double sd_4 = 0.1;

	arma::mat sd = { {sd_1*sd_1,sd_1*sd_2, sd_1*sd_3, sd_1*sd_4},
	 {sd_2*sd_1, sd_2*sd_2, sd_2*sd_3, sd_2*sd_4},
	 {sd_3*sd_1,sd_3*sd_2,sd_3*sd_3,sd_3*sd_4},
	 {sd_4*sd_1,sd_4*sd_2,sd_4*sd_3, sd_4*sd_4} };


	// arma::mat sd = { {sd_2*sd_2,sd_2*sd_3, sd_2*sd_4},
	//  				{sd_3*sd_2,sd_3*sd_3, sd_3*sd_4},
	//  				{sd_4*sd_2,sd_4*sd_3, sd_4*sd_4} };

	arma::mat sigma = (0.01*sd) % cor; // element wise multiplication

	NumericVector gamma_vec_new = clone(gamma_vec);
	// gamma_vec_new[0] = gamma_vec[0]; // this wont change

	gamma_vec_new[1] = -9999;
	while(gamma_vec_new[0]<0  | gamma_vec_new[1]<0 | gamma_vec_new[2]<0 |gamma_vec_new[3]<0  ){

		arma::mat z_arma_mat = mvrnormArma(1, mu, sigma); 

		NumericMatrix z_mat = Rcpp::as<Rcpp::NumericMatrix>(wrap(z_arma_mat));

		NumericVector z = z_mat(0,_);

		gamma_vec_new[0] = z[0];
		gamma_vec_new[1] = z[1];
		gamma_vec_new[2] = z[2];
		gamma_vec_new[3] = z[3];


	}


	Rcpp::List para_vec_age_new =clone(para_vec_age);
	para_vec_age_new["gamma_vec_2"] = gamma_vec_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar, para_vec_age_new, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


// Rcout <<  log_lh_current << "," << log_lh_new << endl;


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

// Rcout <<  acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1)){

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				para_vec_age["gamma_vec_2"] = gamma_vec_new ;
			break;
			}
			case 0: {
			break;
			}
		}
	}
}





void mcmc_beta(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	NumericVector beta_sym_vec = para_vec_age["beta_sym_vec"];

	NumericVector beta_sym_vec_new = clone(beta_sym_vec);

	beta_sym_vec_new[0] = -9999;
	// while(beta_sym_vec_new[0]<0){
		 // beta_sym_vec_new[0] = beta_sym_vec[0] + 0.01*R::rnorm(0.0,1.0);
		 beta_sym_vec_new[0] = exp(log(beta_sym_vec[0]) + 0.1*R::rnorm(0.0,1.0)); //propose on log-scale

	// }


	beta_sym_vec_new[1] = beta_sym_vec_new[2] = beta_sym_vec_new[3] = beta_sym_vec_new[0];

	Rcpp::List para_vec_age_new =clone(para_vec_age);
	para_vec_age_new["beta_sym_vec"] = beta_sym_vec_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar, para_vec_age_new, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


	double log_proposal = 0.0;
	log_proposal = log(beta_sym_vec_new[0]) - log(beta_sym_vec[0]); // if proposed on log-scale

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current + log_proposal));

// Rcout <<  acp_pr << endl;

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_vec_age["beta_sym_vec"] = beta_sym_vec_new ;
		break;
		}
		case 0: {
		break;
		}
	}

}









void mcmc_b_a(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){
// jointly propose b and prop_asym//


	arma::vec mu = {para_scalar["b"], para_scalar["a"]};

	double cor_12 = -0.2;


	// arma::mat cor = { {1,cor_12, cor_13, cor_14},
	//  {cor_12, 1, cor_23, cor_24},
	//  {cor_13,cor_23,1,cor_34},
	//  {cor_14,cor_24,cor_34, 1} };

	arma::mat cor = { {1,cor_12},
	 {cor_12, 1} };

	double sd_1 = 0.001;
	double sd_2 = 0.02;
	

	arma::mat sd = { {sd_1*sd_1,sd_1*sd_2},
	 {sd_2*sd_1, sd_2*sd_2} };

	arma::mat sigma = sd % cor; // element wise multiplication


	double b_new = -9999;
	double a_new = -9999;

	double f_mt = -9999; 
	double p_i =  0.0; 


	while(b_new<0){
	// while(b_new<0 | a_new>0){

		arma::mat z_arma_mat = mvrnormArma(1, mu, sigma); 

		NumericMatrix z_mat = Rcpp::as<Rcpp::NumericMatrix>(wrap(z_arma_mat));

		NumericVector z = z_mat(0,_);

		b_new = z[0];
		a_new = z[1];


		// f_mt =  a_new + b_new*90; 
		// p_i =  exp(f_mt)/(1+exp(f_mt)); 


	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["b"] = b_new;
	para_scalar_new["a"] = a_new;


	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


// Rcout <<  log_lh_current << "," << log_lh_new << endl;


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

// Rcout <<  acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1)){

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				para_scalar["b"] = b_new ;
				para_scalar["a"] = a_new ;
			break;
			}
			case 0: {
			break;
			}
		}
	}

}


void mcmc_b_2_a_2(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){
// jointly propose b and prop_asym//


	arma::vec mu = {para_scalar["b_2"], para_scalar["a_2"]};

	double cor_12 = -0.2;


	// arma::mat cor = { {1,cor_12, cor_13, cor_14},
	//  {cor_12, 1, cor_23, cor_24},
	//  {cor_13,cor_23,1,cor_34},
	//  {cor_14,cor_24,cor_34, 1} };

	arma::mat cor = { {1,cor_12},
	 {cor_12, 1} };

	double sd_1 = 0.001;
	double sd_2 = 0.02;
	

	arma::mat sd = { {sd_1*sd_1,sd_1*sd_2},
	 {sd_2*sd_1, sd_2*sd_2} };

	arma::mat sigma = sd % cor; // element wise multiplication


	double b_new = -9999;
	double a_new = -9999;

	double f_mt = -9999; 
	double p_i =  0.0; 


	while(b_new<0){
	// while(b_new<0 | a_new>0){

		arma::mat z_arma_mat = mvrnormArma(1, mu, sigma); 

		NumericMatrix z_mat = Rcpp::as<Rcpp::NumericMatrix>(wrap(z_arma_mat));

		NumericVector z = z_mat(0,_);

		b_new = z[0];
		a_new = z[1];


		// f_mt =  a_new + b_new*90; 
		// p_i =  exp(f_mt)/(1+exp(f_mt)); 


	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["b_2"] = b_new;
	para_scalar_new["a_2"] = a_new;


	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


// Rcout <<  log_lh_current << "," << log_lh_new << endl;


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

// Rcout <<  acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1)){

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				para_scalar["b_2"] = b_new ;
				para_scalar["a_2"] = a_new ;
			break;
			}
			case 0: {
			break;
			}
		}
	}

}



void mcmc_b_3_a_3(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){
// jointly propose b and prop_asym//


	arma::vec mu = {para_scalar["b_3"], para_scalar["a_3"]};

	double cor_12 = -0.2;


	// arma::mat cor = { {1,cor_12, cor_13, cor_14},
	//  {cor_12, 1, cor_23, cor_24},
	//  {cor_13,cor_23,1,cor_34},
	//  {cor_14,cor_24,cor_34, 1} };

	arma::mat cor = { {1,cor_12},
	 {cor_12, 1} };

	double sd_1 = 0.001;
	double sd_2 = 0.02;
	

	arma::mat sd = { {sd_1*sd_1,sd_1*sd_2},
	 {sd_2*sd_1, sd_2*sd_2} };

	arma::mat sigma = sd % cor; // element wise multiplication


	double b_new = -9999;
	double a_new = -9999;

	double f_mt = -9999; 
	double p_i =  0.0; 


	while(b_new<0){
	// while(b_new<0 | a_new>0){

		arma::mat z_arma_mat = mvrnormArma(1, mu, sigma); 

		NumericMatrix z_mat = Rcpp::as<Rcpp::NumericMatrix>(wrap(z_arma_mat));

		NumericVector z = z_mat(0,_);

		b_new = z[0];
		a_new = z[1];


		// f_mt =  a_new + b_new*90; 
		// p_i =  exp(f_mt)/(1+exp(f_mt)); 


	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["b_3"] = b_new;
	para_scalar_new["a_3"] = a_new;


	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


// Rcout <<  log_lh_current << "," << log_lh_new << endl;


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

// Rcout <<  acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1)){

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				para_scalar["b_3"] = b_new ;
				para_scalar["a_3"] = a_new ;
			break;
			}
			case 0: {
			break;
			}
		}
	}

}


void mcmc_b_4_a_4(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){
// jointly propose b and prop_asym//


	arma::vec mu = {para_scalar["b_4"], para_scalar["a_4"]};

	double cor_12 = -0.2;


	// arma::mat cor = { {1,cor_12, cor_13, cor_14},
	//  {cor_12, 1, cor_23, cor_24},
	//  {cor_13,cor_23,1,cor_34},
	//  {cor_14,cor_24,cor_34, 1} };

	arma::mat cor = { {1,cor_12},
	 {cor_12, 1} };

	double sd_1 = 0.001;
	double sd_2 = 0.02;
	

	arma::mat sd = { {sd_1*sd_1,sd_1*sd_2},
	 {sd_2*sd_1, sd_2*sd_2} };

	arma::mat sigma = sd % cor; // element wise multiplication


	double b_new = -9999;
	double a_new = -9999;

	double f_mt = -9999; 
	double p_i =  0.0; 


	while(b_new<0){
	// while(b_new<0 | a_new>0){

		arma::mat z_arma_mat = mvrnormArma(1, mu, sigma); 

		NumericMatrix z_mat = Rcpp::as<Rcpp::NumericMatrix>(wrap(z_arma_mat));

		NumericVector z = z_mat(0,_);

		b_new = z[0];
		a_new = z[1];


		// f_mt =  a_new + b_new*90; 
		// p_i =  exp(f_mt)/(1+exp(f_mt)); 


	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["b_4"] = b_new;
	para_scalar_new["a_4"] = a_new;


	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


// Rcout <<  log_lh_current << "," << log_lh_new << endl;


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

// Rcout <<  acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1)){

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				para_scalar["b_4"] = b_new ;
				para_scalar["a_4"] = a_new ;
			break;
			}
			case 0: {
			break;
			}
		}
	}

}


void mcmc_omega_1(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat,NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	double omega_1_new = -9999;

	// while(omega_1_new<0 | omega_1_new>1){
	while(omega_1_new<0){

		omega_1_new = Rcpp::as<double>(para_scalar["omega_1"])+ 0.1*R::rnorm(0.0,1.0);
	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["omega_1"] = omega_1_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A,contact_mat_B, dir_results);

// Rcout <<  log_lh_current << "," << log_lh_new << endl;

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_scalar["omega_1"] = omega_1_new ;
		break;
		}
		case 0: {
		break;
		}
	}

}

void mcmc_omega_2(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	double omega_2_new = -9999;

	while(omega_2_new<0){
		omega_2_new = Rcpp::as<double>(para_scalar["omega_2"])+ 0.1*R::rnorm(0.0,1.0);
	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["omega_2"] = omega_2_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);

// Rcout <<  log_lh_current << "," << log_lh_new << endl;

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_scalar["omega_2"] = omega_2_new ;
		break;
		}
		case 0: {
		break;
		}
	}

}


void mcmc_omega_3(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat,NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	double omega_3_new = -9999;

	while(omega_3_new<0){
		omega_3_new = Rcpp::as<double>(para_scalar["omega_3"])+ 0.1*R::rnorm(0.0,1.0);
	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["omega_3"] = omega_3_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);

// Rcout <<  log_lh_current << "," << log_lh_new << endl;

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_scalar["omega_3"] = omega_3_new ;
		break;
		}
		case 0: {
		break;
		}
	}

}

void mcmc_t_intervention(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat,NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	// IntegerVector steps = Rcpp::seq(-2,2);

	double t_intervention_new = -9999;

	// while(t_intervention_new<0 | t_intervention_new> Rcpp::as<double>(para_scalar["t_max"])){
	while(t_intervention_new<0){

		// IntegerVector rw_vec = Rcpp::sample(steps, 1, false, R_NilValue);
		// t_intervention_new = Rcpp::as<double>(para_scalar["t_intervention"])+ rw_vec[0];

		t_intervention_new = Rcpp::as<double>(para_scalar["t_intervention"])+ 1*R::rnorm(0.0,1.0);

	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["t_intervention"] = t_intervention_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat,contact_mat_A, contact_mat_B, dir_results);

// Rcout <<  log_lh_current << "," << log_lh_new << endl;

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_scalar["t_intervention"] = t_intervention_new ;
// para_scalar["t_mask"] = t_intervention_new ;

		break;
		}
		case 0: {
		break;
		}
	}

}

void mcmc_t_intervention_2(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	// IntegerVector steps = Rcpp::seq(-2,2);

	double t_intervention_1 = Rcpp::as<double>(para_scalar["t_intervention"]);
	double t_intervention_3 = Rcpp::as<double>(para_scalar["t_intervention_3"]);

	double t_max = Rcpp::as<double>(para_scalar["t_max"]);

	double t_intervention_new = -9999;

	// while(t_intervention_new<0 | t_intervention_new> Rcpp::as<double>(para_scalar["t_max"])){
	while(t_intervention_new<=t_intervention_1 | t_intervention_new>=t_intervention_3 | t_intervention_new>t_max){

		// IntegerVector rw_vec = Rcpp::sample(steps, 1, false, R_NilValue);
		// t_intervention_new = Rcpp::as<double>(para_scalar["t_intervention"])+ rw_vec[0];

		t_intervention_new = Rcpp::as<double>(para_scalar["t_intervention_2"])+ 1*R::rnorm(0.0,1.0);

	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["t_intervention_2"] = t_intervention_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);

// Rcout <<  log_lh_current << "," << log_lh_new << endl;

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_scalar["t_intervention_2"] = t_intervention_new ;
		break;
		}
		case 0: {
		break;
		}
	}

}


void mcmc_t_intervention_3(Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B,std::string dir_results, double& log_lh_current){


	// IntegerVector steps = Rcpp::seq(-2,2);

	double t_intervention_1 = Rcpp::as<double>(para_scalar["t_intervention"]);
	double t_intervention_2 = Rcpp::as<double>(para_scalar["t_intervention_2"]);

	double t_max = Rcpp::as<double>(para_scalar["t_max"]);

	double t_intervention_new = -9999;

	while(t_intervention_new<=t_intervention_2 | t_intervention_new>t_max){

		// IntegerVector rw_vec = Rcpp::sample(steps, 1, false, R_NilValue);
		// t_intervention_new = Rcpp::as<double>(para_scalar["t_intervention"])+ rw_vec[0];

		t_intervention_new = Rcpp::as<double>(para_scalar["t_intervention_3"])+ 1*R::rnorm(0.0,1.0);

	}

	Rcpp::List para_scalar_new = clone(para_scalar);
	para_scalar_new["t_intervention_3"] = t_intervention_new;

	double log_lh_new =  llh_func(epi_list, sero_list,  para_scalar_new, para_vec_age, para_vec_sero,  mt_vec, contact_mat,contact_mat_A, contact_mat_B, dir_results);

// Rcout <<  log_lh_current << "," << log_lh_new << endl;

	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current));

	double uniform_rv = R::runif(0.0,1.0);

	switch(uniform_rv<=acp_pr){
		case 1: {
			log_lh_current = log_lh_new;
			para_scalar["t_intervention_3"] = t_intervention_new ;
		break;
		}
		case 0: {
		break;
		}
	}

}




void mcmc_n_SE(int i_age, int t, Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat, NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B,std::string dir_results, double& log_lh_current){
// update n_SE (and other downstream variables), at one age and time //

	// double log_lh_new = (log_lh_current);

	int t_max =  Rcpp::as<int>(para_scalar["t_max"]);

	int n_age =  Rcpp::as<int>(para_scalar["n_age"]); // number of age gps

	int EU_p1 = Rcpp::as<int>(para_scalar["EU_p1"]); // d_EU
	int EI_p1 = Rcpp::as<int>(para_scalar["EI_p1"]); // d_EI
	int UR_p1 = Rcpp::as<int>(para_scalar["UR_p1"]); // d_UR
	int IR_p1 = Rcpp::as<int>(para_scalar["IR_p1"]); //d_IR


	NumericVector N_age =  para_vec_age["N_age"]; // pop in each age gp

	NumericVector beta_sym_vec = para_vec_age["beta_sym_vec"];
	NumericVector gamma_vec = para_vec_age["gamma_vec"];

	NumericVector t_sero = para_vec_sero["t_sero"];
	NumericVector n_sero = para_vec_sero["n_sero"];

	//---

	// epi data of a particular age gp: columns age_gp, time, S, E, I, U, R, n_SE, n_EU, n_EI, n_UR, n_IR //
	Rcpp::NumericMatrix epi = Rcpp::as<Rcpp::NumericMatrix>(epi_list[i_age]);

	// double gamma = gamma_vec[i_age];

	NumericVector S = epi(_,2);
	NumericVector E = epi(_,3);
	NumericVector I = epi(_,4);
	NumericVector U = epi(_,5);
	NumericVector R = epi(_,6);

	NumericVector n_SE = epi(_,7);
	NumericVector n_EU = epi(_,8);
	NumericVector n_EI = epi(_,9);
	NumericVector n_UR = epi(_,10);
	NumericVector n_IR = epi(_,11);

	NumericVector S_new = clone(S);
	NumericVector E_new = clone(E);
	NumericVector I_new = clone(I);
	NumericVector U_new = clone(U);
	NumericVector R_new = clone(R);

	NumericVector n_SE_new = clone(n_SE);
	NumericVector n_EU_new = clone(n_EU);
	NumericVector n_EI_new = clone(n_EI);
	NumericVector n_UR_new = clone(n_UR);
	NumericVector n_IR_new = clone(n_IR);



	double low_n_SE_t = 0.0;
	if ((t+EI_p1)<=(t_max-1)){
		low_n_SE_t = n_EI_new[t+EI_p1]; // note that n_EI is assumed not chaning 
	}

	double up_n_SE_t = S_new[t-1]; // cant really use this may be changing due to updates of n_SE at other time t

	IntegerVector times = Rcpp::seq(-10,10);

	double n_SE_t = -999.0;


		while(n_SE_t<low_n_SE_t | n_SE_t>up_n_SE_t){

	 		IntegerVector rw = Rcpp::sample(times, 1, false, R_NilValue);
			n_SE_t = n_SE[t] + rw[0];

	 	}




	//----


	n_SE_new[t] = n_SE_t;

	//-- change the n_EU at (t+d_EU) according to the newly drawn n_SE_t ---//

	if((t+EU_p1)<=(t_max-1)) {

		double n_EU_t_plus_d_EU = n_SE_new[t] - n_EI_new[t+EI_p1] ; 

		n_EU_new[t+EU_p1] = n_EU_t_plus_d_EU; // update n_EU at t+d_EU accordingly
	}

	//-- update n_UR(t+d_EU+d_UR) --//
	if((t+EU_p1+UR_p1)<=(t_max-1)) {
		n_UR_new[t+EU_p1+UR_p1] = n_EU_new[t+EU_p1];
	}

	//-- update S(t+1), S(t+2), .. S(t_max). --//

	for (int t_ahead=0;t_ahead<(t_max-1);t_ahead++){
			S_new[t_ahead+1] = S_new[t_ahead] - n_SE_new[t_ahead+1]; // note that in this sampler, S(0) will not be change as new n_SE(0) does not affect S(0)

			if(S_new[t_ahead+1]<0) {
				S_new[t_ahead+1] = 0.0;
			}

// if (S_new[t_ahead+1]<0) Rcout << "negative S_new[t_ahead+1]" << "," << S_new[t_ahead+1] << "," <<  S_new[t_ahead]<<  "," << n_SE_new[t_ahead] << ","  << i_age << "," << t << "," <<  t_ahead << endl;


	}


	//-- update E(t), E(t+1), ... --//

	for (int t_ahead=1;t_ahead<=(t_max-1);t_ahead++){

		// double E_prev = 3.0;
		// if ((t_ahead-1)>=0) {
		double E_prev = E_new[t_ahead-1];
		// }

		E_new[t_ahead] =  E_prev + n_SE_new[t_ahead] - n_EI_new[t_ahead] - n_EU_new[t_ahead] ;


	}


	//- also update U(t+d_EU+1), U(t+d_EU+2), .... --//
	for (int t_ahead=1;t_ahead<=(t_max-1);t_ahead++){
			U_new[t_ahead] = U_new[t_ahead-1] + n_EU_new[t_ahead] - n_UR_new[t_ahead];

// if (U_new[t_ahead]<0) Rcout << "negative U_new[t_ahead]" << "," << U_new[t_ahead] <<  "," << t_ahead << endl;

	}	



	//- also update R(t+d_EI+d_IR+1), R(t+d_EI+d_IR+2), .... --//
	for (int t_ahead=1;t_ahead<=(t_max-1);t_ahead++){
			R_new[t_ahead] = R_new[t_ahead-1] + n_IR_new[t_ahead] + n_UR_new[t_ahead];

// if (R_new[t_ahead]<0) Rcout << "negative R_new[t_ahead]" << "," << R_new[t_ahead] << endl;

	}

	//-- accept or reject ---//

	Rcpp::List epi_list_new = clone(epi_list);
	Rcpp::NumericMatrix epi_new = clone(epi);

	epi_new(_,2) = S_new; // this works as updating the epi_list
	epi_new(_,3) = E_new; // 
	epi_new(_,4) = I_new; // 
	epi_new(_,5) = U_new; // this works as updating the epi_list
	epi_new(_,6) = R_new; // this works as updating the epi_list
	epi_new(_,7) = n_SE_new; // this works as updating the epi_list
	epi_new(_,8) = n_EU_new; // this works as updating the epi_list
	epi_new(_,9) = n_EI_new; // this works as updating the epi_list
	epi_new(_,10) = n_UR_new; // this works as updating the epi_list
	epi_new(_,11) = n_IR_new; // this works as updating the epi_list

	epi_list_new[i_age] = epi_new;





	double log_proposal_diff = 0.0;
	double  log_proposal_diff_1 = 0.0;
	double log_proposal_diff_2 = 0.0;

	if ((t+EU_p1)<=(t_max-1)){


		log_proposal_diff = log_proposal_diff_1 - log_proposal_diff_2;



	}




	double log_lh_new =  llh_func(epi_list_new, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current + log_proposal_diff));

// Rcout  <<  acp_pr << "," << log_lh_current << "," << log_lh_new  << endl;

// Rcout  << i_age << "," << t << "," << log_lh_current << "," << log_lh_new <<  "," << log_partial_old  << "," << log_partial_new << endl;

// Rcout << exp(log_lh_new - log_lh_current + log_proposal_diff)<< "," << acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1) & (arma::is_finite(log_proposal_diff)==1) ){

	// acp_pr = 1.0;

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				epi_list[i_age] = epi_new ;
				// epi_list = epi_list_new;

	// Rcout << "yes, a change!" << endl;

			break;
			}
			case 0: {
				// log_lh_current = log_lh_current;
				// epi_list = epi_list ;
			break;
			}
		}

		
	}



}





void mcmc_R0(int i_age, Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat,NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){
// update R(0) //

	// double log_lh_new = (log_lh_current);

	int t_max =  Rcpp::as<int>(para_scalar["t_max"]);

	int n_age =  Rcpp::as<int>(para_scalar["n_age"]); // number of age gps

	int EU_p1 = Rcpp::as<int>(para_scalar["EU_p1"]); // d_EU
	int EI_p1 = Rcpp::as<int>(para_scalar["EI_p1"]); // d_EI
	int UR_p1 = Rcpp::as<int>(para_scalar["UR_p1"]); // d_UR
	int IR_p1 = Rcpp::as<int>(para_scalar["IR_p1"]); //d_IR


	NumericVector N_age =  para_vec_age["N_age"]; // pop in each age gp

	NumericVector beta_sym_vec = para_vec_age["beta_sym_vec"];
	NumericVector gamma_vec = para_vec_age["gamma_vec"];

	NumericVector t_sero = para_vec_sero["t_sero"];
	NumericVector n_sero = para_vec_sero["n_sero"];

	//---

	// epi data of a particular age gp: columns age_gp, time, S, E, I, U, R, n_SE, n_EU, n_EI, n_UR, n_IR //
	Rcpp::NumericMatrix epi = Rcpp::as<Rcpp::NumericMatrix>(epi_list[i_age]);

	// double gamma = gamma_vec[i_age];

	NumericVector S = epi(_,2);
	NumericVector E = epi(_,3);
	NumericVector I = epi(_,4);
	NumericVector U = epi(_,5);
	NumericVector R = epi(_,6);

	NumericVector n_SE = epi(_,7);
	NumericVector n_EU = epi(_,8);
	NumericVector n_EI = epi(_,9);
	NumericVector n_UR = epi(_,10);
	NumericVector n_IR = epi(_,11);

	NumericVector S_new = clone(S);
	NumericVector E_new = clone(E);
	NumericVector I_new = clone(I);
	NumericVector U_new = clone(U);
	NumericVector R_new = clone(R);

	NumericVector n_SE_new = clone(n_SE);
	NumericVector n_EU_new = clone(n_EU);
	NumericVector n_EI_new = clone(n_EI);
	NumericVector n_UR_new = clone(n_UR);
	NumericVector n_IR_new = clone(n_IR);


	double low_R_0 = 0.0;
	
	// double up_R_0 = N_age[i_age];
	double up_R_0 = N_age[i_age] - (S[0]+E[0]+ I[0]+U[0]) ;



	NumericVector cul_rate_0 = {0.023,0.023,0.0285,0.022}; // at March 25 2020

	NumericVector cul_rate_t_max = {0.092,0.092,0.113,0.094}; // at June 30, t_max

	double up_R_0_cul_cases = N_age[i_age]*cul_rate_0[i_age];
	double up_R_t_max_cul_cases = N_age[i_age]*cul_rate_t_max[i_age];

	IntegerVector times = Rcpp::seq(-10,10);

	double new_R_0 = -999.0;
	while(new_R_0<low_R_0 | new_R_0>up_R_0 | R_new[t_max-1]>up_R_t_max_cul_cases | S_new[0]<0){

 		IntegerVector rw = Rcpp::sample(times, 1, false, R_NilValue);
		new_R_0 = R[0] + rw[0];

		double old_R_0 = R[0];

		R_new[0] = new_R_0;


		//- update R(t) --//
		for (int t_ahead=1;t_ahead<=(t_max-1);t_ahead++){
				R_new[t_ahead] = R_new[t_ahead-1] + n_UR_new[t_ahead] + n_IR_new[t_ahead];
		}


		//--update S(0)--//
		S_new[0] = S[0] - (new_R_0 - old_R_0);

	}

	//--update S(t)--//
	for (int t_ahead=0;t_ahead<(t_max-1);t_ahead++){
			S_new[t_ahead+1] = S_new[t_ahead] - n_SE_new[t_ahead+1];

	}




	//-- accept or reject ---//

	Rcpp::List epi_list_new = clone(epi_list);
	Rcpp::NumericMatrix epi_new = clone(epi);

	epi_new(_,2) = S_new; // this works as updating the epi_list
	epi_new(_,3) = E_new; // 
	epi_new(_,4) = I_new; // 
	epi_new(_,5) = U_new; // this works as updating the epi_list
	epi_new(_,6) = R_new; // this works as updating the epi_list
	epi_new(_,7) = n_SE_new; // this works as updating the epi_list
	epi_new(_,8) = n_EU_new; // this works as updating the epi_list
	epi_new(_,9) = n_EI_new; // this works as updating the epi_list
	epi_new(_,10) = n_UR_new; // this works as updating the epi_list
	epi_new(_,11) = n_IR_new; // this works as updating the epi_list

	epi_list_new[i_age] = epi_new;




	double log_proposal_diff = 0.0;
	double  log_proposal_diff_1 = 0.0;
	double log_proposal_diff_2 = 0.0;


	log_proposal_diff = log_proposal_diff_1 - log_proposal_diff_2;


	double log_lh_new =  llh_func(epi_list_new, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B, dir_results);


	double acp_pr = min(1.0,exp(log_lh_new - log_lh_current + log_proposal_diff));

// Rcout << R[0] << "," << R_new[0] << "," << acp_pr << endl;


// Rcout  <<  acp_pr << "," << log_lh_current << "," << log_lh_new  << endl;

// Rcout  << i_age << "," << t << "," << log_lh_current << "," << log_lh_new <<  "," << log_partial_old  << "," << log_partial_new << endl;

// Rcout << exp(log_lh_new - log_lh_current + log_proposal_diff)<< "," << acp_pr << endl;

	if ((arma::is_finite(log_lh_new)==1) & (arma::is_finite(log_proposal_diff)==1) ){

	// acp_pr = 1.0;

		double uniform_rv = R::runif(0.0,1.0);

		switch(uniform_rv<=acp_pr){
			case 1: {
				log_lh_current = log_lh_new;
				epi_list[i_age] = epi_new ;
				// epi_list = epi_list_new;

			break;
			}
			case 0: {
				// log_lh_current = log_lh_current;
				// epi_list = epi_list ;
			break;
			}
		}

		
	}



}





// [[Rcpp::export]]
void mcmc_HH(int& n_iter, Rcpp::List& epi_list, Rcpp::List& sero_list,  Rcpp::List& para_scalar, Rcpp::List& para_vec_age, Rcpp::List& para_vec_sero,  NumericVector& mt_vec, NumericMatrix& contact_mat,NumericMatrix& contact_mat_A,NumericMatrix& contact_mat_B, std::string dir_results, double& log_lh_current){


	int t_max =  Rcpp::as<int>(para_scalar["t_max"]);
	int n_age =  Rcpp::as<int>(para_scalar["n_age"]); // number of age gps


	int EU_p1 = Rcpp::as<int>(para_scalar["EU_p1"]); // d_EU

 	// IntegerVector times = Rcpp::seq(1,t_max-1);
 	IntegerVector times = Rcpp::seq(1,t_max-EU_p1-1);

 	int n_time_to_update = 1;

	NumericVector beta_sym_vec = para_vec_age["beta_sym_vec"];
	NumericVector gamma_vec = para_vec_age["gamma_vec"];

	std::ofstream myfile_para;
	std::ofstream myfile_lh;
	std::ofstream myfile_epi_0;
	std::ofstream myfile_epi_1;
	std::ofstream myfile_epi_2;
	std::ofstream myfile_epi_3;
	// std::ofstream myfile_epi_4;

	myfile_para.open(dir_results.c_str() + std::string("para_current.txt"));
	myfile_lh.open(dir_results.c_str() + std::string("log_lh_current.txt"));

	myfile_epi_0.open(dir_results.c_str() + std::string("epi_0_current.txt"));
	myfile_epi_1.open(dir_results.c_str() + std::string("epi_1_current.txt"));
	myfile_epi_2.open(dir_results.c_str() + std::string("epi_2_current.txt"));
	myfile_epi_3.open(dir_results.c_str() + std::string("epi_3_current.txt"));


	myfile_para << "a" <<  "," << "b" << "," << "a_2" <<  "," << "b_2" << "," << "a_3" <<  "," << "b_3" << "," << "a_4" <<  "," << "b_4" << "," << "prop_asym" << "," << "beta_1" << "," << "beta_2" <<  "," <<  "beta_3" << "," << "beta_4" << "," << "gamma_1" << "," << "gamma_2" << "," << "gamma_3" << "," << "gamma_4" << "," << "gamma_mask_1" << "," << "gamma_mask_2" << "," << "gamma_mask_3" << "," << "gamma_mask_4" << "," << "omega_1" <<"," << "omega_2" <<"," << "omega_3" << "," << "t_intervention" << "," << "t_intervention_2" << "," << "t_intervention_3" << "\n";


	myfile_epi_0 << "i_iter_epi" << "," <<  "age_gp" << "," << "time "<< "," << "S" << "," << "E"<< ","  << "I" << "," << "U" << "," << "R" << ","  << "n_SE" << "," << "n_EU" << "," << "n_EI" << "," << "n_UR" << "," << "n_IR" << "\n";

	myfile_epi_1 << "i_iter_epi" << "," <<  "age_gp" << "," << "time "<< "," << "S" << "," << "E"<< ","  << "I" << "," << "U" << "," << "R" << ","  << "n_SE" << "," << "n_EU" << "," << "n_EI" << "," << "n_UR" << "," << "n_IR" << "\n";

	myfile_epi_2 << "i_iter_epi" << "," <<  "age_gp" << "," << "time "<< "," << "S" << "," << "E"<< ","  << "I" << "," << "U" << "," << "R" << ","  << "n_SE" << "," << "n_EU" << "," << "n_EI" << "," << "n_UR" << "," << "n_IR" << "\n";

	myfile_epi_3 << "i_iter_epi" << "," <<  "age_gp" << "," << "time "<< "," << "S" << "," << "E"<< ","  << "I" << "," << "U" << "," << "R" << ","  << "n_SE" << "," << "n_EU" << "," << "n_EI" << "," << "n_UR" << "," << "n_IR" << "\n";

	int i_iter_epi = 0;

	for (int i_iter=0;i_iter<n_iter; i_iter++){


 		IntegerVector t_to_update = Rcpp::sample(times, n_time_to_update, false, R_NilValue);

		mcmc_gammas(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // estimate susceptibilities before lifting the lockdown

		mcmc_gammas_mask_full(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // estimate susceptibilities after lifting the lockdown

		mcmc_beta(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // estimate population level transmissibility before first change point 

	
		mcmc_b_a(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // parameters governing reporting rate for 1st age gp


		mcmc_b_2_a_2(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); //.. 2nd ..

		mcmc_b_3_a_3(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // ..3rd ...

		mcmc_b_4_a_4(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // ..4th..

		mcmc_omega_1(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // parameter governs the transmissibility between 1st and 2nd change points
		mcmc_omega_2(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); //.. after second change point 

		mcmc_t_intervention(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current);  // estimate first change point
		mcmc_t_intervention_2(epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat,contact_mat_A, contact_mat_B, dir_results, log_lh_current); //.. 2nd...


		for (int i_age=0;i_age<=(n_age-1); i_age++){
		for (int k=0;k<=(n_time_to_update-1); k++){
			mcmc_n_SE(i_age, t_to_update[k], epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // impute number of infections at t
		}
		}



		for (int i_age=0;i_age<=(n_age-1); i_age++){
		 	mcmc_R0(i_age,  epi_list, sero_list,  para_scalar, para_vec_age, para_vec_sero,  mt_vec, contact_mat, contact_mat_A, contact_mat_B,dir_results, log_lh_current); // update R(0) and R(t)

		}
		

		myfile_para << Rcpp::as<double>(para_scalar["a"]) <<  "," <<  Rcpp::as<double>(para_scalar["b"]) << ","<< Rcpp::as<double>(para_scalar["a_2"]) <<  "," <<  Rcpp::as<double>(para_scalar["b_2"]) <<  "," << Rcpp::as<double>(para_scalar["a_3"]) <<  "," <<  Rcpp::as<double>(para_scalar["b_3"]) <<   ","<< Rcpp::as<double>(para_scalar["a_4"]) <<  "," <<  Rcpp::as<double>(para_scalar["b_4"]) <<  ","  << Rcpp::as<double>(para_scalar["prop_asym"])  << "," <<  Rcpp::as<NumericVector>(para_vec_age["beta_sym_vec"])[0] <<"," <<  Rcpp::as<NumericVector>(para_vec_age["beta_sym_vec"])[1] << "," <<  Rcpp::as<NumericVector>(para_vec_age["beta_sym_vec"])[2] << "," <<  Rcpp::as<NumericVector>(para_vec_age["beta_sym_vec"])[3]  << "," <<  Rcpp::as<NumericVector>(para_vec_age["gamma_vec"])[0] <<"," << Rcpp::as<NumericVector>(para_vec_age["gamma_vec"])[1] << "," <<  Rcpp::as<NumericVector>(para_vec_age["gamma_vec"])[2]  << "," << Rcpp::as<NumericVector>(para_vec_age["gamma_vec"])[3] <<  "," <<  Rcpp::as<NumericVector>(para_vec_age["gamma_vec_2"])[0] <<"," << Rcpp::as<NumericVector>(para_vec_age["gamma_vec_2"])[1] << "," <<  Rcpp::as<NumericVector>(para_vec_age["gamma_vec_2"])[2]  << "," << Rcpp::as<NumericVector>(para_vec_age["gamma_vec_2"])[3] << "," << Rcpp::as<double>(para_scalar["omega_1"]) << "," << Rcpp::as<double>(para_scalar["omega_2"]) << "," << Rcpp::as<double>(para_scalar["omega_3"]) << "," << Rcpp::as<double>(para_scalar["t_intervention"])  << "," << Rcpp::as<double>(para_scalar["t_intervention_2"])  << "," << Rcpp::as<double>(para_scalar["t_intervention_3"])  << "\n";


		myfile_lh << log_lh_current << "\n";



		if((i_iter%1000)==0){


			i_iter_epi = i_iter_epi + 1;


			for (int i_age=0;i_age<=(n_age-1); i_age++){
				Rcpp::NumericMatrix epi = epi_list[i_age];


				NumericVector age_gp = epi(_,0);
				NumericVector time = epi(_,1);

				NumericVector S = epi(_,2);
				NumericVector E = epi(_,3);
				NumericVector I = epi(_,4);
				NumericVector U = epi(_,5);
				NumericVector R = epi(_,6);

				NumericVector n_SE = epi(_,7);
				NumericVector n_EU = epi(_,8);
				NumericVector n_EI = epi(_,9);
				NumericVector n_UR = epi(_,10);
				NumericVector n_IR = epi(_,11);



				for (int t=0;t<=(t_max-1); t++){

					if (i_age==0) {
						myfile_epi_0 << i_iter_epi << "," <<  age_gp[t] << "," << time[t] << "," << S[t] << "," << E[t]  << ","  << I[t] << "," << U[t] << "," << R[t] << ","  << n_SE[t] << "," << n_EU[t] << "," << n_EI[t] << "," << n_UR[t] << "," << n_IR[t] << "\n";
					}

					if (i_age==1) {
						myfile_epi_1 << i_iter_epi << "," << age_gp[t] << "," << time[t] << "," << S[t] << "," << E[t]  << ","  << I[t] << "," << U[t] << "," << R[t] << ","  << n_SE[t] << "," << n_EU[t] << "," << n_EI[t] << "," << n_UR[t] << "," << n_IR[t] << "\n";
					}


					if (i_age==2) {
						myfile_epi_2 << i_iter_epi << "," << age_gp[t] << "," << time[t] << "," << S[t] << "," << E[t]  << ","  << I[t] << "," << U[t] << "," << R[t] << ","  << n_SE[t] << "," << n_EU[t] << "," << n_EI[t] << "," << n_UR[t] << "," << n_IR[t] << "\n";
					}

					if (i_age==3) {
						myfile_epi_3 << i_iter_epi << "," << age_gp[t] << "," << time[t] << "," << S[t] << "," << E[t]  << ","  << I[t] << "," << U[t] << "," << R[t] << ","  << n_SE[t] << "," << n_EU[t] << "," << n_EI[t] << "," << n_UR[t] << "," << n_IR[t] << "\n";
					}


				}

			}
		}


	}

	myfile_para.close();
	myfile_lh.close();
	myfile_epi_0.close();
	myfile_epi_1.close();
	myfile_epi_2.close();
	myfile_epi_3.close();

}






