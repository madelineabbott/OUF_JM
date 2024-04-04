
// This code builds on code originally written by Trung Dung Tran; this code
// is available on the author's Github: https://github.com/tdt01/LOUmodels.

functions {
  // adapted from: https://groups.google.com/g/stan-users/c/quCUh0sHbgA?pli=1
  vector r_style_subset(matrix x, int col_id, int[] row_ids) {
    vector[size(row_ids)] y;
    vector[rows(x)] x_col;
    int pos;
    pos = 1;
    x_col = col(x, col_id);
    for (i in 1:size(row_ids)) {
      y[pos] = x_col[row_ids[i]];
      pos = pos + 1;
    }
    return y;
  }
}

data {
	int Nall; // Sample size (meas occ + grid)
	int Nlong; // Sample size (meas occ only)
	int N; // Number of subjects

	int K; // Number of measured longitudinal outcomes
	int P; // Number of latent factors

  // For longitudinal measurement occasions
	int meas_occ_rows[Nlong];
	int cumu_meas[N]; 
	int repme_meas[N];

  // For longitudinal measurement occasions + grid points
	int cumu[N];  
	int repme[N]; 

	matrix[Nlong, K] Y; // Measured longitudinal outcome

	vector[Nall] deltat; // Time gaps between meas occ + grid points
	
	int status[N]; // Survival outcome status indicator
}

parameters {

  // Measurement submodel parameters
	vector<lower=0.000001>[K] lambda; 
	real<lower=0.000001> sigma_lambda;
	vector<lower=0.000001>[K] sigma_u;
	vector<lower=0.000001>[K] sigma_e;

  // Latent factors (OU process)
	matrix[Nall, P] eta;
	
	// Structural submodel parameters (OU process)
	matrix[P, P] theta_ou;
	real<lower=-0.999999, upper=0.999999> rho;
	
	// Survival submodel parameters (assuming simple constant baseline haz.)
	real beta_0; // baseline hazard
	real beta_1; // for eta1
	real beta_2; // for eta2
}

transformed parameters {

  // constraint on theta_ou's eigenvalues (from Biometrics 2021 by Tran et al)
	real<lower=0.000001> constraint1;
	real<lower=0.000001> constraint2;

	constraint1 = theta_ou[1, 1] + theta_ou[2, 2];
	constraint2 = theta_ou[1, 1] * theta_ou[2, 2] - theta_ou[1, 2] * theta_ou[2, 1];

}

model{
  
  vector[Nall] haz; // Hazard function at ALL times (meas occ + grid pts + events)
  vector[N] cumul_haz; // Cumulative hazard function at event times
  vector[N] haz_etimes; // Hazard function at event times
	
	matrix[P, P] V; // Stationary correlation matrix for OU process
	
	// Priors
	lambda ~ normal(1, sigma_lambda);
	sigma_lambda ~ cauchy(0, 5);
	sigma_u ~ cauchy(0, 5); 
	sigma_e ~ cauchy(0, 5);
	to_vector(theta_ou) ~ normal(0, 10);
	rho ~ uniform(-0.999999, 0.999999);
	beta_0 ~ normal(0, 5);
	beta_1 ~ normal(0, 5);
	beta_2 ~ normal(0, 5);

  // Structural submodel part
  // j = 1 (First time point)
	V = [[1, rho], [rho, 1]];
	for (i in 1 : N){
		int k;
		k = cumu[i] - repme[i] + 1;
		eta[k] ~ multi_normal([0, 0]', V);
	}
	// Second time points and later
	for (i in 1 : N){
	  if (repme[i] > 1){
	    for (j in 2 : repme[i]){
			  int k;
			  vector[P] cond_mean;
			  matrix[P, P] cond_covar;
			  k = cumu[i] - repme[i] + j;
			  cond_mean = matrix_exp(-deltat[k] * theta_ou) * to_vector(eta[k-1]);
			  cond_covar = V - matrix_exp(-deltat[k] * theta_ou) * V * matrix_exp(-deltat[k] * theta_ou');
			  eta[k] ~ multi_normal(cond_mean, cond_covar);
		  }
	  }
	}

  // Measurement submodel part
	for (i in 1 : N){
	  int k;
	  int k_meas_start;
	  int k_meas_stop;
	  int ni = repme_meas[i];
	  
	  vector[ni] eta1;
	  vector[ni] eta2;
	  
	  matrix[ni, ni] Y1_var;
	  matrix[ni, ni] Y2_var;
	  matrix[ni, ni] Y3_var;
	  matrix[ni, ni] Y4_var;
	  
	  k_meas_start = cumu_meas[i] - repme_meas[i] + 1;
	  k_meas_stop = cumu_meas[i];
	  
	  eta1 = r_style_subset(eta, 1, meas_occ_rows[k_meas_start:k_meas_stop]);
	  eta2 = r_style_subset(eta, 2, meas_occ_rows[k_meas_start:k_meas_stop]);
	  
	  Y1_var = add_diag(rep_matrix(sigma_u[1], ni, ni), sigma_e[1]);
	  Y2_var = add_diag(rep_matrix(sigma_u[2], ni, ni), sigma_e[2]);
	  Y3_var = add_diag(rep_matrix(sigma_u[3], ni, ni), sigma_e[3]);
	  Y4_var = add_diag(rep_matrix(sigma_u[4], ni, ni), sigma_e[4]);
	  
    Y[k_meas_start:k_meas_stop, 1] ~ multi_normal(lambda[1]*eta1, Y1_var);
    Y[k_meas_start:k_meas_stop, 2] ~ multi_normal(lambda[2]*eta1, Y2_var);
    Y[k_meas_start:k_meas_stop, 3] ~ multi_normal(lambda[3]*eta2, Y3_var);
    Y[k_meas_start:k_meas_stop, 4] ~ multi_normal(lambda[4]*eta2, Y4_var);
	}
  
	
  // Survival submodel part
  for (j in 1 : Nall){ // Calculate haz across all pts (meas occ + grid pts)
    haz[j] = exp(beta_0 + beta_1 * eta[j, 1] + beta_2 * eta[j, 2]);
  }

  for (i in 1 : N){
    int k;
    cumul_haz[i] = 0;
    if (repme[i] > 1){
      for (j in 2 : repme[i]){
        k = cumu[i] - repme[i] + j;
        cumul_haz[i] = cumul_haz[i] + (haz[k] + haz[k-1])/2 * deltat[k];
      }
      haz_etimes[i] = haz[k]; // Hazard at the event time
    }else {
      k = cumu[i] - repme[i] + 1;
      haz_etimes[i] = haz[k]; // Hazard at the event time
    }

  }

  for (i in 1 : N){
    target += status[i] * log(haz_etimes[i]) - cumul_haz[i];
  }

}

