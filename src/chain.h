#ifndef _CHAIN_H
#define _CHAIN_H

double *alpha_r_t_UO(int nt, double dt,double* para);
double *alpha_r_t_expcos(int nt, double dt, double *para);
double spectral_UO(int n_para, double *para, double o);
double spectral_OM(int npara, double *para, double omega);
double spectral_test(int n_para, double *para, double o);
double *gamma_t_from_spectral_test(int N, double dt, int n_para, double *para, double (*spect)(int, double*, double));
double *gamma_t_from_spectral(int nt, double dt, double tc, int n_para, double *para, double (*spect)(int, double*, double));
double *gamma_t_from_spectral_truncate(int nt, double dt, double tc, int n_para, double *para, double (*spect)(int, double*, double));
int  draw_spectral(char *name, double *para,double (*spect)(int, double*, double));
double *alpha_t_from_gamma_t(int nt, double dt, int n_para, double *para, double *gamma_t);


extern double T_l;//2.0;//Kelvin;//The temperature of left bath.
extern double T_r;//1.0;//Kelvin;//The temperature of right bath.
extern double tau_l;//fs//The
extern double tau_r;//fs
extern double epselon_l;//fs^-1
extern double epselon_r;//fs^-1
extern double mass;//12.0;//amu
extern double k_b;
extern int num_of_atoms;
extern double x_eq;//pm
extern double D;//0.5*0.2700;//3.678e+02/(0.01*0.01);//unit of energy
extern double a;//1.0;//1.875e-2*(0.01);//pm^-1
//double k_force_constant=2*D*a^2=0.2700;//only for harmonic case
extern double *memory_of_p1;
extern double *memory_of_pN;
extern int N_memory_step_left;
extern int N_memory_step_right;
extern double gamma_para_b_l;
extern double gamma_para_b_r;
int chain(double *dt_of_rlt, int *nt_of_rlt,
double** rlt_mpi, int *good_orb, int *nobvs,
char * output_file_name, int argc, char* argv[]);
/**
 * This function is the same as chain(...), except the response function \gamma(t) is determined by spectral function \hat{\gamma}(\omega).
 * @param dt_of_rlt
 * @param nt_of_rlt
 * @param rlt_mpi
 * @param good_orb
 * @param nobvs
 * @param output_file_name
 * @return
 */
int chain_spectral(double *dt_of_rlt,
int *nt_of_rlt, double** rlt_mpi,
int *good_orb, int *nobvs, char * output_file_name,
int argc, char* argv[]);
/**
 * This function is the same as chain(..), except the response function \gamma(t)=c*exp(-a*|t|)*cos(b*t)
 * @param dt_of_rlt
 * @param nt_of_rlt
 * @param rlt_mpi
 * @param good_orb
 * @param n_obvs
 * @param output_file_name
 * @return
 */
int chain_expcos(double *dt_of_rlt, int *nt_of_rlt,
		 double** rlt_mpi, int *good_orb, int *n_obvs, char* output_file_name,
		 int argc, char *argv[]);

int chain_ana(double *dt_of_rlt, int *nt_of_rlt,
		double** rlt_mpi, int *good_orb, int *n_obvs, char* output_file_name, int argc, char *argv[]);

extern double *corr_rand_gener_2bath_mem;
int corr_rand_gener_2bath_malloc(int nt);
int corr_rand_gener_2bath(const int nt, double** rand, double *lamda_r, double *lamda_l, int flag);
int corr_rand_gener_2bath_free();
int white_rand_gener_2bath(const int nt, double dt, double T_right, double T_left,  double** rand, int flag);

/**
 * chain driven by uo noises
 * @param ny number of coordinates that is propagated
 * @param x time
 * @param rdm random force
 * @param y pointer to the array of coordinates that is propagated
 * @param dydx dy/dx calculated from the parameters above with potential function
 */
void propagator_derive_chain_UO_mors(int ny, double x, double *rdm,
			     double* y, double* dydx);
void propagator_derive_chain_expcos_mors(int ny, double x, double *rdm,
					 double* y, double* dydx);
void propagator_derive_chain_UO_harm(int ny, double x, double *rdm,
				 double* y, double* dydx);
void propagator_derive_chain_expcos_harm(int ny, double x, double *rdm,
					 double* y, double* dydx);
void propagator_derive_chain_white_harm(int ny, double x, double *rdm,
					double* y, double* dydx);
void propagator_derive_chain_UO_ana_harm(int ny, double x, double *rdm,
				     double* y, double* dydx);
void propagator_derive_chain_UO_ana_mors(int ny, double x, double *rdm,
			     double* y, double* dydx);
void propagator_derive_chain_general_spectral(int ny, double x, double *rdm,
					      double* y, double* dydx);
void dump_position_of_1st_and_Nth_atoms(int ny, int xn, double *y, double *data);
void dump_momentum_of_1st_and_Nth_atoms(int ny, int xn, double *y, double *data);
void dump_random_force_of_1st_and_Nth_atoms(int ny, int xn, double *y, double *data);

int truncation_time_step_number_estimate(double dt, int nt, double *alpha, double time_scale);
inline double propagator_auxiliary_integration(double * integration_kernal_l,int truncation_time_step_number, double *momentum, int momentum_for_current_time);
inline double propagator_auxiliary_integration_simpson(double * integration_kernal,int truncation_time_step_number, double *momentum, int momentum_for_current_time);
int initial_chain(double *y,int num_of_atoms);
int initial_chain_expcos(double *y,int num_of_atoms);
int initial_chain_UO_an(double *y,int num_of_atoms);

void data_analyse_chain(int nvar,double x,
			double *y, double* rlt, int nobv);
void data_analyse_chain_kenergy(int nvar,double x,
				double *y, double* rlt, int nobv);
void data_analyse_chain_penergy_mors(int nvar,double x,
				     double *y, double* rlt, int nobv);
void data_analyse_chain_penergy_harm(int nvar,double x,
				     double *y, double* rlt, int nobv);
void data_analyse_chain_total_energy_mors(int nvar,double x,
					  double *y, double* rlt, int nobv);
void data_analyse_chain_total_energy_harm(int nvar,double x,
					  double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_mors(int nvar,double x,
				       double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_mors_expcos(int nvar,double x,
					   double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_harm(int nvar,double x,
				       double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_mors_an(int nvar,double x,
				     double *y, double* rlt, int nobv);
void data_analyse_chain_UO_noise_an(int nvar,double x,
				     double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_harm_an(int nvar,double x,
				       double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_harm_expcos(int nvar,double x,
					      double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_harm_mult(int nvar,double x,
					    double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_fast(int nvar,double x,
				       double *y, double* rlt, int nobv);
void data_analyse_chain_heat_flux_fast_mult(int nvar,double x,
					    double *y, double* rlt, int nobv);
void data_analyse_chain_temperature_profile(int nvar,double x,
					    double *y, double* rlt, int nobv);
void data_analyse_chain_position_mult(int nvar,double x,
				      double *y, double* rlt, int nobv);
int memory_truncate(double *alpha_l, int nt);

#endif
