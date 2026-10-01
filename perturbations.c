// DYNAMIC PHANTOM METRIC: SUPERLUMINAL INJECTION & TRANSIENT BRAKING
double current_a = pvecback[pba->index_bg_a];
double damping_factor = pba->a_eq / current_a;
double K_phantom = 0.66688 * (a_prime_over_a * a_prime_over_a);
double Gamma_against_dynamic = 0.0397 * damping_factor;
double Gamma_with = 0.0;

double active_gamma = (current_theta_b > 0.0) ? Gamma_against_dynamic : Gamma_with;
dy[pv->index_pt_theta_b] += - (active_gamma * current_theta_b) - (K_phantom * current_delta_b);