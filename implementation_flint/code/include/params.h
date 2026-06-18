#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
* Domain separators for XOF expansion
**************************************************/
#define DOMAIN_SEPARATOR_F 0
#define DOMAIN_SEPARATOR_E 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3

/*************************************************
* NTRU-TSampler parameters
**************************************************/
// Ring degree for the sampler
#define PARAM_N 1024
// Log of ring degree for the sampler
#define PARAM_LOG_N 10
// Ring degree for the proof
#define PARAM_N_PI 64
// Ring degree gap between R and R_pi
#define PARAM_K_PI 16
// Modulus for the sampler
#define PARAM_Q 123877L
// Low modulus
#define PARAM_QL 733L
// High modulus
#define PARAM_QH 169L
// Modulus inverse
#define PARAM_QL_MUL_QL_INVMOD_QH 63038L
#define PARAM_QH_MUL_QH_INVMOD_QL 60840L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 17
// G_H dimension
#define PARAM_KH 2
// G_H base
#define PARAM_BH 13
// Hamming weight of the tags
#define PARAM_W 8
// Bound on square spectral norm of e^T
#define PARAM_BESQ 5968.30937574009840318467
// Bound on spectral norm of the f_i
#define PARAM_BF 64.00000000000000000000
// Difference s_1² - s_L²
#define PARAM_S1SQ_SLSQ 23283790.35054315999150276184
// Squared Gaussian parameter s_2^2
#define PARAM_S2SQ 15979459.46655039116740226746
// Gaussian parameter s_2
#define PARAM_S2 3997.43160873959050150006
// Negated squared Gaussian parameter s_H^2
#define PARAM_NEG_SHSQ -1938.94565414150201831944
// Gaussian parameter s_L/q_L for Z-sampling of zL
#define PARAM_SL_DIV_QL 3.35446204446931162124
// Gaussian parameter s_H/b_H for Z-sampling of zH
#define PARAM_SH_DIV_BH 3.38718927843394679300
// Squared verification bound on v_1
#define PARAM_B1SQ 5980214451L
// Squared verification bound on [v_2 | v_3]
#define PARAM_B2SQ 8909931057L

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

#endif /* PARAMS_H */
