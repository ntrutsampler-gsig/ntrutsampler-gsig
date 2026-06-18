#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <fcntl.h>
#include <unistd.h>
#include <x86intrin.h>

#include "random.h"
#include "randombytes.h"
#include "params.h"
#include "fpr.h"

static void random_bytes(uint8_t * restrict);

/*************************************************
* Name:        uniform_int_distribution
*
* Description: Sample an integer uniformly in [0,n] using the uniform
*              distribution over [0, 2^32 - 1] provided by random_bytes
*
* Arguments:   - uint32_t n: upper bound on uniform element
* 
* Returns an integer uniformly distributed in [0,n]
**************************************************/
static uint32_t uniform_int_distribution(uint32_t n) {
  uint32_t scaling = (UINT32_MAX) / (n + 1);
  uint32_t past = (n + 1) * scaling;
  uint32_t r_data[4];

  while(true) {
    random_bytes((uint8_t *) r_data);
    for(int i = 0 ; i < 4 ; ++i) {
      uint32_t r = r_data[i];
      if(r < past) {
        return r / scaling;
      }
    }
  }
}

/*
  Code from random_aesni.c
*/

//macros
#define DO_ENC_BLOCK(m,k) \
  do{\
        m = _mm_xor_si128       (m, k[ 0]); \
        m = _mm_aesenc_si128    (m, k[ 1]); \
        m = _mm_aesenc_si128    (m, k[ 2]); \
        m = _mm_aesenc_si128    (m, k[ 3]); \
        m = _mm_aesenc_si128    (m, k[ 4]); \
        m = _mm_aesenc_si128    (m, k[ 5]); \
        __auto_type m5 = m;\
        m = _mm_aesenc_si128    (m, k[ 6]); \
        m = _mm_aesenc_si128    (m, k[ 7]); \
        m = _mm_aesenc_si128    (m, k[ 8]); \
        m = _mm_aesenc_si128    (m, k[ 9]); \
        m = _mm_aesenclast_si128(m, k[10]);\
        m = _mm_xor_si128(m, m5);\
    }while(0)

#define DO_DEC_BLOCK(m,k) \
  do{\
        m = _mm_xor_si128       (m, k[10+0]); \
        m = _mm_aesdec_si128    (m, k[10+1]); \
        m = _mm_aesdec_si128    (m, k[10+2]); \
        m = _mm_aesdec_si128    (m, k[10+3]); \
        m = _mm_aesdec_si128    (m, k[10+4]); \
        m = _mm_aesdec_si128    (m, k[10+5]); \
        m = _mm_aesdec_si128    (m, k[10+6]); \
        m = _mm_aesdec_si128    (m, k[10+7]); \
        m = _mm_aesdec_si128    (m, k[10+8]); \
        m = _mm_aesdec_si128    (m, k[10+9]); \
        m = _mm_aesdeclast_si128(m, k[0]);\
    }while(0)

#define AES_128_key_exp(k, rcon) aes_128_key_expansion(k, _mm_aeskeygenassist_si128(k, rcon))

//the expanded key
static __m128i key_schedule[20];
static uint64_t ctr = 0;

static __m128i aes_128_key_expansion(__m128i key, __m128i keygened){
  keygened = _mm_shuffle_epi32(keygened, _MM_SHUFFLE(3,3,3,3));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  key = _mm_xor_si128(key, _mm_slli_si128(key, 4));
  return _mm_xor_si128(key, keygened);
}

static void aes128_load_key(const int8_t * enc_key){
    key_schedule[0] = _mm_loadu_si128((const __m128i*) enc_key);
  key_schedule[1]  = AES_128_key_exp(key_schedule[0], 0x01);
  key_schedule[2]  = AES_128_key_exp(key_schedule[1], 0x02);
  key_schedule[3]  = AES_128_key_exp(key_schedule[2], 0x04);
  key_schedule[4]  = AES_128_key_exp(key_schedule[3], 0x08);
  key_schedule[5]  = AES_128_key_exp(key_schedule[4], 0x10);
  key_schedule[6]  = AES_128_key_exp(key_schedule[5], 0x20);
  key_schedule[7]  = AES_128_key_exp(key_schedule[6], 0x40);
  key_schedule[8]  = AES_128_key_exp(key_schedule[7], 0x80);
  key_schedule[9]  = AES_128_key_exp(key_schedule[8], 0x1B);
  key_schedule[10] = AES_128_key_exp(key_schedule[9], 0x36);

  // generate decryption keys in reverse order.
    // k[10] is shared by last encryption and first decryption rounds
    // k[0] is shared by first encryption round and last decryption round (and is the original user key)
    // For some implementation reasons, decryption key schedule is NOT the encryption key schedule in reverse order
  key_schedule[11] = _mm_aesimc_si128(key_schedule[9]);
  key_schedule[12] = _mm_aesimc_si128(key_schedule[8]);
  key_schedule[13] = _mm_aesimc_si128(key_schedule[7]);
  key_schedule[14] = _mm_aesimc_si128(key_schedule[6]);
  key_schedule[15] = _mm_aesimc_si128(key_schedule[5]);
  key_schedule[16] = _mm_aesimc_si128(key_schedule[4]);
  key_schedule[17] = _mm_aesimc_si128(key_schedule[3]);
  key_schedule[18] = _mm_aesimc_si128(key_schedule[2]);
  key_schedule[19] = _mm_aesimc_si128(key_schedule[1]);
}

static void aes128_enc(const int8_t * restrict plainText, int8_t * restrict cipherText){
    __m128i m = _mm_loadu_si128((const __m128i *) plainText);

    DO_ENC_BLOCK(m,key_schedule);

    _mm_storeu_si128((__m128i *) cipherText, m);
}

//public API
void random_init(void) {
  // Seed the AES key
    uint8_t seed[16];
    randombytes(seed, sizeof(uint8_t)*16);
    aes128_load_key((int8_t *) seed); // don't know why the AES key is signed but whatever
    ctr = 0;
}

static void random_bytes(uint8_t * restrict data) {
    uint8_t plaintext[16] = {0};
    *((uint64_t *) plaintext) = ctr++;
    aes128_enc((int8_t *) plaintext, (int8_t *) data); // let's cast them to signed because why not
}

/* #############################################
 * # Integer Sampler for Large Variable Widths #
 * #############################################
 *
 * For the perturbation sampling, we  need to sample integers from 
 * potentially large and variable widths, with also variable centers. 
 * For now, we use the following sampler.
 */

/*
  Code from exp_aes.cpp adapted to custom FPR
*/
static fpr algorithm_EA(uint64_t * n) {
  const fpr fpr_ln2 = 0.6931471805599453; // ln(2)              
  const fpr fpr_a = 5.713363152645423;   // (4 + 3*sqrt(2)) * ln(2)
  const fpr fpr_b = 3.414213562373095;   // 2 + sqrt(2)
  const fpr fpr_c = -1.6734053240284923; // -(1 + sqrt(2)) * ln(2)
  const fpr fpr_p = 0.9802581434685472;   // sqrt(2) * ln(2)
  const fpr fpr_A = 5.6005707569738075;   // (4 + 3*sqrt(2)) * sqrt(2) * ln(2)^2 = a * p       
  const fpr fpr_B = 3.3468106480569846;   // (2 + sqrt(2)) * sqrt(2) * ln(2)     = b * p
  //const fpr fpr_z = 0.32659467597150770857;   // -(1 + sqrt(2)) * ln(2) + 2          = c + 2
  //const fpr fpr_r = 1.01008958200144838280;   //                                     = 8*exp(-z)*b*(b-1)*log(2)/a^2  
  //const fpr fpr_h = 0.02983143853290115555;   //                                     = r - p
  const fpr fpr_D = 0.08578643762690495;   // 1/(2 + sqrt(2))^2                   = 1 / (b*b)
  const fpr fpr_H = 0.0026106723602095233;   //                                     = h * D / p           

  uint8_t data[16];
  random_bytes(data);

  const unsigned long long int II = *((uint64_t *) data);
  const int j = __builtin_ctzll(II);
  const fpr fpr_j = fpr_of((int64_t)j);
  const fpr fpr_G = fpr_add(fpr_c, fpr_mul(fpr_j, fpr_ln2));                                            // c + j*ln2

  const fpr fpr_U = fpr_div(fpr_of((int64_t)(*((uint16_t *) (data + 8)))), PARAM_UINT16_MAX);             // *((uint16_t *) (data + 8)) / (UINT16_MAX*1.0)  [double]
  random_bytes(data);

  (*n)++;
  if (fpr_leq(fpr_U, fpr_p)) {                                                                          // if U <= p
    return fpr_add(fpr_G, fpr_div(fpr_A, fpr_sub(fpr_B, fpr_U)));                                       //    return G + A/(B - U)
  }

  while(true) {
    for(unsigned i = 0; i < 4; i++) {
      (*n)++;
      const fpr fpr_U1 = fpr_div(fpr_of((int64_t)(*((uint16_t *) (data + i*2)))), PARAM_UINT16_MAX);      // *((uint16_t *) (data + i*2)) / (UINT16_MAX*1.0)  [double]
      const fpr fpr_U2 = fpr_div(fpr_of((int64_t)(*((uint16_t *) (data + i*2 + 8)))), PARAM_UINT16_MAX);  // *((uint16_t *) (data + i*2 + 8)) / (UINT16_MAX*1.0)  [double]

      const fpr fpr_bU1 = fpr_sub(fpr_b, fpr_U1);                                                       // b - U1
      const fpr fpr_Y = fpr_div(fpr_a, fpr_bU1);                                                        // a / bU1
      const fpr fpr_L = fpr_mul(fpr_add(fpr_D, fpr_mul(fpr_H, fpr_U2)), fpr_sqr(fpr_bU1));              // (D + H*U2) * bU1 * bU1
      const fpr fpr_Z = fpr_add(fpr_Y, fpr_c);                                                          // Y + c

      const fpr fpr_LZ = fpr_sub(fpr_add(fpr_L, fpr_Z), PARAM_ONE);                                       // L + Z - 1
      const fpr fpr_ZZ = fpr_sqr(fpr_Z);                                                                // Z * Z

      const fpr fpr_eval_cond2 = fpr_sub(fpr_add(fpr_LZ, fpr_div(fpr_mul(fpr_ZZ, fpr_Z), PARAM_SIX)), fpr_half(fpr_ZZ)); // LZ + ZZ*Z/6 - ZZ/2

      const bool cond1 = (bool)(fpr_leq(fpr_LZ, PARAM_ZERO));                                             // cond1 = LZ <= 0
      const bool cond2 = (bool)(fpr_leq(fpr_eval_cond2, PARAM_ZERO));                                     // cond2 = LZ + ZZ*Z/6 - ZZ/2 <= 0
      const bool cond3 = (bool)(fpr_leq(fpr_sub(fpr_LZ, fpr_ZZ), PARAM_ZERO));                            // cond3 = LZ - ZZ <= 0
      const bool cond4 = (bool)(fpr_leq(fpr_L, fpr_expm(fpr_Z)));                                       // cond4 = L <= exp(-Z)  (we can show Z < ln(2))

      if (cond1 || cond2 || (cond3 && cond4)) {
        return fpr_add(fpr_G, fpr_Y);
      }
    }
    random_bytes(data);
  }
}

/*
  Code from algoF_aes.cpp
*/

static int algorithmF(fpr mu, fpr sigma) {
  const int64_t sigma_floor = fpr_floor(sigma);
  const fpr fpr_sigma_floor = fpr_of(sigma_floor);
  const unsigned sigma_max = (unsigned) (((fpr_sigma_floor == sigma) && (sigma_floor != 0)) ? (sigma_floor - 1) : (sigma_floor));

  uint8_t randomData[16];
  random_bytes(randomData);
  while(true) {
    for(unsigned i = 0; i < 16; i++) {
      uint64_t n = 0;
      fpr fpr_k = fpr_of(fpr_ceil(fpr_double(algorithm_EA(&n))) - 1);

      int s = (randomData[i] > 127) ? -1 : 1;
      fpr fpr_s = fpr_of((int64_t)s);

      unsigned j = uniform_int_distribution(sigma_max);
      fpr fpr_j = fpr_of((int64_t)j);

      fpr tmp = fpr_add(fpr_mul(sigma, fpr_k), fpr_mul(fpr_s, mu)); 
      int64_t i0 = (int64_t) fpr_ceil(tmp);
      fpr fpr_i0 = fpr_of(i0);

      fpr fpr_x0 = fpr_div(fpr_sub(fpr_i0, tmp), sigma);
      fpr fpr_x = fpr_add(fpr_x0, fpr_div(fpr_j, sigma));

      fpr fpr_k1 = fpr_of(fpr_ceil(fpr_double(algorithm_EA(&n))) - 1);
      fpr fpr_z = algorithm_EA(&n);

      bool cond1 = (bool)(fpr_leq(fpr_mul(fpr_k, fpr_sub(fpr_k, PARAM_ONE)), fpr_k1));  // k1 >= k*(k-1)
      bool cond2 = (bool)(fpr_lt(fpr_x, PARAM_ONE));                                    // x < 1
      bool cond3a = fpr_x != PARAM_ZERO;                                                // x != 0
      bool cond3b = fpr_k != PARAM_ZERO;                                                // k != 0
      bool cond3c = s == 1;                                                           // s == 1
      bool cond3 = cond3a || cond3b || cond3c;                                        
      bool cond4 = (bool)(fpr_lt(fpr_mul(fpr_x, fpr_add(fpr_k, fpr_half(fpr_x))),fpr_z)); // z > 0.5*x*(2*k + x)

      if (cond1 && cond2 && cond3 && cond4) {
        return s*(((int)i0) + ((int) j));
      }
    }
    random_bytes(randomData);
  }
}

/*************************************************
* Name:        sampleZ
*
* Description: Sample an integer from the discrete Gaussian
*              distribution with density proportional to 
*              exp(-pi * |x-c|^2 / sigma^2)
*
* Arguments:   - fpr c: center for the Gaussian distribution
*              - fpr sigma: Gaussian parameter (std dev = sigma / sqrt(2*pi))
* 
* Returns a Gaussian sample from D_{Z, sigma, c}
**************************************************/
int64_t SampleZ(fpr c, fpr sigma) {
  return algorithmF(c, fpr_mul(sigma, PARAM_INV_SQRT_2PI));
}

#if 0
/* ####################################
 * # Integer Sampler for Small Widths #
 * ####################################
 *
 * For the gadget sampling, we only need to sample integers from small
 * fixed widths, i.e., roughly at the smoothing parameter. We however 
 * need to sample with variable centers. We thus use the sampler from
 * Falcon, which satisfies all the same constraints. 
 */

/*************************************************
* Name:        base_sampler
*
* Description: Sample an integer along a half-Gaussian distribution
*              centered on zero and standard deviation 1.8205, with
*              72 bits of precision.
**************************************************/
static int base_sampler(void) {
  static const uint32_t dist[] = {
    10745844u,  3068844u,  3741698u,
     5559083u,  1580863u,  8248194u,
     2260429u, 13669192u,  2736639u,
      708981u,  4421575u, 10046180u,
      169348u,  7122675u,  4136815u,
       30538u, 13063405u,  7650655u,
        4132u, 14505003u,  7826148u,
         417u, 16768101u, 11363290u,
          31u,  8444042u,  8086568u,
           1u, 12844466u,   265321u,
           0u,  1232676u, 13644283u,
           0u,    38047u,  9111839u,
           0u,      870u,  6138264u,
           0u,       14u, 12545723u,
           0u,        0u,  3104126u,
           0u,        0u,    28824u,
           0u,        0u,      198u,
           0u,        0u,        1u
  };

  size_t u;
  uint8_t data[16];
  uint32_t v0, v1, v2;
  int z;

  /*
   * Get a random 72-bit value, into three 24-bit limbs v0..v2.
   */
  random_bytes(data);
  v0 = (uint32_t)data[0] | ((uint32_t)data[1] << 8) | ((uint32_t)data[2] << 16);
  v1 = (uint32_t)data[3] | ((uint32_t)data[4] << 8) | ((uint32_t)data[5] << 16);
  v2 = (uint32_t)data[6] | ((uint32_t)data[7] << 8) | ((uint32_t)data[8] << 16);

  /*
   * Sampled value is z, such that v0..v2 is lower than the first
   * z elements of the table
   */
  z = 0;
  for (u = 0; u < (sizeof dist) / sizeof(dist[0]); u += 3) {
    uint32_t w0, w1, w2, cc;
    w0 = dist[u + 2];
    w1 = dist[u + 1];
    w2 = dist[u + 0];
    cc = (v0 - w0) >> 31;
    cc = (v1 - w1 - cc) >> 31;
    cc = (v2 - w2 - cc) >> 31;
    z += (int)cc;
  }

  return z;
}

/*************************************************
* Name:        ber_exp
*
* Description: Sample a bit with probability ccs.exp(-x) for some x >= 0
* 
* Arguments:   - fpr x: exponent
*              - fpr ccs: scaling factor
**************************************************/
static int ber_exp(fpr x, fpr ccs) {
  int s, i, off;
  fpr r;
  uint32_t sw, w;
  uint64_t z;
  uint8_t data[16];
  /*
   * Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
   * and 0 <= r < log(2). Since x >= 0, we can use fpr_trunc().
   */
  s = (int)fpr_trunc(fpr_mul(x, PARAM_INV_LOG2));
  r = fpr_sub(x, fpr_mul(fpr_of(s), PARAM_LOG2));

  /*
   * It may happen (quite rarely) that s >= 64; if sigma = 1.2
   * (the minimum value for sigma), r = 0 and b = 1, then we get
   * s >= 64 if the half-Gaussian produced a z >= 13, which happens
   * with probability about 0.000000000230383991, which is
   * approximatively equal to 2^(-32). In any case, if s >= 64,
   * then BerExp will be non-zero with probability less than
   * 2^(-64), so we can simply saturate s at 63.
   */
  sw = (uint32_t)s;
  sw ^= (sw ^ 63) & -((63 - sw) >> 31);
  s = (int)sw;

  /*
   * Compute exp(-r); we know that 0 <= r < log(2) at this point, so
   * we can use fpr_expm_p63(), which yields a result scaled to 2^63.
   * We scale it up to 2^64, then right-shift it by s bits because
   * we really want exp(-x) = 2^(-s)*exp(-r).
   *
   * The "-1" operation makes sure that the value fits on 64 bits
   * (i.e. if r = 0, we may get 2^64, and we prefer 2^64-1 in that
   * case). The bias is negligible since fpr_expm_p63() only computes
   * with 51 bits of precision or so.
   */
  z = ((fpr_expm_p63(r, ccs) << 1) - 1) >> s;

  /*
   * Sample a bit with probability exp(-x). Since x = s*log(2) + r,
   * exp(-x) = 2^-s * exp(-r), we compare lazily exp(-x) with the
   * PRNG output to limit its consumption, the sign of the difference
   * yields the expected result.
   */
  random_bytes(data);
  off = 0;
  i = 64;
  do {
    i -= 8;
    w = data[off] - ((uint32_t)(z >> i) & 0xFF);
    off++;
  } while (!w && i > 0);
  return (int)(w >> 31);
}

/*************************************************
* Name:        sampleZ_small_wdth
*
* Description: Sample a discrete Gaussian integer for a 
*              variable center, and small width in [sigma_min, sigma_max]
* 
* Arguments:   - fpr x: exponent
*              - fpr ccs: scaling factor
**************************************************/
int sampleZ_small_wdth(fpr mu, fpr isigma) {
  int s;
  fpr r, dss, ccs;
  uint8_t byte;

  /*
   * Center is mu. We compute mu = s + r where s is an integer
   * and 0 <= r < 1.
   */
  s = (int)fpr_floor(mu);
  r = fpr_sub(mu, fpr_of(s));

  /*
   * dss = 1/(2*sigma^2) = 0.5*(isigma^2).
   */
  dss = fpr_half(fpr_sqr(isigma));

  /*
   * ccs = sigma_min / sigma = sigma_min * isigma.
   */
  ccs = fpr_mul(isigma, PARAM_SIGMAMIN);

  /*
   * We now need to sample on center r.
   */
  for (;;) {
    int z0, z, b;
    fpr x;

    /*
     * Sample z for a Gaussian distribution. Then get a
     * random bit b to turn the sampling into a bimodal
     * distribution: if b = 1, we use z+1, otherwise we
     * use -z. We thus have two situations:
     *
     *  - b = 1: z >= 1 and sampled against a Gaussian
     *    centered on 1.
     *  - b = 0: z <= 0 and sampled against a Gaussian
     *    centered on 0.
     */
    z0 = base_sampler();
    randombytes(byte, 1);
    b = (int)byte & 1;
    z = b + ((b << 1) - 1) * z0;

    /*
     * Rejection sampling. We want a Gaussian centered on r;
     * but we sampled against a Gaussian centered on b (0 or
     * 1). But we know that z is always in the range where
     * our sampling distribution is greater than the Gaussian
     * distribution, so rejection works.
     *
     * We got z with distribution:
     *    G(z) = exp(-((z-b)^2)/(2*sigma0^2))
     * We target distribution:
     *    S(z) = exp(-((z-r)^2)/(2*sigma^2))
     * Rejection sampling works by keeping the value z with
     * probability S(z)/G(z), and starting again otherwise.
     * This requires S(z) <= G(z), which is the case here.
     * Thus, we simply need to keep our z with probability:
     *    P = exp(-x)
     * where:
     *    x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
     *
     * Here, we scale up the Bernouilli distribution, which
     * makes rejection more probable, but makes rejection
     * rate sufficiently decorrelated from the Gaussian
     * center and standard deviation that the whole sampler
     * can be said to be constant-time.
     */
    x = fpr_mul(fpr_sqr(fpr_sub(fpr_of(z), r)), dss);
    x = fpr_sub(x, fpr_mul(fpr_of(z0 * z0), PARAM_INV_2_SIGMAMAX_SQ)); // fpr_inv_2sqrsigma0));
    if (ber_exp(x, ccs)) {
      /*
       * Rejection sampling was centered on r, but the
       * actual center is mu = s + r.
       */
      return s + z;
    }
  }
}
#endif
