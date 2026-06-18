# -*- coding: utf-8 -*-
"""
    Brief: Parameter Selection and Security Estimation Script
"""
#-- Import --#
from math import sqrt, exp, log, log2, floor, ceil, pi, e, gcd
from scipy.special import comb
from scipy.optimize import root
from sympy import isprime, prevprime, nextprime, divisors
from estimate_SIS_LWE import *
#-- End Import --#

#-- Global Parameters --#
COST_MODEL      = 'realistic_sieving'
LOG2_EPS        = -40
Q_SIGN          = 2**32 # Maximal number of group signature queries per user (hardcoded)
#-- End Global Parameters --#

def c_star(dim:int, sec:float):
    """
    Find the Gaussian tailcut rate c so that the upper bound c*s*sqrt(dim) is
    verified with probability at least 1 - 2^(-sec)
    - input:
        (int)   dim     -- Dimension
        (int)   sec     -- Security parameter
    - output:
        (flt)   c_star  -- Tailcut rate
    """
    f = lambda c: sec + dim * (log2(c * sqrt(2*pi)) + (1/2 - pi*c**2)*log2(e))
    return root(f, 1)['x'][0]

def t_star(dim:int, sec:float):
    """
    Find the Gaussian tailcut rate t so that the upper bound t*s in
    Infinity norm is verified with probability at least 1 - 2^(-sec)
    - input:
        (int)   dim     -- Dimension
        (int)   sec     -- Security parameter
    - output:
        (flt)   t_star  -- Tailcut rate
    """
    return sqrt((sec + 1 + log2(dim)) / (pi * log2(e)))

def centered_mod(a:int, q:int):
    r = a % q
    return (r if r <= (q-1)/2 else r - q)

class NTRUSampler_Parameters:
    """
    Main class containing all the parameters for the sampler
    """

    def __init__(self, target_bitsec:int, n:int, n_pi:int, b_H:int, k_H:int, q_L:int, N_min:int):
        """
        Computing all the group signature parameters as in GS.Setup
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   n_pi            -- Ring degree for the proof system
            (int)   b_H             -- Gadget base for G_H
            (int)   k_H             -- Gadget dimension for G_H
            (int)   q_L             -- Low modulus factor
            (int)   N_min           -- Minimal number of group members 
        """
        self.sec = target_bitsec

        # Degrees of the ring R = Z[X]/<X^n + 1> and R_π = Z[X]/<X^n_π + 1>
        self.n = n
        self.n_pi = n_pi

        # Number of splitting factors for b_H is hardcoded: kappa = 2
        assert isprime(b_H) and (b_H % 8 == 5) and b_H >= 2, "b_H does not satisfy the correct properties"
        self.b_H = b_H
        self.k_H = k_H 
        self.q_H = self.b_H ** self.k_H

        assert isprime(q_L) and (q_L % 8 == 5) and q_L >= 2, "q_L does not satisfy the correct properties"
        self.q_L = q_L

        self.q = self.q_L * self.q_H

        # Finding minimum Hamming weight for the tag space
        w = 1
        while comb(self.n_pi, w) < N_min:
            w += 1
        self.w = w
        self.N = comb(self.n_pi, self.w)

        # Smoothing parameters
        self.smoothing_sL = sqrt((log(2 * self.n) - LOG2_EPS * log(2)) / pi)
        self.smoothing_sH = sqrt((log(2 * self.n * self.k_H) - LOG2_EPS * log(2)) / pi)
        self.eta = sqrt((log(2 * self.n * (1 + self.k_H)) - LOG2_EPS * log(2)) / pi)

        # Sampler key bounds (heuristically determined)
        self.B_e = sqrt(self.n) + sqrt(self.n*self.k_H) 
        self.B_f = 2*sqrt(self.n)

        # Computing Gaussian widths for NTRU-Sampler
        self.s_L = self.smoothing_sL * self.q_L
        self.s_H = self.smoothing_sH * self.b_H
        self.gamma = sqrt(2)
        self.s_1 = sqrt(self.s_L**4/(self.s_L**2 - self.eta**2) + self.gamma**2/(self.gamma**2-1) * self.s_H**4/(self.s_H**2 - self.eta**2) * self.B_e**2)
        self.s_2 = self.gamma * self.s_H**2/sqrt(self.s_H**2 - self.eta**2) * self.B_f

        # User secret key verification bounds
        tail_cut_prob_log = 20
        tail_cut_prob = 2**(-tail_cut_prob_log)

        self.B_1_s = floor(c_star(self.n, tail_cut_prob_log)**2 * self.s_1**2 * (self.n))
        self.B_2_s = floor(c_star(self.n*(1 + self.k_H), tail_cut_prob_log)**2 * self.s_2**2 * (self.n*(1 + self.k_H)))
        self.B_s = self.B_1_s + self.B_2_s
        
        self.B_1 = sqrt(self.B_1_s)
        self.B_2 = sqrt(self.B_2_s)
        self.B = sqrt(self.B_s)

        ### Efficiency
        
        self.ipk_bitsize    = self.n * self.k_H * ceil(log2(self.q)) + 256   # issuer pk = {h, seed}
        self.isk_bitsize    = 2 * self.n * self.k_H * ceil(log2(3))          # issuer sk = {e, f}

        # gsk[i] = {id_i, v_2i, v_3i}
        self.id_bitsize    = self.n_pi
        self.v1_bitsize = ceil(self.n * (1/2 + log2(self.s_1)))
        self.v2_v3_bitsize = ceil((1 + self.k_H) * self.n * (1/2 + log2(self.s_2)))

        self.gsk_i_bitsize  = self.id_bitsize + self.v2_v3_bitsize # v_1 recovered by user

    def __repr__(self):
        """
        Printing a NTRUSampler_Parameters object
        """
        tmp = '\n[*] GROUP SIGNATURE PARAMETERS \n'
        tmp += 105 * '=' + '\n'
        tmp += '| {:60s} | {:^15s} | 2^{:<18.5f} |\n'.format('Group Size', 'N', log2(self.N))
        tmp += '[+] TSampler Parameters\n'
        tmp += 105 * '-' + '\n'
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Sampler ring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Low modulus', 'q_L', self.q_L)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('High modulus', 'q_H', self.q_H)        
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('G_H base', 'b_H', self.b_H)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('G_H dimension', 'k_H', self.k_H)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('z_L sampling width', 's_L', self.s_L)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('z_H sampling width', 's_H', self.s_H)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Preimage sampling width 1', 's_1', self.s_1)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Preimage sampling width 2', 's_2', self.s_2)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Hamming weights of tags', 'w', self.w)
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Smoothing loss', 'ε', LOG2_EPS)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Smoothing of Z^n', 'η(n)', self.smoothing_sL)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Smoothing of Z^(nk_H)', 'η(nk_H)', self.smoothing_sH)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Smoothing of Z^(n(1+k_H))', 'η', self.eta)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Key verification bound 1', 'B_1', self.B_1)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Key verification bound 2', 'B_2', self.B_2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Key verification bound', 'B', self.B)
        tmp += 105 * '=' + '\n'
        tmp += '[+] Group Signature Detailed Estimated Performance (KB)\n'
        tmp += 105 * '=' + '\n'
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Issuer public key size (B)', '|ipk|', self.ipk_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Issuer secret key size (B)', '|gmsk|', self.isk_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Identity size (B)', '|id|', self.id_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('v_1 size (B)', '|v_1|', self.v1_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('v_2, v_3 size (B)', '|v_2,v_3|', self.v2_v3_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Overall User Secret Key size (B)', '|gsk[i]|', self.gsk_i_bitsize / 2 ** 3.)
        return tmp

class CPA_PKE_Parameters:
    """
    Main class containing all the parameters for the 
    verifiable encryption scheme (CPA-anonymous case). 
    """

    def __init__(self, target_bitsec:int, n_op:int, d:int):
        """
        Computing all the signature parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n_op            -- Ring degree
            (int)   d               -- Module rank
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R_op = Z[X]/<X^n_op + 1>
        self.n_op = n_op

        # Encryption module rank 
        self.d = d

        # Computing opening parameters
        self.xi = 1.15
        self.B_op = sqrt(floor(self.xi**2 * self.n_op * self.d))

        # Computing proof bounds on (s,e) and e1
        self.B_se_s = floor(self.xi**2 * self.n_op * self.d)       
        self.B_e1_s = floor(self.xi**2 * self.n_op / 2)       
        self.B_se = sqrt(self.B_se_s)
        self.B_e1 = sqrt(self.B_e1_s)

        # Search proper modulus
        self.p = nextprime(ceil(4 * (self.B_e1 + self.B_op * self.B_se))) # B_e1 + B_op*B_se < p/4
        # while p % self.n != self.n/2 + 1: # splits into n/4 factors %(2*self.n) not in [1, self.n+1]: # splits into n/2 or n factors
        #     p = nextprime(p)
        # self.p = p

        # ct_0, ct_1, ct_2
        self.ct0_bitsize    = self.d * self.n_op * ceil(log2(self.p))
        self.ct1_bitsize    = self.n_op * ceil(log2(self.p))
        self.ct2_bitsize    = self.n_op * ceil(log2(self.p))
        self.ct_bitsize     = self.ct0_bitsize + self.ct1_bitsize + self.ct2_bitsize

        self.opk_bitsize    = self.d * self.n_op * ceil(log2(self.p)) # opener pk = {b_{op,1}}
        self.osk_bitsize    = self.d * self.n_op * ceil(log2(3))          # opener sk = s_{op,1}

    def __repr__(self):
        """
        Printing a PKE_Parameters object
        """
        tmp = '\n[+] VERIFIABLE PUBLIC KEY ENCRYPTION PARAMETERS\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Ring degree of Z[X]/(X^n_op + 1)', 'n_op', self.n_op)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'p', self.p)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on opening key', 'B_op', self.B_op)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Proof bound on (s,e)', 'B_se', self.B_se)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Proof bound on e_1, e_2', 'B_e1', self.B_e1)
        tmp += '\n[+] PKE Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Public key size (KB)', '|opk|', self.opk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Secret key size (KB)', '|osk|', self.osk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Ciphertext size (B)', '|ct|', (self.ct_bitsize - self.ct2_bitsize) / 2 ** 3.)
        return tmp    

class CCA_PKE_Parameters:
    """
    Main class containing all the parameters for the 
    verifiable encryption scheme (CCA-anonymous case). 
    """

    def __init__(self, target_bitsec:int, n_op:int, d:int):
        """
        Computing all the encryption parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n_op            -- Ring degree
            (int)   d               -- Module rank
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R_op = Z[X]/<X^n_op + 1>
        self.n_op = n_op

        # Encryption module rank 
        self.d = d

        # Computing opening parameters
        self.xi = 1.15
        self.beta_op = 3/4*(sqrt(self.n_op) + sqrt(2*self.n_op*self.d) + 6) 
        self.B_op = sqrt(floor(self.xi**2 * self.n_op * self.d))

        self.log2_eps_op = -(self.sec + 3)
        self.eps_op = 2**self.log2_eps_op
        self.eta_op = sqrt((log(4 * self.n_op * self.d) - self.log2_eps_op * log(2)) / pi)
        self.s_MLWE = sqrt(2) * self.eta_op 

        self.alpha = 1.88 # 1.65 
        b = sqrt(2) * self.alpha / sqrt(self.alpha**2 - 2)
        self.s_se = self.alpha * self.s_MLWE 
        self.s_e1 = b * self.s_MLWE * self.beta_op 
        assert abs(1/sqrt(2 * (1/self.s_se**2 + self.beta_op**2/self.s_e1**2)) - self.s_MLWE) < 1e-13 # necessary condition for reduction

        # Computing proof bounds on (s,e) and e1
        tail_cut_prob_log = 20
        tail_cut_prob = 2**(-tail_cut_prob_log)
        self.B_se_s = floor(c_star(2*self.d*self.n_op, tail_cut_prob_log)**2 * self.s_se**2 * (2*self.d*self.n_op))       
        self.B_e1_s = floor(c_star(self.n_op, tail_cut_prob_log)**2 * self.s_e1**2 * self.n_op)       
        self.B_se = sqrt(self.B_se_s)
        self.B_e1 = sqrt(self.B_e1_s)

        # Search proper modulus
        self.p = nextprime(ceil(4 * (self.B_e1 + self.B_op * self.B_se))) # B_e1 + B_op*B_se < p/4
        # while p % self.n != self.n/2 + 1: # splits into n/4 factors %(2*self.n) not in [1, self.n+1]: # splits into n/2 or n factors
        #     p = nextprime(p)
        # self.p = p

        # ct_0, ct_1, ct_2
        self.ct0_bitsize    = self.d * self.n_op * ceil(log2(self.p))
        self.ct1_bitsize    = self.n_op * ceil(log2(self.p))
        self.ct2_bitsize    = self.n_op * ceil(log2(self.p))
        self.ct_bitsize     = self.ct0_bitsize + self.ct1_bitsize + self.ct2_bitsize

        self.opk_bitsize    = 2 * self.d * self.n_op * ceil(log2(self.p))   # opener pk = {b_{op,1}, b_{op,2}}
        self.osk_bitsize    = self.d * self.n_op * ceil(log2(3))            # opener sk = s_{op,1}

    def __repr__(self):
        """
        Printing a PKE_Parameters object
        """
        tmp = '\n[+] VERIFIABLE PUBLIC KEY ENCRYPTION PARAMETERS\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Encryption ring degree of Z[X]/(X^n_op + 1)', 'n_op', self.n_op)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'p', self.p)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Spectral norm bound on opening key', 'β_op', self.beta_op)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on opening key', 'B_op', self.B_op)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Matrix-Hint reduction smoothing loss', 'ε_op', self.log2_eps_op)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(2nd)', 'η_op', self.eta_op)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian parameter for (s,e)', 's_se', self.s_se)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian parameter for e_1, e_2 (hint mask)', 's_e1', self.s_e1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Proof bound on (s,e)', 'B_se', self.B_se)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Proof bound on e_1, e_2', 'B_e1', self.B_e1)
        tmp += '\n[+] PKE Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Public key size (B)', '|opk|', self.opk_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Secret key size (B)', '|osk|', self.osk_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Ciphertext size (B)', '|ct|', self.ct_bitsize / 2 ** 3.)
        return tmp     

class ZKP_Parameters:
    """
    Main class containing all the parameters for the 
    zero-knowledge proof system for the group signature. 
    """

    def __init__(self, target_bitsec:int, n_pi:int, d_pi:int, m_2:int, q_start:int, gamma:int, D:int, samp, pke, bimodal=True, compression=True, garbage=True, CCA=False, infARP=False):
        """
        Computing all the zero-knowledge argument parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m_2             -- Commitment randomness dimension
            (int)   q_1_start       -- Starting search modulus (q_1 largest adequately splitted prime below q_1_start)
            (int)   gamma           -- Compression parameter for commitment w
            (int)   D               -- Compression parameter for commitment t_A
            (NTRUSampler_Parameters) samp -- TSampler parameters
            (PKE_Parameters) pke    -- Encryption parameters
        """
        ### Parameters
        self.infARP = infARP

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R'' = Z[X]/<X^n' + 1>
        self.n_pi = n_pi

        # M-SIS module rank 
        self.d_pi = d_pi

        # # Number of splitting factors for modulus (hardcoded)
        self.kappa = 2

        if compression:
            # Finding the largest prime modulus splitting in kappa factors below q_start and such that a divisor of q-1 is close to gamma
            found_q_gamma = False
            q_pi = q_start + 1
            while not(found_q_gamma):
                q_pi = prevprime(q_pi)
                while q_pi % 8 != 5:
                    q_pi = prevprime(q_pi)
                divs = divisors(q_pi-1)
                for div in divs:
                    if (gamma <= div <= 5*gamma/4) and (div % 2 == 0): # find even divisor closest to gamma (not larger than 5.gamma/4)
                        self.gamma = div
                        found_q_gamma = True 
                        break
            # Commitment compression
            self.D = D 
        else:
            self.gamma = 0
            self.D = 0
            q_pi = prevprime(q_start + 1)
            while (q_pi%(4*self.kappa) != 2*self.kappa + 1) or (q_pi < (2*sqrt(self.kappa)) ** self.kappa):
                q_pi = prevprime(q_pi)
        self.q_pi = q_pi
        self.q_min = q_pi

        if garbage:
            # Repetition for soundness amplification
            l = ceil(self.sec / log2(self.q_min))
            self.l = ceil(l/2)
        else:
            self.l = ceil(self.sec / log2(self.q_min))

        # Infinity norm bound on the challenges
        self.rho = ceil(1/2 * (2 ** (2*(self.sec + 1)/self.n_pi) - 1))

        # Manhattan-like norm bound on the challenges (hardcoded)
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n_pi]

        # Size of challenge space
        self.challenge_space_size = (2 * self.rho + 1) ** (self.n_pi // 2) / 2

        # Subring gap
        self.k_pi = samp.n // self.n_pi
        self.k_op_pi = pke.n_op // self.n_pi

        # Witness dimension
        if CCA:
            self.m_11 = (2 * pke.d) * self.k_op_pi + 1  # s,e (+ one element for four-square norm proof)
            self.m_11 += 1                              # id 
            if bimodal:
                self.m_11 += 1                          # bimodal bit
            
            self.m_12 = (2 + samp.k_H) * self.k_pi + 1  # v (+ one element for four-square norm proof)
            self.m_12 += self.k_op_pi + 1               # e1 (+ one element for four-square norm proof)
            self.m_12 += self.k_op_pi + 1               # e2 (+ one element for four-square norm proof)
        else:
            self.m_11 = (2 * pke.d) * self.k_op_pi + 1  # s,e (+ one element for four-square norm proof)
            self.m_11 += 1                              # tag 
            if bimodal:
                self.m_11 += 1                          # bimodal bit
            self.m_11 += self.k_op_pi + 1               # e1 (+ one element for four-square norm proof)

            self.m_12 = (2 + samp.k_H) * self.k_pi + 1  # v (+ one element for four-square norm proof)
        self.m_1 = self.m_11 + self.m_12

        # Commitment randomness dimension and infinity norm bound (hardcoded)
        self.m_2 = m_2
        self.xi_s2 = 1

        # Bounds on Euclidean norm of the witness
        if CCA:
            self.bound_witness_1 = sqrt(pke.B_se_s + samp.w + (1 if bimodal else 0))
            self.bound_witness_2 = sqrt(samp.B_s + 2*pke.B_e1_s)
        else:
            self.bound_witness_1 = sqrt(pke.B_se_s + samp.w + (1 if bimodal else 0) + pke.B_e1_s)
            self.bound_witness_2 = samp.B
        self.bound_witness = sqrt(self.bound_witness_1**2 + self.bound_witness_2**2)

        # Bound on j = 1/q * (v1 + [-h^T | a3]*[v2 | v3] + q_L*id*gH^T*v2 - u)
        self.B_j = 1/samp.q * ( \
                        samp.B_1 + \
                        samp.q/2 * sqrt(samp.n * (samp.n * (1 + samp.k_H))) * samp.B_2 + \
                        samp.q_L * samp.w * sqrt((samp.q_H**2 - 1)/(samp.b_H**2 - 1)) * samp.B_2 + \
                        samp.q/2 * sqrt(samp.n)
                    )

        # Bound on j_0 = 1/p * ([A_op | I_d]*[s | e] - ct_0)
        self.B_j_0 = 1/pke.p * (
                        sqrt(pke.n_op*pke.d) * sqrt(pke.p**2/4 * (pke.n_op*pke.d)**2 + 1) * pke.B_se + \
                        pke.p/2 * sqrt(pke.n_op*pke.d)
                    )

        # Bound on j_1 = 1/p * (b_op1.s + e1 + round(p/2).t - ct1)
        # and on j_2 = 1/p * (b_op2.s + e2 + round(p/2).t - ct2)
        self.B_j_1 = 1/pke.p * ( \
                            pke.p/2 * pke.n_op*sqrt(pke.d) * pke.B_se + \
                            pke.B_e1 + \
                            round(pke.p/2) * sqrt(samp.w) + \
                            pke.p/2 * sqrt(pke.n_op)
                        )

        # Rejection sampling parameters (hardcoded)
        self.M_11 = sqrt(2)
        self.M_12 = sqrt(2)
        self.M_2 = sqrt(2)
        self.M_3_l2 = sqrt(2)
        self.M_3_linf = sqrt(2)
        if bimodal:
            self.eps_11 = 0
            self.eps_12 = 0
            self.eps_2 = 0
            self.eps_3_l2 = 0
            self.eps_3_linf = 0
            self.alpha_11 = sqrt(pi/log(self.M_11))
            self.alpha_12 = sqrt(pi/log(self.M_12))
            self.alpha_2 = sqrt(pi/log(self.M_2))
            self.alpha_3_l2 = sqrt(pi/log(self.M_3_l2))
            self.alpha_3_linf = sqrt(pi/log(self.M_3_linf))
        else:
            self.eps_11 = 2**(-130)
            self.eps_12 = 2**(-130)
            self.eps_2 = 2**(-130)
            self.eps_3_l2 = 2**(-130)
            self.eps_3_linf = 2**(-130)
            self.alpha_11 = sqrt(pi)/log(self.M_11) * (sqrt(log(1/self.eps_11) + log(self.M_11)) + sqrt(log(1/self.eps_11)))
            self.alpha_12 = sqrt(pi)/log(self.M_12) * (sqrt(log(1/self.eps_12) + log(self.M_12)) + sqrt(log(1/self.eps_12)))
            self.alpha_2 = sqrt(pi)/log(self.M_2) * (sqrt(log(1/self.eps_2) + log(self.M_2)) + sqrt(log(1/self.eps_2)))
            self.alpha_3_l2 = sqrt(pi)/log(self.M_3_l2) * (sqrt(log(1/self.eps_3_l2) + log(self.M_3_l2)) + sqrt(log(1/self.eps_3_l2)))
            self.alpha_3_linf = sqrt(pi)/log(self.M_3_linf) * (sqrt(log(1/self.eps_3_linf) + log(self.M_3_linf)) + sqrt(log(1/self.eps_3_linf)))

        # Gaussian width
        self.s_11 = self.alpha_11 * self.eta * self.bound_witness_1
        self.s_12 = self.alpha_12 * self.eta * self.bound_witness_2
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n_pi * self.m_2)
        if infARP:
            self.s_3_l2 = self.alpha_3_l2 * sqrt(337) * sqrt(self.bound_witness**2)
            self.s_3_linf = self.alpha_3_linf * sqrt(337) * sqrt(self.B_j**2 + self.B_j_0**2 + (2 if CCA else 1) * self.B_j_1**2)
        else:
            self.s_3_l2 = self.alpha_3_l2 * sqrt(337) * sqrt(self.bound_witness**2 + self.B_j**2 + self.B_j_0**2 + (2 if CCA else 1) * self.B_j_1**2)
            self.s_3_linf = 1

        # Checking approximate range proofs bounds
        self.B_z3_l2 = sqrt(floor(c_star(256, self.sec + 3) ** 2 * self.s_3_l2 ** 2 * 256))
        self.B_z3_linf = ceil(t_star(256, self.sec + 3) * self.s_3_linf)

        self.B_arp_l2 = self.B_z3_l2 * 2/sqrt(26)
        self.B_arp_linf = self.B_z3_linf * 2

        self.cond_1_bound = max(samp.B_s, pke.B_se_s, pke.B_e1_s, samp.w) # Lower bound -q < -B^2
        self.cond_2_bound = self.B_arp_l2**2 - min(samp.B_s, pke.B_se_s, pke.B_e1_s, samp.w) # Upper bound (2/root(26) * B_e)^2 - B^2 < q
        self.cond_3_bound = 41 * self.n_pi * (self.m_11 + self.m_12 + (0 if infARP else self.k_pi + self.k_op_pi*(pke.d + (2 if CCA else 1)))) * self.B_arp_l2 # Condition for modular JL bound for ARP bound
        self.cond_4_bound = samp.w + sqrt(samp.w*self.n_pi) # bound for proving tag binary

        # Bound on INFINITY norm of extracted lifting equation [1 | -h^T | a3]v* + qL.id*.gH^T.v2* - u - q.j* = 0 mod q_proof
        self.cond_j_bound   = sqrt(1 + samp.q**2/4 * samp.n * (1+samp.k_H)) * samp.B + \
                                sqrt(samp.w) * samp.q_L * sqrt((samp.q_H**2 - 1)/(samp.b_H**2 - 1)) * samp.B + \
                                samp.q/2 + \
                                samp.q * (self.B_arp_linf if infARP else self.B_arp_l2)

        # Bound on INFINITY norm of extracted lifting equation e* + A.s* - ct0 - p.j_0* = 0 mod q_proof
        self.cond_j_0_bound = sqrt(1 + (pke.p/2)**2 * pke.n_op*pke.d) * pke.B_se + \
                            pke.p / 2 + \
                            pke.p * (self.B_arp_linf if infARP else self.B_arp_l2)

        # Bound on INFINITY norm of extracted lifting equation b_op1.s + e1 + round(p/2).t - ct1 - p.j_1* = 0 mod q_proof
        # and on extracted lifting equation b_op2.s + e2 + round(p/2).t - ct2 - p.j_2* = 0 mod q_proof
        self.cond_j_1_bound = pke.p/2 * sqrt(pke.n_op*pke.d) * pke.B_se + \
                            pke.B_e1 + \
                            round(pke.p/2) + \
                            pke.p/2 + \
                            pke.p * (self.B_arp_linf if infARP else self.B_arp_l2)

        assert (self.q_pi > self.cond_1_bound), "ZKP modulus too small"
        assert (self.q_pi > self.cond_2_bound), "ZKP modulus too small"
        assert (self.q_pi > self.cond_3_bound), "ZKP modulus too small"
        assert (self.q_pi > self.cond_4_bound), "ZKP modulus too small"
        assert (self.q_pi > self.cond_j_bound), "ZKP modulus too small"
        assert (self.q_pi > self.cond_j_0_bound), "ZKP modulus too small"
        assert (self.q_pi > self.cond_j_1_bound), "ZKP modulus too small"

        # Square Verification bounds
        self.B_11_s = 4 * floor(c_star(self.m_11 * self.n_pi, self.sec + 3) ** 2 * self.s_11 ** 2 * (self.m_11 * self.n_pi))
        self.B_12_s = 4 * floor(c_star(self.m_12 * self.n_pi, self.sec + 3) ** 2 * self.s_12 ** 2 * (self.m_12 * self.n_pi))
        self.B_1_s = self.B_11_s + self.B_12_s
        if self.D != 0 and self.gamma != 0:
            self.B_2_s = floor((2 * sqrt(floor(c_star(self.m_2 * self.n_pi, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n_pi))) + (2**self.D * self.eta + self.gamma)*sqrt(self.n_pi*self.d_pi)) ** 2)
        else:
            self.B_2_s = 4 * floor(c_star(self.m_2 * self.n_pi, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n_pi))

        # M-SIS bound
        self.msis_beta = 4 * self.eta * sqrt(self.B_1_s + self.B_2_s)

        ### Security
        self.pre_soundness_error = self.q_min ** (-self.n_pi/2) + 2 / self.challenge_space_size 
        if garbage:
            self.pre_soundness_error += self.q_min ** (-2*self.l)
        else:
            self.pre_soundness_error += self.q_min ** (-self.l)

        ### Efficiency

        # CRS size
        ajtai_crs_size = self.d_pi * (self.m_1 + self.m_2) * (self.n_pi * ceil(log2(self.q_pi)))
        bdlop_crs_size = ((2 if self.infARP else 1) * 256 / self.n_pi + self.l + 1 + 1) * self.m_2 * (self.n_pi * ceil(log2(self.q_pi)))
        self.crs_size  = ajtai_crs_size + bdlop_crs_size

        # Proof size
        self.size_z11 = ceil(self.n_pi * self.m_11 * (1/2 + log2(self.s_11)))
        self.size_z12 = ceil(self.n_pi * self.m_12 * (1/2 + log2(self.s_12)))
        self.size_z2 = ceil(self.n_pi * (self.m_2-(self.d_pi if self.D != 0 else 0)) * (1/2 + log2(self.s_2)))
        self.size_z3_l2 = ceil(256 * (1/2 + log2(self.s_3_l2)))
        self.size_z3_linf = ceil(256 * (1/2 + log2(self.s_3_linf)))
        self.size_c  = self.n_pi * ceil(log2(2 * self.rho + 1))
        self.size_tA = self.n_pi * self.d_pi * (ceil(log2(self.q_pi)) - self.D)
        self.size_hints = ceil(self.n_pi * self.d_pi * (2.25 if self.D != 0 else 0))
        self.size_tB = self.n_pi * (256 / self.n_pi + self.l) * ceil(log2(self.q_pi))
        self.size_h  = self.n_pi * self.l * ceil(log2(self.q_pi))
        self.size_t1 = self.n_pi * ceil(log2(self.q_pi))

        # Size of incompressible elements (those uniform in R_q)
        self.incompressible_bitsize = self.size_tA + self.size_hints + self.size_tB + self.size_h + self.size_t1

        # Size of compressible elements (those Gaussians): 
        self.compressible_bitsize = self.size_z11 + self.size_z12 + self.size_z2 + self.size_z3_l2 + (self.size_z3_linf if infARP else 0)

        # Size of proof
        self.proof_bitsize = self.incompressible_bitsize + self.compressible_bitsize + self.size_c
    
    def __repr__(self):
        """
        Printing a ZKP_Parameters object
        """
        tmp = '\n[+] ZERO-KNOWLEDGE PROOF PARAMETERS\n'
        tmp += 105 * '=' + '\n'
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Subring degree of Z[X]/(X^n_π + 1)', 'n_π', self.n_pi)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Subring gap factor (n / n_π)', 'k_π', self.k_pi)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Subring gap factor (n_op / n_π)', 'k_{op,π}', self.k_op_pi)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Module rank', 'd_2', self.d_pi)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Smallest modulus factor', 'q_min', self.q_min)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Modulus', 'q_pi', self.q_pi)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Parameter for soundness amplification', 'l', self.l)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Dimension of witness', 'm_1', self.m_1)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Dimension of witness (first part)', 'm_11', self.m_11)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Dimension of witness (second part)', 'm_12', self.m_12)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Dimension of ABDLOP commitment randomness', 'm_2', self.m_2)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Infinity norm of commitment randomness', 'ξ', self.xi_s2)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Infinity norm of challenges', 'ρ', self.rho)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Manhattan-like norm of challenges', 'η', self.eta)
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Size of challenge space', '|C|', floor(log2(self.challenge_space_size)))        
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Compression parameter 1', 'gamma', self.gamma)
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('Compression parameter 2', 'D', self.D)        
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 11', 'M_11', self.M_11)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 12', 'M_12', self.M_12)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 2', 'M_2', self.M_2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 3 (l2 ARP)', 'M_3_l2', self.M_3_l2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 3 (linf ARP)', 'M_3_linf', self.M_3_linf) if self.infARP else ''
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Overall repetition rate', 'M', self.M_11*self.M_12*self.M_2*self.M_3_l2*(self.M_3_linf if self.infARP else 1))
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Rejection sampling loss 11', 'ε_11', 0 if self.eps_11 == 0 else int(log2(self.eps_11)))
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Rejection sampling loss 12', 'ε_12', 0 if self.eps_12 == 0 else int(log2(self.eps_12)))
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Rejection sampling loss 2', 'ε_2', 0 if self.eps_2 == 0 else int(log2(self.eps_2)))
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Rejection sampling loss 3 (l2)', 'ε_3_l2', 0 if self.eps_3_l2 == 0 else int(log2(self.eps_3_l2)))
        tmp += '| {:60s} | {:^15s} | 2^{:<18d} |\n'.format('Rejection sampling loss 3 (linf)', 'ε_3_linf', 0 if self.eps_3_linf == 0 else int(log2(self.eps_3_linf))) if self.infARP else ''
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling slack 11', 'α_11', self.alpha_11)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling slack 12', 'α_12', self.alpha_12)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling slack 2', 'α_2', self.alpha_2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling slack 3 (l2)', 'α_3_l2', self.alpha_3_l2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Rejection sampling slack 3 (linf)', 'α_3_linf', self.alpha_3_linf) if self.infARP else ''
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Gaussian width for y_11', 'σ_11', self.s_11)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Gaussian width for y_12', 'σ_12', self.s_12)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Gaussian width for y_2', 'σ_2', self.s_2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Gaussian width for y_3 (l2)', 'σ_3_l2', self.s_3_l2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Gaussian width for y_3 (linf)', 'σ_3_linf', self.s_3_linf) if self.infARP else ''
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness 1', 'B1', self.bound_witness_1)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness 2', 'B2', self.bound_witness_2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness', 'B', self.bound_witness)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Euclidean norm bound on the lifting quotient j', 'B_j', self.B_j)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Euclidean norm bound on the lifting quotient j_0', 'B_j_0', self.B_j_0)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Euclidean norm bound on the lifting quotient j_1/j_2', 'B_j_1', self.B_j_1)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Norm bound from approximate range proof (l2)', 'B_arp_l2', self.B_arp_l2)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Norm bound from approximate range proof (linf)', 'B_arp_linf', self.B_arp_linf) if self.infARP else ''
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 1: B_i^2', 'Cond. 1', floor(self.cond_1_bound))
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 2: B_arp^2 - B_i^2', 'Cond. 2', floor(self.cond_2_bound))
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 3: 41*n*m_1*B_arp', 'Cond. 3', floor(self.cond_3_bound))
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 4: for binary tag', 'Cond. 4', floor(self.cond_4_bound))
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 6: for mod-q lifting', 'j', floor(self.cond_j_bound))
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 7: for mod-p lifting (ct0)', 'j_0', floor(self.cond_j_0_bound))
        tmp += '| {:60s} | {:^15s} | {:<20d} |\n'.format('ARP bound condition 7: for mod-p lifting (ct1/ct2)', 'j_1', floor(self.cond_j_1_bound))
        tmp += '| {:60s} | {:^15s} | 2^{:<18.5f} |\n'.format('(Pre) Soundness error', 'δ_s', log2(self.pre_soundness_error))
        tmp += '\n[+] Proof Estimated Performance (KB)\n'
        tmp += 105 * '=' + '\n'
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Common Random String (KB)', '|crs|', self.crs_size / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Incompressible elements (KB)', '|π_1|', self.incompressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Compressible elements (KB)', '|π_2|', self.compressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Challenge (KB)', '|π_3|', self.size_c / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Total proof bitsize (KB)', '|π|', self.proof_bitsize / 2 ** 13.)
        tmp += 105 * '=' + '\n'
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of z11 (KB)', '|z_11|', self.size_z11 / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of z12 (KB)', '|z_12|', self.size_z12 / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of z11,z12 (KB)', '|z_1|', (self.size_z11 + self.size_z12) / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of z2 (KB)', '|z_2|', self.size_z2 / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of z3_l2 (KB)', '|z_3_l2|', self.size_z3_l2 / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of z3_linf (KB)', '|z_3_linf|', self.size_z3_linf / 2 ** 13.) if self.infARP else ''
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of c (KB)', '|c|', self.size_c / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of tA (KB)', '|tA|', self.size_tA / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of compression hints (KB)', '|hints|', self.size_hints / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of tB (KB)', '|tB|', self.size_tB / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of h (KB)', '|h|', self.size_h / 2 ** 13.)
        tmp += '| {:60s} | {:^15s} | {:<20.5f} |\n'.format('Size of t1 (KB)', '|t1|', self.size_t1 / 2 ** 13.)

        return tmp

def estimate_group_signature(samp, pke, zkp, no_guessing=False, CCA=True):
    """
    Estimate concrete security and efficiency of the group signature based on chosen parameters
    - input:
        (NTRUSampler_Parameters) samp -- Sampler parameters
        (PKE_Parameters) pke -- Verifiable encryption parameters
        (ZKP_Parameters) zkp -- Group signature zero-knowledge proof parameters
    """
    # Print all parameters
    print(samp)
    print(pke)
    print(zkp)

    # Compute group signature sizes
    byte = lambda x: ceil(x / 2**3)
    kbyte = lambda x: x / 2**13.
    gpk_sz = samp.ipk_bitsize + pke.opk_bitsize
    gmsk_sz = pke.osk_bitsize
    gs_sz = pke.ct_bitsize + zkp.proof_bitsize - (0 if CCA else pke.ct2_bitsize)

    tmp = '\n[KEY SIZES]\n'
    tmp += '\n'
    tmp += '    |gpk|    =  %4d B    ( %1.3f KB)  [group public key]\n' % (byte(gpk_sz), kbyte(gpk_sz))
    tmp += '    |gmsk|   =  %4d B    ( %1.3f KB)  [group manager secret key]\n' % (byte(gmsk_sz), kbyte(gmsk_sz))
    tmp += '    |gsk[i]| =  %4d B    ( %1.3f KB)  [user secret key]\n' % (byte(samp.gsk_i_bitsize), kbyte(samp.gsk_i_bitsize))
    tmp += '\n'
    tmp += '\n[GROUP SIGNATURE SIZE]\n'
    tmp += '\n'
    tmp += '  + |ct_0|  =  %5d B    ( %2.3f KB)\n' % (byte(pke.ct0_bitsize), kbyte(pke.ct0_bitsize))
    tmp += '  + |ct_1|  =  %5d B    ( %2.3f KB)\n' % (byte(pke.ct1_bitsize), kbyte(pke.ct1_bitsize))
    tmp += ('  + |ct_2|  =  %5d B    ( %2.3f KB)\n' % (byte(pke.ct2_bitsize), kbyte(pke.ct2_bitsize)) if CCA else '')
    tmp += '  + |π|     =  %5d B    (%2.3f KB)\n' % (byte(zkp.proof_bitsize), kbyte(zkp.proof_bitsize))
    tmp += '  _________________________\n'
    tmp += '  =  %5d B    (%2.3f KB)\n' % (byte(gs_sz), kbyte(gs_sz))
    print(tmp)

    # Compute security of all assumptions
    print("\n[+] START LOG [+]")
    ## M-LWE assumptions
    mlwe_op = estimate_LWE(n=pke.n_op*pke.d, q=pke.p, Xs=ND.CenteredBinomial(1), Xe=ND.CenteredBinomial(1), m=pke.n_op*pke.d, cost_model=COST_MODEL, rough=ROUGH)
    
    if CCA:
        mlwe_for_matrix_hint = estimate_LWE(n=pke.n_op*pke.d, q=pke.p, Xs=ND.DiscreteGaussian(pke.s_MLWE/sqrt(2*pi)), Xe=ND.DiscreteGaussian(pke.s_MLWE/sqrt(2*pi)), m=pke.n_op*(pke.d+1), cost_model=COST_MODEL, rough=ROUGH)
        mlwe_ct_an = estimate_LWE(n=pke.n_op*pke.d, q=pke.p, Xs=ND.DiscreteGaussian(pke.s_se/sqrt(2*pi)), Xe=ND.DiscreteGaussian(pke.s_se/sqrt(2*pi)), m=pke.n_op*(pke.d+1), cost_model=COST_MODEL, rough=ROUGH)
    else:
        mlwe_ct_an = estimate_LWE(n=pke.n_op*pke.d, q=pke.p, Xs=ND.CenteredBinomial(1), Xe=ND.CenteredBinomial(1), m=pke.n_op*(pke.d+1), cost_model=COST_MODEL, rough=ROUGH)

    ## NTRU assumptions
    ntru = estimate_NTRU(n=samp.n, q=samp.q, Xs=ND.CenteredBinomial(1), Xe=ND.CenteredBinomial(1), m=samp.n, ntru_type="circulant", cost_model=COST_MODEL, rough=ROUGH)

    ## RISIS assumption for T-NTRU-ISIS
    risis = estimate_ISIS(n=samp.n, m=samp.n*(samp.k_H + 2), q=samp.q, beta=samp.B, cost_model=COST_MODEL, estimator=False, rough=ROUGH)

    ## ZKP assumptions
    mlwe_zk = estimate_LWE(n=zkp.n_pi * (zkp.m_2 - (zkp.d_pi + floor(256/zkp.n_pi) + (floor(256/zkp.n_pi) if zkp.infARP else 0) + zkp.l + 1)), q=zkp.q_pi, Xs=ND.CenteredBinomial(zkp.xi_s2), Xe=ND.CenteredBinomial(zkp.xi_s2), m=zkp.n_pi * zkp.m_2, cost_model=COST_MODEL, rough=ROUGH)
    
    msis_zk = estimate_SIS(n=zkp.n_pi * zkp.d_pi, m=zkp.n_pi * (zkp.m_1 + zkp.m_2), q=zkp.q_pi, beta=zkp.msis_beta, cost_model=COST_MODEL, estimator=False, rough=ROUGH)
    print("[+] END LOG [+]\n")

    # Print hardness
    print_sis = lambda name,res,beta,smaller: print("\t[*] Hardness of %s \t>>\tBKZ-β = %d \t CSec = %.3f \t QSec = %.3f \t Bound β = %.2f \t[β < q: %s]" % (name,res[2],res[3],res[4],beta,smaller))
    print_assumption = lambda name,res: print("\t[*] Hardness of %s \t>>\tBKZ-β = %d \t CSec = %.3f \t QSec = %.3f" % (name,res[0],res[2],res[3]))

    print_assumption("M-LWE (op)                ", mlwe_op)
    print_assumption("M-LWE (ct, anonymity)     ", mlwe_ct_an)
    if CCA: print_assumption("M-LWE (for matrix hint)   ", mlwe_for_matrix_hint)
    print_assumption("NTRU (trapdoor switching)", ntru)
    print_sis("R-ISIS (for T-NTRU-ISIS)  ", risis, samp.B, samp.B<samp.q)
    print_sis("M-SIS (π)                 ", msis_zk, zkp.msis_beta, zkp.msis_beta<zkp.q_pi)
    print_assumption("M-LWE (π)                 ", mlwe_zk)
    print()

    # Computing hardness bounds
    hardness_bound  = lambda res: 2**(-res[-2]) # Classical bit security
    e_mlwe_op       = hardness_bound(mlwe_op)
    e_mlwe_ct_an    = hardness_bound(mlwe_ct_an)
    if CCA:
        e_mlwe_for_matrix_hint = hardness_bound(mlwe_for_matrix_hint)
        e_mhmlwe      = e_mlwe_for_matrix_hint + pke.eps_op/(1 - pke.eps_op)
        print("\t[o] Hardness of MH-M-LWE                  \t >> Csec: %.3f" % (-log2(e_mhmlwe)))
    else:
        e_mhmlwe      = 2**(-200)
    e_ntru          = hardness_bound(ntru)
    e_risis         = hardness_bound(risis)
    e_tntruisis     = samp.k_H * e_ntru + e_risis
    print("\t[o] Hardness of T-NTRU-ISIS               \t >> Csec: %.3f" % (-log2(e_tntruisis)))
    e_mlwe_zk       = hardness_bound(mlwe_zk)
    e_msis_zk       = hardness_bound(msis_zk)

    e_snd = zkp.pre_soundness_error + e_msis_zk
    e_zk = zkp.eps_11/zkp.M_11 + zkp.eps_12/zkp.M_12 + zkp.eps_2/zkp.M_2 + zkp.eps_3_l2/zkp.M_3_l2 + (zkp.eps_3_linf/zkp.M_3_linf if zkp.infARP else 0) + e_mlwe_zk

    # Computing loss factors
    a = 1 - 1/(2 * samp.sec)
    o = 2 * samp.sec

    ## Trapdoor switching loss terms
    eps = 2 ** LOG2_EPS
    delta_gpk = ((1 + eps/(samp.k_H))/(1 - eps/(samp.k_H)))**(2*samp.k_H) * ((1+eps)/(1-eps))**(4*(1+samp.k_H)*(samp.n-1) + 8) - 1
    mu = (1 + o*(o-1)*delta_gpk**2/(2*(1-delta_gpk)**(o+1)) )**((samp.N-1)/o)

    if CCA:
        e_anonymity = e_snd + e_zk + 3*e_mlwe_op + 2*e_mhmlwe + e_mlwe_ct_an
    else:
        e_anonymity = e_zk + e_mlwe_op + e_mlwe_ct_an
    security_anonymity = -log2(e_anonymity)

    # Full Traceability
    h = lambda Adv: e_ntru + mu * (2*e_ntru + mu * (e_ntru + Adv)**a)**a 

    e_tmp = e_tntruisis
    for _ in range(samp.k_H):
        e_tmp = h(e_tmp)
    e_traceability = e_snd + (1 if no_guessing else samp.N) * (Q_SIGN*e_zk + e_tmp)
    security_traceability = -log2(e_traceability)

    tmp = '\n'
    tmp += ('[%s-ANONYMITY CONCRETE SECURITY] Achieved: %.3f\n' % ('CCA' if CCA else 'CPA', security_anonymity))
    tmp += '[TRACEABILITY CONCRETE SECURITY] Achieved: %.3f\n' % (security_traceability)

    print(tmp)

ROUGH            = False
ZK_OPTIMIZATIONS = True
CCA_version      = True

if CCA_version:
    samp_pms = NTRUSampler_Parameters(target_bitsec=128, n=1024, n_pi=64, b_H=13, k_H=2, q_L=733, N_min=2**32)
    pke_pms = CCA_PKE_Parameters(target_bitsec=128, n_op=256, d=3)
    # ZKP with ARP only in L2 norm
    zkp_pms = ZKP_Parameters(target_bitsec=128, n_pi=64, d_pi=16, m_2=78, q_start=round(2**(68.31)), gamma=2**26, D=18, samp=samp_pms, pke=pke_pms, bimodal=ZK_OPTIMIZATIONS, compression=ZK_OPTIMIZATIONS, garbage=ZK_OPTIMIZATIONS, CCA=CCA_version, infARP=False)
    # ZKP with ARP in both L2 and Linf norms
    # zkp_pms = ZKP_Parameters(target_bitsec=128, n_pi=64, d_pi=20, m_2=75, q_start=round(2**(52.51)), gamma=2**24, D=16, samp=samp_pms, pke=pke_pms, bimodal=ZK_OPTIMIZATIONS, compression=ZK_OPTIMIZATIONS, garbage=ZK_OPTIMIZATIONS, CCA=CCA_version, infARP=True)
else:
    samp_pms = NTRUSampler_Parameters(target_bitsec=128, n=1024, n_pi=64, b_H=13, k_H=2, q_L=733, N_min=2**32)
    pke_pms = CPA_PKE_Parameters(target_bitsec=128, n_op=256, d=3)
    # ZKP with ARP only in L2 norm
    # zkp_pms = ZKP_Parameters(target_bitsec=128, n_pi=64, d_pi=16, m_2=78, q_start=round(2**(68.31)), gamma=2**26, D=18, samp=samp_pms, pke=pke_pms, bimodal=ZK_OPTIMIZATIONS, compression=ZK_OPTIMIZATIONS, garbage=ZK_OPTIMIZATIONS, CCA=CCA_version, infARP=False)
    # ZKP with ARP in both L2 and Linf norms
    zkp_pms = ZKP_Parameters(target_bitsec=128, n_pi=64, d_pi=20, m_2=75, q_start=round(2**(52.51)), gamma=2**24, D=16, samp=samp_pms, pke=pke_pms, bimodal=ZK_OPTIMIZATIONS, compression=ZK_OPTIMIZATIONS, garbage=ZK_OPTIMIZATIONS, CCA=CCA_version, infARP=True)

estimate_group_signature(samp_pms, pke_pms, zkp_pms, no_guessing=True, CCA=CCA_version)