// RLWE
#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <hexl/hexl.hpp>
#include <cstdint>
#include <iostream>
#include <vector>
#include <pthread.h>
#ifndef PORTABLE_BUILD
#include <x86intrin.h>
#endif

// Note: RNS_* and RNSc_* must have the same memory layout

typedef struct _incNTT{
  intel::hexl::NTT ** ntt = NULL;
  uint64_t split_degree;
  uint64_t ** w;
  uint64_t *** D, *** Dhat;
} * incNTT;
 

typedef struct _RNS_Polynomial {
  uint64_t ** coeffs;
  incNTT ntt; 
  uint64_t N, l;
} * RNS_Polynomial;

/* RNS polynomial in coefficient representation*/
typedef struct _RNSc_Polynomial {
  uint64_t ** coeffs;
  incNTT ntt;
  uint64_t N, l;
} * RNSc_Polynomial;

typedef struct _ZqVector {
  uint64_t ** elements;
  uint64_t n, *q, l;
} * ZqVector;


/* RLWE RNS */

typedef struct _IntPolynomial
{
  uint64_t * coeffs;
  uint64_t N;
} * IntPolynomial;

typedef struct _RNS_RLWE_Key {
  IntPolynomial s;
  RNS_Polynomial s_RNS;
  uint64_t N, l, ** delta;
  double sigma;
} * RNS_RLWE_Key;

typedef struct _RNS_RLWE {
  RNS_Polynomial a, b;
} * RNS_RLWE;

typedef struct _RNSc_RLWE {
  RNSc_Polynomial a, b;
} * RNSc_RLWE;

typedef struct _RNS_RLWE_KS_Key {
  RNS_RLWE * s;
  uint64_t l;
} * RNS_RLWE_KS_Key;


#ifdef __cplusplus
extern "C" {
#endif

// polynomial
IntPolynomial polynomial_new_int_polynomial(uint64_t N);
RNS_Polynomial polynomial_new_RNS_polynomial(uint64_t N, uint64_t l, incNTT ntt);
void polynomial_RNS_zero(RNS_Polynomial p);
RNS_Polynomial * polynomial_new_array_of_RNS_polynomials(uint64_t N, uint64_t l, uint64_t size, incNTT ntt);
void polynomial_to_RNS(RNS_Polynomial out, IntPolynomial in);
void polynomial_gen_random_RNSc_polynomial(RNSc_Polynomial out);
void polynomial_mul_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2);
void polynomial_sub_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2);
void polynomial_sub_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2);
void polynomial_RNSc_to_RNS(RNS_Polynomial out, RNSc_Polynomial in);
void polynomial_RNS_to_RNSc(RNSc_Polynomial out, RNS_Polynomial in);
void polynomial_RNSc_add_noise(RNSc_Polynomial out, RNSc_Polynomial in, double sigma);
void polynomial_base_reduce_RNSc(RNSc_Polynomial out);
void polynomial_base_reduce_RNSc_wo_free(RNSc_Polynomial out);
void polynomial_base_extend_RNSc(RNSc_Polynomial out, uint64_t * in, uint64_t p);
// void polynomial_RNSc_reduce_mod_XQ_minus_1(RNSc_Polynomial out);
void polynomial_RNSc_permute(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t gen);
void free_RNS_polynomial(void * p);
void polynomial_RNSc_negate(RNSc_Polynomial out, RNSc_Polynomial in);
void polynomial_add_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2);
void polynomial_add_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2);
void polynomial_int_permute_mod_Q(IntPolynomial out, IntPolynomial in, uint64_t gen);
void polynomial_copy_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in);
void polynomial_reduce_mod_Xd_minus_1_mod_q(uint64_t * p, int64_t d, uint64_t q, size_t size);
void polynomial_reduce_mod_Xd_minus_1(uint64_t * p, int64_t d, size_t size);
void polynomial_RNSc_mul_by_xai(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t a);
void polynomial_int_decompose_i(IntPolynomial out, IntPolynomial in, uint64_t Bg_bit, uint64_t l, uint64_t q, uint64_t bit_size, uint64_t i);

void free_polynomial(void * p);
void array_to_RNS(RNS_Polynomial out, uint64_t ** in);
void polynomial_RNS_get_hash(uint64_t * out, RNS_Polynomial p);
uint64_t * polynomial_RNS_get_hash_p(RNS_Polynomial p);
void array_conj_to_RNS(RNS_Polynomial out, uint64_t * in);
void polynomial_base_reduce_RNSc_and_scale(RNSc_Polynomial out, uint64_t p);
RNS_Polynomial * polynomial_new_RNS_polynomial_array(uint64_t size, uint64_t N, uint64_t l, incNTT ntt);
void free_RNS_polynomial_array(uint64_t size, RNS_Polynomial * p);
void polynomial_scale_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t scale);
void polynomial_scale_RNS_polynomial_RNS(RNS_Polynomial out, RNS_Polynomial in1, uint64_t * scale);
void polynomial_base_reduce_round_RNSc_wo_free(RNSc_Polynomial out);
bool polynomial_eq(RNS_Polynomial a, RNS_Polynomial b);
void polynomial_base_extend_RNSc_2(RNSc_Polynomial out, RNSc_Polynomial in);
void polynomial_multo_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in);
void polynomial_RNSc_mod_reduce(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t idx);
void polynomial_RNS_broadcast_slot(RNS_Polynomial out, RNS_Polynomial in, uint64_t slot_idx);
void polynomial_RNS_rotate_slot(RNS_Polynomial out, RNS_Polynomial in, uint64_t rot);
void polynomial_RNS_copy_slot(RNS_Polynomial out, uint64_t dst, RNS_Polynomial in, uint64_t src);
void polynomial_gen_gaussian_RNSc_polynomial(RNSc_Polynomial out, double sigma);
void int_array_to_RNS(RNS_Polynomial out, uint64_t * in);
void polynomial_RNSc_add_integer(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t in2);
void polynomial_scale_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, uint64_t scale);
void polynomial_RNS_negate(RNS_Polynomial out, RNS_Polynomial in);
void polynomial_RNS_add_integer(RNS_Polynomial out, RNS_Polynomial in1, uint64_t in2);
void polynomial_scale_addto_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t scale);
void polynomial_scale_addto_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, uint64_t scale);

// rlwe rns
RNS_RLWE_Key rlwe_get_RNS_key_from_array(uint64_t N, uint64_t l, uint64_t * array, incNTT ntt, double sigma);
RNS_RLWE_Key rlwe_alloc_RNS_key(uint64_t N, uint64_t l, incNTT ntt, double sigma);
void free_RNS_rlwe_sample(RNS_RLWE c);
RNS_RLWE_Key rlwe_new_RNS_gaussian_key(uint64_t N, uint64_t l, double key_sigma, incNTT ntt, double sigma);
RNS_RLWE rlwe_alloc_RNS_sample(uint64_t N, uint64_t l, incNTT ntt);
RNSc_RLWE rlwe_alloc_RNSc_sample(uint64_t N, uint64_t l, incNTT ntt);
RNS_RLWE rlwe_new_RNS_sample(RNS_RLWE_Key key, uint64_t * m, uint64_t p);
void rlwe_RNS_sample_of_zero(RNS_RLWE out, RNS_RLWE_Key key);
void rlwe_RNSc_sample_of_zero(RNSc_RLWE out, RNS_RLWE_Key key);
RNS_RLWE rlwe_new_RNS_sample_of_zero(RNS_RLWE_Key key);
RNSc_RLWE rlwe_new_RNSc_sample_of_zero(RNS_RLWE_Key key);
RNS_RLWE rlwe_new_RNS_trivial_sample_of_zero(uint64_t N, uint64_t l, incNTT ntt);
void rlwe_RNS_phase(RNS_Polynomial out, RNS_RLWE in, RNS_RLWE_Key key);
void rlwe_RNSc_to_RNS(RNS_RLWE out, RNSc_RLWE in);
void rlwe_RNS_to_RNSc(RNSc_RLWE out, RNS_RLWE in);
RNSc_RLWE rlwe_new_RNSc_sample(RNS_RLWE_Key key, uint64_t m);
void rlwe_RNS_trivial_sample_of_zero(RNS_RLWE out);
void rlwe_copy_RNS_sample(RNS_RLWE out, RNS_RLWE in);
void rlwe_copy_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in);
void rlwe_shrink_RNSc_sample(RNSc_RLWE c);
void rlwe_RNSc_mul_by_xai(RNSc_RLWE out, RNSc_RLWE in, uint64_t a);
void rlwe_RNSc_encrypt_RNS(RNSc_RLWE out, RNS_RLWE_Key key, RNSc_Polynomial m);
void rlwe_RNSc_decrypt_RNS(RNSc_Polynomial out, RNSc_Polynomial tmp, RNS_RLWE in, RNS_RLWE_Key key);
void rlwe_scale_RNS_rlwe_RNS(RNS_RLWE c, uint64_t * scale);
void rlwe_add_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in1, RNSc_RLWE in2);
void rlwe_sub_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in1, RNSc_RLWE in2);

void rlwe_scale_RNSc_rlwe(RNSc_RLWE c, uint64_t scale);
void rlwe_RNSc_mod_switch(RNSc_RLWE c, uint64_t q);
void rlwe_addto_RNSc_sample(RNSc_RLWE out, RNSc_RLWE in);
RNS_RLWE * rlwe_alloc_RNS_sample_array(uint64_t size, uint64_t N, uint64_t l, incNTT ntt);
RNS_RLWE * rlwe_alloc_RNS_sample_array2(uint64_t size, RNS_RLWE c);
void free_RNS_rlwe_array(uint64_t size, RNS_RLWE * v);
void free_rlwe_RNS_sample(void * p);
void rlwe_scale_RNS_rlwe_addto(RNS_RLWE out, RNS_RLWE in, uint64_t scale);
void rlwe_RNS_mul_by_poly(RNS_RLWE out, RNS_RLWE in, RNS_Polynomial poly);
void rlwe_RNSc_extract_lwe(uint64_t * out, RNSc_RLWE in, uint64_t idx);
void rlwe_add_RNSc_polynomial(RNSc_RLWE out, RNSc_RLWE in1, RNSc_Polynomial in2);
void rlwe_sub_RNSc_polynomial(RNSc_RLWE out, RNSc_RLWE in1, RNSc_Polynomial in2);
void rlwe_RNSc_encrypt_int(RNSc_RLWE out, RNS_RLWE_Key key, uint64_t m, uint64_t l);

// rlwe ks
RNS_RLWE_KS_Key rlwe_new_RNS_ks_key(RNS_RLWE_Key out_key, RNS_RLWE_Key in_key);
void rlwe_RNSc_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key ksk);
RNS_RLWE_KS_Key rlwe_new_RNS_automorphism_key(RNS_RLWE_Key key, uint64_t gen);
void rlwe_automorphism_RNSc(RNSc_RLWE out, RNSc_RLWE in, uint64_t gen, RNS_RLWE_KS_Key ksk);
void rlwe_RNSc_priv_keyswitch(RNSc_RLWE out, RNSc_RLWE in, RNS_RLWE_KS_Key kska, RNS_RLWE_KS_Key kskb);

// misc
void gen_sparse_ternary_array_modq(uint64_t * out, uint64_t size, uint64_t h, uint64_t q);
uint64_t next_power_of_2(uint64_t x);
void array_reduce_mod_N(uint64_t * out,  uint64_t * in, uint64_t size, uint64_t p);
void array_mod_switch(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n);
void array_mod_switch_from_2k(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n);
uint64_t int_mod_switch(uint64_t in, uint64_t p, uint64_t q);
intel::hexl::NTT ** new_ntt_list(uint64_t * primes, uint64_t N, uint64_t l);
incNTT new_incomplete_ntt_list(uint64_t * primes, uint64_t split_degree, uint64_t N, uint64_t l);
uint64_t ** incNTT_get_rou_matrix(incNTT ntt);
uint64_t double2int(double x);
uint64_t __debug_get_exp_message_from_noisy_RNS(RNS_Polynomial in, uint64_t * p);
void compute_RNS_Qhat_array(uint64_t * out, uint64_t * p, uint64_t l);
uint64_t __debug_get_exp_message_from_noisy_RNSc(RNSc_Polynomial in, uint64_t * p);
uint64_t __debug_get_exp_message_from_noisy_RNS_2(RNS_Polynomial in, uint64_t * p, uint64_t * Q_hat);
void array_additive_inverse_mod_switch(uint64_t * out, uint64_t * in, uint64_t p, uint64_t q, uint64_t n);
uint64_t mod_dist(uint64_t a, uint64_t b, uint64_t q);
void print_array(const char * msg, uint64_t * v, size_t size);
uint64_t mod_switch(uint64_t v, uint64_t p, uint64_t q);
unsigned char char_rev(unsigned char b);
uint32_t int_rev(uint32_t b);
void bit_rev(uint64_t * out, uint64_t * in, uint64_t n, uint64_t log_n);

// Misc from third party
void generate_random_bytes(uint64_t amount, uint8_t * pointer);
double generate_normal_random(double sigma);
void * safe_malloc(size_t size);
void * safe_aligned_malloc(size_t size);
#ifdef __cplusplus
}
#endif