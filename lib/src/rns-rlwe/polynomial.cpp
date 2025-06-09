#include "rns-rlwe.h"
#include <blake3.h>


incNTT new_incomplete_ntt_list(uint64_t * primes, uint64_t split_degree, uint64_t N, uint64_t l){
  const uint64_t poly_size = N/split_degree;
  const uint64_t log_poly_size = (uint64_t) log2(poly_size);
  uint64_t w_p[2*poly_size];
  incNTT ntt = (incNTT) safe_malloc(sizeof(*ntt));
  ntt->ntt = new_ntt_list(primes, poly_size, l);
  ntt->split_degree = split_degree;
  ntt->w = (uint64_t **) safe_malloc(sizeof(uint64_t *)*l);
  for (size_t i = 0; i < l; i++){
    ntt->w[i] = (uint64_t *) safe_aligned_malloc(poly_size*sizeof(uint64_t));
    uint64_t w1 = ntt->ntt[i]->GetMinimalRootOfUnity();
    uint64_t p = ntt->ntt[i]->GetModulus();
    w_p[0] = w1;
    for (size_t j = 1; j < 2*poly_size; j++){
      w_p[j] = intel::hexl::MultiplyMod(w_p[j - 1], w1, p);
    }
    bit_rev(ntt->w[i], w_p, poly_size, log_poly_size + 1);
  }
  // compute D and Dhat
  ntt->Dhat = (uint64_t***) safe_malloc(sizeof(uint64_t*)*l);
  ntt->D = (uint64_t***) safe_malloc(sizeof(uint64_t*)*l);
  for (size_t i = 0; i < l; i++){
    ntt->Dhat[i] = (uint64_t**) safe_aligned_malloc(sizeof(uint64_t*)*l);
    ntt->D[i] = (uint64_t**) safe_aligned_malloc(sizeof(uint64_t*)*l);
    for (size_t i2 = 0; i2 < l; i2++){
      ntt->Dhat[i][i2] = (uint64_t*) safe_aligned_malloc(sizeof(uint64_t)*l);
      ntt->D[i][i2] = (uint64_t*) safe_aligned_malloc(sizeof(uint64_t)*l);
      for (size_t j = 0; j < l; j++){
        const uint64_t q = ntt->ntt[j]->GetModulus();
      
        ntt->D[i][i2][j] = 1;
        for (size_t k = 0; k < i; k++){
          const uint64_t p = ntt->ntt[k]->GetModulus();
          if(i2 != k){
            ntt->D[i][i2][j] = intel::hexl::MultiplyMod(ntt->D[i][i2][j], p, q);
          }
        }
        ntt->Dhat[i][i2][j] = intel::hexl::InverseMod(ntt->D[i][i2][j], q);
      }
    }
    
  }
  return ntt;
}

uint64_t ** incNTT_get_rou_matrix(incNTT ntt){
  return ntt->w;
}

RNS_Polynomial polynomial_new_RNS_polynomial(uint64_t N, uint64_t l, incNTT ntt){
  RNS_Polynomial res;
  res = (RNS_Polynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t **) safe_malloc(sizeof(uint64_t*) * l);
  for (size_t i = 0; i < l; i++){
    res->coeffs[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t) * N);
  }
  res->N = N;
  res->l = l;
  res->ntt = ntt;
  return res;
}

RNS_Polynomial * polynomial_new_RNS_polynomial_array(uint64_t size, uint64_t N, uint64_t l, incNTT ntt){
  RNS_Polynomial * res;
  res = (RNS_Polynomial *) safe_malloc(sizeof(RNS_Polynomial)*size);
  for (size_t i = 0; i < size; i++){
    res[i] = polynomial_new_RNS_polynomial(N, l, ntt);
  }
  return res;
}

bool polynomial_eq(RNS_Polynomial a, RNS_Polynomial b){
  for (size_t i = 0; i < a->l; i++){
    if(memcmp(a->coeffs[i], b->coeffs[i], a->N*sizeof(uint64_t)) != 0){
      return false;
    }
  }
  return true;
}

void polynomial_copy_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in){
  for (size_t i = 0; i < in->l; i++){
    memcpy(out->coeffs[i], in->coeffs[i], sizeof(uint64_t)*out->N);
  }
}

void polynomial_RNS_zero(RNS_Polynomial p){
  for (size_t i = 0; i < p->l; i++){
    memset(p->coeffs[i], 0, sizeof(uint64_t)*p->N);
  }
}

void free_RNS_polynomial(void * p){
  RNS_Polynomial pp = (RNS_Polynomial) p;
  for (size_t i = 0; i < pp->l; i++){
    free(pp->coeffs[i]);
  }
  free(pp->coeffs);
  free(pp);
}

void free_RNS_polynomial_array(uint64_t size, RNS_Polynomial * p){
  for (size_t i = 0; i < size; i++){
    free_RNS_polynomial(p[i]);
  }
  free(p);
}

RNS_Polynomial * polynomial_new_array_of_RNS_polynomials(uint64_t N, uint64_t l, uint64_t size, incNTT ntt){
  RNS_Polynomial * res = (RNS_Polynomial *) safe_malloc(sizeof(RNS_Polynomial)*size);
  for (size_t i = 0; i < size; i++) res[i] = polynomial_new_RNS_polynomial(N, l, ntt);
  return res;
}

// out = RNS(in)
// Assumes ||in||_inf < min(p)^2
void polynomial_to_RNS(RNS_Polynomial out, IntPolynomial in){
  const uint64_t modMask = out->ntt->split_degree - 1, poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][(j&modMask)*poly_size + j/out->ntt->split_degree] = (p + in->coeffs[j])%p;
    }
  }
  polynomial_RNSc_to_RNS(out, (RNSc_Polynomial) out);
}

void int_array_to_RNS(RNS_Polynomial out, uint64_t * in){
  const uint64_t modMask = out->ntt->split_degree - 1, poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][(j&modMask)*poly_size + j/out->ntt->split_degree] = (p + in[j])%p;
    }
  }
  polynomial_RNSc_to_RNS(out, (RNSc_Polynomial) out);
}

void array_to_RNS(RNS_Polynomial out, uint64_t ** in){
  const uint64_t modMask = out->ntt->split_degree - 1, poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < out->l; i++){
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][(j&modMask)*poly_size + j/out->ntt->split_degree] = in[i][j];
    }
  }
  polynomial_RNSc_to_RNS(out, (RNSc_Polynomial) out);
}

// void array_conj_to_RNS(RNS_Polynomial out, uint64_t * in){
//   for (size_t i = 0; i < out->l; i++){
//     const uint64_t p = out->ntt[i]->GetModulus();
//     out->coeffs[i][0] = (p + in[0])%p;
//     for (size_t j = 1; j < out->N; j++){
//       out->coeffs[i][j] = (p - in[out->N - j])%p;
//     }
//     out->ntt[i]->ComputeForward(out->coeffs[i], out->coeffs[i], 1, 1);
//   }
// }

void polynomial_gen_random_RNSc_polynomial(RNSc_Polynomial out){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt->ntt[i]->GetModulus();
    generate_random_bytes(sizeof(uint64_t)*out->N, (uint8_t *) out->coeffs[i]);
    array_mod_switch_from_2k(out->coeffs[i], out->coeffs[i], p, p, out->N);
  }
}

void polynomial_gen_gaussian_RNSc_polynomial(RNSc_Polynomial out, double sigma){
  for (size_t j = 0; j < out->N; j++){
    const double real_noise = generate_normal_random(sigma);
    const int64_t rounded_noise = (int64_t) round(real_noise);
    for (size_t i = 0; i < out->l; i++){
      const uint64_t p = out->ntt->ntt[i]->GetModulus();
      const uint64_t noise = rounded_noise < 0? rounded_noise + p : rounded_noise;
      out->coeffs[i][j] = noise;
    }
  }
}


void polynomial_gen_random_small_RNSc_polynomial(RNSc_Polynomial out, uint64_t p_idx){
  const uint64_t p = out->ntt->ntt[p_idx]->GetModulus();
  generate_random_bytes(sizeof(uint64_t)*out->N, (uint8_t *) out->coeffs[p_idx]);
  array_mod_switch_from_2k(out->coeffs[p_idx], out->coeffs[p_idx], p, p, out->N);
  for (size_t i = 0; i < out->l; i++){
    const uint64_t p = out->ntt->ntt[i]->GetModulus();
    generate_random_bytes(sizeof(uint64_t)*out->N, (uint8_t *) out->coeffs[i]);
    array_mod_switch_from_2k(out->coeffs[i], out->coeffs[i], p, p, out->N);
  }
}


/* out = in1*in2 mod (X^N + 1) */
void polynomial_multo_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in){
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  const uint64_t l = in->l;
  assert(out->l >= l);
  uint64_t * tmp = (uint64_t *) safe_aligned_malloc(poly_size*sizeof(uint64_t));
  uint64_t * tmp2 = (uint64_t *) safe_aligned_malloc(out->N*sizeof(uint64_t));
  for (size_t i = 0; i < l; i++){
    memcpy(tmp2, out->coeffs[i], sizeof(uint64_t)*out->N);
    memset(out->coeffs[i], 0, sizeof(uint64_t)*out->N);
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->ntt->split_degree; j++){
      for (size_t k = 0; k < out->ntt->split_degree - j; k++){
        intel::hexl::EltwiseMultMod(tmp, &in->coeffs[i][j*poly_size], &tmp2[k*poly_size], poly_size, q, 1);
        intel::hexl::EltwiseAddMod(&out->coeffs[i][(j+k)* poly_size], &out->coeffs[i][(j+k)* poly_size], tmp, poly_size, q);
      }
      for (size_t k = out->ntt->split_degree - j; k < out->ntt->split_degree; k++){
        intel::hexl::EltwiseMultMod(tmp, &in->coeffs[i][j*poly_size], &tmp2[k*poly_size], poly_size, q, 1);
        intel::hexl::EltwiseMultMod(tmp, tmp, out->ntt->w[i], poly_size, q, 1);
        intel::hexl::EltwiseAddMod(&out->coeffs[i][(j+k - out->ntt->split_degree)* poly_size], &out->coeffs[i][(j+k - out->ntt->split_degree)* poly_size], tmp, poly_size, q);
      }
    }
  }
  free(tmp);
  free(tmp2);
}


/* out = in1*in2 mod (X^N + 1) */
void polynomial_mul_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  assert(out != in1);
  assert(out != in2);
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  uint64_t * tmp = (uint64_t *) safe_aligned_malloc(poly_size*sizeof(uint64_t));
  for (size_t i = 0; i < l; i++){
    memset(out->coeffs[i], 0, sizeof(uint64_t)*out->N);
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->ntt->split_degree; j++){
      for (size_t k = 0; k < out->ntt->split_degree - j; k++){
        intel::hexl::EltwiseMultMod(tmp, &in1->coeffs[i][j*poly_size], &in2->coeffs[i][k*poly_size], poly_size, q, 1);
        intel::hexl::EltwiseAddMod(&out->coeffs[i][(j+k)* poly_size], &out->coeffs[i][(j+k)* poly_size], tmp, poly_size, q);
      }
      for (size_t k = out->ntt->split_degree - j; k < out->ntt->split_degree; k++){
        intel::hexl::EltwiseMultMod(tmp, &in1->coeffs[i][j*poly_size], &in2->coeffs[i][k*poly_size], poly_size, q, 1);
        intel::hexl::EltwiseMultMod(tmp, tmp, out->ntt->w[i], poly_size, q, 1);
        intel::hexl::EltwiseAddMod(&out->coeffs[i][(j+k - out->ntt->split_degree)* poly_size], &out->coeffs[i][(j+k - out->ntt->split_degree)* poly_size], tmp, poly_size, q);
      }
    }
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
  free(tmp);
}

/* out = in1 - in2 mod (X^N + 1) */
void polynomial_sub_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    intel::hexl::EltwiseSubMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

/* out = in1 - in2 mod (X^N + 1) */
void polynomial_sub_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    intel::hexl::EltwiseSubMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}


void polynomial_add_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, RNSc_Polynomial in2){
  const uint64_t l = in1->l < in2->l ? in1->l : in2->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    intel::hexl::EltwiseAddMod(out->coeffs[i], in1->coeffs[i], in2->coeffs[i], out->N, q);
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

void polynomial_add_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, RNS_Polynomial in2){
  polynomial_add_RNSc_polynomial((RNSc_Polynomial) out, (RNSc_Polynomial) in1, (RNSc_Polynomial) in2);
}

void polynomial_RNSc_add_integer(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t in2){
  const uint64_t l = in1->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    memcpy(out->coeffs[i], in1->coeffs[i], out->N*sizeof(uint64_t));
    out->coeffs[i][0] = intel::hexl::AddUIntMod(in1->coeffs[i][0], q + in2, q);

  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

void polynomial_RNS_add_integer(RNS_Polynomial out, RNS_Polynomial in1, uint64_t in2){
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  const uint64_t l = in1->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    memcpy(out->coeffs[i], in1->coeffs[i], out->N*sizeof(uint64_t));
    for (size_t j = 0; j < poly_size; j++){
      out->coeffs[i][j] = intel::hexl::AddUIntMod(in1->coeffs[i][j], q + in2, q);
    }
  }
  for (size_t i = l; i < out->l; i++){
    memset(out->coeffs[i], 0, out->N*sizeof(uint64_t));
  } 
}

void polynomial_scale_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t scale){
  const uint64_t l = in1->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    intel::hexl::EltwiseFMAMod(out->coeffs[i], in1->coeffs[i], q + scale, nullptr, in1->N, q, 2);
  }
}

void polynomial_scale_addto_RNSc_polynomial(RNSc_Polynomial out, RNSc_Polynomial in1, uint64_t scale){
  const uint64_t l = in1->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    intel::hexl::EltwiseFMAMod(out->coeffs[i], in1->coeffs[i], scale, out->coeffs[i], in1->N, q, 1);
  }
}

void polynomial_scale_addto_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, uint64_t scale){
  polynomial_scale_addto_RNSc_polynomial((RNSc_Polynomial) out, (RNSc_Polynomial) in1, scale);
}


void polynomial_scale_RNS_polynomial(RNS_Polynomial out, RNS_Polynomial in1, uint64_t scale){
  polynomial_scale_RNSc_polynomial((RNSc_Polynomial) out, (RNSc_Polynomial) in1, scale);
}



void polynomial_scale_RNS_polynomial_RNS(RNS_Polynomial out, RNS_Polynomial in1, uint64_t * scale){
  const uint64_t l = in1->l;
  assert(out->l >= l);
  for (size_t i = 0; i < l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    intel::hexl::EltwiseFMAMod(out->coeffs[i], in1->coeffs[i], q + scale[i], nullptr, in1->N, q, 2);
  }
}

void polynomial_RNSc_negate(RNSc_Polynomial out, RNSc_Polynomial in){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    for (size_t j = 0; j < out->N; j++){
      out->coeffs[i][j] = intel::hexl::SubUIntMod(0, in->coeffs[i][j], q);
    }
  }
}

void polynomial_RNS_negate(RNS_Polynomial out, RNS_Polynomial in){
  polynomial_RNSc_negate((RNSc_Polynomial) out, (RNSc_Polynomial) in);
}

void polynomial_RNSc_to_RNS(RNS_Polynomial out, RNSc_Polynomial in){
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < out->l; i++){
    for (size_t k = 0; k < out->ntt->split_degree; k++){
      out->ntt->ntt[i]->ComputeForward(&out->coeffs[i][k*poly_size], &in->coeffs[i][k*poly_size], 1, 1);
    }
  }
}

void polynomial_RNS_to_RNSc(RNSc_Polynomial out, RNS_Polynomial in){
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < out->l; i++){
    for (size_t k = 0; k < out->ntt->split_degree; k++){
      out->ntt->ntt[i]->ComputeInverse(&out->coeffs[i][k*poly_size], &in->coeffs[i][k*poly_size], 1, 1);
    }
  }
}

// out = in + e(sigma)
void polynomial_RNSc_add_noise(RNSc_Polynomial out, RNSc_Polynomial in, double sigma){
  for (size_t j = 0; j < out->N; j++){
    const double real_noise = generate_normal_random(sigma);
    for (size_t i = 0; i < out->l; i++){
      const uint64_t p = out->ntt->ntt[i]->GetModulus();
      const uint64_t noise = ((uint64_t) (p + real_noise))%p;
      out->coeffs[i][j] = intel::hexl::AddUIntMod(in->coeffs[i][j], noise, p);
    }
  }
}

/* Assumes q_i/q_j < 2 for all i, j*/
void polynomial_base_extend_RNSc(RNSc_Polynomial out, uint64_t * in, uint64_t p){
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    if(q >= p){ // copy
      memcpy(out->coeffs[i], in, sizeof(uint64_t)*out->N);
    }else{ // reduce
      intel::hexl::EltwiseReduceMod(out->coeffs[i], in, out->N, q, 2, 1);
    }
  }
}

// reduce mod p_idx
void polynomial_RNSc_mod_reduce(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t idx){
  const uint64_t p = out->ntt->ntt[idx]->GetModulus();
  for (size_t i = 0; i < out->l; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    if(q >= p){ // copy
      memcpy(out->coeffs[i], in->coeffs[idx], sizeof(uint64_t)*out->N);
    }else{ // reduce
      intel::hexl::EltwiseReduceMod(out->coeffs[i], in->coeffs[idx], out->N, q, 2, 1);
    }
  }
}

// digit decompose polynomial supposing ||p|| < p_{l-1}
void polynomial_RNSc_decompose_small(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t log_base, uint64_t level){
  uint64_t tmp[out->N] __attribute__ ((aligned(64)));
  const uint64_t mask = (1ULL << log_base) - 1;
  const uint64_t shift = log_base*level;
  for (size_t i = 0; i < out->N; i++){
    tmp[i] = (in->coeffs[in->l - 1][i] >> shift) & mask;
  }
  // broadcast to Rq
  for (size_t i = 1; i < out->l; i++){
    memcpy(out->coeffs[i], tmp, sizeof(uint64_t)*out->N);
  }
}

void polynomial_RNSc_to_multiprecision(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t log_base, uint64_t level){
  uint64_t tmp[out->N] __attribute__ ((aligned(64)));
  const uint64_t mask = (1ULL << log_base) - 1;
  const uint64_t shift = log_base*level;
  for (size_t i = 0; i < out->N; i++){
    tmp[i] = (in->coeffs[in->l - 1][i] >> shift) & mask;
  }
  // broadcast to Rq
  for (size_t i = 1; i < out->l; i++){
    memcpy(out->coeffs[i], tmp, sizeof(uint64_t)*out->N);
  }
}

uint64_t get_D_j_mod_p_i(uint64_t * D, uint64_t j, uint64_t p_i, uint64_t w){
  uint64_t D_j = 1;
  for (size_t k = 0; k < w; k++){
    if(k != j) D_j = intel::hexl::MultiplyMod(D_j, D[k], p_i);
  }
  return D_j;
}

void polynomial_base_extend_RNSc_2(RNSc_Polynomial out, RNSc_Polynomial in){
  assert(out->l >= in->l);
  polynomial_copy_RNS_polynomial((RNS_Polynomial) out, (RNS_Polynomial) in);
  for (size_t i = in->l; i < out->l; i++){
    const uint64_t p_i = out->ntt->ntt[i]->GetModulus();
    memset(out->coeffs[i], 0, sizeof(uint64_t)*(out->N));
    for (size_t j = 0; j < in->l; j++){
      const uint64_t p_j = out->ntt->ntt[j]->GetModulus();
      const uint64_t D_j_mod_p_i = in->ntt->D[in->l][j][i];
      const uint64_t Dhat_j_mod_p_j = in->ntt->Dhat[in->l][j][j];
      for (size_t k = 0; k < out->N; k++){
        uint64_t tmp = intel::hexl::MultiplyMod(in->coeffs[j][k], Dhat_j_mod_p_j, p_j);
        tmp = intel::hexl::MultiplyMod(tmp, D_j_mod_p_i, p_i);
        out->coeffs[i][k] = intel::hexl::AddUIntMod(out->coeffs[i][k], tmp, p_i);
      }
    }
  }
}

void polynomial_RNS_get_hash(uint64_t * out, RNS_Polynomial p){
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  for (size_t i = 0; i < p->l; i++){
    blake3_hasher_update(&hasher, p->coeffs[i], p->N*sizeof(uint64_t));
  }
  blake3_hasher_finalize(&hasher, (uint8_t *) out, BLAKE3_OUT_LEN);
}

uint64_t * polynomial_RNS_get_hash_p(RNS_Polynomial p){
  uint64_t * out = (uint64_t *) safe_malloc (4*sizeof(uint64_t));
  polynomial_RNS_get_hash(out, p);
  return out;
}



// Removes the last prime of the RNS polynomial
// todo: optimize
void polynomial_base_reduce_RNSc_wo_free(RNSc_Polynomial out){
  const uint64_t l = out->l, p = out->ntt->ntt[l-1]->GetModulus();
  const uint64_t N = out->N;
  for (size_t i = 0; i < l - 1; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    const uint64_t inv_p = intel::hexl::InverseMod(p, q);
    if(q >= p){ // copy
      for (size_t j = 0; j < N; j++){
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }else{ // reduce
      for (size_t j = 0; j < N; j++){
        const uint64_t in_j = intel::hexl::ReduceMod<2>(out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], in_j, q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }
  } 
  out->l -= 1;
}

void polynomial_base_reduce_round_RNSc_wo_free(RNSc_Polynomial out){
  const uint64_t l = out->l, p = out->ntt->ntt[l-1]->GetModulus(), half_p = p/2;
  const uint64_t N = out->N;
  for (size_t i = 0; i < N; i++){
    out->coeffs[l - 1][i] = intel::hexl::AddUIntMod(out->coeffs[l - 1][i], half_p, p);
  }
  
  for (size_t i = 0; i < l - 1; i++){
    const uint64_t q = out->ntt->ntt[i]->GetModulus();
    const uint64_t inv_p = intel::hexl::InverseMod(p, q);
    if(q >= p){ // copy
      for (size_t j = 0; j < N; j++){
        out->coeffs[i][j] = intel::hexl::AddUIntMod(out->coeffs[i][j], half_p, q);
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }else{ // reduce
      const uint64_t half_p_mod_q = intel::hexl::ReduceMod<2>(half_p, q);
      for (size_t j = 0; j < N; j++){
        out->coeffs[i][j] = intel::hexl::AddUIntMod(out->coeffs[i][j], half_p_mod_q, q);
        const uint64_t in_j = intel::hexl::ReduceMod<2>(out->coeffs[l - 1][j], q);
        out->coeffs[i][j] = intel::hexl::SubUIntMod(out->coeffs[i][j], in_j, q);
        out->coeffs[i][j] = intel::hexl::MultiplyMod(out->coeffs[i][j], inv_p, q);
      }
    }
  } 
  out->l -= 1;
}


// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc(RNSc_Polynomial out){
  polynomial_base_reduce_RNSc_wo_free(out);
  free(out->coeffs[out->l]);
}

// Removes the last prime of the RNS polynomial
void polynomial_base_reduce_RNSc_and_scale(RNSc_Polynomial out, uint64_t p){
  polynomial_base_reduce_RNSc_wo_free(out);
  polynomial_scale_RNSc_polynomial(out, out, (out->ntt->ntt[out->l]->GetModulus())%p);
  free(out->coeffs[out->l]);
}



void polynomial_RNSc_permute(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t gen){
  const uint64_t N = out->N;
  assert(gen < N);
  assert(gen > 0);
  uint64_t idx = 0;
  polynomial_RNS_zero((RNS_Polynomial) out);
  for (size_t i = 0; i < N; i++){
    for (size_t j = 0; j < out->l; j++){
      out->coeffs[j][idx] = in->coeffs[j][i];
    }
    idx = intel::hexl::AddUIntMod(idx, gen, N);
  }
}

void polynomial_int_permute_mod_Q(IntPolynomial out, IntPolynomial in, uint64_t gen){
  const uint64_t N = in->N;
  uint64_t idx = 0;
  for (size_t i = 0; i < N; i++){
    out->coeffs[idx] = in->coeffs[i];
    idx = intel::hexl::AddUIntMod(idx, gen, N);
  }
}

// todo: check
void polynomial_RNSc_mul_by_xai(RNSc_Polynomial out, RNSc_Polynomial in, uint64_t a){
  assert(in != out);
  const uint64_t N = out->N;
  assert(a < N);
  polynomial_RNS_zero((RNS_Polynomial) out);
  for (size_t j = 0; j < out->l; j++){
    for (size_t i = 0; i < N-a; i++){
      out->coeffs[j][i+a] = in->coeffs[j][i];
    }
    for (size_t i = N-a; i < N; i++){
      out->coeffs[j][i-N+a] = -in->coeffs[j][i];
    }
  }
}

void polynomial_int_decompose_i(IntPolynomial out, IntPolynomial in, uint64_t Bg_bit, uint64_t l, uint64_t q, uint64_t bit_size, uint64_t i){
  const uint64_t N = in->N;
  const uint64_t h_mask = (1UL << Bg_bit) - 1;
  const uint64_t h_bit = bit_size - (i + 1) * Bg_bit;

  uint64_t offset = 1ULL << (bit_size - l * Bg_bit - 1);

  for (size_t c = 0; c < N; c++){
    const uint64_t coeff_off = in->coeffs[c] + offset;
    out->coeffs[c] = (coeff_off>>h_bit) & h_mask;
  }
}

IntPolynomial polynomial_new_int_polynomial(uint64_t N){
  IntPolynomial res;
  res = (IntPolynomial) safe_malloc(sizeof(*res));
  res->coeffs = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t) * N);
  res->N = N;
  return res;
}

void free_polynomial(void * p){
  free(((IntPolynomial) p)->coeffs);
  free(p);
}

// polynomial ring slot functions
// todo: optimize
void polynomial_RNS_broadcast_slot(RNS_Polynomial out, RNS_Polynomial in, uint64_t slot_idx){
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < in->l; i++){
    for (size_t j = 0; j < out->ntt->split_degree; j++){
      for (size_t k = 0; k < poly_size; k++){
        out->coeffs[i][j*poly_size + k] = in->coeffs[i][j*poly_size + slot_idx];
      }
    }
  }
}

void polynomial_RNS_broadcast_RNS_comp(RNS_Polynomial out, RNS_Polynomial in, uint64_t rns_comp){
  for (size_t i = 0; i < in->l; i++){
    for (size_t k = 0; k < out->N; k++){
      out->coeffs[i][k] = in->coeffs[rns_comp][k];
    }
  }
}

void polynomial_RNS_rotate_slot(RNS_Polynomial out, RNS_Polynomial in, uint64_t rot){
  assert(out != in);
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < in->l; i++){
    for (size_t j = 0; j < out->ntt->split_degree; j++){
      for (size_t k = 0; k < poly_size - rot; k++){
        out->coeffs[i][j*poly_size + k] = in->coeffs[i][j*poly_size + k + rot];
      }
      for (size_t k = poly_size - rot; k < poly_size; k++){
        out->coeffs[i][j*poly_size + k] = in->coeffs[i][j*poly_size + k + rot - poly_size];
      }
    }
  }
}

void polynomial_RNS_copy_slot(RNS_Polynomial out, uint64_t dst, RNS_Polynomial in, uint64_t src){
  const uint64_t poly_size = out->N/out->ntt->split_degree;
  for (size_t i = 0; i < in->l; i++){
    for (size_t j = 0; j < out->ntt->split_degree; j++){
      out->coeffs[i][j*poly_size + dst] = in->coeffs[i][j*poly_size + src];
    }
  }
}
