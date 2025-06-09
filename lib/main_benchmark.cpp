#include <rns-rlwe.h>
#include <benchmark_util.h>


#define EXT_MOD 5
#define BASE_MOD 562949954961409

void field_mul_addto(uint64_t * out, uint64_t * in1, uint64_t * in2, uint64_t size){
  for (size_t i = 0; i < size; i++){
    for (size_t j = 0; j < size - i; j++){
      const uint64_t tmp = intel::hexl::MultiplyMod(in1[i], in2[j], BASE_MOD);
      out[i+j] = intel::hexl::AddUIntMod(out[i+j], tmp, BASE_MOD);
    }
    for (size_t j = size - i; j < size; j++){
      const uint64_t tmp = intel::hexl::MultiplyMod(in1[i], in2[j], BASE_MOD);
      const uint64_t tmp2 = intel::hexl::MultiplyMod(tmp2, EXT_MOD, BASE_MOD);
      out[i+j] = intel::hexl::AddUIntMod(out[i+j], tmp, BASE_MOD);
    }    
  }
}


#define REPS 100

#define MATRIX_y 256*5
#define MATRIX_x ((int)(MATRIX_y*0.178/5))
#define RS_cw_sz 64
#define MATRIX_hw 7
#define FIELD_sz 4
// out is transposed
void matrix_mul(uint64_t * out, uint64_t * in, uint64_t ** A_val,  uint64_t ** A_idx){
  for (size_t i = 0; i < MATRIX_x; i++){
    for (size_t j = 0; j < MATRIX_hw; j++){
      const uint64_t A_idx_ij = A_idx[i][j];
      for (size_t k = 0; k < MATRIX_y; k++){
        field_mul_addto(&out[(i*MATRIX_y + k)*FIELD_sz], &in[(A_idx_ij*MATRIX_y + k)*FIELD_sz], &A_val[i][4*j], FIELD_sz);
      }
    }    
  }
}

void RS_encode(uint64_t * data, intel::hexl::NTT * ntt){
  uint64_t tmp[RS_cw_sz];
  for (size_t i = 0; i < MATRIX_y*FIELD_sz; i++){
    memset(tmp, 0, sizeof(uint64_t)*RS_cw_sz);
    memcpy(tmp, &data[i*MATRIX_x], sizeof(uint64_t)*MATRIX_x);
    ntt->ComputeForward(tmp, tmp, 1, 1);
  }
}


void test_encoding(){
  const uint64_t n = MATRIX_y*MATRIX_y, N = RS_cw_sz;
  uint64_t primes[1] = {BASE_MOD};
  auto ntt = new_ntt_list(primes, N, 1);
  uint64_t * A_idx[MATRIX_x];
  uint64_t * A_val[MATRIX_x];
  // generate matrix A
  for (size_t i = 0; i < MATRIX_x; i++){
    A_idx[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*MATRIX_hw);
    A_val[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*MATRIX_hw*FIELD_sz);
    generate_random_bytes(sizeof(uint64_t)*MATRIX_hw, (uint8_t *) A_idx[i]);
    generate_random_bytes(sizeof(uint64_t)*MATRIX_hw*FIELD_sz, (uint8_t *) A_val[i]);
    for (size_t j = 0; j < MATRIX_hw; j++){
      A_idx[i][j] &= (MATRIX_y - 1);
    }
    for (size_t j = 0; j < MATRIX_hw*FIELD_sz; j++){
      A_val[i][j] &= ((1ULL << 50) - 1);
    }
    intel::hexl::EltwiseReduceMod(A_val[i], A_val[i], MATRIX_hw*FIELD_sz, primes[0], 2, 1);
  }
  // generate input data
  uint64_t * in = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*n*FIELD_sz);
  generate_random_bytes(sizeof(uint64_t)*n*FIELD_sz, (uint8_t *) in);
  for (size_t i = 0; i < n*FIELD_sz; i++) in[i] &= ((1ULL << 50) - 1);
  intel::hexl::EltwiseReduceMod(in, in, n*FIELD_sz, primes[0], 2, 1);
  
  // compute encoding
  uint64_t * out = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*MATRIX_x*MATRIX_y*FIELD_sz);

  MEASURE_TIME("", 10, "Matrix mul",
    matrix_mul(out, in, A_val, A_idx);
  );
  MEASURE_TIME("", REPS, "RS encode",
    RS_encode(out, ntt[0]);
  );

  free(out);
}

// #define REPS 100

void test_arith(){
  const uint64_t L = 3, N = 1ULL<<13, split_degree = 4;
  uint64_t primes[7] = {562949954961409, 562949955551233, 562949957386241, 562949957779457, 562949960335361, 562949960728577, 562949961842689};
  incNTT ntt = new_incomplete_ntt_list(primes, split_degree, N, L);
  RNS_Polynomial p0 = polynomial_new_RNS_polynomial(N, L, ntt);
  RNS_Polynomial p1 = polynomial_new_RNS_polynomial(N, L, ntt);
  RNS_Polynomial p2 = polynomial_new_RNS_polynomial(N, L, ntt);
  polynomial_gen_random_RNSc_polynomial((RNSc_Polynomial) p0);
  polynomial_gen_random_RNSc_polynomial((RNSc_Polynomial) p1);
  MEASURE_TIME("", REPS, "NTT",
    polynomial_RNSc_to_RNS(p0, (RNSc_Polynomial) p0);
  );

  MEASURE_TIME("", REPS, "MUL",
    polynomial_mul_RNS_polynomial(p2, p0, p1)
  );

  MEASURE_TIME("", REPS, "INTT",
    polynomial_RNS_to_RNSc((RNSc_Polynomial) p2, p2);
  );

}
#define THREADS 1

typedef struct _encoding_arg{
  uint64_t * in;
  intel::hexl::NTT * ntt;
} encoding_arg;


void * encoding_thread(void * arg_p){
  encoding_arg * arg = (encoding_arg*) arg_p;
  RS_encode(arg->in, arg->ntt);
  return NULL;
}

void test_encoding_mp(){
  const uint64_t N = RS_cw_sz;
  uint64_t primes[1] = {BASE_MOD};
  pthread_t threads[THREADS];
  encoding_arg args[THREADS];
  auto ntt = new_ntt_list(primes, N, 1);

  for (size_t i = 0; i < THREADS; i++){
    // generate input data
    uint64_t * in = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*MATRIX_x*MATRIX_y*FIELD_sz);
    generate_random_bytes(sizeof(uint64_t)*MATRIX_x*MATRIX_y*FIELD_sz, (uint8_t *) in);
    for (size_t i = 0; i < MATRIX_x*MATRIX_y*FIELD_sz; i++) in[i] &= ((1ULL << 50) - 1);
    intel::hexl::EltwiseReduceMod(in, in, MATRIX_x*MATRIX_y*FIELD_sz, primes[0], 2, 1);
    args[i].in = in;
    args[i].ntt = ntt[0];
  }
  
  MEASURE_TIME("", REPS, "RS encoding",
    if(THREADS == 1){
      encoding_thread((void *) &args[0]);
    }else{
      for (size_t i = 0; i < THREADS; i++){
        pthread_create(&threads[i], NULL, *encoding_thread, (void *) &(args[i]));
      }
      for (size_t i = 0; i < THREADS; i++){
        pthread_join(threads[i], NULL);
      }
    }
  );
}

int main(int argc, char const *argv[])
{
  test_encoding_mp();
  // test_arith();
  // test_encoding();
  return 0;
}
