#include <rns-rlwe.h>

class Libra
{
private:
  uint64_t ell;
public:
  Libra(uint64_t ell);
  ~Libra();
  void functionEvaluations(RNS_Polynomial *** F, RNS_Polynomial * A, RNS_Polynomial * r);
  void sumCheck(RNS_Polynomial ** a, RNS_Polynomial * A, RNS_Polynomial * r);
  void sumCheck_mt(RNS_Polynomial ** a, RNS_Polynomial * A, RNS_Polynomial * r);
};

Libra::Libra(uint64_t ell){
  this->ell = ell;
}

Libra::~Libra()
{
}

void Libra::functionEvaluations(RNS_Polynomial *** F, RNS_Polynomial * A, RNS_Polynomial * r){
  RNS_Polynomial OneMinusR = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);
  for (size_t i = 1; i < this->ell + 1; i++){
    const uint64_t half = 1ULL << (this->ell - i);
    polynomial_RNS_negate(OneMinusR, r[i-1]);
    polynomial_RNS_add_integer(OneMinusR, OneMinusR, 1);
    for (size_t b = 0; b < half; b++){
      ///// Compute F
      // t = 0
      polynomial_copy_RNS_polynomial(F[i][b][0], A[b]);
      // t = 1
      polynomial_copy_RNS_polynomial(F[i][b][1], A[b + half]);
      // t = 2
      polynomial_scale_RNS_polynomial(F[i][b][2], A[b + half], 2);
      polynomial_sub_RNS_polynomial(F[i][b][2], F[i][b][2], A[b]);
      ///// Update A
      polynomial_mul_RNS_polynomial(A[b + half], A[b + half], r[i-1]);
      polynomial_mul_RNS_polynomial(A[b], A[b], OneMinusR);
      polynomial_add_RNS_polynomial(A[b], A[b], A[b + half]);
    }
  }
  free_RNS_polynomial(OneMinusR);
}

#define THREADS 48

typedef struct _sumcheck_arg{
  uint64_t start, end, half, i;
  RNS_Polynomial OneMinusR;
  RNS_Polynomial ** a, * A, * r;
} sumcheck_arg;

void * sumcheck_thread(void * arg_p){
  sumcheck_arg * arg = (sumcheck_arg*) arg_p;
  RNS_Polynomial tmp = polynomial_new_RNS_polynomial(arg->A[0]->N, arg->A[0]->l, arg->A[0]->ntt);
  RNS_Polynomial tmp2 = polynomial_new_RNS_polynomial(arg->A[0]->N, arg->A[0]->l, arg->A[0]->ntt);
  const uint64_t i = arg->i;
  // printf("start %lu, end: %lu half: %lu\n", arg->start, arg->end, arg->half);
  for (size_t b = arg->start; b < arg->end; b++){
    ///// Compute F and sum to a
    // t = 0
    polynomial_add_RNS_polynomial(arg->a[i-1][0], arg->a[i-1][0], arg->A[b]);
    // t = 1
    polynomial_add_RNS_polynomial(arg->a[i-1][1], arg->a[i-1][1], arg->A[b + arg->half]);
    // t = 2
    polynomial_scale_addto_RNS_polynomial(arg->a[i-1][2], arg->A[b + arg->half], 2);
    polynomial_sub_RNS_polynomial(arg->a[i-1][2], arg->a[i-1][2], arg->A[b]);
    ///// Update A
    polynomial_mul_RNS_polynomial(tmp, arg->A[b + arg->half], arg->r[i-1]);
    polynomial_mul_RNS_polynomial(tmp2, arg->A[b], arg->OneMinusR);
    polynomial_add_RNS_polynomial(arg->A[b], tmp2, tmp);
  }
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
  return NULL;
}

void Libra::sumCheck_mt(RNS_Polynomial ** a, RNS_Polynomial * A, RNS_Polynomial * r){
  RNS_Polynomial OneMinusR = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);
  RNS_Polynomial tmp = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);
  RNS_Polynomial tmp2 = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);

  RNS_Polynomial * a_threads[THREADS][this->ell];
  for (size_t i = 0; i < ell; i++){
    polynomial_RNS_zero(a[i][0]);
    polynomial_RNS_zero(a[i][1]);
    polynomial_RNS_zero(a[i][2]);
    for (size_t j = 0; j < THREADS; j++){
      a_threads[j][i] = polynomial_new_array_of_RNS_polynomials(A[0]->N, A[0]->l, 3, A[0]->ntt);
      polynomial_RNS_zero(a_threads[j][i][0]);
      polynomial_RNS_zero(a_threads[j][i][1]);
      polynomial_RNS_zero(a_threads[j][i][2]);
    } 
  }
  pthread_t threads[THREADS];
  sumcheck_arg args[THREADS];

  for (size_t i = 1; i < this->ell + 1; i++){
    const uint64_t half = 1ULL << (this->ell - i);
    polynomial_RNS_negate(OneMinusR, r[i-1]);
    polynomial_RNS_add_integer(OneMinusR, OneMinusR, 1);
    uint64_t num_threads = half/32 > THREADS ? THREADS : half / 32;
    uint64_t started_threads = 0;
    if(num_threads <= 1){
      // printf("halft st: %lu\n", half);
      for (size_t b = 0; b < half; b++){
        ///// Compute F and sum to a
        // t = 0
        polynomial_add_RNS_polynomial(a[i-1][0], a[i-1][0], A[b]);
        // t = 1
        polynomial_add_RNS_polynomial(a[i-1][1], a[i-1][1], A[b + half]);
        // t = 2
        polynomial_scale_addto_RNS_polynomial(a[i-1][2], A[b + half], 2);
        polynomial_sub_RNS_polynomial(a[i-1][2], a[i-1][2], A[b]);
        ///// Update A
        polynomial_mul_RNS_polynomial(tmp, A[b + half], r[i-1]);
        polynomial_mul_RNS_polynomial(tmp2, A[b], OneMinusR);
        polynomial_add_RNS_polynomial(A[b], tmp2, tmp);
      }
    }else{
      const uint64_t n_points = ceil((float)(half)/num_threads);
      for (size_t j = 0; j < num_threads; j++){
        args[j].start = j*n_points; 
        args[j].end = (j+1)*n_points;
        if(args[j].end > half) args[j].end = half;
        args[j].OneMinusR = OneMinusR;
        args[j].i = i;
        args[j].half = half;
        args[j].a = a_threads[j];
        args[j].A = A;
        args[j].r = r;
        pthread_create(&threads[j], NULL, *sumcheck_thread, (void *) &(args[j])); started_threads++;
        if(args[j].end >= half) break;
      }
      for (size_t j = 0; j < started_threads; j++){
        pthread_join(threads[j], NULL);
        polynomial_add_RNS_polynomial(a[i-1][0], a[i-1][0], a_threads[j][i-1][0]);
        polynomial_add_RNS_polynomial(a[i-1][1], a[i-1][1], a_threads[j][i-1][1]);
        polynomial_add_RNS_polynomial(a[i-1][2], a[i-1][2], a_threads[j][i-1][2]);
      }
    }
  }
  free_RNS_polynomial(OneMinusR);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
  for (size_t j = 0; j < THREADS; j++){
    for (size_t i = 0; i < ell; i++){
      free_RNS_polynomial_array(3, a_threads[j][i]);
    }
  }
}


void Libra::sumCheck(RNS_Polynomial ** a, RNS_Polynomial * A, RNS_Polynomial * r){
  RNS_Polynomial OneMinusR = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);
  RNS_Polynomial tmp = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);
  RNS_Polynomial tmp2 = polynomial_new_RNS_polynomial(A[0]->N, A[0]->l, A[0]->ntt);
  for (size_t i = 0; i < ell; i++){
    polynomial_RNS_zero(a[i][0]);
    polynomial_RNS_zero(a[i][1]);
    polynomial_RNS_zero(a[i][2]);
  }

  for (size_t i = 1; i < this->ell + 1; i++){
    const uint64_t half = 1ULL << (this->ell - i);
    polynomial_RNS_negate(OneMinusR, r[i-1]);
    polynomial_RNS_add_integer(OneMinusR, OneMinusR, 1);
    for (size_t b = 0; b < half; b++){
      ///// Compute F and sum to a
      // t = 0
      polynomial_add_RNS_polynomial(a[i-1][0], a[i-1][0], A[b]);
      // t = 1
      polynomial_add_RNS_polynomial(a[i-1][1], a[i-1][1], A[b + half]);
      // t = 2
      polynomial_scale_addto_RNS_polynomial(a[i-1][2], A[b + half], 2);
      polynomial_sub_RNS_polynomial(a[i-1][2], a[i-1][2], A[b]);
      ///// Update A
      polynomial_mul_RNS_polynomial(tmp, A[b + half], r[i-1]);
      polynomial_mul_RNS_polynomial(tmp2, A[b], OneMinusR);
      polynomial_add_RNS_polynomial(A[b], tmp2, tmp);
    }
  }
  free_RNS_polynomial(OneMinusR);
  free_RNS_polynomial(tmp);
  free_RNS_polynomial(tmp2);
}

#ifdef __cplusplus
extern "C" { 
#endif

// C wrapper
void libra_sumcheck(RNS_Polynomial ** a, RNS_Polynomial * A, RNS_Polynomial * r, uint64_t ell){
  Libra libra(ell);
  libra.sumCheck_mt(a, A, r);
}

#ifdef __cplusplus
}
#endif