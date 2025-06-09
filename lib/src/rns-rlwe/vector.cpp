#include "rns-rlwe.h"
// functions for vectors of RNS integers mod q

ZqVector alloc_ZqVector(uint64_t n, uint64_t * q, uint64_t l){
  ZqVector res;
  res = (ZqVector) safe_malloc(sizeof(*res));
  res->n = n;
  res->l = l;
  res->q = q;
  res->elements = (uint64_t **) safe_malloc(sizeof(uint64_t *)*l);
  for (size_t i = 0; i < l; i++){
    res->elements[i] = (uint64_t *) safe_aligned_malloc(sizeof(uint64_t)*n);
  }
  return res;
}

void ZqVector_add(ZqVector out, ZqVector in1, ZqVector in2){
  for (size_t i = 0; i < out->l; i++){
    intel::hexl::EltwiseAddMod(out->elements[i], in1->elements[i], in2->elements[i], out->n, out->q[i]);
  }
}

void ZqVector_sub(ZqVector out, ZqVector in1, ZqVector in2){
  for (size_t i = 0; i < out->l; i++){
    intel::hexl::EltwiseSubMod(out->elements[i], in1->elements[i], in2->elements[i], out->n, out->q[i]);
  }
}

void ZqVector_scale(ZqVector out, ZqVector in1, uint64_t scale){
  for (size_t i = 0; i < out->l; i++){
    intel::hexl::EltwiseFMAMod(out->elements[i], in1->elements[i], scale, nullptr, out->n, out->q[i], 1);
  }
}