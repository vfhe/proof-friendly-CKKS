from __future__ import annotations
from polynomial import *
from librings import *
import random
from math import log2, ceil

def functionEvaluations(ell, A:list[Polynomial], r:list[Polynomial]) -> dict[tuple:Polynomial]:
  F = {}
  for i in range(1, ell + 1):
    for b in range(2**(ell - i)):
      for t in [0,1,2]:
        F[(i, b, t)] = A[b]*(1-t) + A[b + 2**(ell - i)]*t
      A[b] = A[b]*(1 - r[i-1]) + A[b + 2**(ell - i)]*r[i-1]
  return F

def sumCheck(ell, A:list[Polynomial], r:list[Polynomial]) -> list[list[Polynomial]]:
  F = functionEvaluations(A, r)
  a = [0]*ell
  for i in range(1, ell +1):
    a[i-1] = [0,0,0]
    for t in [0,1,2]:
      for b in range(2**(ell - i)):
        a[i-1][t] += F[(i,b,t)]
  return a
  
def fast_sumCheck(self, A:list[Polynomial], r:list[Polynomial]) -> list[list[Polynomial]]:
  ring = r[0].ring
  res = [[Polynomial(ring), Polynomial(ring), Polynomial(ring)] for _ in range(self.ell)]
  ty_res_in = c_void_p*3
  res_obj = [ty_res_in(*[i.obj for i in j]) for j in res]
  A_obj = [i.obj for i in A] 
  r_obj = [i.obj for i in r]
  ty_res = POINTER(c_void_p)*len(res)
  ty_A = c_void_p*len(A)
  ty_r = c_void_p*len(r)
  librings.lib.libra_sumcheck(ty_res(*res_obj), ty_A(*A_obj), ty_r(*r_obj), self.ell)
  return res

import time
import sys
def test_libra():
  REPs = 100
  Rq = Ring(2**13, 120, split_degree=2)
  num_var = int(sys.argv[1])

  A = [Rq.random_element() for _ in range(2**num_var)]
  # A_copy = [copy(i) for i in A]
  r = [Rq.random_element() for _ in range(num_var)]
  # r_copy = [copy(i) for i in r]
  
  start = time.time_ns()
  for _ in range(REPs):
    a0 = fast_sumCheck(A, r)
    # a1 = sumCheck(A_copy, r_copy)
  end = time.time_ns()
  
  a2 = 0
  for i in A:
    a2 += i
    
  print("Fast sumcheck", "pass" if a2 == a0 else "failed")
  print("Sumcheck: %lf ms" % ((end - start)/1000000/REPs))

  

if __name__ == "__main__":
  test_libra()
