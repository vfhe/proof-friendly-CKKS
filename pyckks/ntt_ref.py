from __future__ import annotations
import sage.all
from sage.all import PolynomialRing, ZZ, QQ, RR, CC, GF, cyclotomic_polynomial, random_prime, primitive_root, pi, exp
import random, sys
from math import sqrt, log

class FFT:
  def __init__(self, N:int) -> FFT:
    self.rou = CC(exp(2*pi*1j/(2*N)))
    self.rous = [self.rou**i for i in range(N)]
    # just change to 5*i and see what happens
    self.inv_rous = [i**-1 for i in self.rous]
    self.br_rous = self.bit_reverse_reorder(self.rous)
    self.br_inv_rous = self.bit_reverse_reorder(self.inv_rous)
    self.N = N

  def DFT(self, p, ws):
    res = []
    n = len(p)
    for i in range(n):
      res.append(sum([p[j]*ws[(i*j) % n] for j in range(n)])/n)
    return res
  

# mul = lambda x,y: [i*j for i,j in zip(x,y)]


  # Cooley-Tukey Natural to Reverse
  # ws in bit reversed order
  def CT_NR(self, p, ws):
    A = p.copy()
    n = len(A)
    t = n
    m = 1
    while m < n:
      t >>= 1
      for i in range(m):
        j1 = 2*i*t
        j2 = j1+t
        w = ws[m+i]
        for j in range(j1, j2):
          U = A[j]
          V = A[j + t]*w
          A[j] = U + V
          A[j+t] = U - V
      m <<= 1
    return A

  # Gentleman-Sande GS Reverse to Natural 
  # ws in bit reversed order
  def GS_RN(self, p, ws):
    A = p.copy()
    n = len(A)
    t = 1
    m = n
    while m > 1:
      j1 = 0
      h = m//2
      for i in range(h):
        j2 = j1 + t
        w = ws[h + i]
        for j in range(j1, j2):
          U = A[j]
          V = A[j + t]
          A[j] = U + V
          A[j + t] = (U - V)*w
        j1 += 2*t
      t <<= 1
      m >>= 1
    return A

  def bit_reverse_reorder(self, v, size=None):
    if(not size): size = int(log(len(v), 2))
    brev = lambda x: int(bin(x)[2:].rjust(size, "0")[::-1], 2) 
    return [v[brev(i)] for i in range(2**size)]

## example:
def basic_test():
  q = 769
  N = 16
  fft = FFT(N)

  mul = lambda x,y: [i*j for i,j in zip(x,y)]
  
  a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
  b = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

  # print(fft.CT_NR(a, fft.br_rous))
  # print(fft.CT_NR(b, fft.br_rous))
  print(mul(fft.GS_RN(fft.CT_NR(a, fft.br_rous), fft.br_inv_rous), [1/N]*N))

  print(mul(fft.GS_RN(mul(fft.CT_NR(a, fft.br_rous), fft.CT_NR(b, fft.br_rous)), fft.br_inv_rous), [1/N]*N))

if __name__ == "__main__":
  basic_test()