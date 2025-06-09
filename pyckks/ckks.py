from __future__ import annotations
from polynomial import *
from ntt_ref import FFT, CC
from math import sqrt


class CKKS:
  def __init__(self, ring:Ring, delta, sigma_s, sigma_e) -> None:
    self.ring = ring
    self.N = ring.N
    self.delta = delta
    self.delta_sq = int(sqrt(delta))
    self.delta_poly = Polynomial(ring).from_array([self.delta_sq] + [0]*(self.N - 1))
    self.sigma_s = sigma_s
    self.sigma_e = sigma_e
    self.fft = FFT(self.N)
    self.inv_p = self.compute_p_inv()
    self.minus_pw = []

  def compute_p_inv(self):
    inv_p = [0]*self.ring.ell
    for j in range(2, self.ring.ell):
      p = self.ring.primes[j - 1]
      inv_p[j] = [int(ZZ.quo(ZZ*self.ring.primes[i])(p)**-1) for i in range(j - 1)] + [0]*(self.ring.ell - j + 1)
    return inv_p

  def gen_key(self):
    sk = self.ring.random_gaussian_element(self.sigma_s)
    sk_sq = sk*sk
    self.evk = []
    for i in range(self.ring.ell):
      scaling_factor = [0]*i + [1] + [0]*(self.ring.ell - 1 - i)
      minus_scaling_factor = [0]*i + [-1] + [0]*(self.ring.ell - 1 - i)
      self.minus_pw.append(minus_scaling_factor)
      self.evk.append(self._encrypt(sk_sq*scaling_factor, sk))
    return sk
  
  def encode(self, x):
    assert(len(x) == self.N//2)
    # inv pi
    x = x + list(reversed(x))
    x_rev = self.fft.bit_reverse_reorder(x)
    fft_res = self.fft.GS_RN(x_rev, self.fft.br_inv_rous)
    float_m = [float(CC(i/self.N).real()) for i in fft_res]
    int_m_h1 = [round(i*self.delta_sq) for i in float_m]
    int_m_h2 = [(j*self.delta_sq - round(i*(self.delta_sq**2))) for i,j in zip(float_m,int_m_h1)]
    m_h1 = Polynomial(self.ring).from_array(int_m_h1)
    m_h2 = Polynomial(self.ring).from_array(int_m_h2)
    return m_h1*self.delta_poly + m_h2
  
  def decode(self, x, delta=None, l=None):
    if not delta: delta = self.delta
    Rq = self.ring.get_sage_ring()
    Q = self.ring.q_l
    Q_l = (Q//self.ring.primes[self.ring.ell - 1])
    emb_m = Rq(x.get_polynomial())
    emb_m = [i%Q_l for i in emb_m]
    int_m = [int(i) if int(i) < Q_l/2 else (int(i) - Q_l) for i in emb_m]
    dec = self.fft.bit_reverse_reorder(self.fft.CT_NR(int_m,self.fft.br_rous))
    return [round(float(CC(i).real())/delta) for i in dec][:self.N//2]
  
  def _encrypt(self, m:list[int], sk:Polynomial):
    a = self.ring.random_element()
    e = self.ring.random_gaussian_element(self.sigma_e)
    return (-a*sk + m + e, a)
  
  def encrypt(self, m:list[int], sk:Polynomial):
    emb_m = self.encode(m)
    # print("emb_m", emb_m)
    return self._encrypt(emb_m, sk)

  def _decrypt(self, ct:list[Polynomial], sk:Polynomial):
    b, a = ct
    return b + a*sk
  
  def decrypt(self, ct:list[Polynomial], sk:Polynomial, delta=None):
    phase = self._decrypt(ct, sk)
    return self.decode(phase, delta=delta)
  
  def scale(self, ct0, poly):
    pass

  def add(self, ct0, ct1):
    ct00, ct01 = ct0
    ct10, ct11 = ct1
    return ct00 + ct10, ct01 + ct11

  def mul(self, ct0, ct1, l):
    ct00, ct01 = ct0
    ct10, ct11 = ct1
    # pre multiply
    d0 = ct00*ct10
    d1 = ct00*ct11 + ct01*ct10
    d2 = ct01*ct11
    # relim
    tildeC0 = d0
    tildeC1 = d1
    for i in range(l):
      tildeC0 += (d2 % self.ring.primes[i])*self.evk[i][0]
      tildeC1 += (d2 % self.ring.primes[i])*self.evk[i][1]
    # rescale
    p = self.ring.primes[l - 1]
    inv_p = [int(ZZ.quo(ZZ*self.ring.primes[i])(p)**-1) for i in range(l - 1)] + [0]*(self.ring.ell - l + 1)
    tildeC0prime = (tildeC0 - (tildeC0 % p))*inv_p
    tildeC1prime = (tildeC1 - (tildeC1 % p))*inv_p
    return (tildeC0prime, tildeC1prime)

import time

def test_ckks():
  REPs = 1000
  Rq = Ring(2**14, 250, split_degree=4)
  print(Rq.primes, Rq.split_degree)
  ckks = CKKS(Rq, 2**49, 1, 1)
  sk = ckks.gen_key()
  a = list(range(1,Rq.N//2 + 1))
  a2 = list(range(1,Rq.N//2 + 1))
  b = ckks.encode(a)
  b2 = ckks.encode(a2)
  bmul = b*b2
  c = ckks.decode(bmul, delta=ckks.delta**2)
  b_dec = ckks.decode(b)
  # print(a, b_dec, c, sep="\n\n")
  print("########")
  b = ckks.encrypt(a, sk)
  b2 = ckks.encrypt(a2, sk)
  start = time.time_ns()
  for _ in range(REPs):
    bmul= ckks.mul(b, b2, len(Rq.primes))
  end = time.time_ns()
  print("Multiplication: %lf ms" % ((end - start)/1000000/REPs))
  # bmul= ckks.add(b, b2)
  # print("b", ckks.decrypt(b, sk, delta=2**20))
  # print("b2", ckks.decrypt(b2, sk, delta=2**20))
  c = ckks.decrypt(bmul, sk, delta=ckks.delta**2/Rq.primes[Rq.ell - 1])
  print(c[:10], sep="\n\n")


if __name__ == "__main__":
  test_ckks()