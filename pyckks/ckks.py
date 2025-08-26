from __future__ import annotations
from .polynomial import *
from .ntt_ref import FFT, CC
from math import sqrt
from copy import copy

class Ciphertext:
  def __init__(self, ckks:CKKS, obj:tuple[Polynomial, Polynomial]=None, ell:int=None):
    self.ckks = ckks
    self.ring = ckks.ring
    self.obj = obj if obj is not None else (self.ring.alloc_polynomial(), self.ring.alloc_polynomial())
    self.ell = ell if ell is not None else self.ring.ell

  def __add__(self, other:Ciphertext|Polynomial):
    if(type(other) is Ciphertext):
      assert(self.ell == other.ell)
      return Ciphertext(self.ckks, self.ckks.add_ciphertexts(self.obj, other.obj), self.ell)
    else:
      return Ciphertext(self.ckks, (self.obj[0]+other, copy(self.obj[1])), self.ell)
  
  def __mul__(self, other:Ciphertext|int|Polynomial):
    if(type(other) is Ciphertext):
      assert(self.ell == other.ell)
      assert(self.ell >= 2), "modulus is too small to multiply" 
      return Ciphertext(self.ckks, self.ckks.mul(self.obj, other.obj, self.ell), self.ell - 1)
    else:
      return Ciphertext(self.ckks, (self.obj[0]*other, self.obj[1]*other), self.ell)
    
  def automomorphism(self, gen:int):
    return Ciphertext(self.ckks, self.ckks.automorphism(self.obj, gen, self.ell), self.ell - 1)

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
  
  def gen_ksk(self, key_in, key_out):
    ksk = []
    for i in range(self.ring.ell):
      scaling_factor = [0]*i + [1] + [0]*(self.ring.ell - 1 - i)
      ksk.append(self._encrypt(key_in*scaling_factor, key_out))
    return ksk

  def gen_ksk_auto(self, key:Polynomial, gen:int):
    key_auto = key.automorphism(gen)
    return self.gen_ksk(key_auto, key)
  
  def gen_ksk_auto_set(self, key:Polynomial, gens:list[int]):
    self.ksk_auto = {}
    for i in gens:
      self.ksk_auto[i] = self.gen_ksk_auto(key, i)

  def gen_key(self, gen_automorpism:bool=False):
    sk = self.ring.random_gaussian_element(self.sigma_s)
    sk_sq = sk*sk
    self.evk = []
    for i in range(self.ring.ell):
      scaling_factor = [0]*i + [1] + [0]*(self.ring.ell - 1 - i)
      minus_scaling_factor = [0]*i + [-1] + [0]*(self.ring.ell - 1 - i)
      self.minus_pw.append(minus_scaling_factor)
      self.evk.append(self._encrypt(sk_sq*scaling_factor, sk))
    if(gen_automorpism):
      self.ksk_auto = {}
      for i in range(1, 2*self.ring.N, 1):
        self.ksk_auto[i] = self.gen_ksk_auto(sk, i)
    return sk
  
  def automorphism(self, ct:tuple[Polynomial, Polynomial], gen, l):
    # TODO: perform automorphism without consuming levels
    # which requires base extension
    assert(gen in self.ksk_auto)
    ct0, ct1 = ct
    a = (ct1*self.delta).automorphism(gen)
    tildeC0 = (ct0*self.delta).automorphism(gen)
    tildeC1 = 0
    ksk_auto = self.ksk_auto[gen]
    for i in range(l):
      tildeC0 += (a % self.ring.primes[i])*ksk_auto[i][0]
      tildeC1 += (a % self.ring.primes[i])*ksk_auto[i][1]
    # rescale
    p = self.ring.primes[l - 1]
    inv_p = [int(ZZ.quo(ZZ*self.ring.primes[i])(p)**-1) for i in range(l - 1)] + [0]*(self.ring.ell - l + 1)
    tildeC0prime = (tildeC0 - (tildeC0 % p))*inv_p
    tildeC1prime = (tildeC1 - (tildeC1 % p))*inv_p
    return (tildeC0prime, tildeC1prime)

  def automorphism5j(self, ct, j, l):
    pass

  def encode(self, x:list[int|float]) -> Polynomial:
    assert(len(x) == self.N//2)
    # inv pi
    x = x + list(reversed(x))
    x_rev = self.fft.bit_reverse_reorder(x)
    x_rev = [i/self.ring.N for i in x_rev]
    fft_res = self.fft.GS_RN(x_rev, self.fft.br_inv_rous)
    float_m = [float(CC(i).real()) for i in fft_res]
    int_m_h1 = [floor(i*self.delta_sq) for i in float_m]
    int_m_h2 = [(round(i*(self.delta_sq**2)) - j*self.delta_sq) for i,j in zip(float_m,int_m_h1)]
    m_h1 = Polynomial(self.ring).from_array(int_m_h1)
    m_h2 = Polynomial(self.ring).from_array(int_m_h2)
    return m_h1*self.delta_poly + m_h2
  
  def decode(self, x, ell):
    delta = self.compute_delta(ell)
    Rq = self.ring.get_sage_ring()
    Q_l = math.prod(self.ring.primes[:ell])
    emb_m = Rq(x.get_polynomial())
    emb_m = [i%Q_l for i in emb_m]
    int_m = [int(i)/delta if int(i) < Q_l/2 else (int(i) - Q_l)/delta for i in emb_m]
    dec = self.fft.bit_reverse_reorder(self.fft.CT_NR(int_m,self.fft.br_rous))
    return [round(float(CC(i).real())) for i in dec][:self.N//2]
  
  def compute_delta(self, ell):
    delta = self.delta
    for i in range(1, self.ring.ell - ell + 1):
      delta = (delta**2)/self.ring.primes[self.ring.ell - i]
    return delta 

  def _encrypt(self, m:Polynomial, sk:Polynomial) -> tuple[Polynomial, Polynomial]:
    a = self.ring.random_element()
    e = self.ring.random_gaussian_element(self.sigma_e)
    return (-a*sk + m + e, a)
  
  def encrypt(self, m:list[int], sk:Polynomial):
    emb_m = self.encode(m)
    # print("emb_m", emb_m)
    return Ciphertext(self, self._encrypt(emb_m, sk))

  def _decrypt(self, ct:list[Polynomial], sk:Polynomial):
    b, a = ct
    return b + a*sk
  
  def decrypt(self, ct:Ciphertext, sk:Polynomial):
    phase = self._decrypt(ct.obj, sk)
    return self.decode(phase, ct.ell)

  def add_ciphertexts(self, ct0, ct1):
    ct00, ct01 = ct0
    ct10, ct11 = ct1
    return ct00 + ct10, ct01 + ct11
  
  def change_level(self, x:Ciphertext|Polynomial, ell):
    pass

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

def bench_ckks():
  REPs = 100
  Rq = Ring(2**14, 250, split_degree=4)
  print(Rq.primes, Rq.split_degree)
  ckks = CKKS(Rq, 2**49, 3.2, 3.2)
  sk = ckks.gen_key()
  a = list(range(1,Rq.N//2 + 1))
  a2 = list(range(1,Rq.N//2 + 1))
  b = ckks.encrypt(a, sk)
  b2 = ckks.encrypt(a2, sk)
  start = time.time_ns()
  for _ in range(REPs):
    bmul = b*b2
  end = time.time_ns()
  print("Multiplication: %lf ms" % ((end - start)/1000000/REPs))
  c = ckks.decrypt(bmul, sk)
  print(c[:10], sep="\n\n")

def test_ckks():
  Rq = Ring(2**14, 250, split_degree=4)
  ckks = CKKS(Rq, 2**49, 3.2, 3.2)
  sk = ckks.gen_key()
  a = [2, 2] + ([0] * (Rq.N//2 - 2))
  a2 = [2, 2] + ([0] * (Rq.N//2 - 2))
  b = ckks.encrypt(a, sk)
  b2 = ckks.encrypt(a2, sk)
  bmul = b*b2
  c = ckks.decrypt(bmul, sk)
  print(c[:10], sep="\n\n")
  for i in range(4):
    bmul = bmul*bmul
    c = ckks.decrypt(bmul, sk)
    print(c[:10], sep="\n\n")

def test_ckks_auto():
  Rq = Ring(16, 250, split_degree=4)
  ckks = CKKS(Rq, 2**49, 3.2, 3.2)
  sk = ckks.gen_key()
  ckks.gen_ksk_auto_set(sk, [1, 3, 5, 7, 9, 11, 13, 15])
  a = [i for i in range(Rq.N//2)]
  b = ckks.encrypt(a, sk)
  c = ckks.decrypt(b, sk)
  print(c[:10], sep="\n\n")
  for i in ckks.ksk_auto:
    print("\ni: %d" % i)
    b2 = b.automomorphism(i)
    c = ckks.decrypt(b2, sk)
    print(c[:10], sep="\n\n")

def test_ckks_add_mul():
  Rq = Ring(2**14, 250, split_degree=4)
  ckks = CKKS(Rq, 2**49, 3.2, 3.2)
  sk = ckks.gen_key()
  a = [2, 2] + ([0] * (Rq.N//2 - 2))
  a2 = [1, 2] + ([0] * (Rq.N//2 - 2))
  b = ckks.encrypt(a, sk)
  b2 = ckks.encrypt(a2, sk)
  b3 = [3, 5, 1] + ([0] * (Rq.N//2 - 3))
  b3_enc = ckks.encode(b3)
  b = b+b3_enc
  c = ckks.decrypt(b, sk)
  print(c[:10], sep="\n\n")
  bmul = b*b2
  c = ckks.decrypt(bmul, sk)
  print(c[:10], sep="\n\n")
  bmul = bmul*2
  c = ckks.decrypt(bmul, sk)
  print(c[:10], sep="\n\n")
  for i in range(4):
    bmul = bmul*bmul
    c = ckks.decrypt(bmul, sk)
    print(c[:10], sep="\n\n")

import random
def test_encode():
  Rq = Ring(2**14,650, split_degree=4)
  print("primes", Rq.primes)
  ckks = CKKS(Rq, 2**49, 3.2, 3.2)
  sk = ckks.gen_key()
  a = [2, 1] + [random.randint(1,2) for i in range(Rq.N//2 - 2)]
  a2 = [2, 1] + [random.randint(1,2) for i in range(Rq.N//2 - 2)]
  b = ckks.encode(a)
  b2 = ckks.encode(a2)
  bmul = b*b2
  p = ckks.ring.primes[Rq.ell - 1]
  bmul = (bmul - (bmul % p)) * ([int(ZZ.quo(ZZ*ckks.ring.primes[i])(p)**-1) for i in range(Rq.ell - 1)] + [0])
  c = ckks.decode(bmul, Rq.ell - 1)
  print(c[:10], sep="\n\n")
  for i in range(2, 6):
    bmul = bmul*bmul
    p = ckks.ring.primes[Rq.ell - i]
    bmul = (bmul - (bmul % p)) * ([int(ZZ.quo(ZZ*ckks.ring.primes[i])(p)**-1) for i in range(Rq.ell - i)] + ([0]*i))
    c = ckks.decode(bmul, Rq.ell - i)
    print("i: %d" % i, c[:10], sep="\n", end="\n\n")

if __name__ == "__main__":
  bench_ckks()
  # test_ckks_auto()
  # test_ckks_add_mul()
  # test_encode()