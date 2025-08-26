from __future__ import annotations
from .librings import *
import math
from sage.all import next_prime, previous_prime, is_prime, ZZ, cyclotomic_polynomial, crt
import numpy as np
from .polynomial import *


# linear only RLWE for the proofs
# CKKS will be implemented using polynomials, because of GKR

class LibRLWE:
  def __init__(self) -> None:

    self.lib = librings.lib

    self.lib.rlwe_new_RNS_gaussian_key.argtypes = (c_uint64, c_uint64, c_double, c_void_p, c_double)
    self.lib.rlwe_new_RNS_gaussian_key.restype = c_void_p
    self.lib.rlwe_alloc_RNS_sample.argtypes = (c_uint64, c_uint64, c_void_p)
    self.lib.rlwe_alloc_RNS_sample.restype = c_void_p

    self.lib.rlwe_RNSc_encrypt_RNS.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_RNSc_decrypt_RNS.argtypes = (c_void_p, c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_RNSc_encrypt_int.argtypes = (c_void_p, c_void_p, c_uint64, c_uint64)

    self.lib.rlwe_RNSc_to_RNS.argtypes = (c_void_p, c_void_p)
    self.lib.rlwe_RNS_to_RNSc.argtypes = (c_void_p, c_void_p) 
    self.lib.free_rlwe_RNS_sample.argtypes = (c_void_p,)
    self.lib.rlwe_copy_RNS_sample.argtypes = (c_void_p, c_void_p) 

    self.lib.rlwe_scale_RNS_rlwe_RNS.argtypes = (c_void_p, POINTER(c_uint64))
    self.lib.rlwe_RNS_mul_by_poly.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_add_RNSc_sample.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_sub_RNSc_sample.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_RNS_phase.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_add_RNSc_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.rlwe_sub_RNSc_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)


lib_rlwe = LibRLWE()
    

class RLWE_Key:
  def __init__(self, ring:Ring, plaintext_ring:Ring, sigma_s:float, sigma_e:float) -> None:
    self.obj = lib_rlwe.lib.rlwe_new_RNS_gaussian_key(ring.N, ring.ell, sigma_s, ring.NTT, sigma_e)
    self.ring = ring
    self.plaintext_ring = plaintext_ring
    self.tmp = Polynomial(self.ring)
    

  def encrypt(self, msg:Polynomial|int, out:RLWE=None):
    if(not out): out = RLWE(self.ring)
    if(type(msg) == int):
      lib_rlwe.lib.rlwe_RNSc_encrypt_int(out.obj, self.obj, msg, self.plaintext_ring.ell)
    else:
      msg.to_coeff()
      lib_rlwe.lib.rlwe_RNSc_encrypt_RNS(out.obj, self.obj, msg.obj)
    out.repr = repr.coeff
    return out
  
  def phase(self, rlwe:RLWE, out:Polynomial=None):
    if(not out):
      out = Polynomial(self.ring)
    rlwe.to_NTT()
    lib_rlwe.lib.rlwe_RNS_phase(out.obj, rlwe.obj, self.obj)
    out.repr = repr.ntt
    return out

  def decrypt(self, rlwe:RLWE, out=None):
    if(type(rlwe) is Polynomial): return rlwe
    if(not out): out = Polynomial(self.plaintext_ring)
    rlwe.to_NTT()
    lib_rlwe.lib.rlwe_RNSc_decrypt_RNS(out.obj, self.tmp.obj, rlwe.obj, self.obj)
    out.repr = repr.coeff
    return out


class RLWE:
  def __init__(self, ring:Ring) -> None:
    self.obj = lib_rlwe.lib.rlwe_alloc_RNS_sample(ring.N, ring.ell, ring.NTT)
    self.ring = ring
    self.repr = repr.empty

  def __del__(self) -> None:
    lib_rlwe.lib.free_rlwe_RNS_sample(self.obj)

  def encrypt(self, msg:Polynomial, key:RLWE_Key):
    key.encrypt(msg, out=self)

  def to_NTT(self):
    if(self.repr == repr.ntt): return
    lib_rlwe.lib.rlwe_RNSc_to_RNS(self.obj, self.obj)
    self.repr = repr.ntt

  def to_coeff(self):
    if(self.repr == repr.coeff): return
    lib_rlwe.lib.rlwe_RNS_to_RNSc(self.obj, self.obj)
    self.repr = repr.coeff

  def multiply_poly(out, in_rlwe, in_poly):
    assert(in_rlwe.repr == in_poly.repr == repr.ntt)
    lib_rlwe.lib.rlwe_RNS_mul_by_poly(out.obj, in_rlwe.obj, in_poly.obj)
    out.repr = repr.ntt

  def multiply_scalar(out, in_rlwe, pointer_to_int_list):
    lib_rlwe.lib.rlwe_scale_RNS_rlwe_RNS(out.obj, in_rlwe.obj, pointer_to_int_list)
    out.repr = in_rlwe.ntt

  def add_RLWE(out, in1, in2):
    assert(in1.repr == in2.repr)
    lib_rlwe.lib.rlwe_add_RNSc_sample(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr

  def add_poly(out:RLWE, in1:RLWE, in2:Polynomial):
    assert(in1.ring == in2.ring), "trying to add things in different rings"
    assert(in1.repr == in2.repr)
    lib_rlwe.lib.rlwe_add_RNSc_polynomial(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr

  def sub_RLWE(out, in1, in2):
    assert(in1.repr == in2.repr)
    lib_rlwe.lib.rlwe_sub_RNSc_sample(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr
  
  def sub_poly(out, in1, in2):
    assert(in1.repr == in2.repr)
    lib_rlwe.lib.rlwe_sub_RNSc_polynomial(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr

  def copy(self):
    res = RLWE(self.ring)
    lib_rlwe.lib.rlwe_copy_RNS_sample(res.obj, self.obj)
    res.repr = self.repr
    return res
  
  def __copy__(self):
    return self.copy()

  def __add__(self, other:RLWE|int|Polynomial):
    if(type(other) == int):
      if(other == 0): return self.copy()
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    print(type(other))
    res = RLWE(self.ring)
    if(type(other) == RLWE): res.add_RLWE(self, other)
    if(type(other) == Polynomial): res.add_poly(self, other)
    return res
  
  def __iadd__(self, other:RLWE|int|Polynomial):
    if(type(other) == int):
      if(other == 0): return self
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    if(type(other) == RLWE): self.add_RLWE(self, other)
    if(type(other) == Polynomial): self.add_poly(self, other)
    return self

  def __sub__(self, other) -> RLWE:
    if(type(other) == int):
      if(other == 0): return self.copy()
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    res = RLWE(self.ring)
    if(type(other) == RLWE): res.sub_RLWE(self, other)
    if(type(other) == Polynomial): res.sub_poly(self, other)
    return res
  
  def __isub__(self, other:RLWE|int|Polynomial):
    if(type(other) == int):
      if(other == 0): return self
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    if(type(other) == RLWE): self.sub_RLWE(self, other)
    if(type(other) == Polynomial): self.sub_poly(self, other)
    return self

  def __mul__(self, other) -> RLWE:
    res = RLWE(self.ring)
    if(type(other) is Polynomial):
      if(other.ring != self.ring):
        assert(other.ring.is_sub_ring(self.ring))
        other_sr = other.base_extend(self.ring)
      else:
        other_sr = other
      self.to_NTT()
      other_sr.to_NTT()
      res.multiply_poly(self, other_sr)
    else: # assuming other is a pointer to int
      res.multiply_scalar(self, other)
    return res
  
  def __rmul__(self, other):
    return self.__mul__(other)
  
  def __radd__(self, other):
    return self.__add__(other)
  
  # def __rsub__(self, other):
  #   self.negate()
  #   return self.__add__(other)

def test_rlwe():
  Rq = Ring(2048, 300)
  Rp = Rq.sub_ring(30)
  key = RLWE_Key(Rq, Rp, 3.2, 3.2)
  Rp_sage = Rp.get_sage_ring()

  m0 = Rp.random_element()
  m1 = Rp.random_element()
  m0s = Rp_sage(m0.get_polynomial())
  m1s = Rp_sage(m1.get_polynomial())
  c0 = key.encrypt(m0)
  c1 = key.encrypt(m1)

  m0dec = key.decrypt(c0)
  m1dec = key.decrypt(c1)
  print("Encrypt 1:", "pass" if m0dec.get_polynomial() == list(m0s) else "fail")
  print("Encrypt 2:", "pass" if m1dec.get_polynomial() == list(m1s) else "fail")

  c2 = c0 + c1
  m2s = m0s + m1s
  m2dec = key.decrypt(c2)
  print("Add:", "pass" if m2dec.get_polynomial() == list(m2s) else "fail")

  c2 = c0 - c1
  m2s = m0s - m1s
  m2dec = key.decrypt(c2)
  print("Sub:", "pass" if m2dec.get_polynomial() == list(m2s) else "fail")

  z = Rp.random_element()
  zs = Rp_sage(z.get_polynomial())
  c2 = c0*z
  m2s = m0s*zs
  m2dec = key.decrypt(c2)
  print("Mul:", "pass" if m2dec.get_polynomial() == list(m2s) else "fail")

  z = Rp.random_element()
  z2 = Rp.random_element()
  z3 = Rp.random_element()
  zs = Rp_sage(z.get_polynomial())
  z2s = Rp_sage(z2.get_polynomial())
  z3s = Rp_sage(z3.get_polynomial())
  c2 = z3
  c2 += c0*z
  c2 += c0*z2
  m2s = z3s + m0s*zs + m0s*z2s
  m2dec = key.decrypt(c2)
  print("Sum mul:", "pass" if m2dec.get_polynomial() == list(m2s) else "fail")



if __name__ == "__main__":
  test_rlwe()