from __future__ import annotations
from librings import *
import math
from sage.all import next_prime, previous_prime, is_prime, ZZ, cyclotomic_polynomial, crt, var
import numpy as np
import rlwe

next_power_of_2 = lambda x: 1 << int(math.ceil(math.log2(x)))

class Ring:
  def __init__(self, N, mod_size, split_degree=None, NTT_proc=None, primes=None, prime_size=49, exceptional_set_size=128) -> None:
    self.prime_size = prime_size
    self.ell = math.ceil(mod_size/self.prime_size)
    if(not split_degree):
      split_degree = next_power_of_2(exceptional_set_size/self.prime_size)
    self.split_degree = split_degree
    self.N = N

    self.lib = librings.lib

    # functions
    self.lib.new_incomplete_ntt_list.argtypes = (ctypes.POINTER(c_uint64), c_uint64, c_uint64, c_uint64)
    self.lib.new_incomplete_ntt_list.restype = c_void_p

    # 
    self.lib.incNTT_get_rou_matrix.argtypes = c_void_p,
    self.lib.incNTT_get_rou_matrix.restype = POINTER(self.ell*POINTER((self.N//self.split_degree)*c_uint64))

    # new polynomial
    self.lib.polynomial_new_RNS_polynomial.argtypes  = (c_uint64, c_uint64, c_void_p)
    self.lib.polynomial_new_RNS_polynomial.restype = c_void_p

    # copy
    self.lib.polynomial_copy_RNS_polynomial.argtypes = (c_void_p, c_void_p)

    # free polynomial
    self.lib.free_RNS_polynomial.argtypes = (c_void_p,)

    # sampling
    self.lib.polynomial_gen_random_RNSc_polynomial.argtypes = (c_void_p,)
    self.lib.polynomial_gen_gaussian_RNSc_polynomial.argtypes = (c_void_p,c_double)

    # arithmetic
    self.lib.polynomial_mul_RNS_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.polynomial_multo_RNS_polynomial.argtypes = (c_void_p, c_void_p)
    self.lib.polynomial_sub_RNS_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.polynomial_sub_RNSc_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.polynomial_add_RNSc_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.polynomial_add_RNS_polynomial.argtypes = (c_void_p, c_void_p, c_void_p)
    self.lib.polynomial_scale_RNSc_polynomial.argtypes = (c_void_p, c_void_p, c_uint64)
    self.lib.polynomial_scale_RNS_polynomial_RNS.argtypes = (c_void_p, c_void_p, POINTER(c_uint64))    
    self.lib.polynomial_RNSc_negate.argtypes = (c_void_p, c_void_p)
    self.lib.polynomial_RNSc_add_integer.argtypes  = (c_void_p, c_void_p, c_uint64)
    self.lib.polynomial_RNS_add_integer.argtypes  = (c_void_p, c_void_p, c_uint64)

    # ntt
    self.lib.polynomial_RNSc_to_RNS.argtypes = (c_void_p, c_void_p)
    self.lib.polynomial_RNS_to_RNSc.argtypes = (c_void_p, c_void_p)

    # base extend
    self.lib.polynomial_base_extend_RNSc_2.argtypes = (c_void_p, c_void_p)
    self.lib.polynomial_RNSc_mod_reduce.argtypes = (c_void_p, c_void_p, c_uint64)

    # logic
    self.lib.polynomial_eq.argtypes = (c_void_p, c_void_p)
    self.lib.polynomial_eq.restype = c_bool

    # slot fuctions
    self.lib.polynomial_RNS_broadcast_slot.argtypes = (c_void_p, c_void_p, c_uint64)
    self.lib.polynomial_RNS_rotate_slot.argtypes =  (c_void_p, c_void_p, c_uint64)
    self.lib.polynomial_RNS_copy_slot.argtypes =  (c_void_p, c_uint64, c_void_p, c_uint64)

    # input
    self.lib.int_array_to_RNS.argtypes = (c_void_p, POINTER(c_uint64))
    self.lib.array_to_RNS.argtypes = (c_void_p, POINTER(POINTER(c_uint64)))

    # hash
    self.lib.polynomial_RNS_get_hash.argtypes = POINTER(c_uint64), c_void_p
    self.lib.polynomial_RNS_get_hash_p.argtypes = c_void_p,
    self.lib.polynomial_RNS_get_hash_p.restype = POINTER(c_uint64*4)

    if(primes):
      assert(len(primes) >= self.ell), "not enough primes for subring for size 2^%d" % mod_size
      self.primes = primes[:self.ell]
    else:
      self.primes = self.gen_primes()

    self.q_l = math.prod(self.primes)

    self.bit_size = math.ceil(math.log2(self.q_l))

    if(NTT_proc):
      self.NTT = NTT_proc
    else:
      array_type = c_uint64*self.ell
      self.NTT = self.lib.new_incomplete_ntt_list(array_type(*self.primes), self.split_degree, self.N, self.ell)

    w = self.lib.incNTT_get_rou_matrix(self.NTT)
    self.rou_matrix = [list(i.contents) for i in w.contents]

    self.Rij = [0]*self.ell
    x = var("x")
    for i in range(self.ell):
      Zq = ZZ.quo(ZZ*self.primes[i])
      self.Rij[i] = []
      for j in range(self.N//self.split_degree):
        R_ij = Zq['x'].quo(x**self.split_degree - self.rou_matrix[i][j])
        self.Rij[i].append(R_ij)

  def sub_ring(self, mod_size):
    res = Ring(self.N, mod_size, self.split_degree, self.NTT, self.primes)
    return res

  def is_sub_ring(self, extension:Ring):
    return (extension.primes[:self.ell] == self.primes)

  #generates the RNS primes of size prime_size
  def gen_primes(self):
    k = 2*self.N//self.split_degree
    num_primes = 0
    primes = []
    a = ((2**self.prime_size - 1)//k) | 1
    while num_primes < self.ell:
      if(is_prime(a*k + 1)):
        # check (a*prim_rou_order + 1)%(2*prim_rou_order) !=1 is not needed because a is odd
        primes.append(a*k + 1)
        num_primes += 1
      a += 2
    return primes
  
  #in this file used to test against the definition of Ring, used in other files when efficiency is not required
  def get_sage_ring(self, prime=None):
    # q_l = math.prod(self.primes)
    mod = self.primes[prime] if prime!=None else self.q_l
    Zq = ZZ.quo(ZZ*mod)
    Rq = Zq["x"].quo(cyclotomic_polynomial(self.N*2))
    return Rq
  
  def alloc_polynomial(self):
    return self.lib.polynomial_new_RNS_polynomial(self.N, self.ell, self.NTT)
  
  def random_element(self, ntt=True):
    return Polynomial(self).sample_uniform(ntt)
  
  def random_gaussian_element(self, sigma, ntt=True):
    return Polynomial(self).sample_gaussian(sigma, ntt)

  
  def random_exceptional(self, size="minimal", ntt=True):
    return Polynomial(self).sample_exceptional(size, ntt)

from enum import Enum

repr = Enum('Polynomial Representation', ['empty', 'ntt', 'coeff'])

class Polynomial:
  def __init__(self, ring:Ring) -> None:
    self.ring = ring
    self.obj = ring.alloc_polynomial()
    self.repr = repr.empty

  def from_array(self, array:list):
    array_type = c_uint64*len(array)
    self.ring.lib.int_array_to_RNS(self.obj, array_type(*array))
    self.repr = repr.ntt
    return self
  
  def from_bigint_array(self, array:list):
    matrix = (POINTER(c_uint64)*self.ring.ell)()
    for i in range(self.ring.ell):
      p = self.ring.primes[i]
      array_type = c_uint64*self.ring.N
      matrix[i] = array_type(*[v%p for v in array])
    self.ring.lib.array_to_RNS(self.obj, matrix)
    self.repr = repr.ntt
    return self

  # inverse using sage
  def slow_inverse(self):
    matrix = self.get_coeff_matrix(repr=repr.ntt)
    d = self.ring.split_degree
    matrix_inv = [0]*self.ring.ell
    for i in range(self.ring.ell):
      matrix_inv[i] = []
      for j in range(self.ring.N//d):
        poly_sage = self.ring.Rij[i][j](matrix[i][j*d:(j+1)*d])
        matrix_inv[i] += list(poly_sage**-1)
    return Polynomial(self.ring).from_coeff_matrix(matrix_inv, repr.ntt)

  def __del__(self) -> None:
    self.ring.lib.free_RNS_polynomial(self.obj)

  # def from_int(self, i:int) -> Polynomial:
  def multiply(out, in1, in2):
    assert(in1.repr == in2.repr == repr.ntt)
    out.ring.lib.polynomial_mul_RNS_polynomial(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr

  def negate(out, in1=None):
    if(not in1): in1 = out
    out.ring.lib.polynomial_RNSc_negate(out.obj, in1.obj)
    out.repr = in1.repr

  def sub(out, in1, in2):
    assert(in1.repr == in2.repr)
    if(in1.repr == repr.ntt):
      out.ring.lib.polynomial_sub_RNS_polynomial(out.obj, in1.obj, in2.obj)
    else:
      out.ring.lib.polynomial_sub_RNSc_polynomial(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr

  def add(out, in1, in2):
    assert(in1.repr == in2.repr)
    if(in1.repr == repr.ntt):
      out.ring.lib.polynomial_add_RNS_polynomial(out.obj, in1.obj, in2.obj)
    else:
      out.ring.lib.polynomial_add_RNSc_polynomial(out.obj, in1.obj, in2.obj)
    out.repr = in1.repr

  def to_NTT(self):
    if(self.repr == repr.ntt): return
    self.ring.lib.polynomial_RNSc_to_RNS(self.obj, self.obj)
    self.repr = repr.ntt
  
  def to_coeff(self):
    if(self.repr == repr.coeff): return
    self.ring.lib.polynomial_RNS_to_RNSc(self.obj, self.obj)
    self.repr = repr.coeff
  
  def to_repr(self, repr:repr):
    if(repr == repr.ntt): self.to_NTT()
    if(repr == repr.coeff): self.to_coeff()

  def sample_uniform(self, ntt=True):
    self.ring.lib.polynomial_gen_random_RNSc_polynomial(self.obj)
    self.repr = repr.ntt if ntt else repr.coeff
    return self
  
  def sample_gaussian(self, sigma, ntt=True):
    self.ring.lib.polynomial_gen_gaussian_RNSc_polynomial(self.obj, sigma)
    self.repr = repr.coeff
    if(ntt): self.to_NTT()
    return self
  
  def sample_exceptional(self, size="minimal", ntt=True):
    self.ring.lib.polynomial_gen_random_RNSc_polynomial(self.obj)
    self.ring.lib.polynomial_RNS_broadcast_slot(self.obj, self.obj, 0)
    self.repr = repr.ntt
    if(not ntt): self.to_coeff()
    return self
  
  def base_extend(self, ring:Ring|None = None, out:Polynomial|None = None):
    self.to_coeff()
    if(out is None and ring is not None): 
      out_ = Polynomial(ring)
    elif type(out) is Polynomial:
      out_:Polynomial = out
    self.ring.lib.polynomial_base_extend_RNSc_2(out_.obj, self.obj)
    out_.repr = repr.coeff
    return out_
  
  def __eq__(self, value: Polynomial|int) -> bool:
    if(type(value) is int):
      self.to_coeff()
      return all([all([value == i[0]] + [j == 0 for j in i[1:]]) for i in self])
    elif(type(value) is list):
      self.to_coeff()
      s_list = list(self)
      return all([all([value[i] == s_list[i][0]] + [j == 0 for j in s_list[i][1:]]) for i in range(self.ring.ell)])
    else:
      return self.ring.lib.polynomial_eq(self.obj, value.obj)
  
  def __repr__(self) -> str:
    return str(list(self))
  
  def get_hash_pointer(self):
    self.to_NTT()
    return self.ring.lib.polynomial_RNS_get_hash_p(self.obj)
  
  def get_hash(self):
    hash_p = self.get_hash_pointer()
    return list(hash_p.contents)

  def __hash__(self):
    assert(False), "not a Python hashable object. Call polynomial.get_hash() for cryptographic hash"

    element_hash = (c_uint64*4)()
    self.ring.lib.polynomial_RNS_get_hash(element_hash, self.obj)
    return hash((int(i) for i in element_hash))
  
  # for compatibility with list, not an efficient iterator
  def __iter__(self):
    self.to_coeff()
    out = self.get_coeff_matrix(repr=repr.coeff)
    return iter(out)
  
  def get_coeff_matrix(self, repr=repr.coeff):
    if(self.repr != repr): self.to_repr(repr)
    pts = ctypes.cast(self.obj, POINTER(POINTER(POINTER(c_uint64*self.ring.N)*self.ring.ell)))
    values = [list(i.contents) for i in list(pts.contents.contents)]
    modMask = self.ring.split_degree - 1
    poly_size = self.ring.N//self.ring.split_degree
    c = values
    out = []
    for i in range(self.ring.ell):
      out_i = [0]*self.ring.N
      for j in range(self.ring.N):
        out_i[j] = c[i][(j&modMask)*poly_size + j//self.ring.split_degree]
      out += [out_i]
    return out

  def from_coeff_matrix(self, matrix, repr=repr.coeff):
    pts = ctypes.cast(self.obj, POINTER(POINTER(POINTER(c_uint64*self.ring.N)*self.ring.ell)))
    modMask = self.ring.split_degree - 1
    poly_size = self.ring.N//self.ring.split_degree
    for i in range(self.ring.ell):
      for j in range(self.ring.N):
        pts.contents.contents[i].contents[(j&modMask)*poly_size + j//self.ring.split_degree] = int(matrix[i][j])
    self.repr = repr
    return self
  
  #returns the list of coefficients of a Polynomial element
  #list of list with RNS coefficients??
  def get_polynomial(self) -> list:
    self.to_coeff()
    rns = list(self)
    if(self.ring.ell == 1):
      return list(sum(rns, start=[]))
    return [crt([rns[i][j] for i in range(self.ring.ell)],self.ring.primes) for j in range(self.ring.N)]
  
  def copy(self) -> Polynomial:
    res = Polynomial(self.ring)
    self.ring.lib.polynomial_copy_RNS_polynomial(res.obj, self.obj)
    res.repr = self.repr
    return res

  def decompose(self, base, small = False):
    res = []
    lst = self.get_polynomial()
    for j in range(self.ring.bit_size//base):
      dec = []
      for i in range(len(lst)):
        dec.append( lst[i] & ( ( 1 << base ) - 1))
        lst[i] >>= base
      res.append(Polynomial(self.ring).from_array(dec))
    return res

  def __mod__(self, value: int):
    self.to_coeff()
    res = Polynomial(self.ring)
    value_idx = self.ring.primes.index(value)
    self.ring.lib.polynomial_RNSc_mod_reduce(res.obj, self.obj, value_idx)
    return res

  def __copy__(self) -> Polynomial:
    return self.copy()
  
  def __mul__(self, other) -> Polynomial:
    if(type(other) is Polynomial):
      self.to_NTT() 
      other.to_NTT() 
      res = Polynomial(self.ring)
      res.multiply(self, other)
    elif(type(other) is int):
      if(other == 0): return 0
      if(other == 1): return self.copy()
      res = Polynomial(self.ring)
      self.ring.lib.polynomial_scale_RNSc_polynomial(res.obj, self.obj, other)
      res.repr = self.repr
    elif(type(other) is list):
      res = Polynomial(self.ring)
      assert(len(other) == self.ring.ell)
      array_type = c_uint64*self.ring.ell
      self.ring.lib.polynomial_scale_RNS_polynomial_RNS(res.obj, self.obj, array_type(*other))
      res.repr = self.repr
    else:
      print(type(other))
      assert(False), "not implemented"
    return res
  
  def __rmul__(self, other):
    return self.__mul__(other)
  
  def __radd__(self, other):
    return self.__add__(other)
  
  def __rsub__(self, other):
    return (-self) + other
  
  def __imul__(self, other):
    if(type(other) is Polynomial):
      self.to_NTT() 
      other.to_NTT()
      self.ring.lib.polynomial_multo_RNS_polynomial(self.obj, other.obj)
    elif(type(other) is int):
      if(other == 0): return 0
      if(other == 1): return self
      assert(other < 2**self.ring.prime_size) # large scaling not implemented
      self.ring.lib.polynomial_scale_RNSc_polynomial(self.obj, self.obj, other)
    else:
      print(type(other))
      assert(False) # not implemented
    return self
  
  def __add__(self, other:Polynomial|int):
    if(type(other) is int):
      if(other == 0): 
        return self.copy()
      else: 
        res = Polynomial(self.ring)
        if(self.repr == repr.coeff):
          self.ring.lib.polynomial_RNSc_add_integer(res.obj, self.obj, other)
        else:
          self.ring.lib.polynomial_RNS_add_integer(res.obj, self.obj, other)
        res.repr = self.repr
        return res
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    res = Polynomial(self.ring)
    res.add(self, other)
    return res
  
  def __iadd__(self, other:Polynomial|rlwe.RLWE|int):
    if(type(other) is int):
      if(other == 0): return self
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    if(type(other) is Polynomial): 
      self.add(self, other)
      return self
    else: # assume it's rlwe
      print(other)
      return other + self
  
  def __isub__(self, other:Polynomial|int):
    if(type(other) is int):
      if(other == 0): return self
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    self.sub(self, other)
    return self
  
  def __neg__(self):
    res = Polynomial(self.ring)
    res.negate(self)
    return res

  def __sub__(self, other:Polynomial|int):
    if(type(other) is int):
      if(other == 0): return self.copy()
      else: assert(False) # not implemented
    if(self.repr != other.repr):
      self.to_NTT()
      other.to_NTT()
    res = Polynomial(self.ring)
    res.sub(self, other)
    return res


def check_primes(primes, rou_order):
  rou_order = int(rou_order)
  for p in primes:
    Zp = ZZ.quo(ZZ*p)
    try:
      Zp.zeta(rou_order)
    except:
      print(Zp, "has no RoU of order", rou_order)
      return False
    
    try:
      Zp.zeta(rou_order*2)
      print(Zp, "has RoU of order", 2*rou_order, ":", Zp.zeta(rou_order*2))
      return False
    except:
      pass
  return True


def test_decomposition(ring):
  f = ring.random_element(ntt=False)
  Rq = ring.get_sage_ring()
  fs = Rq(f.get_polynomial())
  base = 4
  lst = f.decompose(base)
  #lsts = list(map(Rq, lst))
  lsts = list(map(lambda x: Rq(x.get_polynomial()), lst))
  g = 0
  for i in range(ring.bit_size//base):
    g += lsts[i]*(2**(i*base))
  print("Decomposition:", "pass" if g == fs else "fail")

def test_inverse(ring):
  Rq = ring.get_sage_ring()
  p0 = ring.random_element(ntt=False)
  p0_inv = p0.slow_inverse()
  p0_inv_p0 = p0_inv * p0
  print("slow_inverse:", "pass" if Rq(p0_inv_p0.get_polynomial()) == 1 else "fail")

def test_poly():
  r = Ring(2**14, 300, split_degree=4)

  print("Primes:", r.primes)

  r2 = r.sub_ring(100)
  Rq = r.get_sage_ring()

  p0 = r.random_element(ntt=False)
  p1 = r.random_element(ntt=False)

  p0s = Rq(p0.get_polynomial())
  p1s = Rq(p1.get_polynomial())

  p1r = Polynomial(r).from_bigint_array(list(p1s))

  print("Input/Output:", "pass" if p1r.get_polynomial() == list(p1s) else "fail")

  p2 = p0 * p1
  p2s = p0s*p1s
  print("Multiplication:", "pass" if p2.get_polynomial() == list(p2s) else "fail")

  p2 = p0 + p1
  p2s = p0s + p1s
  print("Add:", "pass" if p2.get_polynomial() == list(p2s) else "fail")


  p2 = 1 - p0
  p2s = 1 - p0s
  print("rsub:", "pass" if p2.get_polynomial() == list(p2s) else "fail")

  test_inverse(r)

  p3_r2 = r2.random_element()
  p3_r = p3_r2.base_extend(r)


  p3_diff = list(Rq(p3_r.get_polynomial()) - Rq(p3_r2.get_polynomial()))
  r2_mod = int(r2.get_sage_ring().base_ring().order())
  p3_diff = [int(i)/r2_mod if int(i)%r2_mod == 0 else 100 for i in p3_diff]
  print("Base extend:", "pass" if max(p3_diff) < 3  else "fail")
  test_decomposition(r)
  return
  print("Ring splitting:", "pass" if check_primes(r.primes, 2*r.N/r.split_degree) else "fail")

  

if __name__ == "__main__":
  test_poly()
