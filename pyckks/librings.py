import ctypes, os
from ctypes import c_uint64, c_void_p, c_double, POINTER, c_bool
from math import floor, log2
import subprocess, sys, argparse

class LibRings:
  def __init__(self, full_recompile=False, recompile_librings=False, CC=None, CXX=None, arm=False) -> None:
    self.compiler = f"CC={CC} " if CC else ""
    self.compiler += f"CXX={CXX} " if CXX else ""
    self.compiler += "ARM=true" if arm else ""
    self.multithreaded = False
    self.num_threads = 1
    script_dir = os.path.dirname(__file__)
    lib_path = os.path.join(script_dir,"../lib/lib/librings.so")
    self.lib_make_path = os.path.join(script_dir,"../lib/")
    if full_recompile:
      print("Cleaning up...")
      os.system("make -C %s clean" % self.lib_make_path)
    if not os.path.exists(lib_path):
      self.compile_dependencies()
      self.compile()
    elif recompile_librings:
      self.compile(force_recompile=True)
    self.lib = ctypes.CDLL(lib_path)

  def is_AVX512ifma_supported(self):
    try:
      output = subprocess.check_output("cat /proc/cpuinfo", shell=True)
      return "avx512ifma" in str(output)
    except:
      return False
    
  def compile_dependencies(self):
    print("Compiling HEXL\n\n")
    os.system("make %s -C %s hexl" % (self.compiler, self.lib_make_path))

  def compile(self, force_recompile=False):
    _B = "-B" if force_recompile else ""
    print("\nCompiling Rings C library")
    if(self.is_AVX512ifma_supported()):
      os.system("make %s %s -C %s ENABLE_VAES=true ENABLE_AVX512=true" % (self.compiler, _B, self.lib_make_path))
      print("\nFinished compiling\nResuming python code\n")
    else:
      print("\033[93m" + "AVX512ifma not detected" +
      " - Compiling unoptimized version\n\033[0m" +
      "\033[1m\033[93mDO NOT USE IT FOR BENCHMARKING!\033[0m\033[0m\n" +
      "Press [enter] to continue")
      input()
      os.system("make %s %s -C %s" % (self.compiler, _B, self.lib_make_path))
      print("\nFinished compiling\nResuming python code\n")


parser = argparse.ArgumentParser()
parser.add_argument('--cc')
parser.add_argument('--cxx')
parser.add_argument('--full-recompile', action='store_true') 
parser.add_argument('--recompile', action='store_true')
parser.add_argument('--arm', action='store_true')
args = parser.parse_known_args()[0]

librings = LibRings(recompile_librings=args.recompile, 
                    full_recompile=args.full_recompile,
                    CC=args.cc,
                    CXX=args.cxx,
                    arm=args.arm)
