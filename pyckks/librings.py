import ctypes, os
from ctypes import c_uint64, c_void_p, c_double, POINTER, c_bool
from math import floor, log2
import subprocess

class LibRings:
  def __init__(self) -> None:
    self.multithreaded = False
    self.num_threads = 1
    script_dir = os.path.dirname(__file__)
    lib_path = os.path.join(script_dir,"../lib/lib/librings.so")
    lib_make_path = os.path.join(script_dir,"../lib/")
    if not os.path.exists(lib_path):
      print("Compiling HEXL\n\n")
      os.system("make -C %s hexl" % lib_make_path)
      print("\n\n\n\nCompiling Rings C library")
      output = subprocess.check_output("cat /proc/cpuinfo", shell=True)
      if("avx512ifma" in str(output)):
        os.system("make -C %s ENABLE_VAES=1 ENABLE_AVX512=1" % lib_make_path)
      else:
        print("AVX512ifma not detected - Compiling unoptimized version" % lib_make_path)
        os.system("make -C %s")
    self.lib = ctypes.CDLL(lib_path)

librings = LibRings()