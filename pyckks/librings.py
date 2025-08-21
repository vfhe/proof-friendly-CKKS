import ctypes, os
from ctypes import c_uint64, c_void_p, c_double, POINTER, c_bool
from math import floor, log2
import subprocess, sys

class LibRings:
  def __init__(self, full_recompile=False, recompile_librings=False) -> None:
    self.multithreaded = False
    self.num_threads = 1
    script_dir = os.path.dirname(__file__)
    lib_path = os.path.join(script_dir,"../lib/lib/librings.so")
    lib_make_path = os.path.join(script_dir,"../lib/")
    if full_recompile:
      print("Cleaning up...")
      os.system("make -C %s clean" % lib_make_path)
    if not os.path.exists(lib_path):
      print("Compiling HEXL\n\n")
      os.system("make -C %s hexl" % lib_make_path)
      print("\n\n\n\nCompiling Rings C library")
      output = subprocess.check_output("cat /proc/cpuinfo", shell=True)
      if("avx512ifma" in str(output)):
        os.system("make -C %s ENABLE_VAES=1 ENABLE_AVX512=1" % lib_make_path)
      else:
        print("\033[93m" + "AVX512ifma not detected" +
        " - Compiling unoptimized version\n\033[0m" +
        "\033[1m\033[93mDO NOT USE IT FOR BENCHMARKING!\033[0m\033[0m\n" +
        "Press [enter] to continue")
        input()
        os.system("make -C %s" % lib_make_path)
        print("\nRunning python code\n")
    elif recompile_librings:
      print("\nRecompiling Rings C library")
      output = subprocess.check_output("cat /proc/cpuinfo", shell=True)
      if("avx512ifma" in str(output)):
        os.system("make -B -C %s ENABLE_VAES=1 ENABLE_AVX512=1" % lib_make_path)
      else:
        print("\033[93m" + "AVX512ifma not detected" +
        " - Compiling unoptimized version\n\033[0m" +
        "\033[1m\033[93mDO NOT USE IT FOR BENCHMARKING!\033[0m\033[0m\n" +
        "Press [enter] to continue")
        input()
        os.system("make -B -C %s" % lib_make_path)
        print("\nRunning python code\n")
    self.lib = ctypes.CDLL(lib_path)

if("--full-recompile" in sys.argv):
  librings = LibRings(full_recompile=True)
elif("--recompile" in sys.argv):
  librings = LibRings(recompile_librings=True)
else:
  librings = LibRings()
