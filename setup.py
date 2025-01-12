import multiprocessing as mp
import os
import setuptools
import shutil
import subprocess
import sys

from pybind11.setup_helpers import Pybind11Extension, build_ext

CRC32C_DIR = os.path.join("third_party", "crc32c")
CRC32C_INCLUDE_DIR = os.path.join(CRC32C_DIR, "include")
CRC32C_BUILD_DIR = os.path.join(CRC32C_DIR, "build")

def build_crc32c():
  library_names = [
    "libcrc32c.a",
    "crc32c.lib"
  ]

  exists = [ 
    os.path.exists(os.path.join(CRC32C_BUILD_DIR, name)) 
    for name in library_names
  ]

  if any(exists):
    print("crackle: libcrc32c already built.")
    return

  if os.path.exists(CRC32C_BUILD_DIR):
    shutil.rmtree(CRC32C_BUILD_DIR)
  else:
    os.makedirs(CRC32C_BUILD_DIR)

  res = subprocess.check_call([
    "cmake", 
    f"-DCRC32C_BUILD_TESTS=0",
    f"-DCRC32C_BUILD_BENCHMARKS=0",
    f"-DBUILD_SHARED_LIBS=OFF",
    f"-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
    f"-DCMAKE_OSX_DEPLOYMENT_TARGET=11.0",
    f"-S {CRC32C_DIR}", 
    f"-B {CRC32C_BUILD_DIR}",
    "-DCMAKE_BUILD_TYPE=Release",
  ])

  if res != 0:
    print("crackle: unable to cmake libcrc32c.")
    return

  res = subprocess.check_call([
    "cmake", 
    f"--build",
    f"{CRC32C_BUILD_DIR}",
  ])

  if res != 0:
    print("crackle: unable to build libcrc32c.")
    return

build_crc32c()

extra_compile_args = []
if sys.platform == 'win32':
  extra_compile_args += [
    '/std:c++20', '/O2'
  ]
else:
  extra_compile_args += [
    '-std=c++2a', '-O3'
  ]

setuptools.setup(
  setup_requires=['pbr','pybind11','numpy'],
  cmdclass={"build_ext": build_ext},
  extras_require={
    "remote": [
      "cloud-files",
    ],
  },
  ext_modules=[
    Pybind11Extension(
        "fastcrackle",
        ["src/fastcrackle.cpp"],
        include_dirs=[CRC32C_INCLUDE_DIR],
        libraries=["crc32c"],
        library_dirs=[CRC32C_BUILD_DIR],
        extra_compile_args=extra_compile_args,
        language="c++",
    ),
  ],
  entry_points={
    "console_scripts": [
      "crackle=crackle_cli:main"
    ],
  },
  pbr=True
)
