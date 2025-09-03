import os
import setuptools
import sys
import platform

from pybind11.setup_helpers import Pybind11Extension, build_ext

CRC32C_DIR = os.path.join("third_party", "fastcrc")    

extra_compile_args = []
if sys.platform == 'win32':
  extra_compile_args += [
    '/std:c++20', '/O2', '/arch:AVX2'
  ]
else:
  extra_compile_args += [
    '-std=c++2a', '-O3',
  ]
  # Apple Silicon arm64 machines cross compile for x86_64, 
  # but we want to not include these flags for aarch64
  if (
    platform.machine().upper() in ('X86_64', 'AMD64') 
  ):
    extra_compile_args += [
      '-msse4.2', '-mpclmul' # for x86 and cross compiling x86
    ]
  elif platform.machine().upper() == 'AARCH64':
    extra_compile_args += [
      "-march=armv8-a+crc+simd"  # Enable NEON for aarch64
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
        include_dirs=[CRC32C_DIR],
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
