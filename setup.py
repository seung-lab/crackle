import setuptools
from pybind11.setup_helpers import Pybind11Extension, build_ext


ext_modules = [
    Pybind11Extension(
        "fastcrackle",
        ["src/fastcrackle.cpp"],
        extra_compile_args=["-std=c++17"],
    ),
]

setuptools.setup(
  setup_requires=['pbr'],
  cmdclass={"build_ext": build_ext},
  ext_modules=ext_modules,
  pbr=True
)
