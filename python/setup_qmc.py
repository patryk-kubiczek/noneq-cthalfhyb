from distutils.core import setup, Extension
from Cython.Build import cythonize
import os


# IMPORT ENVIRONMENTAL VARIABLES
qmc_version = os.getenv("QMC", "CTHALFHYB")
build_type = os.getenv("BUILD_TYPE", "Release")
cxx_flags = os.getenv("CXXFLAGS").split()
link_libraries = os.getenv("LINK_LIBRARIES", "ON")

cpp_project = "cthalfhyb_qmc" if qmc_version == "CTHALFHYB" else "cthyb_qmc"

# SET CYTHON AND C++ SOURCES
head_cpp_dir = os.sep.join(os.getcwd().split(os.sep)[:-1])
project_cpp_dir = os.path.join(head_cpp_dir, cpp_project)

def cpp_path(cpp_file):
    return os.path.join(project_cpp_dir, cpp_file)

source_files = [qmc_version + ".pyx", cpp_path("QMC.cpp")]

# SET INCLUDE DIRECTORIES
include_dirs = [project_cpp_dir]

# SET LIBRARY DIRECTORIES
library_dirs = [cpp_path("lib")]
runtime_library_dirs = [] 
if link_libraries == "ON":
    runtime_library_dirs = library_dirs

# SET LIBRARIES
libraries = ["cthalfhyb" if qmc_version == "CTHALFHYB" else "cthyb"]
if os.getenv("USE_MKL", "OFF") == "ON":
    libraries.extend(["mkl_intel_lp64", "mkl_gnu_thread", "mkl_core", "gomp", "mkl_def"])
else:
    libraries.extend(["openblas"])
if link_libraries == "ON":
    libraries.extend(["m", "blas", "lapack"])

# SET COMPILE ARGUMENTS

compile_args = ["-std=c++14",
                "-Wno-sign-compare",
                "-Wno-misleading-indentation",
                "-Wno-unused-local-typedefs",
                "-Wno-unused-variable",
                "-Wno-reorder",
                "-Wno-return-type",
                "-Wno-maybe-uninitialized"]
compile_args.extend(["-DARMA_DONT_USE_WRAPPER", "-DUSE_MPI", "-DUSE_MPI4PY", "-D" + qmc_version + "_QMC"])
compile_args.extend(cxx_flags)

# SET DEBUGGING MODE?
undef_macros = []
if build_type == "Debug":
    undef_macros.append("NDEBUG")

setup(ext_modules =
    cythonize(
        Extension(qmc_version,
                  language = "c++",
                  sources = source_files,
                  include_dirs = include_dirs,
                  library_dirs = library_dirs,
                  runtime_library_dirs = runtime_library_dirs,
                  libraries = libraries,
                  extra_compile_args = compile_args,
                  extra_link_args = compile_args,
                  undef_macros = undef_macros
                  ),
        compiler_directives = {'language_level': 3,
                               'embedsignature': False,
                               'profile': False,
                               'binding': True}
    )
)
