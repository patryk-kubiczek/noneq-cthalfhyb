from distutils.core import setup, Extension
from Cython.Build import cythonize
import os

project_name = "input_generator"

# IMPORT ENVIRONMENTAL VARIABLES
build_type = os.getenv("BUILD_TYPE", "Release")
cxx_flags = os.getenv("CXXFLAGS").split()
link_libraries = os.getenv("LINK_LIBRARIES", "ON")

# SET CYTHON AND C++ SOURCES
head_cpp_dir = os.sep.join(os.getcwd().split(os.sep)[:-1])
project_cpp_dir = os.path.join(head_cpp_dir, project_name)

def cpp_path(cpp_file):
    return os.path.join(project_cpp_dir, cpp_file)

source_files = ["InputGenerator.pyx"]
for file in os.listdir(project_cpp_dir):
    if file.endswith(".cpp") and file != "main.cpp" and "test" not in file:
        source_files.append(cpp_path(file))

# SET INCLUDE DIRECTORIES
include_dirs = [project_cpp_dir]

# SET LIBRARY DIRECTORIES
eom_library_dir = os.path.join(head_cpp_dir, "equations_of_motion", "lib")
library_dirs = [eom_library_dir]
runtime_library_dirs = [] 
if link_libraries == "ON":
    runtime_library_dirs = library_dirs

# SET LIBRARIES
libraries = ["eom"]
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
                "-Wno-unused-variable"]
compile_args.extend(["-DARMA_DONT_USE_WRAPPER"])
compile_args.extend(cxx_flags)

# SET DEBUGGING MODE?
undef_macros = []
if build_type == "Debug":
    undef_macros.append("NDEBUG")

setup(
    ext_modules =
      cythonize(
          Extension("InputGenerator",
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
