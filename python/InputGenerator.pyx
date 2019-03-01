from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp cimport bool
cimport libcpp.complex

import numpy as np

cdef extern from "InputGenerator.h":
    cdef cppclass CPP_InputGenerator "InputGenerator":
        CPP_InputGenerator(vector[vector[vector[double]]] &,
                           vector[vector[vector[libcpp.complex.complex[double]]]] &,
                           vector[vector[double]] &,
                           vector[double] &,
                           double, int, double, double, int, int, string)
        void generate_p()
        void generate_delta()
        void generate_all_for_determinant()
        void generate_all_for_one_spin()
        void generate_benchmark_data()



cdef class InputGenerator:
    cdef CPP_InputGenerator *thisptr
    def __cinit__(self,
                  bath_energies,
                  hybridizations,
                  local_energies,
                  Hubbard_U,
                  float dt, int n_times, float beta, float t_max, int real_grid_size, int imag_grid_size, str params_name):
        self.thisptr = new CPP_InputGenerator(bath_energies,
                                              hybridizations,
                                              local_energies,
                                              Hubbard_U,
                                              dt, n_times, beta, t_max, real_grid_size, imag_grid_size, params_name.encode('utf-8'))
    def __dealloc__(self):
        del self.thisptr

    def generate_p(self):
        self.thisptr.generate_p()
    def generate_delta(self):
        self.thisptr.generate_delta()
    def generate_p(self):
        self.thisptr.generate_p()
    def generate_all_for_determinant(self):
        self.thisptr.generate_all_for_determinant()
    def generate_all_for_one_spin(self):
        self.thisptr.generate_all_for_one_spin()
    def generate_benchmark_data(self):
        self.thisptr.generate_benchmark_data()