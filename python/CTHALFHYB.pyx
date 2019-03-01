from libcpp.vector cimport vector
from libcpp.string cimport string
cimport libcpp.complex
from libcpp cimport bool

cdef extern from "QMC.h":
    cdef cppclass CXX_QMC "QMC":
        CXX_QMC(string, double, double,
                double, int, double, int)
        CXX_QMC(string, double, double,
                double, double, double, double, double, double,
                libcpp.complex.complex[double], libcpp.complex.complex[double], double, double,
                libcpp.complex.complex[double], libcpp.complex.complex[double], double, double)
        CXX_QMC(string, double, double,
                double, int, double, int,
                double, double, double, double, double, double,
                vector[libcpp.complex.complex[double]], vector[libcpp.complex.complex[double]], vector[double], vector[double],
                bool)

        void set_n_imag_times(int)
        void set_n_real_times(int)
        void set_n_matsubara(int)
        void set_block_size(int)
        void set_n_blocks(int)
        void run_measurement(long)
        void run_warm_up(long)
        void run_measurement(double)
        void run_warm_up(double)
        void collect_results()
        void reset_measurements()
        void print_results()
        void save_results(string)
        void save_individual_results(string, bool)

cdef class QMC:
    cdef CXX_QMC *thisptr
    cdef str params_name
    cdef float beta
    cdef float t_max
    def __cinit__(self, str params_name, float beta, float t_max):
        self.params_name, self.beta, self.t_max = params_name, beta, t_max
    def numerical_input(self, float dt, int real_grid_size, float dtau, int imag_grid_size):
        del self.thisptr
        self.thisptr = new CXX_QMC(self.params_name.encode('utf-8'), self.beta, self.t_max,
                               dt, real_grid_size, dtau, imag_grid_size)
    def analytical_input(self, eps_0_up, eps_up, eps_0_down, eps_down, U_0, U,
                         V_0_up, V_up, eps_bath_0_up, eps_bath_up, V_0_down, V_down, eps_bath_0_down, eps_bath_down):
        del self.thisptr
        self.thisptr = new CXX_QMC(self.params_name.encode('utf-8'), self.beta, self.t_max,
                                   eps_0_up, eps_up, eps_0_down, eps_down, U_0, U,
                                   V_0_up, V_up, eps_bath_0_up, eps_bath_up,
                                   V_0_down, V_down, eps_bath_0_down, eps_bath_down)

    def mixed_input(self, float dt, int real_grid_size, float dtau, int imag_grid_size,
                    float eps_0_up, float eps_up, float eps_0_down, float eps_down, float U_0, float U,
                    V_0_up, V_up, eps_bath_0_up, eps_bath_up):
        del self.thisptr
        self.thisptr = new CXX_QMC(self.params_name.encode('utf-8'), self.beta, self.t_max,
                                   dt, real_grid_size, dtau, imag_grid_size,
                                   eps_0_up, eps_up, eps_0_down, eps_down, U_0, U,
                                   V_0_up, V_up, eps_bath_0_up, eps_bath_up,
                                   True)
    def __dealloc__(self):
        del self.thisptr
    def set_n_imag_times(self, int n):
        self.thisptr.set_n_imag_times(n)
    def set_n_real_times(self, int n):
        self.thisptr.set_n_real_times(n)
    def set_n_matsubara(self, int n):
        self.thisptr.set_n_matsubara(n)
    def set_cycle_size(self, int n):
        self.thisptr.set_block_size(n)
    def set_n_blocks(self, int n):
        self.thisptr.set_n_blocks(n)
    def run_measurement(self, long N):
        self.thisptr.run_measurement(N)
    def run_warm_up(self, long N):
        self.thisptr.run_warm_up(N)
    def run_fixed_time_measurement(self, float minutes):
        self.thisptr.run_measurement(minutes)
    def run_fixed_time_warm_up(self, float minutes):
        self.thisptr.run_warm_up(minutes)
    def collect_results(self):
        self.thisptr.collect_results()
    def reset_measurements(self):
        self.thisptr.reset_measurements()
    def print_results(self):
        self.thisptr.print_results()
    def save_results(self, str output_prefix):
        self.thisptr.save_results(output_prefix.encode('utf-8'))
    def save_individual_results(self, str output_prefix, bool truncate = False):
        self.thisptr.save_individual_results(output_prefix.encode('utf-8'), truncate)

