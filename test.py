import sys, os, importlib
from math import cos, pi, sqrt

sys.path.append(os.path.join(os.getenv("HOME"),'Dropbox/Projects/C++/noneq-cthalfhyb/python/'))
sys.path.append(os.path.join(os.getenv("HOME"),'projects/noneq-cthalfhyb/python/'))

# IMPORT QMC

qmc_version = "CTHALFHYB"
if 'QMC' in os.environ:
    qmc_version = os.environ['QMC'] #Get qmc_version from environment variable

qmc_module = importlib.import_module(qmc_version)

from mpi4py import MPI
comm = MPI.COMM_WORLD

# IMPORT INPUT GENERATOR

from InputGenerator import InputGenerator

# PARAMETERS

generate_input = True
if 'GENERATE_INPUT' in os.environ:
    generate_input = bool(int(os.environ['GENERATE_INPUT']))
params_name = "input/test_python"

n_bath = 1 # 10
u = 3 # 0

n_flavor = 2
beta = 1
t_max = 0.5

input_t_max = 1

input_real_grid_size = 100
input_dt = input_t_max / (input_real_grid_size - 1)

real_grid_size = input_real_grid_size
imag_grid_size = int(beta / input_t_max * real_grid_size)
dt = input_t_max / (real_grid_size - 1)
dtau = beta / (imag_grid_size - 1)

voltage = 0
eps_ini = 2
D = 1 # band half-width

def U(t):
    if t == 0:
        return u
    else:
        return u

def eps(a, t):
    if t == 0:
        return -U(t) / 2 + eps_ini
    else:
        return -U(t) / 2

def mu_bath(a, p, t):
    if t == 0:
        return 0
    else:
        if p < n_bath // 2:
            return voltage / 2
        else:
            return -voltage / 2

def Gamma(t):
    if t==0:
        return 1
    else:
        return 1

times = [i * input_dt for i in range(input_real_grid_size)]

if voltage == 0:
    energies = [0]
    if(n_bath > 1):
        energy_spacing = 2 * D / (n_bath - 1)
        energies = [-D + i * energy_spacing for i in range(n_bath)]
else:
    energies = [0, 0]
    if(n_bath > 2):
        energy_spacing = 2 * D / (n_bath / 2 - 1)
        energies = [-D + i * energy_spacing for i in range(int(n_bath / 2))]
        energies += energies


def V(a, p, t):
    return sqrt(2 * D / (pi * n_bath * (1 if voltage == 0 else 0.5)) * Gamma(t))

def energy(a, p):
    return [energies[p] - mu_bath(a, p, t) for t in times]

def hybridization(a, p):
    return [V(a, p, t) for t in times]


bath_energies = [[energy(a, p) for a in range(n_flavor)] for p in range(n_bath)]
bath_hybridizations = [[hybridization(a, p) for a in range(n_flavor)] for p in range(n_bath)]


def local_energy(a):
    return [eps(a, t) for t in times]

local_energies = [local_energy(a) for a in range(n_flavor)]
Hubbard_U = [U(t) for t in times]

# GENERATE INPUT
if comm.Get_rank() == 0 and generate_input:
    os.system("mkdir -p input")
    os.system("mkdir -p output")
    os.system("mkdir -p output/" + qmc_version)

    input_generator = InputGenerator(bath_energies, bath_hybridizations, local_energies, Hubbard_U,
                                     input_dt, input_real_grid_size, beta, input_t_max, real_grid_size, imag_grid_size, params_name)

    input_generator.generate_all_for_determinant()
    input_generator.generate_all_for_one_spin()

    input_generator.generate_benchmark_data()

comm.Barrier()

# QMC PARAMS

output_prefix = "output/" + qmc_version + "/test_python"

n_real = 10
n_imag = 50
n_cycle = 10

minutes_MC = 1

# RUN QMC

qmc_type = "DISCRETE_BATH"
if qmc_version == "CTHYB" and qmc_type == "DISCRETE_BATH":
    qmc_type = "GENERAL"

qmc = qmc_module.QMC(params_name, beta, t_max)
if qmc_type == "GENERAL":
    qmc.numerical_input(dt, real_grid_size, dtau, imag_grid_size)
if qmc_type == "1BATH":
    qmc.analytical_input(eps_0_up=eps(0, 0), eps_up=eps(0, 1.), eps_0_down=eps(1, 0), eps_down=eps(1, 1.), U_0=U(0), U=U(1.),
                         V_0_up=V(0, 0, 0.), V_up=V(0, 0, 1.), eps_bath_0_up=0, eps_bath_up=0,
                         V_0_down=V(1, 0, 0), V_down=V(1, 0, 1.), eps_bath_0_down=0, eps_bath_down=0)
if qmc_type == "DISCRETE_BATH":
    qmc.mixed_input(dt, real_grid_size, dtau, imag_grid_size,
                eps_0_up=eps(0, 0), eps_up=eps(0, 1.), eps_0_down=eps(1, 0), eps_down=eps(1, 1.), U_0=U(0), U=U(1.),
                V_0_up=[V(0, 0, 0.)], V_up=[V(0, 0, 1.)], eps_bath_0_up=[0], eps_bath_up=[0])
    #qmc.set_n_blocks(1)


qmc.set_n_real_times(n_real)
qmc.set_n_imag_times(n_imag)
qmc.set_cycle_size(n_cycle)

qmc.reset_measurements()
qmc.run_fixed_time_warm_up(0.05 * minutes_MC)


n_runs = 10
truncate = True
for i in range(n_runs):
    qmc.reset_measurements()
    qmc.run_fixed_time_measurement(minutes_MC / n_runs)
    qmc.collect_results()
    qmc.print_results()
    qmc.save_results(output_prefix)
    qmc.save_individual_results(output_prefix, truncate=truncate)
    if truncate == True: truncate = False