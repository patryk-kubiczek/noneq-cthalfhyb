from matplotlib import pyplot as plt

output_extension = '.png'
output_dpi = 180
show = True

import sys, os
sys.path.append(os.path.join(sys.path[0],'python'))

from read_data import *
from math import ceil, pi

input_dir = "input"
output_dir = lambda qmc: "output/" + qmc
prefix = "test_python"

beta = 1
t_max = 0.5
input_t_max = 1


# DICTIONARY OF LABELS
text = {"t" : r"$t$", "tau" : r"$\tau$", "iwn" : r"$i\omega_n$",
        "CTHALFHYB" : "CT-1/2-HYB", "CTHYB" : "CT-HYB",
        "occupancies" : r"$n_{\downarrow}(t)$",
        "double_occupancy" : r"$d(t)$",
        "gf_lesser_t_t" : r"$G_{\downarrow}^{<}(t, 0)$",
        "gf_greater_t_t" : r"$G_{\downarrow}^{>}(t, 0)$",
        "gf_tau" : r"$G_{\downarrow}(\tau)$",
        "gf_tau_t" : r"$G_{\downarrow}(\tau = 0, t)$",
        "gf_iwn" : r"$G_{\downarrow}(i\omega_n)$",
        "gf_iwn_t" : r"$G_{\downarrow}(i\omega_n = \pi / \beta, t)$"}

# Function to evaluate weighted average and standard deviation

# Function to transform benchmark
def generate_arguments(data, argument="t"):
    n_args = len(data)
    if argument != "iwn":
        dt = (t_max if argument == "t" else beta) / (n_args - 1)
        return [i * dt for i in range(n_args)]
    else:
        dwn = 2 * pi / beta
        return [dwn * (i + 0.5) for i in range(n_args)]

def real_time_benchmark(data):
    n_times = len(data)
    dt = input_t_max / (n_times - 1)
    i_max = int(ceil(t_max / input_t_max * n_times) + 0.1)
    times = [i * dt for i in range(i_max+1)]
    return times, data[:i_max+1]

def imag_time_benchmark(data):
    n_times = len(data)
    dtau = beta / (n_times - 1)
    times = [i * dtau for i in range(n_times)]
    return times, data

def matsubara_benchmark(data):
    n_freq = len(data)
    dwn = 2 * pi / beta
    freqs = [dwn * (i + 0.5) for i in range(n_freq)]
    return freqs, data

def benchmark(data, argument="t"):
    if argument == "t":
        return real_time_benchmark(data)
    elif argument == "tau":
        return imag_time_benchmark(data)
    else:
        return matsubara_benchmark(data)

def arg_lim(argument="t"):
    if argument == "t":
        return [0, t_max]
    elif argument == "tau":
        return [0, beta]
    else:
        return [0, 50 * pi / beta]


# OCCUPANCIES

def plot_occ(name, qmc_list=("CTHYB", "CTHALFHYB")):
    try:
        times, exact_result = real_time_benchmark(get_benchmark_data(input_dir, prefix, name)[-1])
        plt.plot(times, exact_result)
    except:
        pass
    for data, sign, n_meas, qmc in [get_data(output_dir(QMC), prefix, name) + [QMC] for QMC in qmc_list]:
        result, error = postprocess_data(data[:,:,-1,:], n_meas)
        plt.errorbar(generate_arguments(result), result, fmt='.', yerr=error, label=text[qmc])
    plt.xlabel(r"$t$")
    plt.ylabel(text[name])
    plt.xlim([0, t_max])
    plt.legend(loc='best')
    plt.savefig("plots/" + name + output_extension, dpi=output_dpi)
    if show: plt.show()
    plt.close()


plot_occ("occupancies")
plot_occ("double_occupancy", qmc_list=["CTHYB"])

# GREEN FUNCTIONS

def plot_gf(name, argument="t", use_second_argument=False, qmc_list=("CTHYB", "CTHALFHYB"), benchmark_name=None):
    i, j = [0, slice(None)] if use_second_argument else [slice(None), 0]
    try:
        times, exact_result = benchmark(get_benchmark_data(input_dir, prefix, benchmark_name if benchmark_name else name,
                                                    dtype=complex)[i,j,-1,-1], argument)
        plt.plot(times, exact_result.real, 'o' if argument == "iwn" else '-')
        plt.plot(times, exact_result.imag, 'o' if argument == "iwn" else '-')
    except:
        pass
    for data, sign, n_meas, qmc in [get_data(output_dir(QMC), prefix, name, dtype=complex) + [QMC] for QMC in qmc_list]:
        result, error = postprocess_data(data[:,:,i,j,-1,-1], n_meas)
        plt.errorbar(generate_arguments(result, argument), result.real, fmt='.', yerr=error[0], label=text[qmc] + " (Re)")
        plt.errorbar(generate_arguments(result, argument), result.imag, fmt='.', yerr=error[1], label=text[qmc] + " (Im)")
        plt.xlabel(text[argument])
    plt.ylabel(text[name])
    plt.xlim(arg_lim(argument))
    plt.legend(loc='best')
    plt.savefig("plots/" + name + output_extension, dpi=output_dpi)
    if show: plt.show()
    plt.close()

plot_gf("gf_lesser_t_t")
plot_gf("gf_greater_t_t")
plot_gf("gf_tau", argument="tau")
plot_gf("gf_tau_t", use_second_argument=True)
plot_gf("gf_iwn", argument="iwn", benchmark_name="gf_iwn_t")
plot_gf("gf_iwn_t", use_second_argument=True)
