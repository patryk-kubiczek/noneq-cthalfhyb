import csv
import os
import numpy as np
from math import sqrt

def get_matrix(filename):
    matrix = []
    with open(filename, 'r') as csvfile:
        filereader = csv.reader(csvfile, delimiter=' ', skipinitialspace=True)
        for row in filereader:
            matrix.append([float(x) for x in row])
    return np.array(matrix)

def get_complex_matrix(filename):
    matrix = []
    with open(filename, 'r') as csvfile:
         filereader = csv.reader(csvfile, delimiter=' ')
         for row in filereader:
             new_row = []
             for number in row:
                 if number != "":
                    number = number[1:-1]
                    real_x, imag_x = map(lambda x: float(x), number.split(",") )
                    new_row.append(complex(real_x, imag_x))
             matrix.append(new_row)
    return np.array(matrix, dtype=complex)

def get_column(filename):
    column = []
    with open(filename, 'r') as csvfile:
         filereader = csv.reader(csvfile)
         for row in filereader:
            column.append(float(row[0]))
    return np.array(column)

def get_complex_column(filename):
    column = []
    with open(filename, 'r') as csvfile:
        filereader = csv.reader(csvfile, delimiter=' ')
        for row in filereader:
            number = row[1][1:-1]
            real_x, imag_x = map(lambda x: float(x), number.split(",") )
            column.append(complex(real_x, imag_x))
    return np.array(column, dtype=complex)

def get_number(filename):
    with open(filename, 'r') as csvfile:
        filereader = csv.reader(csvfile)
        for row in filereader:
            return float(row[0])

def get_complex_number(filename):
    with open(filename, 'r') as csvfile:
        filereader = csv.reader(csvfile)
        for row in filereader:
            number = row[0][2:-1]
            real_x, imag_x = map(lambda x: float(x), number.split(","))
            return complex(real_x, imag_x)


def get_dim_and_shape(filename):
    with open(filename, 'r') as file:
        first_row = next(csv.reader(file, delimiter=' '))
        dim = int(first_row[0])
        shape = [int(n) for n in first_row[1:dim+1]]
    return dim, shape

def get_n_procs(filename):
    n_procs = 0
    while(os.path.isfile(filename(n_procs))): n_procs += 1
    return n_procs

def get_n_blocks(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        n_blocks = sum(1 for row in reader)
    return n_blocks


def get_benchmark_data(dir, prefix, suffix, rank=None, dtype=float):
    filename = dir + "/" + prefix + '_benchmark_' + suffix + ('_' + rank if rank else '') +'.txt'
    dim, shape = get_dim_and_shape(filename)
    datacols = np.arange(1 + len(shape), 1 + len(shape) + int(np.prod(shape)))
    data = np.loadtxt(filename, dtype=dtype, usecols=datacols).reshape(shape, order='F')
    return data


def get_data(dir, prefix, suffix, dtype=float):
    filename = lambda i : dir + "/" + prefix + '_' + suffix + '_' + str(i) + '.txt'
    dim, shape = get_dim_and_shape(filename(0))
    n_blocks = get_n_blocks(filename(0))
    n_procs = get_n_procs(filename)
    datacols = np.arange(1 + len(shape), 1 + len(shape) + int(np.prod(shape)))
    signcol = datacols[-1] + 1
    Ncol = signcol + 1
    data = np.empty(shape=[n_procs, n_blocks]+shape, dtype=dtype)
    sign = np.empty(shape=[n_procs, n_blocks], dtype=float)
    N = np.empty(shape=[n_procs, n_blocks], dtype=int)
    for i in range(n_procs):
        data[i] = np.loadtxt(filename(i), dtype=dtype, usecols=datacols).reshape([n_blocks]+shape, order='F')
        sign[i] = np.real(np.loadtxt(filename(i), dtype=dtype, usecols=(signcol,)))
        N[i] = np.real(np.loadtxt(filename(i), dtype=dtype, usecols=(Ncol,)))
    return [data, sign, N]

def get_merged_data(dir, prefix, suffix, dtype=float):
    filename = dir + "/" + prefix + '_' + suffix + '.txt'
    dim, shape = get_dim_and_shape(filename)
    n_blocks = get_n_blocks(filename)
    datacols = np.arange(1 + len(shape), 1 + len(shape) + int(np.prod(shape)))
    signcol = datacols[-1] + 1
    Ncol = signcol + 1
    data = np.empty(shape=[1, n_blocks]+shape, dtype=dtype)
    sign = np.empty(shape=[1, n_blocks], dtype=float)
    N = np.empty(shape=[1, n_blocks], dtype=int)
    data[0] = np.loadtxt(filename, dtype=dtype, usecols=datacols).reshape([n_blocks]+shape, order='F')
    sign[0] = np.real(np.loadtxt(filename, dtype=dtype, usecols=(signcol,)))
    N[0] = np.real(np.loadtxt(filename, dtype=dtype, usecols=(Ncol,)))
    return [data, sign, N]

default_suffices = ['occupancies', 'double_occupancy',
            'gf_greater_t_t', 'gf_lesser_t_t', 'gf_tau_t', 'gf_tau',
            'gf_iwn_t', 'gf_iwn',
            'initial_occupancies',
            'mean_complex_sign',
            'mean_order', 'order_histogram',
            'spin_up_density_matrix_elements',
            'times_statistics']

def merge_data(dir, prefix, n_procs, suffices=default_suffices, truncate=True):
    filename = lambda suffix, i : dir + "/" + prefix + '_' + suffix + '_' + str(i) + '.txt'
    for suffix in suffices:
        master_filename = dir + "/" + prefix + '_' + suffix + ".txt"
        with open(master_filename, 'w' if truncate else 'a') as master_file:
            n = 0
            counter = 0
            while n < n_procs:
                try:
                    file = open(filename(suffix, n), 'r')
                    for line in file:
                        master_file.write(line)
                    os.remove(filename(suffix, n))
                    counter += 1
                except:
                    pass
                n += 1
            print(prefix + "_" + suffix + ": " + str(counter) + " processes")
            if counter == 0:
                print("To remove: ", master_filename)
                os.remove(master_filename)

def postprocess_data(data, weights):
    reduced_weights = weights / max(weights)
    tiled_weights = np.tile(reduced_weights, list(data.shape[2:]) + [1, 1])
    tiled_weights = np.transpose(tiled_weights, [-2, -1] + list(range(data.ndim - 2)))
    average = np.average(data, weights=tiled_weights, axis=(0,1))
    if data.dtype == complex:
        error = np.sqrt(np.average(np.absolute(data.real - average.real) ** 2, weights=tiled_weights, axis=(0,1))) \
                + 1j * np.sqrt(np.average(np.absolute(data.imag - average.imag) ** 2, weights=tiled_weights, axis=(0,1)))
    else:
        error = np.sqrt(np.average(np.absolute(data - average) ** 2, weights=tiled_weights, axis=(0,1)))
    normalization = sqrt(np.sum(reduced_weights**2) / np.sum(reduced_weights)**2)
    return average, error * normalization
