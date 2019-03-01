#pragma once

#include "mpi.h"
#include "../data_types/Tetracube.h"
#include <complex>
#include <vector>
#include <cassert>
#include <armadillo>

class MPI_Initializer{
public:
    MPI_Initializer(){ MPI_Init(nullptr, nullptr); }
    ~MPI_Initializer(){ MPI_Finalize(); }
};

inline int mpi_rank() {
    int id = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    return id;
}

inline int mpi_size() {
    int n_processes = 1;
    MPI_Comm_size(MPI_COMM_WORLD, &n_processes);
    return n_processes;
}

inline int mpi_barrier() {
    return MPI_Barrier(MPI_COMM_WORLD);
}

template <typename T>
auto mpi_datatype();

template<>
inline auto mpi_datatype<int>(){ return MPI_INT; };
template<>
inline auto mpi_datatype<long int>(){ return MPI_LONG; };
template<>
inline auto mpi_datatype<long long>(){ return MPI_LONG_LONG; };
template<>
inline auto mpi_datatype<unsigned int>(){ return MPI_UNSIGNED; };
template<>
inline auto mpi_datatype<long unsigned int>(){ return MPI_UNSIGNED_LONG; };
template<>
inline auto mpi_datatype<double>(){ return MPI_DOUBLE; };
template<>
inline auto mpi_datatype<long double>(){ return MPI_LONG_DOUBLE; };
template<>
inline auto mpi_datatype<std::complex<double>>(){ return MPI_CXX_DOUBLE_COMPLEX; };
template<>
inline auto mpi_datatype<std::complex<long double>>(){ return MPI_CXX_LONG_DOUBLE_COMPLEX; };

template <typename DataType>
void mpi_reduction(const DataType &local_data, DataType &global_data);

template <typename T>
void mpi_reduction(const T &local_data, T &global_data){
    MPI_Reduce(&local_data, &global_data, 1, mpi_datatype<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_reduction(const std::vector<T> &local_data, std::vector<T> &global_data){
    assert(local_data.size() == global_data.size());
    MPI_Reduce(local_data.data(), global_data.data(), global_data.size(), mpi_datatype<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_reduction(const std::vector<std::vector<T>> &local_data, std::vector<std::vector<T>> &global_data){
    assert(local_data.size() == global_data.size());
    for(int i = 0; i < global_data.size(); ++i){
        assert(local_data[i].size() == global_data[i].size());
        MPI_Reduce(local_data[i].data(), global_data[i].data(), global_data[i].size(), mpi_datatype<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

template <typename T>
void mpi_reduction(const arma::Col<T> &local_data, arma::Col<T> &global_data){
    assert(local_data.size() == global_data.size());
    assert(local_data.size() > 0);
    MPI_Reduce(&local_data[0], global_data.memptr(), global_data.size(), mpi_datatype<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_reduction(const arma::Mat<T> &local_data, arma::Mat<T> &global_data){
    assert(local_data.size() == global_data.size());
    assert(local_data.size() > 0);
    MPI_Reduce(&local_data[0], global_data.memptr(), global_data.size(), mpi_datatype<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
}


template <typename T>
void mpi_reduction(const Tetracube<T> &local_data, Tetracube<T> &global_data){
    assert(local_data.n_slices_1() == global_data.n_slices_1());
    assert(local_data.n_slices_2() == global_data.n_slices_2());
    assert(local_data.slice(0, 0).size() == global_data.slice(0, 0).size());
    MPI_Reduce(&local_data()[0], global_data().memptr(), global_data().size(), mpi_datatype<T>(), MPI_SUM, 0, MPI_COMM_WORLD);
}



template <typename T>
void mpi_all_reduction(const T &local_data, T &global_data){
    MPI_Allreduce(&local_data, &global_data, 1, mpi_datatype<T>(), MPI_SUM, MPI_COMM_WORLD);
}

template <typename T>
void mpi_all_reduction(const std::vector<T> &local_data, std::vector<T> &global_data){
    assert(local_data.size() == global_data.size());
    MPI_Allreduce(local_data.data(), global_data.data(), global_data.size(), mpi_datatype<T>(), MPI_SUM, MPI_COMM_WORLD);
}

template <typename T>
void mpi_all_reduction(const std::vector<std::vector<T>> &local_data, std::vector<std::vector<T>> &global_data){
    assert(local_data.size() == global_data.size());
    for(int i = 0; i < global_data.size(); ++i){
        assert(local_data[i].size() == global_data[i].size());
        MPI_Allreduce(local_data[i].data(), global_data[i].data(), global_data[i].size(), mpi_datatype<T>(), MPI_SUM, MPI_COMM_WORLD);
    }
}

template <typename T>
void mpi_all_reduction(const arma::Col<T> &local_data, arma::Col<T> &global_data){
    assert(local_data.size() == global_data.size());
    assert(local_data.size() > 0);
    MPI_Allreduce(&local_data[0], global_data.memptr(), global_data.size(), mpi_datatype<T>(), MPI_SUM, MPI_COMM_WORLD);
}

template <typename T>
void mpi_all_reduction(const arma::Mat<T> &local_data, arma::Mat<T> &global_data){
    assert(local_data.size() == global_data.size());
    assert(local_data.size() > 0);
    MPI_Allreduce(&local_data[0], global_data.memptr(), global_data.size(), mpi_datatype<T>(), MPI_SUM, MPI_COMM_WORLD);
}

template <typename T>
void mpi_all_reduction(const Tetracube<T> &local_data, Tetracube<T> &global_data){
    assert(local_data.n_slices_1() == global_data.n_slices_1());
    assert(local_data.n_slices_2() == global_data.n_slices_2());
    assert(local_data.slice(0, 0).size() == global_data.slice(0, 0).size());
    MPI_Allreduce(&local_data()[0], global_data().memptr(), global_data().size(), mpi_datatype<T>(), MPI_SUM, MPI_COMM_WORLD);
}



template <typename DataType>
void mpi_broadcast(DataType &data);

template <typename T>
void mpi_broadcast(T &data){
    MPI_Bcast(&data, 1, mpi_datatype<T>(), 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_broadcast(std::vector<T> &data){
    MPI_Bcast(data.data(), data.size(), mpi_datatype<T>(), 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_broadcast(arma::Col<T> &data){
    MPI_Bcast(data.memptr(), data.size(), mpi_datatype<T>(), 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_broadcast(arma::Mat<T> &data){
    MPI_Bcast(data.memptr(), data.size(), mpi_datatype<T>(), 0, MPI_COMM_WORLD);
}

template <typename T>
void mpi_broadcast(Tetracube<T> &data){
    MPI_Bcast(data().memptr(), data().size(), mpi_datatype<T>(), 0, MPI_COMM_WORLD);
}