#pragma once

#include <armadillo>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "../data_types/Tetracube.h"



namespace save_to_row {
    const auto delimiter = ' ';

    template <typename T>
    void print_number(std::ostream &os, const T &number){
        os << number << delimiter;
    }

    template <typename T>
    void print_number(std::ostream& os, const std::complex<T> &number){
        os << number.real() << (number.imag() < 0 ? '-' : '+') << std::abs(number.imag()) << 'j' << delimiter;
    }

    // Format:  # of dimensions,
    //          # of elements dim 0, # of elements dim 1, ..., # of elements dim (D - 1),
    //          element 0, ..., element N - 1,
    //          # additional number 1, ...

    template <typename T>
    void print_data(std::ostream& os, const T &data, const std::vector<double>& additional_data){
        os << '0' << delimiter;
        print_number(os, data);
        for(const auto & number : additional_data) print_number(os, number);
        os << std::endl;
    }

    template <typename T>
    void print_data(std::ostream& os, const arma::Col<T> &data, const std::vector<double>& additional_data){
        os << '1' << delimiter << data.n_rows << delimiter;
        for(const auto & number : data) print_number(os, number);
        for(const auto & number : additional_data) print_number(os, number);
        os << std::endl;
    }

    template <typename T>
    void print_data(std::ostream& os, const arma::Mat<T> &data, const std::vector<double>& additional_data){
        os << '2' << delimiter << data.n_rows << delimiter << data.n_cols << delimiter;
        for(const auto & number : data) print_number(os, number);
        for(const auto & number : additional_data) print_number(os, number);
        os << std::endl;
    }

    template <typename T>
    void print_data(std::ostream& os, const arma::Cube<T> &data, const std::vector<double>& additional_data){
        os << '3' << delimiter << data.n_rows << delimiter << data.n_cols << delimiter << data.n_slices << delimiter;
        for(const auto & number : data) print_number(os, number);
        for(const auto & number : additional_data) print_number(os, number);
        os << std::endl;
    }

    template <typename T>
    void print_data(std::ostream& os, const Tetracube<T> &data, const std::vector<double>& additional_data){
        os << '4' << delimiter << data().n_rows << delimiter << data().n_cols << delimiter
           << data.n_slices_1() << delimiter << data.n_slices_2() << delimiter;
        for(const auto & number : data()) print_number(os, number);
        for(const auto & number : additional_data) print_number(os, number);
        os << std::endl;
    }
    
    template <typename T>
    void save(const std::string &filename, const T &data, const std::vector<double>& additional_data, bool truncate = false){
        std::ofstream myfile;
        myfile.open(filename, std::ios::out | (truncate ? std::ios::trunc : std::ios::app));
        print_data(myfile, data, additional_data);
        myfile.close();
    }

}


