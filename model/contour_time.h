#pragma once

#include <cassert>
#include <complex>

enum Branch{PLUS, MINUS, IMAG};

struct contour_time{
    Branch branch;
    double t;

    bool operator<(const contour_time &other) const;
    bool operator>(const contour_time &other) const{ return !(*this < other); }
    bool operator==(const contour_time &other) const{ return branch == other.branch && t == other.t; }
    bool operator!=(const contour_time &other) const{ return !(*this == other); }

    double s(double t_max) const;
    static contour_time make_contour_time(double s, double t_max);

    static double distance(const contour_time &ct_end, const contour_time &ct_start, double t_max, double beta);
    std::complex<double> integration_sign() const;

};


inline bool contour_time::operator<(const contour_time &other) const{
    if(this->branch == PLUS){
        if(other.branch == PLUS)
            return this->t < other.t;
        else
            return true;
    }
    else if(this->branch == MINUS){
        if(other.branch == PLUS)
            return false;
        else if(other.branch == MINUS)
            return this->t > other.t;
        else
            return true;
    }
    else{ //this->branch == IMAG
        if(other.branch == IMAG)
            return this->t < other.t;
        else
            return false;
    }
}

inline double contour_time::s(double t_max) const {
    if(branch == PLUS)
        return t;
    else if(branch == MINUS)
        return 2 * t_max - t;
    else
        return 2 * t_max + t;

}

inline double contour_time::distance(const contour_time &ct_end, const contour_time &ct_start, double t_max, double beta){
    double result = ct_end.s(t_max) - ct_start.s(t_max);
    if(result > 0)
        return result;
    else
        return result + 2 * t_max + beta;
}

inline std::complex<double> contour_time::integration_sign() const {
    switch(branch){
        case PLUS : return {1, 0};
        case MINUS : return {-1, 0};
        case IMAG : return {0, -1};
    }
}

inline contour_time contour_time::make_contour_time(double s, double t_max) {
    assert(s >= 0);
    Branch branch;
    double t;
    if(s < t_max){
        branch = PLUS;
        t = s;
    }
    else if(s < 2 * t_max){
        branch = MINUS;
        t = 2 * t_max - s;
    }
    else{
        branch = IMAG;
        t = s - 2 * t_max;
    }
    return {branch, t};
}



