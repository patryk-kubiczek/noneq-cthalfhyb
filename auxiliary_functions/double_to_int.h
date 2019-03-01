#pragma once

#include <cassert>

inline int int_floor(double x){
    assert(x >= 0);
    return static_cast<int>(x);
}

inline int int_ceil(double x){
    assert(x >= 0);
    return static_cast<int>(x) + 1;
}

inline int int_round(double x){
    assert(x >= 0);
    x = std::round(x);
    return static_cast<int>(x + 0.1);
}
