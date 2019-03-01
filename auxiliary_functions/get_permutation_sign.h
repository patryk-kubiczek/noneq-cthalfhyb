#pragma once

#include <vector>

template <typename PermutationFunction>
int get_permutation_sign(PermutationFunction permutation, size_t n){
    // permutation must be callable - int permutation(int i)
    std::vector<bool> visited(n, false);
    int sgn = 1;
    int next = 0;
    int cycle_length = 0;
    for(int i = 0; i < n; ++i){
        if(!visited[i]) {
            next = i;
            cycle_length = 0;
            while (!visited[next]) {
                ++cycle_length;
                visited[next] = true;
                next = permutation(next);
            }
            if (cycle_length % 2 == 0)
                sgn *= -1;
        }
    }
    return sgn;
}