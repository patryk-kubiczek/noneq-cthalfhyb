#pragma once

template <typename T>
void move_row(arma::Mat<T> &matrix, int old_pos, int new_pos) {
    int swaps = old_pos - new_pos;
    int dir = (swaps > 0 ? -1 : 1);
    int i = 0;
    while(i < abs(swaps)){
        matrix.swap_rows(old_pos + dir * i, old_pos + dir * (i + 1));
        i++;
    }
}

template <typename T>
void move_col(arma::Mat<T> &matrix, int old_pos, int new_pos) {
    int swaps = old_pos - new_pos;
    int dir = (swaps > 0 ? -1 : 1);
    int i = 0;
    while(i < abs(swaps)){
        matrix.swap_cols(old_pos + dir * i, old_pos + dir * (i + 1));
        i++;
    }
}