#pragma once

#include <vector>

template <typename VectorElement>
void insert_to_vector(std::vector<VectorElement> &vec, int position, const VectorElement &el){
    vec.insert(begin(vec) + position, el);
}

template <typename VectorElement>
void erase_from_vector(std::vector<VectorElement> &vec, int position){
    vec.erase(begin(vec) + position);
}

template <typename VectorElement>
auto get_insertion_position(std::vector<VectorElement> &vec, const VectorElement &el){
    return std::distance(vec.begin(), std::upper_bound(vec.begin(), vec.end(), el));
}

template <typename ListElement, typename ListIterator>
void export_to_vector_of_lists(std::vector<std::list<ListElement>> &vec, std::list<ListElement> &list,
                               ListIterator iter){
    vec.emplace_back(std::list<ListElement>());
    vec.back().splice(std::begin(vec.back()), list, iter);
}

template <typename ListElement, typename ListIterator>
void export_to_vector_of_lists(std::vector<std::list<ListElement>> &vec, std::list<ListElement> &list,
                               ListIterator first,  ListIterator last){
    vec.emplace_back(std::list<ListElement>());
    vec.back().splice(std::begin(vec.back()), list, first, last);
}
