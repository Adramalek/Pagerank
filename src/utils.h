//
// Created by Artur on 11/03/2018.
//

#ifndef PAGERANK_UTILS_H
#define PAGERANK_UTILS_H

#include <cstddef>
using std::size_t;

size_t count_set_bits(int);
size_t count_set_bits_range(int, size_t, size_t);
template <class T>
size_t find(T* const, size_t, T const&);

#endif //PAGERANK_UTILS_H
