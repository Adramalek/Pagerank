#include "utils.h"

size_t count_set_bits(int n){
    if (!n) return 0;
    return 1 + count_set_bits(n & (n-1));
}

size_t count_set_bits_range(int n, size_t start, size_t end){
    return count_set_bits(n & ((1 << end) - 1) ^ ((1 << (start-1)) - 1));
}

template <class T>
static void find(T* const arr, size_t size, T const& val, size_t* res){
    if (size == 1){
        if (arr[0] == val)
            *res += 1;
        else
            *res = 0;
    } else {
        size_t half = size / 2;
        if (arr[half-1] == val){
            *res += half-1;
            return;
        } else if (arr[half-1] > val){
            find(arr, half, val, res);
        } else {
            *res += half;
            find(arr+half, half, val, res);
        }
    }
}

template <class T>
size_t find(T* const arr, size_t size, T const& val){
    size_t ind = 0;
    find(arr, size, val, &ind);
    return ind;
}