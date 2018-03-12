#include <algorithm>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include "pagerank.h"
#include "utils.h"

// ------------ Pagerank -----------

Pagerank::Pagerank()
        : pages(0),
          matrix_size(0),
          d(0.0),
          weights(nullptr),
          matrix(nullptr)
{}

explicit Pagerank::Pagerank(size_t pages, double d)
        : pages(pages),
          matrix_size((pages*pages) / INTBITS),
          d(d),
          weights(new double[pages])
{
    size_t additional = (pages*pages) % INTBITS;
    matrix_size +=  (additional ? 1 : 0);
    matrix = new int[matrix_size];
    std::fill(weights, weights+pages, 1.0);
    std::fill(matrix, matrix+matrix_size, 0);
}

virtual Pagerank::~Pagerank(){
    delete [] weights;
    delete [] matrix;
    weights = nullptr;
    matrix = nullptr;
}

Pagerank::Pagerank(Pagerank const& other)
        : pages(other.pages),
          matrix_size(other.matrix_size),
          d(other.d),
          weights(pages ?
                  (double*)memcpy(new double[pages], other.weights, sizeof(double)*pages) :
                  nullptr),
          matrix(matrix_size ?
                 (int*)memcpy(new int[matrix_size], other.matrix, sizeof(int)*matrix_size) :
                 nullptr)
{}

Pagerank& Pagerank::operator=(Pagerank const& other){
    if (this != &other){
        delete [] weights;
        delete [] matrix;
        weights = nullptr;
        matrix = nullptr;

        pages = other.pages;
        matrix_size = other.matrix_size;
        weights = pages ?
                  (double*) memcpy(new double[pages], other.weights, sizeof(double)*pages) :
                  nullptr;
        matrix = matrix_size ?
                 (int*) memcpy(new int[matrix_size], other.matrix, sizeof(int)*matrix_size) :
                 nullptr;
    }
    return *this;
}

Pagerank* Pagerank::shallow_copy(){
    Pagerank* copy = new Pagerank();
    copy->d = d;
    copy->pages = pages;
    copy->matrix_size = matrix_size;
    copy->weights = weights;
    copy->matrix = matrix;
    return copy;
}

bool Pagerank::operator()(size_t i, size_t j) const {
    assert(i < pages && j < pages && matrix != nullptr);
    size_t ind = i*pages+j;
    return matrix[ind/INTBITS] & (1 << (ind % INTBITS));
}

bool Pagerank::operator()(size_t i, size_t j) {
    return const_cast<Pagerank const&>(*this)(i, j);
}

bool Pagerank::operator()(size_t, size_t, bool value) {
    assert(i < pages && j < pages && matrix != nullptr);
    size_t ind = i*pages+j;
    if (value)
        matrix[ind/INTBITS] |= (1 << (ind % INTBITS));
    else
        matrix[ind/INTBITS] &= ~(1 << (ind % INTBITS));
    return true;
}

double const& Pagerank::operator[](size_t i) const {
    assert(i < pages && weights != nullptr);
    return weights[i];
}

double& Pagerank::operator[](size_t i){
    return const_cast<double &>(const_cast<Pagerank const&>(*this)[i]);
}

Pagerank& Pagerank::operator()(size_t cycles){
    assert(weights != nullptr && matrix != nullptr);
    while(cycles--){
        double* ranks = new double[pages];
        for (size_t i = 0; i < pages; ++i){
            ranks[i] = 0.0;
            for (size_t j = 0; j < pages; ++j){
                if (i == j) continue;
                if ((*this)(j, i)){
                    ranks[i] += (*this)[j]/links_from_page(j);
                }
            }
            ranks[i] += (ranks[i] *= d, 1.0 - d);
        }
        for (size_t i = 0; i < pages; ++i){
            weights[i] += ranks[i];
        }
        delete [] ranks;
    }
    return *this;
}

size_t Pagerank::links_from_page(size_t page) const {
    assert(page < pages && matrix != nullptr);
    size_t start_bit = page*pages;
    size_t end_bit = (page+1)*pages;
    size_t i1 = start_bit/INTBITS;
    size_t i2 = end_bit/INTBITS;
    if (i1 == i2){
        return count_set_bits_range(matrix[i1], (start_bit % INTBITS) + 1, end_bit % INTBITS);
    } else {
        return count_set_bits_range(matrix[i1], start_bit % INTBITS, INTBITS) +
               count_set_bits_range(matrix[i2], 1, end_bit % INTBITS);
    }
}

// ------------ RarefiedPagerank -----------

RarefiedPagerank::RarefiedPagerank()
        : pages(0),
          indexes_size(0),
          d(0.0),
          weights(nullptr),
          indexes(nullptr)
{}

explicit RarefiedPagerank::RarefiedPagerank(size_t pages, size_t non_null_count, double d)
        : pages(pages),
          indexes_size(non_null_count),
          d(d),
          weights(new double[pages]),
          indexes(new int[indexes_size])
{
    std::fill(weights, weights+pages, 0.0);
    std::fill(indexes, indexes+indexes_size, 0);
}

RarefiedPagerank::~RarefiedPagerank()
{
    delete [] weights;
    delete [] matrix;
    delete [] indexes;
    weights = nullptr;
    indexes = nullptr;
}

RarefiedPagerank::RarefiedPagerank(RarefiedPagerank const &other)
        : pages(other.pages),
          indexes_size(other.indexes_size),
          d(other.d),
          weights(pages ? (double*)memcpy(new double[pages], other.weights, pages * sizeof(double)) : nullptr),
          indexes(indexes_size ? (int*)memcpy(new int[indexes_size], other.indexes, indexes_size * sizeof(int)) : nullptr)
{}

RarefiedPagerank& RarefiedPagerank::operator=(RarefiedPagerank const & other)
{
    if (this != &other){
        delete [] weights;
        delete [] matrix;
        delete [] indexes;
        weights = nullptr;
        indexes = nullptr;

        pages = other.pages;
        matrix_size = other.matrix_size;
        indexes_size = other.indexes_size;
        if (pages){
            weights = (double*) memcpy(new double[pages], other.weights, pages * sizeof(double));
        }
        if (indexes_size){
            indexes = (int*) memcpy(new int[indexes_size], other.indexes, indexes_size * sizeof(int));
        }
    }
    return *this;
}

RarefiedPagerank* RarefiedPagerank::shallow_copy() {
    RarefiedPagerank* copy = new Pagerank();
    copy->pages = pages;
    copy->indexes_size = indexes_size;
    copy->d = d;
    copy->weights = weights;
    copy->indexes = indexes;
    return copy;
}

bool RarefiedPagerank::operator()(size_t i, size_t j) const {
    assert(i < page && j < page && matrix != nullptr && indexes != nullptr);
    return find(indexes, indexes_size, i * pages + j);
}

bool RarefiedPagerank::operator()(size_t i, size_t j) {
    return const_cast<RarefiedPagerank const&>(*this)(i, j);
}

bool RarefiedPagerank::operator()(size_t i, size_t j, bool val) {
    assert(i < page && j < page && matrix != nullptr && indexes != nullptr);
    size_t ind = i * pages + j;
    size_t elem_ind = find(indexes, indexes_size, ind);
    if (!elem_ind && val){
        int* tmp = (int*) realloc(indexes, (indexes_size+1) * sizeof(int));
        if (tmp == nullptr) return false;
        indexes = tmp;
        indexes[indexes_size++] = ind;
        std::sort(indexes, indexes+indexes_size);
    } else if (elem_ind && !val){
        static size_t null_count = 0;
        if (null_count >= indexes_size){
            int* tmp = (int*) memcpy(new int[indexes_size], indexes, (indexes_size+null_count) * sizeof(int));
            --elem_ind;
            memmove(tmp+elem_ind, tmp+elem_ind+1, (indexes_size+null_count)-elem_ind+1);
            int* tmp_ = (int*) realloc(tmp, indexes_size-1);
            if (tmp_ == nullptr){
                delete [] tmp;
                return false;
            }
            delete [] indexes;
            indexes = tmp_;
            null_count = 0;
            return true;
        } else {
            --elem_ind;
            memmove(indexes+elem_ind, indexes+elem_ind+1, indexes_size-elem_ind+1);
            indexes[indexes_size--] = 0;
            null_count++;
            return true;
        }
    }
}

double const& RarefiedPagerank::operator[](size_t page) const {
    assert(page < pages && weights != nullptr);
    return weights[page];
}

double& RarefiedPagerank::operator[](size_t page) {
    return const_cast<double &>(const_cast<RarefiedPagerank const&>(*this)[page]);
}

size_t RarefiedPagerank::links_from_page(size_t page) const {
    assert(page < pages && indexes != nullptr);
    
}
