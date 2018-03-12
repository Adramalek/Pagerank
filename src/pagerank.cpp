#include <algorithm>
#include <cstring>
#include <cassert>
#include <malloc.h>
#include "pagerank.h"
#include "utils.h"


// ------------ IPagerank -----------

IPagerank::IPagerank() : d_(IPagerank::D_DEFAULT), pages_(0) { }

IPagerank::IPagerank(size_t pages, double d) : d_(d), pages_(pages){}

IPagerank& IPagerank::operator=(IPagerank const & other) {
    if (this != &other){
        d_ = other.d_;
        pages_ = other.pages_;
    }
    return *this;
}

double const& IPagerank::d() const {
    return d_;
}

double& IPagerank::d() {
    return const_cast<double&>(const_cast<IPagerank const&>(*this).d());
}

size_t const& IPagerank::pages() const {
    return pages_;
}

IPagerank& IPagerank::operator()(size_t cycles) {
    while(cycles--){
        double* ranks = new double[this->pages()];
        for (size_t i = 0; i < this->pages(); ++i){
            ranks[i] = 0.0;
            for (size_t j = 0; j < this->pages(); ++j){
                if (i == j) continue;
                if ((*this)(j, i)){
                    ranks[i] += (*this)[j]/links_from_page(j);
                }
            }
            ranks[i] += (ranks[i] *= d(), 1.0 - d());
        }
        for (size_t i = 0; i < this->pages(); ++i){
            (*this)[i] += ranks[i];
        }
        delete [] ranks;
    }
    return *this;
}

// ------------ Pagerank -----------

Pagerank::Pagerank()
        : IPagerank(),
          matrix_size(0),
          weights(nullptr),
          matrix(nullptr)
{}

explicit Pagerank::Pagerank(size_t pages, double d)
        : IPagerank(pages, d),
          matrix_size((pages*pages) / INTBITS),
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
        : IPagerank(other.pages(), other.d()),
          matrix_size(other.matrix_size),
          weights(other.pages() ?
                  (double*)memcpy(new double[other.pages()], other.weights, sizeof(double)*other.pages()) :
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

        static_cast<IPagerank&>(*this) = other;
        matrix_size = other.matrix_size;
        weights = this->pages() ?
                  (double*) memcpy(new double[this->pages()], other.weights, sizeof(double)*this->pages()) :
                  nullptr;
        matrix = matrix_size ?
                 (int*) memcpy(new int[matrix_size], other.matrix, sizeof(int)*matrix_size) :
                 nullptr;
    }
    return *this;
}

Pagerank* Pagerank::shallow_copy(){
    Pagerank* copy = new Pagerank();
    static_cast<IPagerank&>(*copy) = *this;
    copy->matrix_size = matrix_size;
    copy->weights = weights;
    copy->matrix = matrix;
    return copy;
}

bool Pagerank::operator()(size_t i, size_t j) const {
    assert(i < this->pages() && j < this->pages() && matrix != nullptr);
    size_t ind = i*this->pages()+j;
    return bool(matrix[ind/INTBITS] & (1 << (ind % INTBITS)));
}

bool Pagerank::operator()(size_t i, size_t j, bool value) {
    assert(i < this->pages() && j < this->pages() && matrix != nullptr);
    size_t ind = i*this->pages()+j;
    if (value)
        matrix[ind/INTBITS] |= (1 << (ind % INTBITS));
    else
        matrix[ind/INTBITS] &= ~(1 << (ind % INTBITS));
    return true;
}

double const& Pagerank::operator[](size_t i) const {
    assert(i < this->pages() && weights != nullptr);
    return weights[i];
}

double& Pagerank::operator[](size_t i){
    return const_cast<double &>(const_cast<Pagerank const&>(*this)[i]);
}


Pagerank& Pagerank::operator()(size_t cycles){
    assert(weights != nullptr && matrix != nullptr);
    return static_cast<Pagerank&>(static_cast<IPagerank&>(*this)(cycles));
}

size_t Pagerank::links_from_page(size_t page) const {
    assert(page < this->pages() && matrix != nullptr);
    size_t start_bit = page*this->pages();
    size_t end_bit = (page+1)*this->pages();
    size_t i1 = start_bit/INTBITS;
    size_t i2 = end_bit/INTBITS;
    if (i1 == i2){
        return count_set_bits_range(matrix[i1], (start_bit % INTBITS) + 1, end_bit % INTBITS);
    } else {
        return count_set_bits_range(matrix[i1], start_bit % INTBITS, static_cast<size_t>(INTBITS)) +
               count_set_bits_range(matrix[i2], 1, end_bit % INTBITS);
    }
}

// ------------ RarefiedPagerank -----------

RarefiedPagerank::RarefiedPagerank()
        : IPagerank(),
          indexes_size(0),
          weights(nullptr),
          indexes(nullptr)
{}

explicit RarefiedPagerank::RarefiedPagerank(size_t pages, size_t non_null_count, double d)
        : IPagerank(pages, d),
          indexes_size(non_null_count),
          weights(new double[pages]),
          indexes(new size_t[indexes_size])
{
    std::fill(weights, weights+pages, 0.0);
    std::fill(indexes, indexes+indexes_size, 0);
}

RarefiedPagerank::~RarefiedPagerank()
{
    delete [] weights;
    delete [] indexes;
    weights = nullptr;
    indexes = nullptr;
}

RarefiedPagerank::RarefiedPagerank(RarefiedPagerank const &other)
        : IPagerank(other.pages(), other.d()),
          indexes_size(other.indexes_size),
          weights(other.pages() ? (double*)memcpy(new double[other.pages()], other.weights, other.pages() * sizeof(double)) : nullptr),
          indexes(indexes_size ? (size_t*)memcpy(new size_t[indexes_size], other.indexes, indexes_size * sizeof(size_t)) : nullptr)
{}

RarefiedPagerank& RarefiedPagerank::operator=(RarefiedPagerank const & other)
{
    if (this != &other){
        delete [] weights;
        delete [] indexes;
        weights = nullptr;
        indexes = nullptr;

        static_cast<IPagerank&>(*this) = other;
        indexes_size = other.indexes_size;
        if (this->pages()){
            weights = (double*) memcpy(new double[this->pages()], other.weights, this->pages() * sizeof(double));
        }
        if (indexes_size){
            indexes = (size_t*) memcpy(new size_t[indexes_size], other.indexes, indexes_size * sizeof(int));
        }
    }
    return *this;
}

RarefiedPagerank* RarefiedPagerank::shallow_copy() {
    RarefiedPagerank* copy = new RarefiedPagerank();
    copy->pages = pages;
    copy->indexes_size = indexes_size;
    copy->weights = weights;
    copy->indexes = indexes;
    return copy;
}

bool RarefiedPagerank::operator()(size_t i, size_t j) const {
    assert(i < pages() && j < pages() && indexes != nullptr);
    return static_cast<bool>(find(indexes, indexes_size, i * pages() + j));
}

bool RarefiedPagerank::operator()(size_t i, size_t j, bool val) {
    assert(i < pages() && j < pages() && indexes != nullptr);
    size_t ind = i * pages() + j;
    size_t elem_ind = find(indexes, indexes_size, ind);
    static size_t null_count = 0;
    if (!elem_ind && val){
        if (!null_count) {
            size_t *tmp = (size_t *) realloc(indexes, (indexes_size + 1) * sizeof(size_t));
            if (tmp == nullptr) return false;
            indexes = tmp;
        } else {
            --null_count;
        }
        indexes[indexes_size++] = ind;
        std::sort(indexes, indexes + indexes_size);
    } else if (elem_ind && !val){
        if (null_count >= indexes_size){
            size_t* tmp = (size_t*) memcpy(new size_t[indexes_size], indexes, (indexes_size+null_count) * sizeof(size_t));
            --elem_ind;
            memmove(tmp+elem_ind, tmp+elem_ind+1, (indexes_size+null_count)-elem_ind+1);
            size_t* tmp_ = (size_t*) realloc(tmp, (indexes_size-1) * sizeof(size_t));
            if (tmp_ == nullptr){
                delete [] tmp;
                return false;
            }
            delete [] indexes;
            indexes = tmp_;
            null_count = 0;
        } else {
            --elem_ind;
            memmove(indexes+elem_ind, indexes+elem_ind+1, indexes_size-elem_ind+1);
            indexes[indexes_size--] = 0;
            null_count++;
        }
    }
    return true;
}

double const& RarefiedPagerank::operator[](size_t page) const {
    assert(page < pages() && weights != nullptr);
    return weights[page];
}

double& RarefiedPagerank::operator[](size_t page) {
    return const_cast<double &>(const_cast<RarefiedPagerank const&>(*this)[page]);
}

size_t RarefiedPagerank::links_from_page(size_t page) const {
    assert(page < pages() && indexes != nullptr);
    size_t count = 0;
    for (size_t* ptr = indexes; ptr != indexes+indexes_size; ++ptr){
        if (*ptr/pages()!=page) continue;
        ++count;
    }
    return count;
}

RarefiedPagerank& RarefiedPagerank::operator()(size_t cycles) {
    assert(weights != nullptr && indexes != nullptr);
    return static_cast<RarefiedPagerank&>(static_cast<IPagerank&>(*this)(cycles));
}
