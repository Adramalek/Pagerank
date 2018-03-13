#include <algorithm>
#include <cstring>
#include <cassert>
#include <malloc.h>
#include "pagerank.h"
#include "utils.h"


// ------------ IPagerank -----------

template <size_t P>
IPagerank<P>::IPagerank() : d_(IPagerank<P>::D_DEFAULT) { }

template <size_t P>
IPagerank<P>::IPagerank(double d)
        : d_(d),
          weights_(new double[pages])
{
    std::fill(weights_, weights_+P, 1.0);
}

template <size_t P>
IPagerank<P>::IPagerank(IPagerank<P> const &other)
        : d_(other.d_),
          weights_(P?(double*)memcpy(new double[P], other.weights_, P * sizeof(double)): nullptr)
{}

template <size_t P>
IPagerank<P>::~IPagerank() {
    delete [] weights_;
}

template <size_t P>
IPagerank<P>& IPagerank<P>::operator=(IPagerank<P> const & other) {
    if (this != &other){
        d_ = other.d_;
        delete [] weights_;
        weights_ = P ? (double*) memcpy(new double[P], other.weights_, P * sizeof(double)) : nullptr;
    }
    return *this;
}

template <size_t P>
double const& IPagerank<P>::d() const {
    return d_;
}

template <size_t P>
double& IPagerank<P>::d() {
    return const_cast<double&>(const_cast<IPagerank<P> const&>(*this).d());
}

template <size_t P>
size_t const& IPagerank<P>::pages() const {
    return P;
}

template <size_t P>
IPagerank<P>& IPagerank<P>::operator()(size_t cycles) {
    while(cycles--){
        double* ranks = new double[P];
        for (size_t i = 0; i < P; ++i){
            ranks[i] = 0.0;
            for (size_t j = 0; j < P; ++j){
                if (i == j) continue;
                if ((*this)(j, i)){
                    ranks[i] += (*this)[j]/links_from_page(j);
                }
            }
            ranks[i] += (ranks[i] *= d(), 1.0 - d());
        }
        for (size_t i = 0; i < P; ++i){
            (*this)[i] += ranks[i];
        }
        delete [] ranks;
    }
    return *this;
}

template <size_t P>
double const* IPagerank<P>::weights() const {
    return weights_;
}

template <size_t P>
double const& IPagerank<P>::operator[](size_t i) const {
    assert(i < P && weights_ != nullptr);
    return (*this)[i];
}

template <size_t P>
double& IPagerank<P>::operator[](size_t i){
    assert(weights_ != nullptr);
    return const_cast<double &>(const_cast<IPagerank<P> const&>(*this)[i]);
}

// ------------ Pagerank -----------

template <size_t P>
Pagerank<P>::Pagerank()
        : IPagerank<P>(),
          matrix_size(0),
          matrix_(nullptr)
{}

template <size_t P>
explicit Pagerank<P>::Pagerank(double d)
        : IPagerank<P>(d),
          matrix_size((P*P) / INTBITS)
{
    size_t additional = (P*P) % INTBITS;
    matrix_size +=  (additional ? 1 : 0);
    matrix_ = new int[matrix_size];
    std::fill(matrix_, matrix_+matrix_size, 0);
}

template <size_t P>
virtual Pagerank<P>::~Pagerank(){
    delete [] matrix_;
    matrix_ = nullptr;
}

template <size_t P>
Pagerank<P>::Pagerank(Pagerank<P> const& other)
        : IPagerank<P>(other),
          matrix_size(other.matrix_size),
          matrix_(matrix_size ?
                 (int*)memcpy(new int[matrix_size], other.matrix_, sizeof(int)*matrix_size) :
                 nullptr)
{}

template <size_t P>
Pagerank<P>& Pagerank<P>::operator=(Pagerank<P> const& other){
    if (this != &other){
        delete [] matrix_;
        matrix_ = nullptr;

        static_cast<IPagerank<P>&>(*this) = other;
        matrix_size = other.matrix_size;
        matrix_ = matrix_size ?
                 (int*) memcpy(new int[matrix_size], other.matrix_, sizeof(int)*matrix_size) :
                 nullptr;
    }
    return *this;
}

template <size_t P>
Pagerank<P>* Pagerank<P>::shallow_copy(){
    Pagerank<P>* copy = new Pagerank<P>();
    static_cast<IPagerank<P>&>(*copy) = *this;
    copy->matrix_size = matrix_size;
    copy->matrix_ = matrix_;
    return copy;
}

template <size_t P>
bool Pagerank<P>::operator()(size_t i, size_t j) const {
    assert(i < P && j < P && matrix_ != nullptr);
    size_t ind = i*P+j;
    return bool(matrix_[ind/INTBITS] & (1 << (ind % INTBITS)));
}

template <size_t P>
bool Pagerank<P>::operator()(size_t i, size_t j, bool value) {
    assert(i < P && j < P && matrix_ != nullptr);
    size_t ind = i*P+j;
    if (value)
        matrix_[ind/INTBITS] |= (1 << (ind % INTBITS));
    else
        matrix_[ind/INTBITS] &= ~(1 << (ind % INTBITS));
    return true;
}

template <size_t P>
Pagerank<P>& Pagerank<P>::operator()(size_t cycles){
    assert(matrix_ != nullptr);
    return static_cast<Pagerank<P>&>(static_cast<IPagerank<P>&>(*this)(cycles));
}

template <size_t P>
size_t Pagerank<P>::links_from_page(size_t page) const {
    assert(page < P && matrix_ != nullptr);
    size_t start_bit = page*P;
    size_t end_bit = (page+1)*P;
    size_t i1 = start_bit/INTBITS;
    size_t i2 = end_bit/INTBITS;
    if (i1 == i2){
        return count_set_bits_range(matrix_[i1], (start_bit % INTBITS) + 1, end_bit % INTBITS);
    } else {
        return count_set_bits_range(matrix_[i1], start_bit % INTBITS, static_cast<size_t>(INTBITS)) +
               count_set_bits_range(matrix_[i2], 1, end_bit % INTBITS);
    }
}

template <size_t P>
std::bitset<sizeof(int)>** Pagerank<P>::int_matrix() const {
    std::bitset<sizeof(int)> ** result = new std::bitset<sizeof(int)>*[matrix_size];
    for (size_t i = 0; i < matrix_size; ++i){
        result[i] = new std::bitset<sizeof(int)>(static_cast<unsigned long long>(matrix_[i]));
    }
    return result;
}

// ------------ RarefiedPagerank -----------

template <size_t P>
RarefiedPagerank<P>::RarefiedPagerank()
        : IPagerank<P>(),
          indexes_size(0),
          indexes(nullptr)
{}

template <size_t P>
explicit RarefiedPagerank<P>::RarefiedPagerank(size_t non_null_count, double d)
        : IPagerank<P>(d),
          indexes_size(non_null_count),
          indexes(new size_t[indexes_size])
{
    std::fill(indexes, indexes+indexes_size, 0);
}

template <size_t P>
RarefiedPagerank<P>::~RarefiedPagerank()
{
    delete [] indexes;
    indexes = nullptr;
}

template <size_t P>
RarefiedPagerank<P>::RarefiedPagerank(RarefiedPagerank<P> const &other)
        : IPagerank<P>(other),
          indexes_size(other.indexes_size),
          indexes(indexes_size ? (size_t*)memcpy(new size_t[indexes_size], other.indexes, indexes_size * sizeof(size_t)) : nullptr)
{}

template <size_t P>
RarefiedPagerank<P>& RarefiedPagerank<P>::operator=(RarefiedPagerank<P> const & other)
{
    if (this != &other){
        delete [] indexes;
        indexes = nullptr;

        static_cast<IPagerank<P>&>(*this) = other;
        indexes_size = other.indexes_size;
        if (indexes_size){
            indexes = (size_t*) memcpy(new size_t[indexes_size], other.indexes, indexes_size * sizeof(int));
        }
    }
    return *this;
}

template <size_t P>
RarefiedPagerank<P>* RarefiedPagerank<P>::shallow_copy() {
    RarefiedPagerank<P>* copy = new RarefiedPagerank<P>();
    static_cast<IPagerank<P>&>(*copy) = *this;
    copy->indexes_size = indexes_size;
    copy->indexes = indexes;
    return copy;
}

template <size_t P>
bool RarefiedPagerank<P>::operator()(size_t i, size_t j) const {
    assert(i < P && j < P && indexes != nullptr);
    return static_cast<bool>(find(indexes, indexes_size, i * P + j));
}

template <size_t P>
bool RarefiedPagerank<P>::operator()(size_t i, size_t j, bool val) {
    assert(i < P && j < P && indexes != nullptr);
    size_t ind = i * P + j;
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

template <size_t P>
size_t RarefiedPagerank<P>::links_from_page(size_t page) const {
    assert(page < P && indexes != nullptr);
    size_t count = 0;
    for (size_t* ptr = indexes; ptr != indexes+indexes_size; ++ptr){
        if (*ptr/P!=page) continue;
        ++count;
    }
    return count;
}

template <size_t P>
RarefiedPagerank<P>& RarefiedPagerank<P>::operator()(size_t cycles) {
    assert(indexes != nullptr);
    return static_cast<RarefiedPagerank<P>&>(static_cast<IPagerank<P>&>(*this)(cycles));
}

template <size_t P>
std::bitset<sizeof(int)>** RarefiedPagerank<P>::int_matrix() const {
    std::bitset<sizeof(int)> ** result = new std::bitset<sizeof(int)>*[P*P/(8* sizeof(int))+1];
    int buf = 0;

}