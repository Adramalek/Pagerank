#ifndef PAGERANK_PAGERANK_H
#define PAGERANK_PAGERANK_H

#include <cstddef>
#include <bitset>

using std::size_t;

template <size_t Pages>
struct IPagerank {
    virtual IPagerank();
    virtual explicit IPagerank(double);
    virtual IPagerank(IPagerank<Pages> const&);
    virtual ~IPagerank();
    virtual IPagerank<Pages>* shallow_copy() = 0;
    virtual IPagerank<Pages>& operator=(IPagerank<Pages> const&);
    virtual bool operator()(size_t, size_t) const = 0;
    virtual bool operator()(size_t, size_t, bool) = 0;
    virtual double const& operator[](size_t i) const;
    virtual double& operator[](size_t);
    virtual IPagerank<Pages>& operator()(size_t);
    virtual size_t links_from_page(size_t) const = 0;
    double const& d() const;
    double& d();
    size_t const& pages() const;
    double const * weights() const;
    virtual std::bitset<sizeof(int)> ** int_matrix() const = 0;

    static double const D_DEFAULT = 0.85;
private:
    double d_;
    double * weights_;
};

template <size_t Pages>
struct Pagerank : IPagerank<Pages> {
    Pagerank();
    explicit Pagerank(double);
    ~Pagerank();
    Pagerank(Pagerank<Pages> const&);
    Pagerank<Pages>& operator=(Pagerank<Pages> const&);
    Pagerank<Pages>* shallow_copy();
    bool operator()(size_t, size_t) const;
    bool operator()(size_t, size_t, bool);
    Pagerank<Pages>& operator()(size_t);
    size_t links_from_page(size_t) const;
    std::bitset<sizeof(int)> ** int_matrix() const;

private:
    static int const INTBITS = 8 * sizeof(int);
    size_t matrix_size;
    int * matrix_;
};

template <size_t Pages>
struct RarefiedPagerank : IPagerank<Pages> {
    RarefiedPagerank();
    explicit RarefiedPagerank(size_t, double);
    ~RarefiedPagerank();
    RarefiedPagerank(RarefiedPagerank<Pages> const& other);
    RarefiedPagerank<Pages>& operator=(RarefiedPagerank<Pages> const&);
    RarefiedPagerank<Pages>* shallow_copy();
    bool operator()(size_t, size_t) const;
    bool operator()(size_t, size_t, bool);
    RarefiedPagerank<Pages>& operator()(size_t);
    size_t links_from_page(size_t) const;
    std::bitset<sizeof(int)> ** int_matrix() const;

private:
    size_t indexes_size;
    size_t * indexes;
};

#endif //PAGERANK_PAGERANK_H
