#ifndef PAGERANK_PAGERANK_H
#define PAGERANK_PAGERANK_H

#include <cstddef>

using std::size_t;

struct IPagerank {
    virtual IPagerank();
    virtual explicit IPagerank(size_t, double);
    virtual ~IPagerank() = 0;
    virtual IPagerank* shallow_copy() = 0;
    virtual IPagerank& operator=(IPagerank const&);
    virtual bool operator()(size_t, size_t) const = 0;
    virtual bool operator()(size_t, size_t, bool) = 0;
    virtual double const& operator[](size_t i) const = 0;
    virtual double& operator[](size_t) = 0;
    virtual IPagerank& operator()(size_t);
    virtual size_t links_from_page(size_t) const = 0;
    double const& d() const;
    double& d();
    size_t const& pages() const;

    static double const D_DEFAULT = 0.85;
private:
    double d_;
    size_t pages_;
};

struct Pagerank : IPagerank {
    Pagerank();
    explicit Pagerank(size_t, double);
    ~Pagerank();
    Pagerank(Pagerank const&);
    Pagerank& operator=(Pagerank const&);
    Pagerank* shallow_copy();
    bool operator()(size_t, size_t) const;
    bool operator()(size_t, size_t, bool);
    double const& operator[](size_t i) const;
    double& operator[](size_t);
    Pagerank& operator()(size_t);
    size_t links_from_page(size_t) const;

private:
    static int const INTBITS = 8 * sizeof(int);
    size_t matrix_size;
    double * weights;
    int * matrix;
};

struct RarefiedPagerank : IPagerank {
    RarefiedPagerank();
    explicit RarefiedPagerank(size_t, size_t, double);
    ~RarefiedPagerank();
    RarefiedPagerank(RarefiedPagerank const& other);
    RarefiedPagerank& operator=(RarefiedPagerank const&);
    RarefiedPagerank* shallow_copy();
    bool operator()(size_t, size_t) const;
    bool operator()(size_t, size_t, bool);
    double const& operator[](size_t i) const;
    double& operator[](size_t);
    RarefiedPagerank& operator()(size_t);
    size_t links_from_page(size_t) const;

private:
    size_t indexes_size;
    double * weights;
    size_t * indexes;
};

#endif //PAGERANK_PAGERANK_H
