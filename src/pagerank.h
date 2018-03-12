#ifndef PAGERANK_PAGERANK_H
#define PAGERANK_PAGERANK_H

#include <cstddef>

using std::size_t;

struct IPagerank {
    virtual ~IPagerank();
    virtual IPagerank* shallow_copy() = 0;
    virtual bool const& operator()(size_t, size_t) const = 0;
    virtual bool operator()(size_t, size_t) = 0;
    virtual bool operator()(size_t, size_t, bool);
    virtual double const& operator[](size_t i) const = 0;
    virtual double& operator[](size_t) = 0;
    virtual IPagerank& operator()(size_t) = 0;
    virtual size_t links_from_page(size_t) const = 0;
};

struct Pagerank : IPagerank {
    Pagerank();
    explicit Pagerank(size_t, double);
    ~Pagerank();
    Pagerank(Pagerank const&);
    Pagerank& operator=(Pagerank const&);
    Pagerank* shallow_copy();
    bool operator()(size_t, size_t) const;
    bool operator()(size_t, size_t);
    bool operator()(size_t, size_t, bool);
    double const& operator[](size_t i) const;
    double& operator[](size_t);
    Pagerank& operator()(size_t);
    size_t links_from_page(size_t) const;

private:
    static int const INTBITS = 8 * sizeof(int);
    size_t pages;
    size_t matrix_size;
    double d;
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
    bool operator()(size_t, size_t);
    bool operator()(size_t, size_t, bool);
    double const& operator[](size_t i) const;
    double& operator[](size_t);
    RarefiedPagerank& operator()(size_t);
    size_t links_from_page(size_t) const;

private:
    static int const INTBITS = 8 * sizeof(int);
    size_t pages;
    size_t indexes_size;
    double d;
    double * weights;
    int * indexes;
};
#endif //PAGERANK_PAGERANK_H
