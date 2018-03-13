//
// Created by User on 12.03.2018.
//

#ifndef PAGERANK_IPAGERANKMATRIXREADER_H
#define PAGERANK_IPAGERANKMATRIXREADER_H

#include <string>
#include <iostream>
#include <fstream>
#include "pagerank.h"

struct IPagerankIO {
    IPagerankIO(std::string file);
    ~IPagerankIO();
    IPagerankIO& set_file(std::string new_file);
    IPagerankIO& operator<<(IPagerankIO&, IPagerank const&);
    IPagerankIO& operator>>(IPagerankIO&, IPagerank&);

private:
    IPagerankIO& operator=(IPagerankIO const& other) = delete;
    IPagerankIO(IPagerankIO const& other) = delete;
    std::fstream stream;
};

void generate_random_matrix(std::string file, size_t pages);

#endif //PAGERANK_IPAGERANKMATRIXREADER_H
