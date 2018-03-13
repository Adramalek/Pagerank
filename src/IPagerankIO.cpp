#include <bitset>
#include "ipagerankio.h"

IPagerankIO::IPagerankIO(std::string file) : stream(file, std::fstream::in | std::fstream::out)
{}

IPagerankIO::~IPagerankIO() {
    stream.close();
}

IPagerankIO& IPagerankIO::set_file(std::string new_file) {
    stream.close();
    stream.clear();
    stream.open(new_file, std::fstream::in | std::fstream::out);
    return *this;
}

