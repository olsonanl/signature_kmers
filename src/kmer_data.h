#ifndef _kmer_data_h
#define _kmer_data_h

#include <cstdint>
#include <cstddef>
#include <limits>
#include <array>
#include <unordered_map>
#include <iostream>
#include <iterator>
#include <algorithm>


/*!
 * Function indexes are use to reference the entries in function.index
 * which represent function strings assigned to proteins.
 */
typedef uint16_t FunctionIndex;

/**
 * Value representing a missing or undefined function.
 */
const FunctionIndex UndefinedFunction = std::numeric_limits<FunctionIndex>::max();

/**
 * OTU indexes are use to reference the entries in otu.index
 * which represent OTUs associated with kmers.
 */
typedef uint16_t OTUIndex;

/**
 * Value representing a missing or undefined function.
 */
const FunctionIndex UndefinedOTU = std::numeric_limits<FunctionIndex>::max();

template <int K>
using Kmer = std::array<char, K>;

template <size_t K>
std::ostream &operator<<(std::ostream &os, const Kmer<K> &k)
{
    std::copy(k.begin(), k.end(),
	      std::ostream_iterator<char>(os, ""));
    return os;
}

/*! @brief Hash function on kmers.
 */
/*
namespace std {
template<class T, size_t N>
struct hash<std::array<T, N>> {
    auto operator() (const array<T, N>& key) const {
	hash<T> hasher;
	size_t result = 0;
	for(size_t i = 0; i < N; ++i) {
	    result = result * 31 + hasher(key[i]); // ??
	}
	return result;
    }
};
}
*/

template <int K>
struct tbb_hash {
    tbb_hash() {}
    size_t operator()(const Kmer<K>& k) const {
        size_t h = 0;
	for (auto s: k)
	    h = (h * 17) ^ (unsigned int) s;
	return h;
    }
};

template <int N, typename F>
void for_each_kmer(const std::string &str, F cb) {
    const char *ptr = str.c_str();
    const char *end = ptr  + str.length();
    const char *last_kmer= end - N;

    const char *next_ambig = std::find_if(ptr, end, [](char c) -> bool { return c == '*' || c == 'X'; });
    // std::cerr << "next_ambig=" << next_ambig << "\n";
    std::array<char, N> kmer;
    while (ptr <= last_kmer)
    {
	// std::cerr << "ptr=" << ptr << "\n";
	const char *kend = ptr + N;
	// std::cerr << "kend=" << kend << "\n";
	if (next_ambig != end && kend >= next_ambig)
	{
	    // std::cerr << "hit abmig\n";
	    ptr = next_ambig + 1;
	    next_ambig = std::find_if(ptr, end, [](char c) { return c == '*' || c == 'X'; });
	    continue;
	}
	std::copy(ptr, kend, kmer.data());
	// std::cerr << "cb " << kmer << "\n";
	cb(kmer, ptr - str.c_str());
	ptr++;
    }
}


struct KmerAttributes
{
    FunctionIndex func_index;
    OTUIndex otu_index;
    unsigned short offset;
    unsigned int seq_id;
    unsigned int protein_length;
};

struct StoredKmerData
{
    /*!
      @brief Stored form of kmer data.

      These is the data that is stored in the database for each signature kmer.
    */

    uint16_t  avg_from_end = 0;
    FunctionIndex function_index = UndefinedFunction;
    
    uint16_t mean = 0;
    uint16_t median = 0;
    uint16_t var = 0;
};

inline std::ostream &operator<<(std::ostream &os, const StoredKmerData &c)
{
    os << "(" << c.function_index
       << ", " << c.avg_from_end
       << ", " << c.mean
       << ", " << c.median
       << ", " << c.var
       << ")";
    return os;
}


#endif
