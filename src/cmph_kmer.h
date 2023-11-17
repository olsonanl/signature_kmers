#ifndef _cmph_kmer_db_h
#define _cmph_kmer_db_h

/**
 * Kmer database using an flat on-disk mapped data file with cmph-generated perfect hash.
 *
 * Encapsulate all the database ops here.
 */

#include <cmph.h>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>
#include <string>
#include <array>
#include <cstdint>
#include <system_error>
#include <errno.h>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include "kmer_data.h"

namespace fs = boost::filesystem;
namespace ip = boost::interprocess;

template <typename StoredData, int K>
class CmphKmerDb
{
public:
    static constexpr int kmer_size = K;
    static constexpr int KmerSize = K;
    using KData = StoredData;
    using key_type = Kmer<K>;

    CmphKmerDb(const fs::path &file_base)
	: file_base_(file_base)
	, dat_path_(file_base.native() + ".dat")
	, mph_path_(file_base.native() + ".mph")
	, hash_(0)
	{
	    load_hash();
	}

    ~CmphKmerDb() {
	if (hash_)
	{
	    cmph_destroy(hash_);
	}
    }

    void open() {
	map_backing_data();
    }

    /*! Create backing data store for the given hash size.
     */
    void create_backing_data() {
	std::filebuf fbuf;
	size_t sz = hash_size() * sizeof(StoredData);
	fbuf.open(dat_path_, std::ios_base::in | std::ios_base::out
		  | std::ios_base::trunc | std::ios_base::binary);
	fbuf.pubseekoff(sz - 1, std::ios_base::beg);
	fbuf.sputc(0);
    }

    void map_backing_data() {
	mapping_ = ip::file_mapping(dat_path_.native().c_str(), ip::read_write);
	mapped_region_ = ip::mapped_region(mapping_, ip::read_write);
	data_ = (StoredData *) mapped_region_.get_address();

	std::cerr << "Invoke madvise size=" << fs::file_size(dat_path_) << "\n";
	if (madvise(data_, fs::file_size(dat_path_), MADV_POPULATE_READ) != 0)
	{
	    std::cerr << "madvise failed: " << strerror(errno) << "\n";
	}
	std::cerr << "Madvise done\n";
    }

    unsigned int lookup_key(const std::string &key) {
	unsigned int id = cmph_search(hash_, key.c_str(), key.length());
	return id;
    }

    unsigned int lookup_key(const Kmer<K> key) {
	unsigned int id = cmph_search(hash_, key.data(), K);
	return id;
    }

    void load_hash() {
	FILE *fp = fopen(mph_path_.native().c_str(), "rb");
	if (fp == 0)
	{
	    throw std::system_error(errno, std::generic_category(), mph_path_.native());
	}
	hash_ = cmph_load(fp);
	hash_size_ = cmph_size(hash_);
	fclose(fp);
    }
	
    unsigned int hash_size() {
	return hash_size_;
    }


    bool exists() {
	return fs::exists(dat_path_);
    }

    key_type convert_key(const std::string &key) {
	key_type ka;
	if (key.length() != kmer_size)
	    throw std::runtime_error("Invalid kmer size");
	std::copy(key.begin(), key.end(), ka.data());
	return ka;
    }	

    void insert(const std::string &key, const KData &kdata) {
	key_type ka = convert_key(key);;
	int ec;
	insert(ka, kdata, ec);
    }
    void insert(const key_type &key, const KData &kdata, int &ec) {
	unsigned int kidx = lookup_key(key);
	if (kidx < 0 || kidx >= hash_size_)
	{
	    ec = 1;
	    return;
	}
	data_[kidx] = kdata;
    }

    template <typename CB>
    void fetch(const key_type &key, CB cb, int &iec) {
	unsigned int kidx = lookup_key(key);
	if (kidx < 0 || kidx >= hash_size_)
	{
	    iec = 1;
	    return;
	}
	cb(data_[kidx]);
    }
    template <typename CB>
    void fetch(const std::string &key, CB cb, int &iec) {
	fetch(convert_key(key), cb, iec);
    }

private:
    fs::path file_base_;
    fs::path dat_path_, mph_path_;

    cmph_t *hash_;
    unsigned int hash_size_;

    ip::file_mapping mapping_;
    ip::mapped_region mapped_region_;

    StoredData *data_;
};



#endif // _cmph_kmer_db_h
