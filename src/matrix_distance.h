#ifndef _matrix_distance_h
#define _matrix_distance_h

#include "fasta_parser.h"

#include <tbb/global_control.h>
#include <tbb/blocked_range2d.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_map.h>
#include <tbb/concurrent_set.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/asio/streambuf.hpp>

#include <thread>

#include <stdexcept>
#include <vector>


namespace po = boost::program_options;
namespace fs = boost::filesystem;

template <typename Caller>
class MatrixDistance
{
public:


MatrixDistance(Caller &caller, const fs::path &in_file, const fs::path &out_file, bool verbose) 
    : caller_(caller), in_file_(in_file), out_file_(out_file), verbose_(verbose) {
	
    };
    
    struct Counter
    {
	struct value_type { template<typename T> value_type(const T&) { } };
	void push_back(const value_type&) { ++count; }
	size_t count = 0;
    };

    void compute()
    {
	using hit_set = tbb::concurrent_set<Kmer<8>>;
	using id_hit_map = tbb::concurrent_map<std::string, hit_set>;

	id_hit_map id_hit_sets;
	
	auto hit_cb = [&id_hit_sets](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
	    auto iter = id_hit_sets.find(id);
	    if (iter == id_hit_sets.end())
	    {
		auto n = id_hit_sets.emplace(std::make_pair(id, hit_set()));
		n.first->second.insert(kmer);
	    }
	    else
	    {
		iter->second.insert(kmer);
	    }
	};


	tbb::concurrent_unordered_map<std::string, size_t> prot_sizes;
	auto call_cb = [&prot_sizes](const std::string &id, const std::string &func, FunctionIndex func_index, float score, size_t prot_len) {
	    prot_sizes.insert(std::make_pair(id, prot_len));
	};

	fs::ifstream ifstr(in_file_);
	    
	caller_.process_fasta_stream(ifstr, hit_cb, call_cb);
	    
	ifstr.close();

	std::cerr << in_file_ << "Sequence count: " << id_hit_sets.size() << "\n";
	std::cerr << in_file_<< " Start all to all comparison\n";

	auto compare = id_hit_sets.key_comp();

	tbb::concurrent_vector<std::string> keys;
	for (auto ent: id_hit_sets)
	{
	    keys.push_back(ent.first);
	}

	using shared_buf_t = std::shared_ptr<boost::asio::streambuf>;

	tbb::concurrent_vector<shared_buf_t> output_bufs;
	tbb::parallel_for(tbb::blocked_range2d<int,int>(0, keys.size(), 0, keys.size()),
			  [&id_hit_sets, &prot_sizes, &keys, &output_bufs](auto r)
			  {
			      shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
			      std::ostream bufstr(buf.get());
			      
			      for (size_t i=r.rows().begin(); i!=r.rows().end(); ++i)
			      {
				  auto &id1 = keys[i];
				  auto &set1 = id_hit_sets[id1];
				  size_t len1 = prot_sizes[id1];
				  
				  for ( size_t j=r.cols().begin(); j!=r.cols().end(); ++j )
				  {
				      if (i <= j)
					  return;
				      auto &id2 = keys[j];
				      auto &set2 = id_hit_sets[id2];
				      size_t len2 = prot_sizes[id2];

				      Counter c;

				      std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
							    std::back_inserter(c));

				      if (c.count > 0)
				      {
					  float score = static_cast<float>(c.count) / static_cast<float>(len1 + len2);
					  bufstr << id1 << "\t" << id2 << "\t" << c.count << "\t" << score << "\n";
				      }
				  }
			      }
			      if (buf->size() > 0)
			      {
				  output_bufs.push_back(buf);
			      }
			  });
					  
	std::cerr << in_file_ << " all to all done\n";
	fs::ofstream ofstr(out_file_);
	for (auto b: output_bufs)
	{
	    ofstr << b;
	}
    };


private:
    Caller &caller_;
    const fs::path &in_file_;
    const fs::path &out_file_;
    bool verbose_ = false;
};

#endif // _matrix_distance_h
