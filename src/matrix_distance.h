#ifndef _matrix_distance_h
#define _matrix_distance_h

#include "fasta_parser.h"
#include "seq_id_map.h"

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
	: caller_(caller), in_files_{in_file}, out_file_(out_file), verbose_(verbose) {
	
    };
    
    MatrixDistance(Caller &caller, const std::vector<fs::path> &in_files, const fs::path &out_file, bool verbose) 
	: caller_(caller), in_files_(in_files), out_file_(out_file), verbose_(verbose) {
	
    };
    
    void compute()
    {
	/*
	 * kmer_hit_map maps from a kmer to the set of IDs containing that kmer
	 */
	tbb::concurrent_unordered_map<Kmer<8>, tbb::concurrent_unordered_set<int>, tbb_hash<8>> kmer_hit_map;

	auto hit_cb = [&kmer_hit_map, this](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
	    // std::cerr << id << " " << seqlen << " " << kd << "\n";


	    /*
	     * Discard any hit that is outside either 2 standard deviations from the mean
	     * (If we have a variance reported) or outside 20% of the reference sequence length.
	     */

	    int idx = idmap_.lookup_id(id);
	
	    double cutoff_b, cutoff_t;
	    double mean = static_cast<double>(kd.mean);
	    double stddev;
	    if (kd.var == 0)
	    {
		stddev = seqlen * 0.1;
	    }
	    else
	    {
		stddev = std::sqrt(static_cast<double>(kd.var));
	    }
	    cutoff_b = mean - stddev * 2.0;
	    cutoff_t = mean + stddev * 2.0;

	    if (seqlen < cutoff_b || seqlen > cutoff_t)
		return;

	    kmer_hit_map[kmer].insert(idx);
	};

	tbb::concurrent_unordered_map<std::string, size_t> prot_sizes;
	auto call_cb = [&prot_sizes](const std::string &id, const std::string &func, FunctionIndex func_index, float score, size_t prot_len) {
	    prot_sizes.insert(std::make_pair(id, prot_len));
	};

	caller_.ignore_hypothetical(true);
	std::string label;
	for (auto in_file: in_files_)
	{
	    if (!fs::is_regular_file(in_file) || fs::is_empty(in_file))
		continue;
	    
	    fs::ifstream ifstr(in_file);
	    caller_.process_fasta_stream_parallel(ifstr, hit_cb, call_cb, idmap_);
	    ifstr.close();
	    if (label == "")
		label = in_file.string();
	    else
	    {
		label += ",";
		label += in_file.string();
	    }
	}

	if (verbose_)
	{
	    std::cerr << label << " Start all to all comparison\n";
	    std::cerr << "kmer_hit_map size " << kmer_hit_map.size() << "\n";
	}

	/*
	 * If we encountered no files to process, just return.
	 */
	if (label == "")
	{
	    if (verbose_)
		std::cerr << "Skip compute " << in_files_[0] << "\n";
	    return;
	}

	/*
	 * seq_dist is a map from id1 to id2 containing the count between the pair.
	 * We only add an entry when id1 < id2.
	 */

	tbb::concurrent_unordered_map<int, tbb::concurrent_unordered_map<int, int>> seq_dist;
    
	tbb::parallel_for(kmer_hit_map.range(), 
			  [&seq_dist](auto r) {
			      for (auto ent: r)
			      {
				  const auto &kmer = ent.first;
				  tbb::concurrent_unordered_set<int> &s = ent.second;

				  for (auto &id1: s)
				  {
				      for (auto &id2: s)
				      {
					  if (id1 < id2)
					  {
					      seq_dist[id1][id2]++;
					  }
				      }
				  }
			      }
			  });

	if (verbose_)
	    std::cerr << label << " all to all done\n";
	fs::ofstream ofstr(out_file_);
	for (auto &ent1: seq_dist)
	{
	    auto &id1 = ent1.first;
	    auto &l1 = ent1.second;
	    const std::string &seq1 = idmap_.lookup_index(id1);
	    size_t len1 = prot_sizes[seq1];
	    for (auto ent2: l1)
	    {
		auto &id2 = ent2.first;
		const std::string &seq2 = idmap_.lookup_index(id2);
		size_t len2 = prot_sizes[seq2];
		int count = ent2.second;
		float score = static_cast<float>(count) / static_cast<float>(len1 + len2);
		ofstr << seq1 << "\t" << seq2 << "\t" << count << "\t" << score << "\n";
	    }
	}

    };


private:
    Caller &caller_;
    const std::vector<fs::path> &in_files_;
    const fs::path &out_file_;
    bool verbose_ = false;
    SeqIdMap idmap_;
};

#endif // _matrix_distance_h
