#include "nudb_kmer_db.h"
#include "cmph_kmer.h"
#include "call_functions.h"
#include "fasta_parser.h"
#include "seq_id_map.h"
#include "calc_natural_breaks.h"

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
#include <execution>
#include <stdexcept>
#include <vector>

using namespace calc_natural_breaks;

/*!

  @mainpage kmers-matrix-distance

  Compute the all-to-all distance matrix for the given sequences.

  For each sequence, we compute the set of signature kmers.
  Then for each pair of sequences, we compute the size of the intersection of these sets.

*/

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct program_parameters
{
    fs::path data_dir;
    fs::path fasta_file;
    fs::path output_file;
    int min_hits = 3;
    bool debug_hits = false;
    bool verbose = false;
    int n_threads = 1;
};

void process_options(int argc, char **argv, program_parameters &params)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " data-dir input-file\nAllowed options";

    po::options_description desc(x.str());
    desc.add_options()
	("data-dir,d", po::value<fs::path>(&params.data_dir), "Data directory")
	("input-file,i", po::value<fs::path>(&params.fasta_file), "Input fasta file")
	("output-file,o", po::value<fs::path>(&params.output_file), "Output file")
	("min-hits", po::value<int>(&params.min_hits), "Minimum shared kmer hits to emit a match")
	("n-threads,j", po::value<int>(&params.n_threads), "Number of threads")
	("debug-hits", po::bool_switch(&params.debug_hits), "Debug kmer hits")
	("verbose", po::bool_switch(&params.verbose), "Verbose mode")
	("help,h", "show this help message");

    po::positional_options_description pos;
    pos.add("data-dir", 1)
	.add("input-file", -1);
    
    po::variables_map vm;

    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pos).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	exit(0);
    }
}

struct Counter
{
    struct value_type { template<typename T> value_type(const T&) { } };
    void push_back(const value_type&) { ++count; }
    size_t count = 0;
};

int main(int argc, char **argv)
{
    program_parameters params;
    process_options(argc, argv, params);

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, params.n_threads);

    auto db_base = params.data_dir / "kmer_data";

    using DbType = CmphKmerDb<StoredKmerData, 8>;

//    NuDBKmerDb<StoredKmerData, 8> nudb(db_base);
    DbType nudb(db_base);

    if (!nudb.exists())
    {
	std::cerr << "Database " << db_base << " does not exist\n";
	exit(1);
    }
    nudb.open();
    FunctionCaller<DbType> caller(nudb, params.data_dir / "function.index");

    SeqIdMap idmap;

    /*
     * kmer_hit_map maps from a kmer to the set of IDs containing that kmer
     */
    tbb::concurrent_unordered_map<Kmer<8>, tbb::concurrent_unordered_set<int>, tbb_hash<8>> kmer_hit_map;

    auto hit_cb = [&kmer_hit_map, &idmap](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
	// std::cerr << id << " " << seqlen << " " << kd << "\n";


	/*
	 * Discard any hit that is outside either 2 standard deviations from the mean
	 * (If we have a variance reported) or outside 20% of the reference sequence length.
	 */

	int idx = idmap.lookup_id(id);
	
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

    if (params.verbose)
	std::cerr << "Start fasta load\n";

    tbb::concurrent_unordered_map<std::string, size_t> prot_sizes;
    auto call_cb = [&prot_sizes](const std::string &id, const std::string &func, FunctionIndex func_index, float score, size_t prot_len) {
	prot_sizes.insert(std::make_pair(id, prot_len));
    };

    fs::ifstream ifstr(params.fasta_file);

    caller.ignore_hypothetical(true);
    caller.process_fasta_stream_parallel(ifstr, hit_cb, call_cb, idmap);
	    
    ifstr.close();

    std::cerr << "kmer_hit_map size " << kmer_hit_map.size() << "\n";

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

    std::cerr << "write output\n";
    for (auto &ent1: seq_dist)
    {
	auto &id1 = ent1.first;
	auto &l1 = ent1.second;
	const std::string &seq1 = idmap.lookup_index(id1);
	for (auto ent2: l1)
	{
	    auto &id2 = ent2.first;
	    const std::string &seq2 = idmap.lookup_index(id2);
	    int count = ent2.second;
	    std::cout << seq1 << "\t" << seq2 << "\t" << count << "\n";
	}
    }
}

