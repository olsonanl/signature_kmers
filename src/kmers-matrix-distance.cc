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

    using hit_set = tbb::concurrent_set<Kmer<8>>;
    using id_hit_map = tbb::concurrent_map<std::string, hit_set>;

    SeqIdMap idmap;
    
    /*
      id_hit_sets maps from a sequence identifier to the set of kmer hits in that sequence
    */
    id_hit_map id_hit_sets;

    /*
     * kmer_hit_map maps from a kmer to the set of IDs containing that kmer
     */
    tbb::concurrent_unordered_map<Kmer<8>, tbb::concurrent_unordered_set<int>, tbb_hash<8>> kmer_hit_map;

    auto hit_cb = [&id_hit_sets, &kmer_hit_map, &idmap](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
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

    if (params.verbose)
	std::cerr << "Start fasta load\n";

    tbb::concurrent_vector<std::string> keys;
    tbb::concurrent_map<std::string, int> key_to_index;
    tbb::concurrent_unordered_map<std::string, size_t> prot_sizes;
    tbb::concurrent_vector<double> just_sizes;
    auto call_cb = [&prot_sizes, &just_sizes](const std::string &id, const std::string &func, FunctionIndex func_index, float score, size_t prot_len) {
	prot_sizes.insert(std::make_pair(id, prot_len));
	just_sizes.emplace_back(static_cast<double>(prot_len));
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

    exit(0);

    
    tbb::concurrent_vector<std::string> keys_with_hits;
    std::copy_if(keys.begin(), keys.end(), std::back_inserter(keys_with_hits),
		 [&id_hit_sets](auto &key) {
		     return id_hit_sets[key].size() > 5;
		 });


    if (params.verbose)
    {
	std::cerr << "Sequence count " << id_hit_sets.size() << "\n";
	std::cerr << "Filtered Sequence count " << keys_with_hits.size() << "\n";
	std::cerr << "Start all to all comparison\n";
    }
    std::vector<tbb::concurrent_vector<std::string>> partitions;

    /*
     * If we have 500K or fewer sequence, process as a single set.
     * Otherwise partition into length classes and handle those individually.
     */
    
    if (just_sizes.size() < 500000)
    {
	partitions.resize(1);
	std::copy(keys_with_hits.begin(), keys_with_hits.end(), std::back_inserter(partitions[0]));
    }
    else
    {
	int k = std::ceil(just_sizes.size() / 200000);

	std::vector<double> vec;
	std::copy(just_sizes.begin(),just_sizes.end(),
		  std::back_inserter(vec));
	
	ValueCountPairContainer sortedUniqueValueCounts;
	GetValueCountPairs(sortedUniqueValueCounts, &vec[0], vec.size());
	
	LimitsContainer resultingbreaksArray;
	ClassifyJenksFisherFromValueCountPairs(resultingbreaksArray, k, sortedUniqueValueCounts);

	if (params.verbose)
	{
	    std::cerr << "Reporting results..." << std::endl;
	    for (double breakValue: resultingbreaksArray)
		std::cerr << breakValue << std::endl << std::endl;
	}
	
	auto bucket_of = [&resultingbreaksArray](int len) -> int {
	    for (int i = 0; i < resultingbreaksArray.size(); i++)
	    {
		if (len < resultingbreaksArray[i])
		    return i;
	    }
	    return static_cast<int>(resultingbreaksArray.size());
	};
	
	partitions.resize(resultingbreaksArray.size() + 1);
	for (auto key: keys_with_hits)
	{
	    int len = prot_sizes[key];
	    int bucket = bucket_of(len);
	    partitions[bucket].push_back(key);
	}
    }
    
    auto compare = id_hit_sets.key_comp();

    if (params.verbose)
	std::cerr << "done setting up key list\n";

    /*
     * Set up output thread
     */
    std::streambuf *sbuf;
    fs::ofstream ofstr;
    if (params.output_file.empty())
    {
	sbuf = std::cout.rdbuf();
    }
    else
    {
	ofstr.open(params.output_file);
	sbuf = ofstr.rdbuf();
    }
    std::ostream anno_out(sbuf);

    using shared_buf_t = std::shared_ptr<boost::asio::streambuf>;
    tbb::concurrent_bounded_queue<shared_buf_t> output_queue;
    output_queue.set_capacity(100);
    std::thread writer_thread([&output_queue, &anno_out]{
	shared_buf_t buf;
	while (true)
	{
	    output_queue.pop(buf);
	    if (buf->size() == 0)
	    {
		std::cerr << "exiting thread\n";
		break;
	    }
	    anno_out << buf;
	}
    });

    for (auto part: partitions)
    {
	std::cerr << "Start partition size=" << part.size() << "\n";
	tbb::parallel_for(tbb::blocked_range2d<int,int>(0, part.size(), 0, part.size()),
			  [&id_hit_sets, &prot_sizes, &part, &output_queue, &params](auto r)
			      {
				  shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
				  std::ostream bufstr(buf.get());
			      
				  for (size_t i=r.rows().begin(); i!=r.rows().end(); ++i)
				  {
				      auto &id1 = part[i];
				      auto &set1 = id_hit_sets[id1];
				      size_t len1 = prot_sizes[id1];
				  
				      for ( size_t j=r.cols().begin(); j!=r.cols().end(); ++j )
				      {
					  if (i <= j)
					      return;
					  auto &id2 = part[j];
					  auto &set2 = id_hit_sets[id2];
					  size_t len2 = prot_sizes[id2];

					  Counter c;

					  std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
								std::back_inserter(c));

					  if (c.count >= params.min_hits)
					  {
					      float score = static_cast<float>(c.count) / static_cast<float>(len1 + len2);
					      bufstr << id1 << "\t" << id2 << "\t" << c.count << "\t" << score << "\n";
					  }
				      }
				  }
				  if (buf->size() > 0)
				  {
				      output_queue.push(buf);
				  }

			      });
    }
					  
    // std::cerr << "all to all done\n";
    shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
    output_queue.push(buf);
    // std::cerr << "wait for writer thread\n";
    writer_thread.join();
    // std::cerr << "Exiting\n";
}

