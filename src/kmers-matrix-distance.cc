#include "nudb_kmer_db.h"
#include "cmph_kmer.h"
#include "call_functions.h"
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
    bool debug_hits = false;
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
	("output-files,o", po::value<fs::path>(&params.output_file), "Output file")
	("n-threads,j", po::value<int>(&params.n_threads), "Number of threads")
	("debug-hits", po::bool_switch(&params.debug_hits), "Debug kmer hits")
	("help,h", "show this help message");

    po::positional_options_description pos;
    pos.add("data-dir", 1)
	.add("input-files", -1);
    
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

    fs::ifstream ifstr(params.fasta_file);
	    
    caller.process_fasta_stream(ifstr, hit_cb, call_cb);
	    
    ifstr.close();

    std::cerr << "Start all to all comparison\n";

    auto compare = id_hit_sets.key_comp();

    tbb::concurrent_vector<std::string> keys;
    for (auto ent: id_hit_sets)
    {
	keys.push_back(ent.first);
    }

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

    tbb::parallel_for(tbb::blocked_range2d<int,int>(0, keys.size(), 0, keys.size()),
		      [&id_hit_sets, &prot_sizes, &keys, &output_queue](auto r)
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
					  bufstr << id1 << "\t" << id2 << "\t" << score << "\n";
				      }
				  }
			      }
			      if (buf->size() > 0)
			      {
				  output_queue.push(buf);
			      }
			  });
					  
    std::cerr << "all to all done\n";
    shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
    output_queue.push(buf);
    std::cerr << "wait for writer thread\n";
    writer_thread.join();
    std::cerr << "Exiting\n";
/*


    tbb::parallel_for(id_hit_sets.range(), [&id_hit_sets, &prot_sizes, &compare](auto r) {
	for (auto s1: r)
	{
	    auto &id1 = s1.first;
	    auto &set1 = s1.second;
	    size_t len1 = prot_sizes[id1];

	    for (auto s2: id_hit_sets)
	    {
		auto &id2 = s2.first;
		
		if (compare(id1, id2))
		{
		    size_t len2 = prot_sizes[id2];
		    auto &set2 = s2.second;
		    
		    std::vector<Kmer<8>> out;
		    std::set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(),
					  std::back_inserter(out));
		    
		    if (!out.empty())
		    {
			float score = static_cast<float>(out.size()) / static_cast<float>(len1 + len2);
			std::cout << id1 << "\t" << id2 << "\t" << score << "\n";
		    }
		}
	    }
	}
    });

*/
}

