#include "nudb_kmer_db.h"
#include "cmph_kmer.h"
#include "call_functions.h"
#include "fasta_parser.h"
#include "path_utils.h"

#include <tbb/global_control.h>
#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <boost/asio/streambuf.hpp>

#include <thread>

#include <stdexcept>
#include <vector>

/*!

  @mainpage kmers-annotate-seqs

  Dropin replacement for pf-annotate-seqs from the pattyfam compute pipeline.

*/

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct program_parameters
{
    fs::path data_dir;
    fs::path genus_data_dir;
    fs::path sequences_dir;
    fs::path calls_file;
    fs::path uncalled_ids_file;
    int n_threads = 1;
};

void process_options(int argc, char **argv, program_parameters &params)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " kmer-data-dir genus-data-dir sequences-dir calls-file uncalled-ids-file\nAllowed options";

    po::options_description desc(x.str());
    desc.add_options()
	("kmer-data-dir,d", po::value<fs::path>(&params.data_dir), "Kmer data directory")
	("genus-data-dir,g", po::value<fs::path>(&params.genus_data_dir), "Genus data directory")
	("sequences-dir", po::value<fs::path>(&params.sequences_dir), "Sequence directory")
	("calls-file", po::value<fs::path>(&params.calls_file), "Output calls file")
	("uncalled-ids-file", po::value<fs::path>(&params.uncalled_ids_file), "Output uncalled IDs file")
	("parallel,j", po::value<int>(&params.n_threads), "Number of threads")
	("help,h", "show this help message");

    po::positional_options_description pos;
    pos.add("kmer-data-dir", 1)
	.add("genus-data-dir", 1)
	.add("sequences-dir", 1)
	.add("calls-file", 1)
	.add("uncalled-ids-file", 1);
    
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

int main(int argc, char **argv)
{
    program_parameters params;
    process_options(argc, argv, params);

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, params.n_threads);

    auto db_base = params.data_dir / "kmer_data";

    using DbType = CmphKmerDb<StoredKmerData, 8>;

    DbType kdb(db_base);

    if (!kdb.exists())
    {
	std::cerr << "Database " << db_base << " does not exist\n";
	exit(1);
    }
    kdb.open();
    FunctionCaller<DbType> caller(kdb, params.data_dir / "function.index");

    auto hit_cb = [](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
    };

    fs::ofstream anno_out(params.calls_file);

    using shared_buf_t = std::shared_ptr<boost::asio::streambuf>;
    tbb::concurrent_bounded_queue<shared_buf_t> output_queue;
    output_queue.set_capacity(100);
    std::thread writer_thread([&output_queue, &anno_out]{
	shared_buf_t buf;
	while (true)
	{
	    output_queue.pop(buf);
	    if (buf->size() == 0)
		break;
	    anno_out << buf;
	}
    });

    // For this app we keep a vector of uncalled IDs that we'll write when we complete
    //
    tbb::concurrent_vector<std::string> uncalled_ids;

    tbb::concurrent_vector<fs::path> ivec;
    populate_path_list(params.sequences_dir, ivec);
    
    tbb::parallel_for(ivec.range(), [&caller, &output_queue, &db_base, &params, &hit_cb, &uncalled_ids](auto inp)
    {
	for (auto input_path: inp)
	{
	    fs::ifstream ifstr(input_path);
	    shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
	    
	    std::ostream bufstr(buf.get());
	    
	    auto call2_cb = [&bufstr, &uncalled_ids](const std::string &id, const std::string &func, FunctionIndex func_index, float score, size_t seq_len)
		{
		    if (func_index == UndefinedFunction)
		    {
			uncalled_ids.push_back(id);
		    }
		    else
		    {
			bufstr << id << "\t" << func << "\t" << func_index << "\t" << score << "\n";
		    }
		};
	    
	    caller.process_fasta_stream(ifstr, hit_cb, call2_cb);
	    
	    ifstr.close();
	    if (buf->size() > 0)
	    {
		output_queue.push(buf);
	    }
	}
    });

    shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
    output_queue.push(buf);

    writer_thread.join();

    std::ofstream uncalled(params.uncalled_ids_file);
    for (auto id: uncalled_ids)
    {
	uncalled << id << "\n";
    }
}

