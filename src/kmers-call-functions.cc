#include "nudb_kmer_db.h"
#include "cmph_kmer.h"
#include "call_functions.h"
#include "fasta_parser.h"

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

  @mainpage kmers-call-functions

  # Call protein function using signature kmers

  The first parameter is the data directory for the kmer build which
  includes the function.index file which defines the function number
  to function mapping, and the NuDB database files for the saved data.

*/

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct program_parameters
{
    fs::path data_dir;
    std::vector<fs::path> input_files;
    fs::path output_file;
    std::vector<std::string> fasta_dirs;
    bool debug_hits = false;
    bool ignore_hypo = false;
    int n_threads = 1;
};

void process_options(int argc, char **argv, program_parameters &params)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " data-dir input-file [input-file, ...]\nAllowed options";

    po::options_description desc(x.str());
    desc.add_options()
	("data-dir,d", po::value<fs::path>(&params.data_dir), "Data directory")
	("input-files,i", po::value<std::vector<fs::path>>(&params.input_files)->multitoken(), "Input files")
	("output-files,o", po::value<fs::path>(&params.output_file), "Output file")
//	("fasta-dir,F", po::value<std::vector<std::string>>(&params.fasta_dirs)->multitoken(), "Directory of fasta files of protein data")
	("n-threads,j", po::value<int>(&params.n_threads), "Number of threads")
	("ignore-hypo", po::bool_switch(&params.ignore_hypo), "Ignore hypothetical protein kmers when making calls")
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
    if (params.input_files.size() == 0)
    {
	std::cout << desc << "\n";
	exit(1);
    }
}

int main(int argc, char **argv)
{
    program_parameters params;
    process_options(argc, argv, params);

    std::cerr << "Data size " << sizeof(StoredKmerData) << "\n";

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
    caller.ignore_hypothetical(params.ignore_hypo);

    using cbf = std::function<void(const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd)>;

    cbf hit_cb;
    if (params.debug_hits)
    {
	hit_cb = [&caller](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
	    std::cout << kmer << "\t" << offset << "\t" << caller.function_at_index(kd.function_index) << "\t" << kd.median << "\t" << kd.mean << "\t" << kd.var << "\t" << sqrt(kd.var) << "\t" << "\n";
	};
    }
    else
    {
	hit_cb = [](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {};
    }
    
/*
    auto hit_cb = [](const std::string &id, const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
    };
    auto hit_cb = [&caller, &params](const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
	if (params.debug_hits)
	{
	    std::cout << kmer << "\t" << caller.function_at_index(kd.function_index) << "\t" << kd.median << "\t" << kd.mean << "\t" << kd.var << "\t" << sqrt(kd.var) << "\t" << "\n";
	}
    };
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
		break;
	    anno_out << buf;
	}
    });

//    auto call_cb = [&caller, &anno_out](const std::string &id, const std::string &func, FunctionIndex func_index, float score) {
	
//	anno_out << id << "\t" << func << "\t" << func_index << "\t" << score << "\n";
//    };

    tbb::concurrent_vector<fs::path> ivec(params.input_files.begin(), params.input_files.end());
    tbb::parallel_for(ivec.range(), [&caller, &output_queue, &db_base, &params, &hit_cb](auto inp)
    {
	for (auto input_path: inp)
	{
	    fs::ifstream ifstr(input_path);
	    shared_buf_t buf = std::make_shared<boost::asio::streambuf>();
	    
	    std::ostream bufstr(buf.get());
	    
	    auto call2_cb = [&bufstr](const std::string &id, const std::string &func, FunctionIndex func_index, float score, size_t seq_size)
		{
		    bufstr << id << "\t" << func << "\t" << func_index << "\t" << score << "\n";
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
}

