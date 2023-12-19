#include "nudb_kmer_db.h"
#include "cmph_kmer.h"
#include "call_functions.h"
#include "fasta_parser.h"
#include "matrix_distance.h"

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

  @mainpage kmers-matrix-distance-folder

  Compute the all-to-all distance matrix for each fasta in a folder of fasta files. Writes output
  to a file in the output folder with the same name as the file in the input folder.

  For each sequence, we compute the set of signature kmers.
  Then for each pair of sequences, we compute the size of the intersection of these sets.

*/

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct program_parameters
{
    fs::path data_dir;
    fs::path input_dir;
    fs::path output_dir;
    bool debug_hits = false;
    bool verbose = false;
    int n_threads = 1;
};

void process_options(int argc, char **argv, program_parameters &params)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " data-dir input-dir output-dir\nAllowed options";

    po::options_description desc(x.str());
    desc.add_options()
	("data-dir,d", po::value<fs::path>(&params.data_dir), "Data directory")
	("input-dir", po::value<fs::path>(&params.input_dir), "Input directory")
	("output-dir", po::value<fs::path>(&params.output_dir), "Output directory")
	("n-threads,j", po::value<int>(&params.n_threads), "Number of threads")
	("j", po::value<int>(&params.n_threads), "Number of threads")
	("debug-hits", po::bool_switch(&params.debug_hits), "Debug kmer hits")
	("verbose", po::bool_switch(&params.verbose), "Enable verbose mode")
	("help,h", "show this help message");

    po::positional_options_description pos;
    pos.add("data-dir", 1)
	.add("input-dir", 1)
	.add("output-dir", 1);
    
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

//    NuDBKmerDb<StoredKmerData, 8> nudb(db_base);
    DbType nudb(db_base);

    if (!nudb.exists())
    {
	std::cerr << "Database " << db_base << " does not exist\n";
	exit(1);
    }
    nudb.open();
    FunctionCaller<DbType> caller(nudb, params.data_dir / "function.index");

    tbb::concurrent_vector<std::pair<fs::path, fs::path>> work;
    for (auto dit: fs::directory_iterator(params.input_dir))
    {
	if (fs::is_regular_file(dit.path()))
	{
	    fs::path output = params.output_dir / dit.path().filename();
	    if (!fs::exists(output))
		work.emplace_back(dit.path(), output);
	}
    }

    for (auto p: work)
    {
	std::cerr << p.first << " " << p.second << "\n";
    }

    tbb::parallel_for(work.range(), [&caller, &params](auto r) {
	for (auto went: r)
	{
	    const fs::path &input = went.first;
	    const fs::path &output = went.second;

	    MatrixDistance<FunctionCaller<DbType>> md(caller, input, output, params.verbose);
	    md.compute();
	}
    });
}

