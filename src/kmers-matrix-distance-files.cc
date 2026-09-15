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
#include <boost/format.hpp>

#include <thread>

#include <stdexcept>
#include <vector>

/*!

  @mainpage kmers-matrix-distance-folder

  Compute the all-to-all distance matrix for one or more fasta files.

  For each sequence, we compute the set of signature kmers.
  Then for each pair of sequences, we compute the size of the intersection of these sets.

*/

namespace po = boost::program_options;
namespace fs = boost::filesystem;

/*
 * Ugh. cf https://github.com/boostorg/program_options/issues/69
 * It's possible this isn't necessary if we use std::filesystem instead tho.
 */

namespace boost
{
    template <>
    inline fs::path lexical_cast<fs::path, std::basic_string<char>>(const std::basic_string<char> &arg)
    {
	return fs::path(arg);
    }
}

struct program_parameters
{
    fs::path data_dir;
    std::vector<fs::path> input_files;
    fs::path output_file;
    bool debug_hits = false;
    bool verbose = false;
    int n_threads = 1;
    int min_kmers_in_common = 1;
    int min_protein_len = 0;
    bool write_counts = false;
    fs::path kmer_stats;
};

void process_options(int argc, char **argv, program_parameters &params)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " data-dir output-file input-file [input-file...]\nAllowed options";

    fs::path input_file_list;
    
    po::options_description desc(x.str());
    desc.add_options()
	("min-kmers-in-common", po::value<int>(&params.min_kmers_in_common),
	 "Minimum number of signature kmers in common required to report a match")
	("min-protein-len", po::value<int>(&params.min_protein_len),
	 "Minimum size of a protein to include in the computation")
	("data-dir,d", po::value<fs::path>(&params.data_dir), "Data directory")
	("input-file-list", po::value<fs::path>(&input_file_list), "File containing list of input files")
	("input-file", po::value<std::vector<fs::path>>(&params.input_files), "Input file(s)")
	("output-file", po::value<fs::path>(&params.output_file), "Output file")
	("n-threads,j", po::value<int>(&params.n_threads), "Number of threads")
	("j", po::value<int>(&params.n_threads), "Number of threads")
	("debug-hits", po::bool_switch(&params.debug_hits), "Debug kmer hits")
	("write-counts", po::bool_switch(&params.debug_hits), "Write kmer counts. False if piping to MCL")
	("kmer-stats", po::value<fs::path>(&params.kmer_stats), "Write kmer stats here")
	("verbose", po::bool_switch(&params.verbose), "Enable verbose mode")
	("help,h", "show this help message");

    po::positional_options_description pos;
    pos.add("data-dir", 1)
	.add("output-file", 1)
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

    if (!input_file_list.empty())
    {
	fs::ifstream ifstr(input_file_list);
	std::string line;
	while (std::getline(ifstr, line, '\n'))
	{
	    fs::path path{line};
	    if (!fs::is_regular_file(fs::status(path)))
	    {
		throw std::runtime_error(str(boost::format("input file %1% does not exist") % path));
	    }
	    params.input_files.push_back(path);
	}
    }
}

int main(int argc, char **argv)
{
    program_parameters params;
    process_options(argc, argv, params);

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, params.n_threads);

    auto db_base = params.data_dir / "kmer_data";

    using DbType = CmphKmerDb<StoredKmerData, 8>;

    DbType nudb(db_base);

    for (auto x: params.input_files)
    {
	std::cerr << "input: " << x <<"\n";
    }
    std::cerr << "output: " << params.output_file << "\n";
       
    if (!nudb.exists())
    {
	std::cerr << "Database " << db_base << " does not exist\n";
	exit(1);
    }
    nudb.open();
    FunctionCaller<DbType> caller(nudb, params.data_dir / "function.index");

    MatrixDistance<FunctionCaller<DbType>> md(caller, params.input_files, params.output_file, params.verbose,
					      params.min_kmers_in_common, params.min_protein_len);
    size_t n_distinct_kmers = md.compute(params.write_counts);
    if (!params.kmer_stats.empty())
    {
	std::ofstream of(params.kmer_stats);
	of << n_distinct_kmers;
    }
    
}

