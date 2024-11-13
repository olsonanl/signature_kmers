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

  @mainpage kmers-matrix-distance-merge

  Given a base directory (typically the genus.data directory in a family build)
  and an optional list of family ids, compute the kmer distances for the families.
  The data will come from fam-dir/genus-name/fasta_by_function/<family-id>

  If no family id is specified, use the indices from the function.index file.

  For each sequence, we compute the set of signature kmers.
  Then for each pair of sequences, we compute the size of the intersection of these sets.

*/

namespace po = boost::program_options;
namespace fs = boost::filesystem;

struct program_parameters
{
    fs::path data_dir;
    fs::path base_dir;
    std::vector<std::string> family_ids;
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
	("base-dir", po::value<fs::path>(&params.base_dir), "Bae directory")
	("output-dir", po::value<fs::path>(&params.output_dir), "Output directory")
	("family-ids", po::value<std::vector<std::string>>(&params.family_ids), "Family ids")
	("n-threads,j", po::value<int>(&params.n_threads), "Number of threads")
	("debug-hits", po::bool_switch(&params.debug_hits), "Debug kmer hits")
	("verbose", po::bool_switch(&params.verbose), "Enable verbose mode")
	("help,h", "show this help message");

    po::positional_options_description pos;
    pos.add("data-dir", 1)
	.add("base-dir", 1)
	.add("output-dir", 1)
	.add("family-ids", -1)
	;
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

    if (!fs::is_directory(params.base_dir))
    {
	std::cerr << "Base directory " << params.base_dir << " is not a valid directory\n";
	exit(1);
    }

    using DbType = CmphKmerDb<StoredKmerData, 8>;

    DbType nudb(db_base);

    if (!nudb.exists())
    {
	std::cerr << "Database " << db_base << " does not exist\n";
	exit(1);
    }
    nudb.open();
    FunctionCaller<DbType> caller(nudb, params.data_dir / "function.index");


    /*
     * If we don't have a family-ids list, use the indices from the caller.
     * In either event we initialize into a concurrent_vector.
     */

    tbb::concurrent_vector<std::string> work;

    if (params.family_ids.empty())
    {
	auto funcs = caller.function_index();
	work.resize(funcs.size());
	for (int i = 0; i < funcs.size(); i++)
	    work[i] = std::to_string(i);
    }
    else
    {
	std::copy(params.family_ids.begin(), params.family_ids.end(), std::back_inserter(work));
    }
    
    std::vector<fs::path> genus_dirs;
    
    for (auto dit: fs::directory_iterator(params.base_dir))
    {
	fs::path gpath = dit.path();
	if (fs::is_directory(gpath))
	{
	    if (fs::is_regular_file(gpath / "local.family.defs"))
	    {
		genus_dirs.push_back(gpath);
	    }
	}
    }
    if (genus_dirs.empty())
    {
	std::cerr << "No valid genus directories found in " << params.base_dir << "\n";
	exit(1);
    }

    tbb::parallel_for(work.range(), [&caller, &params, &genus_dirs](auto r) {
	for (auto fam: r)
	{
	    std::vector<fs::path> inputs;
	    for (auto g:  genus_dirs)
		inputs.push_back(g / "fasta_by_function" / fam);
	    fs::path output = params.output_dir / fam;
	    MatrixDistance<FunctionCaller<DbType>> md(caller, inputs, output, params.verbose);
	    md.compute();
	}
    });
}

