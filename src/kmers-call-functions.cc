#include "nudb_kmer_db.h"
#include "call_functions.h"
#include "fasta_parser.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>

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
    bool debug_hits = false;
};

void process_options(int argc, char **argv, program_parameters &params)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " data-dir input-file [input-file, ...]\nAllowed options";

    po::options_description desc(x.str());
    desc.add_options()
	("data-dir,d", po::value<fs::path>(&params.data_dir), "Data directory")
	("input-files,i", po::value<std::vector<fs::path>>(&params.input_files), "Input files")
	("output-files,o", po::value<fs::path>(&params.output_file), "Output file")
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

    auto db_base = params.data_dir / "kmer_data";
    NuDBKmerDb<StoredKmerData, 8> nudb(db_base);

    if (!nudb.exists())
    {
	std::cerr << "Database " << db_base << " does not exist\n";
	exit(1);
    }
    nudb.open();

    FunctionCaller<NuDBKmerDb<StoredKmerData, 8>> caller(nudb, params.data_dir / "function.index");

    auto hit_cb = [&caller, &params](const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &kd) {
	if (params.debug_hits)
	{
	    std::cout << kmer << "\t" << caller.function_at_index(kd.function_index) << "\t" << kd.median << "\t" << kd.mean << "\t" << kd.var << "\t" << sqrt(kd.var) << "\t" << "\n";
	}
    };


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
	
    auto call_cb = [&caller, &anno_out](const std::string &id, const std::string &func, FunctionIndex func_index, float score) {
	anno_out << id << "\t" << func << "\t" << func_index << "\t" << score << "\n";
    };

    for (auto input_path:  params.input_files)
    {
	fs::ifstream ifstr(input_path);

	std::cerr << input_path << "\n";
	caller.process_fasta_stream(ifstr, hit_cb, call_cb);

	ifstr.close();
    }
	
}

