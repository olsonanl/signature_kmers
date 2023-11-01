#include "signature_build.h"
#include "path_utils.h"
#include "kept_kmer_db.h"
#include "call_functions.h"
#include "nudb_kmer_db.h"

#include <boost/program_options.hpp>

#include <tbb/global_control.h>
#include <tbb/concurrent_map.h>

namespace po = boost::program_options;
namespace fs = boost::filesystem;

const int K = 8;
const int MaxSequencesPerFile = 100000;

static bool process_command_line_options(int argc, char *argv[],
				  std::vector<fs::path> &function_definitions,
				  std::vector<fs::path> &fasta_data,
				  std::vector<fs::path> &fasta_data_kept_functions,
				  std::vector<std::string> &good_functions,
				  std::vector<std::string> &good_roles,
				  fs::path &deleted_fids_file,
				  int &min_reps_required,
				  fs::path &kmer_data_dir,
				  fs::path &final_kmers,
				  std::string &nudb_file,
				  int &n_threads)
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options]\nAllowed options";
    po::options_description desc(x.str());

    std::vector<std::string> definition_dirs;
    std::vector<std::string> fasta_dirs;
    std::vector<std::string> fasta_keep_dirs;
    std::vector<std::string> good_function_files;
    std::vector<std::string> good_role_files;

    n_threads = 1;

    desc.add_options()
	("definition-dir,D", po::value<std::vector<std::string>>(&definition_dirs)->multitoken(), "Directory of function definition files")
	("fasta-dir,F", po::value<std::vector<std::string>>(&fasta_dirs)->multitoken(), "Directory of fasta files of protein data")
	("fasta-keep-functions-dir,K", po::value<std::vector<std::string>>(&fasta_keep_dirs), "Directory of fasta files of protein data (keep functions defined here)")
	("good-functions", po::value<std::vector<std::string>>(&good_function_files), "File containing list of functions to be kept")
	("good-roles", po::value<std::vector<std::string>>(&good_role_files), "File containing list of roles to be kept")
	("deleted-features-file", po::value<fs::path>(&deleted_fids_file), "File containing list of deleted feature IDs")
	("kmer-data-dir", po::value<fs::path>(&kmer_data_dir), "Write kmer data files to this directory")
	("nudb-file", po::value<std::string>(&nudb_file), "Write saved kmers to this NuDB file base. Should be on a SSD drive.")
	("min-reps-required", po::value<int>(&min_reps_required), "Minimum number of genomes a function must be seen in to be considered for kmers")
	("final-kmers", po::value<fs::path>(&final_kmers), "Write final.kmers file to be consistent with km_build_Data")
	("n-threads", po::value<int>(&n_threads), "Number of threads to use")
	("help,h", "show this help message");

    po::variables_map vm;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n";
	return false;
    }

    po::notify(vm);

    /*
     * Read definition and fasta dirs to populate the path lists.
     */
    populate_path_list(definition_dirs, function_definitions);
    populate_path_list(fasta_dirs, fasta_data);
    populate_path_list(fasta_keep_dirs, fasta_data_kept_functions);

    std::cout << "definitions: ";
    for (auto x: definition_dirs)
	std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "fasta: ";
    for (auto x: fasta_dirs)
	std::cout << x << " ";
    std::cout << std::endl;
    std::cout << "keep: ";
    for (auto x: fasta_keep_dirs)
	std::cout << x << " ";
    std::cout << std::endl;

    load_strings(good_function_files, good_functions);
    load_strings(good_role_files, good_roles);

    return true;
}

void write_nudb_data(const std::string &nudb_file, const KeptKmers<8> &kmers)
{
    typedef NuDBKmerDb<StoredKmerData, 8> KDB;

    KDB db(nudb_file);

    if (!db.exists())
    {
	std::cerr << "creating new db\n";
	db.create();
    }
    db.open();
    
    for (auto ent: kmers)
    {
	auto k = ent.second;
	nudb::error_code ec;
	db.insert(k.kmer, k.stored_data, ec);
	if (ec)
	    std::cerr << "insert error: " << ec.message() << "\n";
    }
}


int main(int argc, char *argv[])
{
    std::vector<fs::path> function_definitions;
    std::vector<fs::path> fasta_data;
    std::vector<fs::path> fasta_data_kept_functions;

    std::vector<std::string> good_functions;
    std::vector<std::string> good_roles;

    fs::path final_kmers;
    fs::path deleted_fids_file;
    fs::path kmer_data_dir;

    int min_reps_required = 3;
    
    int n_threads;

    std::string nudb_file;

    if (!process_command_line_options(argc, argv,
				      function_definitions,
				      fasta_data,
				      fasta_data_kept_functions,
				      good_functions,
				      good_roles,
				      deleted_fids_file,
				      min_reps_required,
				      kmer_data_dir,
				      final_kmers,
				      nudb_file,
				      n_threads))
    {
	return 1;
    }

    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, n_threads);

    SignatureBuilder<K> builder(n_threads, MaxSequencesPerFile);

    builder.load_function_data(good_functions, good_roles, function_definitions);

    std::set<std::string> deleted_fids = load_set_from_file(deleted_fids_file);

    ensure_directory(kmer_data_dir);

    std::cerr << "load fasta\n";
    builder.load_fasta(fasta_data, false, deleted_fids);
    builder.load_fasta(fasta_data_kept_functions, true, deleted_fids);

    builder.process_kept_functions(min_reps_required, kmer_data_dir);

    if (!kmer_data_dir.empty())
    {
	fs::ofstream otu(kmer_data_dir / "otu.index");
	otu.close();
	fs::ofstream genomes(kmer_data_dir / "genomes");
	genomes << "empty genomes\n";
	genomes.close();
    }

    std::cerr << "extract kmers\n";
    builder.extract_kmers(deleted_fids); 
    std::cerr << "process kmers\n";
    builder.process_kmers();
    
    if (!final_kmers.empty())
    {
	std::cerr << "writing kmers to " << final_kmers << "\n";
	fs::ofstream kf(final_kmers);
	std::for_each(builder.kept_kmers().begin(), builder.kept_kmers().end(),
		      [&kf](const auto &k) {
			  const auto &v = k.second;
			  kf <<
			      k.first << "\t" <<
			      v.stored_data.avg_from_end << "\t" <<
			      v.stored_data.function_index << "\t" <<
			      "\n";
	    // kf << "\t" << k.seqs_containing_sig << "\t" << kmer_stats.seqs_with_func[k.function_index] << "\n";
	});
    }

    /*
     * Write kmer_stats.distinct_functions table ; this allows us
     * to reason about the coverage of functions in the database, as well as
     * which fusions have data for the component functions.
     */

    {
	fs::ofstream dfstr(kmer_data_dir / "distinct_functions");
	for (auto ent: builder.kmer_stats().distinct_functions)
	{
	    dfstr << ent.first << "\t" << builder.lookup_function(ent.first) << "\t" << ent.second << "\n";
	}
    }

    KeptKmerDB<K>  kdb(builder.kept_kmers());

    fs::path report_dir = kmer_data_dir / "recall.report.d";
    if (!fs::create_directory(report_dir))
    {
	std::cerr << "mkdir " << report_dir << " failed\n";
    }

    std::string fi_file = (kmer_data_dir / "function.index").string();

    /*
     * Begin recall of source data using newly created kmers.
     */
    
    FunctionCaller<KeptKmerDB<K>> kmer_caller(kdb, fi_file);

    struct call_data
    {
	std::string id;
	std::string old_func;
	std::string old_func_stripped;
	std::string new_func;
	int func_index;
	float score;
    };
    tbb::concurrent_map<std::string, call_data> recall_report;

    struct saver
    {
	const FunctionMap &fm;
	std::map<std::string, call_data> data;
	
	void operator()(const std::string &id, const std::string &func, int func_index, float score) {

	    std::string orig, orig_stripped;
	    fm.lookup_original_assignment(id, orig, orig_stripped);
	    if (orig_stripped != func)
	    {
		data.emplace(id, call_data { id, orig, orig_stripped, func, func_index, score});
		// std::cout << "CALL "  << id << "\t" << orig_call << "\t" << func << "\t" << func_index << "\t" << score << "\n";
	    }
	}
    };

    auto call_cb = [&builder, &recall_report](const std::string &id, const std::string &func, FunctionIndex func_index, float score) {

	std::string orig, orig_stripped;
	builder.function_map().lookup_original_assignment(id, orig, orig_stripped);
	if (orig_stripped != func)
	{
	    recall_report.emplace(id, call_data { id, orig, orig_stripped, func, func_index, score});
	    // std::cout << "CALL "  << id << "\t" << orig << "\t" << func << "\t" << func_index << "\t" << score << "\n";
	}
	
    };

    auto hit_cb = [&builder](const Kmer<8> &kmer, size_t offset, double seqlen, const StoredKmerData &k) {

	if (false)
	{
	    std::cout << kmer << "\t" <<
		builder.function_map().lookup_function(k.function_index) << "\t" <<
		k.median << "\t" <<
		k.mean << "\t" <<
		k.var << "\t" <<
		sqrt(k.var) << "\t" << "\n";
	}
    };

    std::cerr << "Begin recall\n";

    tbb::parallel_for(builder.all_fasta_data().range(), [&report_dir, &builder, &kmer_caller, &hit_cb, &call_cb](auto r) {
	for (auto file: r)
	{

	    fs::path outfile(report_dir / file.filename());

	    saver s { builder.function_map() } ;

	    fs::ifstream ifstr(file);

	    kmer_caller.process_fasta_stream(ifstr, hit_cb, s);

	    ifstr.close();

	    fs::ofstream ofstr(outfile);
	    for (auto ent: s.data)
	    {
		call_data &c = ent.second;
		ofstr << ent.first << "\t" << c.old_func << "\t" << c.old_func_stripped << "\t" << c.new_func << "\t" << c.func_index << "\t" << c.score << "\n";
	    }
	    ofstr.close();
	}
    });
    
    if (!nudb_file.empty())
    {
	std::cerr << "write nudb data " << nudb_file << "\n";
	write_nudb_data(nudb_file, builder.kept_kmers());
    }

    std::cerr << "all done\n";

    return 0;

}
