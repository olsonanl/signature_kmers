#include <iostream>
#include <fstream>
#include <istream>
#include <ostream>
#include <string>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options.hpp>
#include <chrono>

#include "kmer_request_server.h"

#include <unistd.h>

using namespace boost::filesystem;
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
    std::ostringstream x;
    x << "Usage: " << argv[0] << " [options] listen-port kmer-data-dir\nAllowed options";
    po::options_description desc(x.str());

    std::string listen_port;
    std::string listen_port_file;
    std::string kmer_data;
    std::string peg_kmer_data;
    int n_load_threads = 1;
    int n_kmer_threads = 1;
    bool daemonize = false;
    std::string pid_file;
    std::string kmer_family_distribution_file;
    bool no_listen = false;
    std::string family_reps;
    
    desc.add_options()
	("help,h", "show this help message")
	("n-family-file-threads", po::value<int>()->default_value(1), "number of family file reader threads")
	("n-inserter-threads", po::value<int>()->default_value(1), "number of kmer inserter threads")
	("n-load-threads", po::value<int>(&n_load_threads)->default_value(1), "number of NR load threads")
	("n-kmer-threads", po::value<int>(&n_kmer_threads)->default_value(1), "number of kmer processing threads")
	("listen-port-file", po::value<std::string>(&listen_port_file)->default_value("/dev/null"), "save the listen port to this file")
	("peg-kmer-data", po::value<std::string>(&peg_kmer_data), "precomputed PEG/kmer data file")
	("listen-port,l", po::value<std::string>(&listen_port)->required(), "port to listen on. 0 means to choose a random port")
	("kmer-data-dir,d", po::value<std::string>(&kmer_data)->required(), "kmer data directory")
	("kmer-version", po::value<std::string>(), "kmer data version string")
	("families-genus-mapping", po::value<std::string>(), "genus name to taxid mapping file")
	("families-file", po::value<std::string>(), "families file")
	("families-nr", po::value<std::vector<std::string>>()->multitoken(), "families NR data")
	("families-version", po::value<std::string>(), "families data version string")
	("family-reps", po::value<std::string>(&family_reps), "family representative pegs")
	("kmer-family-distribution-file", po::value<std::string>(&kmer_family_distribution_file), "kmer family distribution logfile")
	("reserve-mapping", po::value<unsigned long>(), "Reserve this much space in global mapping table")
	("no-populate-mmap", po::bool_switch(), "Don't populate mmap data at startup")
	("debug-http", po::bool_switch(), "Debug HTTP protocol")
	("no-listen", po::bool_switch(&no_listen), "Don't listen - just load data and quit. For profiling.")
	("daemonize", po::bool_switch(&daemonize), "Run the service in the background")
	("pid-file", po::value<std::string>(&pid_file), "Write the process id to this file")
	;
    po::positional_options_description pd;
    pd.add("listen-port", 1)
	.add("kmer-data-dir", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
	      options(desc).positional(pd).run(), vm);

    if (vm.count("help"))
    {
	std::cout << desc << "\n" <<
	    "If the kmer data directory contains files families.dat and a\n" <<
	    "directory families.nr it will be assumed that these files contain\n" <<
	    "family data files and will be loaded at startup. A file VERSION\n" <<
	    "will set the family version data to the contents of the first line of that file\n";
	return 1;
    }

    try {
	po::notify(vm);
    } catch (po::required_option &e)
    {
	std::cerr << "Invalid command line: " << e.what() << "\n";
	std::cerr << desc << "\n";
	exit(1);
    }

    /*
     * Check the kmer data directory for family and version files. If present,
     * modify the parameters accordingly.
     */
    {
	path fams_file(kmer_data);
	fams_file /= "families.dat";
	
	path fams_nr(kmer_data);
	fams_nr /= "families.nr";

	path version(kmer_data);
	version /= "VERSION";

	path fams_version(kmer_data);
	fams_version /= "families.version";

	path fams_genus_map(kmer_data);
	fams_genus_map /= "families.genus_map";

	char *extra_options[500];
	int extra_count = 0;
	extra_options[extra_count++] = strdup(argv[0]);

	std::string kversion("unknown"), fversion("unknown");

	if (is_regular_file(version))
	{
	    ifstream v(version);
	    std::getline(v, kversion);
	    fversion = kversion;
	}

	if (is_regular_file(fams_version))
	{
	    ifstream v(fams_version);
	    std::getline(v, fversion);
	}

	extra_options[extra_count++] = strdup("--kmer-version");
	extra_options[extra_count++] = strdup(kversion.c_str());
	
	extra_options[extra_count++] = strdup("--families-version");
	extra_options[extra_count++] = strdup(fversion.c_str());

	if (is_regular_file(fams_genus_map))
	{
	    extra_options[extra_count++] = strdup("--families-genus-mapping");
	    extra_options[extra_count++] = strdup(fams_genus_map.c_str());
	}
	
	if (is_regular_file(fams_file))
	{
	    extra_options[extra_count++] = strdup("--families-file");
	    extra_options[extra_count++] = strdup(fams_file.c_str());
	}
	if (is_directory(fams_nr))
	{
	    bool printed_arg = 0;
	    for (auto dit : directory_iterator(fams_nr))
	    {
		if (is_regular_file(dit.path()))
		{
		    if (!printed_arg)
		    {
			printed_arg = true;
			extra_options[extra_count++] = strdup("--families-nr");
		    }
		    extra_options[extra_count++] = strdup(dit.path().string().c_str());
		}
	    }
	}
	for (int i = 0; i < extra_count; i++)
	    std::cerr << i << ": " << extra_options[i] << "\n";

	po::store(po::command_line_parser(extra_count, extra_options).
		  options(desc).positional(pd).run(), vm);

	for (int i = 0; i < extra_count; i++)
	    free(extra_options[i]);
    }

    if (daemonize)
    {
	pid_t child = fork();
	if (child < 0)
	{
	    std::cerr << "fork failed: " << strerror(errno) << "\n";
	    exit(1);
	}
	if (child > 0)
	{
	    if (!pid_file.empty())
	    {
		std::ofstream pf(pid_file);
		pf << child << "\n";
		pf.close();
	    }
	    exit(0);
	}
	pid_t sid = setsid();
	if (sid < 0)
	{
	    std::cerr << "setsid failed: " << strerror(errno) << "\n";
	    exit(1);
	}
    }
    else if (!pid_file.empty())
    {
	std::ofstream pf(pid_file);
	pf << getpid() << "\n";
	pf.close();
    }

    std::ofstream kmer_family_distribution_stream;
    if (!kmer_family_distribution_file.empty())
    {
	kmer_family_distribution_stream.open(kmer_family_distribution_file);
    }


    //
    // Load family reps if in use
    //

    std::shared_ptr<FamilyReps> reps = 0;
    if (family_reps.size())
    {
	reps = std::make_shared<FamilyReps>();
	path reps_path(family_reps);
	if (is_directory(reps_path))
	{
	    reps->load_reps_directory(reps_path);
	}
	else if (is_regular_file(reps_path))
	{
	    reps->load_reps_file(reps_path);
	}
	else
	{
	    std::cerr << "Invalid family reps file " << reps_path << "\n";
	    exit(1);
	}
	
    }

    boost::asio::io_service io_service;

    std::shared_ptr<ThreadPool> tp = std::make_shared<ThreadPool>(kmer_data);

    tp->start(n_load_threads);

    bool family_mode = g_parameters->count("families-file");

    std::chrono::time_point<std::chrono::high_resolution_clock> start, end, end2;
    start = std::chrono::high_resolution_clock::now();
    std::shared_ptr<KmerRequestServer> kserver = std::make_shared<KmerRequestServer>(io_service, listen_port, listen_port_file,
										     tp, family_mode);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cerr << "family load time " << elapsed_seconds.count() << "\n";

    if (reps)
	kserver->set_family_reps(reps);

    if (!kmer_family_distribution_file.empty())
    {
	std::cerr << "write distribution to " << kmer_family_distribution_file << "\n";
	kserver->root_mapping()->write_kmer_distribution(kmer_family_distribution_stream);
	kmer_family_distribution_stream.close();

	end2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = end2 - end;
	std::cerr << "distribution write time " << elapsed_seconds.count() << "\n";
    }

    if (no_listen)
    {
	std::cerr << "Quitting due to --no-listen being set\n";
	exit(0);
    }

    int additional_threads = n_kmer_threads - n_load_threads;
    if (additional_threads > 0)
    {
	std::cerr << "Starting " << additional_threads << " additional threads\n";
	tp->add_threads(additional_threads);
    }
    kserver->startup();

    #ifdef GPROFILER
    std::cerr << "profiler enable\n";
    ProfilerStart("prof.out");
    #endif
    io_service.run();

    tp->stop();

    #ifdef GPROFILER
    ProfilerStop();
    std::cerr << "profiler disable\n";
    #endif

    return 0;
}
