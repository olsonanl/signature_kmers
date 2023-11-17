
#include "kserver.h"

#include <boost/asio.hpp>
#include <boost/filesystem.hpp>
#include <boost/bind.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <memory>
#include "global.h"
#include "kguts.h"

KmerRequestServer::KmerRequestServer(boost::asio::io_service& io_service,
				     const std::string &port,
				     const std::string &port_file,
				     std::shared_ptr<ThreadPool> thread_pool,
				     bool family_mode
    ) :
    family_mode_(family_mode),
    thread_pool_(thread_pool),
    io_service_(io_service),
    acceptor_(io_service_),
    signals_(io_service_),
    port_(port),
    port_file_(port_file),
    family_reps_(0),
    load_state_("master")
    
{
    mapping_map_ = std::make_shared<std::map<std::string, std::shared_ptr<KmerPegMapping> > >();

    auto x = mapping_map_->insert(std::make_pair("", std::make_shared<KmerPegMapping>()));
    // std::cerr << "insert created new: " << x.second << " key='" << x.first->first << "'\n";
    root_mapping_ = x.first->second;
    
    root_mapping_->load_genome_map((*g_parameters)["kmer-data-dir"].as<std::string>() + "/genomes");

    if (g_parameters->count("families-genus-mapping"))
    {
	std::string mapfile = (*g_parameters)["families-genus-mapping"].as<std::string>();
	root_mapping_->load_genus_map(mapfile);
    }

    /*
     * If we are preloading a families file, start that off in the background
     * using the thread pool.
     *
     * If we are in family mode we cannot load in a thread as we need the family
     * ID assignments for loading the NR that we're likely also prelaoding.
     *
     * Legacy test: we have now defined family mode as what we're operating in
     * when the families file has been specified.
     */
    
    if (g_parameters->count("families-file"))
    {
	std::string ff = (*g_parameters)["families-file"].as<std::string>();

	if (family_mode_)
	{
	    std::cerr << "Loading (immediate) families from " << ff << "...\n";
	    root_mapping_->load_families(ff);
	    std::cerr << "Loading families from " << ff << "... done\n";
	}
	else
	{
	    load_state_.pending_inc();
	    thread_pool_->post([this, ff]() {
		    std::cerr << "Loading families from " << ff << "...\n";
		    root_mapping_->load_families(ff);
		    std::cerr << "Loading families from " << ff << "... done\n";
		    load_state_.pending_dec();
		});
	}
    }

    if (g_parameters->count("reserve-mapping"))
    {
	unsigned long count = (*g_parameters)["reserve-mapping"].as<unsigned long>();
	std::cerr << "Reserving " << count << " bytes in mapping table\n";
	root_mapping_->reserve_mapping_space(count);
    }

    int n_inserters = 1;
    if (g_parameters->count("n-inserter-threads"))
    {
	n_inserters = (*g_parameters)["n-inserter-threads"].as<int>();
    }

    KmerInserter inserter(n_inserters, root_mapping_);
    if (family_mode)
    {
	inserter.start();
    }
    std::vector<NRLoader *> active_loaders;

    if (g_parameters->count("families-nr"))
    {
	auto files = (*g_parameters)["families-nr"].as<std::vector<std::string> >();
	for (auto file: files)
	{
	    load_state_.pending_inc();
	    std::cerr << "Queue load NR file " << file << "\n";
	    NRLoader *loader = new NRLoader(load_state_, file, root_mapping_, thread_pool_, files.size(), family_mode_, inserter);
	    active_loaders.push_back(loader);
	    loader->start();
	}
    }
	
    std::cerr << "wait for threads to finish\n";

    load_state_.pending_wait();
    
    std::cerr << "done waiting\n";

    for (auto v: active_loaders)
    {
	// std::cerr << "remove loader for " << v->file_ << "\n";
	delete v;
    }

    /*
     * When we are done, we need to clear the inserters.
     */
    inserter.stop();
}

/*
 * Need to split startup code from constructor due to use of enable_shared_from_this().
 */
void KmerRequestServer::startup()
{

    /*
     * Set up for clean signal handling / termination
     */
    signals_.add(SIGINT);
    signals_.add(SIGTERM);
    signals_.add(SIGQUIT);
    do_await_stop();

    /*
     * Set up listener
     */

    boost::asio::ip::tcp::resolver resolver(io_service_);
    boost::asio::ip::tcp::endpoint endpoint = *resolver.resolve({"0.0.0.0", port_});
    acceptor_.open(endpoint.protocol());
    acceptor_.set_option(boost::asio::ip::tcp::acceptor::reuse_address(true));
    acceptor_.bind(endpoint);
    acceptor_.listen();
    std::cout << "Listening on " << acceptor_.local_endpoint() << "\n";
    if (!port_file_.empty())
    {
	std::ofstream out(port_file_);
	out << acceptor_.local_endpoint().port() << "\n";
	out.close();
    }
	    
    do_accept2();
}

void KmerRequestServer::do_accept2()
{
    std::shared_ptr<KmerRequest2> r = std::make_shared<KmerRequest2>(shared_from_this(), io_service_, mapping_map_, thread_pool_);
    // std::cerr << "create " << r << "\n";
    acceptor_.async_accept(r->socket(),
			   boost::bind(&KmerRequestServer::on_accept2, this,
				       boost::asio::placeholders::error, r));
    //std::cerr << "leaving do_accept2 r use count=" << r.use_count() << "\n";
}

void KmerRequestServer::on_accept2(boost::system::error_code ec, std::shared_ptr<KmerRequest2> r)
{
    // std::cerr << "on-accept2 use=" << r.use_count() << "\n";
    g_timer.start();
    // Check whether the server was stopped by a signal before this
    // completion handler had a chance to run.
    if (!acceptor_.is_open())
    {
	std::cout << "not open\n";
	return;
    }
    
    if (!ec)
    {
	/*
	 * Connection has come in.
	 * Begin parsing the request line and headers.
	 */
	
	active_.insert(r);
	// std::cerr << "start read " << r << "\n";
	r->do_read();
    }
    
    do_accept2();
}

void KmerRequestServer::deactivate(std::shared_ptr<KmerRequest2> x)
{
    active_.erase(x);
}

void KmerRequestServer::do_await_stop()
{
    signals_.async_wait(
	[this](boost::system::error_code ec, int signo)
	{
	    std::cout << "Exiting with signal " << signo << "\n";
	    acceptor_.close();
	});
}


