#ifndef _KSERVER_H
#define _KSERVER_H


/*
 * Kmer-request server.
 *
 * We speak pidgin HTTP here so we can make this available over a proxy.
 * This is NOT a general purpose HTTP server.
 *
 */

#include <boost/asio.hpp>
#include <set>
#include <memory>
#include "kmer.h"
#include "kmer_inserter.h"
#include "krequest2.h"
#include "threadpool.h"
#include "tbb/atomic.h"
#include "family_reps.h"

#include "nr_loader.h"

class KmerRequestServer : public std::enable_shared_from_this<KmerRequestServer>
{
public:
    KmerRequestServer(boost::asio::io_service& io_service,
		      const std::string &port,
		      const std::string &port_file,
		      std::shared_ptr<ThreadPool> thread_pool,
		      bool family_mode = false);

    void load_families_nr(const std::string &file);
    void startup();
    void deactivate(std::shared_ptr<KmerRequest2> x);

    bool family_mode() { return family_mode_; }

    std::shared_ptr<KmerPegMapping> root_mapping() { return root_mapping_; }
    std::shared_ptr<FamilyReps> family_reps() { return family_reps_; }
    void set_family_reps(std::shared_ptr<FamilyReps> reps) { family_reps_ = reps; }

private:
    bool family_mode_;

    void do_accept2();
    void on_accept2(boost::system::error_code ec, std::shared_ptr<KmerRequest2>);
    
    void do_await_stop();
    std::shared_ptr<ThreadPool> thread_pool_;
    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::acceptor acceptor_;
    boost::asio::signal_set signals_;
    std::string port_;
    std::string port_file_;
    std::set<std::shared_ptr<KmerRequest2> > active_;

    std::shared_ptr<std::map<std::string, std::shared_ptr<KmerPegMapping>>> mapping_map_;
    std::shared_ptr<KmerPegMapping> root_mapping_;

    std::shared_ptr<FamilyReps> family_reps_;

private:
    NRLoadState load_state_;

};


#endif
