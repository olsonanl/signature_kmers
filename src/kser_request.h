#ifndef _KREQUEST2_H
#define _KREQUEST2_H

/*
 * A kmer request wraps up the code for parsing out
 * the incoming HTTP request from a client and handling
 * the execution of the request itself.
 */

#include <vector>
#include <string>
#include <set>
#include <memory>
#include <experimental/optional>
#include <boost/asio.hpp>
#include <boost/timer/timer.hpp>
#include <ostream>
#include "kmer.h"
#include "klookup.h"
#include "klookup2.h"
#include "klookup3.h"
#include "compute_request.h"
#include "threadpool.h"

class KmerRequestServer;

class KmerRequest2 : public std::enable_shared_from_this<KmerRequest2>
{
public:
    KmerRequest2(std::shared_ptr<KmerRequestServer> server,
		 boost::asio::io_service &io_service,
		 std::shared_ptr<std::map<std::string, std::shared_ptr<KmerPegMapping>>> mapping_map,
		 std::shared_ptr<ThreadPool> thread_pool
	);
    ~KmerRequest2();

    void do_read();

    boost::asio::io_service &io_service()
    {
	return io_service_;
    }

    std::shared_ptr<ThreadPool> thread_pool()
    {
	return thread_pool_;
    }

    boost::asio::ip::tcp::socket& socket()
    {
	return socket_;
    }

    boost::asio::streambuf &request()
    {
	return request_;
    }

    std::map<std::string, std::string> &parameters()
    {
	return parameters_;
    }

    void write_header(std::ostream &os, int code, const std::string &status);

    void exit_request();

    std::shared_ptr<KmerRequestServer> server() { return server_; }

    std::experimental::optional<std::string> header(const std::string &c) const {
	auto x = headers_.find(c);
	if (x == headers_.end())
	    return {};
	else
	    return x->second;
    }

    void respond(int code, const std::string &status, const std::string &result, std::function<void()> on_done);

private:

    void read_initial_line(boost::system::error_code err, size_t bytes);
    void read_headers(boost::system::error_code err, size_t bytes);

    void process_request();

    std::shared_ptr<KmerRequestServer> server_;

    std::string request_type_;
    std::string path_;
    std::string parameters_raw_;
    std::string fragment_;
    std::string http_version_;
    std::map<std::string, std::string> parameters_;
    std::map<std::string, std::string> headers_;

    boost::asio::io_service &io_service_;
    boost::asio::ip::tcp::socket socket_;
    boost::asio::streambuf request_;
    std::shared_ptr<std::map<std::string, std::shared_ptr<KmerPegMapping>>> mapping_map_;

    boost::asio::streambuf response_;

    std::shared_ptr<ComputeRequest> active_request_;
    std::shared_ptr<ThreadPool> thread_pool_;
};


#endif
