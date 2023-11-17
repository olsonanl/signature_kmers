#include "krequest2.h"
#include <iostream>
#include <sstream>
#include <boost/bind.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <boost/asio/buffer.hpp>
#include <boost/array.hpp>
#ifdef BLCR_SUPPORT
#include "libcr.h"
#endif

#include "klookup.h"
#include "kserver.h"
#include "kguts.h"
#include "global.h"
#include "add_request.h"
#include "matrix_request.h"
#include "query_request.h"
#include "lookup_request.h"
#include "fq_process_request.h"
#include "debug.h"

const boost::regex request_regex("^([A-Z]+) ([^?#]*)(\\?([^#]*))?(#(.*))? HTTP/(\\d+\\.\\d+)");
/* $1 = type
 * $2 = path
 * $4 = parameters
 * $6 = fragment
 */
const boost::regex mapping_path_regex("^/mapping/([^/]+)(/(add|matrix|lookup))$");
const boost::regex genus_lookup_path_regex("^/genus_lookup/([^/]+)$");

KmerRequest2::KmerRequest2(std::shared_ptr<KmerRequestServer> server,
			   boost::asio::io_service &io_service,
			   std::shared_ptr<std::map<std::string, std::shared_ptr<KmerPegMapping> > > mapping_map,
			   std::shared_ptr<ThreadPool> thread_pool) :
    server_(server),
    io_service_(io_service),
    socket_(io_service_),
    request_(1048576),
    mapping_map_(mapping_map),
    thread_pool_(thread_pool)
   
{
    // std::cerr << "construct krequest2 " << this << "\n";
}

KmerRequest2::~KmerRequest2()
{
//     std::cerr << "destroy krequest2 " << this << "\n";
}
			   
/*
 * Start processing our request.
 *
 * This was invoked when the acceptor completed the new connection on our socket_.
 *
 */
void KmerRequest2::do_read()
{
    // std::cerr << "do_read " << this << "\n";
    boost::asio::async_read_until(socket_, request_, "\n",
				  boost::bind(&KmerRequest2::read_initial_line, this,
					      boost::asio::placeholders::error,
					      boost::asio::placeholders::bytes_transferred));
}

			   
std::string make_string(boost::asio::streambuf& streambuf)
{
    return {buffers_begin(streambuf.data()),
	    buffers_end(streambuf.data())};
}

void stringPurifier ( std::string& s )
{
    for ( std::string::iterator it = s.begin(), itEnd = s.end(); it!=itEnd; ++it)
    {
	if ( static_cast<unsigned int>(*it) < 32 || static_cast<unsigned int>(*it) > 127 )
	{
	    (*it) = ' ';
	}
    }
}

void KmerRequest2::read_initial_line(boost::system::error_code err, size_t bytes)
{
    if (!err)
    {
	std::istream is(&request_);
	std::string line;
	if (std::getline(is, line))
	{
	    size_t s = line.find('\r');
	    if (s != std::string::npos)
	    {
		line.erase(s);
	    }
	    if (debug_is_set("debug-http"))
		std::cerr << "Request: " << line << "\n";
	    
	    boost::smatch match;
	    if (boost::regex_match(line, match, request_regex))
	    {
		request_type_ = match[1];
		path_ = match[2];
		parameters_raw_ = match[4];
		fragment_ = match[6];
		http_version_ = match[7];

		if (!parameters_raw_.empty())
		{
		    std::vector<std::string> parts;
		    boost::split(parts, parameters_raw_, boost::is_any_of(";&"));
		    for (auto it : parts)
		    {
			size_t pos = it.find('=');
			if (pos != std::string::npos)
			{
			    parameters_[it.substr(0, pos)] = it.substr(pos+1);
			}
		    }
		}

		// std::cerr << "Got " << request_type_ << " for path " << path_ << "\n";

		if (!is.eof())
		{
		    boost::asio::async_read_until(socket_, request_, "\n",
						  boost::bind(&KmerRequest2::read_headers, this,
							      boost::asio::placeholders::error,
							      boost::asio::placeholders::bytes_transferred));
		}
		else
		{
		    // std::cerr << "eof\n";
		    exit_request();
		}
	    }
	    else
	    {
		std::cerr << "Invalid request '" << line << "'\n";
		exit_request();
	    }
	}
	else
	{
	    std::cerr << "getline failed\n";
	    exit_request();
	}
    }
    else
    {
	std::cerr << "error " << err << "\n";
	exit_request();
    }
    
}
void KmerRequest2::read_headers(boost::system::error_code err, size_t bytes)
{
    // std::cerr << "read_headers " << this << " err=" << err << "\n";
    if (!err)
    {
	// std::string r = make_string(request_);
	// std::cerr << "Read buffer contains: "  << r << std::endl;

	std::istream is(&request_);
	std::string line;
	bool finished = false;
	while (!finished && std::getline(is, line))
	{
	    size_t s = line.find('\r');
	    if (s != std::string::npos)
	    {
		line.erase(s);
	    }
	    // std::cerr << "Got line '" << line << "'\n";

	    if (line != "")
	    {
		size_t x = line.find(':');
		std::string k(line.substr(0,x));
		x++;
		while (line[x] == ' ')
		    x++;
		std::string v(line.substr(x));
		std::transform(k.begin(), k.end(), k.begin(), ::tolower);
		headers_[k] = v;
	    }
	    else
	    {
		finished = true;
	    }
	}
	if (finished)
	{
	    if (debug_is_set("debug-http"))
	    {
		std::cerr << "Headers:\n";
		for (auto x: headers_)
		{
		    std::cerr << x.first << ": " << x.second << "\n";
		}
	    }
	    /*
	     * We don't support chunked encoding.
	     */
	    auto x = headers_.find("transfer-encoding");
	    if (x != headers_.end() && x->second == "chunked")
	    {
		respond(501, "Chunked encoding not implemented",  "Chunked encoding not implemented\n", [this](){ });
		return;
	    }
	    try {
		process_request();
	    }
	    catch (std::exception &e)
	    {
		std::ostringstream o;
		o << "Caught exception " << e.what() << "\n";
		std::cerr << o.str();
		respond(500, "Failed", o.str(), [](){});
	    }
	    catch (...)
	    {
		std::ostringstream o;
		o << "Caught default exception\n";
		std::cerr << o.str();
		respond(500, "Failed", o.str(), [](){});
	    }
		    
	}
	else
	{
	    boost::asio::async_read_until(socket_, request_, "\n",
					  boost::bind(&KmerRequest2::read_headers, this,
						      boost::asio::placeholders::error,
						      boost::asio::placeholders::bytes_transferred));
	}
    } 
    else if (err == boost::asio::error::eof)
    {
	// std::cerr << "eof\n";
	exit_request();
    }
    else
    {
	std::cerr << "error " << err << "\n";
	exit_request();
    }
}

void KmerRequest2::process_request()
{
    // std::cerr << "Process request type " << request_type_ << "\n";

    /*
     * Process Expect: 100-continue
     */

    auto it = headers_.find("expect");
    if (it != headers_.end() && it->second == "100-continue")
    {
	// std::cerr << "handling expect 100-continue\n";
	boost::asio::streambuf s;
	std::ostream os(&s);
	os << "HTTP/" << http_version_ << " 100 Continue\n\n";
	boost::asio::write(socket_, s);
    }


    if (request_type_ == "GET")
    {
	boost::smatch match;
	if (path_ == "/quit")
	{
	    respond(200, "OK", "OK, quitting\n", [this](){
		    std::cerr << "stopping io service\n";
		    io_service_.stop();
		});
	}
	else if (path_ == "/version")
	{
	    std::ostringstream os;

	    if (g_parameters->count("kmer-version"))
	    {
		os << "kmer\t" << (*g_parameters)["kmer-version"].as<std::string>() << "\n";
	    }
	    if (g_parameters->count("families-version"))
	    {
		os << "families\t" << (*g_parameters)["families-version"].as<std::string>() << "\n";
	    }
	    os << "family-mode\t" << (server_->family_mode() ? "1" : "0") << "\n";
	    respond(200, "OK", os.str(), [this](){});

	}
	else if (boost::regex_match(path_, match, genus_lookup_path_regex))
	{
	    std::string genus = match[1];

	    auto xmap = mapping_map_->find("");
	    if (xmap == mapping_map_->end())
	    {
		respond(404, "Not Found", "genus not found\n", [this](){});
	    }
	    else
	    {
		auto root_mapping = xmap->second;
		auto hit = root_mapping->genus_map_.find(genus);
		if (hit == root_mapping->genus_map_.end())
		{
		    respond(404, "Not Found", "genus not found\n", [this](){});
		}
		else
		{
		    respond(200, "OK", hit->second + "\n", [this](){});
		}
	    }
	}
	else if (path_ == "/dump_mapping")
	{
	    auto xmap = mapping_map_->find("");
	    auto map = xmap->second;

	    for (auto it: map->kmer_to_id_)
	    {
	        char kmer[10];
		KmerGuts::decoded_kmer(it.first, kmer);
		std::cout << kmer << "\t";
		//os << it.first << "\t";
		for (auto elt: it.second)
		{
		    std::cout << " " << map->decode_id(elt);
		    //std::cout << " " << elt;
		}
		std::cout << "\n";
	    }

	    /* need to go back and fix; this is debugging code anyway
	    for (auto it: map->family_mapping_)
	    {
		std::cout << it.first << "\t" << it.second.pgf << "\t" << it.second.plf << "\t" << it.second.function << "\n";
	    }
	    */
	    respond(200, "OK", "Mapping dumped\n", [this](){ });
	}
	else if (path_ == "/dump_sizes")
	{
	    std::ostringstream os;
	    os << "memory dump\n";
	    for (auto mapping_it : *mapping_map_)
	    {
		os << "Mapping '" << mapping_it.first << "':\n";
		mapping_it.second->dump_sizes(os);
	    }
	    respond(200, "OK", os.str(), [this](){ });
	    
	}
#ifdef BLCR_SUPPORT
	else if (path_ == "/checkpoint")
	{
	    cr_checkpoint_handle_t handle;
	    cr_checkpoint_args_t args;
	    cr_initialize_checkpoint_args_t(&args);

	    std::ostringstream os;
	    os << "checkpoint." << getpid();
	    std::string file(os.str());
	    int fd;
	    fd = open(file.c_str(), O_WRONLY | O_CREAT, 0644);
	    if (fd < 0)
	    {
		std::cerr << "Error opening " << file << ": " << strerror(errno) << "\n";
	    }
	    
	    args.cr_scope = CR_SCOPE_TREE;
	    args.cr_flags = CR_CHKPT_ASYNC_ERR;
	    args.cr_target = getpid();
	    args.cr_fd = fd;
	    args.cr_signal = 0;
	    args.cr_timeout = 0;

	    thread_pool_->image_->detach();
	    io_service_.stop();
	    thread_pool_->io_service_.stop();

	    std::cerr << "Request checkpoint\n";
	    int rc = cr_request_checkpoint(&args, &handle);
	    std::cerr << "Request checkpoint: rc=" << rc << "\n";
	    if (rc == 0)
	    {
		rc = cr_wait_checkpoint(&handle, 0);
		std::cerr << "Wait returns " << rc << "\n";
		if (rc < 0)
		{
		    std::cerr << "Wait error: " << strerror(errno) << "\n";
		}
	    }
	    else
	    {
		std::cerr << "Request failed: " << strerror(errno) << "\n";
	    }
	    thread_pool_->image_->attach();
	    respond(200, "OK", "OK\n", [this]() {});
	}
#endif
	else
	{
	    respond(404, "Not found", "path not found\n", [this](){ });
	}
    }
    else if (request_type_ == "POST")
    {
	auto x = headers_.find("content-length");
	if (x == headers_.end())
	{
	    respond(500, "Missing content length", "Missing content length header\n", [this](){ });
	    return;
	}

	size_t len = std::stoul(x->second);

	std::shared_ptr<KmerPegMapping> mapping;
	std::string key("");
	std::string action(path_);

	boost::smatch match;
	if (boost::regex_match(path_, match, mapping_path_regex))
	{
	    key = match[1];
	    action = match[2];	    
	    std::cerr << "Got keyed mapping '" << key << "' '" << action << "'\n";
	}

	auto xmap = mapping_map_->find(key);
	if (xmap == mapping_map_->end())
	{
	    auto ymap = mapping_map_->insert(std::make_pair(key, std::make_shared<KmerPegMapping>()));
	    mapping = ymap.first->second;
	}
	else
	{
	    mapping = xmap->second;
	}
	
	if (action == "/add")
	{
	    boost::asio::streambuf s;
	    std::ostream os(&s);
	    os << "HTTP/" << http_version_ << " 200 OK\n";
	    os << "Content-type: text/plain\n";
	    os << "\n";
	    boost::asio::write(socket_, s);

	    auto add_request = std::make_shared<AddRequest>(shared_from_this(), mapping, len);
	    add_request->run();
	    active_request_ = add_request;
	}   
	else if (action == "/matrix")
	{
	    auto matrix_request = std::make_shared<MatrixRequest>(shared_from_this(), mapping, len);
	    matrix_request->run();
	    active_request_ = matrix_request;
	}   
	else if (action == "/lookup")
	{
	    auto lookup_request = std::make_shared<LookupRequest>(shared_from_this(), mapping, server_->family_mode(), len);
	    lookup_request->run();
	    active_request_ = lookup_request;
	}   
	else if (action == "/fq_lookup")
	{
	    auto lookup_request = std::make_shared<FqProcessRequest>(shared_from_this(), mapping, server_->family_mode(), len);
	    lookup_request->run();
	    active_request_ = lookup_request;
	}   
	else if (action == "/query")
	{
	    auto query_request = std::make_shared<QueryRequest>(shared_from_this(), len);
	    query_request->run();
	    active_request_ = query_request;
	}
	else
	{
	    respond(404, "Not found", "path not found\n", [this](){ });
	}
    }
}

void KmerRequest2::write_header(std::ostream &os, int code, const std::string &status)
{
    os << "HTTP/" << http_version_ << " " << code << " " << status << "\n";
    os << "Content-type: text/plain\n";
}

void KmerRequest2::respond(int code, const std::string &status, const std::string &result, std::function<void()> on_done)
{
    std::ostringstream os;
    write_header(os, code, status);
    os << "Content-length: " << result.size() << "\n";
    os << "\n";

    std::vector<boost::asio::const_buffer> bufs; 
    std::string s(os.str());
    bufs.push_back(boost::asio::buffer(s));
    bufs.push_back(boost::asio::buffer(result));

    //
    // Use the shared ptr here to keep this object alive through the chain of lambda invocations.
    //
    auto obj = shared_from_this();
    boost::asio::async_write(socket_, bufs,
			     [obj, on_done](const boost::system::error_code &err, const long unsigned int &bytes){
				 // std::cerr << "write all finished\n";
				 obj->exit_request();
				 on_done();
			     });
}

void KmerRequest2::exit_request()
{
    server_->deactivate(shared_from_this());
}
