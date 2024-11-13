#ifndef _signature_build_h
#define _signature_build_h

#define TBB_PREVIEW_NUMA_SUPPORT 1
#define TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS 1

#include "kmer_data.h"
#include "function_map.h"

#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_unordered_set.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

namespace acc = boost::accumulators;

/*!
  @brief Module encapsulating code for building signature Kmers.

  

*/

/*! @brief Object to keep state for kmers we are keeping.
  We hang onto
  some statistics in order to compute weights later on.
*/
template <int K>
struct KeptKmer
{
    Kmer<K> kmer;
    StoredKmerData stored_data;

//    unsigned int seqs_containing_sig;	// Count of sequences containing this kmer
//    unsigned int seqs_containing_function; // Count of sequences with the kmer that have the function
};

struct KmerStatistics
{
    tbb::atomic<int> distinct_signatures = 0;
    tbb::concurrent_unordered_map<int, int> distinct_functions;
    tbb::concurrent_unordered_map<FunctionIndex, int> seqs_with_func;
    tbb::concurrent_unordered_set<unsigned int> seqs_with_a_signature;
};

template <int K>
using KeptKmers = tbb::concurrent_unordered_map<Kmer<K>, KeptKmer<K>, tbb_hash<K>>;

template <int K>
class SignatureBuilder
{
public:
    SignatureBuilder(int n_threads, int max_seqs_per_file);
    
    using KmerAttributeMap =  tbb::concurrent_unordered_multimap<Kmer<K>, KmerAttributes, tbb_hash<K>>;

    void load_function_data(const std::vector<std::string> &good_functions,
			    const std::vector<std::string> &good_roles,
			    const std::vector<fs::path> &function_definitions);
		   
    void load_fasta(const std::vector<fs::path> &fasta_files, bool keep_functions,
		    const std::set<std::string> &deleted_fids);

    void process_kept_functions(int min_reps_required, const fs::path &function_index_file, std::set<std::string> &ignored_functions);

    void extract_kmers(const std::set<std::string> &deleted_fids);
    void process_kmers();

private:
    void load_kmers_from_fasta(unsigned file_number, const fs::path &file,
			       const std::set<std::string> &deleted_fids);

    void load_kmers_from_sequence(unsigned int &next_sequence_id,
				  const std::string &id, const std::string &def, const std::string &seq);

    struct KmerSet
    {
        KmerSet() : count(0) {}
	void reset() {
	    count = 0;
	    func_count.clear();
	    set.clear();
	}
	// ~KmerSet() { std::cerr << "destroy " << this << "\n"; }
	Kmer<K> kmer;
	std::map<FunctionIndex, int> func_count;
	int count;
	std::vector<KmerAttributes> set;
    };


    KeptKmers<K> kept_kmers_;
    
    void process_kmer_set(KmerSet &set);

    const std::set<unsigned char> ok_prot_ = { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y',
			   'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};

public:
    const KeptKmers<K> &kept_kmers() { return kept_kmers_; }
    const KmerStatistics &kmer_stats() { return kmer_stats_; }
    const std::string lookup_function(FunctionIndex idx) { return fm_.lookup_function(idx); }
    const tbb::concurrent_vector<fs::path> all_fasta_data() { return all_fasta_data_; }
    const FunctionMap &function_map() { return fm_; }

private:

    KmerStatistics kmer_stats_;

    /*! Max allowed sequences per file. Used to assign unique sequence IDs efficiently in parallel.
     */
    int max_seqs_per_file_;
    
    /*! Multimap from a kmer to a set of attributes.
     */
    KmerAttributeMap kmer_attributes_;

    /*! Number of threads to use for processing.
     */
    int n_threads_;
    
    /*! List of all fasta files being processed.
     * Initialized using the load_fasta() method.
     */
    tbb::concurrent_vector<fs::path> all_fasta_data_;
    

    /*! FunctionMap instance used to manage function lists etc.
     */
    FunctionMap fm_;

    /*! Output directory
     */
    std::string kmer_dir_;

    /*! Pathnames to all fasta files being processed.
     */
    tbb::concurrent_vector<fs::path> fasta_data_files_;
    
    
};

#include "signature_build.tcc"

#endif // _signature_build_h
