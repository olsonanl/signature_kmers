#ifndef _call_functions_h
#define _call_functions_h

#include "kmer_data.h"

#include "operators.h"
#include "fasta_parser.h"

#include <map>
#include <boost/math/statistics/univariate_statistics.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_map.h>
#include <tbb/parallel_for.h>
#include "seq_id_map.h"

class KmerCall
{
public:
    unsigned int start;
    unsigned int end;
    int count;
    FunctionIndex function_index;
    unsigned int protein_length_median;
    float protein_length_med_avg_dev;
    bool is_short;
    bool is_long;

KmerCall() : start(0), end(0), count(0), function_index(UndefinedFunction),
	     protein_length_median(0), protein_length_med_avg_dev(0.0), is_short(false), is_long(false) {
    }
    KmerCall(unsigned int s, unsigned int e, int c, FunctionIndex f, 
	     unsigned int med = 0, float mad = 0.0, bool is_short = false, bool is_long = false)  :
    start(s), end(e), count(c), function_index(f), 
    protein_length_median(med), protein_length_med_avg_dev(mad),
    is_short(is_short), is_long(is_long) { }
    KmerCall(KmerCall &&k) :
        start(k.start), end(k.end), count(k.count),
	function_index(k.function_index), 
	protein_length_median(k.protein_length_median), protein_length_med_avg_dev(k.protein_length_med_avg_dev),
	is_short(k.is_short), is_long(k.is_long) { }
    KmerCall(const KmerCall &k) :
        start(k.start), end(k.end), count(k.count), function_index(k.function_index),
	protein_length_median(k.protein_length_median), protein_length_med_avg_dev(k.protein_length_med_avg_dev),
	is_short(k.is_short), is_long(k.is_long) { }

};

inline std::ostream &operator<<(std::ostream &os, const KmerCall &c)
{
    os << "KmerCall(" << c.start << "-" << c.end << ": " << c.count << ", " << c.function_index
       << ", " << c.protein_length_median
       << ", " << c.protein_length_med_avg_dev
       << ")";
    return os;
}


template <class KmerDb>
class FunctionCaller
{
public:
    enum OffLengthCallBehavior {
	OffLengthCall = 0,
	OffLengthNoCall = 1,
	OffLengthMark = 2,
    } offlength_call_mode_;

    typedef KmerDb kmer_db_type;

    struct HitInSequence {
	StoredKmerData kmer_data;
	KmerDb::key_type kmer;
	KmerDb::encoded_key_type kidx;
	size_t offset;
    };
    
    FunctionCaller(const KmerDb &db, const boost::filesystem::path &function_index_file,
		   int min_hits = 5, int max_gap = 200);

    void read_function_index(const boost::filesystem::path &function_index_file);
    const std::vector<std::string> &function_index() { return function_index_;}
    
    /* Fasta processing methods */

    template <typename HitCB, typename CallCB>
    void process_fasta_stream(std::istream &istr, HitCB &hit_cb, CallCB &call_cb) const;

    template <typename HitCB, typename CallCB>
    void process_fasta_stream_parallel(std::istream &istr, HitCB &hit_cb, CallCB &call_cb
				       ,SeqIdMap &idmap
	) const;

    template <typename HitCB>
    void process_aa_seq(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			HitCB hit_cb,
			const std::map<std::string, std::string>& parameters =
			std::map<std::string, std::string>()
	) const;

    void find_best_call(const std::string &id, std::vector<KmerCall> &calls, FunctionIndex &function_index,
			std::string &function, float &score, float &score_offset) const;


    const std::string &function_at_index(int idx) const {
	if (idx == UndefinedFunction)
	    return undefined_function_;
	else
	    return function_index_[idx];
    }

    std::string format_call(const KmerCall &c) const;
    std::string format_hit(const HitInSequence &h) const;

#if 0
    void gather_hits(const std::string &seqstr,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
    void process_seq(const char *id,const char *data,
		     std::shared_ptr<std::vector<KmerCall>> calls,
		     std::function<void(const hit_in_sequence_t<Caller> &)> hit_cb,
		     std::shared_ptr<KmerOtuStats> otu_stats);
	
    void process_aa_seq_hits(const std::string &id, const std::string &seq,
			std::shared_ptr<std::vector<KmerCall>> calls,
			     std::shared_ptr<std::vector<hit_in_sequence_t<Caller>>> hits,
			std::shared_ptr<KmerOtuStats> otu_stats);

    std::string format_otu_stats(const std::string &id, size_t size, KmerOtuStats &otu_stats);

    void find_best_call(const std::string &id, std::vector<KmerCall> &calls, FunctionIndex &function_index,
			std::string &function, float &score, float &score_offset);

#endif

private:

    const KmerDb &kmer_db_;

    bool order_constraint_;
    int min_hits_;
    int max_gap_;

    std::vector<std::string> function_index_;
    std::string undefined_function_;
    
};

#include "call_functions.tcc"

#endif // _call_functions_h
