
#ifndef DEBUG_SCORING
#define DEBUG_SCORING 0
#endif

template <class KmerDb>
class HitSet
{
public:
    struct hit
    {
	StoredKmerData kdata;
	unsigned long pos;
    };
    using HitVector = std::vector<hit>;

    HitSet(int seq_len, int min_hits)
	: seq_len_(seq_len)
	,min_hits_(min_hits)
	{
	}

    void reset();
    template< class... Args >
    void emplace_back( Args&&... args ) {
	hits_.emplace_back(std::forward<Args>(args)...);
    }

    auto rbegin() { return hits_.rbegin(); }

    bool empty() { return hits_.empty(); }
    int count() { return static_cast<int>(hits_.size()); }
    auto clear() { return hits_.clear(); }
    hit& last_hit() { return hits_.back(); }
    void process(const std::string &id, double seqlen, FunctionIndex &current_fI,
		 std::shared_ptr<std::vector<KmerCall>> calls) {

	int fI_count = 0;
	typename HitVector::iterator last_hit;

	std::vector<float> protein_lengths;
	for (auto h_iter = hits_.begin(); h_iter != hits_.end(); h_iter++)
	{
	    if (h_iter->kdata.function_index == current_fI)
	    {
		last_hit = h_iter;
		fI_count++;
		protein_lengths.push_back(static_cast<float>(h_iter->kdata.mean));
	    }
	}
	auto mean_length = boost::math::statistics::mean(protein_lengths);
	auto median_length = boost::math::statistics::median(protein_lengths);
	auto mad_length = boost::math::statistics::median_absolute_deviation(protein_lengths);
	if (mad_length == 0)
	    mad_length = 30;
	auto cutoff_b = mean_length - 2.0 * mad_length;
	auto cutoff_t = mean_length + 2.0 * mad_length;
	// std::cout << "hit stats: " << mean_length << " " << median_length << " " << mad_length << "\n";
	// for (auto x: protein_lengths) { std::cout << x << " " ; }; std::cout << "\n";
	if (fI_count >= min_hits_)
	{
	    if (seqlen < cutoff_b || seqlen > cutoff_t)
	    {
		/*
		std::cerr << "Skip hit" << "\t"
			  << id << "\t"
			  << seqlen << "\t"
			  << cutoff_b << "\t"
			  << cutoff_t << "\t"
			  << current_fI << "\n";
		*/
	    }
	    else
	    {
		if (calls)
		{
		    calls->push_back({
			    static_cast<unsigned int>(hits_[0].pos),
			    static_cast<unsigned int>(last_hit->pos + (KmerDb::KmerSize - 1)),
			    fI_count,
			    current_fI,
			    static_cast<unsigned int>(median_length),
			    mad_length });
		}
	    }
	}
	
	auto end = hits_.rbegin();
	
	if (end[1].kdata.function_index != current_fI &&
	    end[1].kdata.function_index == end[0].kdata.function_index)
	{
	    current_fI = end[1].kdata.function_index;
	    // std::cerr << "reset cur=" << cur_fi << "\n";
	    hits_.erase(hits_.begin(), hits_.end() - 2);
	    // std::cerr << "after erase:\n";
	    //for (auto x: hits)
	    //std::cerr << x << "\n";
	}
	else {
	    hits_.clear();
	}
    }

    HitVector hits_;
    int seq_len_;
    int min_hits_;
};


template <class KmerDb>
FunctionCaller<KmerDb>::FunctionCaller(KmerDb &kmer_db, const fs::path &function_index_file,
		   int min_hits, int max_gap) :
    kmer_db_(kmer_db),
    order_constraint_(false),
    min_hits_(min_hits),
    max_gap_(max_gap),
    ignore_hypothetical_(false)
{
    read_function_index(function_index_file);
}

template <class KmerDb>
void FunctionCaller<KmerDb>::read_function_index(const fs::path &function_index_file)
{
    boost::filesystem::ifstream ifstr(function_index_file);
    std::string line;
    int max_id = 0;
    while (std::getline(ifstr, line, '\n'))
    {
	auto tab = line.find('\t');
	int id = std::stoi(line.substr(0, tab));
	if (id > max_id)
	    max_id = id;
    }
    ifstr.close();
    ifstr.open(function_index_file);

    function_index_.resize(max_id + 1);
    
    while (std::getline(ifstr, line, '\n'))
    {
	auto parts = split(line, "\t");
	int id = std::stoi(parts[0]);
	
	function_index_[id] = parts[1];
    }
}

struct Sequence
{
    std::string id;
    std::string seq;
    Sequence(const std::string &i, const std::string &s) : id(i), seq(s) {}
};

template <class KmerDb>
template <typename HitCB, typename CallCB>
void FunctionCaller<KmerDb>::process_fasta_stream_parallel(std::istream &istr, HitCB &hit_cb, CallCB &call_cb
							   ,SeqIdMap &idmap
    )
{

    try {
	FastaParser parser;
	
	tbb::concurrent_vector<Sequence> seqs;
    
	parser.set_callback([this, &seqs, &idmap ](const std::string &id, const std::string &seq) {

	    if (id.empty())
		return 0;
	    // if (seq.length() > 50)
	    {
		int idx = idmap.lookup_id(id);

		seqs.emplace_back(id, seq);
	    }
	    return 1;
	});
	parser.parse(istr);
	parser.parse_complete();

	tbb::parallel_for(seqs.range(), [this, &hit_cb, &call_cb](auto r) {

	    for (auto entry: r)
	    {
		std::string &id = entry.id;
		std::string &seq = entry.seq;
		
		double slen = static_cast<double>(seq.length());
		auto calls = std::make_shared<std::vector<KmerCall>>();
		
		// std::cerr << "stream " << id << " " << seq << "\n";
		process_aa_seq(id, seq, calls, hit_cb);
		for (auto c: *calls)
		{
		    // std::cout << c << "\n";
		}
		FunctionIndex fi;
		std::string func;
		float score;
		float offset;
		find_best_call(id, *calls, fi, func, score, offset);
		call_cb(id, func, fi, score, seq.size());
		// std::cout << id << "\t" << func << "\t" << fi << "\t" << score << "\n";
	    }
	});;

    }
    catch (std::runtime_error &x)
    {
	std::cerr<< "caught " << x.what() << "\n";
    }
}

template <class KmerDb>
template <typename HitCB, typename CallCB>
void FunctionCaller<KmerDb>::process_fasta_stream(std::istream &istr, HitCB &hit_cb, CallCB &call_cb)
{

    try {
	FastaParser parser;
    
	parser.set_callback([this, &hit_cb, &call_cb](const std::string &id, const std::string &seq) {

	    if (id.empty())
		return 0;
	    double slen = static_cast<double>(seq.length());
	    auto calls = std::make_shared<std::vector<KmerCall>>();

	    // std::cerr << "stream " << id << " " << seq << "\n";
	    process_aa_seq(id, seq, calls, hit_cb);
	    for (auto c: *calls)
	    {
		// std::cout << c << "\n";
	    }
	    FunctionIndex fi;
	    std::string func;
	    float score;
	    float offset;
	    find_best_call(id, *calls, fi, func, score, offset);
	    call_cb(id, func, fi, score, seq.size());
	    // std::cout << id << "\t" << func << "\t" << fi << "\t" << score << "\n";
	    return 0;
	});

	parser.parse(istr);
	parser.parse_complete();
    }
    catch (std::runtime_error &x)
    {
	std::cerr<< "caught " << x.what() << "\n";
    }
}



template <class KmerDb>
template <typename HitCB>
void FunctionCaller<KmerDb>::process_aa_seq(const std::string &idstr, const std::string &seqstr,
					    std::shared_ptr<std::vector<KmerCall>> calls,
					    HitCB hit_cb)
{
    HitSet<KmerDb> hits(seqstr.length(), min_hits_);
    FunctionIndex current_fI = UndefinedFunction;
    double seqlen = static_cast<double>(seqstr.length());

    auto it = std::find(function_index_.begin(), function_index_.end(), "hypothetical protein");
    if (it == function_index_.end())
    {
	std::cerr << "Cannot find hypothetical protein index\n";
	exit(1);
    }
    std::ptrdiff_t hypo_pos = it - function_index_.begin();
    for_each_kmer<KmerDb::KmerSize>(seqstr, [this, &idstr, &calls, &hit_cb, &hits, &current_fI, seqlen, hypo_pos]
				    (const std::array<char, KmerDb::KmerSize> &kmer, size_t offset) {
	// std::cerr << "process " << kmer << "\n";
	
	int ec;
	kmer_db_.fetch(kmer, [this, hit_cb, offset, &idstr, &hits, &calls, &current_fI, &kmer, seqlen, hypo_pos]
		      (const StoredKmerData &kdata) {


	    if (ignore_hypothetical_ && kdata.function_index == hypo_pos)
	    {
		// std::cerr << "Skipping hypo " << kmer << "\t" << offset << "\t" << kdata.function_index << "\n";
		return;
	    }
		
	    hit_cb(idstr, kmer, offset, seqlen, kdata);

	    // std::cerr << kmer << "\t" << offset << "\t" << kdata->function_index << "\n";
	    // Is this hit beyond max_gap_ of the last one?
	    if (!hits.empty() && hits.last_hit().pos + max_gap_ < offset)
	    {
		if (hits.count() >= min_hits_)
		    hits.process(idstr, seqlen, current_fI, calls);
		else
		    hits.clear();
	    }
	    if (hits.empty())
	    {
		current_fI = kdata.function_index;
	    }
	    
	    if (!order_constraint_ ||
		hits.empty() ||
		(kdata.function_index == hits.last_hit().kdata.function_index &&
		 labs((offset - hits.last_hit().pos) -
		      (hits.last_hit().kdata.avg_from_end - kdata.avg_from_end)) <= 20))
	    {
//		using hittype = typename HitSet<KmerDb,K>::hit;
//		hits.emplace_back(hittype {kdata, offset} );
		hits.emplace_back(typename HitSet<KmerDb>::hit{kdata, offset} );
		/*
		 * If we have a pair of new functions, it is time to
		 * process one set and initialize the next.
		 */
		if (hits.count() > 1 && current_fI != kdata.function_index)
		{
		    auto end = hits.rbegin();
		    if (end[1].kdata.function_index == end[0].kdata.function_index)
		    {
			hits.process(idstr, seqlen, current_fI, calls);
		    }
		}
	    }

	}, ec);
	if (ec)
	{
	    // std::cerr << "Error " << ec << " returned from caller\n";
	}
    });
    if (hits.count() >= min_hits_)
	hits.process(idstr, seqlen, current_fI, calls);
}


/*
 * Find the best call from this set of calls.
 *
 * This code replicates the amino acid version of the SEED pipeline
 * km_process_hits_to_regions | km_pick_best_hit_in_peg
 */
template <class KmerDb>
void FunctionCaller<KmerDb>::find_best_call(const std::string &id, std::vector<KmerCall> &calls, FunctionIndex &function_index, std::string &function, float &score, float &score_offset)
{
    function_index = UndefinedFunction;
    function = "";
    score = 0.0;

    if (calls.size() == 0)
    {
	return;
    }
    
#if DEBUG_SCORING
    std::cout << "Initial calls:\n";
    for (auto iter = calls.begin(); iter != calls.end(); iter++)
    {
	std::cout << *iter << "\n";
    }
#endif


    /*
     * First merge adjacent hits that have the same function.
     */
    std::vector<KmerCall> collapsed;

    auto comp = calls.begin();

    while (comp != calls.end())
    {
	// std::cerr << "Inspect start " << * comp << "\n";
	collapsed.push_back(*comp);
	comp++;
	KmerCall &cur = collapsed.back();

	while (comp != calls.end() && cur.function_index == comp->function_index)
	{
	    // std::cerr << "Add call " << *comp << "\n";
	    cur.end = comp->end;
	    cur.count += comp->count;
	    comp++;
	}
    }
#if DEBUG_SCORING
    std::cout << "after collapse:\n";
    for (auto iter = collapsed.begin(); iter != collapsed.end(); iter++)
    {
	std::cout << *iter << "\n";
    }
#endif

    /*
     *
     * Merge hits when we have a case with
     * +------+--+-------+--+-----------+
     * |  F1  |  |   F2  |  |   F1      |
     * +------+--+-------+--+-----------+
     *
     * where the score for F2 is below 5 and the combined scores for
     * the two F1 hits is 10 or more.
     *
     * If that is the case we discard the F2 hit and combine
     * the F1 hits.
     */

    std::vector<KmerCall> merged;

    int merge_interior_thresh = 5;
    int merge_exterior_thresh = 10;

    comp = collapsed.begin();
    while (comp != collapsed.end())
    {
	merged.push_back(*comp);
	comp++;
	auto comp2 = comp + 1;
	KmerCall &cur = merged.back();
	while (comp != collapsed.end() && comp2 != collapsed.end() &&
	       cur.function_index == comp2->function_index &&
	       comp->count < merge_interior_thresh &&
	       (cur.count + comp2->count) >= merge_exterior_thresh)
	{
	    cur.end = comp2->end;
	    cur.count += comp2->count;
	    comp += 2;
	    comp2 = comp + 1;
	}
    }
   
#if DEBUG_SCORING
    std::cerr << "after merge:\n";
    for (auto iter = merged.begin(); iter != merged.end(); iter++)
    {
	std::cerr << *iter << "\n";
    }
#endif

    /*
     * Determine best call.
     *
     * In the current perl kmer_search (km_pick_best_hit_in_peg) code we just take the best
     * function in terms of weighted score. However, that allows tied scores to be called
     * arbitrarily, and close-to-tied scores to be settled based on insignificant differences.
     *
     * In the original kmerv1 code, we required a score threshold difference (typically 5 hits)
     * between the best function and next best function. We resurrect that here.
     *
     */

    /*
     * Fusions
     *
     * We make an attempt to call fusions here.
     */

    if (merged.size() > 1)
    {
	char next_func_key = 'A';
	char next_fusion_key = 'W';
	std::map<std::string, char> func_map;
	std::map<std::string, char> fusion_map;
	std::map<char, std::pair<FunctionIndex, std::string>> key_to_function_info;
	namespace acc = boost::accumulators;
	using accum = acc::accumulator_set<float, acc::stats<acc::tag::mean, acc::tag::variance>>;
	std::map<char, accum> part_stats;
	std::string exp;

	struct fusion_part
	{
	    char ident;
	    const KmerCall call;
	};
	std::vector<fusion_part> fusion_parts;
	int sum_scores = 0;
	for (auto c: merged)
	{
	    /* accumulate scores in the event we call the fusion */
	    sum_scores += c.count;
	       
	    auto func = function_at_index(c.function_index);
	    std::vector<std::string> parts = split(func, " / ");
	    std::string fusion_key;
	    for (auto part: parts)
	    {
		if (func_map.find(part) == func_map.end())
		{
		    char f = next_func_key++;
		    func_map[part] = f;
		}
		fusion_key += func_map[part];
	    }
	    if (parts.size() > 1)
	    {
		if (fusion_map.find(fusion_key) == fusion_map.end())
		{
		    fusion_map[fusion_key] = next_fusion_key++;
		}
		char fkey = fusion_map[fusion_key];
		exp += fkey;
		fusion_parts.emplace_back(fusion_part { fkey, c });
		part_stats[fkey](static_cast<float>(c.protein_length_median));

		key_to_function_info[fkey] = std::make_pair(c.function_index, func);
	    }
	    else
	    {
		exp += func_map[func];
		fusion_parts.emplace_back(fusion_part { func_map[func], c });
		part_stats[func_map[func]](static_cast<float>(c.protein_length_median));
		key_to_function_info[func_map[func]] = std::make_pair(c.function_index, func);
	    }
		
	}
//	std::cout << id << "\t" << exp << "\n";
#if DEBUG_SCORING
	std::cout << "Exp list: " << exp << "\n";
	for (auto x: fusion_parts) {
	    std::cout << x.ident << ": " << x.call << "\n";
	}
#endif

	const boost::regex fusion_re("^W?A[A|W]*W[B|W]*BW?");
	if (regex_match(exp, fusion_re))
	{
	    #if DEBUG_SCORING
	    std::cerr << "potential fusion " << exp << "\n";
	    for (auto p: { 'A', 'W', 'B' })
	    {
		auto stats = part_stats[p];
		std::cerr << p << " " << acc::mean(stats) << " " << acc::variance(stats) << "\n";
	    }
	    #endif
	    auto a_mean = acc::mean(part_stats['A']);
	    auto w_mean = acc::mean(part_stats['W']);
	    auto b_mean = acc::mean(part_stats['B']);
	    auto diff = (a_mean + b_mean) - w_mean;
	    auto frac_dif = abs(diff) / w_mean;
	    if (frac_dif < 0.1)
	    {
#if DEBUG_SCORING
		std::cout << "call fusion " << id << " " << exp << "\n";
		std::cout << a_mean << " " << w_mean << " " << b_mean << " " << diff  << " " << frac_dif << "\n";
		for (auto x: fusion_parts) 
		    std::cout << x.ident << ": " << x.call << "\n";
#endif		
		function_index = key_to_function_info['W'].first;
		function = key_to_function_info['W'].second;
		score = sum_scores;
		score_offset = 0.0;
		return;
	    }
	    else
	    {
		// std::cout << "no-call fusion " << id << " " << exp << "\n";
		// std::cout << a_mean << " " << w_mean << " " << b_mean << " " << diff  << " " << frac_dif << "\n";
	    }
	}
	    
    }

    typedef std::map<int, int> map_t;

    map_t by_func;

    for (auto c: merged)
    {
	auto it = by_func.find(c.function_index);
	if (it == by_func.end())
	{
	    by_func.insert(std::make_pair(c.function_index, c.count));
	}
	else
	{
	    it->second += c.count;
	}
    }

    //typedef map_t::value_type ent_t;
    typedef std::pair<FunctionIndex, int> ent_t;

    std::vector<ent_t> vec;
    for (auto it = by_func.begin(); it != by_func.end(); it++)
	vec.push_back(*it);

//    std::cerr << "vec len " << vec.size() << "\n";
    if (vec.size() > 1)
    {
	std::partial_sort(vec.begin(), vec.begin() +  2, vec.end(), 
			  [](const ent_t& s1, const ent_t& s2) {
			      return (s1.second > s2.second); });
    }
    
#if DEBUG_SCORING
    for (auto x: vec)
    {
	std::cerr << x.first << " " << x.second << " ";
	std::cerr << function_at_index(x.first) << "\n";
    }
#endif
    
    if (vec.size() == 1)
	score_offset = (float) vec[0].second;
    else
	score_offset = (float) (vec[0].second - vec[1].second);

#if DEBUG_SCORING
    std::cerr << "Offset=" << score_offset << "\n";
#endif

    if (score_offset >= 5.0)
    {
	auto best = vec[0];
	function_index = best.first;
	function = function_at_index(function_index);
	score = (float) best.second;
    }
    else
    {
	function_index = UndefinedFunction;
	function = "";
	score = 0.0;

	/*
	 * Try to compute a fallback function naming the two best hits if there are two hits within the
	 * threshold but greater than the next hit.
	 */
#if 1
	if (vec.size() >= 2)
	{
	    std::string f1 = function_at_index(vec[0].first);
	    std::string f2 = function_at_index(vec[1].first);
	    if (f2 > f1)
		std::swap(f1, f2);

	    if (vec.size() == 2)
	    {
		function = f1 + " ?? " + f2;
		score = (float) vec[0].second;
	    }
	    else if (vec.size() > 2)
	    {
		float pair_offset = (float) (vec[1].second - vec[2].second);
		if (pair_offset > 2.0)
		{
		    function = f1 + " ?? " + f2;
		    score = (float) vec[0].second;
		    score_offset = pair_offset;
		}
	    }
	}
#endif
    }
}

