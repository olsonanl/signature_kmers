
template <int K>
SignatureBuilder<K>::SignatureBuilder(int n_threads, int max_seqs_per_file) :
    n_threads_(n_threads),
    max_seqs_per_file_(max_seqs_per_file)
{
}


template <int K>
void SignatureBuilder<K>::load_function_data(const std::vector<std::string> &good_functions,
					     const std::vector<std::string> &good_roles,
					     const std::vector<fs::path> &function_definitions)
{
    fm_.add_good_roles(good_roles);
    fm_.add_good_functions(good_functions);

    for (auto def: function_definitions)
    {
	fm_.load_id_assignments(def);
    }

}

template <int K>
void SignatureBuilder<K>::load_fasta(const std::vector<fs::path> &fasta_files,
				     bool keep_functions,
				     const std::set<std::string> &deleted_fids)
{
    for (auto fasta: fasta_files)
    {
	fm_.load_fasta_file(fasta, false, deleted_fids);
	all_fasta_data_.emplace_back(fasta);
    }
}

template <int K>
void SignatureBuilder<K>::process_kept_functions(int min_reps_required, const fs::path &output_dir)
{
    fm_.process_kept_functions(min_reps_required);
    if (!output_dir.empty())
    {
	fm_.write_function_index(output_dir);
    }
}

template <int K>
void SignatureBuilder<K>::extract_kmers(const std::set<std::string> &deleted_fids)
{
    if (n_threads_ < 2)
    {
	for (unsigned i = 0; i < (unsigned) all_fasta_data_.size(); i++)
	{
	    load_kmers_from_fasta(i, all_fasta_data_[i], deleted_fids);
	}
    }
    else
    {
	size_t n = all_fasta_data_.size();
	tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
			  [this, &deleted_fids](const tbb::blocked_range<size_t> &r) {
			      for (size_t i = r.begin(); i != r.end(); ++i)
			      {
				  auto fasta = all_fasta_data_[i];
				  // std::cout << "load file " << i << " " << fasta << "\n";
				  load_kmers_from_fasta((unsigned) i, fasta, deleted_fids);
			      }
			  });
    }
}

/*!
  @brief Load sequence data.

  Parse the fasta data using FastaParser.

  For each sequence, invoke load_sequence() to create and store kmers.
  
  If a function index is not defined, this is not a function we wish
  to process so skip this sequence.
  
 */
template <int K>
void SignatureBuilder<K>::load_kmers_from_fasta(unsigned file_number, const fs::path &file,
						const std::set<std::string> &deleted_fids)
{
    fs::ifstream ifstr(file);

    FastaParser parser;
    
    unsigned next_sequence_id = file_number * max_seqs_per_file_;

    parser.set_def_callback([this, &next_sequence_id, &deleted_fids](const std::string &id, const std::string &def, const std::string &seq) {
	if (deleted_fids.find(id) == deleted_fids.end())
	{
	    load_kmers_from_sequence(next_sequence_id, id, def, seq);
	}
	return 0;
    });
    parser.parse(ifstr);
    parser.parse_complete();
}

/*!
  @brief Load a single sequence.

  Use the provided FunctionMap to look up the function for the given sequence. If it is not present,
  skip this sequence.

  Update KmerStatistics::seqs_with_func to reflect the additional sequence having this function.

  For each kmer in the sequence,

  - If kmer has no invalid characters, insert it into the @ref KmerAttributeMap. This logs the existence
  of the kmer with the given function, offset from the end of its protein, and length of the source protein.

  
*/

template <int K>
void SignatureBuilder<K>::load_kmers_from_sequence(unsigned int &next_sequence_id,
						   const std::string &id, const std::string &def, const std::string &seq)
{
    if (id.empty())
	return;

    std::string func = fm_.lookup_function(id);
    
    /*
     * Empty means empty (and perhaps deleted feature).
     */
    
    if (func.empty())
    {
	return;
    }
    
    unsigned int seq_id = next_sequence_id++;

    FunctionIndex function_index = fm_.lookup_index(func);

    if (false)
    {
	if (function_index == UndefinedFunction)
	{
	    function_index = fm_.lookup_index("hypothetical protein");
	    if (function_index == UndefinedFunction)
	    {
		std::cerr << "No function defined for hypothetical protein\n";
		exit(1);
	    }
	}
    }

    if (function_index == UndefinedFunction)
    {
    	return;
    }

    kmer_stats_.seqs_with_func[function_index]++;

    for (auto it = seq.begin(); it < seq.end() - K + 1; it++)
    {
        unsigned short n = (unsigned short) std::distance(it, seq.end());
	Kmer<K> kmer;
	bool ok = true;
	std::copy_n(it, K, kmer.begin());
	for (auto x: kmer)
	{
	    if (ok_prot_.find(x) == ok_prot_.end())
	    {
		ok = false;
		break;
	    }
	}
	if (ok)
	{
	    kmer_attributes_.insert({kmer, { function_index, UndefinedOTU, n, seq_id, static_cast<unsigned int>(seq.length())}});
	}
    }
}

template <int K>
void SignatureBuilder<K>::process_kmers()
{
    tbb::parallel_for(kmer_attributes_.range(), [this](auto r) {
	    KmerSet cur_set;
	    Kmer<K> cur { 0 };
	    for (auto ent = r.begin(); ent != r.end(); ent++)
	    {
		const Kmer<K> &kmer = ent->first;
		KmerAttributes &attr = ent->second;

		if (kmer != cur)
		{
		    if (cur_set.count > 0)
			process_kmer_set(cur_set);
		    
		    cur_set.reset();
		    cur_set.kmer = kmer;
		    cur = kmer;
		}
		cur_set.func_count[attr.func_index]++;
		cur_set.count++;
		cur_set.set.emplace_back(attr);
	    }
	    process_kmer_set(cur_set);
	});

    std::cout << "Kept " << kept_kmers_.size() << " kmers\n";
    std::cout << "distinct_signatures=" << kmer_stats_.distinct_signatures << "\n";
    std::cout << "num_seqs_with_a_signature=" << kmer_stats_.seqs_with_a_signature.size() << "\n";
}

/*! @brief Process a set of instances of a given kmer.

 */
template <int K>
void SignatureBuilder<K>::process_kmer_set(KmerSet &set)
{
    FunctionIndex best_func_1 = UndefinedFunction, best_func_2 = UndefinedFunction;
    int best_count_1 = -1, best_count_2 = -1;
    

    // we want the top two elements by value in the map; don't know how
    // with standard STL without copying to a vector, but if we're copying
    // anyway we can just search for them.
    for (auto x: set.func_count)
    {
	if (best_func_1 == UndefinedFunction)
	{
	    best_func_1 = x.first;
	    best_count_1 = x.second;
	}
	else if (x.second > best_count_1)
	{
	    best_func_2 = best_func_1;
	    best_count_2 = best_count_1;

	    best_func_1 = x.first;
	    best_count_1 = x.second;
	}
	else if (x.second > best_count_2)
	{
	    best_func_2 = x.first;
	    best_count_2 = x.second;
	}
    }

    float thresh = float(set.count) * 0.8f;
    int best_count = best_count_1;
    FunctionIndex best_func = best_func_1;

    if ((float) best_count < thresh)
    {
	    return;
    }

//    unsigned int seqs_containing_func = 0;
    std::vector<unsigned short> offsets;

    acc::accumulator_set<unsigned short, acc::stats<acc::tag::mean,
						  acc::tag::median,
						  acc::tag::variance> > acc;

    for (auto item: set.set)
    {
	if (item.func_index == best_func)
	{
//	    seqs_containing_func++;
	    acc(item.protein_length);
	}
	offsets.push_back(item.offset);
	kmer_stats_.seqs_with_a_signature.insert(item.seq_id);
    }

    unsigned short mean = acc::mean(acc);
    unsigned short median = acc::median(acc);
    unsigned short var = acc::variance(acc);

    std::sort(offsets.begin(), offsets.end());
    unsigned short avg_from_end = offsets[offsets.size() / 2];
    // std::cout << seqs_containing_func << " " << avg_from_end<< "\n";

    kmer_stats_.distinct_signatures++;
    kmer_stats_.distinct_functions[best_func]++;

    kept_kmers_.emplace(set.kmer, KeptKmer<K> { set.kmer, { avg_from_end, best_func, mean, median, var }
	    // , (unsigned int) set.set.size()
	    // , seqs_containing_func
	});

}

