#ifndef _kept_kmer_db_h
#define _kept_kmer_db_h

/*
 * This is a little wrapper so we can use the KeptKmers map as
 * a database for kmer calling.
 */

template <int K>
class KeptKmerDB
{
public:
    static const int KmerSize = K;
    
    KeptKmerDB(const KeptKmers<K> &kk) :
	kept_kmers_(kk) {
    }

    template <typename CB>
    void fetch(const Kmer<K> &k, CB cb, int &ec) const
    {
	auto iter = kept_kmers_.find(k);
	if (iter != kept_kmers_.end())
	{
	    cb(iter->second.stored_data);
	}
	ec = 0;
    };

private:
    const KeptKmers<K> &kept_kmers_;
};


#endif // _kept_kmer_db_h
