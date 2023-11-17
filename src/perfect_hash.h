#include <algorithm>
#include <cmph.h>

/*!
  Build perfect hash from signature data in builder.

  We transform the keys of the kept kmers data into a vector, then pass that to the cmph code.
*/

template <int K>
void build_perfect_hash(SignatureBuilder<K> &builder,
			const fs::path &perfect_hash_file,
			const fs::path &data_file)
{
    std::cerr << "build perfect hash into " << perfect_hash_file << " with data in " << data_file << "\n";

    auto &map = builder.kept_kmers();

    std::vector<char *> strvec;
    strvec.reserve(map.size());

    std::transform(map.begin(), map.end(), strvec.begin(), [](auto pair){
	    return strndup(pair.first.data(), K);
	});

    char **d = strvec.data();

    FILE *mphf_fd = fopen(perfect_hash_file.native().c_str(), "wb");
    cmph_io_adapter_t *source = cmph_io_vector_adapter(d, map.size());
    cmph_config_t *config = cmph_config_new(source);
    cmph_config_set_algo(config, CMPH_BDZ);
    cmph_config_set_mphf_fd(config, mphf_fd);
    cmph_t *hash = cmph_new(config);

    /*
     * We have computed our hash. Now write the kmer data.
     */
    
    int hash_size = cmph_size(hash);
    
    StoredKmerData *kd = new StoredKmerData[hash_size];
    std::memset(kd, 0, hash_size);

    std::atomic<int> n = 0;
    tbb::parallel_for(builder.kept_kmers().range(), [kd, hash, &n ](auto r) {
	    for (auto ent = r.begin(); ent != r.end(); ent++)
	    {
		const Kmer<K> &kmer = ent->first;
		const KeptKmer<K> &kept = ent->second;
		unsigned int idx = cmph_search(hash, kmer.data(), K);
		kd[idx] = kept.stored_data;
		n++;
	    }
	});
    std::cerr << "Wrote " << n << " values\n";

    std::FILE *fp = std::fopen(data_file.native().c_str(), "wb");
    if (fp == 0)
    {
	throw std::system_error(errno, std::generic_category(), data_file.native());
    }
    std::fwrite(kd, sizeof(StoredKmerData), hash_size, fp);
    std::fclose(fp);

    cmph_config_destroy(config);
    cmph_dump(hash, mphf_fd);
    cmph_destroy(hash);
    fclose(mphf_fd);
}

