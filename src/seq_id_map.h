#ifndef _seq_id_map_h
#define _seq_id_map_h

#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_map.h>

class SeqIdMap
{
public:
    SeqIdMap() { }

    int lookup_id(const std::string &id) {
	auto iter = id_to_index_.find(id);
	int idx;
	if (iter == id_to_index_.end())
	{
	    auto ent_iter = index_to_id_.push_back(id);
	    idx = std::distance(index_to_id_.begin(), ent_iter);
	    id_to_index_.emplace(id, idx);
	}
	else
	{
	    idx = iter->second;

	}
	return idx;
    }
    const std::string &lookup_index(int index) {
	return index_to_id_[index];
    }

private:
    tbb::concurrent_vector<std::string> index_to_id_;
    tbb::concurrent_map<std::string, int> id_to_index_;
};

#endif // _seq_id_map_h
