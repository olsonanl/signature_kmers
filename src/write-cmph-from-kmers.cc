#include "cmph_kmer.h"
#include <regex>
#include <fstream>

/*! Read a final.kmers file and create the mmap data file.
 */

int main(int argc, char **argv)
{
    if (argc != 3) {
	std::cerr << "usage: tst-cmph basename kmer-file\n";
	exit(1);
    }
    fs::path base = argv[1];
    fs::path kmer_file = argv[2];

    CmphKmerDb<StoredKmerData, 8> kmer_db(base);

    std::ifstream instr(kmer_file);

    std::string line;
    int n = 0;

    int hash_size = kmer_db.hash_size();
    StoredKmerData *kd = new StoredKmerData[hash_size];
    std::memset(kd, 0, hash_size);

    std::regex re("\t");
    while (std::getline(instr, line))
    {
	std::vector<std::string> cols;
	std::copy(std::sregex_token_iterator(line.begin(), line.end(), re, -1), {}, std::back_inserter(cols));

	std::string &kmer = cols[0];

	unsigned int idx = kmer_db.lookup_key(kmer);
	kd[idx].avg_from_end = static_cast<uint16_t>(std::stoul(cols[1]));
	kd[idx].function_index =static_cast<FunctionIndex>(std::stoul(cols[2]));
	if (n++ % 100000 == 0)
	{
	    std::cerr << n << "\n";
	}
	//if (n > 1000000)
//	    break;
    }
    std::cerr << "done initializing\n";

    return 0;
}
