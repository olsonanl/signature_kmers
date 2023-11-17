#include "cmph_kmer.h"
#include <regex>
#include <fstream>

int main(int argc, char **argv)
{
    if (argc != 4) {
	std::cerr << "usage: tst-cmph basename kmer-file [R|W]\n";
	exit(1);
    }
    fs::path base = argv[1];
    fs::path kmer_file = argv[2];
    char what = argv[3][0];
    CmphKmerDb<StoredKmerData, 8> kmer_db(base);

    if (what != 'R' && what != 'W')
    {
	std::cerr << "usage: tst-cmph basename kmer-file [R|W]\n";
	exit(1);
    }

    if (what == 'W')
    {
	std::cerr << "Create store\n";
	kmer_db.create_backing_data();
	std::cerr << "done\n";
    }
    
    std::cerr << "Create mapping\n";
    kmer_db.map_backing_data();
    std::cerr << "done\n";

    std::ifstream instr(kmer_file);

    std::string line;
    int n = 0;

    std::regex re("\t");
    while (std::getline(instr, line))
    {
	std::vector<std::string> cols;
	std::copy(std::sregex_token_iterator(line.begin(), line.end(), re, -1), {}, std::back_inserter(cols));

	std::string &kmer = cols[0];

	unsigned int idx = kmer_db.lookup_key(kmer);
//	std::cerr << kmer << " " << idx << "\n";

	if (what == 'R')
	{
	    int ec;
	    kmer_db.fetch(kmer, [&kmer](const StoredKmerData &k) { std::cout << kmer << " " << k << "\n"; }, ec);
	}
	else
	{
	    StoredKmerData kd { static_cast<uint16_t>(std::stoul(cols[1])), static_cast<FunctionIndex>(std::stoul(cols[2])) };
//	    std::cerr << "Inserting " << kmer << " " << kd  << "\n";
	    kmer_db.insert(kmer, kd);
	}
	    
//	if (n++ > 10)
//	    break;
    }

    
    
    return 0;
}
