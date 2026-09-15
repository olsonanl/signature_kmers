#include "fasta_parser.h"
#include <fstream>
#include <iostream>

int main(int argc, char **argv)
{
    std::ifstream inp(argv[1]);
    FastaParser parser;
    parser.set_def_callback([](const std::string &id, const std::string &def, const std::string &seq) {
	std::cout << "'" << id << "' '" << def << "' '" << seq << "'\n";
    });
    parser.parse(inp);
    parser.parse_complete();
}
