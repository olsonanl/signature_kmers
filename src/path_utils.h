#ifndef _path_utils
#define _path_utils

/*!
  @file path_utils.h
  @brief Some useful utilities for managing paths for command line processing.
*/

/*!
  @brief Populate a list of paths from a list of directories.

  For each path in dirs, add any regular file found in that path to paths.

  @param dirs List of directories to load
  @param paths List of paths to populate
*/
void populate_path_list(const std::vector<std::string> &dirs, std::vector<fs::path> &paths)
{
    for (auto dir: dirs)
    {
	fs::path p(dir);
	for (auto dit: fs::directory_iterator(dir))
	{
	    if (fs::is_regular_file(dit.path()))
	    {
		paths.emplace_back(dit);
	    }
	}
    }
}    

void load_strings(const std::vector<std::string> &files, std::vector<std::string> &strings)
{
    for (auto f: files)
    {
	std::ifstream ifstr(f);
	if (ifstr.good())
	{
	    // std::cout << "load " << f << "\n";
	    std::string line;
	    while (std::getline(ifstr, line, '\n'))
	    {
		strings.emplace_back(line);
	    }
	}
	else
	{
	    std::cerr << "could not open " << f << "\n";
	}
    }
}

std::set<std::string> load_set_from_file(const fs::path &file)
{
    /*
     * Read data if file is present.
     */

    std::set<std::string> set;
    if (!file.empty())
    {
	fs::ifstream ifstr(file);
	std::string line;
	while (std::getline(ifstr, line, '\n'))
	{
	    set.emplace(line);
	}
	
    }
    return set;
}

void ensure_directory(const fs::path &dir)
{

    if (!dir.empty())
    {
	if (!fs::is_directory(dir))
	{
	    if (!fs::create_directory(dir))
	    {
		std::cerr << "Error creating " << dir << "\n";
		exit(1);
	    }
	}
    }
}

#endif

