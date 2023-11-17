#include <thread>
#include <iostream>
#include <array>
#include <cstring>

int main(int argc, char **argv)
{

    std::array<char, 4> a{ 'a', 'b', 'c', 'd'};
    char *s = strndup(a.data(), 4);
    free(s);
}
