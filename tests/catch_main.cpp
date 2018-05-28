#define CATCH_CONFIG_MAIN
#include <catch.hpp>

#include <functional>

/// class to ensure that the buffer is replaced regardless of how the function exits
class CoutRedirect {
    std::streambuf *old;

  public:
    CoutRedirect(std::streambuf *buf)
        : old(std::cout.rdbuf(buf)) {}
    ~CoutRedirect() { std::cout.rdbuf(old); }
};

std::string capture_stdout(std::function<void()> &&func) {
    std::stringstream buffer;
    CoutRedirect tmp{buffer.rdbuf()};

    func();

    std::string text = buffer.str();
    return text;
}
