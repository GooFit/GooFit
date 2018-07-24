#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>

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

#ifdef GOOFIT_MPI
#include <mpi.h>
#endif

int main(int argc, char *argv[]) {
#ifdef GOOFIT_MPI
    printf("MPI_Init!\n");
    MPI_Init(&argc, &argv);
#endif

    int result = Catch::Session().run(argc, argv);

#ifdef GOOFIT_MPI
    printf("MPI_Finalize\n");
    MPI_Finalize();
#endif

    return result;
}
