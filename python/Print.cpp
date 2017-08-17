#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Based originally on https://stackoverflow.com/questions/12826751/c-execute-function-any-time-a-stream-is-written-to

#include <streambuf>
#include <ostream>
#include <functional>
#include <string>
#include <memory>

namespace py = pybind11;
using namespace py::literals;


namespace pybind11 {
    
class pythonbuf : public std::streambuf {
private:
    typedef std::streambuf::traits_type traits_type;
    char d_buffer[1024];
    std::string name_;
    
    int overflow(int c) {
        if (!traits_type::eq_int_type(c, traits_type::eof())) {
            *this->pptr() = traits_type::to_char_type(c);
            this->pbump(1);
        }
        return this->sync() ? traits_type::not_eof(c) : traits_type::eof();
    }
    
    int sync() {
        if (this->pbase() != this->pptr()) {
            std::string line(this->pbase(), this->pptr());
            
            object file;
            
            try {
                file = module::import("sys").attr(name_.c_str());
            } catch (const error_already_set &) {
                /* If print() is called from code that is executed as
                 part of garbage collection during interpreter shutdown,
                 importing 'sys' can fail. Give up rather than crashing the
                 interpreter in this case. */
                return 0;
            }
            
            
            auto write = file.attr("write");
            write(line);
            file.attr("flush")();
            
            this->setp(this->pbase(), this->epptr());
        }
        return 0;
    }
public:
    pythonbuf(std::string name = "stdout") : name_(name) {
        this->setp(this->d_buffer, this->d_buffer + sizeof(this->d_buffer) - 1);
    }
};

class opythonstream : private virtual pythonbuf, public std::ostream {
public:
    opythonstream(std::string name = "stdout")
    : pythonbuf(name), std::ostream(static_cast<std::streambuf*>(this)) {
        // this->flags(std::ios_base::unitbuf);
    }
};
}

// Global output
py::opythonstream buffer;
std::streambuf * old = nullptr;

void init_Print(py::module &m) {
    m.def("redirect_stdout", [](bool redirect=true){
        if(redirect && old == nullptr)
            old = std::cout.rdbuf(buffer.rdbuf());
        else if (!redirect && old != nullptr) {
            std::cout.rdbuf(old);
            old = nullptr;
        }
    });

}
