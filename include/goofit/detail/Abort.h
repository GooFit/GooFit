#pragma once

#include <string>

class PdfBase;

namespace GooFit {

/// Smart abort that includes the file name and location, and prints a stack trace if possible
void abort(std::string file, int line, std::string reason, const PdfBase* pdf = nullptr);

}

/// Classic abort name
[[deprecated("Use GooFit::abort instead")]]
inline void abortWithCudaPrintFlush(std::string file, int line, std::string reason, const PdfBase* pdf = nullptr) {
    GooFit::abort(file, line, reason, pdf);
}
