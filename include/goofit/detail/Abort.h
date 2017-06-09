#pragma once

#include <string>

namespace GooFit {

class PdfBase;

/// Smart abort that includes the file name and location, and prints a stack trace if possible
void abort(std::string file, int line, std::string reason, const PdfBase *pdf = nullptr);

} // namespace GooFit

/// Classic abort name
[[deprecated("Use GooFit::abort instead")]] inline void
abortWithCudaPrintFlush(std::string file, int line, std::string reason, const GooFit::PdfBase *pdf = nullptr) {
    GooFit::abort(file, line, reason, pdf);
}
