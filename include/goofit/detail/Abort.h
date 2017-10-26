#pragma once

#include <string>

namespace GooFit {

class PdfBase;

/// Smart abort that includes the file name and location, and prints a stack trace if possible
void abort(std::string file, int line, std::string reason, const PdfBase *pdf = nullptr);

} // namespace GooFit
