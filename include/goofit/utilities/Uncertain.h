// The MIT License (MIT)
//
// Copyright (c) 2015 Claus Jensby Madsen
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// https://github.com/clausjensbymadsen/uncertain-values

#pragma once

#include <fmt/format.h>
#include <iosfwd>
#include <stdexcept>
#include <string>

namespace GooFit {

class Variable;

class Uncertain {
  private:
    fptype m_value;
    fptype m_uncertainty;

  public:
    // Constructors
    Uncertain(fptype value, fptype uncertainty = 0)
        : m_value(value)
        , m_uncertainty(uncertainty) {
        if(m_uncertainty < 0)
            throw std::invalid_argument("Uncertainty can not be negative.");
    }

    Uncertain(const Variable &val)
        : Uncertain(val.getValue(), val.getError()) {}

    // Class getters
    fptype get_value() const { return m_value; }

    fptype get_uncertainty() const { return m_uncertainty; }

    fptype get_relative_uncertainty() const { return m_uncertainty / m_value; }

    // Equality operators
    bool operator==(const Uncertain &other) const {
        return (m_value == other.m_value) && (m_uncertainty == other.m_uncertainty);
    }

    bool operator!=(const Uncertain &other) const { return !operator==(other); }

    // Arithmetic operators
    Uncertain operator+(const Uncertain &other) const {
        return Uncertain(m_value + other.m_value, m_uncertainty + other.m_uncertainty);
    }

    Uncertain operator-(const Uncertain &other) const {
        return Uncertain(m_value - other.m_value, m_uncertainty + other.m_uncertainty);
    }

    Uncertain operator*(const Uncertain &other) const {
        return Uncertain(m_value * other.m_value, get_relative_uncertainty() + other.get_relative_uncertainty());
    }

    /// Allow int and float multiplies
    Uncertain operator*(fptype other) const { return Uncertain(m_value * other, m_uncertainty * other); }

    Uncertain operator/(const Uncertain &other) const {
        return Uncertain(m_value / other.m_value, get_relative_uncertainty() + other.get_relative_uncertainty());
    }

    Uncertain operator/(const fptype &other) const { return Uncertain(m_value / other, m_uncertainty / other); }
};

/// Allow int and float multiplies
inline Uncertain operator*(fptype other, const Uncertain &self) {
    return Uncertain(self.get_value() * other, self.get_uncertainty() * other);
}

/// Simple << output
inline std::ostream &operator<<(std::ostream &stream, Uncertain value) {
    return stream << value.get_value() << " ± " << value.get_uncertainty();
}

/// fmt support
inline void format_arg(fmt::BasicFormatter<char> &f, const char *&format_str, const Uncertain &s) {
    std::string value{format_str};
    auto end        = value.find('}');
    std::string val = std::string{"{" + value.substr(0, end) + "}"};

    f.writer().write(val.c_str(), s.get_value());
    f.writer().write(" ± ");
    f.writer().write(val.c_str(), s.get_uncertainty());

    format_str += end + 1; // Remove :xxx}
}

} // namespace GooFit
