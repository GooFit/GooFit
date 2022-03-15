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
#include <fmt/ostream.h>

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
    auto get_value() const -> fptype { return m_value; }

    auto get_uncertainty() const -> fptype { return m_uncertainty; }

    auto get_relative_uncertainty() const -> fptype { return m_uncertainty / m_value; }

    // Equality operators
    auto operator==(const Uncertain &other) const -> bool {
        return (m_value == other.m_value) && (m_uncertainty == other.m_uncertainty);
    }

    auto operator!=(const Uncertain &other) const -> bool { return !operator==(other); }

    // Arithmetic operators
    auto operator+(const Uncertain &other) const -> Uncertain {
        return {m_value + other.m_value, m_uncertainty + other.m_uncertainty};
    }

    auto operator-(const Uncertain &other) const -> Uncertain {
        return {m_value - other.m_value, m_uncertainty + other.m_uncertainty};
    }

    auto operator*(const Uncertain &other) const -> Uncertain {
        return {m_value * other.m_value, get_relative_uncertainty() + other.get_relative_uncertainty()};
    }

    /// Allow int and float multiplies
    auto operator*(fptype other) const -> Uncertain { return {m_value * other, m_uncertainty * other}; }

    auto operator/(const Uncertain &other) const -> Uncertain {
        return {m_value / other.m_value, get_relative_uncertainty() + other.get_relative_uncertainty()};
    }

    auto operator/(const fptype &other) const -> Uncertain { return {m_value / other, m_uncertainty / other}; }
};

/// Allow int and float multiplies
inline auto operator*(fptype other, const Uncertain &self) -> Uncertain {
    return {self.get_value() * other, self.get_uncertainty() * other};
}

/// Simple << output, will also support fmt printout
inline auto operator<<(std::ostream &stream, Uncertain value) -> std::ostream & {
    fptype val = value.get_value();
    fptype err = value.get_uncertainty();

    if(err > 0 && err < 2.5) {
        auto numdig = (size_t)std::ceil(-log10(err) + log10(2.5));
        return stream << fmt::format("{0:.{2}f} +/- {1:.{2}f}", val, err, numdig);
    } else if(err >= 2.5) {
        return stream << fmt::format("{0:.0f} +/- {1:.0f}", val, err);
    } else {
        return stream << val << " +/- " << err;
    }
}

} // namespace GooFit
