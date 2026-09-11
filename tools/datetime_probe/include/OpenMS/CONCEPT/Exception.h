// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#pragma once
#include <stdexcept>
namespace OpenMS::Exception {
struct ParseError : std::runtime_error { template<class... T> ParseError(T&&...) : std::runtime_error("ParseError") {} };
struct InvalidValue : std::runtime_error { template<class... T> InvalidValue(T&&...) : std::runtime_error("InvalidValue") {} };
}
