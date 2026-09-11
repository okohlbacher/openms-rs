// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#pragma once
#include <string>
namespace OpenMS::StringUtils {
inline bool has(const std::string& value, char c) { return value.find(c) != std::string::npos; }
inline std::string prefix(const std::string& value, char c) { return value.substr(0, value.find(c)); }
template<class T> std::string toStr(T value) { return std::to_string(value); }
}
