// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#include <ctime>
#include <cerrno>
#include <iostream>
int main(){for(int year : {-1,0,1,400,1800,1899,1900,1901,1969,1970,2147483647}){
 std::tm t{};t.tm_year=year-1900;t.tm_mon=0;t.tm_mday=1;t.tm_isdst=-1;errno=0;
 std::time_t result=timegm(&t);
 std::cout<<year<<'\t'<<result<<'\t'<<errno<<'\n';
}}
