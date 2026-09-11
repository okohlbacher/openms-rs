// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <locale>
#include <clocale>
static std::string hex(const std::string& text) {
 const char* digits="0123456789abcdef"; std::string out;
 for (unsigned char c:text) {out+=digits[c>>4];out+=digits[c&15];}return out;
}
static std::string unhex(const std::string& text) {
 std::string out;for(size_t i=0;i<text.size();i+=2)out+=static_cast<char>(std::stoi(text.substr(i,2),nullptr,16));return out;
}
static std::vector<std::string> split(const std::string& text,char delimiter) {
 std::vector<std::string> out;size_t start=0;
 for(;;) {auto pos=text.find(delimiter,start);out.push_back(text.substr(start,pos-start));if(pos==std::string::npos)return out;start=pos+1;}
}
int main() {
 std::setlocale(LC_ALL,"C");
 std::string line;
 const std::vector<std::string> formats={"yyyy-MM-ddThh:mm:ss","yyyy-MM-ddThh:mm:ss.zzz","yyyy-MM-dd hh:mm:ss","yyyy-MM-dd+hh:mm","yyyy-MM-ddThh:mm:ssZ","yyyy-MM-dd","hh:mm:ss"};
 while(std::getline(std::cin,line)) {
  auto f=split(line,'\t'); if(f.size()!=6)return 2;
  OpenMS::DateTime d;std::string status="ok";
  try {
   if(!f[1].empty()) { auto seed=unhex(f[1]); if(seed.starts_with("time:"))d.setTime(seed.substr(5)); else d.set(seed); }
   auto text=unhex(f[3]);auto format=unhex(f[4]);
   if(f[2]=="set")d.set(text);
   else if(f[2]=="date")d.setDate(text);
   else if(f[2]=="time")d.setTime(text);
   else if(f[2]=="from")d=OpenMS::DateTime::fromString(text,format);
   else if(f[2]=="add")d.addSecs(std::stoi(f[5]));
   else if(f[2]=="clear")d.clear();
   else if(f[2]=="parts"||f[2]=="dateparts"||f[2]=="timeparts") {
    auto n=split(f[5],',');std::vector<unsigned> nums;
    for(auto v:n)nums.push_back(static_cast<unsigned>(std::stoul(v)));
    if(f[2]=="parts")d.set(nums.at(0),nums.at(1),nums.at(2),nums.at(3),nums.at(4),nums.at(5));
    else if(f[2]=="dateparts")d.setDate(nums.at(0),nums.at(1),nums.at(2));
    else d.setTime(nums.at(0),nums.at(1),nums.at(2));
   } else if(f[2]!="none")return 3;
  } catch(const std::exception& e) {status=e.what();}
  OpenMS::UInt m,day,y,h,min,s;d.get(m,day,y,h,min,s);
  std::cout<<line<<'\t'<<status<<'\t'<<d.isValid()<<'\t'<<static_cast<int>(m)<<','<<static_cast<int>(day)<<','<<static_cast<int>(y)<<','<<static_cast<int>(h)<<','<<static_cast<int>(min)<<','<<static_cast<int>(s);
  for(auto format:formats)std::cout<<'\t'<<hex(d.toString(format));
  std::cout<<'\t'<<hex(d.get())<<'\t'<<hex(d.getDate())<<'\t'<<hex(d.getTime())<<'\n';
 }
}
