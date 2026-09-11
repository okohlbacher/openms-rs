// Test driver only: links the unmodified pinned raw C++ source.
#include <OpenMS/FORMAT/MSNUMPRESS/MSNumpress.h>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
namespace np = ms::numpress::MSNumpress;
double value(const std::string& word) { const std::uint64_t bits = std::stoull(word,nullptr,16); double x; std::memcpy(&x,&bits,8); return x; }
std::string bits(double x) { std::uint64_t raw; std::memcpy(&raw,&x,8); std::ostringstream out; out<<std::hex<<std::setw(16)<<std::setfill('0')<<raw; return out.str(); }
int main() {
 std::string mode,fp,word; std::size_t count;
 while(std::cin>>mode>>fp>>count) {
  std::vector<double> data(count); for(auto& x:data) { std::cin>>word; x=value(word); }
  try {
   std::vector<unsigned char> encoded(std::max(count*8+8,count*5+8)); std::size_t n=0;
   if(mode=="linear") n=np::encodeLinear(data.data(),data.size(),encoded.data(),value(fp));
   else if(mode=="pic") n=np::encodePic(data.data(),data.size(),encoded.data());
   else if(mode=="slof") n=np::encodeSlof(data.data(),data.size(),encoded.data(),value(fp));
   else if(mode=="safe") n=np::encodeSafe(data.data(),data.size(),encoded.data());
   else return 2;
   encoded.resize(n); std::vector<double> decoded(n*2+2); std::size_t m=0;
   if(mode=="linear") m=np::decodeLinear(encoded.data(),n,decoded.data());
   else if(mode=="pic") m=np::decodePic(encoded.data(),n,decoded.data());
   else if(mode=="slof") m=np::decodeSlof(encoded.data(),n,decoded.data());
   else if(n!=0) m=np::decodeSafe(encoded.data(),n,decoded.data()); // C++ empty decode is undefined
   for(auto byte:encoded) std::cout<<std::hex<<std::setw(2)<<std::setfill('0')<<static_cast<unsigned>(byte);
   std::cout<<'\t'; for(std::size_t i=0;i<m;++i) { if(i) std::cout<<','; std::cout<<bits(decoded[i]); }
   std::cout<<'\t'<<bits(np::optimalLinearFixedPoint(data.data(),data.size()))
     <<'\t'<<bits(np::optimalLinearFixedPointMass(data.data(),data.size(),0.001))
     <<'\t'<<bits(np::optimalSlofFixedPoint(data.data(),data.size()))<<'\n'<<std::dec;
  } catch(const char* error) { std::cout<<"ERROR:"<<error<<'\n'; }
 }
}
