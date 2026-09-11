// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// Standalone harness appended to exact source AST and writer extraction.
using namespace OpenMS;
using P = ProForma;

P::Modification scored(P::ModificationTag tag, double value)
{
  P::Modification mod;
  mod.alternatives.emplace_back(tag, P::Label{P::Label::Type::AMBIGUOUS, "g1", value});
  return mod;
}
P::Peptidoform peptide()
{
  P::Peptidoform p;
  p.sequence.emplace_back(P::SequenceElement{'X', {}});
  return p;
}
int main()
{
  const uint64_t bits[] = {
    0x0000000000000000ULL, 0x8000000000000000ULL,
    0x3fc0000000000000ULL, 0xbfc0000000000000ULL,
    0x3ff0147ae147ae14ULL, 0xbff0147ae147ae14ULL,
    0x3f1a36e2bd504417ULL, 0x412e847fcccccccdULL,
    0x0000000000000001ULL, 0x8000000000000001ULL,
    0x0010000000000000ULL, 0x8010000000000000ULL,
    0x7fefffffffffffffULL, 0xffefffffffffffffULL,
    0x400921fb54442d18ULL, 0xc00921fb54442d18ULL,
    0x4340000000000000ULL, 0x3f847ae147ae147bULL,
    0x3ee4f8b588e368f1ULL, 0x40c3880000000000ULL,
  };
  for (uint64_t raw : bits)
  {
    double value;
    std::memcpy(&value, &raw, sizeof(value));
    for (int mode = 0; mode < 2; ++mode)
    {
      auto write_mode = mode == 0 ? P::WriteMode::LOSSLESS : P::WriteMode::CANONICAL;
      for (int kind = 0; kind < 4; ++kind)
      {
        auto p = peptide();
        auto& mods = std::get<P::SequenceElement>(p.sequence[0]).modifications;
        P::MassDelta delta{P::MassDelta::Source::OBS, value, kind == 2 ? "+001.2300" : ""};
        if (kind == 1) mods.push_back(scored(P::InfoTag{"before"}, value));
        mods.push_back(scored(delta, value));
        mods.push_back(scored(P::InfoTag{"after"}, value));
        std::string text;
        if (kind == 3)
        {
          auto second = peptide();
          std::get<P::SequenceElement>(second.sequence[0]).modifications.push_back(scored(P::InfoTag{"fresh"}, value));
          P::PeptidoformIon ion;
          ion.name = "omitted";
          ion.chains = {p, second};
          ion.chains[0].charge = P::ChargeState{2};
          ion.chains[1].charge = P::ChargeState{-1};
          ion.is_chimeric = true;
          ion.charge = P::ChargeState{3};
          text = detail::ProFormaWriter::toString(ion, write_mode);
        }
        else text = detail::ProFormaWriter::toString(p, write_mode);
        std::cout << std::hex << raw << std::dec << '\t' << mode << '\t' << kind << '\t' << text << '\n';
      }
    }
  }
}
