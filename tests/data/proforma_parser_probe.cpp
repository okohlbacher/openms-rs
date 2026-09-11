// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// Harness for unchanged source tokenizer/parser/writer extraction.
using P = OpenMS::ProForma;

std::string hex(const std::string& input)
{
  const char* digits = "0123456789abcdef";
  std::string output;
  for (unsigned char c : input) {
    output += digits[c >> 4];
    output += digits[c & 15];
  }
  return output;
}

int main()
{
  std::string input;
  while (std::getline(std::cin, input)) {
    for (int grammar = 0; grammar < 2; ++grammar) {
      std::cout << hex(input) << '\t' << grammar << '\t';
      try {
        OpenMS::detail::ProFormaParserImpl parser(input);
        std::string lossless, canonical;
        if (grammar == 0) {
          auto value = parser.parsePeptidoform();
          lossless = OpenMS::detail::ProFormaWriter::toString(value, P::WriteMode::LOSSLESS);
          canonical = OpenMS::detail::ProFormaWriter::toString(value, P::WriteMode::CANONICAL);
        } else {
          auto value = parser.parsePeptidoformIon();
          lossless = OpenMS::detail::ProFormaWriter::toString(value, P::WriteMode::LOSSLESS);
          canonical = OpenMS::detail::ProFormaWriter::toString(value, P::WriteMode::CANONICAL);
        }
        std::cout << "ok\t" << hex(lossless) << '\t' << hex(canonical) << '\n';
      } catch (const P::ParseError& error) {
        std::cout << "error\t" << static_cast<int>(error.code) << '\t'
                  << error.position << '\t' << hex(error.what()) << '\n';
      }
    }
  }
}
