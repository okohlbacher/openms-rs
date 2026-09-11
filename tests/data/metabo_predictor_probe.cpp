// SPDX-License-Identifier: BSD-3-Clause
// Test-only driver for unmodified LIBSVM v337, not a full OpenMS SDK build.
#include <svm.h>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if (argc != 2) return 1;
  std::array<svm_model*, 2> models{};
  std::array<std::array<double, 4>, 2> centers{}, scales{};
  for (int i = 0; i < 2; ++i)
  {
    const auto stem = std::string(argv[1]) + "/MetaboliteIsoModelNoised" + (i == 0 ? "2" : "5");
    models[i] = svm_load_model((stem + ".svm").c_str());
    if (!models[i] || models[i]->nr_class != 2 || models[i]->param.kernel_type != RBF ||
        models[i]->label[0] != 2 || models[i]->label[1] != 1) return 2;
    std::ifstream stream(stem + ".scale");
    for (int j = 0; j < 4; ++j)
      if (!(stream >> centers[i][j] >> scales[i][j])) return 3;
  }
  int noise;
  while (std::cin >> std::dec >> noise)
  {
    if (noise != 2 && noise != 5) return 4;
    const int index = noise == 2 ? 0 : 1;
    std::array<svm_node, 5> nodes{};
    std::array<std::uint64_t, 4> bits{};
    for (int j = 0; j < 4; ++j)
    {
      if (!(std::cin >> std::hex >> bits[j])) return 5;
      double raw;
      std::memcpy(&raw, &bits[j], sizeof(raw));
      if (!std::isfinite(raw)) return 6;
      nodes[j] = {j + 1, (raw - centers[index][j]) / scales[index][j]};
    }
    nodes[4].index = -1;
    double decision = 0;
    const auto label = svm_predict_values(models[index], nodes.data(), &decision);
    if (!std::isfinite(decision)) return 7;
    std::uint64_t decision_bits;
    std::memcpy(&decision_bits, &decision, sizeof(decision_bits));
    std::cout << std::dec << noise;
    for (const auto raw_bits : bits) std::cout << '\t' << std::hex << std::setw(16) << std::setfill('0') << raw_bits;
    std::cout << '\t' << std::setw(16) << decision_bits << '\t' << std::dec << static_cast<int>(label) << '\n';
  }
  for (auto& model : models) svm_free_and_destroy_model(&model);
}
