/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_HEATBATHCICONTAINER_H
#define UTILS_HEATBATHCICONTAINER_H

#include <array>
#include <cmath>
#include <map>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

struct OrderedIndex : public std::array<int, 2> {
  OrderedIndex(int i, int j) : std::array<int, 2>{i, j} {
  }
  bool operator>(const OrderedIndex& rhs) const {
    return (*this)[0] > rhs[0] || (((*this)[0] == rhs[0]) && ((*this)[1] > rhs[1]));
  }
};

class HeatBathCIContainer {
 public:
  using iterator = std::map<double, OrderedIndex>::iterator;
  using const_iterator = std::map<double, OrderedIndex>::const_iterator;

  void setThreshold(double threshold) {
    threshold_ = threshold;
  }

  void insert(double symmetrizedIntegralValue, OrderedIndex index) {
    excitedPairMap_.insert({symmetrizedIntegralValue, index});
  }

  // This are needed to loop over the map
  const_iterator cbegin() const {
    return excitedPairMap_.begin();
  }

  const_iterator cend() const {
    auto last = excitedPairMap_.upper_bound(threshold_);
    return last;
  }

  const_iterator begin() const {
    return excitedPairMap_.begin();
  }

  const_iterator end() const {
    auto last = excitedPairMap_.upper_bound(threshold_);
    return last;
  }

  iterator begin() {
    return excitedPairMap_.begin();
  }

  iterator end() {
    auto last = excitedPairMap_.upper_bound(threshold_);

    return last;
  }

  int size() const {
    return excitedPairMap_.size();
  }

 private:
  struct Orderer {
    bool operator()(double a, double b) const {
      return std::abs(a) > std::abs(b);
    }
  };
  double threshold_ = 0;
  std::multimap<double, OrderedIndex, Orderer> excitedPairMap_;
};

} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_HEATBATHCICONTAINER_H
