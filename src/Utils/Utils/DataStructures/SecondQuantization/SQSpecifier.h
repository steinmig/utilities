/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.\n
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.\n
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_MO_INTEGRALS_H
#define UTILS_MO_INTEGRALS_H

#include "CasSpecifier.h"
#include <Eigen/Core>
#include <boost/optional.hpp>

namespace Scine {
namespace Utils {
namespace SecondQuantization {

struct MoIntegrals {
  enum class Type { OneBody, TwoBody, ThreeBody };

  auto get(Type type) const -> const Eigen::MatrixXd& {
    if (type == Type::OneBody && oneBodyIntegrals) {
      return *oneBodyIntegrals;
    }
    if (type == Type::TwoBody && twoBodyIntegrals) {
      return *twoBodyIntegrals;
    }
    if (type == Type::ThreeBody && threeBodyIntegrals) {
      return *threeBodyIntegrals;
    }
    throw std::runtime_error("Integrals not implemented");
  }

  boost::optional<Eigen::MatrixXd> oneBodyIntegrals;
  boost::optional<Eigen::MatrixXd> twoBodyIntegrals;
  boost::optional<Eigen::MatrixXd> threeBodyIntegrals;
};

struct MoIntegralsUnrestricted {
  enum class Type { OneBodyAlpha, OneBodyBeta, TwoBodyAlphaAlpha, TwoBodyAlphaBeta, TwoBodyBetaAlpha, TwoBodyBetaBeta };

  auto get(Type type) const -> const Eigen::MatrixXd& {
    if (type == Type::OneBodyAlpha && oneBodyIntegralsAlpha) {
      return *oneBodyIntegralsAlpha;
    }
    if (type == Type::OneBodyBeta && oneBodyIntegralsBeta) {
      return *oneBodyIntegralsBeta;
    }
    if (type == Type::TwoBodyAlphaAlpha && twoBodyIntegralsAlphaAlpha) {
      return *twoBodyIntegralsAlphaAlpha;
    }
    if (type == Type::TwoBodyAlphaBeta && twoBodyIntegralsAlphaBeta) {
      return *twoBodyIntegralsAlphaBeta;
    }
    if (type == Type::TwoBodyBetaAlpha && twoBodyIntegralsBetaAlpha) {
      return *twoBodyIntegralsBetaAlpha;
    }
    if (type == Type::TwoBodyBetaBeta && twoBodyIntegralsBetaBeta) {
      return *twoBodyIntegralsBetaBeta;
    }
    throw std::runtime_error("Integrals not implemented");
  }

  boost::optional<Eigen::MatrixXd> oneBodyIntegralsAlpha;
  boost::optional<Eigen::MatrixXd> oneBodyIntegralsBeta;

  boost::optional<Eigen::MatrixXd> twoBodyIntegralsAlphaAlpha;
  boost::optional<Eigen::MatrixXd> twoBodyIntegralsAlphaBeta;
  boost::optional<Eigen::MatrixXd> twoBodyIntegralsBetaAlpha;
  boost::optional<Eigen::MatrixXd> twoBodyIntegralsBetaBeta;
};

struct Specifier {
  MoIntegrals integrals;
  MoIntegralsUnrestricted unrestrictedIntegrals;
  CasSpecifier<Reference::Restricted> cas;
  CasSpecifier<Reference::Unrestricted> unrestrictedCas;
  Eigen::VectorXd occupations;
  double coreEnergy = 0;

  auto getOneBody() const -> const Eigen::MatrixXd& {
    return integrals.get(MoIntegrals::Type::OneBody);
  }

  auto getOneBodyAlpha() const -> const Eigen::MatrixXd& {
    return unrestrictedIntegrals.get(MoIntegralsUnrestricted::Type::OneBodyAlpha);
  }

  auto getOneBodyBeta() const -> const Eigen::MatrixXd& {
    return unrestrictedIntegrals.get(MoIntegralsUnrestricted::Type::OneBodyBeta);
  }

  auto getTwoBody() const -> const Eigen::MatrixXd& {
    return integrals.get(MoIntegrals::Type::TwoBody);
  }

  auto getTwoBodyAlphaAlpha() const -> const Eigen::MatrixXd& {
    return unrestrictedIntegrals.get(MoIntegralsUnrestricted::Type::TwoBodyAlphaAlpha);
  }

  auto getTwoBodyAlphaBeta() const -> const Eigen::MatrixXd& {
    return unrestrictedIntegrals.get(MoIntegralsUnrestricted::Type::TwoBodyAlphaBeta);
  }

  auto getTwoBodyBetaAlpha() const -> const Eigen::MatrixXd& {
    return unrestrictedIntegrals.get(MoIntegralsUnrestricted::Type::TwoBodyBetaAlpha);
  }
  auto getTwoBodyBetaBeta() const -> const Eigen::MatrixXd& {
    return unrestrictedIntegrals.get(MoIntegralsUnrestricted::Type::TwoBodyBetaBeta);
  }

  auto getThreeBody() const -> const Eigen::MatrixXd& {
    return integrals.get(MoIntegrals::Type::ThreeBody);
  }
};
} // namespace SecondQuantization
} // namespace Utils
} // namespace Scine

#endif // UTILS_MO_INTEGRALS_H
