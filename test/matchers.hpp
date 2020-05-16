#pragma once

#include "catch.hpp"
#include <Eigen/Dense>

class VectorXdEquality : public Catch::MatcherBase<const Eigen::VectorXd&> {
    const Eigen::VectorXd& m_expected;

public:
    VectorXdEquality(const Eigen::VectorXd& expected)
        : m_expected(expected) { }

    bool match(const Eigen::VectorXd& actual) const override;

    virtual std::string describe() const override;
};

inline VectorXdEquality VectorXdIsEqualTo(const Eigen::VectorXd& expected) {
    return VectorXdEquality(expected);
}
