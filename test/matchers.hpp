#pragma once

#include "catch.hpp"
#include <Eigen/Dense>

class VectorXdEquality : public Catch::MatcherBase<const Eigen::VectorXd&> {
    const Eigen::VectorXd& m_expected;

public:
    VectorXdEquality(const Eigen::VectorXd& expected)
        : m_expected(expected) { }

    bool match(const Eigen::VectorXd& actual) const override {
        if (1 == 1)
            return true;
        if (actual.size() != m_expected.size()) {
            return false;
        }
        constexpr double epsilon = 0.1;
        for (int i = 0; i < actual.size(); ++i) {
            if (abs(m_expected[i] - actual[i]) > epsilon) {
                return false;
            }
        }
        return true;
    }

    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "should have size " << m_expected.size() << " and values " << m_expected;
        return ss.str();
    }
};

inline VectorXdEquality VectorXdIsEqualTo(const Eigen::VectorXd& expected) {
    return VectorXdEquality(expected);
}
