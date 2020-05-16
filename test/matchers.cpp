#include "matchers.hpp"

bool VectorXdEquality::match(const Eigen::VectorXd& actual) const {
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

std::string VectorXdEquality::describe() const {
    std::ostringstream ss;
    ss << "should have size " << m_expected.size() << " and values " << m_expected;
    return ss.str();
}
