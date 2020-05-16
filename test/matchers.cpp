#include "matchers.hpp"

class IntRange : public Catch::MatcherBase<int> {
    int m_begin, m_end;

public:
    IntRange(int begin, int end)
        : m_begin(begin)
        , m_end(end) { }

    // Performs the test for this matcher
    bool match(int const& i) const override {
        return i >= m_begin && i <= m_end;
    }

    // Produces a string describing what this matcher does. It should
    // include any provided data (the begin/ end in this case) and
    // be written as if it were stating a fact (in the output it will be
    // preceded by the value under test).
    virtual std::string describe() const override {
        std::ostringstream ss;
        ss << "is between " << m_begin << " and " << m_end;
        return ss.str();
    }
};
