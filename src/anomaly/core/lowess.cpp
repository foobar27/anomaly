#include "lowess.hpp"

namespace anomaly::core::lowess {

using Scalar = Eigen::VectorXd::Scalar;
using Index  = Eigen::Index;

constexpr double square(const Scalar x) {
    return x * x;
}

constexpr double cube(const Scalar x) {
    return x * x * x;
}

constexpr double biSquare(const Scalar x) {
    return square(1.0 - square(x));
}

constexpr double triCube(Scalar x) {
    return cube(1.0 - cube(x));
}

bool LowessAlgorithm::lowest(const Eigen::VectorXd& positions,
                             const Eigen::VectorXd& input,
                             const Scalar           position,
                             Scalar&                fittedValue,
                             const Index            n_left,
                             const Index            n_right,
                             Eigen::VectorXd&       weights,
                             const bool             use_rw) {
    auto n     = positions.size();
    auto range = positions[n - 1] - positions[0];
    auto h     = std::max(position - positions[n_left], positions[n_right - 1] - position);
    auto h9    = 0.999 * h;
    auto h1    = 0.001 * h;

    Index n_rt{};

    Scalar sumOfWeights = 0.0; // sum of weights
    Index  j            = n_left; // The loop might terminate earlier, we need to remember the last loop position.
    for (; j < n; j++) { // compute weights (pick up all ties on right)
        weights[j] = 0.0;
        auto r     = abs(positions[j] - position);
        if (r <= h9) { // small enough for non-zero weight
            if (r > h1)
                weights[j] = triCube(r / h);
            else
                weights[j] = 1.0;
            if (use_rw)
                weights[j] *= m_robustnessWeights[j];
            sumOfWeights += weights[j];
        } else if (positions[j] > position)
            break; // get out at first zero wt on right
    }
    // TODO shift:
    n_rt = j + 1 - 1; // rightmost pt (may be greater than nright because of ties)
    if (sumOfWeights <= 0.0)
        return false;
    // make sum of w(j) == 1
    for (Index j = n_left; j < n_rt; ++j)
        weights[j] /= sumOfWeights; // TODO(sw) simplify scalar multiplication

    // weighted least squares
    if (h > 0.0) { // use linear fit
        double a = 0.0;
        for (Index j = n_left; j < n_rt; ++j) // TODO(sw) simplify via dot product
            a += weights[j] * positions[j]; // weighted center of x values
        auto   b = position - a;
        double c = 0.0;
        for (Index j = n_left; j < n_rt; ++j) {
            c += weights[j] * square(positions[j] - a); // TODO(sw) simplify dot product
        }

        if (sqrt(c) > .001 * range) {
            // points are spread out enough to compute slope
            b /= c;
            for (Index j = n_left; j < n_rt; ++j)
                weights[j] = weights[j] * (1.0 + b * (positions[j] - a));
        }
    }

    // Compute output
    fittedValue = 0.0;
    for (Index j = n_left; j < n_rt; ++j) {
        fittedValue += weights[j] * input[j]; // TODO(sw) simplify dot product
    }
    return true;
}

void LowessAlgorithm::lowess(const Eigen::VectorXd& positions, const Eigen::VectorXd& input, const double f, const Index nSteps, const double delta) {
    auto n = m_numberOfPoints;
    assert(positions.size() == n);
    assert(input.size() == n);
    assert(m_output.size() == n);
    assert(m_robustnessWeights.size() == n);
    assert(m_residuals.size() == n);
    if (n < 2) {
        m_output[0] = input[0];
        return;
    }
    Index ns = std::max(std::min(Index(f * n), n), Index(2)); // at least two, at most n points

    // robustness iterations
    for (Index iter = 1; iter <= nSteps + 1; iter = iter + 1) {
        Index n_left  = 0;
        Index n_right = ns;
        Index last    = 0; // index of prev estimated point
        Index i       = 1; // index of current point
        do {
            while (n_right < n) {
                // move n_left, n_right to right if radius decreases
                auto d1 = positions[i - 1] - positions[n_left];
                auto d2 = positions[n_right] - positions[i - 1];
                // if d1<=d2 with x(nright+1)==x(nright), lowest fixes
                if (d1 <= d2)
                    break;
                // radius will not decrease by move right
                n_left++;
                n_right++;
            }

            // fitted value at x(i)
            if (!lowest(positions, input, positions[i - 1], m_output[i - 1], n_left, n_right, m_residuals, iter > 1))
                m_output[i - 1] = input[i];

            // all weights zero - copy over value (all rw==0)

            if (last < i - 1) { // skipped points -- interpolate
                auto denom = positions[i - 1] - positions[last - 1]; // non-zero - proof?
                for (Index j = last + 1; j < i; j++) {
                    auto alpha      = (positions[j - 1] - positions[last - 1]) / denom;
                    m_output[j - 1] = alpha * m_output[i - 1] + (1.0 - alpha) * m_output[last - 1];
                }
            }
            last     = i; // last point actually estimated
            auto cut = positions[last - 1] + delta; // x coord of close points
            for (i = last + 1; i <= n; i = i + 1) { // find close points
                if (positions[i - 1] > cut)
                    break; // i one beyond last pt within cut
                if (positions[i - 1] == positions[last - 1]) { // exact match in x
                    m_output[i - 1] = m_output[last - 1];
                    last            = i;
                }
            }
            i = std::max(last + 1, i - 1);
            // back 1 point so interpolation within delta, but always go forward
        } while (last < n);
        m_residuals = input - m_output;
        if (iter > nSteps)
            break; // compute robustness weights except last time

        m_robustnessWeights = m_residuals.cwiseAbs();
        std::sort(m_robustnessWeights.data(), m_robustnessWeights.data() + m_robustnessWeights.size());
        auto m1   = 1 + n / 2;
        auto m2   = n - m1 + 1;
        auto cmad = 3.0 * (m_robustnessWeights[m1 - 1] + m_robustnessWeights[m2 - 1]); // 6 median abs resid
        auto c9   = 0.999 * cmad;
        auto c1   = 0.001 * cmad;
        for (Index i = 0; i < n; ++i) {
            auto r = abs(m_residuals[i]);
            if (r <= c1)
                m_robustnessWeights[i] = 1.0; // near 0, avoid underflow
            else if (r > c9)
                m_robustnessWeights[i] = 0.0; // near 1, avoid underflow
            else
                m_robustnessWeights[i] = biSquare(r / cmad);
        }
    }
}

} // end namespace anomaly::core::lowess
