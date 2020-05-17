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

    Index j = n_left; // The loop might terminate earlier, we need to remember the last loop position.
    for (; j < n; j++) { // compute weights (pick up all ties on right)
        weights[j] = 0.0;
        auto r     = abs(positions[j] - position);
        if (r <= h9) { // small enough for non-zero weight
            if (r > h1)
                weights[j] = triCube(r / h);
            else
                weights[j] = 1.0;
        } else if (positions[j] > position)
            break; // get out at first zero wt on right
    }
    // rightmost pt (may be greater than n_right because of ties)
    auto weights_seg            = weights.segment(n_left, j - n_left);
    auto robustness_weights_seg = m_robustnessWeights.segment(n_left, j - n_left);
    auto input_seg              = input.segment(n_left, j - n_left);
    auto positions_seg          = positions.segment(n_left, j - n_left);

    if (use_rw)
        weights_seg.array() *= robustness_weights_seg.array();

    // make sum of w(j) == 1 (if possible)
    auto sum_of_weights = weights_seg.sum();
    if (sum_of_weights <= 0.0)
        return false;
    weights_seg /= sum_of_weights;

    // weighted least squares
    if (h > 0.0) { // use linear fit
        double a = weights_seg.dot(positions_seg); // weighted center of x values
        auto   b = position - a;
        double c = weights_seg.dot((positions_seg.array() - a).square().matrix());
        if (sqrt(c) > .001 * range) {
            // points are spread out enough to compute slope
            b /= c;
            weights_seg.array() *= 1.0 + b * (positions_seg.array() - a);
        }
    }

    // Compute output
    fittedValue = weights_seg.dot(input_seg);
    return true;
}

static void convertResidualsToRobustnessWeights(const Eigen::VectorXd& residuals, Eigen::VectorXd& robustnessWeights) {
    auto n            = residuals.size();
    robustnessWeights = residuals.cwiseAbs();
    std::sort(robustnessWeights.data(), robustnessWeights.data() + robustnessWeights.size());
    auto m1           = n / 2;
    auto m2           = n - m1;
    auto cmad         = 3.0 * (robustnessWeights[m1] + robustnessWeights[m2]); // 6 median abs resid
    auto c9           = 0.999 * cmad;
    auto c1           = 0.001 * cmad;
    robustnessWeights = residuals.cwiseAbs().unaryExpr([c1, c9, cmad](Scalar r) {
        if (r <= c1)
            return 1.0; // near 0, avoid underflow
        else if (r > c9)
            return 0.0; // near 1, avoid underflow
        else
            return biSquare(r / cmad);
    });
}

// Interpolates the given segment from the first and last values (which are fixed).
static void interpolate(Eigen::VectorXd::ConstSegmentReturnType positions, Eigen::VectorXd::SegmentReturnType values) {
    assert(positions.size() == values.size());
    assert(positions.size() > 1);
    auto n           = positions.size();
    auto denom       = positions[n - 1] - positions[0];
    auto left_value  = values[0];
    auto right_value = values[n - 1];
    for (Index i = 1; i < n - 1; ++i) {
        auto alpha = (positions[i] - positions[0]) / denom;
        values[i]  = alpha * right_value + (1.0 - alpha) * left_value;
    }
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
    Index ns = std::clamp(f * n, 2.0, (double)n);

    // robustness iterations
    for (Index iter = 1; iter <= nSteps + 1; iter = iter + 1) {
        Index n_left  = 0;
        Index n_right = ns;
        Index last    = -1; // index of prev estimated point
        Index i       = 0; // index of current point
        do {
            while (n_right < n) {
                // move n_left, n_right to right if radius decreases
                auto d1 = positions[i] - positions[n_left];
                auto d2 = positions[n_right] - positions[i];
                // if d1<=d2 with position[nright]==position[nright-1], lowest fixes
                if (d1 <= d2)
                    break;
                // radius will not decrease by move right
                n_left++;
                n_right++;
            }

            // fitted value at x[i]
            if (!lowest(positions, input, positions[i], m_output[i], n_left, n_right, m_residuals, iter > 1))
                m_output[i] = input[i];

            // all weights zero - copy over value (all robustnessWeights==0)

            if (last + 1 < i) { // skipped points -- interpolate
                interpolate(positions.segment(last, i - last + 1), m_output.segment(last, i - last + 1));
            }
            last     = i; // last point actually estimated
            auto cut = positions[last] + delta; // x coord of close points
            for (i = last + 1; i < n; i++) { // find close points
                if (positions[i] > cut)
                    break; // i one beyond last pt within cut
                if (positions[i] == positions[last]) { // exact match in x
                    m_output[i] = m_output[last];
                    last        = i;
                }
            }
            i = std::max(last + 1, i - 1);
            // back 1 point so interpolation within delta, but always go forward
        } while (last + 1 < n);
        m_residuals = input - m_output;
        if (iter > nSteps)
            break; // compute robustness weights except last time

        convertResidualsToRobustnessWeights(m_residuals, m_robustnessWeights);
    }
}

} // end namespace anomaly::core::lowess
