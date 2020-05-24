#include "stl.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/framework/extractor.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/array.hpp>
#include <eigen3/Eigen/Dense>

#include "utils.hpp"

using Scalar = Eigen::VectorXd::Scalar;
using Index  = Eigen::Index;

namespace anomaly::core::stl {

// void Smoother::repair() {
//    m_length = std::max(size_t(3), m_length);
//    if (m_length % 2 == 0) {
//        m_length++;
//    }
//}

// void StlConfiguration::repair() {
//    // Apply boundaries to parameters
//    m_seasonalSmoother.repair();
//    m_trendSmoother.repair();
//    m_lowPassSmoother.repair();
//    m_period = std::max(size_t(2), m_period);
//}

static constexpr auto repairPeriod(Index period) {
    return std::max(Index(2), period);
}

static constexpr auto tmpSize(const StlConfiguration& config, Index numberOfPoints) {
    return numberOfPoints + 2 * repairPeriod(config.m_period);
}

static constexpr lowess::LowessConfiguration repair(lowess::LowessConfiguration config) {
    config.m_length = std::max(Index(3), config.m_length);
    if (config.m_length % 2 == 0)
        config.m_length++; // make odd
    return config;
}

static constexpr StlConfiguration repair(StlConfiguration config) {
    repair(config.m_trendSmoother);
    repair(config.m_lowPassSmoother);
    repair(config.m_trendSmoother);
    return config;
}

StlAlgorithm::StlAlgorithm(const StlConfiguration& config, Index numberOfPoints)
    : m_config(repair(config))
    , m_season(numberOfPoints)
    , m_trend(numberOfPoints)
    , m_tmp1(tmpSize(config, numberOfPoints))
    , m_tmp2(tmpSize(config, numberOfPoints))
    , m_tmp3(tmpSize(config, numberOfPoints))
    , m_tmp4(tmpSize(config, numberOfPoints))
    , m_tmp5(tmpSize(config, numberOfPoints)) { }

static void
computeRobustnessWeights(const Eigen::VectorXd& input, const Eigen::VectorBlock<Eigen::VectorXd> fit, Eigen::VectorXd& robustnessWeights) {
    using utils::biSquare;
    const auto n = input.size();

    robustnessWeights = (input - fit).cwiseAbs();

    auto m0 = n / 2 + 1;
    auto m1 = n - m0 + 1;
    std::sort(robustnessWeights.data(), robustnessWeights.data() + robustnessWeights.size());
    auto cmad         = 3.0 * (robustnessWeights(m0) + robustnessWeights(m1)); // 6 * median abs resid
    auto c9           = 0.999 * cmad;
    auto c1           = 0.001 * cmad;
    robustnessWeights = (input - fit).cwiseAbs().unaryExpr([c1, c9, cmad](double r) {
        if (r <= c1)
            return 1.0;
        if (r <= c9)
            return biSquare(r / cmad);
        return 0.0;
    });
}

void StlAlgorithm::stl(const Eigen::VectorXd& input) {
    bool use_rw = false;
    m_trend.array().setZero();
    Index k = 0;
    do {
        innerLoop(input, use_rw);
        ++k;
        if (k > m_config.m_nIterationsOuterLoop)
            break;
        m_tmp1.segment(0, m_trend.size()) = m_trend + m_season;
        computeRobustnessWeights(input, m_tmp1.segment(0, input.size()), m_robustnessWeights);
        use_rw = true;
    } while (true);
    if (m_config.m_nIterationsOuterLoop <= 0)
        m_robustnessWeights.array().setOnes();
}

// void StlAlgorithm::stlez(const Eigen::VectorXd& input, long period, Index season_smoother_length, Index season_degree, Index
// trend_degree, bool robust, Index numberOfOuterIterations) {
//    auto y = input;
//    auto n = input.size();
//    auto np = period;
//    // subroutine stlez(y, n, np, ns, isdeg, itdeg, robust, no, rw, season, trend, work)

//    // integer n, i, j, np, ns, no, nt, nl, ni, nsjump, ntjump, nljump, newns, newnp
//    // real maxs, mins, maxt, mint, maxds, maxdt, difs, dift

//    auto low_pass_degree = trend_degree;
//    auto newns = std::max(Index(3), season_smoother_length);
////    if(newns% 2 == 0) newns = newns+1;
//    // newnp = max0(2,np)
//    // nt = (1.5*newnp)/(1 - 1.5/newns) + 0.5
//    // nt = max0(3,nt)
//    // if(mod(nt,2)==0) nt = nt+1
//    // nl = newnp
//    // if(mod(nl,2)==0) nl = nl+1
//    // if(robust) ni = 1
//    // else ni = 2
//    // nsjump = max0(1,int(float(newns)/10 + 0.9))
//    // ntjump = max0(1,int(float(nt)/10 + 0.9))
//    // nljump = max0(1,int(float(nl)/10 + 0.9))
//    // do i = 1,n
//    //	trend(i) = 0.0
//    // call onestp(y,n,newnp,newns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,nljump,ni,
//    //	.false.,rw,season,trend,work)
//    // no = 0
//    // if(robust){
//    //	for(j=1; j<=15; j=j+1){	#robustness iterations
//    //		do i = 1,n{	#initialize for testing
//    //			work(i,6) = season(i)
//    //			work(i,7) = trend(i)
//    //			work(i,1) = trend(i)+season(i)
//    //			}
//    //		call rwts(y,n,work(1,1),rw)
//    //		call onestp(y, n, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump,
//    //			ntjump, nljump, ni, .true., rw, season, trend, work)
//    //		no = no+1
//    //		maxs = work(1,6)
//    //		mins = work(1,6)
//    //		maxt = work(1,7)
//    //		mint = work(1,7)
//    //		maxds = abs(work(1,6) - season(1))
//    //		maxdt = abs(work(1,7) - trend(1))
//    //		do i = 2,n{
//    //			if(maxs<work(i,6)) maxs = work(i,6)
//    //			if(maxt<work(i,7)) maxt = work(i,7)
//    //			if(mins>work(i,6)) mins = work(i,6)
//    //			if(mint>work(i,7)) mint = work(i,7)
//    //			difs = abs(work(i,6) - season(i))
//    //			dift = abs(work(i,7) - trend(i))
//    //			if(maxds<difs) maxds = difs
//    //			if(maxdt<dift) maxdt = dift
//    //			}
//    //		if((maxds/(maxs-mins)<.01) & (maxdt/(maxt-mint)<.01)) break
//    //		}
//    //	}
//    // if(!robust){
//    //	do i = 1,n
//    //		rw(i) = 1.0
//    //	}
//}

template <typename T>
static bool est(const T&         y, // size: n
                Index            len,
                Index            degree,
                double           xs,
                double&          ys,
                Index            n_left,
                Index            n_right, // TODO(sw) shift
                Eigen::VectorXd& weights, // size: n
                bool             use_rw,
                Eigen::VectorXd& robustness_weights) { // size: n
    using utils::square;
    using utils::triCube;
    auto n     = y.size();
    auto range = double(n - 1);
    auto h     = std::max(xs - double(n_left + 1), double(n_right) - xs);
    if (len > n)
        h += double((len - n) / 2);
    double h9 = .999 * h;
    double h1 = .001 * h;
    // compute weights
    auto weights_seg = weights.segment(n_left, n_right - n_left);
    for (Index j = n_left; j < n_right; ++j) {
        weights[j] = 0.0;
        auto r     = abs(double(j + 1) - xs);
        if (r <= h9) {
            if (r <= h1)
                weights[j] = 1.0;
            else
                weights[j] = triCube(r / h);
        }
    }
    if (use_rw)
        weights_seg.array() *= robustness_weights.segment(n_left, n_right - n_left).array();
    // make sum of w[j] == 1
    double a = weights_seg.sum();
    if (a <= 0.0)
        return false;
    weights_seg /= a;
    // weighted least squares
    if ((h > 0.) & (degree > 0)) { // use linear fit
        // weighted center of x values
        auto   positions = Eigen::VectorXd::LinSpaced(n_right - n_left, n_left + 1, n_right);
        double a         = weights_seg.dot(positions);
        //                for (Index j = n_left; j < n_right; ++j)
        //                    a += weights[j] * double(j + 1);
        auto   b = xs - a;
        double c = 0.0;
        for (Index j = n_left; j < n_right; ++j)
            c += weights[j] * utils::square(double(j + 1) - a);
        if (sqrt(c) > .001 * range) {
            b /= c;
            // points are spread out enough to compute slope
            for (Index j = n_left; j < n_right; ++j)
                weights[j] *= (b * (double(j + 1) - a) + 1.0);
        }
    }
    ys = weights_seg.dot(y.segment(n_left, n_right - n_left));
    return true;
}

template <typename Derived1, typename Derived2>
static void ess(const Eigen::MatrixBase<Derived1>& y, // size: n
                const lowess::LowessConfiguration& config,
                bool                               use_rw,
                Eigen::VectorXd&                   robustness_weights,
                Eigen::MatrixBase<Derived2>&       ys,
                Eigen::VectorXd&                   residuals) {

    auto n = y.size();
    if (n < 2) {
        ys[0] = y[0];
        return;
    }

    auto  new_nj = std::min(Index(config.m_delta), n - 1);
    Index n_left{};
    Index n_right{};
    if (config.m_length >= n) {
        n_left  = 0;
        n_right = n;
        for (Index i = 0; i < n; i += new_nj) {
            if (!est(y, config.m_length, config.m_degree, double(i + 1), ys[i], n_left, n_right, residuals, use_rw, robustness_weights))
                ys[i] = y[i];
        }
    } else if (new_nj == 1) { // newnj equal to one, len less than n
        Index n_sh = (config.m_length + 1) / 2;
        n_left     = 0;
        n_right    = config.m_length;
        for (Index i = 0; i < n; ++i) { // fitted value at i
            if (i + 1 > n_sh && n_right != n) {
                n_left++;
                n_right++;
            }
            if (!est(y, config.m_length, config.m_degree, double(i + 1), ys[i], n_left, n_right, residuals, use_rw, robustness_weights))
                ys[i] = y[i];
        }
    } else { // newnj greater than one, len less than n
        Index n_sh = (config.m_length + 1) / 2;
        for (Index i = 0; i < n; i += new_nj) { // fitted value at i
            if (i + 1 < n_sh) {
                n_left  = 0;
                n_right = config.m_length;
            } else if (i > n - n_sh + 1) {
                n_left  = n - config.m_length;
                n_right = n;
            } else {
                n_left  = i - n_sh + 1;
                n_right = config.m_length + i + 1 - n_sh;
            }
            if (!est(y, config.m_length, config.m_degree, double(i + 1), ys[i], n_left, n_right, residuals, use_rw, robustness_weights))
                ys[i] = y[i];
        }
    }
    if (new_nj != 1) {
        // Do the interpolation
        for (Index i = 0; i < n - new_nj; i += new_nj) {
            auto delta = (ys[i + new_nj] - ys[i]) / double(new_nj);
            for (Index j = i + 1; j < i + 1 + new_nj - 1; ++j)
                ys[j] = ys[i] + delta * double(j - i);
        }
        auto k = ((n - 1) / new_nj) * new_nj + 1;
        if (k != n) {
            if (!est(y, config.m_length, config.m_degree, double(n), ys[n - 1], n_left, n_right, residuals, use_rw, robustness_weights))
                ys[n - 1] = y[n - 1];
            if (k != n - 1) {
                auto delta = (ys[n - 1] - ys[k - 1]) / double(n - k);
                for (Index j = k; j < n - 1; ++j)
                    ys[j] = ys[k - 1] + delta * double(j + 1 - k);
            }
        }
    }
}

// If input has size n, output has size n-len+1.
static void movingAverage(const Eigen::VectorXd& input, Index len, Eigen::VectorXd& output) {
    using namespace std;
    using namespace boost::accumulators;
    accumulator_set<double, stats<tag::rolling_mean>> acc(tag::rolling_window::window_size = len);
    for_each(input.data(), input.data() + len, acc);

    // get the first average
    output[0] = rolling_mean(acc);

    auto new_n = input.size() - len + 1;
    if (new_n > 1) {
        auto k = len - 1;
        for (Index j = 1; j < new_n; ++j) {
            k++;
            acc(input[k]);
            output[j] = rolling_mean(acc);
        }
    }
}

static void fts(const Eigen::VectorXd& input, Index period, Eigen::VectorXd& trend, Eigen::VectorXd& tmp) {
    // Three subsequent moving averages: (N=original_n, P=period)
    // - input has length: n = N + 2*P
    // - smoothing with len=period yields a vector of length: n - P + 1 = N + 2*P - P + 1 = N + P + 1
    // - smoothing with len=period yields a vector of length: N + P + 1 - P  + 1 = N + 2
    // - smoothing with len=3 yields a vector of length: N + 2 - 3 + 1 = N (which is the original size)
    auto n = input.size();
    movingAverage(input, period, trend);
    movingAverage(trend.segment(0, n - period + 1), period, tmp);
    movingAverage(tmp.segment(0, n - 2 * period + 2), 3, trend);
}

template <typename DerivedInput>
static void ss(const Eigen::MatrixBase<DerivedInput>& input, // size: n
               Index                                  period,
               const lowess::LowessConfiguration&     smoother,
               bool                                   use_rw,
               Eigen::VectorXd&                       robustness_weights, // size: n
               Eigen::VectorXd&                       season, // size: n + 2*period
               Eigen::VectorXd&                       work1, // size: n
               Eigen::VectorXd&                       work2, // size: n
               Eigen::VectorXd&                       work3, // size: n
               Eigen::VectorXd&                       work4) { // size: n
    auto n = input.size();
    for (Index j = 0; j < period; ++j) {
        auto k = (n - j - 1) / period + 1;

        // Transpose cycle matrix into work1
        for (Index i = 0; i < k; ++i)
            work1[i] = input[i * period + j];

        if (use_rw) {
            // Transpose robustnessWeight series into work3
            for (Index i = 0; i < k; ++i)
                work3[i] = robustness_weights[i * period + j];
        }
        // Lowess smoothing, from work 1 to work2
        // Robustness_weights in work3 (with skips).
        auto work1_seg = work1.segment(0, k);
        auto work2_seg = work2.segment(1, work2.size() - 1);
        ess(work1_seg, smoother, use_rw, work3, work2_seg, work4);
        {
            // Add left margin to smoothing
            double xs      = 0;
            Index  n_right = std::min(smoother.m_length, k);
            if (!est(work1.segment(0, k), smoother.m_length, smoother.m_degree, xs, work2[0], 1, n_right, work4, use_rw, work3))
                work2[0] = work2[1];
        }

        {
            // Add right margin to smoothing
            auto  xs     = double(k + 1);
            Index n_left = std::max(Index(1), k - smoother.m_length + 1);
            if (!est(work1.segment(0, k), smoother.m_length, smoother.m_degree, xs, work2[k + 1], n_left, k, work4, use_rw, work3))
                work2[k + 2] = work2[k + 1];
        }

        // Transpose work2 into season (backwards)
        for (Index m = 0; m < k + 2; ++m)
            season[m * period + j] = work2[m];
    }
}

void StlAlgorithm::innerLoop(const Eigen::VectorXd& input, bool use_rw) {
    auto n        = input.size();
    auto tmp1_seg = m_tmp1.segment(0, n);
    for (Index j = 0; j < m_config.m_nIterationsInnerLoop; ++j) {
        // Step 1: Detrending
        tmp1_seg = input - m_trend;

        // Step 2: Cycle sub-series smoothing
        ss(tmp1_seg, m_config.m_period, m_config.m_seasonalSmoother, use_rw, m_robustnessWeights, m_tmp2, m_tmp3, m_tmp4, m_tmp5, m_season);

        // Step 3: Low-pass filtering of smoothed cycle series
        fts(m_tmp2, m_config.m_period, m_tmp3, m_tmp1);
        auto tmp3_seg = m_tmp3.segment(0, n);
        ess(tmp3_seg, m_config.m_lowPassSmoother, false, m_tmp4, m_tmp1, m_tmp5);

        // Step 4: Detrending of smoothed cycle sub-series
        m_season = m_tmp2.segment(m_config.m_period, n) - m_tmp1.segment(0, n);

        // Step 5: Deseasonalizing
        m_tmp1.segment(0, n) = input - m_season;

        // Step 6: Trend smoothing
        ess(tmp1_seg, m_config.m_trendSmoother, use_rw, m_robustnessWeights, m_trend, m_tmp3);
    }
}

} // end namespace anomaly::core::stl
