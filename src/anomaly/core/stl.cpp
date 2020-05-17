//#include "stl.hpp"

//#include <algorithm>
//#include <cmath>

//#include <boost/array.hpp>
//#include <eigen3/Eigen/Dense>

// using Scalar = Eigen::VectorXd::Scalar;
// using Index  = Eigen::Index;

// namespace anomaly::core::stl {

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

// inline double square(double x) {
//    return x * x;
//}

// inline double biSquare(double x) {
//    square(1.0 - square(x));
//}

// inline double cube(double x) {
//    return x * x * x;
//}

// inline double triCube(double x) {
//    return cube(1.0 - cube(x));
//}

// StlAlgorithm::StlAlgorithm(const StlConfiguration& config, const Eigen::VectorXd& input)
//    : m_config(config)
//    , m_input(input) {
//    m_config.repair(); // TODO ensure via contract
//}

//// TODO(sw) reuse the corresponding method from lowess.cpp
//// (the difference is the this one also computes the residuals)
// void computeRobustnessWeights(const Eigen::VectorXd& input, const Eigen::VectorXd& fit, Eigen::VectorXd& robustnessWeights) {
//    const auto n = input.size();

//    robustnessWeights = (input - fit).cwiseAbs();

//    auto m0 = n / 2 + 1;
//    auto m1 = n - m0 + 1;
//    std::sort(robustnessWeights.data(), robustnessWeights.data() + robustnessWeights.size());
//    auto cmad         = 3.0 * (robustnessWeights(m0) + robustnessWeights(m1)); // 6 * median abs resid
//    auto c9           = 0.999 * cmad;
//    auto c1           = 0.001 * cmad;
//    robustnessWeights = (input - fit).cwiseAbs().unaryExpr([c1, c9, cmad](double r) {
//        if (r <= c1)
//            return 1.0;
//        else if (r <= c9)
//            return biSquare(r / cmad);
//        else
//            return 0.0;
//    });
//}

// void StlAlgorithm::stl() {
//    auto            n = m_input.size();
//    Eigen::VectorXd robustnessWeights(n);
//    Eigen::VectorXd season(n);
//    Eigen::VectorXd trend(n);
//    Eigen::VectorXd work1;
//    Eigen::VectorXd work2;
//    Eigen::VectorXd work3;
//    Eigen::VectorXd work4;
//    Eigen::VectorXd work5;
//    bool            use_rw{false};

//    size_t k = 0;
//    while (true) {
//        innerLoop(use_rw, robustnessWeights, season, trend, work1, work2, work3, work4, work5);
//        k++;
//        if (k > m_config.m_nIterationsOuterLoop) {
//            break;
//        }
//        work1 = trend + season;
//        computeRobustnessWeights(m_input, work1, robustnessWeights);
//        use_rw = true;
//    }
//    if (m_config.m_nIterationsOuterLoop <= 0) {
//        robustnessWeights.setOnes();
//    }
//    // TODO where are the output components?
//}

// bool est(const Eigen::VectorXd& y,
//         const Smoother&        smoother,
//         const double           xs,
//         double&                ys,
//         const size_t           n_left,
//         const size_t           n_right,
//         Eigen::VectorXd&       w,
//         const bool             use_rw,
//         const Eigen::VectorXd& rw) {
//    auto   len    = smoother.m_length;
//    auto   degree = smoother.m_degree;
//    auto   n      = y.size();
//    double range  = double(n) - double(1);
//    double h      = std::max(xs - double(n_left), double(n_right) - xs);
//    if (len > n)
//        h += (len - n) / 2;
//    double h9 = 0.999 * h;
//    double h1 = 0.001 * h;
//    // compute weights
//    double a = 0.0;
//    for (size_t j = n_left; j <= n_right; ++j) {
//        w[j - 1] = 0.0;
//        double r = abs(double(j) - xs);
//        if (r <= h9) {
//            if (r <= h1)
//                w[j - 1] = 1.0;
//            else
//                w[j - 1] = triCube(r / h);
//            if (use_rw)
//                w[j - 1] = rw[j - 1] * w[j - 1];
//            a += w[j - 1];
//        }
//    }
//    if (a <= 0.0)
//        return false;
//    // weighted least squares
//    // make sum of w[j] == 1
//    for (size_t j = n_left; j <= n_right; ++j)
//        w[j - 1] /= a;
//    // use linear fit
//    if ((h > 0.) & (degree > 0)) {
//        // weighted center of x values
//        double a = 0.0;
//        for (size_t j = n_left; j <= n_right; ++j)
//            a += w[j - 1] * double(j);
//        double b = xs - a;
//        double c = 0.0;
//        for (size_t j = n_left; j <= n_right; ++j)
//            c += w[j - 1] * square(double(j) - a);
//        if (sqrt(c) > .001 * range) {
//            b = b / c;
//            // points are spread out enough to compute slope
//            for (size_t j = n_left; j <= n_right; ++j)
//                w[j - 1] = w[j - 1] * (b * (float(j) - a) + 1.0);
//        }
//    }
//    // Populate result
//    // TODO(sw) could be a return value?
//    ys = 0.0;
//    for (size_t j = n_left; j <= n_right; ++j)
//        ys += w[j - 1] * y[j - 1];
//    return true;
//}

// void ess(const Eigen::VectorXd& y, const Smoother& smoother, const bool use_rw, const Eigen::VectorXd& rw, Eigen::VectorXd& ys,
// Eigen::VectorXd& res) {
//    auto       len    = smoother.m_length;
//    auto       degree = smoother.m_degree;
//    auto       n_jump = smoother.m_jump;
//    const auto n      = y.size();
//    if (n < 2) {
//        ys[0] = y[0];
//        return;
//    }
//    size_t newnj = std::min(n_jump, n - 1);
//    size_t n_left, n_right;
//    if (len >= n) { // len> or = n size_t
//        n_left  = 1;
//        n_right = n;
//        for (size_t i = 1; i <= n; i += newnj)
//            if (!est(y, smoother, float(i), ys[i - 1], n_left, n_right, res, use_rw, rw))
//                ys[i - 1] = y[i - 1];
//    } else if (newnj == 1) { // newnj equal to one, len less than n
//        size_t nsh = (len + 1) / 2;
//        n_left     = 1;
//        n_right    = len;
//        for (int i = 1; i <= n; ++i) { // fitted value at i
//            if (i > nsh && n_right != n) {
//                n_left  = n_left + 1;
//                n_right = n_right + 1;
//            }
//            if (!est(y, smoother, float(i), ys[i - 1], n_left, n_right, res, use_rw, rw))
//                ys[i - 1] = y[i - 1];
//        }
//    } else { // newnj greater than one, len less than n
//        size_t n_sh = (len + 1) / 2;
//        for (int i = 1; i <= n; ++i) { // fitted value at i
//            if (i < n_sh) {
//                n_left  = 1;
//                n_right = len;
//            } else if (i >= n - n_sh + 1) {
//                n_left  = n - len + 1;
//                n_right = n;
//            } else {
//                n_left  = i - n_sh + 1;
//                n_right = len + i - n_sh;
//            }
//            if (!est(y, smoother, double(i), ys[i - 1], n_left, n_right, res, use_rw, rw))
//                ys[i - 1] = y[i - 1];
//        }
//    }
//    if (newnj != 1) {
//        for (int i = 1; i <= n - newnj; i += newnj) {
//            double delta = (ys[i + newnj - 1] - ys[i - 1]) / double(newnj);
//            for (int j = i + 1; j <= i + newnj - 1; ++i)
//                ys[j - 1] = ys[i - 1] + delta * float(j - i);
//        }
//        size_t k = ((n - 1) / newnj) * newnj + 1;
//        if (k != n) {
//            if (!est(y, smoother, double(n), ys[n - 1], n_left, n_right, res, use_rw, rw))
//                ys[n - 1] = y[n - 1];
//            if (k != n - 1) {
//                double delta = (ys(n) - ys(k)) / float(n - k);
//                for (int j = k + 1; j <= n - 1; ++j)
//                    ys[j - 1] = ys[k - 1] + delta * float(j - k);
//            }
//        }
//    }
//}

// void movingAverage(const Eigen::VectorXd& x, const size_t len, Eigen::VectorXd& ave) {
//    auto   n    = x.size();
//    auto   newn = n - len + 1;
//    double flen = double(len);
//    double v    = 0.0;
//    // get the first average
//    for (size_t i = 1; i <= len; ++i)
//        v += x[i - 1];
//    ave[0] = v / flen;
//    if (newn > 1) {
//        int k = len;
//        int m = 0;
//        for (size_t j = 2; j <= newn; ++j) {
//            // window down the array
//            k          = k + 1;
//            m          = m + 1;
//            v          = v - x[m - 1] + x[k - 1];
//            ave[j - 1] = v / flen;
//        }
//    }
//}

// void StlAlgorithm::fts(const Eigen::VectorXd& x, Eigen::VectorXd& trend, Eigen::VectorXd& work) {
//    auto n  = x.size();
//    auto np = m_config.m_period;

//    movingAverage(x, np, trend);
//    movingAverage(trend, np, work);
//    movingAverage(work, 3, trend);
//}

// void StlAlgorithm::ss(const Eigen::VectorXd& y,
//                      const bool             use_rw,
//                      const Eigen::VectorXd& rw,
//                      Eigen::VectorXd&       season,
//                      Eigen::VectorXd&       work1,
//                      Eigen::VectorXd&       work2,
//                      Eigen::VectorXd&       work3,
//                      Eigen::VectorXd&       work4) {
//    const auto n  = y.size();
//    auto       np = m_config.m_period;
//    auto       ns = m_config.m_seasonalSmoother.m_length;

//    // integer ns, isdeg, nsjump, nright, nleft, i, j, k
//    // real y(n), rw(n), season(n+2*np), work1(n), work2(n), work3(n), work4(n), xs
//    // logical userw,ok
//    for (size_t j = 1; j <= np; j = j + 1) {
//        int k = (n - j) / np + 1;
//        for (size_t i = 1; i < k; ++i)
//            work1[i - 1] = y[(i - 1) * np + j - 1];
//        if (use_rw)
//            for (size_t i = 1; i <= k; ++i)
//                work3[i - 1] = rw[(i - 1) * np + j - 1];
//        // TODO what is the k-parmeter in the following line? Do we need to provide only a segment of work1?
//        ess(work1, k, m_config.m_seasonalSmoother, use_rw, work3, work2[1], work4);
//        double xs      = 0;
//        size_t n_right = std::min(ns, k);
//        // TODO what is the k-parmeter in the following line? Do we need to provide only a segment of work1?
//        if (!est(work1, k, m_config.m_seasonalSmoother, xs, work2[0], 1, n_right, work4, use_rw, work3))
//            work2[0] = work2[1];
//        xs            = k + 1;
//        size_t n_left = std::max(1, k - ns + 1);
//        // TODO what is the k-parmeter in the following line? Do we need to provide only a segment of work1?
//        if (!est(work1, k, m_config.m_seasonalSmoother, xs, work2[k + 1], n_left, k, work4, use_rw, work3))
//            work2[k + 1] = work2[k + 1];

//        for (size_t m = 1; m <= k + 2; ++m)
//            season[(m - 1) * np + j - 1] = work2[m - 1];
//    }
//}

// void StlAlgorithm::innerLoop(bool             use_rw,
//                             Eigen::VectorXd& rw,
//                             Eigen::VectorXd& season,
//                             Eigen::VectorXd& trend,
//                             Eigen::VectorXd& work1,
//                             Eigen::VectorXd& work2,
//                             Eigen::VectorXd& work3,
//                             Eigen::VectorXd& work4,
//                             Eigen::VectorXd& work5) {
//    auto   n  = m_input.size();
//    size_t np = m_config.m_period;
//    for (size_t j = 1; j <= m_config.m_nIterationsInnerLoop; ++j) {
//        for (size_t i = 1; i <= n; ++i)
//            work1[i - 1] = m_input[i - 1] - trend[i - 1];
//        ss(work1, use_rw, rw, work2, work3, work4, work5, season);
//        fts(work2, work3, work1);
//        ess(work3, m_config.m_lowPassSmoother, false, work4, work1, work5);
//        for (size_t i = 1; i <= n; ++i)
//            season[i - 1] = work2[np + i - 1] - work1[i - 1];
//        for (size_t i = 1; i <= n; ++i)
//            work1[i - 1] = m_input[i] - season(i);
//        ess(work1, m_config.m_trendSmoother, use_rw, rw, trend, work3);
//    }
//}

//// void StlAlgoritm::StlAlgoritm::innerLoop(bool use_rw, Eigen::VectorXd& rw, Eigen::VectorXd& season, Eigen::VectorXd& trend, matrix&
//// work)
//// {
////    // TODO(sw) Fix array indices

////    auto n  = m_input.size();
////    auto ni = m_config.m_nIterationsInnerLoop;
////    auto np = m_config.m_period;

////    for (size_t j = 1; j <= ni; ++j) {
////        for (size_t i = 1; i <= n; ++i) {
////            work(i, 1) = m_input[i] - trend[i];
////        }
////        seasonalSmoothing(
////            work(1, 1), n, np, m_config.m_seasonalSmoother.m_length, isdeg, nsjump, use_rw, rw, work(1, 2), work(1, 3), work(1, 4),
////            work(1, 5), season);
////        fts(work(1, 2), n + 2 * np, np, work(1, 3), work(1, 1));
////        ess(work(1, 3), n, nl, ildeg, nljump, false, work(1, 4), work(1, 1), work(1, 5));
////        for (size_t i = 1; i <= n; ++i) {
////            season[i] = work(np + i, 2) - work(i, 1);
////        }
////        for (size_t i = 1; i <= n; ++i) {
////            work(i, 1) = y(i) - season[i];
////        }
////        ess(work(1, 1), n, nt, itdeg, ntjump, use_rw, rw, trend, work(1, 3));
////    }
////}

//// void StlAlgoritm::ess(const size_t           len,
////                      const int              degree,
////                      const size_t           njump,
////                      const bool             use_rw,
////                      const Eigen::VectorXd& rw,
////                      Eigen::VectorXd&       ys, // TODO return this, since it's the output value?
////                      const Eigen::VectorXd& res) {
////    double delta;

////    auto n = m_input.size();
////    if (n < 2) {
////        ys[0] = m_input[0];
////        return;
////    }
////    auto   newnj = std::min(njump, n - 1);
////    size_t n_left, n_right;
////    if (len >= n) {
////        n_left  = 0;
////        n_right = n - 1;
////        for (size_t i = 0; i < n; i += newnj) {
////            if (!est(len, degree, double(i + 1), ys[i], n_left, n_right, res, use_rw, rw)) {
////                ys[i] = m_input[i];
////            }
////        }
////    } else if (newnj == 1) { // newnj equal to one, len less than n
////        size_t n_sh = (len + 1) / 2;
////        n_left      = 0;
////        n_right     = len - 1;
////        for (size_t i = 0; i < n; ++i) { // fitted value at i
////            if (i + 1 > n_sh && n_right + 1 != n) {
////                n_left++;
////                n_right++;
////            }
////            if (!est(len, degree, double(i + 1), ys[i], n_left, n_right, res, use_rw, rw)) {
////                ys[i] = m_input[i];
////            }
////        }
////    } else { // newnj greater than one, len less than n
////        size_t n_sh = (len + 1) / 2;
////        for (size_t i = 0; i < n; i += newnj) { // fitted value at i
////            if (i + 1 < n_sh) {
////                n_left  = 0;
////                n_right = len - 1;
////            } else if (i + 1 >= n - n_sh + 1) {
////                n_left  = n - len;
////                n_right = n - 1;
////            } else {
////                n_left  = i + 1 - n_sh;
////                n_right = len + i - n_sh;
////            }
////            if (!est(len, degree, double(i + 1), ys[i], n_left, n_right, res, use_rw, rw)) {
////                ys[i] = m_input[i];
////            }
////        }
////    }
////    if (newnj != 1) {
////        for (size_t i = 0; i < n - newnj; i += newnj) {
////            delta = (ys[i + newnj] - ys[i]) / double(newnj);
////            for (size_t j = i + 1; j < i + newnj - 1; ++j) {
////                ys[j] += delta * double(j + 1 - i);
////            }
////        }

////        size_t k = ((n - 1) / newnj) * newnj;
////        if (k + 1 != n) {
////            if (!est(len, dewgreee, double(n), ys[n - 1], n_left, n_right, res, use_rw, rw)) {
////                ys[n - 1] = m_input[n - 1];
////            }
////            if (k + 1 != n - 1) {
////                delta = (ys[n - 1] - ys[k]) / double(n - k - 1);
////                for (size_t j = k + 1; j < n - 1; ++j) {
////                    ys[j] = ys[k] + delta * double(j - k);
////                }
////            }
////        }
////    }
////}

//// template <typename V>
//// bool makeSum1(V& v) {
////    // Make sum of the segment equal to 1.
////    double sum = v.sum();
////    if (sum <= 0.0) {
////        return false;
////    }
////    v /= sum;
////    return true;
////}

//// bool StlAlgoritm::est(const size_t           len,
////                      const size_t           ideg,
////                      const double           xs,
////                      double&                ys,
////                      const size_t           n_left,
////                      const size_t           n_right,
////                      Eigen::VectorXd&       w,
////                      const bool             use_rw, // TODO replace by std::optional<Vector&>
////                      const Eigen::VectorXd& rw) {
////    // TODO ys is actually an output argument => should we return it?
////    // TOOD what about the other arguments?

////    // TODO adjust n_left and n_right (maybe some are already adjusted)
////    size_t n = m_input.size();

////    auto range = double(n) - 1.0;
////    auto h     = std::max(xs - double(n_left + 1), double(n_right + 1) - xs);
////    if (len > n) {
////        h += (len - n) / 2;
////    }
////    double h9 = 0.999 * h;
////    double h1 = 0.001 * h;

////    // compute weights
////    auto w_seg = w.segment(n_left - 1, n_right - n_left + 1); // TODO correct segment boundaries?
////    w_seg.setZero();
////    for (size_t j = n_left - 1; j < n_right; ++j) {
////        double r = std::abs(double(j + 1) - xs);
////        if (r <= h9) {
////            if (r <= h1)
////                w[j] = 1.0;
////            else
////                w[j] = triCube(r / h);
////        }
////    }
////    if (use_rw)
////        w = rw * w;
////    if (!makeSum1(w_seg))
////        return false;

////    // weighted least squares
////    if ((h > 0.0) & (ideg > 0)) { // use linear fit
////        double a = 0.0;
////        for (size_t j = n_left - 1; j < n_right; ++j) { // weighted center of x values
////            a += w[j] * double(j + 1);
////        }
////        double b = xs - a;
////        double c = 0.0;
////        for (size_t j = n_left - 1; j < n_right; ++j) {
////            c += w[j] * square(double(j + 1) - a);
////        }
////        if (sqrt(c) > .001 * range) {
////            b = b / c;
////            // points are spread out enough to compute slope
////            for (size_t j = n_left - 1; j < n_right; ++j) {
////                w[j] *= (b * (double(j + 1) - a) + 1.0);
////            }
////        }
////    }
////    // TODO replace by eigen VectorXd.dot()
////    auto input_seg = m_input.segment(n_left - 1, n_right - n_left + 1);
////    ys             = w_seg.dot(input_seg);
////    return true;
////}

//// static void movingAverage(const Eigen::VectorXd& x, const size_t len, Eigen::VectorXd& output) {
////    size_t n     = x.size();
////    size_t new_n = n - len + 1;
////    double ilen  = 1.0 / len;
////    // get the first average
////    double v  = x.segment(0, len).sum();
////    output[0] = v * ilen;

////    if (new_n > 1) {
////        size_t k = len;
////        size_t m = 0;
////        for (size_t j = 1; j < new_n; ++j) {
////            // window down the array
////            k         = k + 1;
////            m         = m + 1;
////            v         = v - x[m] + x[k];
////            output[j] = v * ilen;
////        }
////    }
////}

//// void StlAlgoritm::fts(const Eigen::VectorXd& x, Eigen::VectorXd& trend, Eigen::VectorXd& work) {
////    auto np = m_config.m_period;
////    movingAverage(x, np, trend);
////    movingAverage(trend, np, work);
////    movingAverage(work, 3, trend);
////}

//// void StlAlgoritm::innerLoop(const bool use_rw, Eigen::VectorXd& rw, Eigen::VectorXd& season, Eigen::VectorXd& trend,
//// StlAlgoritm::matrix& work) {
////    // TODO(sw) Fix array indices
////    auto n  = m_input.size();
////    auto ni = m_config.m_nIterationsInnerLoop;
////    auto np = m_config.m_period;

////    for (size_t j = 1; j <= ni; ++j) {
////        // Step 1: Detrending
////        // TODO replace by work1 = m_input - trend;
////        for (size_t i = 1; i <= n; ++i) {
////            work(i, 1) = m_input[i] - trend[i];
////        }
////        // Cycle-subseries Smoothing
////        seasonalSmoothing(
////            work(1, 1), n, np, m_config.m_seasonalSmoother.m_length, isdeg, nsjump, use_rw, rw, work(1, 2), work(1, 3), work(1, 4),
////            work(1, 5), season);
////        // Low-pass Filtering of smoothed cycle-subseries
////        fts(work(1, 2), n + 2 * np, np, work(1, 3), work(1, 1));
////        ess(work(1, 3), n, nl, ildeg, nljump, false, work(1, 4), work(1, 1), work(1, 5));
////        // Step 6: Detrending of the smoothed cycle-subseries
////        for (size_t i = 1; i <= n; ++i) {
////            season[i] = work(np + i, 2) - work(i, 1);
////        }
////        // Step 5: Deseasonalizing
////        for (size_t i = 1; i <= n; ++i) {
////            work(i, 1) = y(i) - season[i];
////        }
////        // Step 6: Trend smoothing
////        ess(work(1, 1), n, nt, itdeg, ntjump, use_rw, rw, trend, work(1, 3));
////    }
////}

//// void StlAlgoritm::rwts(Eigen::VectorXd& fit, Eigen::VectorXd& rw) {
////    // TODO(sw) Fix array indices
////    size_t          n = m_input.size();
////    Eigen::Vector2d mid;

////    // integer n
////    // real y(n), fit(n), rw(n), cmad, c9, c1, r

////    rw     = (m_input - fit).cwiseAbs();
////    mid[0] = n / 2 + 1;
////    mid[1] = n - mid[0] + 1;
////    psort(rw, mid);
////    double cmad = 3.0 * (rw[mid[1]] + rw[mid[2]]); // 6 * median abs resid
////    double c9   = 0.999 * cmad;
////    double c1   = 0.001 * cmad;
////    for (size_t i = 1; i <= n; ++i) {
////        double r = fabs(m_input[i] - fit[i]);
////        if (r <= c1)
////            rw[i] = 1.0;
////        else if (r <= c9)
////            rw[i] = square(1.0 - square(r / cmad));
////        else
////            rw[i] = 0.0;
////    }
////}

//// void StlAlgoritm::seasonalSmoothing() {
////    // subroutine ss(use_rw,rw,season,work1,work2,work3,work4)
////    // TODO(sw) Fix array indices

////    // integer n, np, ns, isdeg, nsjump, nright, nleft, i, j, k
////    // real y(n), rw(n), season(n+2*np), xs
////    // logical use_rw,ok
////    auto            n  = m_input.size();
////    auto            np = m_config.m_period;
////    auto            ns = m_config.m_seasonalSmoother.m_length;
////    Eigen::VectorXd work1(n);
////    Eigen::VectorXd work2(n);
////    Eigen::VectorXd work3(n);
////    Eigen::VectorXd work4(n);

////    for (size_t j = 1; j <= np; j = j + 1) {
////        size_t k = (n - j) / np + 1;
////        for (size_t i = 1; i <= k; ++i) {
////            work1[i] = m_input[(i - 1) * np + j];
////        }
////        if (use_rw)
////            for (size_t i = 1; i < k; ++i) {
////                work3[i] = rw[(i - 1) * np + j];
////            }
////        ess(work1, k, ns, isdeg, nsjump, use_rw, work3, work2[2], work4);
////        double xs      = 0.0;
////        int    n_right = std::min(ns, k);
////        if (!est(work1, k, ns, isdeg, xs, work2[1], 1, n_right, work4, use_rw, work3)) {
////            work2[1] = work2[2];
////        }
////        xs         = k + 1;
////        int n_left = std::max(1, k - ns + 1);
////        if (!est(work1, k, ns, isdeg, xs, work2[k + 2], n_left, k, work4, use_rw, work3)) {
////            work2[k + 2] = work2[k + 1];
////        }
////        for (size_t m = 1; m <= k + 2; ++m) {
////            season[(m - 1) * np + j] = work2[m];
////        }
////    }
////}

////// TODO(sw) should be a configuration builder?
////// subroutine stlez(y, n, np, ns, isdeg, itdeg, robust, no, rw, season, trend, work)

////// logical robust
////// integer n, i, j, np, ns, no, nt, nl, ni, nsjump, ntjump, nljump, newns, newnp
////// integer isdeg, itdeg, ildeg
////// real y(n), rw(n), season(n), trend(n), work(n+2*np,7)
////// real maxs, mins, maxt, mint, maxds, maxdt, difs, dift

////// ildeg = itdeg
////// newns = max0(3,ns)
////// if(mod(newns,2)==0) newns = newns+1
////// newnp = max0(2,np)
////// nt = (1.5*newnp)/(1 - 1.5/newns) + 0.5
////// nt = max0(3,nt)
////// if(mod(nt,2)==0) nt = nt+1
////// nl = newnp
////// if(mod(nl,2)==0) nl = nl+1
////// if(robust) ni = 1
////// else ni = 2
////// nsjump = max0(1,int(float(newns)/10 + 0.9))
////// ntjump = max0(1,int(float(nt)/10 + 0.9))
////// nljump = max0(1,int(float(nl)/10 + 0.9))
////// do i = 1,n
//////	trend(i) = 0.0
////// call onestp(y,n,newnp,newns,nt,nl,isdeg,itdeg,ildeg,nsjump,ntjump,nljump,ni,
//////	.false.,rw,season,trend,work)
////// no = 0
////// if(robust){
//////	for(j=1; j<=15; j=j+1){	#robustness iterations
//////		do i = 1,n{	#initialize for testing
//////			work(i,6) = season(i)
//////			work(i,7) = trend(i)
//////			work(i,1) = trend(i)+season(i)
//////			}
//////		call rwts(y,n,work(1,1),rw)
//////		call onestp(y, n, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump,
//////			ntjump, nljump, ni, .true., rw, season, trend, work)
//////		no = no+1
//////		maxs = work(1,6)
//////		mins = work(1,6)
//////		maxt = work(1,7)
//////		mint = work(1,7)
//////		maxds = abs(work(1,6) - season(1))
//////		maxdt = abs(work(1,7) - trend(1))
//////		do i = 2,n{
//////			if(maxs<work(i,6)) maxs = work(i,6)
//////			if(maxt<work(i,7)) maxt = work(i,7)
//////			if(mins>work(i,6)) mins = work(i,6)
//////			if(mint>work(i,7)) mint = work(i,7)
//////			difs = abs(work(i,6) - season(i))
//////			dift = abs(work(i,7) - trend(i))
//////			if(maxds<difs) maxds = difs
//////			if(maxdt<dift) maxdt = dift
//////			}
//////		if((maxds/(maxs-mins)<.01) & (maxdt/(maxt-mint)<.01)) break
//////		}
//////	}
////// if(!robust){
//////	do i = 1,n
//////		rw(i) = 1.0
//////	}

//} // end namespace anomaly::core::timeseries
