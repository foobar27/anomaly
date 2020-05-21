#pragma once

#include "anomaly/core/lowess.hpp"

namespace anomaly::core::stl {

struct StlConfiguration {
    using Scalar = Eigen::VectorXd::Scalar;
    using Index  = Eigen::Index;

    /**
     * @brief The period of the seasonal component.
     *
     * For example, if the time series is monthly with a yearly cycle,
     * then the period is 12.
     */
    Index m_period;

    /**
     * @brief The configuration of the seasonal smoother.
     *
     * Smoothes the seasonal component at a given point in the seasonal cycle
     * (e.g., January values of a monthly series with a yearly cycle).
     *
     * It is recommended that the length of the smoother is at least 7.
     */
    lowess::LowessConfiguration m_seasonalSmoother;

    /**
     * @brief m_trendSmoother
     *
     * It is recommended that the length of the smoother is between 1.5 * period and 2 * period.
     */
    lowess::LowessConfiguration m_trendSmoother;
    lowess::LowessConfiguration m_lowPassSmoother;

    /**
     * @brief Number of iterations for updating the seasonal and trend components.
     *
     * The value should be a positive integer.
     * See the number of outer iterations for advice on the choice of the number of inner iterations.
     */
    Index m_nIterationsInnerLoop;

    /**
     * @brief number of iterations of robust fitting.
     *
     * The value should be a non-negative integer.
     *
     * If the data are well behaved without outliers, then robustness iterations are not needed.
     * In this case set no=0, and set ni=2 to 5 depending on how much security you want that
     * the seasonal-trend looping converges.
     *
     * If outliers are present then no=3 is a very secure value unless the outliers are radical,
     * in which case the number of outer iterations 5 or even 10 might be better.
     * If the number of outer iterations is 1 or bigger, then set the number of inner loops to 1 or 2.
     */
    Index m_nIterationsOuterLoop;
};

struct StlAlgorithm {
    using Scalar = Eigen::VectorXd::Scalar;
    using Index  = Eigen::Index;

    StlAlgorithm(const StlConfiguration& config, Index numberOfPoints);

    const Eigen::VectorXd& season() const {
        return m_season;
    }

    const Eigen::VectorXd& trend() const {
        return m_trend;
    }

    void stl(const Eigen::VectorXd& input);

    // void stlez(const Eigen::VectorXd& input, long period, long ns);

private:
    void innerLoop(const Eigen::VectorXd& input, bool use_rw);

    const StlConfiguration m_config;
    Eigen::VectorXd        m_season;
    Eigen::VectorXd        m_trend;
    Eigen::VectorXd        m_robustnessWeights; // TODO(sw) move to lowess
    Eigen::VectorXd        m_tmp1;
    Eigen::VectorXd        m_tmp2;
    Eigen::VectorXd        m_tmp3;
    Eigen::VectorXd        m_tmp4;
    Eigen::VectorXd        m_tmp5;
};

} // end namespace anomaly::core::timeseries
