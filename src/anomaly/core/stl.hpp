#pragma once

#include <eigen3/Eigen/Dense>

namespace anomaly::core::stl {

struct Smoother {

    void repair();

    /**
     * @brief The length of the smoother.
     *
     * The value should be an odd integer greater than or equal to 3.
     *
     * As the value increases the values of the times series become smoother
     */
    size_t m_length;

    /**
     * @brief Degree of locally-fitted polynomial in seasonal smoothing.
     *
     * The value can be 0 or 1.
     */
    size_t m_degree;

    /*
     * @brief Skipping value for seasonal smoothing
     *
     * The smoother skips ahead this amount of points and then linearly interpolates in between.
     *
     * The value should be a positive integer.
     * If the value if 1, a seasonal smooth is calculated at all n points.
     * To make the procedure run faster, a reasonable choice is
     * 10%-20% of m_seasonalSmoother.
     */
    size_t m_jump;
};

struct StlConfiguration {

    // TODO ensure via ctor and setters in configuration
    void repair();

    /**
     * @brief The period of the seasonal component.
     *
     * For example, if the time series is monthly with a yearly cycle,
     * then the period is 12.
     */
    size_t m_period;

    /**
     * @brief The configuration of the seasonal smoother.
     *
     * Smoothes the seasonal component at a given point in the seasonal cycle
     * (e.g., January values of a monthly series with a yearly cycle).
     *
     * It is recommended that the length of the smoother is at least 7.
     */
    Smoother m_seasonalSmoother;

    /**
     * @brief m_trendSmoother
     *
     * It is recommended that the length of the smoother is between 1.5 * period and 2 * period.
     */
    Smoother m_trendSmoother;
    Smoother m_lowPassSmoother;

    /**
     * @brief Number of iterations for updating the seasonal and trend components.
     *
     * The value should be a positive integer.
     * See the number of outer iterations for advice on the choice of the number of inner iterations.
     */
    size_t m_nIterationsInnerLoop;

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
    size_t m_nIterationsOuterLoop;
};

struct StlAlgorithm {
    StlAlgorithm(const StlConfiguration& config, const Eigen::VectorXd& input);

    void stl();

private:
    using matrix = boost::numeric::ublas::matrix<double>;

    void innerLoop(bool             use_rw,
                   Eigen::VectorXd& rw,
                   Eigen::VectorXd& season,
                   Eigen::VectorXd& trend,
                   Eigen::VectorXd& work1,
                   Eigen::VectorXd& work2,
                   Eigen::VectorXd& work3,
                   Eigen::VectorXd& work4,
                   Eigen::VectorXd& work5);

    void fts(const Eigen::VectorXd& x, Eigen::VectorXd& trend, Eigen::VectorXd& work);

    void ss(const Eigen::VectorXd& y,
            const bool             use_rw,
            const Eigen::VectorXd& rw,
            Eigen::VectorXd&       season,
            Eigen::VectorXd&       work1,
            Eigen::VectorXd&       work2,
            Eigen::VectorXd&       work3,
            Eigen::VectorXd&       work4);

    const Eigen::VectorXd&  m_input;
    const StlConfiguration& m_config;
};

} // end namespace anomaly::core::timeseries
