#pragma once

#include <eigen3/Eigen/Dense>

namespace anomaly::core::lowess {

struct LowessConfiguration {
    using Scalar = Eigen::VectorXd::Scalar;
    using Index  = Eigen::Index;

    /*@param f */

    /**
     * @brief m_ratio specifies the amount of smoothing
     *
     * Ratio is the fraction of points used to compute each fitted value.
     * As ratio increases the smoothed values become smoother.
     * Choosing ratio in the range 0.2 to 0.8 usually results in a good fit.
     * If you have no idea which value to use, try ratio=0.5
     */
    Scalar m_ratio;

    /* @param nSteps the number of iterations in the robust fit; if nSteps=0, the nonrobust fit is returned; setting nSteps equal to 2
     * should serve most purposes.
     * */
    Index m_numberOfSteps;
    /*
     * @brief nonnegative parameter which may be used to save computations
     *
     * If the number of points is less than 100, set delta equal to 0.0.
     *
     * Very roughly the algorithm is this:
     *
     * On the initial fit and on each of the nSteps iterations locally weighted regression fitted values
     * are computed at points in X which are spaced, roughly, delta apart;
     * then the fitted values at the remaining points are computed using linear interpolation.
     *
     * The first locally weighted regression (l.w.r.) computation is carried out at position[0] and the last is carried out at
     * position[N-1].
     * Suppose the l.w.r. computation is carried out at position[I].
     * If position[I+1] is greater than or equal to position[I]+DELTA, the next l.w.r. computation is carried out at position[I+1].
     * If position[I+1] is less than position[I]+DELTA, the next l.w.r. computation is carried out at
     * the largest position[J] which is greater than or equal to position[I] but is not greater than X[I]+DELTA.
     *
     * Then the fitted values for X[K] between X[I] and X[J], if there are any, are computed by linear interpolation of the fitted values at
     * X[I] and X[J].
     *
     * If the number of points is less than 100 then delta can be set to 0.0 since the computation time will not be too great.
     * For larger numbers of points it is typically not necessary to carry out the l.w.r. computation for all points,
     * so that much computation time can be saved by taking DELTA to be greater than 0.0.
     *
     * If DELTA = Range(positions)/k then, if the values in X were uniformly scattered over the range, the full l.w.r. computation would be
     * carried out at approximately k points. Taking k to be 50 often works well.

     */
    Scalar m_delta;
};

class LowessAlgorithm {
    using Scalar = Eigen::VectorXd::Scalar;
    using Index  = Eigen::Index;

public:
    LowessAlgorithm(Index numberOfPoints)
        : m_numberOfPoints(numberOfPoints)
        , m_output(numberOfPoints)
        , m_robustnessWeights(numberOfPoints)
        , m_residuals(numberOfPoints) { }

    /**
     * @brief lowess computess the smooth of a scatterplot of input against positions using robust locally weighted regression.
     *
     * Fitted values, m_ouptut, are computed at each of the values of the horizontal axis in X.
     *
     * @param config some configuration parameters for the algorithm
     * @param positions abscissas of the points on the scatterplot; the values must be ordered from smallest to largest.
     * @param input ordinates of the points on the scatterplot
     *
     * Method
     * ======
     *
     * The fitted values are computed by using the nearest neighbor routine and robust locally weighted regression of degree 1
     * with the tricube weight function. A few additional features have been added. Suppose r is \f$f\cdot N\f$ truncated to an integer.
     * Let h be the distance to the r-th nearest neighbor from positions[I].
     * All points within h of X[I] are used. Thus if the r-th nearest neighbor is exactly the same distance as other  points,
     * more than r points can possibly be used for the smooth at positions[I].
     * There are two cases where robust locally weighted regression of degree 0 is actually used at positions[I].
     *    1. One case occurs when h is 0.0.
     *    2. The second case occurs when the weighted standard error of the positions [I] with respect to the weights w[j] is less than .001
     *       times the range of the X[I], where w[j] is the weight assigned to the j-th point of positions (the tricube weight times the
     *       robustness weight) divided by the sum of all of the weights.
     *
     * Finally, if the w[j] are all zero for the smooth at positions[I], the  fitted
     *        value is taken to be Y(I).
     */
    void lowess(const LowessConfiguration& config, const Eigen::VectorXd& positions, const Eigen::VectorXd& input);

    const Eigen::VectorXd& output() const {
        return m_output;
    }

    /**
     * Robustness weights for each point.
     * If nSteps==0, the robusness weights are not used.
     */
    const Eigen::VectorXd& robustnessWeights() const {
        return m_robustnessWeights;
    }

    /**
     * @return the residuals (=input-output).
     */
    const Eigen::VectorXd& residuals() const {
        return m_residuals;
    }

private:
    /**
     * @brief Support routine for lowess.
     *
     * The fitted value, m_output, is computed at the value, xs, of the horizontal axis.
     *
     * Robustness weights can be employed in computing the fit.
     *
     * @param positions abscissas of the points on the scatterplot; the values in X must be ordered from smallest to largest, strictly
     * increasing.
     * @param input ordinates of the points on the scatterplot.
     * @param position value of the horizontal axis at which the smooth is computed.
     * @param fittedValue output parameter, fitted value at xs
     * @param n_left index of the first point which should be considered in computing the fitted value.
     * @param n_right index of the last point which should be considered in computing the fitted value.
     * @param weights Weight for each input value used in the expression of ys, which is the sum from I = n_left to n_right  of
     * w[I]*input[I]. Only defined at locations from n_left to n_right.
     * @param use_rw If true, a robust fit is carried out using the weights in m_robutsnessWeights, else the robustness weights are ignored.
     * @return If the weights for the smooth are all 0, the fitted value, ys, is not computed and false is returned.
     *
     * Method
     * ======
     * The smooth at XS is computed using (robust) locally weighted regression of degree 1.
     * The tricube weight function is used with h equal to the maximum of position-positions[NLEFT] and position[n_right] - position.
     * Two cases where the program reverts to locally weighted regression of degree 0 are described in the documentation
     * for lowess.
     */
    bool lowest(const Eigen::VectorXd& positions,
                const Eigen::VectorXd& input,
                const Scalar           position,
                Scalar&                fittedValue,
                const Index            n_left,
                const Index            n_right,
                Eigen::VectorXd&       weights,
                const bool             use_rw);

    Index           m_numberOfPoints;
    Eigen::VectorXd m_output;
    Eigen::VectorXd m_robustnessWeights;
    Eigen::VectorXd m_residuals;
};

} // end namespace anomaly::core::lowess
