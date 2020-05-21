#include <eigen3/Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

static void fooVector(VectorXd& vector) {
    vector[0] = 42.0;
}

static void fooSegment(VectorBlock<VectorXd> segment) {
    segment[0] = 42.0;
}

// Recommendation from: https://eigen.tuxfamily.org/dox/TopicFunctionTakingEigenTypes.html

template <typename Derived>
static void fooGeneric(MatrixBase<Derived>& vector) {
    vector[0] = 42.0;
}

int main() {
    {
        VectorXd v(5);
        v << 1, 2, 3, 4, 5;
        fooVector(v);
        cout << "fooVector: " << v[0] << endl;
    }
    {
        VectorXd v(5);
        v << 1, 2, 3, 4, 5;
        fooSegment(v.segment(0, 1));
        cout << "fooSegment: " << v[0] << endl;
    }
    {
        VectorXd v(5);
        v << 1, 2, 3, 4, 5;
        fooGeneric(v);
        cout << "fooGeneric<Vector>: " << v[0] << endl;
    }
    {
        VectorXd v(5);
        v << 1, 2, 3, 4, 5;
        Eigen::DenseBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>::SegmentReturnType seg = v.segment(0, 1);
        fooGeneric(seg);
        cout << "fooGeneric<Segment>: " << v[0] << endl;
    }
    return 0;
}
