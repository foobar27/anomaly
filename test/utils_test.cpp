#include <anomaly/core/utils.hpp>
#include <eigen3/Eigen/Dense>

#include "catch.hpp"

using namespace anomaly::core::utils;
using Eigen::Index;

template <Index Size>
requires StaticOrDynamicSize<Size> struct SimpleExample {
    constexpr SimpleExample() requires StaticSize<Size> = default;

    explicit SimpleExample(Index size) requires DynamicSize<Size> : m_size(size) { }

    // Solution 1: rely on auto
    // Problem: Return type is not constexpr
    //    auto size() const {
    //        return m_size();
    //    }

    // Solution 2: dispatch via concepts (they are exclusive)
    // Problem: call to member function 'size' is ambiguous
    //    constexpr auto size() const requires StaticSize<Size> {
    //        return m_size();
    //    }
    //    auto size() const requires DynamicSize<Size> {
    //        return m_size();
    //    }

    // Solution 3: Dispatch via enable_if
    // Problem: class member cannot be redeclared.
    //    template <typename = std::enable_if<StaticOrDynamicSizeContainer<Size>::isStatic>>
    //    constexpr auto size() const {
    //        return m_size();
    //    }
    //    template <typename = std::enable_if<!StaticOrDynamicSizeContainer<Size>::isStatic>>
    //    auto size() const {
    //        return m_size();
    //    }

    // Solution 4: Like Solution 2, but with different function names
    // Problem: Shifts responsbility to caller.
    [[nodiscard]] constexpr auto staticSize() const requires StaticSize<Size> {
        return m_size();
    }

    [[nodiscard]] auto dynamicSize() const requires DynamicSize<Size> {
        return m_size();
    }

    // Solution 5: Introduce StaticOrDynamicContainer<Size>::ReturnType
    // Problem: unfortunately I cannot define the type constexpr in the static case

    // Solution 6: Return StaticOrDynamicSizeContainer
    // Problem: constexpr is lost when executing size()().
    //    [[nodiscard]] auto size() const {
    //        return m_size;
    //    }

private:
    StaticOrDynamicSizeContainer<Size> m_size;
};

// TODO(sw) simplify by overloading operator
template <int Size>
requires StaticOrDynamicSize<Size> struct DerivedExample {
    constexpr DerivedExample() requires StaticSize<Size> = default;
    explicit DerivedExample(Eigen::Index size) requires DynamicSize<Size>
        : m_originalSize(size), m_derivedSize(2 * m_originalSize() + 1), m_values(m_derivedSize()) { }

    // private: // TODO(sw) make private again
    StaticOrDynamicSizeContainer<Size>         m_originalSize;
    StaticOrDynamicSizeContainer<2 * Size + 1> m_derivedSize;

private:
    Eigen::Matrix<double, Size, 1>             m_values;
};

TEST_CASE("static or dynamic size") {
    SECTION("simple use case") {
        SECTION("static") {
            constexpr SimpleExample<5> example;
            constexpr auto             size = example.staticSize();
            REQUIRE(size == 5);
        }

        SECTION("dynamic") {
            const SimpleExample<Eigen::Dynamic> example(5);
            const auto                          size = example.dynamicSize();
            REQUIRE(size == 5);
        }
    }

    SECTION("derived use case") {
        SECTION("static") {
            const DerivedExample<5> example;
            constexpr auto          original_size = example.m_originalSize();
            constexpr auto          derived_size  = example.m_derivedSize();
            REQUIRE(original_size == 5);
            REQUIRE(derived_size == 11);
        }

        SECTION("dynamic") {
            const DerivedExample<Eigen::Dynamic> example(5);
            const auto                           original_size = example.m_originalSize();
            const auto                           derived_size  = example.m_derivedSize();
            REQUIRE(original_size == 5);
            REQUIRE(derived_size == 11);
        }
    }
}
