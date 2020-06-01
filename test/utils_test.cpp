#include <anomaly/core/utils.hpp>
#include <eigen3/Eigen/Dense>

#include "catch.hpp"

using namespace anomaly::core::utils;
using Eigen::Index;

template <Index Size>
requires StaticOrDynamicSize<Size> struct SimpleExample {
    constexpr SimpleExample() requires StaticSize<Size> = default;

    explicit SimpleExample(Index size) requires DynamicSize<Size> : m_size(size) { }

    [[nodiscard]] constexpr auto size() const noexcept {
        return m_size();
    }

private:
    StaticOrDynamicSizeContainer<Size> m_size;
};

// TODO(sw) simplify by overloading operator
template <int Size>
requires StaticOrDynamicSize<Size> struct DerivedExample {
    constexpr DerivedExample() requires StaticSize<Size> = default;
    explicit DerivedExample(Eigen::Index size) requires DynamicSize<Size>
        : m_originalSize(size), m_derivedSize(2 * m_originalSize() + 1), m_values(m_derivedSize()) { }

    [[nodiscard]] constexpr auto originalSize() const noexcept {
        return m_originalSize();
    }

    [[nodiscard]] constexpr auto derivedSize() const noexcept {
        return m_derivedSize();
    }

private:
    StaticOrDynamicSizeContainer<Size>         m_originalSize;
    StaticOrDynamicSizeContainer<2 * Size + 1> m_derivedSize;
    Eigen::Matrix<double, Size, 1>             m_values;
};

TEST_CASE("static or dynamic size") {
    SECTION("simple example") {
        SECTION("static") {
            constexpr SimpleExample<5> example;
            constexpr auto             size = example.size();
            REQUIRE(size == 5);
        }

        SECTION("dynamic") {
            const SimpleExample<Eigen::Dynamic> example(5);
            const auto                          size = example.size();
            REQUIRE(size == 5);
        }
    }

    SECTION("operator overloading") {
        SECTION("static") {
            constexpr StaticOrDynamicSizeContainer<1> one;
            constexpr StaticOrDynamicSizeContainer<2> two;
            constexpr auto                            three = one + two;
            static_assert(three() == 3);
        }
//        SECTION("both dynamic") {
//            constexpr StaticOrDynamicSizeContainer<Eigen::Dynamic> one(1);
//            constexpr StaticOrDynamicSizeContainer<Eigen::Dynamic> two(2);
//            auto                                                   three = one + two;
//            REQUIRE(three() == 3);
//        }
//        SECTION("first dynamic") {
//            constexpr StaticOrDynamicSizeContainer<Eigen::Dynamic> one(1);
//            constexpr StaticOrDynamicSizeContainer<2> two;
//            auto                                                   three = one + two;
//            REQUIRE(three() == 3);
//        }
        SECTION("second dynamic") {
            constexpr StaticOrDynamicSizeContainer<1> one;
            constexpr StaticOrDynamicSizeContainer<Eigen::Dynamic> two(2);
            auto                                                   three = one + two;
            REQUIRE(three() == 3);
        }
    }

    SECTION("derived example") {
        SECTION("static") {
            const DerivedExample<5> example;
            constexpr auto          original_size = example.originalSize();
            constexpr auto          derived_size  = example.derivedSize();
            REQUIRE(original_size == 5);
            REQUIRE(derived_size == 11);
        }

        SECTION("dynamic") {
            const DerivedExample<Eigen::Dynamic> example(5);
            const auto                           original_size = example.originalSize();
            const auto                           derived_size  = example.derivedSize();
            REQUIRE(original_size == 5);
            REQUIRE(derived_size == 11);
        }
    }
}
