#ifndef BOOST_HELPERS_HPP
#define BOOST_HELPERS_HPP

/** Compare two iterable containers with tolerance.
 *
 * @param aa First container, e.g. of type `std::vector<double>`.
 * @param bb Second container.
 * @param tolerance Desired maximum deviation between values in first and second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_SMALL for comparing the container contents with given
 * tolerance.
 */
#define REQUIRE_CLOSE_COLLECTIONS(aa, bb, tolerance, containerName)                          \
    {                                                                                        \
        using std::distance;                                                                 \
        using std::begin;                                                                    \
        using std::end;                                                                      \
        auto a = begin(aa), a0 = begin(aa), ae = end(aa);                                    \
        auto b = begin(bb);                                                                  \
        BOOST_TEST_INFO("Vector \"" << containerName << "\" length comparison");             \
        BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb)));                          \
        for (; a != ae; ++a, ++b)                                                            \
        {                                                                                    \
            BOOST_TEST_INFO("Vector \"" << containerName << "\" element " << distance(a0, a) \
                                        << " comparison: " << *a << " vs. " << *b);          \
            BOOST_REQUIRE_SMALL(*a - *b, tolerance);                                         \
        }                                                                                    \
    }

/** Compare two two-dimensional iterable containers with tolerance.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<double>>`.
 * @param bb Second container.
 * @param tolerance Desired maximum deviation between values in first and second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_SMALL for comparing the container contents with given
 * tolerance.
 */
#define REQUIRE_CLOSE_COLLECTIONS_2D(aaa, bbb, tolerance, containerName)                          \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aa = begin(aaa), aa0 = begin(aaa), aae = end(aaa);                                   \
        auto bb = begin(bbb);                                                                     \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(bbb)));                           \
        for (; aa != aae; ++aa, ++bb)                                                             \
        {                                                                                         \
            auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                                  \
            auto b = begin(*bb);                                                                  \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa)     \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                          \
            for (; a != ae; ++a, ++b)                                                             \
            {                                                                                     \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa) \
                                           << "," << distance(a0, a) << " comparison: " << *a     \
                                           << " vs. " << *b);                                     \
                BOOST_REQUIRE_SMALL(*a - *b, tolerance);                                          \
            }                                                                                     \
        }                                                                                         \
    }

/** Compare two three-dimensional iterable containers with tolerance.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<std::vector<double>>>`.
 * @param bb Second container.
 * @param tolerance Desired maximum deviation between values in first and second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_SMALL for comparing the container contents with given
 * tolerance.
 */
#define REQUIRE_CLOSE_COLLECTIONS_3D(aaaa, bbbb, tolerance, containerName)                        \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aaa = begin(aaaa), aaa0 = begin(aaaa), aaae = end(aaaa);                             \
        auto bbb = begin(bbbb);                                                                   \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aaa, aaae), distance(bbb, end(bbbb)));                       \
        for (; aaa != aaae; ++aaa, ++bbb)                                                         \
        {                                                                                         \
            auto aa = begin(*aaa), aa0 = begin(*aaa), aae = end(*aaa);                            \
            auto bb = begin(*bbb);                                                                \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aaa0, aaa)   \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(*bbb)));                      \
            for (; aa != aae; ++aa, ++bb)                                                         \
            {                                                                                     \
                auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                              \
                auto b = begin(*bb);                                                              \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element "                      \
                                           << distance(aaa0, aaa) << "," << distance(aa0, aa)     \
                                           << " 3rd dimension comparison");                       \
                BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                      \
                for (; a != ae; ++a, ++b)                                                         \
                {                                                                                 \
                    BOOST_TEST_INFO("Array \"" << containerName << "\" element "                  \
                                               << distance(aaa0, aaa) << "," << distance(aa0, aa) \
                                               << "," << distance(a0, a) << " comparison: " << *a \
                                               << " vs. " << *b);                                 \
                    BOOST_REQUIRE_SMALL(*a - *b, tolerance);                                      \
                }                                                                                 \
            }                                                                                     \
        }                                                                                         \
    }

/** Compare equality of two iterable containers.
 *
 * @param aa First container, e.g. of type `std::vector<int>`.
 * @param bb Second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_EQUAL for comparing the container contents.
 */
#define REQUIRE_EQUAL_COLLECTIONS(aa, bb, containerName)                                     \
    {                                                                                        \
        using std::distance;                                                                 \
        using std::begin;                                                                    \
        using std::end;                                                                      \
        auto a = begin(aa), a0 = begin(aa), ae = end(aa);                                    \
        auto b = begin(bb);                                                                  \
        BOOST_TEST_INFO("Vector \"" << containerName << "\" length comparison");             \
        BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(bb)));                          \
        for (; a != ae; ++a, ++b)                                                            \
        {                                                                                    \
            BOOST_TEST_INFO("Vector \"" << containerName << "\" element " << distance(a0, a) \
                                        << " comparison: " << *a << " vs. " << *b);          \
            BOOST_REQUIRE_EQUAL(*a, *b);                                                     \
        }                                                                                    \
    }

/** Compare equality of two two-dimensional iterable containers.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<int>>`.
 * @param bb Second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_EQUAL for comparing the container contents.
 */
#define REQUIRE_EQUAL_COLLECTIONS_2D(aaa, bbb, containerName)                                     \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aa = begin(aaa), aa0 = begin(aaa), aae = end(aaa);                                   \
        auto bb = begin(bbb);                                                                     \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(bbb)));                           \
        for (; aa != aae; ++aa, ++bb)                                                             \
        {                                                                                         \
            auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                                  \
            auto b = begin(*bb);                                                                  \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa)     \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                          \
            for (; a != ae; ++a, ++b)                                                             \
            {                                                                                     \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aa0, aa) \
                                           << "," << distance(a0, a) << " comparison: " << *a     \
                                           << " vs. " << *b);                                     \
                BOOST_REQUIRE_EQUAL(*a, *b);                                                      \
            }                                                                                     \
        }                                                                                         \
    }

/** Compare equality of two three-dimensional iterable containers.
 *
 * @param aa First container, e.g. of type `std::vector<std::vector<std::vector<int>>>`.
 * @param bb Second container.
 * @param containerName Descriptive name of the containers' contents.
 *
 * Internally calls BOOST_REQUIRE_EQUAL for size comparison and
 * BOOST_REQUIRE_EQUAL for comparing the container contents.
 */
#define REQUIRE_EQUAL_COLLECTIONS_3D(aaaa, bbbb, containerName)                                   \
    {                                                                                             \
        using std::distance;                                                                      \
        using std::begin;                                                                         \
        using std::end;                                                                           \
        auto aaa = begin(aaaa), aaa0 = begin(aaaa), aaae = end(aaaa);                             \
        auto bbb = begin(bbbb);                                                                   \
        BOOST_TEST_INFO("Array \"" << containerName << "\" 1st dimension comparison");            \
        BOOST_REQUIRE_EQUAL(distance(aaa, aaae), distance(bbb, end(bbbb)));                       \
        for (; aaa != aaae; ++aaa, ++bbb)                                                         \
        {                                                                                         \
            auto aa = begin(*aaa), aa0 = begin(*aaa), aae = end(*aaa);                            \
            auto bb = begin(*bbb);                                                                \
            BOOST_TEST_INFO("Array \"" << containerName << "\" element " << distance(aaa0, aaa)   \
                                       << " 2nd dimension comparison");                           \
            BOOST_REQUIRE_EQUAL(distance(aa, aae), distance(bb, end(*bbb)));                      \
            for (; aa != aae; ++aa, ++bb)                                                         \
            {                                                                                     \
                auto a = begin(*aa), a0 = begin(*aa), ae = end(*aa);                              \
                auto b = begin(*bb);                                                              \
                BOOST_TEST_INFO("Array \"" << containerName << "\" element "                      \
                                           << distance(aaa0, aaa) << "," << distance(aa0, aa)     \
                                           << " 3rd dimension comparison");                       \
                BOOST_REQUIRE_EQUAL(distance(a, ae), distance(b, end(*bb)));                      \
                for (; a != ae; ++a, ++b)                                                         \
                {                                                                                 \
                    BOOST_TEST_INFO("Array \"" << containerName << "\" element "                  \
                                               << distance(aaa0, aaa) << "," << distance(aa0, aa) \
                                               << "," << distance(a0, a) << " comparison: " << *a \
                                               << " vs. " << *b);                                 \
                    BOOST_REQUIRE_EQUAL(*a, *b);                                                  \
                }                                                                                 \
            }                                                                                     \
        }                                                                                         \
    }

/// Check equality with context message.
#define CHECK_EQUAL_MESSAGE(L, R, M) \
    {                                \
        BOOST_TEST_INFO(M);          \
        BOOST_CHECK_EQUAL(L, R);     \
    }
/// Check equality with context message (warning).
#define WARN_EQUAL_MESSAGE(L, R, M) \
    {                               \
        BOOST_TEST_INFO(M);         \
        BOOST_WARN_EQUAL(L, R);     \
    }
/// Check equality with context message (required).
#define REQUIRE_EQUAL_MESSAGE(L, R, M) \
    {                                  \
        BOOST_TEST_INFO(M);            \
        BOOST_REQUIRE_EQUAL(L, R);     \
    }
#endif
