#include <vector>

#include <catch.hpp>

#include <hungarian_algorithm.h>


using rharel::hungarian_algorithm::solve_for_minimum_cost_assignment;


/// Requires that the specified solution is the anti diagonal of an nxn matrix.
void require_solution_is_anti_diagonal(const unsigned int* const solution, 
                                       const unsigned int n)
{
    for (unsigned int i = 0; i < n; ++i) 
    { 
        REQUIRE(solution[i] == n - 1 - i); 
    }
}
/// Solves for an nxn cost matrix C where C(i, j) = (i + 1) * (j + 1).
/// Note: i, j are zero-based.
void test(const unsigned int n)
{
    std::vector<unsigned int> solution(n);
    solve_for_minimum_cost_assignment(
        /* problem_size:  */ n,
        /* cost function: */ [](const unsigned int i, 
                                const unsigned int j) 
                                { 
                                    return (i + 1) * (j + 1); 
                                }, 
        solution.data()
    );
    require_solution_is_anti_diagonal(solution.data(), n);
}
TEST_CASE("Solves worst case [n =  0]") { test(0);  }
TEST_CASE("Solves worst case [n =  1]") { test(1);  }
TEST_CASE("Solves worst case [n =  2]") { test(2);  }
TEST_CASE("Solves worst case [n =  3]") { test(3);  }
TEST_CASE("Solves worst case [n =  4]") { test(4);  }
TEST_CASE("Solves worst case [n =  5]") { test(5);  }
TEST_CASE("Solves worst case [n = 10]") { test(10); }
TEST_CASE("Solves worst case [n = 50]") { test(50); }
