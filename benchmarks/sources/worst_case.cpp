#include <hayai/hayai.hpp>

#include <hungarian_algorithm.h>


using rharel::hungarian_algorithm::solve_for_minimum_cost_assignment;


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
}
BENCHMARK(Hungarian_Algorithm_Worst_Case,   n_5, 1000, 1000) { test(5);   }
BENCHMARK(Hungarian_Algorithm_Worst_Case,  n_10,  500, 1000) { test(10);  }
BENCHMARK(Hungarian_Algorithm_Worst_Case,  n_50,   10,  100) { test(50);  }
BENCHMARK(Hungarian_Algorithm_Worst_Case, n_100,   10,   10) { test(100); }
