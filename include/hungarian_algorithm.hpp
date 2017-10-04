#include "hungarian_algorithm.h"


/// Contains private implementation details.
namespace rharel::hungarian_algorithm::detail
{
    /// Computes the cost matrix.
    template <typename CostComputer>
    inline void compute_cost_matrix(/* problem size:  */ const unsigned int  n,
                                    /* cost function: */ const CostComputer& c,
                                    /* cost matrix:   */ unsigned int**      C)
    {
        for (unsigned int i = 0; i < n; ++i)
        {
            for (unsigned int j = 0; j < n; ++j)
            {
                C[i][j] = static_cast<unsigned int>(c(i, j));
            }
        }
    }
}
namespace rharel::hungarian_algorithm
{
    template <typename CostComputer>
    void solve_for_minimum_cost_assignment(
        /* problem size:  */ const unsigned int  n,
        /* cost function: */ const CostComputer& c,
        /* assignment:    */ unsigned int*       A)

    {
        using std::vector;
        using detail::compute_cost_matrix;

        vector<vector<unsigned int>> C(n);
        vector<unsigned int*>        C_row_pointers(n);
        for (unsigned int i = 0; i < n; ++i)
        {
            C[i].resize(n);
            C_row_pointers[i] = C[i].data();
        }
        unsigned int** C_raw = C_row_pointers.data();
        compute_cost_matrix(n, c, C_raw); 
        solve_for_minimum_cost_assignment(n, C_raw, A); 
    }
    template <class Worker, class Task, typename CostComputer>
    void solve_for_minimum_cost_assignment(
        /* problem size:  */ const unsigned int  n,
        /* workers:       */ const Worker*       W,
        /* tasks:         */ const Task*         T,
        /* cost function: */ const CostComputer& c,
        /* assignment:    */ unsigned int*       A)

    {
        solve_for_minimum_cost_assignment(
            /* problem_size:  */ n, 
            /* cost function: */ [&c, &W, &T] (const unsigned int i, 
                                               const unsigned int j) 
                                 { 
                                    return c(i, j, W[i], T[j]); 
                                 }, 
            /* assignment:    */ A
        );
    }
}
