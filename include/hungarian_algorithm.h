#pragma once

#include <vector>

#include <eigen/Sparse>


/// Contains methods for solving the assignment problem [1] using the Hungarian
/// algorithm [2].
/// 
/// The problem statement: Given n > 0 workers, n tasks, and the cost matrix C 
/// whose member C(i, j) represents the cost of assigning the ith worker to the 
/// jth task, find an assignment of workers to tasks with minimal total cost.
///
/// # References
/// 1. https://en.wikipedia.org/wiki/Assignment_problem
/// 2. https://en.wikipedia.org/wiki/Hungarian_algorithm
namespace rharel::hungarian_algorithm
{
    /// Solves the assignment problem from a given cost matrix.
    ///
    /// @param problem_size
    ///     The number of workers/tasks.
    /// @param cost_matrix
    ///     A square matrix with problem_size rows and columns. Member (i, j)
    ///     represents the cost of assigning worker i to task j.
    /// @param[out] assignment
    ///     An output buffer for the minimum cost assignment.
    void solve_for_minimum_cost_assignment(unsigned int   problem_size,
                                           unsigned int** cost_matrix,
                                           unsigned int*  assignment);

    /// Builds the cost matrix and solves the assignment problem.
    ///
    /// @tparam CostComputer
    ///     The type of a function-like object:
    ///     unsigned int (*)(unsigned int i, unsigned int j);
    ///
    /// @param problem_size
    ///     The number of workers/tasks.
    /// @param compute_cost
    ///     Computes the cost of assigning worker i to task j.
    /// @param[out] assignment
    ///     An output buffer for the minimum cost assignment.
    template <typename CostComputer>
    void solve_for_minimum_cost_assignment(unsigned int        problem_size, 
                                           const CostComputer& compute_cost,
                                           unsigned int*       assignment);

    /// Builds the cost matrix and solves the assignment problem.
    ///
    /// @tparam Worker
    ///     The worker type.
    /// @tparam Task
    ///     The task type.
    /// @tparam CostComputer
    ///     The type of a function-like object:
    ///     unsigned int (*)(unsigned int i, unsigned int j, 
    ///                      const Worker& worker, const Task& task);
    ///
    /// @param problem_size
    ///     The number of workers/tasks.
    /// @param workers
    ///     A list of workers.
    /// @param tasks
    ///     A list of tasks.
    /// @param compute_cost
    ///     Computes the cost of assigning worker i to task j.
    /// @param[out] assignment
    ///     An output buffer for the minimum cost assignment.
    template <class Worker, class Task, typename CostComputer>
    void solve_for_minimum_cost_assignment(unsigned int        problem_size, 
                                           const Worker*       workers, 
                                           const Task*         tasks,
                                           const CostComputer& compute_cost,
                                           unsigned int*       assignment);

    /// Represents an assignment problem instance.
    class Problem
    {
        public:
        /// Creates a new problem with the specified number of workers/tasks 
        /// and assignment cost matrix.
        explicit Problem(unsigned int size, 
                         unsigned int *const *const cost_matrix);

        /// Performs one step towards a solution.
        ///
        /// That is, performs one step out of the 6-step procedure from [1].
        ///
        /// Returns true iff the solution is ready. It can be retrieved by 
        /// invoking solve().
        ///
        /// # References
        /// 1. http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
        bool step();

        /// Solves the problem and writes the solution onto the specified 
        /// output buffer.
        void solve(unsigned int* assignment);

        private:
        /// Enumerates steps of the algorithm.
        ///
        /// We use the 6-step procedure from [1] to arrive at a solution.
        ///
        /// # References
        /// 1. http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
        enum class Step : int
        {
            Done = 0,
            One, Two, Three, Four, Five, Six
        };

        /// For each row in C, subtracts its members by the minimum amongst 
        /// them. Proceeds to step 2.
        Step step_1();
        /// Finds an unstarred zero z = C(i, j), if there is no starred zero 
        /// in either C(i, *) or C(*, j), star z. Repeats for each member of C.
        /// Proceeds to step 3.
        Step step_2();
        /// Covers all columns with a starred zero. If n columns were covered,
        /// we are done. Otherwise, proceeds to step 4.
        Step step_3();
        /// Finds a non covered zero z_p = C(i, j) and primes it. If there is 
        /// no starred zero z_s = C(i, k) in the same row, proceeds to step 5. 
        /// Otherwise, covers row i and uncovers column k. Repeats until there 
        /// are no uncovered zeros left, in which case proceeds to step 6.
        Step step_4();
        /// Traverses a sequence of alternating primed and starred zeros: 
        /// Let z_0 represent the uncovered primed zero found in step 4. 
        /// Let z_1 denote the starred zero in the column of z_0 (if any). 
        /// Let z_2 denote the primed zero in the row of z_1 (there will always
        /// be one). The sequence continues until a primed zero that has no 
        /// starred zero in its column is reached. Unstars each starred zero 
        /// and stars each primed zero of the sequence. Erases all primes and 
        /// uncovers all rows and columns in C. Proceeds to step 3.
        Step step_5();
        /// Finds the smallest uncovered member m of C. Adds m to covered rows
        /// and subtracts it from uncovered columns. Proceeds to step 4.
        Step step_6();
        /// Writes the current assignment to the specified buffer.
        void output_solution(unsigned int* assignment) const;

        /// Finds the minimum member of C(i, *).
        unsigned int minimum_in_row(unsigned int i) const;
        /// Finds the minimum uncovered member of C.
        /// If C does not contain uncovered members, returns the maximum 
        /// unsigned value instead.
        unsigned int minimum_uncovered() const;

        /// Adds the specified value to C(i, *).
        void add_to_row(unsigned int i, int value);
        /// Adds the specified value to C(*, j).
        void add_to_column(unsigned int j, int value);

        /// Finds a non-covered zero and reports its location.
        /// Returns true iff one was found.
        bool find_uncovered_zero(unsigned int& i, 
                                 unsigned int& j) const;

        /// Finds a starred zero in the specified row and reports its column.
        /// Returns true iff one was found.
        bool find_starred_zero_in_row(unsigned int  i, 
                                      unsigned int& j) const;
        /// Finds a starred zero in the specified column and reports its row.
        /// Returns true iff one was found.
        bool find_starred_zero_in_column(unsigned int  j, 
                                         unsigned int& i) const;
        /// Finds a primed zero in the specified row and reports its column.
        /// Returns true iff one was found.
        bool find_primed_zero_in_row(unsigned int  i, 
                                     unsigned int& j) const;

        /// An iterator over members of a column in a sparse mask matrix.
        typedef Eigen::SparseMatrix<bool>::InnerIterator 
                ColumnMemberIterator;

        const unsigned int n;          // Problem size.
        unsigned int *const *const C;  // Cost matrix.
        Eigen::SparseMatrix<bool> S;   // Stars matrix.
        Eigen::SparseMatrix<bool> P;   // Primes matrix.

        std::vector<bool> is_covered_row, 
                          is_covered_column;

        unsigned int uncovered_prime_zero[2];  // Output of step 4.

        Step current_step = Step::One;
    };
}

#include "hungarian_algorithm.hpp"
