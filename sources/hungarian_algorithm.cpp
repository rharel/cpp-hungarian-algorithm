#include <limits>

#include "../include/hungarian_algorithm.h"


using namespace rharel::hungarian_algorithm;


void rharel::hungarian_algorithm::solve_for_minimum_cost_assignment(
    const unsigned int problem_size,
    unsigned int**     cost_matrix,
    unsigned int*      assignment)
{
    if (problem_size == 0) {                    return; }
    if (problem_size == 1) { assignment[0] = 0; return; }

    Problem(problem_size, cost_matrix).solve(assignment);
}

Problem::Problem(const unsigned int size, 
                 unsigned int *const *const cost_matrix)
    : n(size), C(cost_matrix), S(n, n), P(n, n),
      is_covered_row(n, false), is_covered_column(n, false)
{}
bool Problem::step()
{
    switch (current_step)
    {
        case Step::One:   { current_step = step_1(); break; }
        case Step::Two:   { current_step = step_2(); break; }
        case Step::Three: { current_step = step_3(); break; }
        case Step::Four:  { current_step = step_4(); break; }
        case Step::Five:  { current_step = step_5(); break; }
        case Step::Six:   { current_step = step_6(); break; }
    }
    return current_step == Step::Done;
}
void Problem::solve(unsigned int *const assignment)
{
    while (current_step != Step::Done) { step(); }
    output_solution(assignment);
}
Problem::Step Problem::step_1()
{
    S.reserve(Eigen::VectorXi::Constant(n, 4));
    P.reserve(Eigen::VectorXi::Constant(n, 4));

    for (unsigned int i = 0; i < n; ++i)
    {
        add_to_row(i, /* value: */ -1 * minimum_in_row(i));
    }
    return Step::Two;
}
Problem::Step Problem::step_2()
{
    std::vector<bool> is_starred_column(n);
    for (unsigned int i = 0; i < n; ++i)
    {
        for (unsigned int j = 0; j < n; ++j)
        {
            if  (!is_starred_column[j] && C[i][j] == 0)
            {
                S.insert(i, j)       = true;
                is_starred_column[j] = true;
                break;  // Since there is already a starred zero in this row.
            }
        }
    }
    return Step::Three;
}
Problem::Step Problem::step_3()
{
    unsigned int covered_column_count = 0;
    for (unsigned int j = 0; j < n; ++j)
    {
        unsigned int i;  // Dummy.
        if (find_starred_zero_in_column(j, i))
        {
            is_covered_column[j] = true; 
            ++ covered_column_count;
        }
    }
    if (covered_column_count == n) { return Step::Done; }
    else                           { return Step::Four; }
}
Problem::Step Problem::step_4()
{
    unsigned int i, j;
    while (find_uncovered_zero(i, j))
    {
        P.coeffRef(i, j) = true;

        if (find_starred_zero_in_row(i, j))
        {
            is_covered_row[i]    = true;
            is_covered_column[j] = false;
        }
        else
        {
            uncovered_prime_zero[0] = i;
            uncovered_prime_zero[1] = j;
            return Step::Five;
        }
    }
    return Step::Six;
}
Problem::Step Problem::step_5()
{
    typedef std::pair<unsigned int, unsigned int> Location;
    std::vector<Location> primed_zeroes;
    std::vector<Location> starred_zeroes;

    unsigned int i = uncovered_prime_zero[0],
                 j = uncovered_prime_zero[1];

    primed_zeroes.push_back(Location(i, j));
    while (find_starred_zero_in_column(j, i))
    {
        starred_zeroes.push_back(Location(i, j));
        find_primed_zero_in_row(i, j);
        primed_zeroes.push_back(Location(i, j));
    }
    
    for (const auto& location : primed_zeroes) 
    { 
        S.coeffRef(location.first, location.second) = true; 
    }
    for (const auto& location : starred_zeroes)
    {
        S.coeffRef(location.first, location.second) = false;
    }

    P.setZero();

    std::fill(is_covered_row.begin(), 
              is_covered_row.end(), false);
    std::fill(is_covered_column.begin(), 
              is_covered_column.end(), false);

    return Step::Three;
}
Problem::Step Problem::step_6() 
{
    const unsigned int m = minimum_uncovered();
    for (unsigned int i = 0; i < n; ++i)
    {
        if (is_covered_row[i]) { add_to_row(i, m); }
    }
    for (unsigned int j = 0; j < n; ++j)
    {
        if (!is_covered_column[j]) { add_to_column(j, -1 * m); }
    }
    return Step::Four;
}
void Problem::output_solution(unsigned int* assignment) const
{
    for (unsigned int j = 0; j < n; ++j)
    {
        for (ColumnMemberIterator cmi(S, j); cmi; ++cmi)
        {
            if (cmi.value()) { assignment[cmi.row()] = j; break; }
        }
    }
}

unsigned int Problem::minimum_in_row(const unsigned int i) const
{
    unsigned int minimum = C[i][0];
    for (unsigned int j = 1; j < n; ++j)
    {
        if (C[i][j] < minimum) { minimum = C[i][j]; }
    }
    return minimum;
}
unsigned int Problem::minimum_uncovered() const
{
    unsigned int minimum = std::numeric_limits<unsigned int>::max();
    for (unsigned int i = 0; i < n; ++i)
    {
        if (is_covered_row[i]) { continue; }
        for (unsigned int j = 0; j < n; ++j)
        {
            if (is_covered_column[j]) { continue; }

            if (C[i][j] < minimum) { minimum = C[i][j]; }
        }
    }
    return minimum;
}

void Problem::add_to_row(const unsigned int i, const int value)
{
    for (unsigned int j = 0; j < n; ++j) { C[i][j] += value; }
}
void Problem::add_to_column(const unsigned int j, const int value)
{
    for (unsigned int i = 0; i < n; ++i) { C[i][j] += value; }
}

bool Problem::find_uncovered_zero(unsigned int& row_index,
                                  unsigned int& column_index) const
{
    for (unsigned int i = 0; i < n; ++i)
    {
        if (is_covered_row[i]) { continue; }
        for (unsigned int j = 0; j < n; ++j)
        {
            if (is_covered_column[j]) { continue; }

            if (C[i][j] == 0)
            {
                row_index    = i;
                column_index = j;
                return true;
            }
        }
    }
    return false;
}

bool Problem::find_starred_zero_in_row(
    const unsigned int row_index,
    unsigned int&      column_index) const
{
    for (unsigned int j = 0; j < n; ++j)
    {
        for (ColumnMemberIterator cmi(S, j); cmi; ++cmi)
        {
            const unsigned int i = cmi.row();
            
            if (i > row_index) { break; }
            if (i == row_index && cmi.value()) 
            { 
                column_index = j; 
                return true; 
            }
        }
    }
    return false;
}
bool Problem::find_starred_zero_in_column(
    const unsigned int column_index,
    unsigned int&      row_index) const
{
    for (ColumnMemberIterator cmi(S, column_index); cmi; ++cmi)
    {
        if (cmi.value()) 
        { 
            row_index = cmi.row(); 
            return true; 
        }
    }
    return false;
}
bool Problem::find_primed_zero_in_row(
    const unsigned int row_index,
    unsigned int&      column_index) const
{
    for (unsigned int j = 0; j < n; ++j)
    {
        for (ColumnMemberIterator cmi(P, j); cmi; ++cmi)
        {
            const unsigned int i = cmi.row();

            if (i > row_index) { break; }
            if (i == row_index && cmi.value()) 
            { 
                column_index = j; 
                return true; 
            }
        }
    }
    return false;
}
