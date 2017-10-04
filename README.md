# The Hungarian Algorithm
The [Hungarian algorithm](https://en.wikipedia.org/wiki/Hungarian_algorithm) is used to solve the following [problem](https://en.wikipedia.org/wiki/Assignment_problem):

> Given `n > 0` workers, `n` tasks, and the cost matrix `C` whose member `C(i, j)` represents the cost of assigning the `i`-th worker to the `j`-th task, find an assignment of workers to tasks with minimal total cost.

The Hungarian algorithm solves this problem in O(n<sup>3</sup>) time.

# Directory structure

| Directory        | Description               |
| ---------        | -----------               |
|`visual_studio/`  | Solution files (VS 2017). |
| `include/`       | Headers.                  |
| `sources/`       | *.cpp sources.            |
| `libraries/`     | 3rd party dependencies.   |
| `documentation/` | [Doxygen](http://www.stack.nl/~dimitri/doxygen/) configuration. |
| `tests/`         | Test project files.       |
| `benchmarks/`    | Benchmark project files.  |
