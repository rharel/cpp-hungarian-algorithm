#include <hayai/hayai.hpp>


int main()
{ 
    hayai::ConsoleOutputter console_outputter;

    hayai::Benchmarker::AddOutputter(console_outputter);
    hayai::Benchmarker::RunAllTests();

    return 0;
}
