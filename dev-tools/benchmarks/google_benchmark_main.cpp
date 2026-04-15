
#include <iostream>

#include <benchmark/benchmark.h>

namespace {

void print_stats() {
    std::cout << "MorphoTreeAdjust Google Benchmark\n\n";
}

} // namespace

int main(int argc, char **argv) {
    print_stats();
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv)) {
        return 1;
    }
    benchmark::RunSpecifiedBenchmarks();
    return 0;
}
