#include <cstdio>
#include <exception>
void RhsCpuRunTests();
int main()
{
    try { RhsCpuRunTests(); printf("PASS: RHS exact operators, maps, convolution, G/T order, repeated tasks, diagnostics, nonfinite, thread failure, shutdown; threads=40\n"); return 0; }
    catch (const std::exception& e) { fprintf(stderr, "FAIL: %s\n", e.what()); return 2; }
}
