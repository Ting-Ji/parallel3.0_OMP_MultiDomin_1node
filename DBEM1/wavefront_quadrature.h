#ifndef DBEM_WAVEFRONT_QUADRATURE_H
#define DBEM_WAVEFRONT_QUADRATURE_H
#include <array>
#include <cstddef>
#include <memory>
#include <string>
class DSquareElement;
class DynaMat;
namespace Wavefront {
using Vec = std::array<double,3>;
enum class Kernel { U, T1, T2 };
enum class Status { Outside, Converged, Failed };
struct Options {
 double relativeTolerance, scaledAbsoluteTolerance;
 unsigned maxDepth; std::size_t maxEvaluations;
 Options();
};
struct Statistics {
 std::size_t evaluations=0, leaves=0, skipped=0;
 unsigned maxDepth=0; bool planar=false; double missedSupportError=0;
};
struct Result {
 std::array<double,9> block{};
 double error=0, scale=0;
 // Sampled norm, full-panel characteristic scale, and actual stopping target.
 double characteristicScale=0, absoluteScale=0, tolerance=0;
 Status status=Status::Failed;
 Statistics statistics;
 std::string message;
};
struct Geometry {
 double coefficient[3][8]{};
 std::array<Vec,4> vertices{};
 Vec origin{}, normal{}, axis1{}, axis2{};
 double diameter=0, boundingRadius=0;
 bool planar=false, valid=false;
 std::string reason;
};
std::shared_ptr<const Geometry> BuildGeometry(DSquareElement& element);
Result Integrate(DSquareElement& element, Kernel kernel, const Vec& source,
 long lag, const DynaMat& material, const Options& options=Options());
int Classify(DSquareElement& element, const Vec& source, long lag, bool secondTimeWindow, const DynaMat& material);
const char* Name(Kernel kernel);
void RecordWorkerFailure(const char* message) noexcept;
void ThrowWorkerFailure();
}
#endif
