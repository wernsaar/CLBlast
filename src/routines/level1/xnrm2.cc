
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This file implements the Xnrm2 class (see the header for information about the class).
//
// =================================================================================================

#include "internal/routines/level1/xnrm2.h"

#include <string>
#include <vector>

namespace clblast {
// =================================================================================================

// Specific implementations to get the memory-type based on a template argument
template <> const Precision Xnrm2<float>::precision_ = Precision::kSingle;
template <> const Precision Xnrm2<double>::precision_ = Precision::kDouble;
template <> const Precision Xnrm2<float2>::precision_ = Precision::kComplexSingle;
template <> const Precision Xnrm2<double2>::precision_ = Precision::kComplexDouble;

// =================================================================================================

// Constructor: forwards to base class constructor
template <typename T>
Xnrm2<T>::Xnrm2(Queue &queue, EventPointer event, const std::string &name):
    Routine<T>(queue, event, name, {"Xdot"}, precision_) {
  source_string_ =
    #include "../../kernels/level1/xnrm2.opencl"
  ;
}

// =================================================================================================

// The main routine
template <typename T>
StatusCode Xnrm2<T>::DoNrm2(const size_t n,
                            const Buffer<T> &nrm2_buffer, const size_t nrm2_offset,
                            const Buffer<T> &x_buffer, const size_t x_offset, const size_t x_inc) {

  // Makes sure all dimensions are larger than zero
  if (n == 0) { return StatusCode::kInvalidDimension; }

  // Tests the vectors for validity
  auto status = TestVectorX(n, x_buffer, x_offset, x_inc, sizeof(T));
  if (ErrorIn(status)) { return status; }
  status = TestVectorDot(1, nrm2_buffer, nrm2_offset, sizeof(T));
  if (ErrorIn(status)) { return status; }

  // Retrieves the Xnrm2 kernels from the compiled binary
  try {
    const auto program = GetProgramFromCache();
    auto kernel1 = Kernel(program, "Xnrm2");
    auto kernel2 = Kernel(program, "Xnrm2Epilogue");

    // Creates the buffer for intermediate values
    auto temp_size = 2*db_["WGS2"];
    auto temp_buffer = Buffer<T>(context_, temp_size);

    // Sets the kernel arguments
    kernel1.SetArgument(0, static_cast<int>(n));
    kernel1.SetArgument(1, x_buffer());
    kernel1.SetArgument(2, static_cast<int>(x_offset));
    kernel1.SetArgument(3, static_cast<int>(x_inc));
    kernel1.SetArgument(4, temp_buffer());

    // Event waiting list
    auto eventWaitList = std::vector<Event>();

    // Launches the main kernel
    auto global1 = std::vector<size_t>{db_["WGS1"]*temp_size};
    auto local1 = std::vector<size_t>{db_["WGS1"]};
    auto kernelEvent = Event();
    status = RunKernel(kernel1, global1, local1, kernelEvent.pointer());
    if (ErrorIn(status)) { return status; }
    eventWaitList.push_back(kernelEvent);

    // Sets the arguments for the epilogue kernel
    kernel2.SetArgument(0, temp_buffer());
    kernel2.SetArgument(1, nrm2_buffer());
    kernel2.SetArgument(2, static_cast<int>(nrm2_offset));

    // Launches the epilogue kernel
    auto global2 = std::vector<size_t>{db_["WGS2"]};
    auto local2 = std::vector<size_t>{db_["WGS2"]};
    status = RunKernel(kernel2, global2, local2, event_, eventWaitList);
    if (ErrorIn(status)) { return status; }

    // Succesfully finished the computation
    return StatusCode::kSuccess;
  } catch (...) { return StatusCode::kInvalidKernel; }
}

// =================================================================================================

// Compiles the templated class
template class Xnrm2<float>;
template class Xnrm2<double>;
template class Xnrm2<float2>;
template class Xnrm2<double2>;

// =================================================================================================
} // namespace clblast
