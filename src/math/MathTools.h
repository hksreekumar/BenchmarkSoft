#pragma once

#include <fstream>
#include <map>
#include <vector>

#ifdef USE_INTEL_MKL
#define MKL_DIRECT_CALL 1
#include <mkl.h>
#endif


namespace MathTools {
    void testfunction();

}
