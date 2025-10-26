#pragma once

#include "Point.h"

#include <optional>

namespace ZXing {

class BitMatrix;

namespace OneD {

struct OrientationEstimate
{
        PointF scanDir;   // direction perpendicular to bars, used for sampling
        PointF offsetDir; // direction along the bars
        double eccentricity = 0; // relative difference between eigenvalues
};

std::optional<OrientationEstimate> EstimateOrientation(const BitMatrix& matrix);

} // namespace OneD
} // namespace ZXing
