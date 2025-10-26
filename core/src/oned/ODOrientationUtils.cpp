#include "ODOrientationUtils.h"

#include "BitMatrix.h"

#include <algorithm>
#include <cmath>

namespace ZXing::OneD {

namespace {
constexpr double kEpsilon = 1e-6;
}

std::optional<OrientationEstimate> EstimateOrientation(const BitMatrix& matrix)
{
        const int width = matrix.width();
        const int height = matrix.height();

        long long count = 0;
        double meanX = 0;
        double meanY = 0;

        for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                        if (matrix.get(x, y)) {
                                meanX += x;
                                meanY += y;
                                ++count;
                        }
                }
        }

        if (count == 0)
                return std::nullopt;

        meanX /= count;
        meanY /= count;

        double sxx = 0;
        double syy = 0;
        double sxy = 0;
        for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                        if (!matrix.get(x, y))
                                continue;
                        double dx = x - meanX;
                        double dy = y - meanY;
                        sxx += dx * dx;
                        syy += dy * dy;
                        sxy += dx * dy;
                }
        }

        if (sxx + syy < kEpsilon)
                return std::nullopt;

        double angle = 0.5 * std::atan2(2 * sxy, sxx - syy);
        double trace = sxx + syy;
        double diff = sxx - syy;
        double root = std::sqrt(std::max(0.0, diff * diff + 4 * sxy * sxy));
        double lambda1 = (trace + root) / 2;
        double lambda2 = (trace - root) / 2;

        if (lambda1 <= kEpsilon)
                return std::nullopt;

        double eccentricity = lambda1 > 0 ? (lambda1 - lambda2) / lambda1 : 0;
        if (eccentricity < 0.05)
                return std::nullopt;

        PointF major = PointF(std::cos(angle), std::sin(angle));
        if (length(major) < kEpsilon)
                return std::nullopt;
        major = normalized(major);

        PointF scanDir = major;
        PointF offsetDir = normalized(PointF(-major.y, major.x));

        return OrientationEstimate{scanDir, offsetDir, eccentricity};
}

} // namespace ZXing::OneD

