#include "ReadBarcode.h"
#include "ReaderOptions.h"
#include "BitMatrix.h"
#include "Matrix.h"
#include "ImageView.h"
#include "oned/ODCode128Writer.h"

#include "gtest/gtest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

using namespace ZXing;

namespace {
std::vector<uint8_t> RotateImage(const std::vector<uint8_t>& src, int width, int height, double angleDeg, int& outW, int& outH)
{
        const double radians = angleDeg * M_PI / 180.0;
        const double cosA = std::cos(radians);
        const double sinA = std::sin(radians);
        const double cx = (width - 1) / 2.0;
        const double cy = (height - 1) / 2.0;

        auto rotatePoint = [&](double x, double y) {
                double rx = cosA * (x - cx) - sinA * (y - cy);
                double ry = sinA * (x - cx) + cosA * (y - cy);
                return std::pair<double, double>{rx, ry};
        };

        std::array<std::pair<double, double>, 4> corners = {
                rotatePoint(0, 0), rotatePoint(width - 1, 0), rotatePoint(width - 1, height - 1), rotatePoint(0, height - 1)};

        double minX = std::numeric_limits<double>::infinity();
        double maxX = -std::numeric_limits<double>::infinity();
        double minY = std::numeric_limits<double>::infinity();
        double maxY = -std::numeric_limits<double>::infinity();

        for (const auto& [rx, ry] : corners) {
                minX = std::min(minX, rx);
                maxX = std::max(maxX, rx);
                minY = std::min(minY, ry);
                maxY = std::max(maxY, ry);
        }

        outW = static_cast<int>(std::ceil(maxX - minX)) + 1;
        outH = static_cast<int>(std::ceil(maxY - minY)) + 1;

        std::vector<uint8_t> dst(outW * outH, 0xff);

        for (int y = 0; y < outH; ++y) {
                double ry = minY + y;
                for (int x = 0; x < outW; ++x) {
                        double rx = minX + x;
                        double srcX = cosA * rx + sinA * ry + cx;
                        double srcY = -sinA * rx + cosA * ry + cy;
                        if (0 <= srcX && srcX < width && 0 <= srcY && srcY < height) {
                                int sx = static_cast<int>(std::round(srcX));
                                int sy = static_cast<int>(std::round(srcY));
                                sx = std::clamp(sx, 0, width - 1);
                                sy = std::clamp(sy, 0, height - 1);
                                dst[y * outW + x] = src[sy * width + sx];
                        }
                }
        }

        return dst;
}
}

TEST(RotatedOneDTest, Code128ArbitraryAngles)
{
        OneD::Code128Writer writer;
        auto bits = writer.encode("0123456789", 0, 80);
        auto matrix = ToMatrix<uint8_t>(bits, uint8_t(0), uint8_t(255));
        std::vector<uint8_t> buffer(matrix.begin(), matrix.end());
        int width = matrix.width();
        int height = matrix.height();

        auto runAngle = [&](double angle) {
                int rotatedW = 0;
                int rotatedH = 0;
                auto rotated = RotateImage(buffer, width, height, angle, rotatedW, rotatedH);
                ImageView view(rotated.data(), rotatedW, rotatedH, ImageFormat::Lum, rotatedW, 1);

                ReaderOptions opts;
                opts.setFormats(BarcodeFormat::Code128);
                opts.setMinLineCount(1);
                auto result = ReadBarcode(view, opts);

                ASSERT_TRUE(result.isValid());
                EXPECT_EQ(result.text(), "0123456789");

                int orientation = result.orientation();
                int roundedAngle = static_cast<int>(std::lround(angle));
                int diff = std::min({std::abs(orientation - roundedAngle),
                                                        std::abs(orientation - (roundedAngle + 180)),
                                                        std::abs(orientation - (roundedAngle - 180))});
                EXPECT_LE(diff, 40);
        };

        runAngle(27.0);
        runAngle(-33.0);
}

