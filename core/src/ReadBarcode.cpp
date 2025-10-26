/*
* Copyright 2019 Axel Waggershauser
*/
// SPDX-License-Identifier: Apache-2.0

#include "ReadBarcode.h"

#if !defined(ZXING_READERS) && !defined(ZXING_WRITERS)
#include "Version.h"
#endif

#ifdef ZXING_READERS
#include "GlobalHistogramBinarizer.h"
#include "HybridBinarizer.h"
#include "oned/ODOrientationUtils.h"
#include "MultiFormatReader.h"
#include "Pattern.h"
#include "ThresholdBinarizer.h"
#include "Point.h"
#endif

#include <climits>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ZXing {

#ifdef ZXING_READERS

class LumImage : public Image
{
public:
        using Image::Image;

        uint8_t* data() { return const_cast<uint8_t*>(Image::data()); }
};

struct RotatedView
{
        LumImage image;
        double cosA = 0;
        double sinA = 0;
        double minX = 0;
        double minY = 0;
        double cx = 0;
        double cy = 0;
};

static RotatedView RotateLuminance(const ImageView& iv, double angle)
{
        RotatedView rot;
        rot.cosA = std::cos(angle);
        rot.sinA = std::sin(angle);
        rot.cx = (iv.width() - 1) / 2.0;
        rot.cy = (iv.height() - 1) / 2.0;

        auto rotatePoint = [&](double x, double y) {
                double rx = rot.cosA * (x - rot.cx) - rot.sinA * (y - rot.cy);
                double ry = rot.sinA * (x - rot.cx) + rot.cosA * (y - rot.cy);
                return std::pair<double, double>{rx, ry};
        };

        auto corners = {rotatePoint(0, 0), rotatePoint(iv.width() - 1, 0),
                                        rotatePoint(iv.width() - 1, iv.height() - 1), rotatePoint(0, iv.height() - 1)};

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

        int outW = static_cast<int>(std::ceil(maxX - minX)) + 1;
        int outH = static_cast<int>(std::ceil(maxY - minY)) + 1;
        rot.minX = minX;
        rot.minY = minY;
        rot.image = LumImage(outW, outH);

        auto* dst = rot.image.data();
        std::fill(dst, dst + outW * outH, uint8_t(0xff));

        for (int y = 0; y < outH; ++y) {
                double ry = minY + y;
                for (int x = 0; x < outW; ++x) {
                        double rx = minX + x;
                        double srcX = rot.cosA * rx + rot.sinA * ry + rot.cx;
                        double srcY = -rot.sinA * rx + rot.cosA * ry + rot.cy;
                        if (srcX < 0 || srcX > iv.width() - 1 || srcY < 0 || srcY > iv.height() - 1)
                                continue;

                        int x0 = static_cast<int>(std::floor(srcX));
                        int y0 = static_cast<int>(std::floor(srcY));
                        int x1 = std::min(x0 + 1, iv.width() - 1);
                        int y1 = std::min(y0 + 1, iv.height() - 1);

                        double fx = srcX - x0;
                        double fy = srcY - y0;

                        auto sample = [&](int sx, int sy) {
                                return static_cast<double>(iv.data(sx, sy)[0]);
                        };

                        double v00 = sample(x0, y0);
                        double v10 = sample(x1, y0);
                        double v01 = sample(x0, y1);
                        double v11 = sample(x1, y1);

                        double interp = (1 - fx) * (1 - fy) * v00 + fx * (1 - fy) * v10 + (1 - fx) * fy * v01 + fx * fy * v11;
                        dst[y * outW + x] = static_cast<uint8_t>(std::clamp(std::lround(interp), 0l, 255l));
                }
        }

        return rot;
}

static PointF MapRotatedPoint(const RotatedView& rot, PointF p)
{
        double rx = rot.minX + p.x;
        double ry = rot.minY + p.y;
        double x = rot.cosA * rx + rot.sinA * ry + rot.cx;
        double y = -rot.sinA * rx + rot.cosA * ry + rot.cy;
        return {x, y};
}

static bool IsLikelyDuplicateLinear(const Barcode& existing, const Barcode& candidate)
{
        if (!IsLinearBarcode(existing.format()) || !IsLinearBarcode(candidate.format()))
                return false;

        if (existing.format() != candidate.format() || existing.text() != candidate.text())
                return false;

        auto orientationDiff = [](int a, int b) {
                int diff = std::abs(a - b) % 360;
                if (diff > 180)
                        diff = 360 - diff;
                if (diff > 90)
                        diff = 180 - diff;
                return diff;
        }(existing.orientation(), candidate.orientation());

        constexpr int kMaxOrientationDiff = 10;
        if (orientationDiff > kMaxOrientationDiff)
                return false;

        auto centerOf = [](const QuadrilateralI& pos) {
                PointF acc(0.0, 0.0);
                for (const auto& pt : pos)
                        acc += PointF(pt.x, pt.y);
                return acc / 4.0;
        };

        constexpr double kMaxCenterDistance = 40.0;
        return distance(centerOf(existing.position()), centerOf(candidate.position())) <= kMaxCenterDistance;
}

static bool HasLinearFormat(const ReaderOptions& opts)
{
        auto formats = opts.formats();
        if (formats.empty())
                return true;
        for (auto format : formats) {
                if (IsLinearBarcode(format))
                        return true;
        }
        return false;
}

template<typename P>
static LumImage ExtractLum(const ImageView& iv, P projection)
{
	LumImage res(iv.width(), iv.height());

	auto* dst = res.data();
	for(int y = 0; y < iv.height(); ++y)
		for(int x = 0, w = iv.width(); x < w; ++x)
			*dst++ = projection(iv.data(x, y));

	return res;
}

class LumImagePyramid
{
	std::vector<LumImage> buffers;

	template<int N>
	void addLayer()
	{
		auto siv = layers.back();
		buffers.emplace_back(siv.width() / N, siv.height() / N);
		layers.push_back(buffers.back());
		auto& div = buffers.back();
		auto* d   = div.data();

		for (int dy = 0; dy < div.height(); ++dy)
			for (int dx = 0; dx < div.width(); ++dx) {
				int sum = (N * N) / 2;
				for (int ty = 0; ty < N; ++ty)
					for (int tx = 0; tx < N; ++tx)
						sum += *siv.data(dx * N + tx, dy * N + ty);
				*d++ = sum / (N * N);
			}
	}

	void addLayer(int factor)
	{
		// help the compiler's auto-vectorizer by hard-coding the scale factor
		switch (factor) {
		case 2: addLayer<2>(); break;
		case 3: addLayer<3>(); break;
		case 4: addLayer<4>(); break;
		default: throw std::invalid_argument("Invalid ReaderOptions::downscaleFactor"); break;
		}
	}

public:
	std::vector<ImageView> layers;

	LumImagePyramid(const ImageView& iv, int threshold, int factor)
	{
		if (factor < 2)
			throw std::invalid_argument("Invalid ReaderOptions::downscaleFactor");

		layers.push_back(iv);
		// TODO: if only matrix codes were considered, then using std::min would be sufficient (see #425)
		while (threshold > 0 && std::max(layers.back().width(), layers.back().height()) > threshold &&
			   std::min(layers.back().width(), layers.back().height()) >= factor)
			addLayer(factor);
#if 0
		// Reversing the layers means we'd start with the smallest. that can make sense if we are only looking for a
		// single symbol. If we start with the higher resolution, we get better (high res) position information.
		// TODO: see if masking out higher res layers based on found symbols in lower res helps overall performance.
		std::reverse(layers.begin(), layers.end());
#endif
	}
};

ImageView SetupLumImageView(ImageView iv, LumImage& lum, const ReaderOptions& opts)
{
	if (iv.format() == ImageFormat::None)
		throw std::invalid_argument("Invalid image format");

	if (opts.binarizer() == Binarizer::GlobalHistogram || opts.binarizer() == Binarizer::LocalAverage) {
		// manually spell out the 3 most common pixel formats to get at least gcc to vectorize the code
		if (iv.format() == ImageFormat::RGB && iv.pixStride() == 3) {
			lum = ExtractLum(iv, [](const uint8_t* src) { return RGBToLum(src[0], src[1], src[2]); });
		} else if (iv.format() == ImageFormat::RGBA && iv.pixStride() == 4) {
			lum = ExtractLum(iv, [](const uint8_t* src) { return RGBToLum(src[0], src[1], src[2]); });
		} else if (iv.format() == ImageFormat::BGR && iv.pixStride() == 3) {
			lum = ExtractLum(iv, [](const uint8_t* src) { return RGBToLum(src[2], src[1], src[0]); });
		} else if (iv.format() != ImageFormat::Lum) {
			lum = ExtractLum(iv, [r = RedIndex(iv.format()), g = GreenIndex(iv.format()), b = BlueIndex(iv.format())](
									 const uint8_t* src) { return RGBToLum(src[r], src[g], src[b]); });
		} else if (iv.pixStride() != 1) {
			// GlobalHistogram and LocalAverage need dense line memory layout
			lum = ExtractLum(iv, [](const uint8_t* src) { return *src; });
		}
		if (lum.data())
			return lum;
	}
	return iv;
}

std::unique_ptr<BinaryBitmap> CreateBitmap(ZXing::Binarizer binarizer, const ImageView& iv)
{
	switch (binarizer) {
	case Binarizer::BoolCast: return std::make_unique<ThresholdBinarizer>(iv, 0);
	case Binarizer::FixedThreshold: return std::make_unique<ThresholdBinarizer>(iv, 127);
	case Binarizer::GlobalHistogram: return std::make_unique<GlobalHistogramBinarizer>(iv);
	case Binarizer::LocalAverage: return std::make_unique<HybridBinarizer>(iv);
	}
	return {}; // silence gcc warning
}

Barcode ReadBarcode(const ImageView& _iv, const ReaderOptions& opts)
{
	return FirstOrDefault(ReadBarcodes(_iv, ReaderOptions(opts).setMaxNumberOfSymbols(1)));
}

Barcodes ReadBarcodes(const ImageView& _iv, const ReaderOptions& opts)
{
	if (sizeof(PatternType) < 4 && (_iv.width() > 0xffff || _iv.height() > 0xffff))
		throw std::invalid_argument("Maximum image width/height is 65535");

	if (!_iv.data() || _iv.width() * _iv.height() == 0)
		throw std::invalid_argument("ImageView is null/empty");

	LumImage lum;
	ImageView iv = SetupLumImageView(_iv, lum, opts);
	MultiFormatReader reader(opts);

	if (opts.isPure())
		return {reader.read(*CreateBitmap(opts.binarizer(), iv)).setReaderOptions(opts)};

	std::unique_ptr<MultiFormatReader> closedReader;
#ifdef ZXING_EXPERIMENTAL_API
	auto formatsBenefittingFromClosing = BarcodeFormat::Aztec | BarcodeFormat::DataMatrix | BarcodeFormat::QRCode | BarcodeFormat::MicroQRCode;
	ReaderOptions closedOptions = opts;
	if (opts.tryDenoise() && opts.hasFormat(formatsBenefittingFromClosing) && _iv.height() >= 3) {
		closedOptions.setFormats((opts.formats().empty() ? BarcodeFormat::Any : opts.formats()) & formatsBenefittingFromClosing);
		closedReader = std::make_unique<MultiFormatReader>(closedOptions);
	}
#endif
	LumImagePyramid pyramid(iv, opts.downscaleThreshold() * opts.tryDownscale(), opts.downscaleFactor());

        Barcodes res;
        auto isDuplicateLinear = [&](const Barcode& candidate) {
                for (const auto& existing : res) {
                        if (IsLikelyDuplicateLinear(existing, candidate))
                                return true;
                }
                return false;
        };
        int maxSymbols = opts.maxNumberOfSymbols() ? opts.maxNumberOfSymbols() : INT_MAX;
        std::shared_ptr<const BitMatrix> orientationBits;
	for (auto&& iv : pyramid.layers) {
		auto bitmap = CreateBitmap(opts.binarizer(), iv);
		for (int close = 0; close <= (closedReader ? 1 : 0); ++close) {
			if (close) {
				// if we already inverted the image in the first round, we need to undo that first
				if (bitmap->inverted())
					bitmap->invert();
				bitmap->close();
			}

			// TODO: check if closing after invert would be beneficial
                        for (int invert = 0; invert <= static_cast<int>(opts.tryInvert() && !close); ++invert) {
                                if (invert)
                                        bitmap->invert();
                                if (!orientationBits && &iv == &pyramid.layers.front() && close == 0 && invert == 0) {
                                        if (const BitMatrix* bits = bitmap->getBitMatrix())
                                                orientationBits = std::make_shared<BitMatrix>(bits->copy());
                                }
                                auto rs = (close ? *closedReader : reader).readMultiple(*bitmap, maxSymbols);
				for (auto& r : rs) {
					if (iv.width() != _iv.width())
						r.setPosition(Scale(r.position(), _iv.width() / iv.width()));
                                        if (!Contains(res, r) && !isDuplicateLinear(r)) {
                                                r.setReaderOptions(opts);
                                                r.setIsInverted(bitmap->inverted());
                                                res.push_back(std::move(r));
                                                --maxSymbols;
                                        }
				}
				if (maxSymbols <= 0)
					return res;
			}
		}
        }

        if (res.empty() && orientationBits && HasLinearFormat(opts) && opts.tryRotate()) {
                if (auto orientation = OneD::EstimateOrientation(*orientationBits)) {
                        double axisAlignment = std::max(std::abs(dot(orientation->scanDir, PointF{1, 0})),
                                                        std::abs(dot(orientation->scanDir, PointF{0, 1})));
                        if (axisAlignment < 0.98) {
                                double angle = std::atan2(orientation->scanDir.y, orientation->scanDir.x);
                                auto rotated = RotateLuminance(pyramid.layers.front(), -angle);
                                ReaderOptions rotatedOpts(opts);
                                rotatedOpts.setTryRotate(false);
                                if (opts.maxNumberOfSymbols())
                                        rotatedOpts.setMaxNumberOfSymbols(static_cast<uint8_t>(std::clamp(maxSymbols, 0, 255)));
                                auto rotatedResults = ReadBarcodes(rotated.image, rotatedOpts);
                                for (auto& r : rotatedResults) {
                                        auto pos = r.position();
                                        for (auto& pt : pos) {
                                                PointF mapped = MapRotatedPoint(rotated, PointF(pt.x, pt.y));
                                                pt = {static_cast<int>(std::lround(mapped.x)),
                                                      static_cast<int>(std::lround(mapped.y))};
                                        }
                                        r.setPosition(std::move(pos));
                                        r.setReaderOptions(opts);
                                        if (!Contains(res, r) && !isDuplicateLinear(r)) {
                                                r.setReaderOptions(opts);
                                                res.push_back(std::move(r));
                                                if (--maxSymbols <= 0)
                                                        break;
                                        }
                                }
                        }
                }
        }

        if (res.size() < static_cast<size_t>(maxSymbols) && HasLinearFormat(opts) && opts.tryRotate()) {
                const std::vector<int> coarseAngles = {-75, -72, -69, -66, -63, -60, -57, -54, -51, -48, -45, -42, -39,
                                                      -36, -33, -30, -27, -24, -21, -18, -15, -12, -9, -6, -3, 3, 6, 9,
                                                      12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57,
                                                      60, 63, 66, 69, 72, 75};

                for (int deg : coarseAngles) {
                        if (maxSymbols <= 0)
                                break;

                        double angle = deg * M_PI / 180.0;
                        auto rotated = RotateLuminance(pyramid.layers.front(), -angle);
                        ReaderOptions rotatedOpts(opts);
                        rotatedOpts.setTryRotate(false);
                        if (opts.maxNumberOfSymbols())
                                rotatedOpts.setMaxNumberOfSymbols(static_cast<uint8_t>(std::clamp(maxSymbols, 0, 255)));

                        auto rotatedResults = ReadBarcodes(rotated.image, rotatedOpts);
                        for (auto& r : rotatedResults) {
                                auto pos = r.position();
                                for (auto& pt : pos) {
                                        PointF mapped = MapRotatedPoint(rotated, PointF(pt.x, pt.y));
                                        pt = {static_cast<int>(std::lround(mapped.x)),
                                              static_cast<int>(std::lround(mapped.y))};
                                }
                                r.setPosition(std::move(pos));
                                r.setReaderOptions(opts);
                                if (!Contains(res, r) && !isDuplicateLinear(r)) {
                                        res.push_back(std::move(r));
                                        if (--maxSymbols <= 0)
                                                break;
                                }
                        }
                }
        }

        return res;
}

#else // ZXING_READERS

Barcode ReadBarcode(const ImageView&, const ReaderOptions&)
{
	throw std::runtime_error("This build of zxing-cpp does not support reading barcodes.");
}

Barcodes ReadBarcodes(const ImageView&, const ReaderOptions&)
{
	throw std::runtime_error("This build of zxing-cpp does not support reading barcodes.");
}

#endif // ZXING_READERS

} // ZXing
