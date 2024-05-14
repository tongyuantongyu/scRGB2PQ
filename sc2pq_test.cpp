#include <gtest/gtest.h>
#include <hwy/targets.h>

#include <array>
#include <chrono>
#include <cmath>
#include <numeric>
#include <random>

#include "sc2pq.h"

static float clip(float i) {
  i = i < 0 ? 0 : i;
  i = i > 1 ? 1 : i;
  return i;
}

static std::array<float, 3> bt709_bt2100(float r, float g, float b) {
  return {clip(0.627403895934699f * r + 0.3292830383778837f * g + 0.04331306568741723f * b),
          clip(0.06909728935823207f * r + 0.919540395075459f * g + 0.01136231556630918f * b),
          clip(0.01639143887515028f * r + 0.0880133078772257f * g + 0.895595253247624f * b)};
}

static float toGammaPQ(float linear) {
  if (linear > 0.0f) {
    // Scale from extended SDR range to [0.0, 1.0].
    const float powLinear = powf(linear, 0.1593017578125f);
    const float num = 0.1640625f * powLinear - 0.1640625f;
    const float den = 1.0f + 18.6875f * powLinear;
    return powf(1.0f + num / den, 78.84375f);
  }
  else {
    return 0.0f;
  }
}

static uint16_t quant(float i, float bound) {
  i = clip(i) * bound;
  i = roundf(i);
  return uint16_t(i);
}

constexpr float scale_f2h = 1.925929944387236e-34f;
constexpr float scale_h2f = 5.192296858534828e+33f;

static float h2f(uint16_t v) {
  auto bits = uint32_t(v) << 13;
  auto f = std::bit_cast<float>(bits) * scale_h2f;
  return f;
}

static uint16_t f2h(float v) {
  auto scaled = v * scale_f2h;
  auto bits = std::bit_cast<uint32_t>(scaled);
  return uint16_t(bits >> 13);
}

template<class Derived>
struct TestBase : testing::Test {
  void SetUp() override {
    auto arch = get<0>(static_cast<Derived *>(this)->GetParam());
    if (!(hwy::SupportedTargets() & arch)) {
      GTEST_SKIP() << "Current CPU doesn't support this ISA";
    }
    else {
      hwy::SetSupportedTargetsForTest(arch);
    }
  }

  void TearDown() override { hwy::SetSupportedTargetsForTest(0); }
};

// Value

struct ValueTest : TestBase<ValueTest>,
                   testing::WithParamInterface<std::tuple<int64_t, float, uint32_t, PrecisionMode>> {
  std::vector<uint16_t> values;
  std::vector<uint16_t> results;

  void CompareResult() {
    auto unit = get<1>(GetParam());
    auto depth = uint8_t(get<2>(GetParam()));

    ASSERT_EQ(values.size() % 4, 0);
    ASSERT_EQ(results.size() % 4, 0);
    ASSERT_EQ(values.size(), results.size());

    size_t count = values.size() / 4;
    float scale = unit / 10000.0f;

    for (size_t i = 0; i < count; ++i) {
      auto [r, g, b] =
          bt709_bt2100(h2f(values[4 * i + 0]) * scale, h2f(values[4 * i + 1]) * scale, h2f(values[4 * i + 2]) * scale);
      auto r_ref = quant(toGammaPQ(r), float((1 << depth) - 1));
      auto g_ref = quant(toGammaPQ(g), float((1 << depth) - 1));
      auto b_ref = quant(toGammaPQ(b), float((1 << depth) - 1));
      auto a_ref = quant(h2f(values[4 * i + 3]), float((1 << depth) - 1));

      EXPECT_NEAR(r_ref, results[4 * i + 0], 1) << "Mismatched at R value " << i;
      EXPECT_NEAR(g_ref, results[4 * i + 1], 1) << "Mismatched at G value " << i;
      EXPECT_NEAR(b_ref, results[4 * i + 2], 1) << "Mismatched at B value " << i;
      EXPECT_NEAR(a_ref, results[4 * i + 3], 1) << "Mismatched at A value " << i;
    }
  }

  static std::string Name(const testing::TestParamInfo<ParamType> &info) {
    std::string name;
    name += hwy::TargetName(std::get<0>(info.param));
    name += "_range";
    name += std::to_string(uint32_t(std::get<1>(info.param)));
    name += "_depth";
    name += std::to_string(std::get<2>(info.param));
    name += "_prec";
    switch (std::get<3>(info.param)) {
      case UNorm12Bits: name += "12"; break;
      case UNorm16Bits: name += "16";
    }
    return name;
  }
};

TEST_P(ValueTest, AllGray) {
  auto unit = get<1>(GetParam());
  auto depth = uint8_t(get<2>(GetParam()));
  auto prec = get<3>(GetParam());

  if (depth == 16 && prec == UNorm12Bits) {
    SUCCEED() << "Skip, accuracy target mismatch";
    return;
  }

  auto max_value = f2h(10000.0f / unit);
  auto one = f2h(1.0f);
  for (uint16_t i = 0; i < max_value + 1; ++i) {
    values.push_back(i);
    values.push_back(i);
    values.push_back(i);
    values.push_back(std::min(i, one));
  }

  results.resize(size_t(4 * (max_value + 1)), 0xffff);

  scRGB2PQParams params {values.data(),
                         static_cast<ptrdiff_t>(values.size()),
                         results.data(),
                         static_cast<ptrdiff_t>(results.size()),
                         static_cast<size_t>(max_value + 1),
                         1,
                         unit,
                         depth,
                         prec};

  auto result = scRGB2PQ(&params);

  if (result == -2) {
    GTEST_SKIP() << "Routine not compiled";
  }

  ASSERT_EQ(result, 0);

  CompareResult();
}

TEST_P(ValueTest, RandomValues) {
  auto unit = get<1>(GetParam());
  auto depth = uint8_t(get<2>(GetParam()));
  auto prec = get<3>(GetParam());

  if (depth == 16 && prec == UNorm12Bits) {
    SUCCEED() << "Skip, accuracy target mismatch";
    return;
  }

  size_t test_count = 1048576;

  values.reserve(test_count * 4);
  std::mt19937 rng {0u};
  std::uniform_real_distribution<float> rgb {0.0f, 4000.0f / unit};
  std::uniform_real_distribution<float> a {0.0f, 1.0f};

  for (size_t i = 0; i < test_count; ++i) {
    values.push_back(f2h(rgb(rng)));
    values.push_back(f2h(rgb(rng)));
    values.push_back(f2h(rgb(rng)));
    values.push_back(f2h(a(rng)));
  }

  results.resize(test_count * 4, 0xffff);

  scRGB2PQParams params {values.data(),
                         static_cast<ptrdiff_t>(values.size()),
                         results.data(),
                         static_cast<ptrdiff_t>(results.size()),
                         test_count,
                         1,
                         unit,
                         depth,
                         prec};

  auto result = scRGB2PQ(&params);

  if (result == -2) {
    GTEST_SKIP() << "Routine not compiled";
  }

  ASSERT_EQ(result, 0);
  ASSERT_NEAR(params.max_nits, 3985, 1);
  ASSERT_EQ(params.avg_nits, 2001);

  CompareResult();
}

#if HWY_ARCH_X86_64

INSTANTIATE_TEST_SUITE_P(X86_64, ValueTest,
                         testing::Combine(testing::Values(HWY_SSE2, HWY_AVX2, HWY_AVX3), testing::Values(1.0f, 80.0f),
                                          testing::Values(10, 12, 16), testing::Values(UNorm12Bits, UNorm16Bits)),
                         ValueTest::Name);

#elif HWY_ARCH_ARM

INSTANTIATE_TEST_SUITE_P(ARM, ValueTest,
                         testing::Combine(testing::Values(HWY_EMU128, HWY_NEON), testing::Values(1.0f, 80.0f),
                                          testing::Values(10, 12, 16), testing::Values(UNorm12Bits, UNorm16Bits)),
                         ValueTest::Name);

#endif

// Safety

struct CorrectTest : TestBase<CorrectTest>, testing::WithParamInterface<std::tuple<int64_t, size_t>> {
  static std::string Name(const testing::TestParamInfo<ParamType> &info) {
    std::string name = hwy::TargetName(std::get<0>(info.param));
    name += +"_extra";
    name += std::to_string(std::get<1>(info.param));
    return name;
  }
};

TEST_P(CorrectTest, PreciseWrite) {
  auto extra = std::get<1>(GetParam());
  auto width = 64 + extra;
  size_t pad_width = 64 + 64;
  size_t height = 4;

  std::vector<uint16_t> values(size_t(pad_width * height * 4), 0);
  std::vector<uint16_t> results(size_t(pad_width * height * 4), 0xfefe);

  scRGB2PQParams params {values.data(),  static_cast<ptrdiff_t>(4 * pad_width),
                         results.data(), static_cast<ptrdiff_t>(4 * pad_width),
                         width,          height,
                         80.0f,          10,
                         UNorm12Bits};
  auto result = scRGB2PQ(&params);

  if (result == -2) {
    GTEST_SKIP() << "Routine not compiled";
  }

  ASSERT_EQ(result, 0);

  for (int y = 0; y < height; ++y) {
    for (size_t x = 0; x < width; ++x) {
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 0], 0);
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 1], 0);
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 2], 0);
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 3], 0);
    }

    for (size_t x = width; x < pad_width; ++x) {
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 0], 0xfefe);
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 1], 0xfefe);
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 2], 0xfefe);
      EXPECT_EQ(results[4 * pad_width * y + 4 * x + 3], 0xfefe);
    }
  }
}

#if HWY_ARCH_X86_64

INSTANTIATE_TEST_SUITE_P(X86_64, CorrectTest,
                         testing::Combine(testing::Values(HWY_SSE2, HWY_AVX2, HWY_AVX3),
                                          testing::Range<size_t>(0, 16 + 1)),
                         CorrectTest::Name);

#elif HWY_ARCH_ARM

INSTANTIATE_TEST_SUITE_P(ARM, CorrectTest,
                         testing::Combine(testing::Values(HWY_EMU128, HWY_NEON), testing::Range<size_t>(0, 16 + 1)),
                         CorrectTest::Name);

#endif

// Performance

using hw = std::pair<size_t, size_t>;

#define BenchmarkTest DISABLED_BenchmarkTest

struct BenchmarkTest : TestBase<BenchmarkTest>, testing::WithParamInterface<std::tuple<int64_t, PrecisionMode, hw>> {
  using clock = std::chrono::high_resolution_clock;
  using duration = std::chrono::duration<double, std::chrono::milliseconds::period>;

  std::mt19937 rng {0u};
  std::uniform_real_distribution<float> rgb {0.0f, 10000.0f / 80.0f};
  std::uniform_real_distribution<float> a {0.0f, 1.0f};

  constexpr static std::align_val_t align {32};
  std::unique_ptr<uint16_t[], decltype([](auto *p) { operator delete[](p, align); })> src, dst;

  int n = 100;

  void SetUp() override {

    TestBase::SetUp();

    auto [width, height] = get<2>(GetParam());

    auto padded_width = width + 16;

    // MSVC quirk
    src.reset(static_cast<uint16_t *>(operator new[](padded_width * height * 4 * sizeof(uint16_t), align)));
    dst.reset(static_cast<uint16_t *>(operator new[](padded_width * height * 4 * sizeof(uint16_t), align)));

    for (size_t y = 0; y < height; ++y) {
      for (size_t x = 0; x < padded_width; ++x) {
        src[4 * padded_width * y + 4 * x + 0] = f2h(rgb(rng));
        src[4 * padded_width * y + 4 * x + 1] = f2h(rgb(rng));
        src[4 * padded_width * y + 4 * x + 2] = f2h(rgb(rng));
        src[4 * padded_width * y + 4 * x + 3] = f2h(a(rng));
      }
    }
  }

  static std::string Name(const testing::TestParamInfo<ParamType> &info) {
    std::string name;
    name += hwy::TargetName(std::get<0>(info.param));
    name += "_prec";
    switch (std::get<1>(info.param)) {
      case UNorm12Bits: name += "12"; break;
      case UNorm16Bits: name += "16";
    }
    name += "_size";
    auto [w, h] = std::get<2>(info.param);
    name += std::to_string(w);
    name += "x";
    name += std::to_string(h);
    return name;
  }

  void Measure(uint16_t *psrc, uint16_t *pdst, ptrdiff_t src_stride, ptrdiff_t dst_stride, size_t width,
               PrecisionMode prec) const {
    auto [_, height] = get<2>(GetParam());

    auto begin = clock::now();

    scRGB2PQParams params {psrc, src_stride, pdst, dst_stride, width, height, 80.0f, 10, prec};

    for (int i = 0; i < n; ++i) {
      auto result = scRGB2PQ(&params);

      if (result == -2) {
        GTEST_SKIP() << "Routine not compiled";
      }

      ASSERT_EQ(result, 0);
    }

    auto speed = static_cast<duration>(clock::now() - begin) / n;

    std::cout << Name({GetParam(), 0}) << "\tSpeed: " << speed << "/it" << std::endl;
  }
};

TEST_P(BenchmarkTest, Contiguous) {
  auto [width, height] = get<2>(GetParam());

  Measure(src.get(), dst.get(), static_cast<ptrdiff_t>(width * 4), static_cast<ptrdiff_t>(width * 4), width,
          get<1>(GetParam()));
}

TEST_P(BenchmarkTest, Strided) {
  auto [width, height] = get<2>(GetParam());
  auto padded_width = width + 16;

  Measure(src.get(), dst.get(), static_cast<ptrdiff_t>(padded_width * 4), static_cast<ptrdiff_t>(padded_width * 4),
          width, get<1>(GetParam()));
}

TEST_P(BenchmarkTest, Unaligned) {
  auto [width, height] = get<2>(GetParam());
  auto padded_width = width + 1;

  Measure(src.get(), dst.get(), static_cast<ptrdiff_t>(padded_width * 4), static_cast<ptrdiff_t>(padded_width * 4),
          width, get<1>(GetParam()));
}

TEST_P(BenchmarkTest, Odd) {
  auto [width, height] = get<2>(GetParam());

  Measure(src.get(), dst.get(), static_cast<ptrdiff_t>(width * 4), static_cast<ptrdiff_t>(width * 4), width - 1,
          get<1>(GetParam()));
}

TEST_P(BenchmarkTest, OddUnaligned) {
  auto [width, height] = get<2>(GetParam());
  auto padded_width = width + 1;

  Measure(src.get(), dst.get(), static_cast<ptrdiff_t>(padded_width * 4), static_cast<ptrdiff_t>(padded_width * 4),
          width - 1, get<1>(GetParam()));
}

TEST_P(BenchmarkTest, Flipped) {
  auto [width, height] = get<2>(GetParam());
  auto dst_last = dst.get() + (height - 1) * width * 4;

  Measure(src.get(), dst_last, static_cast<ptrdiff_t>(width * 4), -static_cast<ptrdiff_t>(width * 4), width,
          get<1>(GetParam()));
}

#if HWY_ARCH_X86_64

INSTANTIATE_TEST_SUITE_P(X86_64, BenchmarkTest,
                         testing::Combine(testing::Values(HWY_SSE2, HWY_AVX2, HWY_AVX3),
                                          testing::Values(UNorm12Bits, UNorm16Bits),
                                          testing::Values(hw {1280, 720}, hw {1920, 1080}, hw {3840, 2160})),
                         BenchmarkTest::Name);

#elif HWY_ARCH_ARM

INSTANTIATE_TEST_SUITE_P(ARM, BenchmarkTest,
                         testing::Combine(testing::Values(HWY_EMU128, HWY_NEON),
                                          testing::Values(UNorm12Bits, UNorm16Bits),
                                          testing::Values(hw {1280, 720}, hw {1920, 1080}, hw {3840, 2160})),
                         BenchmarkTest::Name);

#endif