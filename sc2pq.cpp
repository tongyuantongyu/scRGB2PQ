#undef HWY_TARGET_INCLUDE
#define HWY_TARGET_INCLUDE "sc2pq.cpp"// this file

#include <hwy/detect_compiler_arch.h>

#if HWY_ARCH_X86_64
#define HWY_BASELINE_TARGETS HWY_SSE2
#else
#define HWY_BASELINE_TARGETS HWY_EMU128
#endif

#define HWY_DISABLED_TARGETS ~(HWY_AVX2 | HWY_AVX3 | HWY_SSE2 | HWY_NEON | HWY_EMU128)

#include <array>
#include <cstring>

#include <hwy/foreach_target.h>// must come before highway.h

#include <hwy/contrib/math/math-inl.h>
#include <hwy/highway.h>

#include "sc2pq.h"

#define SC2PQ_INTRIN_FUNC HWY_ATTR HWY_INLINE

namespace hwy {
namespace HWY_NAMESPACE {

template<class D, class V>
SC2PQ_INTRIN_FUNC V LogR(const D d, V x) {
  // http://git.musl-libc.org/cgit/musl/tree/src/math/log.c for more info.
  using T = TFromD<D>;
  impl::LogImpl<T> impl;

  constexpr bool kIsF32 = (sizeof(T) == 4);

  // Float Constants
  const V kLn2Hi = Set(d, kIsF32 ? static_cast<T>(0.69313812256f) : static_cast<T>(0.693147180369123816490));
  const V kLn2Lo = Set(d, kIsF32 ? static_cast<T>(9.0580006145e-6f) : static_cast<T>(1.90821492927058770002e-10));
  const V kOne = Set(d, static_cast<T>(+1.0));
  const V kMinNormal = Set(d, kIsF32 ? static_cast<T>(1.175494351e-38f) : static_cast<T>(2.2250738585072014e-308));
  const V kScale = Set(d, kIsF32 ? static_cast<T>(3.355443200e+7f) : static_cast<T>(1.8014398509481984e+16));

  // Integer Constants
  using TI = MakeSigned<T>;
  const Rebind<TI, D> di;
  using VI = decltype(Zero(di));
  const VI kLowerBits = Set(di, kIsF32 ? static_cast<TI>(0x00000000L) : static_cast<TI>(0xFFFFFFFFLL));
  const VI kMagic = Set(di, kIsF32 ? static_cast<TI>(0x3F3504F3L) : static_cast<TI>(0x3FE6A09E00000000LL));
  const VI kExpMask = Set(di, kIsF32 ? static_cast<TI>(0x3F800000L) : static_cast<TI>(0x3FF0000000000000LL));
  const VI kExpScale = Set(di, kIsF32 ? static_cast<TI>(-25) : static_cast<TI>(-54));
  const VI kManMask = Set(di, kIsF32 ? static_cast<TI>(0x7FFFFFL) : static_cast<TI>(0xFFFFF00000000LL));

  VI exp_bits = Add(BitCast(di, x), Sub(kExpMask, kMagic));
  V exp = ConvertTo(d, impl.Log2p1NoSubnormal(d, BitCast(d, exp_bits)));

  // Renormalize.
  const V y = Or(And(x, BitCast(d, kLowerBits)), BitCast(d, Add(And(exp_bits, kManMask), kMagic)));

  // Approximate and reconstruct.
  const V ym1 = Sub(y, kOne);
  const V z = Div(ym1, Add(y, kOne));

  return MulSub(exp, kLn2Hi, Sub(MulSub(z, Sub(ym1, impl.LogPoly(d, z)), Mul(exp, kLn2Lo)), ym1));
}

template<class D, class V>
SC2PQ_INTRIN_FUNC V ExpR(const D d, V x) {
  using T = TFromD<D>;

  const V kNegHalf = Set(d, static_cast<T>(-0.5));
  const V kNegZero = Set(d, static_cast<T>(-0.0));
  const V kOne = Set(d, static_cast<T>(+1.0));
  const V kOneOverLog2 = Set(d, static_cast<T>(+1.442695040888963407359924681));

  impl::ExpImpl<T> impl;

  // q = static_cast<int32>((x / log(2)) + ((x < 0) ? -0.5 : +0.5))
  const auto q = impl.ToInt32(d, MulAdd(x, kOneOverLog2, kNegHalf));

  // Reduce, approximate, and then reconstruct.
  const V y = impl.LoadExpShortRange(d, Add(impl.ExpPoly(d, impl.ExpReduce(d, x, q)), kOne), q);
  return y;
}

}// namespace HWY_NAMESPACE
}// namespace hwy

namespace {
namespace HWY_NAMESPACE {

namespace hn = hwy::HWY_NAMESPACE;

template<class T>
SC2PQ_INTRIN_FUNC auto PerceptualQuantizationPedantic(T x) {
  hn::ScalableTag<float> d;

  auto pow_linear =
      IfThenElseZero(hn::Ge(x, hn::Set(d, 1.596992334170612e-13)), ExpR(d, hn::Set(d, 0.1593017578125f) * LogR(d, x)));
  auto num = hn::MulAdd(hn::Set(d, 0.1640625f), pow_linear, hn::Set(d, -0.1640625f));
  auto den = hn::MulAdd(hn::Set(d, 18.6875f), pow_linear, hn::Set(d, 1.0f));
  auto fra = hn::MulAdd(num, hn::ApproximateReciprocal(den), hn::Set(d, 1.0f));
  return ExpR(d, hn::Set(d, 78.84375f) * LogR(d, fra));
}

template<class T>
SC2PQ_INTRIN_FUNC auto PerceptualQuantizationD12(T x) {
  // Error < 0.5/4095 for [0, 1] input: suitable for 12bit data

  // Most time this should be sufficient
  hn::ScalableTag<float> d;
  auto log = LogR(d, x);
  // 3,2 order Pade approximation of Log(PQ(Exp(x))) at -4.15
  auto num = hn::impl::Estrin(log, hn::Set(d, 0.005396294078783058f), hn::Set(d, 78.29447758409737f),
                              hn::Set(d, -1.476062552035652f), hn::Set(d, 0.04729362818571236f));
  auto den = hn::impl::Estrin(log, hn::Set(d, 747.9678350466744f), hn::Set(d, 39.87717455527674f), hn::Set(d, 1.0f));
  auto fra = hn::Div(num, den);
  return IfThenElseZero(hn::Ge(x, hn::Set(d, 1.596992334170612e-13)), ExpR(d, fra));
}

template<class T>
SC2PQ_INTRIN_FUNC auto PerceptualQuantizationD16(T x) {
  // Error < 0.5/65535 for [0, 1] input: suitable for 16bit data

  hn::ScalableTag<float> d;
  auto log = LogR(d, x);
  // 3,3 order Pade approximation of Log(PQ(Exp(x))) at -6.3
  auto num = hn::impl::Estrin(log, hn::Set(d, 0.4048046902710638f), hn::Set(d, 6037.836117036402f),
                              hn::Set(d, -47.289523694165744f), hn::Set(d, 2.9584270989160584f));
  auto den = hn::impl::Estrin(log, hn::Set(d, 57682.5383969059f), hn::Set(d, 3712.4393477306107f),
                              hn::Set(d, 117.07215624733027f), hn::Set(d, 1.0f));
  auto fra = hn::Div(num, den);
  return IfThenElseZero(hn::Ge(x, hn::Set(d, 1.596992334170612e-13)), ExpR(d, fra));
}

template<class T>
SC2PQ_INTRIN_FUNC auto BT709ToBT2100Clip(T r, T g, T b, T a) {
  hn::ScalableTag<float> d;
  auto ro =
      hn::MulAdd(r, hn::Set(d, 0.627409f), hn::MulAdd(g, hn::Set(d, 0.0691248f), hn::Mul(b, hn::Set(d, 0.0164234f))));
  auto go =
      hn::MulAdd(g, hn::Set(d, 0.919549f), hn::MulAdd(r, hn::Set(d, 0.32926f), hn::Mul(b, hn::Set(d, 0.0880478f))));
  auto bo =
      hn::MulAdd(b, hn::Set(d, 0.895617f), hn::MulAdd(r, hn::Set(d, 0.0432719f), hn::Mul(g, hn::Set(d, 0.0113208f))));

  auto zero = hn::Zero(d);
  auto one = hn::Set(d, 1.0f);

  auto rc = hn::Clamp(ro, zero, one);
  auto gc = hn::Clamp(go, zero, one);
  auto bc = hn::Clamp(bo, zero, one);
  auto ac = hn::Clamp(a, zero, one);

  return std::array {rc, gc, bc, ac};
}

template<PrecisionMode mode, class T>
SC2PQ_INTRIN_FUNC auto PerceptualQuantization(T x) {
  if constexpr (mode == PrecisionMode::UNorm12Bits) {
    return PerceptualQuantizationD12(x);
  }
  else if constexpr (mode == PrecisionMode::UNorm16Bits) {
    return PerceptualQuantizationD16(x);
  }
  else {
    static_assert(false, "Unknown precision mode");
  }
}

template<class T>
SC2PQ_INTRIN_FUNC auto QuantNb(T x, float up) {
  auto rx = hn::Round(Mul(x, hn::Set(hn::DFromV<T> {}, up)));
  return hn::BitCast(hn::Rebind<uint32_t, hn::DFromV<T>>(), hn::ConvertTo(hn::Rebind<int32_t, hn::DFromV<T>>(), rx));
}

#if HWY_TARGET == HWY_SSE2 || HWY_TARGET == HWY_AVX2 || HWY_TARGET == HWY_AVX3

SC2PQ_INTRIN_FUNC auto LoadRGBAF16(const uint16_t *src, float unit_scale) {
  hn::ScalableTag<uint16_t> u16;
  hn::ScalableTag<uint32_t> u32;
  hn::ScalableTag<float> f32;

  const auto x1 = hn::LoadU(u16, src);
  const auto x2 = hn::LoadU(u16, src + hn::Lanes(u16));

  auto even = hn::BitCast(u32, hn::InterleaveLower(u16, x1, x2));
  auto odd = hn::BitCast(u32, hn::InterleaveUpper(u16, x1, x2));

  auto rgu = hn::BitCast(u16, hn::InterleaveLower(u32, even, odd));
  auto bau = hn::BitCast(u16, hn::InterleaveUpper(u32, even, odd));

  auto zero = hn::Zero(u16);
  auto ru = hn::ShiftLeft<13>(hn::BitCast(u32, hn::InterleaveLower(u16, rgu, zero)));
  auto gu = hn::ShiftLeft<13>(hn::BitCast(u32, hn::InterleaveUpper(u16, rgu, zero)));
  auto bu = hn::ShiftLeft<13>(hn::BitCast(u32, hn::InterleaveLower(u16, bau, zero)));
  auto au = hn::ShiftLeft<13>(hn::BitCast(u32, hn::InterleaveUpper(u16, bau, zero)));

  auto scale = hn::Set(f32, 5.192296858534828e+33f * unit_scale);
  auto rf = hn::Mul(hn::BitCast(f32, ru), scale);
  auto gf = hn::Mul(hn::BitCast(f32, gu), scale);
  auto bf = hn::Mul(hn::BitCast(f32, bu), scale);

  auto a_scale = hn::Set(f32, 5.192296858534828e+33f);
  auto af = hn::Mul(hn::BitCast(f32, au), a_scale);

  return std::array {rf, gf, bf, af};
}

template<class T>
SC2PQ_INTRIN_FUNC auto StoreRGBAU16(uint16_t *dst, T r, T g, T b, T a) {
  hn::ScalableTag<uint16_t> u16;
  hn::ScalableTag<uint32_t> u32;
  hn::ScalableTag<uint64_t> u64;

  auto rg = hn::Add(r, hn::ShiftLeft<16>(g));
  auto ba = hn::Add(b, hn::ShiftLeft<16>(a));

  auto even = hn::BitCast(u64, hn::InterleaveLower(u32, rg, ba));
  auto odd = hn::BitCast(u64, hn::InterleaveUpper(u32, rg, ba));

  auto y1 = hn::BitCast(u16, hn::InterleaveLower(u64, even, odd));
  auto y2 = hn::BitCast(u16, hn::InterleaveUpper(u64, even, odd));

  hn::StoreU(y1, u16, dst);
  hn::StoreU(y2, u16, dst + hn::Lanes(u16));
}

#endif

#if HWY_TARGET == HWY_EMU128 || HWY_TARGET == HWY_NEON

SC2PQ_INTRIN_FUNC auto LoadRGBAF16(const uint16_t *src, float unit_scale) {
  hn::ScalableTag<uint16_t>::Half u16_h;
  hn::ScalableTag<hwy::float16_t>::Half f16_h;
  hn::ScalableTag<float> f32;

  hn::VFromD<decltype(u16_h)> r, g, b, a;
  hn::LoadInterleaved4(u16_h, src, r, g, b, a);

  auto scale = hn::Set(f32, unit_scale);

  auto r1 = hn::Mul(hn::PromoteTo(f32, hn::BitCast(f16_h, r)), scale);
  auto g1 = hn::Mul(hn::PromoteTo(f32, hn::BitCast(f16_h, g)), scale);
  auto b1 = hn::Mul(hn::PromoteTo(f32, hn::BitCast(f16_h, b)), scale);
  auto a1 = hn::PromoteTo(f32, hn::BitCast(f16_h, a));

  return std::array {r1, g1, b1, a1};
}

template<class T>
SC2PQ_INTRIN_FUNC auto StoreRGBAU16(uint16_t *dst, T r, T g, T b, T a) {
  hn::ScalableTag<uint16_t>::Half u16_h;
  hn::StoreInterleaved4(hn::DemoteTo(u16_h, r), hn::DemoteTo(u16_h, g), hn::DemoteTo(u16_h, b), hn::DemoteTo(u16_h, a),
                        u16_h, dst);
}

#endif

template<PrecisionMode mode>
SC2PQ_INTRIN_FUNC void scRGB2PQOnce(const uint16_t *src, uint16_t *dst, scRGB2PQParams params) {
  const auto [r, g, b, a] = LoadRGBAF16(src, params.unit_scale);
  const auto [rc, gc, bc, ac] = BT709ToBT2100Clip(r, g, b, a);
  const auto rq = QuantNb(PerceptualQuantization<mode>(rc), params.bound);
  const auto gq = QuantNb(PerceptualQuantization<mode>(gc), params.bound);
  const auto bq = QuantNb(PerceptualQuantization<mode>(bc), params.bound);
  const auto aq = QuantNb(ac, params.bound);
  StoreRGBAU16(dst, rq, gq, bq, aq);
};

template<PrecisionMode mode>
HWY_ATTR void scRGB2PQImpl(scRGB2PQParams params) {
  hn::ScalableTag<uint16_t> u16;

  for (size_t y = 0; y < params.height; ++y) {
    const uint16_t *src = params.src + y * params.src_stride;
    uint16_t *dst = params.dst + y * params.dst_stride;
    size_t w = 4 * params.width;

    constexpr size_t buf1 = hn::Lanes(u16);
    constexpr size_t step = 2 * buf1;

    size_t i = 0;
    for (; i + step <= w; i += step) {
      scRGB2PQOnce<mode>(src + i, dst + i, params);
    }

    if (i != w) {
      struct alignas(2 * buf1 * sizeof(uint16_t)) buf {
        uint16_t d[2 * buf1];
      };

      buf s_buf {}, d_buf {};
      memcpy(s_buf.d, src + i, (w - i) * sizeof(uint16_t));
      scRGB2PQOnce<mode>(s_buf.d, d_buf.d, params);
      memcpy(dst + i, d_buf.d, (w - i) * sizeof(uint16_t));
    }
  }
}

int scRGB2PQDispatch(scRGB2PQParams params) {
  switch (params.precision) {
    case UNorm12Bits: scRGB2PQImpl<UNorm12Bits>(params); return 0;
#if !defined(SC2PQ_LEAN)
    case UNorm16Bits: scRGB2PQImpl<UNorm16Bits>(params); return 0;
#endif
    default: return -2;
  }
}

}// namespace HWY_NAMESPACE
}// namespace

#if HWY_ONCE

namespace {

HWY_EXPORT(scRGB2PQDispatch);

}

int scRGB2PQ(scRGB2PQParams params) {
  if (params.src_stride == params.dst_stride && params.src_stride == params.width * 4) {
    params.src_stride = params.dst_stride = 0;
    params.width *= params.height;
    params.height = 1;
  }

  return HWY_DYNAMIC_DISPATCH(scRGB2PQDispatch)(params);
}

#endif
