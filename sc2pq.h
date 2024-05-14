#ifndef NUM_EXP_SC2PQ_H
#define NUM_EXP_SC2PQ_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum PrecisionMode {
  UNorm12Bits,// Suitable for 12bit quantization
  UNorm16Bits,// Suitable for 16bit quantization
};

typedef struct scRGB2PQParams {
  const uint16_t *src;
  ptrdiff_t src_stride;
  uint16_t *dst;
  ptrdiff_t dst_stride;

  size_t width, height;

  float unit_scale;
  float bound;
  PrecisionMode precision;

#ifdef __cplusplus
  scRGB2PQParams(const uint16_t *src, ptrdiff_t src_stride, uint16_t *dst, ptrdiff_t dst_stride, size_t width,
                 size_t height, float unit = 80.0f, uint8_t depth = 10, PrecisionMode precision = UNorm12Bits)
      : src(src), src_stride(src_stride), dst(dst), dst_stride(dst_stride), width(width), height(height),
        unit_scale(unit / 10000.0f), bound(float((1 << depth) - 1)), precision(precision) {}
#endif
} scRGB2PQParams;

inline void scRGB2PQParamsSetUpScaleDepth(scRGB2PQParams *params, float unit, uint8_t depth) {
  params->unit_scale = unit / 10000.0f;
  params->bound = float((1 << depth) - 1);
}

inline void scRGB2PQParamsDefaultScaleDepth(scRGB2PQParams *params) {
  scRGB2PQParamsSetUpScaleDepth(params, 80.0f, 10);
}

// Return 0 for success, otherwise failure
int scRGB2PQ(scRGB2PQParams params);

#ifdef __cplusplus
};
#endif

#endif//NUM_EXP_SC2PQ_H
