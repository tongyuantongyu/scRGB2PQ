# scRGB to PQ Converter

This library does one thing: convert pixels [scRGB](https://en.wikipedia.org/wiki/ScRGB) stored as
[Half float](https://en.wikipedia.org/wiki/Half-precision_floating-point_format) to
[PQ / SMPTE ST2084](https://en.wikipedia.org/wiki/Perceptual_quantizer) stored as
unsigned integer up to 16 bit.

The conversion is implemented in SSE2, AVX2, AVX512 and NEON SIMD, and using a fast approximate
formula to speed up conversion.


Benchmark result on Ryzen 7950x:

```
SSE2_prec12_size1280x720	Speed: 11.1056ms/it
SSE2_prec12_size1920x1080	Speed: 25.0222ms/it
SSE2_prec12_size3840x2160	Speed: 100.323ms/it
SSE2_prec16_size1280x720	Speed: 11.1167ms/it
SSE2_prec16_size1920x1080	Speed: 25.0862ms/it
SSE2_prec16_size3840x2160	Speed: 100.514ms/it
AVX2_prec12_size1280x720	Speed: 3.58655ms/it
AVX2_prec12_size1920x1080	Speed: 8.12563ms/it
AVX2_prec12_size3840x2160	Speed: 33.8692ms/it
AVX2_prec16_size1280x720	Speed: 3.58141ms/it
AVX2_prec16_size1920x1080	Speed: 8.06523ms/it
AVX2_prec16_size3840x2160	Speed: 34.285ms/it
AVX3_prec12_size1280x720	Speed: 2.50691ms/it
AVX3_prec12_size1920x1080	Speed: 5.64231ms/it
AVX3_prec12_size3840x2160	Speed: 23.0271ms/it
AVX3_prec16_size1280x720	Speed: 2.52674ms/it
AVX3_prec16_size1920x1080	Speed: 5.6914ms/it
AVX3_prec16_size3840x2160	Speed: 23.2591ms/it
```