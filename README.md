# scRGB to PQ Converter

This library does one thing: convert pixels [scRGB](https://en.wikipedia.org/wiki/ScRGB) stored as
[Half float](https://en.wikipedia.org/wiki/Half-precision_floating-point_format) to
[PQ / SMPTE ST2084](https://en.wikipedia.org/wiki/Perceptual_quantizer) stored as
unsigned integer up to 16 bit. 

The conversion is implemented in SSE2, AVX2, AVX512 and NEON SIMD, and using a fast approximate
formula to speed up conversion.


Benchmark result on Ryzen 7950x:

```
SSE2_prec12_size1280x720	Speed: 11.9748ms/it
SSE2_prec12_size1920x1080	Speed: 26.9152ms/it
SSE2_prec12_size3840x2160	Speed: 108.035ms/it
SSE2_prec16_size1280x720	Speed: 12.0663ms/it
SSE2_prec16_size1920x1080	Speed: 27.1404ms/it
SSE2_prec16_size3840x2160	Speed: 108.456ms/it
AVX2_prec12_size1280x720	Speed: 4.08695ms/it
AVX2_prec12_size1920x1080	Speed: 9.22785ms/it
AVX2_prec12_size3840x2160	Speed: 38.4899ms/it
AVX2_prec16_size1280x720	Speed: 3.9446ms/it
AVX2_prec16_size1920x1080	Speed: 8.87443ms/it
AVX2_prec16_size3840x2160	Speed: 36.2998ms/it
AVX3_prec12_size1280x720	Speed: 2.44766ms/it
AVX3_prec12_size1920x1080	Speed: 5.49054ms/it
AVX3_prec12_size3840x2160	Speed: 22.2641ms/it
AVX3_prec16_size1280x720	Speed: 2.49396ms/it
AVX3_prec16_size1920x1080	Speed: 5.61471ms/it
AVX3_prec16_size3840x2160	Speed: 22.7128ms/it
```