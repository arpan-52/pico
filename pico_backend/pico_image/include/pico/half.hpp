#pragma once
// IEEE-754 binary16 ↔ binary32 decode. The SPOTLIGHT raw files store each
// (re, im) of a visibility as two 16-bit half-floats packed into one float
// slot — exactly what svfits's half_to_float() handles.
#include <cstdint>

namespace pico {

inline uint32_t as_uint(float x) noexcept {
    uint32_t u; __builtin_memcpy(&u, &x, sizeof(u)); return u;
}
inline float as_float(uint32_t u) noexcept {
    float x; __builtin_memcpy(&x, &u, sizeof(x)); return x;
}

// EXACT port of svfits/utils.c half_to_float: a 1-5-10 half format *WITHOUT
// infinity* — the top exponent (0x7C00, e=31) is a normal large value, not
// Inf/NaN (range ±131008). The GMRT correlator's float_to_half saturates large
// visibilities into this range, so we must decode them as finite (svfits does,
// and lets the MAD clip remove them). The standard IEEE decoder would turn
// these into NaN/Inf and the isfinite guards would silently drop the very
// (large-amplitude RFI) samples svfits keeps — desyncing every downstream stat.
inline float half_to_float(uint16_t h) noexcept {
    const uint32_t e = (uint32_t(h) & 0x7C00u) >> 10;   // exponent
    const uint32_t m = (uint32_t(h) & 0x03FFu) << 13;   // mantissa
    const uint32_t v = as_uint(float(m)) >> 23;         // leading-zero count hack
    return as_float(
        (uint32_t(h) & 0x8000u) << 16                                   // sign
        | (e != 0) * ((e + 112) << 23 | m)                              // normalized
        | ((e == 0) & (m != 0)) * ((v - 37) << 23                       // denormalized
                                   | ((m << (150 - v)) & 0x007FE000u)));
}

} // namespace pico
