#pragma once
// IEEE-754 binary16 ↔ binary32 decode. The SPOTLIGHT raw files store each
// (re, im) of a visibility as two 16-bit half-floats packed into one float
// slot — exactly what svfits's half_to_float() handles.
#include <cstdint>

namespace pico {

// Returns 32-bit float bit pattern equivalent of the given half. Ported from
// svfits/svsubs.c half_to_float (the canonical bit-twiddling).
inline float half_to_float(uint16_t h) noexcept {
    const uint32_t s  = (uint32_t(h) & 0x8000u) << 16;
    uint32_t      e  =  uint32_t(h) & 0x7C00u;
    uint32_t      m  =  uint32_t(h) & 0x03FFu;
    uint32_t      f;
    if (e == 0) {
        if (m == 0) { f = s; }
        else {
            // subnormal → normalize
            while ((m & 0x0400u) == 0) { m <<= 1; e -= 0x0400u; }
            e += 0x1C400u; m &= 0x03FFu;
            f = s | (e << 13) | (m << 13);
        }
    } else if (e == 0x7C00u) {
        // inf / NaN
        f = s | 0x7F800000u | (m << 13);
    } else {
        // normal
        f = s | ((e + 0x1C000u) << 13) | (m << 13);
    }
    float out;
    __builtin_memcpy(&out, &f, sizeof(out));
    return out;
}

} // namespace pico
