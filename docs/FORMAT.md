# PO-33 transfer format

Observed by demodulating and decoding two official transfer WAVs and three independent
controlled captures. Every statement below is either a checked protocol constant or is
marked inferred. No manufacturer specification for this format was located; the values
come from the captures themselves.

## Signal

Stereo differential quadrature phase modulation with four phase states, carrying two bits
per channel per symbol.

| Parameter | Value |
|---|---|
| Nominal symbol rate | 7,812.5 Hz |
| Bits per second across both channels | 31,250 |
| Phase differences | Taken modulo four |
| Observed auto-calibrated rates | 7,812.485 Hz (one official transfer), 7,812.185 Hz (another) |

Preamble: 4,092 bytes of the repeating `00` to `FF` ramp, followed by the four sync bytes
`4D B8 16 5D`. Both matched exactly in every capture.

A phase-only correction does not hold lock across a full-length capture. The receiver
tracks phase and clock rate together, adjusting every 512 symbols, and re-checks the exact
ramp to resolve the several phase minima that fit the same signal.

## Byte packing

Two stereo symbols form one byte, right channel first:

```text
byte[i] = R[2i] | L[2i]<<2 | R[2i+1]<<4 | L[2i+1]<<6
```

## Whitening and checksums

Byte 4,096 onward is whitened with the same additive transform and initial state `1D0F`
that PO-32 uses. The state carries across packet boundaries and is not reset per packet.

```python
from binascii import crc_hqx
state = 0x1d0f
for coded in body:
    raw = (coded - (state & 255)) & 255
    state = crc_hqx(bytes([raw]), state)
```

Subtraction, not XOR, then the CRC-16 update. Removing the preamble and subtracting bytes
is not sufficient, because the trailer bytes are whitened in the same way.

## Packets

After byte 4,096 the stream is a sequence of packets:

| Field | Size |
|---|---|
| Tag | 2 bytes, little-endian |
| Length | 1 byte, unsigned |
| Payload | Length bytes |
| CRC | 2 bytes, little-endian |

The CRC is the running state immediately after the header and payload. The trailer's own
two bytes are whitened and advance the state for the next packet.

| Tag | Observed role |
|---|---|
| `985B` | Start or reset, zero payload |
| `7D5C` | Configuration, 255 bytes |
| `A041` | Select pattern 0 to 15, one-byte payload |
| `C498` | Pattern-data chunk |
| `68B5` | Select sample object 0 to 135, one-byte payload |
| `EF78` | Sample-object chunk |
| `826F` | End, zero payload |

A complete transfer contains one start packet, two configuration packets of 255 bytes each,
16 pattern selectors with 1,032 bytes of pattern data each, 136 sample selectors, and the
end packet. Decoding stops at the end packet, before recording-tail noise.

The first configuration payload begins with a high-entropy-looking arithmetic sequence.
Entropy alone did not establish encryption, and no encryption was found.

## Audio objects

Concatenate `EF78` payloads under each sample selector.

A populated object starts with a nine-byte `<IIB` header: start frame, end frame, source
object index. An owner object carries audio bytes after the header and uses its own index
as the source. An alias carries only the header and points at another object's audio. Empty
objects have no payload. Owner audio buffers are padded to 4,096-byte boundaries.

Objects 0 to 7 map to melodic slots 1 to 8. Objects 8 to 135 map to slots 9 to 16 with 16
slices per slot.

### Trim interpretation, inferred

Trim ends are treated as inclusive, inferred from factory ranges such as 0 to 4,095. Equal
endpoints are preserved in the manifest with no WAV. In one official transfer, ten slices
reference a range longer than the recovered source audio while every packet CRC still
passes; the cause is unresolved and those slices get no WAV.

### Audio conversion, inferred

Audition WAVs assume a native sampling rate of 23,437.5 Hz and a monotonic signed ordering
of mu-law codes. For signed byte `q`, let `m = q` if nonnegative, otherwise `-q-1`, expand
the magnitude as `((132 + ((m & 15) << 3)) << (m >> 4)) - 132`, then apply the sign. `00`
and `FF` both map to silence. Output is resampled to 44.1 kHz mono 16-bit.

Neither the mu-law ordering nor the sample rate has been verified against hardware. The
`--codec signed8` option supplies a linear interpretation for comparison.

## Unmapped

Pattern bytes are recovered exactly but their note and effect layout is not mapped. These
exports do not restore the original songs.
