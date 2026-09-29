# po-backup-decode

Single-file Python tool that demodulates a Teenage Engineering PO-33 sound and pattern
transfer WAV, then exports the decoded packet data, raw sample bytes, patterns, BPM and
audition WAVs.

The module is `po33_decode.py` and the installed command is `po33-decode`. Only the PO-33
transfer format is decoded here; see [Known limits](#known-limits).

This repository contains no device material. You supply the transfer recording, either
captured from your own unit or downloaded from the manufacturer's PO-33 page. The tool
reads that file and writes plain output; it never contacts a network.

The protocol description is in [docs/FORMAT.md](docs/FORMAT.md).

## Requirements

Python 3.9 or newer, NumPy, and a stereo 8-, 16- or 24-bit PCM WAV of the transfer. No
other runtime dependency.

```sh
pip install numpy
```

## Usage

```sh
python3 po33_decode.py                       # self-checks only, no input file needed
python3 po33_decode.py transfer.wav out/     # decode
python3 po33_decode.py transfer.wav out/ --symbol-rate 7812.485
python3 po33_decode.py transfer.wav out/ --codec signed8
```

The output folder must not exist; the tool refuses to overwrite one. It validates the
whole transfer before writing anything, so a bad packet CRC, a missing end marker or a
truncation leaves no partial output.

Install as a command with `pipx install .` after cloning, or run the file directly.

### Options

| Option | Effect |
|---|---|
| `--symbol-rate HZ` | Fix the symbol clock instead of auto-calibrating it. Accepted range 7500 to 8100 Hz, nominal 7812.5. |
| `--codec mulaw` | Default. Treat each audio byte as a signed G.711 mu-law code. |
| `--codec signed8` | Treat each audio byte as linear signed 8-bit, for comparison. |

Auto-calibration is the default because the two official transfers observed here run at
approximately 7,812.485 Hz and approximately 7,812.185 Hz. Both are close to nominal but
differ enough that a fixed clock fails on the longer capture.

## Output files

| File | Contents |
|---|---|
| `manifest.json` | Source file name and SHA-256, detected input format, symbol rate, packet counts, BPM, per-object trim and source references, unresolved trims, codec, counts. |
| `transfer.bin` | Packed, whitened wire bytes through the validated end packet. |
| `packets.bin` | Dewhitened packet headers, payloads and CRCs, concatenated. |
| `configuration-0.bin`, `configuration-1.bin` | The two 255-byte configuration payloads. Byte 129 of the second payload is the BPM value. |
| `pattern-01.bin` to `pattern-16.bin` | The 16 pattern payloads, 1,032 bytes each. |
| `audio-NNN.pcm8` | Raw owner audio bytes before waveform expansion, including 4,096-byte padding. |
| `slot-NN.wav`, `slot-NN-pad-NN.wav` | Individually trimmed audition sounds, resampled to 44.1 kHz mono 16-bit. |

Melodic slots 1 to 8 come from objects 0 to 7. Sliced slots 9 to 16 come from objects 8
to 135, with 16 pads per slot. Object and WAV counts are not counts of independent
original recordings, because several objects refer to different slices of one recording.

## Verified and inferred

The split matters when you use the output.

| Result | Status |
|---|---|
| Transport demodulation | Verified: both official transfers and three independent controlled captures decode with every packet CRC and the end marker accepted. |
| Raw audio bytes | Verified: exact recovered bytes, written to `audio-*.pcm8` before any conversion. |
| Trim endpoints and source references | Recovered bytes. Inclusive-end interpretation is inferred from factory lengths such as 0 to 4,095. |
| Pattern payloads | Verified bytes, 1,032 per pattern. The note and effect layout is not mapped, so patterns do not reconstruct the original songs. |
| BPM | Verified by a controlled round trip: only byte 129 of the second configuration payload changes between BPM 60 and 61. |
| Audition WAV conversion | Inferred. Assumes a native rate of 23,437.5 Hz and a monotonic signed ordering of mu-law codes. Not compared with hardware. |
| Hardware sound | Not established. Nothing here has been verified against a physical unit's output. |

## Known limits

- Compatibility with PO-133, PO-35 and PO-137 is not established. Their guides document a
  backup and restore facility, but no decoder for their formats was verified here. Do not
  assume this format applies to them.
- Ten slice references in one observed transfer point past the end of the recovered source
  audio even though every CRC passes. Those entries are flagged in `manifest.json` and get
  no WAV. No audio is fabricated or clamped for them.
- The tool decodes sound and pattern data. It does not write transfers back to hardware and
  does not reconstruct songs.

## Verify

```sh
python3 po33_decode.py
```

The self-check needs no input and covers stereo packing under all four phase rotations,
the captured official reset packet, CRC-rejection and truncation-rejection paths, sample
alias handling, an out-of-range alias retained without a WAV, duplicate-selector rejection
and the audio conversion endpoints. Continuous integration runs the same check on Python
3.9 and 3.12.

## Legal

This tool ships no manufacturer material. Device factory audio and any decoded output are
Teenage Engineering's material, and third-party samples may be present in a transfer you
did not create. Teenage Engineering's published terms permit use of that material in the
creation, performance, recording and distribution of your own original creation, and
prohibit redistributing, reselling, leasing, sharing, reproducing or publicly displaying
it in its original or modified form unless it is substantially transformed into an
original creation. See
[teenage.engineering/ip-rights-te-products](https://teenage.engineering/ip-rights-te-products).

Publishing decoded samples is your decision and your risk; this repository only provides
the decoder. This is not legal advice.

Not affiliated with, authorised by or endorsed by Teenage Engineering AB. Teenage
Engineering, Pocket Operator, PO-33 and PO-133 are trademarks of Teenage Engineering AB,
used here only to identify the hardware this tool reads.

MIT licensed; see [LICENSE](LICENSE).

## Prior public work

The format was worked out independently here, and the following public investigations and
implementations were consulted. None of their code is included.

- [cspotcode/PO-33-KO-backup-re](https://github.com/cspotcode/PO-33-KO-backup-re): original
  research and captured examples. Its linked manufacturer specification is the source of
  the mu-law and native sample-rate lead.
- [charlesvestal/po33-backup-re](https://github.com/charlesvestal/po33-backup-re) and
  [sZXZ/po33-experements](https://github.com/sZXZ/po33-experements): earlier partial
  payload interpretations.
- [ericlewis/libpo32](https://github.com/ericlewis/libpo32): PO-32 protocol, including the
  whitening transform and CRC convention that PO-33 shares.
