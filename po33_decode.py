"""Decode PO-33 transfer WAVs; export checksummed data and audition WAVs.

Requires NumPy. Run with no arguments for self-checks. See README.md and
docs/FORMAT.md. Raw audio and pattern bytes are verified against the packet
CRCs; the signed mu-law playback mapping and inclusive trim ends are inferred.
"""
import argparse
from binascii import crc_hqx
from collections import Counter
import hashlib
import json
from pathlib import Path
import struct
import wave

import numpy as np

PREAMBLE = (bytes(range(256)) * 16)[:4092] + bytes.fromhex('4db8165d')
NATIVE_RATE = 23437.5
TABLE = [crc_hqx(bytes([i]), 0) for i in range(256)]


def crc(state, value):
    return ((state << 8) & 65535) ^ TABLE[(state >> 8) ^ value]


def unpack_phases(left, right):
    left, right = np.diff(left.astype(int)) % 4, np.diff(right.astype(int)) % 4
    n = min(len(left), len(right)) // 2 * 2
    return (right[:n:2] | left[:n:2] << 2 | right[1:n:2] << 4
            | left[1:n:2] << 6).astype(np.uint8).tobytes()


def interpolate(signal, times):
    indices = np.floor(times).astype(int)[:, None] + np.arange(-8, 9)
    if indices.min() < 0 or indices.max() >= len(signal):
        raise ValueError('Transfer is too short for timing recovery')
    d = times[:, None] - indices
    kernel = np.sinc(d) * np.where(abs(d) < 9, .5 + .5 * np.cos(np.pi*d/9), 0)
    return np.sum(signal[indices] * kernel, axis=1) / np.sum(kernel, axis=1)


def quadrants(signal, times, period):
    x = interpolate(signal, times + period*.125)
    y = interpolate(signal, times + period*.375)
    q = np.array([2, 3, 1, 0])[(x > 0).astype(int)*2 + (y > 0)]
    cost = float(np.mean((abs(x)-abs(y))**2) / max(1e-9, np.mean(x*x+y*y)))
    return q.astype(np.uint8), cost


def read_audio(path):
    with wave.open(str(path)) as wav:
        fs, width, frames = wav.getframerate(), wav.getsampwidth(), wav.getnframes()
        if wav.getnchannels() != 2 or width not in (1, 2, 3):
            raise ValueError('Expected stereo 8-, 16- or 24-bit PCM WAV')
        if not 32000 <= fs <= 192000 or frames*2*width > 160*1024*1024:
            raise ValueError('Supported input: 32–192 kHz, at most 160 MiB of PCM')
        raw = wav.readframes(frames)
        if len(raw) != frames*2*width:
            raise ValueError('Truncated WAV audio')
    if width == 1:
        audio = (np.frombuffer(raw, np.uint8).astype(float)-128)*256
    elif width == 2:
        audio = np.frombuffer(raw, '<i2').astype(float)
    else:
        b = np.frombuffer(raw, np.uint8).reshape(-1, 3).astype(np.int32)
        v = b[:, 0] | b[:, 1] << 8 | b[:, 2] << 16
        audio = ((v ^ 0x800000)-0x800000)/256
    audio = np.pad(audio.reshape(-1, 2), ((1000, 1000), (0, 0)))
    # Centered moving average removes capture DC without SciPy.
    for ch in (0, 1):
        total = np.r_[0, np.cumsum(np.pad(audio[:, ch], (250, 250)))]
        audio[:, ch] -= (total[501:]-total[:-501])/501
    return audio, fs, width


def demodulate(path, symbol_rate=None):
    audio, fs, width = read_audio(path)
    starts = np.flatnonzero(abs(audio[:fs*6, 0]) > 1000)
    if not len(starts):
        raise ValueError('No sufficiently strong signal in the first six seconds')
    start = int(starts[0])
    if symbol_rate is not None and not 7500 <= symbol_rate <= 8100:
        raise ValueError('Symbol rate must be between 7500 and 8100 Hz')
    period = fs/(symbol_rate or 7812.5)
    if symbol_rate is None:
        lag = round(4096*period)
        span = int(np.ceil(3*fs/44100))
        lags = np.arange(lag-span, lag+span+1)
        reference = audio[start+2000:start+14000, 0]
        scores = []
        for lag in lags:
            shifted = audio[start+2000+lag:start+14000+lag, 0]
            if len(shifted) != len(reference) or not len(reference):
                raise ValueError('Transfer is too short for clock calibration')
            scores.append(np.dot(reference, shifted) / max(1e-9, np.linalg.norm(shifted)))
        i = int(np.argmax(scores))
        denominator = 2*(scores[i-1]-2*scores[i]+scores[i+1]) if 0 < i < len(lags)-1 else 0
        fraction = (scores[i-1]-scores[i+1])/denominator if denominator else 0
        period = (lags[i]+fraction)/4096
    initial_rate = fs/period
    times = start + np.arange(8300)*period
    best = None
    # The exact ramp/sync disambiguates the several possible phase minima.
    for offset in np.linspace(-period, period, 81):
        candidates = [quadrants(audio[:, ch], times+offset, period) for ch in (0, 1)]
        if unpack_phases(candidates[0][0], candidates[1][0])[:4096] != PREAMBLE:
            continue
        cost = sum(c[1] for c in candidates)
        if best is None or cost < best[0]:
            best = cost, offset
    if best is None:
        raise ValueError('Preamble lock failed; check stereo channels and --symbol-rate')
    position = start+best[1]
    phases = [[], []]
    # Track both phase and clock rate; phase-only correction fails on Escalate.
    block = 512
    while position+block*period < len(audio)-20:
        times = position+np.arange(block)*period
        delta = period/12
        costs = [quadrants(audio[:, 0], times[64::3]+off, period)[1]
                 for off in (-delta, 0, delta)]
        denominator = costs[0]-2*costs[1]+costs[2]
        correction = (np.clip(delta*(costs[0]-costs[2])/(2*denominator), -delta, delta)
                      if denominator > 0 else (-delta, 0, delta)[int(np.argmin(costs))])
        for ch in (0, 1):
            phases[ch].append(quadrants(audio[:, ch], times+correction, period)[0])
        position += block*period+correction
        period += correction/block*.3
        if not fs/8500 < period < fs/7000:
            raise ValueError('Lost symbol-clock lock')
    if not phases[0]:
        raise ValueError('No symbols recovered')
    data = unpack_phases(*(np.concatenate(ch) for ch in phases))
    return data, {'sample_rate': fs, 'pcm_bits': width*8, 'initial_symbol_rate': initial_rate}


def parse_packets(data):
    if data[:4096] != PREAMBLE:
        raise ValueError('Invalid PO-33 preamble/sync')
    state, position = 0x1d0f, 4096
    raw = bytearray()

    def take(size):
        nonlocal state, position
        if position+size > len(data):
            raise ValueError('Truncated transfer; end packet not reached')
        out = bytearray()
        for coded in data[position:position+size]:
            value = (coded-(state & 255)) & 255
            state = crc(state, value)
            out.append(value)
        position += size
        raw.extend(out)
        return bytes(out)

    packets = []
    while True:
        offset = position
        header = take(3)
        tag, size = struct.unpack('<HB', header)
        payload = take(size)
        expected = state.to_bytes(2, 'little')
        if take(2) != expected:
            raise ValueError(f'Checksum failed at wire byte {offset}; no export written')
        if not packets and (tag != 0x985b or size != 0):
            raise ValueError('Unrecognized start packet')
        packets.append((tag, payload))
        if tag == 0x826f:
            if size:
                raise ValueError('Invalid end packet')
            return packets, bytes(raw), position


def collect_objects(packets):
    patterns, objects, configs = {}, {}, []
    pattern = sample = None
    for tag, payload in packets:
        if tag in (0xa041, 0x68b5):
            limit = 16 if tag == 0xa041 else 136
            if len(payload) != 1 or payload[0] >= limit:
                raise ValueError('Invalid pattern/sample selector')
            target = patterns if tag == 0xa041 else objects
            if payload[0] in target:
                raise ValueError('Duplicate selector')
            target[payload[0]] = bytearray()
            if tag == 0xa041:
                pattern = payload[0]
            else:
                sample = payload[0]
        elif tag == 0xc498:
            if pattern is None:
                raise ValueError('Pattern data without selector')
            patterns[pattern].extend(payload)
        elif tag == 0xef78:
            if sample is None:
                raise ValueError('Sample data without selector')
            objects[sample].extend(payload)
        elif tag == 0x7d5c:
            configs.append(payload)
        elif tag not in (0x985b, 0x826f):
            raise ValueError(f'Unsupported packet tag {tag:04x}')
    if len(patterns) != 16 or any(len(b) != 1032 for b in patterns.values()):
        raise ValueError('Incomplete or unsupported pattern layout')
    if len(objects) != 136 or len(configs) != 2 or any(len(b) != 255 for b in configs):
        raise ValueError('Incomplete or unsupported sound/configuration layout')
    roots, entries = {}, []
    for index, data in objects.items():
        if not data:
            continue
        if len(data) < 9:
            raise ValueError('Truncated sample metadata')
        start, end, source = struct.unpack('<IIB', data[:9])
        if len(data) > 9:
            if source != index:
                raise ValueError('Unexpected audio owner')
            roots[index] = bytes(data[9:])
        entries.append({'object': index, 'source': source, 'start': start, 'end': end,
                        'slot': index+1 if index < 8 else 9+(index-8)//16,
                        'pad': None if index < 8 else 1+(index-8) % 16})
    for entry in entries:
        if entry['source'] not in roots or not 0 <= entry['start'] <= entry['end']:
            raise ValueError('Invalid sample source or trim order')
        if entry['end'] >= len(roots[entry['source']]):
            # Escalate's last bank has ten such references despite valid packet CRCs.
            entry['unresolved'] = 'Trim exceeds recovered source audio; no WAV exported'
        elif entry['end'] == entry['start']:
            entry['unresolved'] = 'Equal trim endpoints; no WAV exported'
    return patterns, roots, entries, configs


def expand_audio(data, codec):
    values = np.frombuffer(data, np.int8).astype(np.int32)
    if codec == 'signed8':
        return values/128
    # Inferred monotonic signed ordering of G.711 mu-law codes; both zero codes map to silence.
    magnitude = np.where(values >= 0, values, -values-1)
    expanded = ((132+((magnitude & 15) << 3)) << (magnitude >> 4))-132
    return np.where(values >= 0, expanded, -expanded)/32768


def audition_wav(path, data, codec):
    audio = expand_audio(data, codec)
    # WAV rates are integers; resample the fractional native clock to 44.1 kHz.
    times = np.arange(round(len(audio)*44100/NATIVE_RATE))*NATIVE_RATE/44100
    pcm = np.round(np.interp(times, np.arange(len(audio)), audio)*32767).astype('<i2')
    with wave.open(str(path), 'wb') as wav:
        wav.setparams((1, 2, 44100, len(pcm), 'NONE', 'not compressed'))
        wav.writeframes(pcm.tobytes())


def export(path, output, codec, symbol_rate):
    if output.exists():
        raise ValueError('Output folder already exists; choose a new folder')
    data, signal_info = demodulate(path, symbol_rate)
    packets, raw, consumed = parse_packets(data)
    patterns, roots, entries, configs = collect_objects(packets)
    # Validate the complete transfer first; preserve unsupported trims without guessing.
    output.mkdir(parents=True)
    (output/'transfer.bin').write_bytes(data[:consumed])
    (output/'packets.bin').write_bytes(raw)
    for i, payload in enumerate(configs):
        (output/f'configuration-{i}.bin').write_bytes(payload)
    for i, payload in patterns.items():
        (output/f'pattern-{i+1:02}.bin').write_bytes(payload)
    for i, payload in roots.items():
        (output/f'audio-{i:03}.pcm8').write_bytes(payload)
    for entry in entries:
        if 'unresolved' in entry:
            continue
        label = f"slot-{entry['slot']:02}" + (f"-pad-{entry['pad']:02}" if entry['pad'] else '')
        entry['wav'] = label+'.wav'
        # Inclusive end is inferred from factory lengths (e.g. 0..4095).
        chunk = roots[entry['source']][entry['start']:entry['end']+1]
        audition_wav(output/entry['wav'], chunk, codec)
    with path.open('rb') as source:
        source_hash = hashlib.file_digest(source, 'sha256').hexdigest()
    manifest = {
        'format': 'po33-decoded', 'source_file': path.name,
        'source_sha256': source_hash,
        **signal_info, 'native_sample_rate': NATIVE_RATE, 'codec': codec,
        'audio_conversion_verified_against_hardware': False,
        'packet_count': len(packets), 'all_checksums_valid': True, 'end_marker_valid': True,
        'tags': {f'{tag:04x}': count for tag, count in Counter(t for t, _ in packets).items()},
        'pattern_count': len(patterns), 'root_audio_count': len(roots),
        'audition_wav_count': sum('wav' in e for e in entries),
        'unresolved_trims': sum('unresolved' in e for e in entries),
        'bpm': configs[1][129], 'samples': entries,
        'limits': 'Raw audio/pattern bytes are verified; mu-law ordering and inclusive trim-end interpretation are inferred. Pattern event semantics are not yet decoded.',
    }
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return {k: manifest[k] for k in ('packet_count', 'all_checksums_valid', 'end_marker_valid', 'root_audio_count', 'audition_wav_count', 'unresolved_trims', 'bpm')}


def self_check():
    ramp = bytes(range(256))*16
    b = np.frombuffer(ramp, np.uint8)
    left = np.column_stack(((b >> 2) & 3, b >> 6)).ravel()
    right = np.column_stack((b & 3, (b >> 4) & 3)).ravel()
    for rotation in range(4):
        phases = [np.r_[rotation, (rotation+np.cumsum(d)) % 4] for d in (left, right)]
        assert unpack_phases(*phases) == ramp
    state, coded = 0x1d0f, bytearray(PREAMBLE)
    for tag in (0x985b, 0x826f):
        header = struct.pack('<HB', tag, 0)
        trailer = crc_hqx(header, state).to_bytes(2, 'little')
        for raw in header+trailer:
            coded.append((raw+(state & 255)) & 255)
            state = crc(state, raw)
    assert coded[4096:4101].hex() == '6a9a34e248'  # Captured official reset packet.
    assert parse_packets(coded)[0] == [(0x985b, b''), (0x826f, b'')]
    for bad in (coded[:-1], coded[:4100], coded[:4103]+bytes([coded[4103] ^ 1])+coded[4104:]):
        try:
            parse_packets(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Corrupt/truncated transfer was accepted')
    assert expand_audio(bytes([0, 255, 127, 128]), 'mulaw').tolist() == [0, 0, 32124/32768, -32124/32768]
    fixture = [(0x7d5c, bytes(255))]*2
    for i in range(16):
        fixture.extend([(0xa041, bytes([i])), (0xc498, bytes(1032))])
    for i in range(136):
        fixture.append((0x68b5, bytes([i])))
        if i < 3:
            # Audio owner, valid alias, and an out-of-range alias retained without a WAV.
            blob = struct.pack('<IIB', 0, 3 if i < 2 else 8, 0)
            fixture.append((0xef78, blob+(bytes(4) if i == 0 else b'')))
    patterns, roots, entries, _ = collect_objects(fixture)
    assert len(patterns) == 16 and roots == {0: bytes(4)}
    assert 'unresolved' not in entries[1] and 'unresolved' in entries[2]
    try:
        collect_objects(fixture+[(0x68b5, b'\x00')])
    except ValueError:
        pass
    else:
        raise AssertionError('Duplicate sample selector was accepted')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wav', nargs='?', type=Path)
    parser.add_argument('output', nargs='?', type=Path)
    parser.add_argument('--symbol-rate', type=float)
    parser.add_argument('--codec', choices=['mulaw', 'signed8'], default='mulaw')
    args = parser.parse_args(argv)
    self_check()
    if args.wav:
        if not args.output:
            parser.error('Choose a new output folder')
        try:
            print(json.dumps(export(args.wav, args.output, args.codec, args.symbol_rate), indent=2))
        except (ValueError, OSError, wave.Error) as error:
            parser.exit(1, str(error)+'\n')
    else:
        print('PASS: stereo packing, captured reset packet, CRC rejection, truncation, sample references, audio conversion')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
