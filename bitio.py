from __future__ import print_function

import re
from pathlib import Path


def _byte_to_bits(byte, bigendian=True):
    if bigendian:
        return [1 if (byte & (0x80 >> i)) else 0 for i in range(8)]
    return [(byte >> i) & 1 for i in range(8)]


def _looks_like_text(data):
    if not data:
        return True
    sample = data[:4096]
    bad = 0
    for b in sample:
        if b in (9, 10, 13):
            continue
        if 32 <= b <= 126:
            continue
        bad += 1
    return (bad / float(len(sample))) < 0.05


def _parse_group_to_bits(text, input_format="auto", bigendian=True):
    s = text.strip()
    if not s:
        return []

    if input_format in ("auto", "bit"):
        bit_str = "".join(ch for ch in s if ch in ("0", "1"))
        if bit_str:
            return [1 if ch == "1" else 0 for ch in bit_str]
        if input_format == "bit":
            return []

    if input_format in ("auto", "hex"):
        hex_str = re.sub(r"[^0-9a-fA-F]", "", s)
        if hex_str and (len(hex_str) % 2 == 0):
            try:
                data = bytes.fromhex(hex_str)
            except ValueError:
                data = None
            if data is not None:
                bits = []
                for b in data:
                    bits.extend(_byte_to_bits(b, bigendian=bigendian))
                return bits
        if input_format == "hex":
            return []

    return []


def _read_groups_from_text_file(file_path):
    content = Path(file_path).read_text(encoding="utf-8", errors="ignore")
    normalized = content.replace("\r\n", "\n").replace("\r", "\n")
    groups = re.split(r"\n\s*\n+", normalized)
    return [g for g in groups if g.strip()]


def _read_bits_from_binary_file(file_path, bigendian=True):
    bits = []
    with open(file_path, "rb") as f:
        while True:
            data = f.read(16384)
            if not data:
                break
            for b in data:
                bits.extend(_byte_to_bits(b, bigendian=bigendian))
    return bits


def iter_bit_sequences(key_path, input_format="auto", bigendian=True, recursive=False):
    p = Path(key_path)
    if p.is_file():
        data = p.read_bytes()
        if _looks_like_text(data):
            groups = _read_groups_from_text_file(p)
            for idx, group in enumerate(groups, start=1):
                bits = _parse_group_to_bits(group, input_format=input_format, bigendian=bigendian)
                yield {"source": str(p), "index": idx, "bits": bits}
        else:
            bits = _read_bits_from_binary_file(p, bigendian=bigendian)
            yield {"source": str(p), "index": 1, "bits": bits}
        return

    if p.is_dir():
        it = p.rglob("*") if recursive else p.glob("*")
        for fp in it:
            if not fp.is_file():
                continue
            for item in iter_bit_sequences(str(fp), input_format=input_format, bigendian=bigendian, recursive=False):
                yield item
        return

    raise FileNotFoundError("找不到这个路径：" + str(p))

