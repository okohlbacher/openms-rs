#!/usr/bin/env python3
"""Synthesise the F7b stream oracle inputs (deterministic: mtime 0, fixed levels)."""
import bz2
import gzip
import os
import struct
import zlib

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "data"))


def gz(data):
    return gzip.compress(data, compresslevel=9, mtime=0)


def raw_deflate(data):
    co = zlib.compressobj(9, zlib.DEFLATED, -15)
    return co.compress(data) + co.flush()


def trailer(data):
    return struct.pack("<II", zlib.crc32(data), len(data) & 0xFFFFFFFF)


payload = b"".join(struct.pack("<I", i) for i in range(20000))
files = {}
files["plain.txt"] = b"plain text, not gzip\n"
files["one.txt"] = b"x"
files["empty.bin"] = b""
files["members.gz"] = gz(b"first member\n") + gz(b"second member\n") + b"\0\0\0 trailing garbage"
full = gz(payload)
files["payload.gz"] = full
files["trunc_trailer.gz"] = full[:-5]
files["trunc_body.gz"] = full[: len(full) // 2]
b = bytearray(full)
b[-8] ^= 1
files["bad_crc.gz"] = bytes(b)
b = bytearray(full)
b[-1] ^= 1
files["bad_len.gz"] = bytes(b)
hdr = b"\x1f\x8b\x08\x02" + b"\0\0\0\0" + b"\x00\x03"
good = zlib.crc32(hdr) & 0xFFFF
files["hcrc_good.gz"] = hdr + struct.pack("<H", good) + raw_deflate(b"header crc ok\n") + trailer(b"header crc ok\n")
files["hcrc_bad.gz"] = hdr + struct.pack("<H", good ^ 1) + raw_deflate(b"header crc bad\n") + trailer(b"header crc bad\n")
ehdr = (b"\x1f\x8b\x08" + bytes([0x04 | 0x08 | 0x10]) + b"\0\0\0\0" + b"\x00\x03"
        + struct.pack("<H", 3) + b"abc" + b"name.txt\0" + b"a comment\0")
files["extra.gz"] = ehdr + raw_deflate(b"extra fields\n") + trailer(b"extra fields\n")
files["bad_second.gz"] = gz(b"good\n") + b"\x1f\x8b\x09\x00" + b"\0" * 20
files["reserved_flag.gz"] = b"\x1f\x8b\x08\x20" + gz(b"x")[4:]
bpayload = b"".join(struct.pack("<I", i * 7) for i in range(30000))
bb = bz2.compress(bpayload, 9)
files["payload.bz2"] = bb
files["two.bz2"] = bb + bz2.compress(b"second stream", 9)
files["trunc.bz2"] = bb[:-3]
files["plain.bz2"] = b"plain text, not bzip2\n"
import io
import zipfile

buf = io.BytesIO()
with zipfile.ZipFile(buf, "w", zipfile.ZIP_STORED) as z:
    info = zipfile.ZipInfo("../escape.txt", date_time=(2020, 1, 1, 0, 0, 0))
    z.writestr(info, b"escaped")
files["traversal.zip"] = buf.getvalue()
for name, data in files.items():
    with open(name, "wb") as f:
        f.write(data)
print(len(payload), len(full), len(bpayload), len(bb))
