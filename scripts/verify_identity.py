#!/usr/bin/env python
"""Recompute the two code fingerprints recorded in every result record.

gate tree  (repository root as shipped)                      -> bedca2849bb2
main tree  (root with reference/main_tree/* overlaid)        -> 08e28b83407c

The fingerprint is the first 12 hex digits of a SHA-1 over the names and SHA-1 digests of the
files listed in turb/runner.py:code_fingerprint().  Usage:  python scripts/verify_identity.py
"""
import hashlib, os, re, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPECT = {"gate": "bedca2849bb2", "main": "08e28b83407c"}


def file_list():
    src = open(os.path.join(ROOT, "turb", "runner.py"), encoding="utf-8").read()
    body = src[src.index("def code_fingerprint"):]
    block = body[body.index("files = ["):body.index("]", body.index("files = [")) + 1]
    return re.findall(r'"([^"]+)"', block)


def fp(overlay=None):
    h = hashlib.sha1()
    for fn in file_list():
        p = os.path.join(ROOT, fn)
        if overlay and os.path.exists(os.path.join(overlay, fn)):
            p = os.path.join(overlay, fn)
        try:
            with open(p, "rb") as fh:
                h.update(fn.encode())
                h.update(hashlib.sha1(fh.read()).digest())
        except OSError:
            h.update(b"MISSING:" + fn.encode())
    return h.hexdigest()[:12]


if __name__ == "__main__":
    got = {"gate": fp(), "main": fp(os.path.join(ROOT, "reference", "main_tree"))}
    ok = True
    for k in ("gate", "main"):
        s = "OK" if got[k] == EXPECT[k] else "MISMATCH"
        ok &= got[k] == EXPECT[k]
        print(f"{k} tree: {got[k]} (expected {EXPECT[k]}) {s}")
    sys.exit(0 if ok else 1)
