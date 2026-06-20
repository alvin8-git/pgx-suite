#!/usr/bin/env python3
"""Assess whether an input BAM/CRAM was aligned ALT-aware (GRCh38).

Non-ALT-aware aligners (minimap2 default, bowtie2) or alignment to a no-ALT
reference mis-distribute reads across paralogs + ALT contigs, which silently
breaks paralog/SV genes — e.g. CYP2D6 called as a false *5/*5 deletion. This
module derives an ALT-awareness verdict from signals recorded in the file:

  1. @PG aligner identity (PN/ID/CL)        — bwa/bwa-mem2/DRAGEN vs minimap2/bowtie2
  2. @SQ ALT-contig count (SN ending _alt)  — 0 => ALT-aware was impossible
  3. fraction of reads on _alt contigs       — ~0 with ALT contigs present => suspect

Verdict feeds pgx-compare.py --alt-aware {yes,no,unknown}; 'no'/'unknown'
forces NO_CALL for paralog/SV genes.

Usage:
    pgx_altcheck.py <bam|cram>            # prints JSON verdict
    pgx_altcheck.py --selftest
"""
import json
import sys

# bwa-mem reads the GRCh38 .alt file at index time; DRAGEN/Isaac/Novoalign are
# ALT-aware. minimap2 (default) and bowtie2 are not.
_ALT_CAPABLE = ("bwa", "dragen", "isaac", "novoalign")
_NOT_ALT = ("minimap", "bowtie")
_MIN_ALT_FRAC = 0.0002  # 0.02%: bwa GRCh38+ALT leaves ~0.1-0.2% reads on ALT


def _aligner_from_header(header: str) -> tuple[str, str]:
    """Return (aligner_name_lower, command_line) of the first real aligner @PG."""
    for ln in header.splitlines():
        if not ln.startswith("@PG"):
            continue
        fields = {}
        for f in ln.split("\t")[1:]:
            k, _, v = f.partition(":")
            fields[k] = v
        name = (fields.get("PN") or fields.get("ID") or "").lower()
        if any(a in name for a in _ALT_CAPABLE + _NOT_ALT):
            return name, fields.get("CL", "")
    return "", ""


def _count_alt_contigs(header: str) -> int:
    n = 0
    for ln in header.splitlines():
        if not ln.startswith("@SQ"):
            continue
        for f in ln.split("\t"):
            if f.startswith("SN:") and f.endswith("_alt"):
                n += 1
    return n


def _alt_read_fraction(idxstats: str) -> float:
    tot = alt = 0
    for ln in idxstats.splitlines():
        p = ln.split("\t")
        if len(p) < 3:
            continue
        try:
            mapped = int(p[2])
        except ValueError:
            continue
        tot += mapped
        if p[0].endswith("_alt"):
            alt += mapped
    return (alt / tot) if tot else 0.0


def assess(header: str, idxstats: str) -> dict:
    aligner, cl = _aligner_from_header(header)
    n_alt = _count_alt_contigs(header)
    frac = _alt_read_fraction(idxstats)
    capable = any(a in aligner for a in _ALT_CAPABLE)
    not_alt = any(a in aligner for a in _NOT_ALT)

    if n_alt == 0:
        verdict, reason = "no", "reference has no ALT contigs (no-ALT/primary assembly) — ALT-aware alignment was impossible"
    elif not_alt:
        verdict, reason = "no", f"aligner '{aligner}' is not ALT-aware"
    elif capable and frac >= _MIN_ALT_FRAC:
        verdict, reason = "yes", f"{aligner}; {n_alt} ALT contigs; {frac*100:.3f}% reads on ALT"
    elif capable:
        verdict, reason = "unknown", (f"{aligner} but only {frac*100:.4f}% reads on ALT — "
                                      "the GRCh38 .alt file may not have been used at index time")
    else:
        verdict, reason = "unknown", f"aligner '{aligner or 'unrecognised'}' — cannot confirm ALT-awareness"

    return {"aligner": aligner or "unknown", "cmdline": cl[:200],
            "n_alt_contigs": n_alt, "alt_read_frac": round(frac, 5),
            "alt_aware": verdict, "reason": reason}


def _selftest() -> None:
    bwa_hdr = ("@HD\tVN:1.6\n"
               "@SQ\tSN:chr1\tLN:1\n@SQ\tSN:chr1_KI270706v1_alt\tLN:1\n"
               "@PG\tID:bwa\tPN:bwa\tCL:bwa mem hg38.fa -Y reads.fq\n")
    mm_hdr = ("@SQ\tSN:chr1\tLN:1\n@SQ\tSN:chr1_KI270706v1_alt\tLN:1\n"
              "@PG\tID:minimap2\tPN:minimap2\tCL:minimap2 -ax sr hg38.fa\n")
    noalt_hdr = "@SQ\tSN:chr1\tLN:1\n@PG\tID:bwa\tPN:bwa\tCL:bwa mem primary.fa\n"
    idx_alt = "chr1\t1\t1000\t0\nchr1_KI270706v1_alt\t1\t2\t0\n"   # ~0.2% on ALT
    idx_noalt = "chr1\t1\t1000\t0\nchr1_KI270706v1_alt\t1\t0\t0\n"  # 0% on ALT

    assert assess(bwa_hdr, idx_alt)["alt_aware"] == "yes"
    assert assess(mm_hdr, idx_alt)["alt_aware"] == "no"          # minimap2 -> no
    assert assess(noalt_hdr, idx_alt)["alt_aware"] == "no"       # no ALT contigs -> no
    assert assess(bwa_hdr, idx_noalt)["alt_aware"] == "unknown"  # bwa but ~0 ALT reads
    print("all pgx_altcheck self-checks passed")


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        _selftest()
        sys.exit(0)
    import subprocess
    bam = sys.argv[1]
    hdr = subprocess.run(["samtools", "view", "-H", bam], capture_output=True, text=True).stdout
    idx = subprocess.run(["samtools", "idxstats", bam], capture_output=True, text=True).stdout
    print(json.dumps(assess(hdr, idx), indent=2))
