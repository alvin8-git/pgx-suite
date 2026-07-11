#!/usr/bin/env python3
"""Minimal Markdown->LaTeX converter tuned to this manuscript.

Not a general converter — handles the constructs used in manuscript.md:
headings, **bold**, `code`, pipe tables (->longtable), '- ' lists,
numbered reference list, blockquote notes, and the unicode/specials in the text.

Run:  python3 paper/md2tex.py   # writes manuscript.tex
Then: pdflatex manuscript.tex   (needs longtable, booktabs, geometry)
"""
import re, os

SRC = os.path.join(os.path.dirname(__file__), "manuscript.md")
OUT = os.path.join(os.path.dirname(__file__), "manuscript.tex")

UNI = {"≥": r"\ensuremath{\geq}", "≤": r"\ensuremath{\leq}", "×": r"\ensuremath{\times}",
       "→": r"\ensuremath{\rightarrow}", "≈": r"\ensuremath{\approx}", "±": r"\ensuremath{\pm}",
       "—": "---", "–": "--", "‰": r"\textperthousand{}",
       "“": "``", "”": "''", "‘": "`", "’": "'", "⟦": "[[", "⟧": "]]", "×": r"\ensuremath{\times}"}
SPECIAL = {"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_", "$": r"\$",
           "{": r"\{", "}": r"\}", "<": r"\textless{}", ">": r"\textgreater{}", "~": r"\textasciitilde{}",
           "^": r"\textasciicircum{}"}


def esc(t):
    for k, v in UNI.items():
        t = t.replace(k, v)
    return "".join(SPECIAL.get(c, c) for c in t)


def inline(t):
    """Escape then re-inject code/bold with literal braces via placeholders."""
    spans = []
    def stash(m, wrap):
        spans.append((wrap, m.group(1)))
        return f"\x00{len(spans)-1}\x00"
    t = re.sub(r"`([^`]+)`", lambda m: stash(m, "texttt"), t)
    t = re.sub(r"\*\*([^*]+)\*\*", lambda m: stash(m, "textbf"), t)
    t = esc(t)
    for i, (wrap, inner) in enumerate(spans):
        t = t.replace(f"\x00{i}\x00", "\\%s{%s}" % (wrap, esc(inner)))
    return t


def main():
    lines = open(SRC, encoding="utf-8").read().split("\n")
    body, i, title = [], 0, "Manuscript"
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("# ") and title == "Manuscript":
            title = inline(ln[2:]); i += 1; continue
        if ln.strip() == "---":
            i += 1; continue
        if ln.startswith("> "):  # blockquote note -> small italic
            body.append(r"\begin{quote}\small\itshape " + inline(ln[2:]) + r"\end{quote}"); i += 1; continue
        if ln.startswith("### "):
            body.append(r"\subsection*{%s}" % inline(ln[4:])); i += 1; continue
        if ln.startswith("## "):
            body.append(r"\section*{%s}" % inline(ln[3:])); i += 1; continue
        if ln.startswith("|"):  # pipe table block
            tbl = []
            while i < len(lines) and lines[i].startswith("|"):
                tbl.append(lines[i]); i += 1
            body.append(render_table(tbl)); continue
        if re.match(r"^\s*-\s+", ln):  # itemize block
            body.append(r"\begin{itemize}")
            while i < len(lines) and re.match(r"^\s*-\s+", lines[i]):
                body.append(r"\item " + inline(re.sub(r"^\s*-\s+", "", lines[i]))); i += 1
            body.append(r"\end{itemize}"); continue
        if re.match(r"^\d+\.\s", ln):  # numbered list (references)
            body.append(r"\begin{enumerate}")
            while i < len(lines) and re.match(r"^\d+\.\s", lines[i]):
                body.append(r"\item " + inline(re.sub(r"^\d+\.\s", "", lines[i]))); i += 1
            body.append(r"\end{enumerate}"); continue
        if ln.strip() == "":
            body.append(""); i += 1; continue
        body.append(inline(ln)); i += 1

    cells0 = "l"
    tex = PREAMBLE % {"title": title} + "\n".join(body) + "\n\\end{document}\n"
    open(OUT, "w", encoding="utf-8").write(tex)
    print(f"Wrote {OUT} ({tex.count(chr(10))} lines). Compile: pdflatex manuscript.tex")


def render_table(rows):
    cells = [[c.strip() for c in r.strip().strip("|").split("|")] for r in rows]
    header = cells[0]
    data = [r for r in cells[2:]]  # skip the |---| separator row
    ncol = len(header)
    spec = "p{%.2f\\linewidth}" % (0.92 / ncol)
    out = [r"\begin{longtable}{%s}" % ("".join(spec for _ in range(ncol))), r"\hline"]
    out.append(" & ".join(r"\textbf{%s}" % inline(h) for h in header) + r" \\ \hline\endhead")
    for r in data:
        r = (r + [""] * ncol)[:ncol]
        out.append(" & ".join(inline(c) for c in r) + r" \\")
    out += [r"\hline", r"\end{longtable}"]
    return "\n".join(out)


PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[margin=1in]{geometry}
\usepackage{longtable,booktabs,textcomp,hyperref,amsmath}
\usepackage{parskip}
\setlength{\LTleft}{0pt}\setlength{\LTright}{\fill}
\title{%(title)s}
\author{\textit{[[author list --- to be completed]]}}
\date{}
\begin{document}
\maketitle
"""

if __name__ == "__main__":
    main()
