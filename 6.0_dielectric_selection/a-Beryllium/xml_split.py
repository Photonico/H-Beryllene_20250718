#!/usr/bin/env python3
# split_vasprun.py
# Split vasprun.xml into N well-formed XML files by <calculation> blocks.
# usage: python3 split_vasprun.py vasprun.xml --parts 4 --outdir . --prefix vasprun_

import argparse
import os
import sys
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

XML_DECL = '<?xml version="1.0" encoding="UTF-8"?>\n'

def _root_open_tag(tag: str, attrib: dict) -> str:
    if not attrib:
        return f"<{tag}>\n"
    attrs = " ".join(f'{k}="{escape(v)}"' for k, v in attrib.items())
    return f"<{tag} {attrs}>\n"

def _count_calcs_and_header(inpath: str):
    """Return (root_tag, root_attrib, header_children_xml:list[str], calc_count:int)"""
    header_xml = []
    calc_count = 0
    root_tag = None
    root_attrib = None

    stack = []
    saw_first_calc_start = False

    # iterparse is streaming; we clear elements to keep memory low
    for event, elem in ET.iterparse(inpath, events=("start", "end")):
        if event == "start":
            stack.append(elem)
            if len(stack) == 1:
                root_tag = elem.tag
                root_attrib = dict(elem.attrib)
            # first calculation starts once we see <calculation> as a direct child of root
            if len(stack) == 2 and elem.tag == "calculation":
                saw_first_calc_start = True

        else:  # end
            depth = len(stack)
            if depth == 2 and elem.tag == "calculation":
                calc_count += 1
                elem.clear()
            elif depth == 2 and (not saw_first_calc_start) and elem.tag != "calculation":
                # header child under root before any calculation
                header_xml.append(ET.tostring(elem, encoding="unicode"))
                elem.clear()

            stack.pop()

    if root_tag is None:
        raise RuntimeError("Failed to detect root tag; is this a valid XML?")
    return root_tag, root_attrib, header_xml, calc_count

def split_vasprun(inpath: str, outdir: str, prefix: str, parts: int):
    if parts < 2:
        raise ValueError("--parts must be >= 2")

    root_tag, root_attrib, header_xml, ncalc = _count_calcs_and_header(inpath)
    if ncalc == 0:
        raise RuntimeError("No <calculation> blocks found; cannot split.")

    # distribute calculations across parts as evenly as possible
    base = ncalc // parts
    rem = ncalc % parts
    sizes = [base + (1 if i < rem else 0) for i in range(parts)]
    boundaries = []
    s = 0
    for sz in sizes:
        boundaries.append((s + 1, s + sz))  # inclusive calc indices (1-based)
        s += sz

    os.makedirs(outdir, exist_ok=True)

    def open_part(i: int):
        outpath = os.path.join(outdir, f"{prefix}{i}.xml")
        f = open(outpath, "w", encoding="utf-8", newline="\n")
        f.write(XML_DECL)
        f.write(_root_open_tag(root_tag, root_attrib))
        for hx in header_xml:
            f.write(hx)
            if not hx.endswith("\n"):
                f.write("\n")
        return f, outpath

    current_part = 1
    part_file, part_path = open_part(current_part)
    start_i, end_i = boundaries[current_part - 1]
    written = 0
    calc_idx = 0

    stack = []

    for event, elem in ET.iterparse(inpath, events=("start", "end")):
        if event == "start":
            stack.append(elem)
        else:
            depth = len(stack)
            if depth == 2 and elem.tag == "calculation":
                calc_idx += 1
                if start_i <= calc_idx <= end_i:
                    part_file.write(ET.tostring(elem, encoding="unicode"))
                    part_file.write("\n")
                    written += 1

                elem.clear()

                if calc_idx == end_i and current_part < parts:
                    # close current part
                    part_file.write(f"</{root_tag}>\n")
                    part_file.close()

                    # open next part
                    current_part += 1
                    part_file, part_path = open_part(current_part)
                    start_i, end_i = boundaries[current_part - 1]

            # clear non-root elements aggressively to keep memory down
            if depth >= 2 and elem.tag != "calculation":
                elem.clear()

            stack.pop()

    # close last file
    part_file.write(f"</{root_tag}>\n")
    part_file.close()

    return ncalc, sizes

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Path to vasprun.xml")
    ap.add_argument("--parts", type=int, default=4, help="Number of parts (default: 4)")
    ap.add_argument("--outdir", default=".", help="Output directory (default: .)")
    ap.add_argument("--prefix", default="vasprun_", help="Output prefix (default: vasprun_)")
    args = ap.parse_args()

    ncalc, sizes = split_vasprun(args.input, args.outdir, args.prefix, args.parts)
    print(f"OK: found {ncalc} <calculation> blocks.")
    for i, sz in enumerate(sizes, 1):
        print(f"  wrote {args.prefix}{i}.xml with {sz} calculations")

if __name__ == "__main__":
    main()
