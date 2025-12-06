#!/usr/bin/env python3
# merge_vasprun.py
# Merge split vasprun parts back into one XML by concatenating <calculation> blocks.

import argparse
import xml.etree.ElementTree as ET
from xml.sax.saxutils import escape

XML_DECL = '<?xml version="1.0" encoding="UTF-8"?>\n'

def _root_open_tag(tag: str, attrib: dict) -> str:
    if not attrib:
        return f"<{tag}>\n"
    attrs = " ".join(f'{k}="{escape(v)}"' for k, v in attrib.items())
    return f"<{tag} {attrs}>\n"

def merge_parts(outpath: str, parts: list[str]):
    if not parts:
        raise ValueError("No part files provided.")

    # Use the first part to determine root tag/attrib
    root_tag = None
    root_attrib = None

    # Open output early
    with open(outpath, "w", encoding="utf-8", newline="\n") as out:
        out.write(XML_DECL)

        # ---- pass part1: write root open tag + header children + its calculations ----
        first = parts[0]
        stack = []
        saw_first_calc_start = False
        opened_root = False

        for event, elem in ET.iterparse(first, events=("start", "end")):
            if event == "start":
                stack.append(elem)
                if len(stack) == 1:  # root
                    root_tag = elem.tag
                    root_attrib = dict(elem.attrib)
                    out.write(_root_open_tag(root_tag, root_attrib))
                    opened_root = True
                if len(stack) == 2 and elem.tag == "calculation":
                    saw_first_calc_start = True

            else:
                depth = len(stack)
                if depth == 2 and elem.tag == "calculation":
                    out.write(ET.tostring(elem, encoding="unicode"))
                    out.write("\n")
                    elem.clear()
                elif depth == 2 and (not saw_first_calc_start) and elem.tag != "calculation":
                    out.write(ET.tostring(elem, encoding="unicode"))
                    out.write("\n")
                    elem.clear()

                if depth >= 2 and elem.tag != "calculation":
                    elem.clear()
                stack.pop()

        if not opened_root:
            raise RuntimeError(f"Failed to parse root from {first}")

        # ---- pass remaining parts: only write calculations ----
        for p in parts[1:]:
            stack = []
            for event, elem in ET.iterparse(p, events=("start", "end")):
                if event == "start":
                    stack.append(elem)
                else:
                    depth = len(stack)
                    if depth == 1:
                        # sanity-check root tag matches
                        if elem.tag != root_tag:
                            raise RuntimeError(f"Root tag mismatch: {p} has <{elem.tag}>, expected <{root_tag}>")
                    if depth == 2 and elem.tag == "calculation":
                        out.write(ET.tostring(elem, encoding="unicode"))
                        out.write("\n")
                        elem.clear()

                    if depth >= 2 and elem.tag != "calculation":
                        elem.clear()
                    stack.pop()

        out.write(f"</{root_tag}>\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="vasprun_merged.xml", help="Output file")
    ap.add_argument("parts", nargs="+", help="Part files in order, e.g. vasprun_1.xml vasprun_2.xml ...")
    args = ap.parse_args()

    merge_parts(args.output, args.parts)
    print(f"OK: merged into {args.output}")

if __name__ == "__main__":
    main()
