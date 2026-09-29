#!/usr/bin/env python
"""Extract section bodies from the two source builders into `_sv_burden_kept_sections.json`, which
`_build_sv_burden_model_nb.py` reads.

`_build_sv_burden_influence_nb.py` and `_build_sv_burden_origin_climate_control_nb.py` produce no live
notebook any more (theirs are in archive/burden_superseded_2026-09-29/), but they hold the tested code for
sections the consolidated notebook still uses. Extracting rather than copying means the text cannot drift
into two versions. Re-run this whenever either source builder changes. `basic` env.

The scan is delimiter-aware on purpose. An earlier regex version looked for a closing triple quote alone on
its own line, which a one-line-ending assignment such as

    md_scale = \"\"\"## 4 ...
    ... All values are from the fitted model.\"\"\"

does not have -- so `md_scale` swallowed the whole of `code_scale` after it, and that Python source was
rendered as markdown in the notebook. This walks forward from the opening quote to its true match.
"""
import os
import json

H = os.path.dirname(os.path.abspath(__file__))
WANT = {
    "_build_sv_burden_influence_nb.py":
        ["code_load", "code_avp", "code_infl", "code_deld", "md_avp", "md_infl", "md_deld"],
    "_build_sv_burden_origin_climate_control_nb.py":
        ["code_scale", "md_scale"],
}
TRIPLES = ('"' * 3, "'" * 3)


def grab(src, name):
    """Body of `name = r?<triple-quote> ... <triple-quote>`, ending at the FIRST closing delimiter."""
    i = src.index("\n" + name + " = ") + 1
    j = i + len(name + " = ")
    if src[j] == "r":
        j += 1
    q = src[j:j + 3]
    assert q in TRIPLES, name + ": unexpected delimiter " + repr(q)
    start = j + 3
    if src[start] == "\n":
        start += 1
    return src[start:src.index(q, start)]


def main():
    out = {}
    for mod, names in WANT.items():
        src = open(os.path.join(H, mod)).read()
        for n in names:
            body = grab(src, n)
            # a body that still contains an assignment opening a triple quote means the scan ran past
            # the end of its own block -- exactly the failure this rewrite exists to prevent
            for t in TRIPLES:
                assert (" = r" + t) not in body and (" = " + t) not in body, \
                    n + " swallowed the following block -- extraction is wrong"
            out[n] = body
    f = os.path.join(H, "_sv_burden_kept_sections.json")
    json.dump(out, open(f, "w"), indent=1)
    print("[wrote] " + f)
    for n in sorted(out):
        print("  %-12s %6d chars" % (n, len(out[n])))


if __name__ == "__main__":
    main()
