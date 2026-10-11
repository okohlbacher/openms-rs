#!/usr/bin/env python3
"""Generate the XQuestResultXMLFile oracle inputs (lane F3a).

Upstream inputs (the class-test fixture and the TOPP xquest.xml/idXML files)
are copied into inputs/ by hand; this writes the generated edge cases x*.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "inputs")
os.makedirs(OUT, exist_ok=True)


def attrs(d):
    return " ".join(f'{k}="{v}"' for k, v in d.items())


OPXL_RESULTS = {
    "xquest_version": "OpenPepXL 1.0", "date": "2018-03-06 11:52:13", "author": "Eugen Netz",
    "tolerancemeasure_ms1": "ppm", "tolerancemeasure_ms2": "Da", "ms1tolerance": "10",
    "ms2tolerance": "0.2", "xlink_ms2tolerance": "0.3", "crosslinkername": "DSS",
    "xlinkermw": "138.0680796", "monolinkmw": "156.07864431, 155.094628715", "database": "db.fasta",
    "database_dc": "decoy.fasta", "xlinktypes": "1111", "AArequired1": "K", "AArequired2": "K",
    "cp_isotopediff": "12.075321", "enzyme_name": "Trypsin", "outputpath": "", "missed_cleavages": "2",
    "ntermxlinkable": "1", "CID_match2ndisotope": "1", "variable_mod": "Oxidation (M)",
    "fixed_mod": "Carbamidomethyl (C),Oxidation (M)", "decoy_prefix": "1", "decoy_string": "decoy_",
    "charges": "4,2,3", "run_path": "a.raw,b.raw", "nocutatxlink": "1",
}

XQ_RESULTS = {
    "xquest_version": "xquest 2.1.1", "date": "Fri Dec 18 12:28:42 2015", "enzyme_num": "1",
    "missed_cleavages": "2", "database": "db.fasta", "ms1tolerance": "10", "tolerancemeasure": "ppm",
    "ms2tolerance": "0.2", "tolerancemeasure_ms2": "Da", "xlink_ms2tolerance": "0.3",
    "variable_mod": "M,15.99491", "fixed_mod": "", "decoy_string": "decoy_", "database_dc": "",
    "monolinkmw": "156.07864431, 155.094628715", "xlinkermw": "138.0680796", "crosslinkername": "DSS",
    "cp_isotopediff": "12.075321", "ntermxlinkable": "1", "AArequired": "K",
}


def hit(**over):
    base = {
        "search_hit_rank": "1", "id": "X-Y-a1-b1", "type": "xlink", "structure": "X-Y",
        "seq1": "PEPKTIDE", "seq2": "LLKQR", "prot1": "P1,decoy_P2", "prot2": "P3",
        "topology": "a4-b3", "xlinkposition": "4,3", "Mr": "1", "mz": "1", "charge": "3",
        "xlinkermass": "138.0680796", "measured_mass": "1", "error": "0", "error_rel": "1.5",
        "num_of_matched_common_ions_alpha": "3", "num_of_matched_common_ions_beta": "2",
        "num_of_matched_xlink_ions_alpha": "1", "num_of_matched_xlink_ions_beta": "4",
        "wTIC": "0.1", "TIC": "0.2", "intsum": "3.5", "match_odds": "2.25", "score": "10.5",
    }
    for key, value in over.items():
        key = key.replace("__", ":")
        if value is None:
            base.pop(key, None)
        else:
            base[key] = value
    return f"<search_hit {attrs(base)} >\n</search_hit>"


def spectrum(hits, **over):
    base = {"spectrum": "run.light.7_run.heavy.9", "mz_precursor": "500.25", "scantype": "light_heavy",
            "charge_precursor": "3", "Mr_precursor": "1497.7", "rtsecscans": "100.5:101.25",
            "mzscans": "500.25:504.26", "scan_index_light": "7", "scan_index_heavy": "9"}
    for key, value in over.items():
        if value is None:
            base.pop(key, None)
        else:
            base[key] = value
    return f"<spectrum_search {attrs(base)} >\n" + "\n".join(hits) + "\n</spectrum_search>"


def doc(results, spectra, **over):
    r = dict(results)
    for key, value in over.items():
        if value is None:
            r.pop(key, None)
        else:
            r[key] = value
    return ('<?xml version="1.0" encoding="UTF-8"?>\n<xquest_results ' + attrs(r) + ">\n"
            + "\n".join(spectra) + "\n</xquest_results>\n")


XQ_SPEC = "aleitner_M1012_006.c.02942.02942.3_aleitner_M1012_006.c.02913.02913.3"
CASES = {
    # OpenPepXL: every link type, heavy and light scans, decoys, meta values.
    "x01_opxl_types": doc(OPXL_RESULTS, [
        spectrum([
            hit(),
            hit(search_hit_rank="2", type="intralink", seq1="KPEPKTIDE", seq2="", prot1="P1",
                prot2="", xlinkposition="1,5", score="7.25", fdr="0.01", xprophet_f="1", fdr_type="q-value"),
            hit(search_hit_rank="3", type="monolink", seq1="PEPKTIDE", seq2="", prot1="decoy_P9",
                prot2="", xlinkposition="4", xlinkermass="156.0786", score="-2.5"),
        ]),
        spectrum([hit(seq1="KLLR", seq2="AAAK", xlinkposition="1,4", score="0")],
                 spectrum="run.light.12_run.heavy.12", scan_index_light="12", scan_index_heavy="12",
                 rtsecscans="200:200", mzscans="", mz_precursor="612.5"),
    ]),
    # A cross-linker name no monolink modification contains: the source
    # prints "Set default mono-link" and takes the first candidate.
    "x02_opxl_monolink_default": doc(dict(OPXL_RESULTS, crosslinkername="NOSUCH"), [
        spectrum([hit(type="monolink", seq1="PEPKTIDE", seq2="", prot2="", xlinkposition="4",
                      xlinkermass="156.0786")]),
    ]),
    # Original xQuest: textual date, enzyme number, spectrum names split in
    # the middle, X for oxidised methionine, AArequired with ntermxlinkable.
    "x03_xquest": doc(XQ_RESULTS, [
        spectrum([hit(seq1="PEPXKTIDE", seq2="LLKXR", xlinkposition="5,3"),
                  hit(search_hit_rank="2", type="monolink", seq1="PEPKTIDE", seq2="", prot2="",
                      xlinkposition="4", xlinkermass="156.0786", score="3")],
                 spectrum=XQ_SPEC, charge_precursor="4", scan_index_light=None, scan_index_heavy=None),
        spectrum([hit(score="1")], spectrum="f_x.00100.00100.2_f_x.00100.00100.2", charge_precursor="2"),
    ]),
    "x04_xquest_single_digit_day": doc(dict(XQ_RESULTS, date="Tue Dec  8 02:03:04 2015"), [
        spectrum([hit()], spectrum=XQ_SPEC),
    ]),
    "x05_xquest_even_underscores": doc(XQ_RESULTS, [spectrum([hit()], spectrum="a_b_c.1.1.2")]),
    "x06_unsupported_type": doc(OPXL_RESULTS, [spectrum([hit(type="foo")])]),
    "x07_missing_seq1": doc(OPXL_RESULTS, [spectrum([hit(seq1=None)])]),
    "x08_unknown_enzyme": doc(dict(OPXL_RESULTS, enzyme_name="NoSuchEnzyme"), [spectrum([hit()])]),
    # Link positions at the termini decide the term specificity.
    "x09_term_spec": doc(dict(OPXL_RESULTS, AArequired1="K", AArequired2="K,N-term"), [
        spectrum([
            hit(seq1="APEPK", seq2="GLLKR", xlinkposition="1,5", score="5"),
            hit(search_hit_rank="2", seq1="KPEPA", seq2="KLLRG", xlinkposition="5,1", score="4"),
        ]),
    ]),
    # xQuest's U+2212 minus in the cross-linker mass.
    "x10_unicode_minus": doc(dict(OPXL_RESULTS, xlinkermw="−" + "1.5"), [
        spectrum([hit(xlinkermass="−" + "1.5")]),
    ]),
    # Two identical hits: std::unique keeps one.
    "x11_duplicate_hits": doc(OPXL_RESULTS, [spectrum([hit(), hit()])]),
    # More than 16 identifications of one light spectrum with tied scores:
    # the order std::sort leaves is observable.
    "x12_many_ties": doc(OPXL_RESULTS, [
        spectrum([hit(search_hit_rank=str(i + 1), id=f"H{i}", score=str([5, 3, 5, 1, 3][i % 5]))
                  for i in range(21)]),
    ]),
    # Hits of pairs sharing one light spectrum are combined and re-ranked.
    "x13_pairs_share_light": doc(OPXL_RESULTS, [
        spectrum([hit(id="A", score="2")]),
        spectrum([hit(id="B", score="9")], spectrum="run.light.7_run.heavy.11", scan_index_heavy="11",
                 rtsecscans="100.5:103"),
    ]),
    "x14_bad_xprophet_f": doc(OPXL_RESULTS, [spectrum([hit(xprophet_f="abc")])]),
    # Optional attributes absent or empty; decoy_prefix="0".
    "x15_optional_absent": doc(dict(OPXL_RESULTS, variable_mod=None, fixed_mod="", decoy_prefix="0",
                                   decoy_string=None, database_dc="", monolinkmw=""), [
        spectrum([hit(wTIC=None, TIC="", intsum=None, match_odds="")]),
    ]),
    "x16_no_hits": doc(OPXL_RESULTS, []),
}

for name, text in CASES.items():
    with open(os.path.join(OUT, name + ".xquest.xml"), "w", encoding="utf-8") as handle:
        handle.write(text)
print(len(CASES), "cases")
