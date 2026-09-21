// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The FileInfo summary of FASTA sequence databases
//! (`FORMAT/FileInfo.cpp:853-1076`, `:2001-2004`, `:2112-2114`, `:2377-2379`).
//!
//! The FASTA branch of the report. The entries are loaded through
//! [`FASTAFile::load`](crate::format::fasta::FASTAFile::load), and the branch then writes, in the source order:
//!
//! 1. the number of sequences;
//! 2. the sequence-length distribution, five lines, but only when there is at
//!    least one sequence. Below three sequences the source does not ask for
//!    quartiles and prints the minimum and the maximum in their place;
//! 3. the number of sequences carrying at least one ambiguous residue, and the
//!    duplicate-header and duplicate-sequence counts, each with a percentage
//!    of the entry count at two decimal places;
//! 4. the total residue count and one line per residue byte with its count;
//! 5. two ambiguity totals, which alphabet they use depending on whether the
//!    file was recognised as nucleic acid.
//!
//! Whether the file is nucleic acid is decided first, over every byte of every
//! sequence: one byte outside the IUPAC nucleotide alphabet makes the whole
//! file amino acid, and that choice selects the labels, the ambiguity
//! alphabets and the two trailing totals.
//!
//! The branch writes nothing at all to the TSV report, and the `-m`, `-p` and
//! `-s` arms for FASTA are empty in the source (`:2001-2004`, `:2112-2114`,
//! `:2377-2379`), so `-m` and `-s` contribute only their titles and `-p` its
//! title and the no-information line.
//!
//! The structured [`FastaInfo`](crate::format::file_info::model::FastaInfo) is filled alongside, as the source fills its
//! `Result`; the duplicate warnings the source writes with `OPENMS_LOG_WARN`
//! go into [`FileInfoResult::warnings`](crate::format::file_info::model::FileInfoResult::warnings), which is where this port collects the
//! messages the source keeps out of the reports.
//!
//! See `docs/FILE_INFO_A7_SUPPORT.md` for the evidence and the native
//! differences.

use super::model::{FastaInfo, FileInfoResult, Options};
use super::report::{
    ReportStream, write_meta_title, write_processing, write_processing_title,
    write_statistics_title,
};
use crate::concept::math_functions::percent_of;
use crate::format::fasta::{FASTAEntry, FASTAFile};
use crate::math::statistic_functions::{
    SummaryStatistics, mean, median_sorted, quantile1st_sorted, quantile3rd_sorted,
    variance_with_mean,
};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// `AA_AMBIGUOUS_BXZJ` (`FORMAT/FileInfo.cpp:890`): B = Asx, Z = Glx, X = unknown,
/// J = Leu/Ile, in both cases.
const AA_AMBIGUOUS_BXZJ: &[u8] = b"BZXbzxJj";

/// `AA_AMBIGUOUS_BXZ` (`FORMAT/FileInfo.cpp:891`), the same without J.
const AA_AMBIGUOUS_BXZ: &[u8] = b"BZXbzx";

/// `NUCLEOTIDE_CHARS` (`FORMAT/FileInfo.cpp:895`): the standard codes A, C, G, T and
/// U and the ambiguity codes, in both cases. A sequence byte outside this set
/// makes the whole file amino acid.
const NUCLEOTIDE_CHARS: &[u8] = b"ACGTUNacgtunRYSWKMBDHVryswkmbdhv";

/// `NA_AMBIGUOUS` (`FORMAT/FileInfo.cpp:896`): every IUPAC nucleotide ambiguity code,
/// in both cases.
const NA_AMBIGUOUS: &[u8] = b"NRYSWKMBDHVnryswkmbdhv";

/// The two residue bytes the nucleic-acid `N` total adds up.
const NA_AMBIGUOUS_N: &[u8] = b"Nn";

/// Load the entries and write the FASTA branch, then the `-m`, `-p` and `-s`
/// titles the source still emits for a FASTA input.
pub(crate) fn report(
    path: &Path,
    options: &Options,
    os: &mut ReportStream,
    os_tsv: &mut ReportStream,
    result: &mut FileInfoResult,
) -> Result<()> {
    let mut file = FASTAFile::new();
    let entries = file.load(path)?;
    reject_non_ascii(&entries)?;

    let is_nucleic_acid = is_nucleic_acid(&entries);
    let ambiguous = if is_nucleic_acid {
        NA_AMBIGUOUS
    } else {
        AA_AMBIGUOUS_BXZJ
    };

    // std::map<char, int>, keyed by the raw sequence byte.
    let mut residue_counts: BTreeMap<u8, u64> = BTreeMap::new();
    let mut number_of_residues = 0_u64;
    let mut dup_header = 0_u64;
    let mut dup_seq = 0_u64;
    let mut seq_has_ambiguous = 0_u64;
    let mut sequence_lengths: Vec<u64> = Vec::new();
    sequence_lengths
        .try_reserve_exact(entries.len())
        .map_err(|_| overflow("cannot allocate the FASTA sequence-length buffer"))?;

    // The source's two SHashmaps. Each bucket is ASSIGNED a one-element vector
    // rather than appended to (`m_headers[id_hash] = { index };`), so only the
    // last index with a given hash survives; see `last_by_hash`.
    let mut headers_by_hash: BTreeMap<u64, usize> = BTreeMap::new();
    let mut sequences_by_hash: BTreeMap<u64, usize> = BTreeMap::new();

    for (index, entry) in entries.iter().enumerate() {
        if let Some(previous) =
            headers_by_hash.insert(string_hash(entry.identifier.as_bytes()), index)
        {
            if entries[previous].header_matches(entry) {
                dup_header += 1;
                result.warnings.push(format!(
                    "Warning: Duplicate header, #{index}, ID: {} = #{previous}, ID: {}",
                    entry.identifier, entries[previous].identifier
                ));
            }
        }
        if let Some(previous) =
            sequences_by_hash.insert(string_hash(entry.sequence.as_bytes()), index)
        {
            if entries[previous].sequence_matches(entry) {
                dup_seq += 1;
                result.warnings.push(format!(
                    "Warning: Duplicate sequence, #{index}, ID: {} == #{previous}, ID: {}",
                    entry.identifier, entries[previous].identifier
                ));
            }
        }

        let sequence = entry.sequence.as_bytes();
        sequence_lengths.push(count(sequence.len())?);

        // The source compares the cumulative ambiguous total before and after
        // this sequence's residues are added, so a sequence is counted once
        // however many ambiguous residues it has.
        let before = count_residues(&residue_counts, ambiguous);
        for &byte in sequence {
            *residue_counts.entry(byte).or_insert(0) += 1;
        }
        let after = count_residues(&residue_counts, ambiguous);
        if before != after {
            seq_has_ambiguous += 1;
        }
        number_of_residues = number_of_residues
            .checked_add(count(sequence.len())?)
            .ok_or_else(|| overflow("FASTA residue count overflows 64 bits"))?;
    }

    let (residue_type, residue_type_cap) = if is_nucleic_acid {
        ("nucleotide", "Nucleotide")
    } else {
        ("amino acid", "Amino acid")
    };
    let entry_count = count(entries.len())?;

    os.text("\nNumber of sequences   : ")
        .value(entry_count)
        .text("\n");

    if !sequence_lengths.is_empty() {
        sequence_lengths.sort_unstable();
        let lengths: Vec<f64> = sequence_lengths.iter().map(|&n| as_double(n)).collect();
        let len_min = sequence_lengths[0];
        let len_max = sequence_lengths[sequence_lengths.len() - 1];
        let len_median = median_sorted(&lengths)?;
        let (len_q1, len_q3) = if sequence_lengths.len() >= 3 {
            (quantile1st_sorted(&lengths)?, quantile3rd_sorted(&lengths)?)
        } else {
            (as_double(len_min), as_double(len_max))
        };
        os.text("Sequence length distribution:\n");
        os.text("  Minimum : ").value(len_min).text("\n");
        os.text("  25%ile  : ").double(len_q1).text("\n");
        os.text("  Median  : ").double(len_median).text("\n");
        os.text("  75%ile  : ").double(len_q3).text("\n");
        os.text("  Maximum : ").value(len_max).text("\n");
    }

    let total = as_double(entry_count);
    os.text("Number of sequences with ambiguous ")
        .text(residue_type)
        .text("s: ")
        .value(seq_has_ambiguous)
        .text(" (")
        .double(percent_of(as_double(seq_has_ambiguous), total, 2)?)
        .text("%)\n");
    os.text("# duplicated headers  : ")
        .value(dup_header)
        .text(" (")
        .double(percent_of(as_double(dup_header), total, 2)?)
        .text("%)\n");
    os.text("# duplicated sequences: ")
        .value(dup_seq)
        .text(" (")
        .double(percent_of(as_double(dup_seq), total, 2)?)
        .text("%) [by exact string matching]\n");
    os.text("Total ")
        .text(residue_type)
        .text("s     : ")
        .value(number_of_residues)
        .text("\n\n");
    os.text(residue_type_cap).text(" counts:\n");
    for (byte, count) in signed_char_order(&residue_counts) {
        // The residue is printed as the single character the byte is; every
        // byte here is ASCII, because `reject_non_ascii` refused the rest.
        os.text("  ")
            .value(char::from(byte))
            .text(":\t")
            .value(count)
            .text("\n");
    }

    let residues = as_double(number_of_residues);
    if is_nucleic_acid {
        let amb_n = count_residues(&residue_counts, NA_AMBIGUOUS_N);
        let amb_all = count_residues(&residue_counts, NA_AMBIGUOUS);
        os.text("Ambiguous nucleotides (N)      : ")
            .value(amb_n)
            .text(" (")
            .double(percent_of(as_double(amb_n), residues, 2)?)
            .text("%)\n");
        os.text("All IUPAC ambiguity codes      : ")
            .value(amb_all)
            .text(" (")
            .double(percent_of(as_double(amb_all), residues, 2)?)
            .text("%)\n\n");
    } else {
        let amb = count_residues(&residue_counts, AA_AMBIGUOUS_BXZ);
        let amb_j = count_residues(&residue_counts, AA_AMBIGUOUS_BXZJ);
        os.text("Ambiguous amino acids (B/Z/X)  : ")
            .value(amb)
            .text(" (")
            .double(percent_of(as_double(amb), residues, 2)?)
            .text("%)\n");
        os.text("                      (B/Z/X/J): ")
            .value(amb_j)
            .text(" (")
            .double(percent_of(as_double(amb_j), residues, 2)?)
            .text("%)\n\n");
    }

    let mut ambiguity_counts = BTreeMap::new();
    if is_nucleic_acid {
        ambiguity_counts.insert(
            "N".to_owned(),
            count_residues(&residue_counts, NA_AMBIGUOUS_N),
        );
        ambiguity_counts.insert(
            "IUPAC".to_owned(),
            count_residues(&residue_counts, NA_AMBIGUOUS),
        );
    } else {
        ambiguity_counts.insert(
            "BZX".to_owned(),
            count_residues(&residue_counts, AA_AMBIGUOUS_BXZ),
        );
        ambiguity_counts.insert(
            "BZXJ".to_owned(),
            count_residues(&residue_counts, AA_AMBIGUOUS_BXZJ),
        );
    }
    result.fasta = Some(FastaInfo {
        num_sequences: entry_count,
        total_residues: number_of_residues,
        is_nucleic_acid,
        length_stats: length_stats(&sequence_lengths)?,
        // The source skips zero-count keys here, because the verbatim print
        // code above may already have operator[]-inserted 'N' and 'n'.
        residue_counts: residue_counts
            .iter()
            .filter(|&(_, &count)| count != 0)
            .map(|(&byte, &count)| (byte, count))
            .collect(),
        seq_with_ambiguous: seq_has_ambiguous,
        dup_headers: dup_header,
        dup_sequences: dup_seq,
        ambiguity_counts,
    });

    if options.meta {
        // FORMAT/FileInfo.cpp:2001-2004: the FASTA arm of the -m block is empty.
        write_meta_title(os);
    }
    if options.processing {
        // FORMAT/FileInfo.cpp:2112-2114: the FASTA arm leaves `dp` empty, so only the
        // no-information line follows the title.
        write_processing_title(os);
        write_processing(os, os_tsv, &[], result);
    }
    if options.statistics {
        // FORMAT/FileInfo.cpp:2377-2379: the FASTA arm of the -s block is empty.
        write_statistics_title(os);
    }
    Ok(())
}

/// `is_nucleic_acid` (`FORMAT/FileInfo.cpp:900-912`): true until one sequence byte
/// falls outside [`NUCLEOTIDE_CHARS`]. An empty entry list leaves it true, as
/// the source's initial value does.
fn is_nucleic_acid(entries: &[FASTAEntry]) -> bool {
    entries.iter().all(|entry| {
        entry
            .sequence
            .bytes()
            .all(|byte| NUCLEOTIDE_CHARS.contains(&byte))
    })
}

/// `count_residues` (`FORMAT/FileInfo.cpp:877-885`): the counts of `which`'s bytes
/// that the table holds, absent keys contributing nothing. The source looks
/// the keys up with `find`, so this never inserts.
fn count_residues(residue_counts: &BTreeMap<u8, u64>, which: &[u8]) -> u64 {
    which
        .iter()
        .filter_map(|byte| residue_counts.get(byte))
        .sum()
}

/// The residue table in the order the source's `std::map<char, int>` iterates
/// it. `char` is signed on the reference build's x86_64 Linux target, so a
/// byte at or above `0x80` sorts *before* `A`; the model keys the table by
/// `u8`, which orders it the other way round.
///
/// [`reject_non_ascii`] refuses every input that could tell the two orders
/// apart, so this only has to agree with the model's order on ASCII, where it
/// does. It is written as the signed comparison anyway, so that the rule the
/// source follows is the rule in the code.
fn signed_char_order(residue_counts: &BTreeMap<u8, u64>) -> Vec<(u8, u64)> {
    let mut rows: Vec<(u8, u64)> = residue_counts.iter().map(|(&b, &c)| (b, c)).collect();
    rows.sort_by_key(|&(byte, _)| byte as i8);
    rows
}

/// Refuse a sequence byte outside ASCII, before anything is written.
///
/// The source counts and prints raw `char`s, so a byte at or above `0x80` ends
/// up in the report as that single byte. This port's reports are Rust strings
/// and cannot hold a lone continuation byte, and the reader only ever yields
/// such a byte as part of a multi-byte UTF-8 sequence, so there is no faithful
/// rendering available. Refusing is the documented native difference; the
/// source reports such a file.
///
/// # Errors
///
/// [`Error::Unsupported`] naming the first offending byte and its entry.
fn reject_non_ascii(entries: &[FASTAEntry]) -> Result<()> {
    for (index, entry) in entries.iter().enumerate() {
        if let Some(byte) = entry.sequence.bytes().find(|byte| !byte.is_ascii()) {
            return Err(Error::Unsupported(format!(
                "FileInfo FASTA branch: sequence #{index} ({}) holds the non-ASCII byte {byte:#04x}, \
                 which this port's report cannot render as the single character the source prints",
                entry.identifier
            )));
        }
    }
    Ok(())
}

/// The structured length statistics (`FORMAT/FileInfo.cpp:1046-1065`): three or more
/// sequences go through the full [`SummaryStatistics`], one or two are filled
/// field by field with the quartiles falling back to the extremes, and none
/// leaves the default.
///
/// `lengths` must already be ascending, as the print code above sorted it.
fn length_stats(lengths: &[u64]) -> Result<SummaryStatistics> {
    let values: Vec<f64> = lengths.iter().map(|&n| as_double(n)).collect();
    if values.len() >= 3 {
        let mut values = values;
        return SummaryStatistics::new(&mut values);
    }
    if values.is_empty() {
        return Ok(SummaryStatistics::default());
    }
    let mean_of_lengths = mean(&values)?;
    Ok(SummaryStatistics {
        count: values.len(),
        min: values[0],
        max: values[values.len() - 1],
        mean: mean_of_lengths,
        median: median_sorted(&values)?,
        lowerq: values[0],
        upperq: values[values.len() - 1],
        variance: if values.len() >= 2 {
            variance_with_mean(&values, mean_of_lengths)?
        } else {
            0.0
        },
    })
}

/// `std::hash<std::string>` as the reference build's libstdc++ computes it:
/// `_Hash_bytes(data, length, 0xc70f6907)`, the 64-bit Murmur variant of
/// `libstdc++-v3/libsupc++/hash_bytes.cc`.
///
/// The hash is observable. FileInfo's duplicate detection ASSIGNS each bucket
/// a one-element vector instead of appending to it (`FORMAT/FileInfo.cpp:931` and
/// `:949`), so a bucket only ever remembers the last index with that hash. Two
/// different strings that collide therefore hide a duplicate that a
/// collision-free hash would have reported, and reproducing the reference
/// build needs its hash function rather than merely some hash function. No
/// crate implements this exact seed and finalisation, so it is written out
/// here; `oracle/a7-fileinfo/scripts/probe_std_hash.cpp` checks it against the
/// reference toolchain over the empty string, every tail length, bytes above
/// `0x7f` and an embedded NUL.
fn string_hash(data: &[u8]) -> u64 {
    /// `_Hash_impl::hash`'s default seed.
    const SEED: u64 = 0xc70f_6907;
    /// `(0xc6a4a793 << 32) + 0x5bd1e995`.
    const MUL: u64 = 0xc6a4_a793_5bd1_e995;

    fn shift_mix(value: u64) -> u64 {
        value ^ (value >> 47)
    }

    let length = data.len() as u64;
    let aligned = data.len() & !7;
    let mut hash = SEED ^ length.wrapping_mul(MUL);
    for chunk in data[..aligned].chunks_exact(8) {
        let mut word = [0_u8; 8];
        word.copy_from_slice(chunk);
        let value = shift_mix(u64::from_le_bytes(word).wrapping_mul(MUL)).wrapping_mul(MUL);
        hash = (hash ^ value).wrapping_mul(MUL);
    }
    if data.len() & 7 != 0 {
        // `load_bytes` reads the remaining bytes least significant first.
        let mut word = [0_u8; 8];
        word[..data.len() - aligned].copy_from_slice(&data[aligned..]);
        hash = (hash ^ u64::from_le_bytes(word)).wrapping_mul(MUL);
    }
    shift_mix(shift_mix(hash).wrapping_mul(MUL))
}

/// The `size_t` to `double` conversion the source performs when it multiplies
/// a count by `100.0` or hands it to `Math::mean`.
#[expect(
    clippy::cast_precision_loss,
    reason = "the source converts its size_t counts to double in exactly this place"
)]
fn as_double(value: u64) -> f64 {
    value as f64
}

fn count(value: usize) -> Result<u64> {
    u64::try_from(value).map_err(|_| overflow("FASTA count overflows 64 bits"))
}

fn overflow(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::string_hash;

    /// 28 of the 60 values `oracle/a7-fileinfo/scripts/probe_std_hash.cpp`
    /// printed on ibminode06 (g++ 13.3.0, the reference build's toolchain):
    /// the empty string, every `length % 8` tail case, one exactly-aligned
    /// input, the identifiers the A7 fixtures use, two bytes above `0x7f` and
    /// an embedded NUL. The probe's own record is
    /// `oracle/a7-fileinfo/scripts/probe_std_hash.cpp`; the hash matters
    /// because the source overwrites each duplicate-detection bucket instead of
    /// appending to it, so a collision hides a duplicate.
    #[test]
    fn string_hash_matches_the_reference_libstdcxx() {
        const PROBE: [(&[u8], u64); 28] = [
            (b"", 0x553e_9390_1e46_2a6e),
            (b"A", 0x6006_68de_4345_e18e),
            (b"AC", 0x365a_fb19_3177_08eb),
            (b"ACD", 0x5266_fbbb_b357_8c7a),
            (b"ACDE", 0x9e97_a9eb_486c_cae7),
            (b"ACDEF", 0x1538_0b27_6091_8a57),
            (b"ACDEFG", 0x4105_bb2d_29df_d0d6),
            (b"ACDEFGH", 0x36bf_acb5_e10c_dcac),
            (b"ACDEFGHI", 0xc326_f1c7_b6b6_7380),
            (b"ACDEFGHIK", 0x959b_3450_528c_4632),
            (b"ACDEFGHIKLMNPQRS", 0xc39f_dc69_f5e7_653d),
            (b"DUP", 0xddae_ab53_8cf9_1be7),
            (b"SAME", 0x52bb_604c_c6c3_24b4),
            (b"SEQ1", 0x6c95_a25e_3c3c_74c1),
            (b"SEQ2", 0x90bc_5527_6222_8415),
            (b"SEQ3", 0xdaeb_6802_1175_e753),
            (b"DNA1", 0x7b14_bebf_849f_067a),
            (b"DNA2", 0xaeb8_3a32_871d_d15b),
            (b"DNA3", 0x9cc1_18f8_987b_8732),
            (b"AA1", 0x682d_3f78_6dc9_2510),
            (b"AA2", 0x3acb_bcac_5774_31f9),
            (b"AA3", 0x6509_6d05_3ab0_0b09),
            (b"MIX", 0xc0a1_9812_dd37_6ede),
            (b"TINY", 0xc3d2_7402_1ab6_8fc0),
            (b"M", 0x5043_d991_72e3_9881),
            (b"MKVLmkvlATTLattl", 0x4d66_1ddb_623a_2463),
            (b"\xc3\xa9\xc3\xa9", 0xb319_0c67_004b_b667),
            (b"a\x00b", 0xad5e_8d18_a187_2c2f),
        ];
        for (input, expected) in PROBE {
            assert_eq!(
                string_hash(input),
                expected,
                "std::hash of {input:?} (length {})",
                input.len()
            );
        }
    }

    /// The hash of a long input is not the hash of its prefix: the aligned loop
    /// runs more than once.
    #[test]
    fn string_hash_consumes_every_eight_byte_block() {
        let long = b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy";
        assert_ne!(string_hash(long), string_hash(&long[..16]));
        assert_ne!(string_hash(long), string_hash(&long[..32]));
    }
}
