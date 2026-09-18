// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The build flag and the message that tells a user how to undo it must agree.
//!
//! `.cargo/config.toml` builds x86_64 with `-C target-feature=+fma`, and
//! `openms::system::cpu_features::FMA_UNSUPPORTED_MESSAGE` tells a user with an
//! older processor to rebuild with `RUSTFLAGS="-C target-feature=-fma"`. That
//! command is correct *only* because `RUSTFLAGS` replaces the config's
//! `rustflags` rather than adding to them, and because `+fma` is the only entry
//! there: a second entry would be silently dropped by the very command this
//! crate prints, and the user would get a build missing a flag nobody told them
//! about.
//!
//! These tests pin that. They read the config file rather than trusting a
//! comment, and they skip when it is absent, which is what a consumer of the
//! packaged crate sees. See `docs/FMA_BUILD_FLAG.md`.

use openms::system::cpu_features::FMA_UNSUPPORTED_MESSAGE;

/// The `rustflags` line the message's opt-out command is written against.
const EXPECTED_RUSTFLAGS: &str = r#"rustflags = ["-C", "target-feature=+fma"]"#;

/// The command the message prints, which must undo exactly that.
const EXPECTED_OPT_OUT: &str =
    r#"RUSTFLAGS="-C target-feature=-fma" cargo build --release --locked"#;

/// The repository's Cargo configuration, or `None` outside a checkout.
fn cargo_config() -> Option<String> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(".cargo/config.toml");
    std::fs::read_to_string(path).ok()
}

/// Lines that are not comments and not blank, with their indentation removed.
fn significant_lines(text: &str) -> Vec<&str> {
    text.lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect()
}

#[test]
fn the_configuration_sets_exactly_one_rustflags_entry_scoped_to_x86_64() {
    let Some(text) = cargo_config() else {
        return; // Not a checkout: nothing to check, and nothing wrong.
    };
    let lines = significant_lines(&text);

    let tables: Vec<&&str> = lines.iter().filter(|line| line.starts_with('[')).collect();
    assert_eq!(
        tables,
        vec![&r#"[target.'cfg(target_arch = "x86_64")']"#],
        "the flag must stay scoped to x86_64, and nothing else may be configured \
         here without revisiting FMA_UNSUPPORTED_MESSAGE"
    );

    let rustflags: Vec<&&str> = lines
        .iter()
        .filter(|line| line.starts_with("rustflags"))
        .collect();
    assert_eq!(
        rustflags,
        vec![&EXPECTED_RUSTFLAGS],
        "RUSTFLAGS replaces these flags rather than adding to them, so a second \
         entry would be dropped by the opt-out command this crate prints"
    );
}

#[test]
fn the_message_prints_the_command_that_undoes_that_entry() {
    assert!(
        FMA_UNSUPPORTED_MESSAGE.contains(EXPECTED_OPT_OUT),
        "the message must carry the exact rebuild command"
    );
    // The opt-out is the configured flag with its sign turned round.
    assert!(EXPECTED_RUSTFLAGS.contains("target-feature=+fma"));
    assert!(EXPECTED_OPT_OUT.contains("target-feature=-fma"));
}

/// A user reading the message must be able to act on it without the source.
#[test]
fn the_message_says_which_processors_satisfy_the_requirement() {
    assert!(FMA_UNSUPPORTED_MESSAGE.contains("FMA3"));
    assert!(FMA_UNSUPPORTED_MESSAGE.contains("Intel Haswell, AMD Piledriver or newer"));
    assert!(FMA_UNSUPPORTED_MESSAGE.contains(".cargo/config.toml"));
}
