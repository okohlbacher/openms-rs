// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! `ParamCWLFile` (`OpenMS/FORMAT/ParamCWLFile.h`, core `bc9cc12`) as the
//! Release build has it, without TDL: both members throw the
//! `std::runtime_error` of `ParamCWLFile.cpp:331`. The Release build's text
//! is the one its `-write_cwl` and `-write_nested_cwl` print after
//! `Unable to initialize or run <tool>: ` (`../oracle/toppbase-completion`,
//! cases `write_cwl` and `write_nested_cwl`, retained in
//! `tests/data/topp_cli_completion` and compared through the command line in
//! `tests/topp_cli_completion.rs`). There is no upstream `ParamCWLFile_test.cpp`
//! at the pin.
#![cfg(feature = "paramxml")]

use openms::Error;
use openms::data_structures::ToolInfo;
use openms::format::param_cwl::{ParamCwlFile, TDL_UNAVAILABLE};
use openms::param::{Param, ParamValue};
use openms::system::file::TempDir;
use std::fs;
use std::path::Path;

/// The retained Release standard error of `-write_cwl` and
/// `-write_nested_cwl`, whose only line carries the writer's message.
fn release_message(case: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/topp_cli_completion")
        .join(case)
        .join("stderr.txt");
    let text = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let line = text.lines().last().unwrap().to_owned();
    let (_, message) = line
        .split_once(": ")
        .unwrap_or_else(|| panic!("{case}: {line}"));
    message.to_owned()
}

#[test]
fn both_members_refuse_as_the_release_build_without_tdl() {
    let mut param = Param::new();
    param
        .set_value("tool:1:value", ParamValue::Float(1.5), "a value", &[])
        .unwrap();
    let info = ToolInfo {
        name: "tool".into(),
        ..ToolInfo::default()
    };
    let dir = TempDir::new_in(std::env::temp_dir(), false).unwrap();
    let target = dir.path().join("tool.cwl");
    for (file, case) in [
        (ParamCwlFile::default(), "write_nested_cwl"),
        (
            ParamCwlFile {
                flat_hierarchy: true,
            },
            "write_cwl",
        ),
    ] {
        assert_eq!(release_message(case), TDL_UNAVAILABLE, "{case}");
        match file.store(&target.to_string_lossy(), &param, &info) {
            Err(Error::Unsupported(message)) => assert_eq!(message, TDL_UNAVAILABLE),
            other => panic!("{other:?}"),
        }
        let mut sink = Vec::new();
        match file.write_cwl_to_stream(&mut sink, &param, &info) {
            Err(Error::Unsupported(message)) => assert_eq!(message, TDL_UNAVAILABLE),
            other => panic!("{other:?}"),
        }
        assert!(sink.is_empty());
    }
    // The Release build opens (and so creates) the target first (CPP-350);
    // this port leaves it alone.
    assert!(!target.exists());
    assert_eq!(
        TDL_UNAVAILABLE,
        "TDL support is not available. Rebuild with -DENABLE_TDL=ON to enable this feature."
    );
}
