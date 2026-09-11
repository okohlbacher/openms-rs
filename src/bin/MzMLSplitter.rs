// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! The `MzMLSplitter` executable. The tool itself is `openms::cli::tools::MzMLSplitter`.

fn main() -> std::process::ExitCode {
    std::process::ExitCode::from(
        openms::cli::run::<openms::cli::tools::MzMLSplitter>().as_i32() as u8,
    )
}
