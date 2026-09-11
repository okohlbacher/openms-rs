// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
fn main() {
    // libxml0.3.14 predates the upstream MSVC static-link fix. These Windows
    // SDK libraries are needed by libxml2; no dependency source is patched.
    if std::env::var_os("CARGO_FEATURE_MZML_SCHEMA").is_some()
        && std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("windows")
        && std::env::var("CARGO_CFG_TARGET_ENV").as_deref() == Ok("msvc")
    {
        println!("cargo:rustc-link-lib=bcrypt");
        println!("cargo:rustc-link-lib=ws2_32");
    }
}
