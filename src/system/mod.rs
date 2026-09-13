// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native filesystem operations, runtime paths, timing, and owned temporary resources.
/// Build and runtime platform identity, from `SYSTEM/BuildInfo.h`.
pub mod build_info;
/// Filesystem queries, copies, temporaries and resource lookup, from `SYSTEM/File.h`.
pub mod file;
/// One-shot URL download, from `SYSTEM/Network.h`.
#[cfg(feature = "network")]
pub mod network;
/// Synchronous HTTP GET, from `SYSTEM/NetworkGetRequest.h` and `SYSTEM/CurlInit.h`.
#[cfg(feature = "network")]
pub mod network_get_request;
/// Lexical basename and path conversion, from `SYSTEM/PathUtils.h`.
pub mod path_utils;
/// Wall-clock and process CPU timing, from `SYSTEM/StopWatch.h`.
pub mod stop_watch;
/// Process and system memory reporting, from `SYSTEM/SysInfo.h`.
pub mod sys_info;
/// Rate-limited version query against the OpenMS REST server, from `SYSTEM/UpdateCheck.h`.
#[cfg(feature = "network")]
pub mod update_check;
