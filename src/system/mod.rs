// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native filesystem operations, runtime paths, timing, external processes and
//! owned temporary resources.
/// Build and runtime platform identity, from `SYSTEM/BuildInfo.h`.
pub mod build_info;
/// Starting an external program and forwarding its output, from `SYSTEM/ExternalProcess.h`.
pub mod external_process;
/// Filesystem queries, copies, temporaries and resource lookup, from `SYSTEM/File.h`.
pub mod file;
/// Detecting a Java installation, from `SYSTEM/JavaInfo.h`.
pub mod java_info;
/// Lexical basename and path conversion, from `SYSTEM/PathUtils.h`.
pub mod path_utils;
/// Detecting a Python installation, from `SYSTEM/PythonInfo.h`.
pub mod python_info;
/// Calling R scripts, from `SYSTEM/RWrapper.h`.
pub mod r_wrapper;
/// Wall-clock and process CPU timing, from `SYSTEM/StopWatch.h`.
pub mod stop_watch;
/// Process and system memory reporting, from `SYSTEM/SysInfo.h`.
pub mod sys_info;
