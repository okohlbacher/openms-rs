// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ports of `ANALYSIS/ID`: peptide and protein identification post-processing: consensus identification, score handling, FDR and Percolator support, search engines and metabolite identification. Older ID ports (`false_discovery_rate`, `peptide_indexing`, `id_filter`, …) predate this module and stay where they are.
//!
//! Submodules are registered below, one per ported header or cohesive group of
//! headers. Parallel porting lanes each append their own `pub mod` lines here.
