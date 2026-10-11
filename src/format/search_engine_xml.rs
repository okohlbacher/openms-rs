// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Resource ceilings of the X! Tandem and OMSSA XML readers.
//!
//! The source readers (`XTandemXMLFile`, `OMSSAXMLFile`) are Xerces SAX2
//! handlers derived from `Internal::XMLHandler` (`FORMAT/HANDLERS/XMLHandler.h`).
//! They are driven by the crate's one Xerces-faithful SAX driver,
//! `format::sax_xml`, which delivers character data where Xerces splits it
//! (at every reference) and white space inside the root, as several of their
//! rules depend on it (OMSSA resets the current modification on the white
//! space inside `<MSHits_mods>`; both act on the first `characters` call of a
//! text node only). Lane F3b wrote its own driver here first; wave-4
//! integration folded it into `sax_xml`, keeping this public limits type.

use crate::format::sax_xml::SaxLimits;
use crate::{Error, Result};

/// Resource ceilings for one search-engine XML document.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct XmlReadLimits {
    /// Maximum input size in bytes.
    pub max_bytes: usize,
    /// Maximum XML events, text and markup alike.
    pub max_events: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum bytes of one run of character data after reference expansion.
    pub max_text_bytes: usize,
    /// Maximum attributes on one element.
    pub max_attributes: usize,
}

impl Default for XmlReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 1024 * 1024 * 1024,
            max_events: 200_000_000,
            max_depth: 256,
            max_text_bytes: 64 * 1024 * 1024,
            max_attributes: 1024,
        }
    }
}

impl XmlReadLimits {
    /// The driver's limits.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidValue`] when a limit is zero.
    pub(crate) fn sax(&self) -> Result<SaxLimits> {
        if self.max_bytes == 0
            || self.max_events == 0
            || self.max_depth == 0
            || self.max_text_bytes == 0
            || self.max_attributes == 0
        {
            return Err(Error::InvalidValue(
                "search-engine XML limits must be positive".into(),
            ));
        }
        Ok(SaxLimits {
            max_input_bytes: u64::try_from(self.max_bytes).unwrap_or(u64::MAX),
            max_depth: self.max_depth,
            max_events: u64::try_from(self.max_events).unwrap_or(u64::MAX),
            max_attributes: self.max_attributes,
            max_text_bytes: self.max_text_bytes,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_limits_are_refused() {
        assert!(XmlReadLimits::default().sax().is_ok());
        let limits = XmlReadLimits {
            max_depth: 0,
            ..XmlReadLimits::default()
        };
        assert!(matches!(limits.sax(), Err(Error::InvalidValue(_))));
    }
}
