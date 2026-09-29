// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A lightweight view over the attribute list of an XML element
//! (`FORMAT/HANDLERS/XMLAttributes.h`).
//!
//! The source wraps a live `xercesc::Attributes`, owned by the SAX parser and
//! valid only during the `startElement` callback, behind an opaque pointer so
//! that no Xerces type appears in a public header. The Rust XML layer parses
//! with `quick-xml` and gives each handler its attributes as decoded
//! `(qualified name, value)` pairs; this type is the same small, read-only
//! surface over such a list. Its lifetime ties it to the list it views, which
//! is what the source's "must not outlive the callback" rule asks for, now
//! checked by the compiler.
//!
//! Values are UTF-8 `&str` rather than UTF-16 `const char16_t*`; where UTF-16
//! is needed, [`StringManager::convert`](super::string_manager::StringManager::convert)
//! transcodes. The source's "`nullptr` if absent" is `None`, and its negative
//! index for a missing attribute is `None` as well.

/// Non-owning view over an element's attributes: value by name, index by
/// name, value by position.
///
/// Copying the view copies the reference, as the source's defaulted copy
/// constructor and assignment copy the pointer.
///
/// # Example
///
/// ```
/// use openms::format::xml_attributes::XMLAttributes;
///
/// let list = vec![
///     ("id".to_owned(), "spectrum=1".to_owned()),
///     ("defaultArrayLength".to_owned(), "42".to_owned()),
/// ];
/// let attributes = XMLAttributes::new(&list);
/// assert_eq!(attributes.value("defaultArrayLength"), Some("42"));
/// assert_eq!(attributes.index("id"), Some(0));
/// assert_eq!(attributes.value("missing"), None);
/// assert_eq!(attributes.value_by_index(1), Some("42"));
/// assert_eq!(attributes.value_by_index(2), None);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct XMLAttributes<'a> {
    attributes: &'a [(String, String)],
}

impl<'a> XMLAttributes<'a> {
    /// Wrap an attribute list of `(qualified name, value)` pairs in document
    /// order (non-owning).
    pub fn new(attributes: &'a [(String, String)]) -> Self {
        Self { attributes }
    }

    /// The value of the attribute whose qualified name is `qname`, or `None`
    /// when absent.
    ///
    /// Names are compared exactly, prefix included, as `Attributes::getValue`
    /// compares qualified names. Should a list name one attribute twice, which
    /// a well-formed document cannot, the first is found. This is the
    /// source's `value(const char*)`; the UTF-16 overload is
    /// [`XMLAttributes::value_utf16`].
    pub fn value(&self, qname: &str) -> Option<&'a str> {
        self.index(qname).map(|i| self.attributes[i].1.as_str())
    }

    /// [`XMLAttributes::value`] for a name given as UTF-16 code units, up to
    /// the first `0` unit: the source's `value(const char16_t*)`. A name
    /// that is not valid UTF-16 matches nothing.
    pub fn value_utf16(&self, qname: &[u16]) -> Option<&'a str> {
        let name = super::string_manager::StringManager::convert_utf16(qname).ok()?;
        self.value(&name)
    }

    /// The position of the attribute whose qualified name is `qname`, or
    /// `None` where the source returns a negative index.
    pub fn index(&self, qname: &str) -> Option<usize> {
        self.attributes.iter().position(|(name, _)| name == qname)
    }

    /// [`XMLAttributes::index`] for a name given as UTF-16 code units, the
    /// source's only `index` overload.
    pub fn index_utf16(&self, qname: &[u16]) -> Option<usize> {
        let name = super::string_manager::StringManager::convert_utf16(qname).ok()?;
        self.index(&name)
    }

    /// The value of the attribute at position `i`, or `None` when `i` is out
    /// of range, where `Attributes::getValue(XMLSize_t)` returns a null
    /// pointer.
    pub fn value_by_index(&self, i: usize) -> Option<&'a str> {
        self.attributes.get(i).map(|(_, value)| value.as_str())
    }

    /// The underlying attribute list: the source's opaque `handle()`.
    pub fn handle(&self) -> &'a [(String, String)] {
        self.attributes
    }

    /// The number of attributes. Native; the source's view exposes no count.
    pub fn len(&self) -> usize {
        self.attributes.len()
    }

    /// Whether the element has no attributes. Native.
    pub fn is_empty(&self) -> bool {
        self.attributes.is_empty()
    }
}
