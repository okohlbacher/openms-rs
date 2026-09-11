// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Caller-owned source-compatible ID generation and the common ID value interface.

use crate::chemistry::decoy_generator::DecoyRandom;
use crate::{Error, Result};
use std::time::{SystemTime, UNIX_EPOCH};

/// Deterministic MT19937-64 draws. Instances replace the source process singleton;
/// callers can share one behind a standard Mutex when a shared sequence is needed.
/// These identifiers are not cryptographic randomness or a uniqueness guarantee.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UniqueIdGenerator {
    seed: u64,
    random: DecoyRandom,
}
impl Default for UniqueIdGenerator {
    fn default() -> Self {
        let micros = match SystemTime::now().duration_since(UNIX_EPOCH) {
            Ok(duration) => duration.as_micros() as u64,
            Err(error) => 0u64.wrapping_sub(error.duration().as_micros() as u64),
        };
        Self::from_seed(micros ^ (u64::from(std::process::id()) << 32))
    }
}
impl UniqueIdGenerator {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn from_seed(seed: u64) -> Self {
        Self {
            seed,
            random: DecoyRandom::seeded(seed),
        }
    }
    pub fn seed(&self) -> u64 {
        self.seed
    }
    pub fn set_seed(&mut self, seed: u64) {
        self.random.reseed(seed);
        self.seed = seed;
    }
    /// One raw engine word, including zero if the engine happens to emit it.
    pub fn get_unique_id(&mut self) -> u64 {
        self.random.next_u64()
    }
    /// Source UUIDv4 layout from two raw words in native byte order. The raw ID
    /// stream is portable; UUID byte order follows the source host's endianness.
    pub fn get_uuid(&mut self) -> String {
        let mut bytes = [0u8; 16];
        bytes[..8].copy_from_slice(&self.get_unique_id().to_ne_bytes());
        bytes[8..].copy_from_slice(&self.get_unique_id().to_ne_bytes());
        bytes[6] = (bytes[6] & 0x0f) | 0x40;
        bytes[8] = (bytes[8] & 0x3f) | 0x80;
        const HEX: &[u8; 16] = b"0123456789abcdef";
        let mut result = String::with_capacity(36);
        for (i, byte) in bytes.into_iter().enumerate() {
            if matches!(i, 4 | 6 | 8 | 10) {
                result.push('-');
            }
            result.push(HEX[usize::from(byte >> 4)] as char);
            result.push(HEX[usize::from(byte & 15)] as char);
        }
        result
    }
}

/// Native standalone counterpart of the concrete source UniqueIdInterface value.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct UniqueId(pub u64);
impl UniqueId {
    pub const INVALID: u64 = 0;
    pub const fn is_valid(value: u64) -> bool {
        value != Self::INVALID
    }
}

/// Implements source ID operations over an existing u64 field, without another
/// wrapper allocation or a second stored ID. Generator ownership is explicit.
pub trait HasUniqueId {
    fn unique_id(&self) -> u64;
    fn unique_id_mut(&mut self) -> &mut u64;

    fn has_valid_unique_id(&self) -> bool {
        UniqueId::is_valid(self.unique_id())
    }
    fn has_invalid_unique_id(&self) -> bool {
        !self.has_valid_unique_id()
    }
    fn set_unique_id(&mut self, value: u64) {
        *self.unique_id_mut() = value;
    }
    fn clear_unique_id(&mut self) -> usize {
        let changed = usize::from(self.has_valid_unique_id());
        self.set_unique_id(0);
        changed
    }
    fn swap_unique_id<T: HasUniqueId>(&mut self, other: &mut T) {
        std::mem::swap(self.unique_id_mut(), other.unique_id_mut());
    }
    /// Always reports one assignment; as in C++, zero is not redrawn.
    fn assign_new_unique_id(&mut self, generator: &mut UniqueIdGenerator) -> usize {
        self.set_unique_id(generator.get_unique_id());
        1
    }
    fn ensure_unique_id(&mut self, generator: &mut UniqueIdGenerator) -> usize {
        if self.has_invalid_unique_id() {
            self.assign_new_unique_id(generator)
        } else {
            0
        }
    }
    /// Last underscore suffix, ASCII digits only, with defined u64 wrapping.
    /// Invalid/empty suffix clears the ID. Inputs over 1 MiB return an error
    /// before changing it; all other text follows the source parser exactly.
    fn set_unique_id_from_str(&mut self, text: &str) -> Result<()> {
        if text.len() > 1024 * 1024 {
            return Err(Error::InvalidValue("unique ID text exceeds 1 MiB".into()));
        }
        let suffix = text.rsplit('_').next().unwrap_or("");
        let mut value = 0u64;
        for byte in suffix.bytes() {
            if !byte.is_ascii_digit() {
                self.set_unique_id(0);
                return Ok(());
            }
            value = value.wrapping_mul(10).wrapping_add(u64::from(byte - b'0'));
        }
        self.set_unique_id(value);
        Ok(())
    }
}
impl HasUniqueId for UniqueId {
    fn unique_id(&self) -> u64 {
        self.0
    }
    fn unique_id_mut(&mut self) -> &mut u64 {
        &mut self.0
    }
}
impl HasUniqueId for u64 {
    fn unique_id(&self) -> u64 {
        *self
    }
    fn unique_id_mut(&mut self) -> &mut u64 {
        self
    }
}
