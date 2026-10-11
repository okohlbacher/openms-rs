// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Mass shift patterns of labelled peptide multiplets.
//!
//! Native equivalent of `FEATUREFINDER/MultiplexDeltaMasses.h`. See
//! `docs/FEATUREFINDER_MULTIPLEX_SUPPORT.md`.

use std::cmp::Ordering;

/// Set of labels associated with a mass shift (source
/// `MultiplexDeltaMasses::LabelSet`, a `std::multiset<std::string>`).
///
/// For example, a set of SILAC labels `[Lys8, Lys8, Arg10]` results in a
/// +26 Da mass shift. As in a multiset, a label may occur several times, and
/// the labels are kept in ascending byte order, which is the order in which
/// `std::string`'s `operator<` (and so the source's iteration) visits them.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct LabelSet {
    labels: Vec<String>,
}

impl LabelSet {
    /// An empty label set.
    pub fn new() -> Self {
        Self::default()
    }

    /// A label set holding the single label `label`.
    pub fn single(label: impl Into<String>) -> Self {
        let mut set = Self::new();
        set.insert(label);
        set
    }

    /// Insert `label` (source `std::multiset::insert`).
    ///
    /// Equal labels are kept, after the ones already present, as the multiset
    /// inserts at the upper bound of the equal range.
    pub fn insert(&mut self, label: impl Into<String>) {
        let label = label.into();
        let at = self
            .labels
            .partition_point(|present| present.as_str() <= label.as_str());
        self.labels.insert(at, label);
    }

    /// The labels in ascending order, repeated labels repeated.
    pub fn iter(&self) -> impl Iterator<Item = &str> {
        self.labels.iter().map(String::as_str)
    }

    /// The labels in ascending order.
    pub fn as_slice(&self) -> &[String] {
        &self.labels
    }

    /// The number of labels, repeated labels counted repeatedly.
    pub fn len(&self) -> usize {
        self.labels.len()
    }

    /// Whether the set holds no label.
    pub fn is_empty(&self) -> bool {
        self.labels.is_empty()
    }

    /// Number of occurrences of `label`.
    pub fn count(&self, label: &str) -> usize {
        self.labels.iter().filter(|l| l.as_str() == label).count()
    }
}

impl<S: Into<String>> FromIterator<S> for LabelSet {
    fn from_iter<I: IntoIterator<Item = S>>(iter: I) -> Self {
        let mut set = Self::new();
        for label in iter {
            set.insert(label);
        }
        set
    }
}

/// A mass shift with its corresponding label set (source
/// `MultiplexDeltaMasses::DeltaMass`).
#[derive(Clone, Debug, PartialEq)]
pub struct DeltaMass {
    /// The mass shift in Da (source `delta_mass`).
    pub delta_mass: f64,
    /// The labels causing the mass shift (source `label_set`).
    pub label_set: LabelSet,
}

impl DeltaMass {
    /// A delta mass with the given label set (source
    /// `DeltaMass(double dm, LabelSet ls)`).
    pub fn new(delta_mass: f64, label_set: LabelSet) -> Self {
        Self {
            delta_mass,
            label_set,
        }
    }

    /// A delta mass whose label set holds the single label `label` (source
    /// `DeltaMass(double dm, const std::string& l)`).
    pub fn with_label(delta_mass: f64, label: impl Into<String>) -> Self {
        Self::new(delta_mass, LabelSet::single(label))
    }
}

/// Mass shift pattern of a group of labelled peptides (source
/// `MultiplexDeltaMasses`).
///
/// Groups of labelled peptides appear with characteristic mass shifts. For an
/// Arg6-labelled SILAC peptide pair we expect mass shifts of 0 and 6 Da. For a
/// dimethyl-labelled pair with a single lysine we expect 56 Da and 64 Da:
/// 28 Da (N-term) + 28 Da (K) and 34 Da (N-term) + 34 Da (K) for the light and
/// heavy partners respectively.
///
/// The structure stores the mass shifts and corresponding labels of a group of
/// matching peptide features, including the zero mass shift of the first
/// peptide.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MultiplexDeltaMasses {
    /// The mass shifts between peptides, including the zero mass shift of the
    /// first peptide (source `getDeltaMasses`, whose non-const overload hands
    /// out the mutable vector; the public field is that access).
    pub delta_masses: Vec<DeltaMass>,
}

impl MultiplexDeltaMasses {
    /// An empty pattern (source default constructor).
    pub fn new() -> Self {
        Self::default()
    }

    /// A pattern holding `delta_masses` (source
    /// `MultiplexDeltaMasses(const std::vector<DeltaMass>&)`).
    pub fn from_delta_masses(delta_masses: Vec<DeltaMass>) -> Self {
        Self { delta_masses }
    }

    /// The mass shifts (source `getDeltaMasses() const`).
    pub fn delta_masses(&self) -> &[DeltaMass] {
        &self.delta_masses
    }

    /// Converts a label set to a string (source `labelSetToString`): the labels
    /// in ascending order, separated by single spaces.
    pub fn label_set_to_string(label_set: &LabelSet) -> String {
        let mut text = String::new();
        for (i, label) in label_set.iter().enumerate() {
            if i > 0 {
                text.push(' ');
            }
            text.push_str(label);
        }
        text
    }

    /// The source's free `operator<(const MultiplexDeltaMasses&, const
    /// MultiplexDeltaMasses&)`, the order in which patterns are searched.
    ///
    /// A pattern with more mass shifts comes first ("complete multiplets, then
    /// knock-out cases"). Between equally long patterns, the first position at
    /// which the shifts relative to each pattern's first shift differ decides,
    /// the smaller relative shift first ("cases without miscleavages first").
    /// Equal patterns are not less than each other. The label sets play no
    /// part.
    ///
    /// This is deliberately not [`PartialOrd`]: two patterns that differ only in
    /// their labels are equivalent here but not equal under [`PartialEq`], and a
    /// relative shift of NaN makes the source's `!=` test true and its `<` test
    /// false, which this reproduces.
    pub fn source_less(&self, other: &Self) -> bool {
        let (a, b) = (&self.delta_masses, &other.delta_masses);
        if a.len() != b.len() {
            return a.len() > b.len();
        }
        for i in 0..a.len() {
            let ms1 = a[i].delta_mass - a[0].delta_mass;
            let ms2 = b[i].delta_mass - b[0].delta_mass;
            if ms1 != ms2 {
                return ms1 < ms2;
            }
        }
        false
    }

    /// [`Self::source_less`] as an [`Ordering`], equivalent patterns comparing
    /// [`Ordering::Equal`]. Stable sorts by this key keep equivalent patterns in
    /// their input order; the source's `std::sort` does not promise that, which
    /// is why the generator sorts with [`crate::math::source_sort`] instead.
    pub fn source_cmp(&self, other: &Self) -> Ordering {
        if self.source_less(other) {
            Ordering::Less
        } else if other.source_less(self) {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }
}
