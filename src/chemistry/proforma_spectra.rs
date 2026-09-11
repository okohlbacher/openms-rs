// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Complete source ProForma spectrum wrapper operations; see support document.
use super::conversion::{ConversionEvaluation, ConversionWarning};
use super::resolution::{Budget, Resolver};
use super::*;
use crate::chemistry::{
    ModificationsDB, ProteinProteinCrossLink, TheoreticalIonSeries, TheoreticalSpectrumGenerator,
    TheoreticalSpectrumGeneratorXLMS,
};
use crate::{MSSpectrum, Peak1D};
use std::mem::size_of;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpectrumGenerationOptions {
    pub min_charge: i32,
    pub max_charge: i32,
    pub ion_types: String,
    pub add_losses: bool,
    pub add_metainfo: bool,
}
impl Default for SpectrumGenerationOptions {
    fn default() -> Self {
        Self {
            min_charge: 1,
            max_charge: 1,
            ion_types: "by".into(),
            add_losses: false,
            add_metainfo: true,
        }
    }
}
#[derive(Clone, Copy)]
enum Input<'a> {
    Chain(&'a Peptidoform),
    Ion(&'a PeptidoformIon),
}
struct Output {
    issues: Vec<ConversionIssue>,
    spectrum: MSSpectrum,
}
impl Peptidoform {
    pub fn spectrum_generation_issues(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<Vec<ConversionIssue>>> {
        let result = run(Input::Chain(self), None, registry, &mut Budget::default())?;
        Ok(ConversionEvaluation {
            value: result.value.issues,
            warnings: result.warnings,
        })
    }
    pub fn can_generate_spectrum(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<bool>> {
        let result = self.spectrum_generation_issues(registry)?;
        Ok(ConversionEvaluation {
            value: result.value.is_empty(),
            warnings: result.warnings,
        })
    }
    pub fn generate_spectrum(
        &self,
        options: &SpectrumGenerationOptions,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<MSSpectrum>> {
        let result = run(
            Input::Chain(self),
            Some(options),
            registry,
            &mut Budget::default(),
        )?;
        Ok(ConversionEvaluation {
            value: result.value.spectrum,
            warnings: result.warnings,
        })
    }
}
impl PeptidoformIon {
    pub fn spectrum_generation_issues(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<Vec<ConversionIssue>>> {
        let result = run(Input::Ion(self), None, registry, &mut Budget::default())?;
        Ok(ConversionEvaluation {
            value: result.value.issues,
            warnings: result.warnings,
        })
    }
    pub fn can_generate_spectrum(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<bool>> {
        let result = self.spectrum_generation_issues(registry)?;
        Ok(ConversionEvaluation {
            value: result.value.is_empty(),
            warnings: result.warnings,
        })
    }
    pub fn generate_spectrum(
        &self,
        options: &SpectrumGenerationOptions,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<MSSpectrum>> {
        let result = run(
            Input::Ion(self),
            Some(options),
            registry,
            &mut Budget::default(),
        )?;
        Ok(ConversionEvaluation {
            value: result.value.spectrum,
            warnings: result.warnings,
        })
    }
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn add_warnings(
    out: &mut Vec<ConversionWarning>,
    resolver: &mut Resolver<'_>,
    extra: Vec<ConversionWarning>,
) -> Result<()> {
    let count = resolver.warning_count().saturating_add(extra.len());
    resolver.budget.items(count)?;
    for warning in resolver
        .take_warnings()
        .into_iter()
        .map(ConversionWarning::Resolution)
        .chain(extra)
    {
        resolver.budget.reserve(out)?;
        out.push(warning);
    }
    Ok(())
}
fn chain_issues(
    chain: &Peptidoform,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<Vec<ConversionIssue>> {
    let first = conversion::issues_with_session(chain, resolver, false)?;
    add_warnings(warnings, resolver, Vec::new())?;
    if first.is_empty() {
        return Ok(first);
    }
    // Source re-runs issue collection after the separate advisory predicate.
    let second = conversion::issues_with_session(chain, resolver, false)?;
    add_warnings(warnings, resolver, Vec::new())?;
    Ok(second)
}
fn issue(text: &str, budget: &mut Budget) -> Result<Vec<ConversionIssue>> {
    budget.items(1)?;
    budget.text(text)?;
    budget.allocate(text.len().saturating_add(size_of::<ConversionIssue>()))?;
    Ok(vec![ConversionIssue {
        issue_type: ConversionIssueType::UnsupportedFeature,
        description: text.into(),
        position: Some(0),
    }])
}
struct Link<'a> {
    position: usize,
    mass: f64,
    id: &'a str,
}
fn first_link<'a>(chain: &'a Peptidoform, budget: &mut Budget) -> Result<Option<Link<'a>>> {
    budget.items(chain.sequence.len())?;
    let mut position = 0usize;
    for section in &chain.sequence {
        if let SequenceSection::Element(element) = section {
            budget.items(element.modifications.len())?;
            for modification in &element.modifications {
                if let Some((tag, Some(label))) = modification.alternatives.first() {
                    if label.label_type == LabelType::Crosslink {
                        budget.text(&label.identifier)?;
                        let mass = if let ModificationTag::MassDelta(delta) = tag {
                            delta.mass
                        } else {
                            modification
                                .resolved_mod
                                .as_ref()
                                .map_or(0., |m| m.diff_mono_mass())
                        };
                        return Ok(Some(Link {
                            position,
                            mass,
                            id: &label.identifier,
                        }));
                    }
                }
            }
            // CPP-047: source excludes range/ambiguous-region letters here.
            position = position
                .checked_add(1)
                .ok_or_else(|| invalid("ProForma crosslink position overflows"))?;
        }
    }
    Ok(None)
}
fn collect(
    input: Input<'_>,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<Vec<ConversionIssue>> {
    match input {
        Input::Chain(chain) => chain_issues(chain, resolver, warnings),
        Input::Ion(ion) => {
            resolver.budget.consume(1)?;
            if ion.chains.is_empty() {
                return issue("No peptide chains to fragment", resolver.budget);
            }
            if ion.is_chimeric {
                return issue(
                    "Theoretical spectrum generation not supported for chimeric spectra.",
                    resolver.budget,
                );
            }
            if ion.chains.len() == 1 {
                return chain_issues(&ion.chains[0], resolver, warnings);
            }
            if ion.chains.len() != 2 {
                return issue(
                    "Only two-chain cross-links are currently supported for spectrum generation",
                    resolver.budget,
                );
            }
            let a = first_link(&ion.chains[0], resolver.budget)?;
            let b = first_link(&ion.chains[1], resolver.budget)?;
            let (Some(a), Some(b)) = (a, b) else {
                return issue("Cross-link label not found in both chains", resolver.budget);
            };
            resolver
                .budget
                .consume(a.id.len().max(b.id.len()).saturating_add(1))?;
            if a.id != b.id {
                return issue(
                    "Cross-link labels don't match between chains",
                    resolver.budget,
                );
            }
            Ok(Vec::new())
        }
    }
}
fn reject(issues: &[ConversionIssue], budget: &mut Budget) -> Result<()> {
    if issues.is_empty() {
        return Ok(());
    }
    budget.items(issues.len())?;
    let length = issues.iter().fold(28usize, |n, i| {
        n.saturating_add(i.description.len()).saturating_add(2)
    });
    if length > MAX_PROFORMA_CONVERSION_TEXT_BYTES {
        return Err(invalid("ProForma spectrum diagnostic text limit exceeded"));
    }
    budget.consume(length)?;
    budget.allocate(length)?;
    let mut text = String::new();
    text.try_reserve_exact(length)
        .map_err(|_| invalid("ProForma spectrum allocation failed"))?;
    text.push_str("Spectrum generation failed: ");
    for issue in issues {
        text.push_str(&issue.description);
        text.push_str("; ")
    }
    Err(invalid(text))
}
fn flags(options: &SpectrumGenerationOptions, budget: &mut Budget) -> Result<[bool; 8]> {
    budget.text(&options.ion_types)?;
    budget.consume(options.ion_types.len().saturating_mul(8))?;
    Ok(['a', 'b', 'c', 'x', 'y', 'z', 'M', 'I'].map(|c| options.ion_types.contains(c)))
}
fn convert(
    chain: &Peptidoform,
    policy: ConversionPolicy,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<crate::chemistry::AASequence> {
    let (sequence, extra) = conversion::convert_with_session(chain, policy, resolver, false)?;
    add_warnings(warnings, resolver, extra)?;
    Ok(sequence)
}
fn generate_chain(
    chain: &Peptidoform,
    options: &SpectrumGenerationOptions,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<MSSpectrum> {
    let issues = chain_issues(chain, resolver, warnings)?;
    reject(&issues, resolver.budget)?;
    generate_chain_after_check(chain, options, resolver, warnings)
}
fn generate_chain_after_check(
    chain: &Peptidoform,
    options: &SpectrumGenerationOptions,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<MSSpectrum> {
    let sequence = convert(chain, ConversionPolicy::FailOnLoss, resolver, warnings)?;
    let f = flags(options, resolver.budget)?;
    resolver.budget.allocate(
        size_of::<TheoreticalSpectrumGenerator>()
            .saturating_add(8 * size_of::<TheoreticalIonSeries>()),
    )?;
    let generator = TheoreticalSpectrumGenerator {
        ion_series: [
            TheoreticalIonSeries::A,
            TheoreticalIonSeries::B,
            TheoreticalIonSeries::C,
            TheoreticalIonSeries::X,
            TheoreticalIonSeries::Y,
            TheoreticalIonSeries::Z,
        ]
        .into_iter()
        .zip(f)
        .filter_map(|(series, yes)| yes.then_some(series))
        .collect(),
        add_precursor_peaks: f[6],
        add_abundant_immonium_ions: f[7],
        add_losses: options.add_losses,
        add_metainfo: options.add_metainfo,
        ..Default::default()
    };
    let min = u8::try_from(options.min_charge)
        .map_err(|_| invalid("ordinary ProForma fragment charges require 1..=255"))?;
    let max = u8::try_from(options.max_charge)
        .map_err(|_| invalid("ordinary ProForma fragment charges require 1..=255"))?;
    generator.generate_for_proforma(
        &sequence,
        min,
        max,
        &mut resolver.budget.work,
        &mut resolver.budget.bytes,
    )
}
fn generate_ion(
    ion: &PeptidoformIon,
    options: &SpectrumGenerationOptions,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<MSSpectrum> {
    if ion.chains.len() == 1 {
        return generate_chain(&ion.chains[0], options, resolver, warnings);
    }
    let a = first_link(&ion.chains[0], resolver.budget)?.expect("checked first link");
    let b = first_link(&ion.chains[1], resolver.budget)?.expect("checked first link");
    let alpha = convert(
        &ion.chains[0],
        ConversionPolicy::BestEffort,
        resolver,
        warnings,
    )?;
    let beta = convert(
        &ion.chains[1],
        ConversionPolicy::BestEffort,
        resolver,
        warnings,
    )?;
    // CPP-038: attach converted endpoint chemistry and add the separate source
    // linker again. Do not compensate for this finite source double counting.
    let mut link = ProteinProteinCrossLink::new(if a.mass > 0.001 { a.mass } else { b.mass })?;
    resolver
        .budget
        .allocate(2 * size_of::<crate::chemistry::AASequence>() + 128)?;
    link.alpha = Some(Arc::new(alpha));
    link.beta = Some(Arc::new(beta));
    link.cross_link_position = (
        isize::try_from(a.position).map_err(|_| invalid("crosslink index overflows"))?,
        isize::try_from(b.position).map_err(|_| invalid("crosslink index overflows"))?,
    );
    resolver.budget.text(a.id)?;
    resolver.budget.allocate(a.id.len())?;
    link.cross_linker_name = a.id.into();
    let f = flags(options, resolver.budget)?;
    let mut generator = TheoreticalSpectrumGeneratorXLMS::default();
    generator.options.add_a_ions = f[0];
    generator.options.add_b_ions = f[1];
    generator.options.add_c_ions = f[2];
    generator.options.add_x_ions = f[3];
    generator.options.add_y_ions = f[4];
    generator.options.add_z_ions = f[5];
    generator.options.add_precursor_peaks = f[6];
    generator.options.add_losses = options.add_losses;
    generator.options.add_metainfo = options.add_metainfo;
    let mut spectrum = MSSpectrum::default();
    for alpha in [true, false] {
        generator.crosslink_with_budget(
            &mut spectrum,
            &link,
            alpha,
            options.min_charge,
            options.max_charge,
            &mut resolver.budget.work,
            &mut resolver.budget.bytes,
        )?;
    }
    let n = spectrum.len();
    resolver.budget.items(n)?;
    let levels = (usize::BITS - n.max(1).leading_zeros()) as usize;
    let names = spectrum.string_data_arrays.first().map_or(0, |array| {
        array
            .data
            .iter()
            .fold(0usize, |n, s| n.saturating_add(s.len()))
    });
    resolver.budget.consume(
        n.saturating_mul(levels.saturating_add(16))
            .saturating_add(names),
    )?;
    resolver.budget.allocate(
        n.saturating_mul((size_of::<Peak1D>() + size_of::<String>() + 3 * size_of::<usize>()) * 4)
            .saturating_add(names.saturating_mul(2))
            .saturating_add(4096),
    )?;
    spectrum.sort_by_position()?;
    Ok(spectrum)
}
fn run(
    input: Input<'_>,
    options: Option<&SpectrumGenerationOptions>,
    registry: &mut ModificationsDB,
    budget: &mut Budget,
) -> Result<ConversionEvaluation<Output>> {
    let mut resolver = Resolver::new(registry, budget);
    let mut warnings = Vec::new();
    let issues = collect(input, &mut resolver, &mut warnings)?;
    let spectrum = if let Some(options) = options {
        reject(&issues, resolver.budget)?;
        match input {
            Input::Chain(chain) => {
                generate_chain_after_check(chain, options, &mut resolver, &mut warnings)?
            }
            Input::Ion(ion) => generate_ion(ion, options, &mut resolver, &mut warnings)?,
        }
    } else {
        MSSpectrum::default()
    };
    let (result, staged) = conversion::finish(Output { issues, spectrum }, resolver, warnings)?;
    if let Some(next) = staged {
        *registry = next;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn db() -> ModificationsDB {
        ModificationsDB::from_records(vec![]).unwrap()
    }
    #[test]
    fn late_shared_counter_failures_never_publish_formula_definitions() {
        let chain = Peptidoform::parse("AM[Formula:O]A").unwrap();
        let ion = PeptidoformIon::parse("AM[Formula:O#XL1]A//AM[#XL1]A").unwrap();
        let options = SpectrumGenerationOptions::default();
        for input in [Input::Chain(&chain), Input::Ion(&ion)] {
            let mut budget = Budget::default();
            let mut registry = db();
            let good = run(input, Some(&options), &mut registry, &mut budget).unwrap();
            assert!(!good.value.spectrum.is_empty());
            assert_eq!(registry.len(), 1);
            let used = [
                MAX_PROFORMA_RESOLUTION_WORK - budget.work,
                MAX_PROFORMA_RESOLUTION_BYTES - budget.bytes,
                MAX_PROFORMA_RESOLUTION_ITEMS - budget.items,
            ];
            for (index, consumed) in used.into_iter().enumerate() {
                let mut budget = Budget::default();
                match index {
                    0 => budget.work = consumed - 1,
                    1 => budget.bytes = consumed - 1,
                    _ => budget.items = consumed - 1,
                }
                let before = (budget.work, budget.bytes, budget.items);
                let mut registry = db();
                assert!(run(input, Some(&options), &mut registry, &mut budget).is_err());
                assert!(registry.is_empty());
                let after = (budget.work, budget.bytes, budget.items);
                assert!(after.0 < before.0 && after.1 < before.1 && after.2 < before.2);
            }
        }
        assert_eq!(chain, Peptidoform::parse("AM[Formula:O]A").unwrap());
        assert_eq!(
            ion,
            PeptidoformIon::parse("AM[Formula:O#XL1]A//AM[#XL1]A").unwrap()
        );
    }
}
