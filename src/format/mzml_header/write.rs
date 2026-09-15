// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
use super::{xml::Xml, *};
use std::hash::{Hash, Hasher};

pub(super) const EMPTY_HISTORY: &str = "openms-rust:empty-processing-history";
pub(super) const EMPTY_ACTIONS: &str = "openms-rust:empty-processing-actions";
const MARKER_VALUE: &str = "1";
/// Name of the software that the processing method of an empty history refers
/// to. It is written, as [`PLACEHOLDER_SOFTWARE`], only when some record has
/// no processing history; `so_default` itself is the source's empty
/// `Software()`.
pub(super) const FALLBACK_SOFTWARE: &str = "OpenMS Rust mandatory mzML processing placeholder";
/// Identifier of the empty-history placeholder software.
const PLACEHOLDER_SOFTWARE: &str = "so_default_empty_history";
/// Source identifier of the default instrument's software (`so_in_0`).
const INSTRUMENT_SOFTWARE: &str = "so_in_0";
/// Source identifier of the run (`ru_0`).
const RUN_ID: &str = "ru_0";

pub(crate) struct ArrayHeader {
    pub attrs: String,
    pub params: String,
}
pub(crate) struct Plan {
    pub prefix: String,
    pub spectra: Vec<String>,
    pub chromatograms: Vec<String>,
    pub arrays: Vec<ArrayHeader>,
    pub spectrum_metadata: Vec<String>,
    pub chromatogram_metadata: Vec<String>,
}
struct RecordRef {
    source: Option<usize>,
    processing: usize,
}
struct Build<'a> {
    sources: Vec<&'a SourceFile>,
    histories: Vec<&'a [Arc<DataProcessing>]>,
    history_hashes: BTreeMap<u64, Vec<usize>>,
    records: Vec<RecordRef>,
    array_refs: Vec<(Option<usize>, &'a MetaInfo)>,
}
impl<'a> Build<'a> {
    fn history(&mut self, history: &'a [Arc<DataProcessing>], work: &mut Work) -> Result<usize> {
        measure_history(history, work)?;
        let hash = history_hash(history);
        if let Some(indices) = self.history_hashes.get(&hash) {
            for &index in indices {
                measure_history(history, work)?;
                if self.histories[index] == history {
                    return Ok(index);
                }
            }
        }
        work.slots::<&[Arc<DataProcessing>]>(4)?;
        work.meter().tree::<(u64, Vec<usize>)>(1)?;
        work.slots::<usize>(4)?;
        let index = self.histories.len();
        self.histories.push(history);
        self.history_hashes.entry(hash).or_default().push(index);
        Ok(index)
    }
    fn record(
        &mut self,
        source: &'a SourceFile,
        history: &'a [Arc<DataProcessing>],
        work: &mut Work,
    ) -> Result<()> {
        let processing = self.history(history, work)?;
        let source = if source_nondefault(source) {
            crate::kernel::acquisition_fields::source(&mut work.meter(), source)?;
            work.slots::<&SourceFile>(4)?;
            let index = self.sources.len();
            self.sources.push(source);
            Some(index)
        } else {
            None
        };
        work.slots::<RecordRef>(4)?;
        self.records.push(RecordRef { source, processing });
        Ok(())
    }
    fn arrays<T>(&mut self, arrays: &'a [DataArray<T>], work: &mut Work) -> Result<()> {
        work.slots::<(Option<usize>, &MetaInfo)>(arrays.len().saturating_mul(4))?;
        for a in arrays {
            a.description_with_budget(&mut work.remaining, &mut work.bytes)?;
            validate_scalar_metadata(&a.metadata)?;
            let history = if a.data_processing.is_empty() {
                None
            } else {
                Some(self.history(&a.data_processing, work)?)
            };
            self.array_refs.push((history, &a.metadata));
        }
        Ok(())
    }
}
/// Plan the header and per-record reference text for `experiment`.
///
/// The allowance is the reader's fixed header allowance ([`Work::default`])
/// multiplied by `writer_shares`: once for the experiment-level header and
/// once for every spectrum and chromatogram, whose metadata, processing
/// history, array descriptions and controlled-vocabulary lookups are charged
/// against the same pool. A single fixed allowance refused realistic runs after
/// about 650 records.
pub(crate) fn prepare(experiment: &MSExperiment) -> Result<Plan> {
    let shares = writer_shares(experiment);
    prepare_with_work(
        experiment,
        Work {
            remaining: MAX_WORK.saturating_mul(shares),
            bytes: MAX_BYTES.saturating_mul(shares),
        },
    )
}
fn prepare_with_work(experiment: &MSExperiment, mut work: Work) -> Result<Plan> {
    guard(experiment)?;
    let settings = &experiment.settings;
    settings.with_budget(&mut work.remaining, &mut work.bytes)?;
    work.slots::<MSSpectrum>(experiment.spectra.len())?;
    work.slots::<MSChromatogram>(experiment.chromatograms.len())?;
    let mut build = Build {
        sources: Vec::new(),
        histories: Vec::new(),
        history_hashes: BTreeMap::new(),
        records: Vec::new(),
        array_refs: Vec::new(),
    };
    work.slots::<&SourceFile>(settings.source_files.len().saturating_mul(4))?;
    build.sources.extend(&settings.source_files);
    for spectrum in &experiment.spectra {
        work.meter().meta(&spectrum.metadata)?;
        record_transport::validate(&spectrum.metadata, false)?;
        build.record(&spectrum.source_file, &spectrum.data_processing, &mut work)?;
        build.arrays(&spectrum.float_data_arrays, &mut work)?;
        build.arrays(&spectrum.integer_data_arrays, &mut work)?;
        build.arrays(&spectrum.string_data_arrays, &mut work)?;
        work.slots::<Acquisition>(spectrum.acquisition_info.acquisitions.len())?;
        for scan in &spectrum.acquisition_info.acquisitions {
            work.meter().meta(&scan.metadata)?;
            if let Some(value) = scan.metadata.get("instrument_configuration_ref") {
                if value.unit().is_some()
                    || !settings
                        .instrument_configurations
                        .contains_key(value.as_str()?)
                {
                    return Err(invalid(
                        "scan references an absent additional instrument configuration",
                    ));
                }
            }
        }
    }
    for c in &experiment.chromatograms {
        work.meter().meta(&c.metadata)?;
        record_transport::validate(&c.metadata, true)?;
        if source_nondefault(&c.source_file) {
            return Err(Error::Unsupported(
                "chromatogram source-file has no mzML schema representation".into(),
            ));
        }
        build.record(&c.source_file, &c.data_processing, &mut work)?;
        build.arrays(&c.float_data_arrays, &mut work)?;
        build.arrays(&c.integer_data_arrays, &mut work)?;
        build.arrays(&c.string_data_arrays, &mut work)?;
    }
    if build.histories.is_empty() {
        build.history(&[], &mut work)?;
    }
    let mut x = Xml::new(&mut work);
    let mut default_instrument = "ic_0".to_owned();
    while settings
        .instrument_configurations
        .contains_key(&default_instrument)
    {
        x.work.charge(
            default_instrument.len().saturating_mul(64),
            default_instrument.len().saturating_mul(2),
        )?;
        default_instrument.push('_');
    }
    // Source root (`MzMLHandler.cpp:4851`): schema location, the document
    // accession (written even when empty) and the version. A read `mzML@id`
    // lives in the run metadata as `mzml_id` and is written back as the run
    // userParam the source writes, not as a root attribute.
    x.raw("<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\"")?;
    x.attribute("accession", &settings.document.identifier)?;
    if let Some(value) = settings.metadata.get("mzml_id") {
        if value.unit().is_some() {
            return Err(invalid("mzml_id cannot have a unit"));
        }
        value.as_str()?;
    }
    x.raw(" version=\"1.1.0\">\n<cvList count=\"5\">\n")?;
    // The source's fixed controlled-vocabulary list (`MzMLHandler.cpp:4855-4861`).
    for (id, name, version, uri) in [
        (
            "MS",
            "Proteomics Standards Initiative Mass Spectrometry Ontology",
            None,
            "http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo",
        ),
        (
            "UO",
            "Unit Ontology",
            None,
            "http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo",
        ),
        (
            "BTO",
            "BrendaTissue545",
            Some("unknown"),
            "http://www.brenda-enzymes.info/ontology/tissue/tree/update/update_files/BrendaTissueOBO",
        ),
        (
            "GO",
            "Gene Ontology - Slim Versions",
            Some("unknown"),
            "http://www.geneontology.org/GO_slims/goslim_goa.obo",
        ),
        (
            "PATO",
            "Quality ontology",
            Some("unknown"),
            "http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/quality.obo",
        ),
    ] {
        x.define_id(id)?;
        x.raw("<cv")?;
        x.attribute("id", id)?;
        x.attribute("fullName", name)?;
        if let Some(version) = version {
            x.attribute("version", version)?;
        }
        x.attribute("URI", uri)?;
        x.raw("/>\n")?;
    }
    x.end("cvList")?;
    x.start("fileDescription", &[])?;
    x.start("fileContent", &[])?;
    // Reuse the already source-audited scan-mode aggregation. Its at-most-15
    // fixed terms fit this fixed precharge; the temporary cannot scale with peaks.
    x.work
        .charge(experiment.spectra.len().saturating_add(8192), 8192)?;
    let mut content = Vec::new();
    settings_metadata::write_file_content(&mut content, experiment)?;
    x.raw(std::str::from_utf8(&content).map_err(|_| invalid("invalid generated file content"))?)?;
    x.end("fileContent")?;
    if !build.sources.is_empty() {
        let count = x.integer(build.sources.len())?;
        x.start("sourceFileList", &[("count", &count)])?;
        for (index, source) in build.sources.iter().enumerate() {
            let id = identifier("sf", index, x.work)?;
            write_source(&mut x, &id, source)?;
        }
        x.end("sourceFileList")?;
    }
    for contact in &settings.contacts {
        write_contact(&mut x, contact)?;
    }
    x.end("fileDescription")?;
    x.start("sampleList", &[("count", "1")])?;
    write_sample(&mut x, &settings.sample)?;
    x.end("sampleList")?;
    let method_count = build
        .histories
        .iter()
        .try_fold(0usize, |n, h| n.checked_add(h.len()).ok_or_else(resource))?;
    // An empty history has no software of its own; its processing method
    // refers to one extra placeholder entry, which the reader recognizes.
    let placeholder = build.histories.iter().any(|history| history.is_empty());
    let software_count = method_count
        .checked_add(settings.instrument_configurations.len())
        .and_then(|n| n.checked_add(2 + usize::from(placeholder)))
        .ok_or_else(resource)?;
    let count = x.integer(software_count)?;
    x.start("softwareList", &[("count", &count)])?;
    x.work.charge(128, 128)?;
    // Source order and identifiers (`MzMLHandler.cpp:5107-5125`): the default
    // instrument's software, one per additional configuration, then the
    // fallback software, which is `Software()` - empty name and version.
    write_software(&mut x, INSTRUMENT_SOFTWARE, &settings.instrument.software)?;
    for (index, value) in settings.instrument_configurations.values().enumerate() {
        let id = configuration_software(index, x.work)?;
        write_software(&mut x, &id, &value.software)?;
    }
    write_software(&mut x, "so_default", &Software::default())?;
    if placeholder {
        x.work.charge(128, 128)?;
        write_software(
            &mut x,
            PLACEHOLDER_SOFTWARE,
            &Software {
                name: FALLBACK_SOFTWARE.into(),
                version: MARKER_VALUE.into(),
                ..Default::default()
            },
        )?;
    }
    for (index, history) in build.histories.iter().enumerate() {
        for (method, value) in history.iter().enumerate() {
            let id = method_id(index, method, x.work)?;
            write_software(&mut x, &id, &value.software)?;
        }
    }
    x.end("softwareList")?;
    let count = x.integer(
        settings
            .instrument_configurations
            .len()
            .checked_add(1)
            .ok_or_else(resource)?,
    )?;
    x.start("instrumentConfigurationList", &[("count", &count)])?;
    write_instrument(
        &mut x,
        &default_instrument,
        INSTRUMENT_SOFTWARE,
        &settings.instrument,
    )?;
    for (index, (id, value)) in settings.instrument_configurations.iter().enumerate() {
        parameter_id(id)?;
        let software = configuration_software(index, x.work)?;
        write_instrument(&mut x, id, &software, value)?;
    }
    x.end("instrumentConfigurationList")?;
    let count = x.integer(build.histories.len())?;
    x.start("dataProcessingList", &[("count", &count)])?;
    for (index, history) in build.histories.iter().enumerate() {
        write_history(&mut x, index, history)?;
    }
    x.end("dataProcessingList")?;
    // The source writes the fixed run identifier `ru_0` (`MzMLHandler.cpp:5211`),
    // whatever the input's run was called; the input slices of the benchmark
    // carry `ru_0` too.
    x.define_id(RUN_ID)?;
    x.raw("<run")?;
    x.attribute("id", RUN_ID)?;
    x.attribute("defaultInstrumentConfigurationRef", &default_instrument)?;
    x.attribute("sampleRef", "sa_0")?;
    if settings.date_time != DateTime::default() {
        if !settings.date_time.is_valid() {
            return Err(invalid("invalid experiment date-time cannot be written"));
        }
        x.work.charge(256, 256)?;
        let time = timestamp_text(settings.date_time)?;
        let preserved = settings
            .metadata
            .get("mzml_start_time_stamp")
            .map(|value| {
                if value.unit().is_some() {
                    return Err(invalid("raw timestamp cannot have a unit"));
                }
                value.as_str()
            })
            .transpose()?;
        let text = preserved
            .filter(|raw| DateTime::parse(raw).is_ok_and(|parsed| parsed == settings.date_time))
            .unwrap_or(&time);
        x.attribute("startTimeStamp", text)?;
    } else if settings.metadata.contains_key("mzml_start_time_stamp") {
        return Err(invalid(
            "raw timestamp metadata requires an initialized date-time",
        ));
    }
    // Source `MzMLHandler.cpp:5219-5222`: the run's first source file is its
    // default source file. `build.sources` starts with the run's own files.
    if !settings.source_files.is_empty() {
        let source = identifier("sf", 0, x.work)?;
        x.attribute("defaultSourceFileRef", &source)?;
    }
    x.raw(">\n")?;
    if !settings.fraction_identifier.is_empty() {
        x.cv_text("MS:1000858", &settings.fraction_identifier)?;
    }
    x.metadata("run", &settings.metadata, &["mzml_start_time_stamp"])?;
    let prefix = std::mem::take(&mut x.text);
    let mut refs = Vec::new();
    x.work.slots::<String>(build.records.len())?;
    refs.try_reserve_exact(build.records.len())
        .map_err(|_| resource())?;
    for record in &build.records {
        // A record whose history is the first one inherits it from the list's
        // `defaultDataProcessingRef`, so the attribute carries nothing; the
        // source omits it for the same reason on every spectrum after the
        // first (`MzMLHandler.cpp:5257-5271`) and on every chromatogram.
        if record.processing != 0 {
            let processing = identifier("dp", record.processing, x.work)?;
            x.attribute("dataProcessingRef", &processing)?;
        }
        if let Some(index) = record.source {
            let source = identifier("sf", index, x.work)?;
            x.attribute("sourceFileRef", &source)?;
        }
        refs.push(std::mem::take(&mut x.text));
    }
    x.work.slots::<String>(experiment.chromatograms.len())?;
    let chromatograms = refs.split_off(experiment.spectra.len());
    let mut arrays = Vec::new();
    x.work.slots::<ArrayHeader>(build.array_refs.len())?;
    arrays
        .try_reserve_exact(build.array_refs.len())
        .map_err(|_| resource())?;
    for (history, metadata) in build.array_refs {
        if let Some(index) = history {
            let id = identifier("dp", index, x.work)?;
            x.attribute("dataProcessingRef", &id)?;
        }
        let attrs = std::mem::take(&mut x.text);
        x.metadata("binaryDataArray", metadata, &[])?;
        arrays.push(ArrayHeader {
            attrs,
            params: std::mem::take(&mut x.text),
        });
    }
    x.work.slots::<String>(
        experiment
            .spectra
            .len()
            .saturating_add(experiment.chromatograms.len()),
    )?;
    let mut spectrum_metadata = Vec::new();
    let mut chromatogram_metadata = Vec::new();
    spectrum_metadata
        .try_reserve_exact(experiment.spectra.len())
        .map_err(|_| resource())?;
    chromatogram_metadata
        .try_reserve_exact(experiment.chromatograms.len())
        .map_err(|_| resource())?;
    for (meta, name, chrom) in experiment
        .spectra
        .iter()
        .map(|s| (&s.metadata, &s.name, false))
        .chain(
            experiment
                .chromatograms
                .iter()
                .map(|c| (&c.metadata, &c.name, true)),
        )
    {
        x.metadata(
            if chrom { "chromatogram" } else { "spectrum" },
            meta,
            if chrom {
                &record_transport::CHROMATOGRAM_SKIP
            } else {
                &record_transport::SPECTRUM_SKIP
            },
        )?;
        if !name.is_empty() {
            let value = MetaValue::from(x.work.copy(name)?);
            x.user(NAME_KEY, &value)?;
        }
        if chrom {
            chromatogram_metadata.push(std::mem::take(&mut x.text));
        } else {
            spectrum_metadata.push(std::mem::take(&mut x.text));
        }
    }
    Ok(Plan {
        spectrum_metadata,
        chromatogram_metadata,
        prefix,
        spectra: refs,
        chromatograms,
        arrays,
    })
}
fn identifier(prefix: &str, index: usize, work: &mut Work) -> Result<String> {
    work.charge(64, 64)?;
    Ok(format!("{prefix}_{index:020}"))
}
/// Source identifier of the software of the `index`-th additional instrument
/// configuration, `so_configuration_<index>` (`MzMLHandler.cpp:5112`).
fn configuration_software(index: usize, work: &mut Work) -> Result<String> {
    work.charge(64, 64)?;
    Ok(format!("so_configuration_{index}"))
}
fn method_id(history: usize, method: usize, work: &mut Work) -> Result<String> {
    work.charge(96, 96)?;
    Ok(format!("so_dp_{history:020}_{method:020}"))
}
fn measure_history(history: &[Arc<DataProcessing>], work: &mut Work) -> Result<()> {
    work.slots::<Arc<DataProcessing>>(history.len())?;
    for value in history {
        if !value.software.cv_terms.is_empty() {
            return Err(Error::Unsupported(
                "arbitrary software CVTermList entries are not represented".into(),
            ));
        }
        work.slots::<DataProcessing>(1)?;
        work.meter().text(&value.software.name)?;
        work.meter().text(&value.software.version)?;
        work.meter().cv(&value.software.cv_terms)?;
        work.meter().tree::<ProcessingAction>(value.actions.len())?;
        work.meter().meta(&value.metadata)?;
        validate_scalar_metadata(&value.metadata)?;
        validate_scalar_metadata(&value.software.cv_terms.metadata)?;
        if value.metadata.contains_key(EMPTY_HISTORY) || value.metadata.contains_key(EMPTY_ACTIONS)
        {
            return Err(invalid("reserved processing marker key"));
        }
        for action in &value.actions {
            if !read::PROCESSING_ACTIONS.iter().any(|r| r.1 == *action) {
                return Err(Error::Unsupported(
                    "processing action has no source mzML mapping".into(),
                ));
            }
        }
    }
    Ok(())
}
fn history_hash(history: &[Arc<DataProcessing>]) -> u64 {
    let mut h = std::collections::hash_map::DefaultHasher::new();
    history.len().hash(&mut h);
    for value in history {
        value.software.name.hash(&mut h);
        value.software.version.hash(&mut h);
        value.software.cv_terms.hash(&mut h);
        value.actions.len().hash(&mut h);
        for action in &value.actions {
            (*action as u32).hash(&mut h);
        }
        value.metadata.hash(&mut h);
        value
            .completion_time
            .map(|t| (t.components(), t.millisecond(), t.is_valid()))
            .hash(&mut h);
    }
    h.finish()
}
pub(super) fn source_nondefault(s: &SourceFile) -> bool {
    !s.name.is_empty()
        || !s.path.is_empty()
        || s.size_mb.to_bits() != 0
        || !s.file_type.is_empty()
        || !s.checksum.is_empty()
        || s.checksum_type != ChecksumType::Unknown
        || !s.native_id_type.is_empty()
        || !s.native_id_type_accession.is_empty()
        || !s.cv_terms.is_empty()
        || !s.cv_terms.metadata.is_empty()
}
pub(crate) fn guard(experiment: &MSExperiment) -> Result<()> {
    let settings = &experiment.settings;
    let s = &settings.sample;
    let h = &settings.hplc;
    if !settings.comment.is_empty()
        || !s.organism.is_empty()
        || !s.subsamples.is_empty()
        || !h.instrument.is_empty()
        || !h.column.is_empty()
        || h.temperature != 21
        || h.pressure != 0
        || h.flux != 0
        || !h.comment.is_empty()
        || !h.gradient.eluents().is_empty()
        || !h.gradient.timepoints().is_empty()
        || !h.gradient.percentages().is_empty()
    {
        return Err(Error::Unsupported(
            "experiment comment/HPLC/sample organism/subsamples have no source mzML representation"
                .into(),
        ));
    }
    if !settings.instrument.software.cv_terms.is_empty() {
        return Err(Error::Unsupported(
            "arbitrary instrument software CV terms are not represented".into(),
        ));
    }
    Ok(())
}
fn write_software(x: &mut Xml<'_>, id: &str, value: &Software) -> Result<()> {
    if !value.cv_terms.is_empty() {
        return Err(Error::Unsupported(
            "arbitrary software CV terms are not represented".into(),
        ));
    }
    validate_scalar_metadata(&value.cv_terms.metadata)?;
    x.start("software", &[("id", id), ("version", &value.version)])?;
    let term = ControlledVocabulary::psi_ms()?.first_child_with_name_with_budget(
        "MS:1000531",
        &value.name,
        &mut x.work.remaining,
        &mut x.work.bytes,
    )?;
    // Exact name transport avoids source's lossy fallback aliases (appending
    // ' software' / 'TOPP '). An unrecognized name remains a custom software.
    if let Some(term) = term.filter(|t| t.id != "MS:1000799") {
        x.cv(&term.id, None)?;
    } else {
        x.cv_text("MS:1000799", &value.name)?;
    }
    x.metadata("software", &value.cv_terms.metadata, &[])?;
    x.end("software")
}
fn write_source(x: &mut Xml<'_>, id: &str, value: &SourceFile) -> Result<()> {
    if value.size_mb.to_bits() != 0 || !value.cv_terms.is_empty() {
        return Err(Error::Unsupported(
            "source-file size/arbitrary CV terms are not represented".into(),
        ));
    }
    value.validate()?;
    validate_scalar_metadata(&value.cv_terms.metadata)?;
    x.start(
        "sourceFile",
        &[("id", id), ("name", &value.name), ("location", &value.path)],
    )?;
    // Unlike source forced fallback, absent metadata stays absent; these CVs
    // are semantic recommendations, not required by XML Schema ParamGroupType.
    match value.checksum_type {
        ChecksumType::Sha1 => x.cv_text("MS:1000569", &value.checksum)?,
        ChecksumType::Md5 => x.cv_text("MS:1000568", &value.checksum)?,
        ChecksumType::Unknown if !value.checksum.is_empty() => {
            return Err(invalid("source checksum requires an algorithm"));
        }
        _ => {}
    }
    for (name, parent) in [
        (&value.file_type, "MS:1000560"),
        (&value.native_id_type, "MS:1000767"),
    ] {
        if name.is_empty() {
            continue;
        }
        let term = ControlledVocabulary::psi_ms()?
            .first_child_with_name_with_budget(
                parent,
                name,
                &mut x.work.remaining,
                &mut x.work.bytes,
            )?
            .ok_or_else(|| invalid("source file type/name has no pinned CV term"))?;
        if parent == "MS:1000767"
            && !value.native_id_type_accession.is_empty()
            && value.native_id_type_accession != term.id
        {
            return Err(invalid("native-ID name/accession mismatch"));
        }
        x.cv(&term.id, None)?;
    }
    if value.native_id_type.is_empty() && !value.native_id_type_accession.is_empty() {
        return Err(invalid("native-ID accession has no corresponding name"));
    }
    x.metadata("sourceFile", &value.cv_terms.metadata, &[])?;
    x.end("sourceFile")
}
fn write_contact(x: &mut Xml<'_>, value: &ContactPerson) -> Result<()> {
    if value.metadata.contains_key("contact_info") {
        return Err(invalid(
            "contact_info metadata collides with dedicated field",
        ));
    }
    // Source comma splitting trims both names and ignores later commas. Reject
    // names that cannot survive that exact source convention.
    if value.first_name.contains(',')
        || value.last_name.contains(',')
        || crate::data_structures::list::trim(&value.first_name) != value.first_name
        || crate::data_structures::list::trim(&value.last_name) != value.last_name
    {
        return Err(Error::Unsupported(
            "contact names cannot roundtrip source comma splitting".into(),
        ));
    }
    x.work.charge(
        value.first_name.len().saturating_add(value.last_name.len()),
        value
            .first_name
            .len()
            .saturating_add(value.last_name.len())
            .saturating_add(2),
    )?;
    let name = format!("{}, {}", value.last_name, value.first_name);
    x.start("contact", &[])?;
    x.cv_text("MS:1000586", &name)?;
    x.cv_text("MS:1000590", &value.institution)?;
    for (id, text) in [
        ("MS:1000587", &value.address),
        ("MS:1000588", &value.url),
        ("MS:1000589", &value.email),
    ] {
        if !text.is_empty() {
            x.cv_text(id, text)?;
        }
    }
    // Dedicated contact_info follows promoted CVs, avoiding source's invalid
    // CV-after-user ordering when ordinary metadata is promotable.
    x.metadata("contact", &value.metadata, &[])?;
    if !value.contact_info.is_empty() {
        let text = x.work.copy(&value.contact_info)?;
        x.user("contact_info", &text.into())?;
    }
    x.end("contact")
}
fn write_sample(x: &mut Xml<'_>, value: &Sample) -> Result<()> {
    if value.metadata.contains_key("comment") {
        return Err(invalid(
            "sample comment metadata collides with dedicated field",
        ));
    }
    x.start("sample", &[("id", "sa_0"), ("name", &value.name)])?;
    if !value.number.is_empty() {
        x.cv_text("MS:1000001", &value.number)?;
    }
    x.cv_float("MS:1000004", value.mass, Some("UO:0000021"))?;
    x.cv_float("MS:1000005", value.volume, Some("UO:0000098"))?;
    x.cv_float("MS:1000006", value.concentration, Some("UO:0000175"))?;
    let state = match value.state {
        SampleState::Unknown => None,
        SampleState::Emulsion => Some("MS:1000047"),
        SampleState::Gas => Some("MS:1000048"),
        SampleState::Liquid => Some("MS:1000049"),
        SampleState::Solid => Some("MS:1000050"),
        SampleState::Solution => Some("MS:1000051"),
        SampleState::Suspension => Some("MS:1000052"),
    };
    if let Some(id) = state {
        x.cv(id, None)?;
    }
    x.metadata("sample", &value.metadata, &[])?;
    if !value.comment.is_empty() {
        let text = x.work.copy(&value.comment)?;
        x.user("comment", &text.into())?;
    }
    x.end("sample")
}
fn write_instrument(x: &mut Xml<'_>, id: &str, software: &str, value: &Instrument) -> Result<()> {
    instrument::validate(value)?;
    x.start("instrumentConfiguration", &[("id", id)])?;
    if value.name.is_empty() {
        x.cv("MS:1000031", None)?;
    } else {
        let term = ControlledVocabulary::psi_ms()?
            .first_child_with_name_with_budget(
                "MS:1000031",
                &value.name,
                &mut x.work.remaining,
                &mut x.work.bytes,
            )?
            .ok_or_else(|| invalid("instrument name has no pinned model CV term"))?;
        x.cv(&term.id, None)?;
    }
    if !value.customizations.is_empty() {
        x.cv_text("MS:1000032", &value.customizations)?;
    }
    if let Some(id) = instrument_terms::write_ion_optics(value.ion_optics)? {
        x.cv(id, None)?;
    }
    x.metadata("instrumentConfiguration", &value.metadata, &[])?;
    let count = value
        .ion_sources
        .len()
        .checked_add(value.mass_analyzers.len())
        .and_then(|n| n.checked_add(value.ion_detectors.len()))
        .ok_or_else(resource)?;
    if count != 0 {
        let count = x.integer(count)?;
        x.start("componentList", &[("count", &count)])?;
        for value in &value.ion_sources {
            let order = x.signed(value.order)?;
            x.start("source", &[("order", &order)])?;
            if let Some(id) = instrument_terms::write_inlet(value.inlet_type)? {
                x.cv(id, None)?;
            }
            let override_id = value
                .metadata
                .get("ionization accession")
                .map(|v| {
                    if v.unit().is_some() {
                        return Err(invalid("instrument accession override cannot have a unit"));
                    }
                    v.as_str()
                })
                .transpose()?;
            let id = override_id.or(instrument_terms::write_ionization(value.ionization_method)?);
            if let Some(id) = id {
                if override_id.is_some()
                    && (!x.work.child(id, "MS:1000008")?
                        || instrument_terms::ionization(id).is_some())
                {
                    return Err(invalid(
                        "ionization override does not retain unknown-method metadata",
                    ));
                }
                x.cv(id, None)?;
            }
            x.metadata("source", &value.metadata, &["ionization accession"])?;
            x.end("source")?;
        }
        for value in &value.mass_analyzers {
            let order = x.signed(value.order)?;
            x.start("analyzer", &[("order", &order)])?;
            let override_id = value
                .metadata
                .get("mass analyzer accession")
                .map(|v| {
                    if v.unit().is_some() {
                        return Err(invalid("instrument accession override cannot have a unit"));
                    }
                    v.as_str()
                })
                .transpose()?;
            let id = override_id.or(instrument_terms::write_analyzer(value.analyzer_type)?);
            if let Some(id) = id {
                if override_id.is_some()
                    && (!x.work.child(id, "MS:1000443")?
                        || instrument_terms::analyzer(id).is_some())
                {
                    return Err(invalid(
                        "analyzer override does not retain unknown-type metadata",
                    ));
                }
                x.cv(id, None)?;
            }
            if let Some(id) = instrument_terms::write_reflectron(value.reflectron_state)? {
                x.cv(id, None)?;
            }
            x.cv_float("MS:1000014", value.accuracy, Some("UO:0000169"))?;
            x.cv_float(
                "MS:1000022",
                value.tof_total_path_length,
                Some("UO:0000008"),
            )?;
            x.cv(
                "MS:1000024",
                Some(&MetaValue::from(i64::from(value.final_ms_exponent))),
            )?;
            x.cv_float(
                "MS:1000025",
                value.magnetic_field_strength,
                Some("UO:0000228"),
            )?;
            x.metadata("analyzer", &value.metadata, &["mass analyzer accession"])?;
            x.end("analyzer")?;
        }
        for value in &value.ion_detectors {
            let order = x.signed(value.order)?;
            x.start("detector", &[("order", &order)])?;
            if let Some(id) = instrument_terms::write_detector(value.detector_type)? {
                x.cv(id, None)?;
            }
            if let Some(id) = instrument_terms::write_acquisition(value.acquisition_mode)? {
                x.cv(id, None)?;
            }
            x.cv_float("MS:1000028", value.resolution, None)?;
            x.cv_float(
                "MS:1000029",
                value.adc_sampling_frequency,
                Some("UO:0000106"),
            )?;
            x.metadata("detector", &value.metadata, &[])?;
            x.end("detector")?;
        }
        x.end("componentList")?;
    }
    x.raw("<softwareRef")?;
    x.attribute("ref", software)?;
    x.raw("/>\n")?;
    x.end("instrumentConfiguration")
}
fn write_history(x: &mut Xml<'_>, index: usize, history: &[Arc<DataProcessing>]) -> Result<()> {
    let id = identifier("dp", index, x.work)?;
    x.start("dataProcessing", &[("id", &id)])?;
    if history.is_empty() {
        x.start(
            "processingMethod",
            &[("order", "0"), ("softwareRef", PLACEHOLDER_SOFTWARE)],
        )?;
        x.cv("MS:1000544", None)?;
        x.user(EMPTY_HISTORY, &MARKER_VALUE.into())?;
        x.end("processingMethod")?;
    }
    for (method, value) in history.iter().enumerate() {
        let software = method_id(index, method, x.work)?;
        // The source writes `order="0"` for every method (`MzMLHandler.cpp:3852`)
        // and both readers keep document order, so the value carries nothing.
        x.start(
            "processingMethod",
            &[("order", "0"), ("softwareRef", &software)],
        )?;
        if value.actions.is_empty() {
            x.cv("MS:1000543", None)?;
        } else {
            for &(id, action, _) in read::PROCESSING_ACTIONS {
                if value.actions.contains(&action) {
                    x.cv(id, None)?;
                }
            }
        }
        if let Some(time) = value.completion_time {
            x.work.charge(128, 128)?;
            x.cv_text("MS:1000747", &timestamp_text(time)?)?;
        }
        x.metadata("processingMethod", &value.metadata, &[])?;
        if value.actions.is_empty() {
            x.user(EMPTY_ACTIONS, &MARKER_VALUE.into())?;
        }
        x.end("processingMethod")?;
    }
    x.end("dataProcessing")
}

pub(super) fn normalize_methods(methods: &mut Vec<Arc<DataProcessing>>) -> Result<()> {
    if methods
        .iter()
        .any(|method| method.metadata.contains_key(EMPTY_HISTORY))
    {
        let exact = methods.len() == 1 && {
            let value = &methods[0];
            value.software.name == FALLBACK_SOFTWARE
                && value.software.version == MARKER_VALUE
                && value.software.cv_terms.is_empty()
                && value.software.cv_terms.metadata.is_empty()
                && value.actions.len() == 1
                && value.actions.contains(&ProcessingAction::ConversionMzML)
                && value.completion_time.is_none()
                && value.metadata.len() == 1
                && marker(value.metadata.get(EMPTY_HISTORY))
        };
        if !exact {
            return Err(invalid(
                "empty-history marker does not match the exact generated placeholder",
            ));
        }
        methods.clear();
        return Ok(());
    }
    for method in methods {
        if method.metadata.contains_key(EMPTY_ACTIONS) {
            if !marker(method.metadata.get(EMPTY_ACTIONS))
                || method.actions.len() != 1
                || !method.actions.contains(&ProcessingAction::DataProcessing)
            {
                return Err(invalid(
                    "empty-action marker does not match its exact generic-action payload",
                ));
            }
            let value = Arc::get_mut(method).ok_or_else(|| {
                invalid("processing marker normalized after references were attached")
            })?;
            value.metadata.remove(EMPTY_ACTIONS);
            value.actions.clear();
        }
    }
    Ok(())
}
fn marker(value: Option<&MetaValue>) -> bool {
    value.is_some_and(|value| {
        value.unit().is_none() && matches!(value.data(),MetaValueData::String(s) if s==MARKER_VALUE)
    })
}

fn timestamp_text(time: DateTime) -> Result<String> {
    if !time.is_valid() || !(1..=9999).contains(&time.date_components().2) {
        return Err(invalid(
            "timestamp is invalid or outside supported XML years",
        ));
    }
    let text = time.format(if time.millisecond() == 0 {
        "yyyy-MM-ddThh:mm:ss"
    } else {
        "yyyy-MM-ddThh:mm:ss.zzz"
    })?;
    if DateTime::parse(&text)? != time {
        return Err(invalid("timestamp has inconsistent partial fields"));
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn aggregate_preflight_precedes_scanning_or_copying_owned_metadata() {
        let mut e = MSExperiment::default();
        e.settings
            .metadata
            .insert("large".into(), "x".repeat(4096).into());
        e.settings.instrument.name = "not a source model".into();
        let error = prepare_with_work(
            &e,
            Work {
                remaining: 0,
                bytes: usize::MAX,
            },
        )
        .err()
        .unwrap();
        assert!(error.to_string().contains("limit") || error.to_string().contains("allowance"));
        let error = prepare_with_work(
            &e,
            Work {
                remaining: usize::MAX,
                bytes: 1,
            },
        )
        .err()
        .unwrap();
        assert!(error.to_string().contains("limit") || error.to_string().contains("allowance"));
        assert_eq!(e.settings.metadata["large"].as_str().unwrap().len(), 4096);
    }
}
