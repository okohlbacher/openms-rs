// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! mzML writer scale and C++ output parity (lane `mzml-writer-scale-parity`).
//!
//! The synthetic experiments are built from one MS1 and one MS2 record in the
//! layout of the benchmark input `sub_centroid_uk222_picked_first600.mzML` (a
//! Q Exactive run peak-picked by OpenMS, cut from PXD-staged `UK222_picked`):
//! the same header, processing history, record metadata, scan and precursor
//! terms, with three peaks per spectrum. Before the per-record writer budgets,
//! every mzML writer entry point refused such an experiment after about 650
//! records with "controlled vocabulary resource limit exceeded"; the smoke
//! benchmark of 2026-09-14 measured 647 passing and 648 failing on the real
//! `UK222_picked` prefixes.
//!
//! The `#[ignore]` tests at the end read the staged benchmark inputs under
//! `/ceph/ibmi/abi/oliver/bench/openms4/inputs` and only run on the IBMI HPC
//! nodes that mount that path: `cargo test --release --test mzml_writer_scale
//! -- --ignored`.
#![cfg(feature = "mzml")]

use openms::format::mzml;
use openms::format::peak_options::PeakFileOptions;
use openms::{MSExperiment, MSSpectrum};
use std::io::{self, Cursor};

/// One MS1 and one MS2 spectrum with the benchmark input's header and
/// metadata; the binary arrays hold three peaks each.
const TEMPLATE: &str = r#"<?xml version="1.0" encoding="ISO-8859-1"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd" accession="" version="1.1.0">
	<cvList count="5">
		<cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>
		<cv id="UO" fullName="Unit Ontology" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>
		<cv id="BTO" fullName="BrendaTissue545" version="unknown" URI="http://www.brenda-enzymes.info/ontology/tissue/tree/update/update_files/BrendaTissueOBO"/>
		<cv id="GO" fullName="Gene Ontology - Slim Versions" version="unknown" URI="http://www.geneontology.org/GO_slims/goslim_goa.obo"/>
		<cv id="PATO" fullName="Quality ontology" version="unknown" URI="http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/quality.obo"/>
	</cvList>
	<fileDescription>
		<fileContent>
			<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" />
			<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" />
		</fileContent>
		<sourceFileList count="1">
			<sourceFile id="sf_ru_0" name="UK222.raw" location="file://">
				<cvParam cvRef="MS" accession="MS:1000569" name="SHA-1" value="781876ca750eb41bde61816d8a8765c730b3cc29" />
				<cvParam cvRef="MS" accession="MS:1000563" name="Thermo RAW format" />
				<cvParam cvRef="MS" accession="MS:1000768" name="Thermo nativeID format" />
			</sourceFile>
		</sourceFileList>
	</fileDescription>
	<sampleList count="1">
		<sample id="sa_0" name="">
			<cvParam cvRef="MS" accession="MS:1000004" name="sample mass" value="0" unitAccession="UO:0000021" unitName="gram" unitCvRef="UO" />
			<cvParam cvRef="MS" accession="MS:1000005" name="sample volume" value="0" unitAccession="UO:0000098" unitName="milliliter" unitCvRef="UO" />
			<cvParam cvRef="MS" accession="MS:1000006" name="sample concentration" value="0" unitAccession="UO:0000175" unitName="gram per liter" unitCvRef="UO" />
		</sample>
	</sampleList>
	<softwareList count="4">
		<software id="so_in_0" version="2.8-280502/2.8.1.2806" >
			<cvParam cvRef="MS" accession="MS:1000532" name="Xcalibur" />
		</software>
		<software id="so_default" version="" >
			<cvParam cvRef="MS" accession="MS:1000799" name="custom unreleased software tool" value="" />
		</software>
		<software id="so_dp_sp_0_pm_0" version="3.0.9987" >
			<cvParam cvRef="MS" accession="MS:1000615" name="ProteoWizard software" />
		</software>
		<software id="so_dp_sp_0_pm_1" version="2.4.0" >
			<cvParam cvRef="MS" accession="MS:1002135" name="TOPP PeakPickerHiRes" />
		</software>
	</softwareList>
	<instrumentConfigurationList count="1">
		<instrumentConfiguration id="ic_0">
			<cvParam cvRef="MS" accession="MS:1001911" name="Q Exactive" />
			<cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="Exactive Series slot #0012"/>
			<componentList count="3">
				<source order="1">
					<cvParam cvRef="MS" accession="MS:1000485" name="nanospray inlet" />
					<cvParam cvRef="MS" accession="MS:1000398" name="nanoelectrospray" />
				</source>
				<analyzer order="2">
					<cvParam cvRef="MS" accession="MS:1000484" name="orbitrap" />
				</analyzer>
				<detector order="3">
					<cvParam cvRef="MS" accession="MS:1000624" name="inductive detector" />
				</detector>
			</componentList>
			<softwareRef ref="so_in_0" />
		</instrumentConfiguration>
	</instrumentConfigurationList>
	<dataProcessingList count="1">
		<dataProcessing id="dp_sp_0">
			<processingMethod order="0" softwareRef="so_dp_sp_0_pm_0">
				<cvParam cvRef="MS" accession="MS:1000544" name="Conversion to mzML" />
			</processingMethod>
			<processingMethod order="0" softwareRef="so_dp_sp_0_pm_1">
				<cvParam cvRef="MS" accession="MS:1000035" name="peak picking" />
				<cvParam cvRef="MS" accession="MS:1000747" name="completion time" value="2019-07-22+10:54" />
				<userParam name="parameter: in" type="xsd:string" value="UK222.mzML"/>
				<userParam name="parameter: out" type="xsd:string" value="UK222_picked.mzML"/>
				<userParam name="parameter: threads" type="xsd:integer" value="1"/>
				<userParam name="parameter: algorithm:signal_to_noise" type="xsd:double" value="0.0"/>
				<userParam name="parameter: algorithm:spacing_difference_gap" type="xsd:double" value="4.0"/>
				<userParam name="parameter: algorithm:ms_levels" type="xsd:string" value="[]"/>
				<userParam name="parameter: algorithm:SignalToNoise:win_len" type="xsd:double" value="200.0"/>
				<userParam name="parameter: algorithm:SignalToNoise:noise_for_empty_window" type="xsd:double" value="1.0e20"/>
			</processingMethod>
		</dataProcessing>
	</dataProcessingList>
	<run id="ru_0" defaultInstrumentConfigurationRef="ic_0" sampleRef="sa_0" startTimeStamp="2016-11-18T23:31:16" defaultSourceFileRef="sf_ru_0">
		<userParam name="mzml_id" type="xsd:string" value="UK222"/>
		<spectrumList count="2" defaultDataProcessingRef="dp_sp_0">
			<spectrum id="controllerType=0 controllerNumber=1 scan=1" index="0" defaultArrayLength="3" dataProcessingRef="dp_sp_0">
				<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
				<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="1" />
				<cvParam cvRef="MS" accession="MS:1000579" name="MS1 spectrum" />
				<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
				<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="371.101997346613018" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="2.027586375e06" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="5.761997102219e07" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="297.014933353315996" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="2020.207901882424039" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
				<userParam name="filter string" type="xsd:string" value="FTMS + p NSI Full lock ms [300.0000-2000.0000]"/>
				<userParam name="preset scan configuration" type="xsd:string" value="1"/>
				<scanList count="1">
					<cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
					<scan >
						<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="60.118806" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
						<userParam name="MS:1000927" type="xsd:double" value="9.999999776483" unitAccession="UO:0000028" unitName="millisecond" unitCvRef="UO"/>
						<scanWindowList count="1">
							<scanWindow>
								<cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="300" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
								<cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="2000" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
							</scanWindow>
						</scanWindowList>
					</scan>
				</scanList>
				<binaryDataArrayList count="2">
					<binaryDataArray encodedLength="32">
						<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
						<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
						<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
						<binary>tT/4x6Exd0BSuB6F69F7QCcnO+TUkJ9A</binary>
					</binaryDataArray>
					<binaryDataArray encodedLength="16">
						<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
						<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
						<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
						<binary>E4L3SQAgekQAAKJB</binary>
					</binaryDataArray>
				</binaryDataArrayList>
			</spectrum>
			<spectrum id="controllerType=0 controllerNumber=1 scan=4" index="1" defaultArrayLength="3">
				<cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
				<cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
				<cvParam cvRef="MS" accession="MS:1000580" name="MSn spectrum" />
				<cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
				<cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="107.96688079834" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="5531.5947265625" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="19911.6002197265625" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="107.96688079834" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
				<cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="704.372009277343977" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS"/>
				<userParam name="filter string" type="xsd:string" value="FTMS + c NSI d Full ms2 684.2034@hcd25.00 [94.3333-1415.0000]"/>
				<userParam name="preset scan configuration" type="xsd:string" value="2"/>
				<scanList count="1">
					<cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
					<scan >
						<cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="61.007892" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
						<userParam name="MS:1000927" type="xsd:double" value="79.999998211860998" unitAccession="UO:0000028" unitName="millisecond" unitCvRef="UO"/>
						<userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="0.0"/>
						<scanWindowList count="1">
							<scanWindow>
								<cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="94.333335876465" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
								<cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="1415" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
							</scanWindow>
						</scanWindowList>
					</scan>
				</scanList>
				<precursorList count="1">
					<precursor>
						<isolationWindow>
							<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="684.203369140625" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
							<cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
							<cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="1" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
						</isolationWindow>
						<selectedIonList count="1">
							<selectedIon>
								<cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="684.203369140625" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
								<cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
							</selectedIon>
						</selectedIonList>
						<activation>
							<cvParam cvRef="MS" accession="MS:1000422" name="beam-type collision-induced dissociation" />
							<cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="25.0" unitAccession="UO:0000266" unitName="electronvolt" unitCvRef="UO"/>
						</activation>
					</precursor>
				</precursorList>
				<binaryDataArrayList count="2">
					<binaryDataArray encodedLength="32">
						<cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
						<cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
						<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
						<binary>CwAAYOH9WkCOj7ut+dFyQAIAAOD5AoZA</binary>
					</binaryDataArray>
					<binaryDataArray encodedLength="16">
						<cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
						<cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
						<cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
						<binary>wtysRZrZXEQAoFxE</binary>
					</binaryDataArray>
				</binaryDataArrayList>
			</spectrum>
		</spectrumList>
	</run>
</mzML>
"#;

/// The template's two records repeated to `count` spectra, alternating MS1 and
/// MS2, with unique native IDs and increasing retention times.
fn scaled(count: usize) -> MSExperiment {
    let mut experiment = mzml::read(Cursor::new(TEMPLATE)).expect("template reads");
    let template: Vec<MSSpectrum> = std::mem::take(&mut experiment.spectra);
    experiment.spectra = (0..count)
        .map(|index| {
            let mut spectrum = template[index % 2].clone();
            spectrum.native_id = format!("controllerType=0 controllerNumber=1 scan={}", index + 1);
            spectrum.rt = 60.0 + index as f64 * 0.25;
            spectrum
        })
        .collect();
    experiment
}

#[test]
fn fifty_thousand_realistic_records_write_through_every_default_entry_point() {
    let experiment = scaled(50_000);
    mzml::write(io::sink(), &experiment).expect("mzml::write");
    mzml::write_with_options(io::sink(), &experiment, &mzml::WriteOptions::default())
        .expect("mzml::write_with_options");
    for write_index in [false, true] {
        let mut options = PeakFileOptions::default();
        options.write_index = write_index;
        mzml::write_with_peak_options(io::sink(), &experiment, &options)
            .expect("mzml::write_with_peak_options");
    }
}

#[test]
fn explicit_whole_document_ceilings_still_refuse_and_derived_ones_scale() {
    use mzml::PeakWriteLimits;
    let experiment = scaled(20_000);
    // The former fixed default is still a valid explicit whole-document
    // ceiling, and it still refuses this experiment before writing anything.
    let mut untouched = b"existing output".to_vec();
    assert!(
        mzml::write_with_peak_options_and_limits(
            &mut untouched,
            &experiment,
            &PeakFileOptions::default(),
            &PeakWriteLimits::default(),
        )
        .is_err()
    );
    assert_eq!(untouched, b"existing output");
    // Derived ceilings equal the fixed ones for an empty experiment and grow
    // with the record count and the stored values.
    let empty = PeakWriteLimits::for_experiment(&MSExperiment::default());
    let base = PeakWriteLimits::default();
    assert_eq!(empty.max_xml_bytes, base.max_xml_bytes);
    assert_eq!(empty.max_work, base.max_work);
    assert_eq!(empty.max_bytes, base.max_bytes);
    let derived = PeakWriteLimits::for_experiment(&experiment);
    assert!(derived.max_work > base.max_work.saturating_mul(20_000));
    assert!(derived.max_bytes > base.max_bytes.saturating_mul(20_000));
    assert!(derived.max_xml_bytes <= u64::MAX / 8);
}

#[test]
fn twenty_thousand_record_file_roundtrips_through_the_file_handler() {
    use openms::format::file_handler::FileHandler;
    use openms::format::file_types::FileType;
    let experiment = scaled(20_000);
    let directory =
        openms::system::file::TempDir::new_in(std::env::temp_dir(), false).expect("temp dir");
    let path = directory.path().join("scaled.mzML");
    FileHandler::store_experiment(&path, &experiment, Some(FileType::MzMl)).expect("store");
    let file = io::BufReader::new(std::fs::File::open(&path).expect("open"));
    let back = mzml::read_with_options(file, &generous_read_options()).expect("read");
    assert_eq!(back.spectra.len(), experiment.spectra.len());
    assert_eq!(back.spectra, experiment.spectra);
}

#[test]
fn default_writer_emits_indexed_mzml_with_verified_offsets_and_sha1() {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let mut experiment = scaled(2_000);
    experiment.chromatograms.push(openms::MSChromatogram {
        native_id: "TIC".into(),
        ..Default::default()
    });
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &experiment).expect("write");
    let text = std::str::from_utf8(&bytes).expect("UTF-8");
    assert!(text.contains("<indexedmzML "));
    assert!(text.ends_with("</indexedmzML>\n"));
    // Independent Python oracle: SHA-1 of the prefix through `<fileChecksum>`,
    // every offset addressing its `<spectrum`/`<chromatogram` tag, in order.
    let script =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tools/mzml_writing/check_output.py");
    let mut child = Command::new("python3")
        .arg(script)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("python3 independent oracle");
    child
        .stdin
        .take()
        .expect("stdin")
        .write_all(&bytes)
        .expect("pipe");
    let result = child.wait_with_output().expect("oracle");
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let back = mzml::read(Cursor::new(&bytes)).expect("read back");
    assert_eq!(back.spectra, experiment.spectra);
    // Nothing to index: an empty experiment stays plain mzML (CPP-050).
    let mut empty = Vec::new();
    mzml::write(&mut empty, &MSExperiment::default()).expect("empty");
    assert!(!std::str::from_utf8(&empty).unwrap().contains("indexedmzML"));
}

/// Header parity with the C++ Release writer on the benchmark input's own
/// metadata. Every literal here is a byte of the C++ Release output of
/// MapNormalizer on `inputs/derived/sub_centroid_uk222_picked_first600.mzML`
/// (smoke-run INI, 2026-09-14 prefix), whose header is that of the input.
#[test]
fn header_matches_the_cpp_release_writer_on_the_benchmark_input_metadata() {
    let experiment = scaled(2);
    let mut bytes = Vec::new();
    mzml::write(&mut bytes, &experiment).expect("write");
    let xml = String::from_utf8(bytes).expect("UTF-8");
    for expected in [
        // Root element: schema location and accession, never an `id`.
        "<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" \
         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" \
         xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml \
         http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" accession=\"\" version=\"1.1.0\">",
        // The pinned cvList.
        "<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" \
         URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\"/>",
        "<cv id=\"BTO\" fullName=\"BrendaTissue545\" version=\"unknown\" \
         URI=\"http://www.brenda-enzymes.info/ontology/tissue/tree/update/update_files/BrendaTissueOBO\"/>",
        // Valueless cvParams carry no `value` attribute.
        "<cvParam cvRef=\"MS\" accession=\"MS:1000579\" name=\"MS1 spectrum\"/>",
        "<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\"/>",
        // The instrument's software comes first and is `so_in_0`; the fallback
        // `so_default` is the source's empty Software().
        "<software id=\"so_in_0\" version=\"2.8-280502/2.8.1.2806\">",
        "<software id=\"so_default\" version=\"\">\n<cvParam cvRef=\"MS\" \
         accession=\"MS:1000799\" name=\"custom unreleased software tool\" value=\"\"/>",
        "<softwareRef ref=\"so_in_0\"/>",
        // Every processing method is order="0", as the source writes it.
        "<processingMethod order=\"0\"",
        // Inherited processing parameters keep the source's number text.
        "<userParam name=\"parameter: algorithm:signal_to_noise\" type=\"xsd:double\" value=\"0.0\"/>",
        "<userParam name=\"parameter: algorithm:SignalToNoise:win_len\" type=\"xsd:double\" \
         value=\"200.0\"/>",
        "<userParam name=\"parameter: algorithm:SignalToNoise:noise_for_empty_window\" \
         type=\"xsd:double\" value=\"1.0e20\"/>",
        // The run keeps the source identifier, its default source file, and
        // carries the document id as the userParam the source writes.
        "<run id=\"ru_0\" defaultInstrumentConfigurationRef=\"ic_0\" sampleRef=\"sa_0\" \
         startTimeStamp=\"2016-11-18T23:31:16\" defaultSourceFileRef=\"sf_00000000000000000000\">",
        "<userParam name=\"mzml_id\" type=\"xsd:string\" value=\"UK222\"/>",
        // The intensity array carries the source's counts unit.
        "<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\" unitCvRef=\"MS\" \
         unitAccession=\"MS:1000131\" unitName=\"number of detector counts\"/>",
    ] {
        let expected: String = expected.split_whitespace().collect::<Vec<_>>().join(" ");
        assert!(xml.replace('\n', " ").contains(&expected), "{expected}");
    }
    // The document `id` is not a root attribute any more.
    assert!(!xml.contains("<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" id="));
    // Records inherit the list's default processing reference; only the list
    // declares one. (The source repeats it on the first spectrum.)
    assert_eq!(xml.matches(" dataProcessingRef=").count(), 0);
    assert_eq!(xml.matches("defaultDataProcessingRef=").count(), 1);
    assert_eq!(
        mzml::read(Cursor::new(xml.as_bytes()))
            .expect("read back")
            .settings,
        experiment.settings
    );
}

/// Reader options with every size ceiling lifted. The reader's own default
/// ceilings belong to a separate lane; these tests are about the writer.
fn generous_read_options() -> mzml::ReadOptions {
    mzml::ReadOptions {
        max_xml_bytes: u64::MAX / 8,
        max_array_bytes: usize::MAX / 8,
        max_total_peaks: usize::MAX / 8,
        max_total_array_bytes: usize::MAX / 8,
        max_total_array_elements: usize::MAX / 8,
        max_total_params: usize::MAX / 8,
        max_param_bytes: usize::MAX / 8,
        ..Default::default()
    }
}

/// Staged benchmark inputs; present only on the IBMI HPC nodes.
const BENCH_INPUTS: &str = "/ceph/ibmi/abi/oliver/bench/openms4/inputs";

/// Reads a staged benchmark input, stores it through the tool path
/// (`FileHandler::store_experiment`) and reads the output back.
fn store_and_reload_bench_input(relative: &str, expected_spectra: usize) {
    use openms::format::file_handler::FileHandler;
    use openms::format::file_types::FileType;
    let options = generous_read_options();
    let input = std::path::Path::new(BENCH_INPUTS).join(relative);
    let file = io::BufReader::new(std::fs::File::open(&input).expect("staged benchmark input"));
    let started = std::time::Instant::now();
    let experiment = mzml::read_with_options(file, &options).expect("read input");
    let read_time = started.elapsed();
    assert_eq!(experiment.spectra.len(), expected_spectra);
    let directory =
        openms::system::file::TempDir::new_in(std::env::temp_dir(), false).expect("temp dir");
    let output = directory.path().join("stored.mzML");
    let started = std::time::Instant::now();
    FileHandler::store_experiment(&output, &experiment, Some(FileType::MzMl)).expect("store");
    let store_time = started.elapsed();
    // What the index and the checksum cost: the same document written plain.
    let plain = directory.path().join("plain.mzML");
    let started = std::time::Instant::now();
    let mut file = io::BufWriter::new(std::fs::File::create(&plain).expect("create"));
    mzml::write_with_options(&mut file, &experiment, &mzml::WriteOptions::default())
        .expect("plain store");
    drop(file);
    let plain_time = started.elapsed();
    println!(
        "{relative}: read {:.2} s, indexed store {:.2} s, plain store {:.2} s, \
         {} bytes in, {} bytes out",
        read_time.as_secs_f64(),
        store_time.as_secs_f64(),
        plain_time.as_secs_f64(),
        std::fs::metadata(&input).map(|m| m.len()).unwrap_or(0),
        std::fs::metadata(&output).map(|m| m.len()).unwrap_or(0),
    );
    std::fs::remove_file(&plain).ok();
    let file = io::BufReader::new(std::fs::File::open(&output).expect("open output"));
    let back = mzml::read_with_options(file, &options).expect("read output");
    assert_eq!(back.spectra.len(), expected_spectra);
    for (stored, original) in back.spectra.iter().zip(&experiment.spectra) {
        assert_eq!(stored.native_id, original.native_id);
        assert_eq!(stored.peaks, original.peaks);
    }
}

#[test]
#[ignore = "HPC only: reads the 547 MB UK222_picked benchmark input from /ceph"]
fn hpc_benchmark_centroid_uk222_picked_stores_through_the_tool_path() {
    store_and_reload_bench_input(
        "centroid_lcms_qe_silac_uk222_picked/UK222_picked.mzML",
        40_856,
    );
}

#[test]
#[ignore = "HPC only: reads the 2.3 GB UK222 profile benchmark input from /ceph"]
fn hpc_benchmark_profile_uk222_stores_through_the_tool_path() {
    store_and_reload_bench_input("profile_hr_qe_silac_uk222/UK222.mzML", 40_856);
}
