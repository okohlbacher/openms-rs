#!/usr/bin/env python3
"""Project only the SemanticValidator::locateTerm membership query, not validation."""
from pathlib import Path
import xml.etree.ElementTree as ET
import json,hashlib,sys
sdk=Path(sys.argv[1]);dest=Path(sys.argv[2]);src=sdk/'share/OpenMS/MAPPING/ms-mapping.xml'
root=ET.parse(src).getroot();owners={'run','contact','sourceFile','sample','software','instrumentConfiguration','source','analyzer','detector','processingMethod','binaryDataArray'}
terms={}
for rule in root.iter('CvMappingRule'):
 path=rule.attrib['cvElementPath'];owner=path.split('/')[-3]
 if owner not in owners:continue
 for term in rule.findall('CvTerm'):
  a=term.attrib;row=(a['termAccession'],a['useTerm']=='true',a['allowChildren']=='true')
  terms.setdefault(owner,[]).append(row)
fixture={'source_commit':'82ce5b373c97f934ffd9b1ffd80215ca66473d0b','source_path':str(src.relative_to(sdk)),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'query':'union of matching path rules, direct useTerm or strict allowChildren ancestry; does not validate MUST/cardinality/value/units','owners':terms}
(dest/'tests/data/mzml_header/mapping_rules.json').write_text(json.dumps(fixture,indent=2)+'\n')
lines=['// Copyright (c) 2002-present, OpenMS Inc.','// SPDX-License-Identifier: BSD-3-Clause','// $Maintainer: OpenMS Rust contributors $','//! Pinned source mapping-path membership only; not full semantic validation.','use super::*;','pub(super) fn permitted(owner: &str,id: &str,work: &mut Work)->Result<bool> {','    let rules: &[(&str,bool,bool)] = match owner {']
for owner,rows in sorted(terms.items()):
 lines.append('        '+json.dumps(owner)+' => &[')
 for accession,direct,children in rows:
  lines.append(f'            ("{accession}",{str(direct).lower()},{str(children).lower()}),')
 lines.append('        ],')
lines+=['        _ => return Ok(false),','    };','    for &(term,direct,children) in rules {','        work.charge(1,0)?;','        if (direct && id==term) || (children && work.child(id,term)?) { return Ok(true); }','    }','    Ok(false)','}']
(dest/'src/format/mzml_header/mapping.rs').write_text('\n'.join(lines)+'\n')
