# B8a lane notes (BLOCK-06a: MzMLFile, FeatureXMLFile, StringUtils)

Branch: worktree-wf_ee239125-5be-34. Prereq worktree-wf_ee239125-5be-1 merged (fast-forward to 9fcb962).

## Gaps (from docs/core-sdk-coverage.json scopes)
- MzMLFile.h (partial): IMPeakType (MS:1003441) ion-mobility peak type; selected-ion drift time
  propagation to spectrum; remaining precursor activation metadata/unit routes; progress of
  loadSize, loadBuffer, storeBuffer and both transform overloads.
- FeatureXMLFile.h (partial): "source ZIP input" unsupported (check what C++ actually does).
- StringUtils.h (partial, only FileInfo scope in src/format/file_info/text_format.rs): not ported:
  toStr(float/double,false), appendToStrLowP, long double overloads, toStr(char), toStr(DataValue),
  toStr(ParamValue), numberLength, and the rest of the header (string utilities/conversions).

## Next step
1. Read skill ~/.claude/skills/openms-port-header/SKILL.md
2. Read C++ StringUtils.h at core bc9cc12 member by member; find existing Rust equivalents (grep).
3. MzMLFile.h member-by-member for the named gaps.
4. FeatureXMLFile ZIP input.
