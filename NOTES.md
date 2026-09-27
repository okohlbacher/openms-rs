# A3a NOTES (DEBT-04: refusals resting on "order unspecified")

Tooling: `crate::math::source_sort::{source_sort_permutation, source_sort_by,
source_sort_reversed_by, source_stable_sort_permutation(len, less, TemporaryBuffer::Allocate)}`.

Sites (from brief), status:
1. src/kernel/range_utils.rs:712-722 (total_cmp + refuse inf) -- TODO
2. src/kernel/spectrum_mobility.rs:581-593 sort_by_ion_mobility, :624 presorted -- TODO
3. src/processing/iterative.rs:440/490/505-518 -- TODO
4. src/kernel/faims_helper.rs:57/113 total_cmp order -- TODO
5. src/analysis/feature_finder_picked/debug.rs:186/588 -- TODO
6. src/kernel/chromatogram_tools.rs:194/333 -- TODO
7. src/math/statistic_functions.rs compute_rank / rank_correlation_coefficient (w_idx.sort_by :1388) -- TODO; needs oracle row on NaN keys for pair.second comparator
8. src/processing/feature_overlap_filter.rs:1102-1112 stable_sort emulation -> source_stable_sort_permutation -- TODO

Oracle: ibminode06, ../oracle/A3a-<topic>/, Release install /ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576
Cargo only via gate: OPENMS_GATE_HOST=spock OPENMS_GATE_JOBS=64 ~/.local/bin/openms-kim-gate.sh A3a <args> (background!)

## Next step
1. Read each site + its C++ source call; record exact std::sort/stable_sort call here.
