# MRM feature selection by linear programming

Rust ports of `ANALYSIS/OPENSWATH/MRMFeatureSelector.h` and
`ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h` (pin `bc9cc12`):

| Header | Rust |
| --- | --- |
| `MRMFeatureSelector.h` | [`src/analysis/openswath/mrm_feature_selector.rs`](../src/analysis/openswath/mrm_feature_selector.rs) |
| `MRMBatchFeatureSelector.h` | [`src/analysis/openswath/mrm_batch_feature_selector.rs`](../src/analysis/openswath/mrm_batch_feature_selector.rs) |

Tests: [`tests/mrm_feature_selector.rs`](../tests/mrm_feature_selector.rs) (12),
[`tests/mrm_oracle.rs`](../tests/mrm_oracle.rs) (executed comparison).
Manifest: [`tests/data/mrm_feature_provenance.json`](../tests/data/mrm_feature_provenance.json).

## The solver gap

The source builds its programs in `LPWrapper` (`DATASTRUCTURES/LPWrapper.h`) and
solves them with GLPK, COIN-OR or HiGHS, whichever the build links. That header
is not ported: the solver choice for `LPWrapper` is a recorded backlog item of
the project (`GATED-01`). What is ported is everything the selectors do around
the solve. The port builds the program into [`LinearProgram`], a minimal model
of the problem-building half of `LPWrapper` (columns with bounds, type and
objective, rows with coefficients and bounds), and hands it to an [`LpSolver`].

- [`PartitionSolver`] is exact for the programs of
  `MRMFeatureSelectorScore`: rows `sum(x) = 1` over disjoint columns bounded to
  `[0, 1]`, for which the continuous and the integer optimum is the column of
  least objective in each row. Ties choose the **last** tied column: the source
  leaves the choice to the solver's pivoting, and the executed Release build's
  GLPK took the last of two (`adn`) and of three (`adpglc`) tied candidates on
  the upstream 210-feature file in six of the seven parameter sets. In the
  seventh (continuous variables, windows of 8 stepping 4) the two overlapping
  segments of a group got different tied columns from GLPK, which a solver that
  looks at one row at a time cannot reproduce; the oracle test compares that run
  group by group (same groups, every pick optimal, differences only in tied
  groups). That is a measurement of one solver on one data set, not a rule GLPK
  states.
- `MRMFeatureSelectorQMIP` builds its program in full (the `qp` and `ABS`
  columns and the five constraints per pair), but needs a general mixed integer
  solver, which this crate does not ship. A caller supplies one through
  `LpSolver`.

Both headers therefore stay **partial**: `MRMFeatureSelector.h` until a solver
is chosen for `LPWrapper`, `MRMBatchFeatureSelector.h` because its
`batchMRMFeaturesQMIP` runs the QMIP selector.

## API mapping

| Source member | Rust |
| --- | --- |
| `MRMFeatureSelector` (abstract) | trait `MRMFeatureSelector` |
| `VariableType` (`INTEGER = 1`, `CONTINUOUS`) | `VariableType::{Integer, Continuous}` |
| `LambdaScore` (`LINEAR = 1`, `INVERSE`, `LOG`, `INVERSE_LOG`, `INVERSE_LOG10`) | `LambdaScore::{Linear, Inverse, Log, InverseLog, InverseLog10}` |
| `SelectorParameters` and its defaults | `SelectorParameters` and `Default` |
| `optimize` (pure virtual) | `MRMFeatureSelector::optimize` |
| `selectMRMFeature` | `MRMFeatureSelector::select_mrm_feature` (provided method) |
| `addVariable_` | `add_variable` |
| `addConstraint_` | `add_constraint` |
| `computeScore_` | `compute_score` |
| `constructTargTransList_` | `construct_targ_trans_list` |
| `weightScore_` | `weight_score` |
| `removeSpaces_` | `remove_spaces` |
| `MRMFeatureSelectorQMIP::optimize` | `MRMFeatureSelectorQMIP::optimize`, `build_problem` |
| `MRMFeatureSelectorScore::optimize` | `MRMFeatureSelectorScore::optimize`, `build_problem` |
| `MRMFeatureSelector_test` (test facade) | not part of the API; the free functions above are public |
| `MRMBatchFeatureSelector::batchMRMFeatures` | `batch_mrm_features` |
| `MRMBatchFeatureSelector::batchMRMFeaturesScore` | `batch_mrm_features_score` |
| `MRMBatchFeatureSelector::batchMRMFeaturesQMIP` | `batch_mrm_features_qmip` |
| `LPWrapper` (the used members) | `LinearProgram`, `Column`, `Row`, `BoundType`, `ColumnType`, `Sense`, `LpSolver`, `LpSolution`, `SolverStatus` |

## Preserved source conventions

- **Variable names.** A candidate is `"<name>_<unique id>"`; a pair variable is
  `"<name1>_<i>-<name2>_<j>"` with the positions in the candidate lists, and
  its absolute value variable adds `-ABS`.
- **`constructTargTransList_` tests `names.contains(component_name)` where it
  means the negation.** With `select_transition_group = false` a component name
  enters `time_to_name` only when it was already seen, and then once per further
  occurrence; a first occurrence is never added. Kept (defect candidate).
- **Selection thresholds.** The Score selector keeps a column at or above
  `optimal_threshold`, the QMIP selector strictly above it and only candidate
  columns.
- **`segment_window_length`/`segment_step_length` of `-1`** mean one segment of
  everything. Other negative values convert to huge `Size` values, as in the
  source.
- **`std::sort` of `(retention time, name)` pairs** is the total order on the
  pair; equal pairs are indistinguishable.
- **Scores multiply in key order** of the weight map, a missing metavalue is
  skipped and a transformed value that is not positive, NaN or infinite is left
  out; more than one weight takes the `n`-th root in the QMIP program.

## Native differences

| Where | Source | Here |
| --- | --- | --- |
| solver | `LPWrapper` with GLPK / COIN-OR / HiGHS | the `LpSolver` trait; `PartitionSolver` for the Score program |
| failed solve | `Exception::FailedAPICall` | `Error::InvalidValue` naming the status |
| `segment_step_length == 0` | integer division by zero | `Error::InvalidValue` |
| unknown name in `time_to_name` | `std::map::at` throws | `Error::InvalidValue` |
| row names | LPWrapper passes duplicate row names to GLPK, which **aborts** the process (`** duplicate name ..._constraint`, exit 134) | no abort; the second row of a name has no column left (its candidate exists already), so the program is infeasible and the selection is refused. With `select_transition_group = false` the source reaches this whenever a component name repeats a group name (retained: `tests/data/mrm_oracle/selector_dupname.*`) |
| `getColumnIndex` of a repeated name | solver dependent | the first column of that name |
| `LambdaScore` outside the enumeration | `Exception::IllegalArgument` | unreachable (a closed enum) |
| `OPENMS_LOG_WARN` for a missing metavalue | a warning | none |
| `MRMBatchFeatureSelector` | a class with deleted constructor and destructor | three free functions |

## Checked boundaries and evidence

| Behaviour | Tier | Source |
| --- | --- | --- |
| `selectMRMFeature` of the Score selector on the 210-feature file (37 selected, the first and tenth subordinate), `removeSpaces_`, `constructTargTransList_`, `weightScore_`, `computeScore_` | 3 | `MRMFeatureSelector_test.cpp` in `tests/mrm_feature_selector.rs` |
| the Score selector's selection on seven parameter sets over both upstream files (integer and continuous, windows, weight sets, `select_transition_group` false, no weights) | 1 | `tests/mrm_oracle.rs` against `../oracle/O3b-mrm` |
| the features the QMIP batches keep (all of the small file, one candidate per group) | 1 | same |
| the abort on a repeated row name | 1 | same, exit code and output retained |
| the QMIP program's columns, rows, bounds and coefficients for a small instance, derived from reading `MRMFeatureSelectorQMIP::optimize` | 4 | `tests/mrm_feature_selector.rs` |
| segment arithmetic, failed solve, the inverted `contains` | 4 | same |

The class test's QMIP sections are not run through a solver here: the small
file has one candidate per component group, so the optimum is forced, and the
large continuous one is documented upstream as flaky between compilers. The
program the port builds for them is checked structurally.
