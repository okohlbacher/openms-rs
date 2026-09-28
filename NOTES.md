# FF4 lane notes (Biosaur2)

## Scope
- Header: OpenMS/FEATUREFINDER/Biosaur2Algorithm.h (+ source/FEATUREFINDER/Biosaur2Algorithm.cpp) @ core bc9cc12
- Class test: src/tests/class_tests/openms/source/Biosaur2Algorithm_test.cpp @ core bc9cc12
- TOPP consumer: topp 174b576 src/FeatureFinderLFQ.cpp (222 lines, wraps Biosaur2Algorithm)
- TOPP test: test-data 0cb15f2 topp/CMakeLists.txt:429-431 TOPP_FeatureFinderLFQ_1
  (-in FeatureFinderCentroided_1_input.mzML -mini 20 -htol 100 -itol 100 -minlh 2 -cmin 2 -cmax 2;
   expected topp/FeatureFinderLFQ_1_output.featureXML, FuzzyDiff whitelist "id=")
- Also used by ProteomicsLFQ -Seeding:algorithm biosaur2 (not this lane)
- Rust target: src/analysis/featurefinder/biosaur2.rs (append `pub mod` in src/analysis/featurefinder.rs)

## Git access to pins
Worktree guard refuses `git -C`; use `cd <pkg> && git show <sha>:<path>` as separate simple commands.
Pins root: /Users/kohlbach/Claude/OpenMS/OpenMS4-Exploration/OpenMS4-tests/packages

## Next step
1. Read the porting skill, the header, the .cpp and the class test.
2. Survey existing Rust MSExperiment / FeatureMap / Param / featureXML writer APIs.
3. Port the algorithm, then the class tests, then the TOPP differential test.
