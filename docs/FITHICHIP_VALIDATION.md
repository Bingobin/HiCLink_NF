# FitHiChIP integration validation — 2026-11-05

- Official base image: `aylab/fithichip@sha256:bbb2cbcfc4bc61123d5100b111b2a648100da5fe8d24236ca0e4a3b7d8162fb4`.
- Official source commit: `0ea1ac21be870908c672316ffbb630189dc6fae2`; scientific scripts unchanged.
- Installed local image: `hiclink/fithichip:11-0ea1ac2`.
- Installed image ID: `sha256:df3b160c92efa6b7e4e6b584e84529e9a85d5bcc56d8a1504baf3f4270279f1b`.
- Runtime checked on Apple Silicon using `linux/amd64`, Python 3.10.16, cooler 0.10.4 and R 3.6.1 in the container. All required R packages, including data.table, loaded; bedToBigBed ran successfully.
- 21 Python regression tests passed; Nextflow lint and Bash syntax checks passed.
- A real Docker/Nextflow smoke run used three **synthetic** libraries: one directed comparison and one unpaired library. All 12 tasks completed with exit code 0.
- The synthetic differential test contained 43,841 tested contacts. Nine injected treatment-enhanced contacts passed the official differential/FitHiChIP criteria; the annotated table matched the nine native filtered records and their positive direction.
- Each container task retained the actual Docker command and image ID. Container analysis used no network; inputs were mounted read-only.
- Production input discovery found five HMGA2 libraries, all five peaks and the two expected RAD21 comparisons. A production-config preview completed without analysis tasks.
- **No real HMGA2 FitHiChIP caller or differential analysis was started.** Runtime/RAM requirements and biological recovery for those data remain to be evaluated. Current one-library-per-condition contrasts use an assumed BCV, not estimated biological-replicate variance.

The weak-effect fixture also completed successfully with zero significant differential contacts; valid empty results are not converted into synthetic positives. Full native tables, logs and outputs were retained in isolated temporary test directories, separate from the user's existing results.
