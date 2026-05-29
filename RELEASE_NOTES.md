# V.0.3.0

## New Features
- Change Mono exon transcript representative TSS/TES selection to guide mode.
- Add new output table for tracking the representative TSS/TES selection for each source transcript.
- Change the output option for merge mode to output prefix instead of output gtf. This allows isomatch output multiple files with one output prefix.
- Use version system for isomx and isoms file. The indexes with outdated version will not work with the new version of isomatch. I made this change because re-build index is cheap. 


## Bug Fixes
- Fixed an issue when counting guide tss/tes hits.