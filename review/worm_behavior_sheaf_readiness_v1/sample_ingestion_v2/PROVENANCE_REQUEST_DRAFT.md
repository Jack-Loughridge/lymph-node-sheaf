# Small, specific provenance request, if native metadata does not resolve these fields

No message has been sent. First inspect the native dictionary's full metadata-field list and the actual candidate-specific export code/version if available.

The three downloaded HDF5s are numerically intact, and their neural/NIR clock means reconcile. To distinguish genuine observations from reconstruction without reprocessing microscopy, the smallest useful author exports would be:

- A per-final-trace-row/per-time boolean detection mask for both the activity and reference channels, before interpolation, including an index mapping to the public HDF5 rows and frames.
- If those masks are not directly available, the sparse observed `activity_traces` and `marker_traces` dictionaries and final ROI-to-row map used by the exported recordings.
- The exact processing/export code revision and run parameters, including interpolation, bleaching, background/reference division and any denoising settings. Bleach parameters/curves or pre-correction traces would allow assessment of causal preprocessing without whole-test-record normalization.
- The head/camera/stage transformation used to compute `dist_from_cuCenter`, the circle-coordinate convention, and original confocal-to-NIR bin indices or DAQ-derived exposure intervals if readily exportable.

The dictionary's manual identity confidence is distinct from the continuous registration-score array. We will validate native row/class views and left/right labels locally. We are not requesting raw microscopy, new annotations, model fitting or conclusions from the authors.
