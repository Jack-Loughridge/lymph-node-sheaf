# Minimal author exports — drafts only; no messages sent

These requests address precise unresolved fields. They are not requests for terabytes of microscopy and cannot repair a variable never recorded. No contact was sent or form submitted.

## Elvers / Gottschalk: direct coupling benchmark

Request processed QuasAr2 traces and timestamps from the cOVC cohort, genotype/control-strain labels, an anonymous stable animal/session key, trial boundaries, cell1/2 ROI masks or coordinates, clamp setpoint, achieved fluorescence and both actuator channels, calibration/bleaching parameters and any missing/artefact flags. Request the corresponding three-cell spontaneous ensemble arrays with positions and their actual animal IDs so experiments are not incorrectly treated as matched. A CSV/HDF5 export of a few exemplar animals per condition should be small; verify format/size before transfer.

The decisive question: was a third cell recorded with usable voltage during the SAME locally targeted cOVC command? Could it be exported from the existing camera field? If not, simultaneous three-or-more-cell perturbational data require a new recording, not another download. Confirm whether muscle-voltage and behaviour measurements belong to the same animals or merely genotype-matched cohorts. Clarify muscle-specific versus organism-wide innexin loss and available excitability/light controls. Cx36 overexpression is not automatically a rescue of unc-9 loss.

The published future-tense data statement and two example movies are not a structured numerical cohort. Ask for a public processed-data deposit or stable small export, not permission to reinterpret figure means as recordings.

## Cahill / Poskanzer: spatial local-to-network experiment

Request ONE complete representative FOV in WT, genetic Cx43 reduction and CBX (and a matched control FOV), with continuous cell/compartment traces or enough raw signal/ROI information for a small extraction; original and corrected signal definitions; full ROI/cell roster including zero-event cells; cell masks/centroids and physical scale; uncaged cell/landmark; measured laser/Pockels timestamps; anonymous mouse→slice→FOV→paired-NT mapping; genotype/Cre target and neighbour status where measured; original exclusion/ramping/artefact flags. Confirm whether acquisition date is a unique mouse key. Seek explicit cell-count/coverage before moving beyond this FOV.

The compact event exporter is useful but has no continuous waveform or complete cell catalog, its coordinates are event centroids rather than fixed cell ROI geometry, and it uses preset stimulus frames. Please specify whether the12.01GB MATLAB object is7.3/HDF5 and which fields contain independent cell fluorescence and observation masks. A selective HDF5/MAT export of these fields is preferable to transferring the full object. If only event-derived average fluorescence exists, label it as such rather than raw optical measurement. Confirm independent mouse counts and whether controls are littermate/batch matched; FOVs and paired neurotransmitters do not increase biological N.

## Ruhl / Heinemann: voltage and knockout/rescue

Inspect the public Source Data ZIP first. If its traces lack metadata, request a small table joining every trace to dish/culture, recording/session, genotype/construct/rescue, illumination/clock, ROI geometry and simultaneously acquired cell group. Fifty plotted traces must not be presumed a fifty-cell simultaneous network. Identify whether a local voltage/current/optical perturbation with responses outside the targeted cell was recorded. Provide those command/readout traces if they exist; do not substitute global ionomycin treatment for a selective local drive.

Clarify which source traces are raw, background-corrected, cropped, binned or bleach-corrected, whether any samples were interpolated, and how many independent cultures and clone lines contribute per condition. Raw imaging is only needed if processed geometry/simultaneous grouping truly cannot be exported; avoid requesting it initially. Relative GEVI voltage is acceptable, and absolute potential is not required for the pilot.

## Yale / other candidates: stop unless missing primary fields exist

For the endothelial study, ask specifically whether public quantification workbooks contain actual cell×time arrays; request original arrays/ROIs/session clocks with mouse and treatment joins if not. Determine how vascular function assays match calcium recordings and distinguish intrinsic Ca-channel compensation from restored Cx43 coupling.

For Choi request matched control/innexin cell-resolved time series (or a small original movie with segmentation and fly/retina IDs), not more morphology summaries. For Ammer ask for a repository manifest and simultaneous ROI/session grouping before any acquisition. For Gao the missing requirement is independent mechanism-specific coupling confirmation/controls and replicated records, not simply extracting more example-video frames. For beta/delta islets request the canonical processed-data location and mouse/islet/ROI/clock/condition map; a redirecting algorithm repository is not verified physiological data access.
