# Conditional design: perturbational influence, internal compatibility and gluing

This is a proposed construction, not a fitted model or a DATA GO. It applies only after simultaneous physiological arrays and the intervention/replicate joins pass audit. No H0, regional energy or pushforward was computed.

The cleanest mechanistic pilot is locally driven voltage in a small coupled ensemble. Elvers supplies a published example of this design, but its numerical cohort is not yet public/inspected and the documented local-clamp assay observes only two cells; three-cell spontaneous recordings are a different assay. Do not manufacture a jointly observed three-cell perturbation experiment by joining separate animals. A two-cell pilot is an allowed smaller-system exception: it can establish influence, intrinsic dynamics and compatibility across one boundary. It cannot by itself establish a distinct multiport/hierarchical integration effect beyond a single-edge conventional model. A third observed port is needed only for that specific multicellular/multi-input extension, not for DATA GO or a valuable local physiological benchmark. Simultaneous third-cell clamp data are not universally required if a different genuinely independent validation target and jointly observed regional states are available. Whether the source camera includes such cells must be established before acquisition.

## Measured states and observation operators

Let V_i(t) be membrane potential only when that is the underlying modality. For relative GEVI measurements write y_i(t)=s_i V_i(t)+o_i+noise, or a calibrated nonlinear sensor response. Relative voltage is sufficient; absolute millivolts are not mandatory. Unknown cell-specific gains/offsets prevent treating equal fluorescence as equal voltage or identifying conductance uniquely. Use an effective optical dynamical model if gains cannot be independently calibrated and label it accordingly.

Keep original fluorescence, background, drift-corrected signal and available detection/observation masks as distinct representations. A separate binary mask is not mandatory if complete frame acquisition and absence of interpolation are established by the source pipeline; the gate is relevant provenance, not a demand for uncollected raw videos. The observation map M_a(t) encodes actual recorded/artefact-valid cells in biological unit a. Never infer a genuine-observation mask from finiteness, fill a silent cell roster from maximum event ID, or count interpolation as measurement. Time-varying optical transforms must be fixed from calibration/baseline, with uncertainty, not learned from held-out responses.

For astrocytes replace V by the actually measured calcium-related state and include indicator kinetics if supported. Calcium/IP3 propagation is not membrane-voltage diffusion. Event onsets/marks are a legitimate different observation modality, but an event-only model needs the entire observation/cell roster, explicit exposure time and detection model. The compact Cahill exporter omits those full-state fields; no continuous waveform may be fabricated from peak/decay summaries.

## Topology and local equations

Candidate edges come from independently measured cell contacts, anatomical identities and perturbation-supported coupling. Contact is support, not proof of a conducting junction. Separate uncertain, directed chemical and candidate junctional edges. Genotype/protein expression alone is not a measured edge weight. A direct clamp-response transfer can support an influence relation, but indirect paths remain possible.

For short, subthreshold voltage responses, a plausible physical model is

C_i dV_i/dt = −g_i(V_i−E_i) − sum_j g_ij(V_i−V_j) + I_i^ext + eta_i.

C_i, intrinsic leak/excitability and junctional conductance must be distinguished. Innexin deletion can affect excitability as well as coupling. If junction current/calibration is not measured, fit predictive effective coefficients; do not report them in siemens. APs, nonlinear channels or rectification may invalidate a passive model, in which case use a separately declared nonlinear/state-space constraint and compare standard baselines fairly. Closed-loop optical clamp inputs depend on recorded state; keep command, achieved state and actuator output separate and do not estimate an open-loop causal kernel from feedback logs as though inputs were random.

A finite-history effective factor can be

r_it = x_i(t+1) − A_i h_i(t) − sum_j B_ij h_j(t) − U_i u(t) − b_i,

where h contains only pre-target samples, u includes known stimulation/context, and all maps are trained independently of evaluation units. Multi-input factors retain the full receiving-cell relation rather than splitting a common residual into fictitious pairwise errors.

## Stalks, restrictions and compatibility operator

Variable stalks contain observed cell histories and their next state, plus known input coordinates or explicitly marked latent states. History-copy factors equate overlapping samples. A factor stalk is its affine solution space {z_f:A_f z_f=b_f}; incidence-edge stalks are its cell/input ports. Factor restrictions project solutions to ports; variable restrictions are identities/copies. The affine coboundary equates shared ports; eliminating factor copies gives D_theta z=b_theta. Homogenisation, if used, fixes the added coordinate to one so zero cannot satisfy every affine relation trivially.

A fitted forward dynamical system normally has mathematically compatible trajectories for arbitrary initial states. Therefore a zero set or minimised unanchored residual is not biological validation. Anchor to independent observations and report constraint rank, supported degrees and missing coverage. Latent unconstrained ports must not absorb all disagreement.

## Three distinct empirical quantities

1. **Perturbational influence:** on held-out biological units measure the impulse/step response in B to a controlled input in A relative to sham and pre-input history. Contrast control, communication disruption and rescue, accounting for actual achieved A response and intrinsic excitability. Keep a response-amplitude/transfer statistic separate from compatibility.
2. **Within-region compatibility:** use control-trained, independently noise-calibrated internal factors in each anatomical region. For held-out residuals e_f, E_R=sum observed whitened residual squares / supported scalar constraint degrees. Zero supported degrees is NOT EVALUABLE. Report the biological-unit distribution and measurement support; do not treat low activity or few observed edges as high compatibility.
3. **Cross-region gluing:** maintain boundary ports and compare anchored regional sections with a jointly constrained section. J_sep is the minimum observation-distance cost with boundary copies allowed to differ; J_joint imposes equality of those copies. Allocate each observation's total weight once across regions, including any overlapping boundary. Then J_joint−J_sep>=0 by nested feasible sets. Normalise by declared boundary degrees and use coverage/noise-matched controls. This difference can vanish because no boundary is observed or constraints are underidentified; that is not integration.

Weak or absent B_ij need not create incompatibility. An uncoupled network can satisfy its own coherent dynamics. Evaluate both (i) frozen **control** boundary maps after disruption and (ii) independently recalibrated disrupted maps. The first tests loss of control-like communication; the second distinguishes coherent reorganisation from failure of the new local dynamics. A fall in influence with preserved internal coherence is a possible signature, not an assumed outcome. Increased coupling need not improve function; an excessively synchronised preparation can be dysfunctional.

## Falsification and conventional equivalences

Require held-out perturbation propagation beyond receiver-only/context predictors before treating effective maps as biological constraints. A linear factor energy can be exactly the innovation loss of a VAR/state-space model; an identity-restriction sheaf yields a graph Laplacian. Neither equivalence is a new contribution. Compare with dynamic correlations, functional connectivity, graph components/modularity/Laplacians, standard state-space/VAR, and ordinary perturbation transfer maps on matched observations and tuning budgets.

The distinct prospective claim is a multi-region boundary/gluing statistic that explains or predicts independent perturbation/outcome variation beyond internal innovations, pairwise correlations, activity amplitude and graph structure. A purely relabelled conventional residual is a negative methodological result. A three-cell/two-port benchmark can validate bookkeeping and influence but cannot establish extensive multiscale organisation. Larger coverage, ageing, complete RNA or absolute voltage are not prerequisites for that first benchmark; simultaneously observed boundary dynamics are.
