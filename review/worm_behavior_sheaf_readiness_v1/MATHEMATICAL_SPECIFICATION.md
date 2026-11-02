# Prospective dynamical factor sheaf: specification, not a fitted model

This design is conditional on the new A–F readiness gates. No H0, sheaf energy, fitted restrictions or pushforwards are computed in this task. Calcium observations are fluorescence states, not measured voltage or membrane conductance. Anatomical support is a prior, not proof of direct causal action in these animals.

## Variables, factors and restriction maps

Let i index a side-resolved individual neuron, t an actual acquired frame, and a an independent animal. An elementary neuron-time stalk is X_it=R, containing a training-scaled calcium state c_it. A history stalk H_it=R^L contains L causal samples ending at t; initial L=10 after verifying timestamp units and continuity. Overlapping histories have explicit copy/shift consistency constraints. Class-only labels cannot silently instantiate two bilateral neuron stalks; a class-level model would require a separately justified aggregation/latent observation model.

Each receiver update is one multiport factor

r_jt = c_j,t+1 - alpha_j^T h_jt - sum_{i in N(j)} beta_ji^T h_it - gamma_j^T u_t - b_j.

u_t contains known sensory/context and permitted behavioural histories available by the prediction cutoff. N(j) is an independently versioned anatomical adjacency intersected with reliably observed identities, or an explicitly labelled predictive sparse graph selected using training data only. Electrical, chemical and uncertain edges remain separately identified. No adjacency is inferred merely from proximity/correlation. Graph omissions in an incompletely annotated reference are not measured absences.

Stack update, history-copy and declared boundary equations as D_theta z=b_theta. D is the affine compatibility operator. For a conventional cellular sheaf on the incidence graph, variable vertices have their state stalks; a factor vertex has the local solution space {v:A_f v=0}. Incidence-edge stalks are the relevant variable spaces. Factor-to-edge restrictions project to ports, variable-to-edge restrictions are identities, and the coboundary equates copies. For b≠0 use affine factor stalks or homogenize with a coordinate fixed to one. An unconstrained extra coordinate must not let every affine equation collapse to zero.

A global section assigns shared state values satisfying all local factors and copy constraints. Its existence is a mathematical property of the chosen maps, not evidence that those maps are biologically correct. A forward-time state-space system generally has compatible trajectories for any initial condition; its unobserved zero set is therefore often nonempty by construction. Nontrivial empirical discrepancies require anchoring to observations and independently calibrated noise. Do not claim a physiological discovery because an unconstrained dynamical system has zero residual.

## Observation masks and latent states

Write y=M H z+eta. M encodes actual availability, identity uncertainty and verified artefact masks; reconstructed/interpolated samples have separate provenance and are not counted as observations. H maps latent variables to measured calcium if such a model is used. Timestamps, clock offsets and idx_splits semantics must be verified before constructing histories. No window crosses a missing interval or session boundary. All scaling, factor coefficients, noise models and latent transformations use training animals only.

Primary diagnostics omit unsupported histories rather than invent values. A causal filter can integrate missing states only if the same filter/information set is used in comparators and its conditional uncertainty is propagated. Future-data smoothing is a different task. Unobserved neurons can create effective interactions between observed neurons; they cannot be treated as observed direct partners. Free latent-state models also have similarity/gauge ambiguities, so their coordinates and couplings are not uniquely physiological.

## Coupling and independent prediction

A group S→j predictive-coupling statistic is the animal-weighted held-out change in forecast loss between the frozen full model and a matched model omitting that group's history, with both trained using the same inner folds. Report signed loss differences and whole-animal intervals. A magnitude or nonzero fitted coefficient alone is not coupling evidence; common inputs, feedback and calcium filtering can explain apparent influence. Joint versus individual group ablations can test predictive redundancy/synergy, but remain noncausal without interventions. Report instability and collinearity rather than select a favourable edge attribution.

Calcium at approximately 1.67 Hz cannot identify fast junction conductances. Even direct optical voltage needs sensor gain/offset and current/capacitance information for physical parameters. For known calibrated voltages and independent g_ij, the physical alternatives would include I_ij=g_ij(V_j−V_i) plus membrane-current balance; neither equal voltages nor an anatomical junction count is an estimate of g_ij. Rectification/heterotypic channels may require nonlinear directional laws. These data do not supply a matched whole-network conductance experiment.

## Regional observational compatibility and joint sections

For held-out, causally predicted residuals e_f with training/calibration covariance S_f, define

E_R = [sum_{f in R, observed} e_f^T S_f^(-1)e_f] / [sum_{f in R, observed} dim(e_f)].

This is a noise-normalized regional innovation energy with expected value near one under a calibrated Gaussian model. Missing parents require marginalization/forecast uncertainty or exclusion with coverage reporting. Zero supported degrees means NOT EVALUABLE, not zero energy. Weights/covariance must not be fitted on the region being tested. Compare calibration across animals and report the observation fraction alongside every regional value.

An anchored section discrepancy is instead

J_R = min_{z:D_R z=b_R} ||W_R^(1/2)(H_R z-y_R)||^2 / n_R.

W comes from independent training noise and n_R counts observed scalar degrees. Specify initial states and boundary ports; otherwise free latent states can trivialize the projection. This is an observation-to-section projection, not held-out forecast performance. If exact equations are inappropriate under process noise, declare a penalized residual likelihood rather than calling the observations an exact section.

To test whether local regional solutions glue jointly, retain boundary ports. Compare the sum of separately minimized anchored costs with the joint minimum enforcing equality of shared boundary states, using the same observations once and a declared boundary allocation. Joint cost cannot decrease when equality constraints are added to the same objective; any claim of positive joint excess must control duplicated observations, effective degrees and noise. A fitted state-space trajectory can produce the same quantity: this is not automatically a new biological invariant.

## Compatibility filtration and H0

Only after independent predictive validation, freeze factor uncertainty, observation support and graph selection. For each observed multiport factor, q_f=e_f^T S_f^(-1)e_f and score s_f=exp(−q_f/(2 d_f)), d_f=observed dimension. Increasing allowed mismatch τ retains factors with q_f/d_f≤τ, yielding a nested incidence/hypergraph filtration. Retain a factor's whole supported port relation, not arbitrary pairwise copies of its total residual. No unvalidated near-zero coefficient may act as a spurious connection solely because it appears among the factor's ports.

H0 gives components of locally compatible constraints; it does not prove the component admits a jointly anchored section. Residual accumulation and cyclic/boundary mismatch need J_R separately. Temporal persistence alone must not count as cross-neuron integration: report temporal-only connectivity separately, and declare any projection from neuron-time components to neurons/windows. Missing/uncertain ports do not supply edges. Changes to observation masks can change components independently of biology; report coverage-matched nulls and downsampling.

## Distance and long-range integration

If actual stable positions and units are unavailable, use independently mapped anatomical hop distance, a prespecified anatomical path-length measure where available, or declared circuit membership. Do not label any of these micrometres. Reference positions are not measured positions in a moving animal. Use a reference distance independent of the inferred compatibility graph, so a component does not determine its own notion of being long-range.

A prospective statistic is the fraction of sufficiently observed reference-neuron pairs at distance≥d whose neuron-time ports belong to the same retained cross-neuron component within a fixed observation window. Denominator, windows, pair support and distance definition are frozen before comparison. Repeated neuron pairs are not independent animals. Threshold sweeps are descriptive until held-out behavioural added value and mask stability are established; more integration is not assumed healthier or more avoidant.

## Functional circuits and pushforwards

Choose circuit regions independently of held-out activity/behaviour, from named anatomical/functional references. A region map π from the original incidence domain to a region domain induces the direct-image sheaf (π_*F)(U)=F(π^(-1)(U)); restrictions follow inclusion. In a finite-poset implementation, region stalks retain compatible sections over the appropriate inverse-image neighbourhood, including shared boundary ports. They are not merely averages of neuronal calcium or a relabelled block adjacency. Overlapping circuits require a cover/incidence representation rather than silently assigning a neuron twice as two independent measurements. No pushforward is constructed before local validation.

## Added value and falsification

An identity-restriction scalar sheaf yields a graph Laplacian. Linear update-factor energy is a conventional VAR/state-space innovation loss. These equivalences are baselines, not evidence of added value. Match causal inputs, parameter counts, tuning budgets, masks and noise calibration when comparing receiver-only autoregression, regularised multivariate models, sensory/behaviour models and low-rank/connectome-constrained dynamics.

Potential value is explicit multiport/state consistency, interpretable boundary gluing and a systematic circuit-level representation across scales. Demonstration would require a distinct anchored gluing or integration summary with robust independent behavioural value beyond conventional activity, correlations and regional innovation statistics, or a prespecified interpretability test with independent circuit perturbations. If the sheaf merely renames the same loss or cannot improve explanation/generalisation, stop the sheaf-specific claim. A decisive future mechanistic experiment perturbs/rescues an independently supported innexin connection or circuit while recording its endpoints and behavioural response in independent animals; current correlations cannot supply that causal evidence.
