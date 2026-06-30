RCTD-Fast Doublet: Exact-Preserving C++ Backend for RCTD Doublet Mode
1. Executive Summary
Problem statement

RCTD doublet mode is too slow for production-scale Xenium cell annotation. Your current baseline is approximately:

Input modality: Xenium cells
Slide canvas: ~100k × 100k pixels
Number of cells/locations: millions
Reference cell types: ~20
Hardware: 8 cores
RCTD setting: MAX_CORES = 8
Runtime: ~10 hours / slide
Primary bottleneck: runtime, not crash/failure
Target: ~50× speedup
Use case: internal production tool
Accuracy requirement: exact identity with RCTD first; speed second

The current RCTD doublet code is structurally expensive. Per bead/cell, it runs a full all-type decomposition, thresholds candidate cell types at all_weights > 0.01, computes singlet scores, scores all candidate pairs, applies doublet/singlet/reject logic, then refits the selected pair. The pair-scoring loop calls decompose_sparse(..., score_mode = TRUE), which uses 25 IRWLS iterations; the final selected doublet refit uses 50 iterations.

The main implementation goal is therefore:

Build run.RCTD.fast.doublet(): a drop-in, exact-preserving replacement for RCTD doublet mode, called from R, with hot kernels implemented in C++ and parallelized over Xenium cells.

RCTD’s own README says doublet-mode outputs are stored in @results$results_df and @results$weights_doublet, with spot_class, first_type, and second_type as key outputs. Those must remain API-compatible.

Non-goals for v1

These are tempting, but they break exact identity:

Non-goal	Why not v1
Extra top-K candidate pruning	RCTD uses all_weights > 0.01; additional pruning can change labels.
True 1D MLE for two-type mixtures	Faster, but not identical to RCTD’s IRWLS/QP procedure.
GPU backend	Harder to keep numerically identical and likely not necessary for the first 10–50×.
Full rewrite of spacexr	Overkill. We only need doublet mode.
Changing platform normalization / sigma estimation	Would change outputs.
Key innovations table
Innovation	Exact-preserving?	ROI	Complexity	v1 priority
New R wrapper run.RCTD.fast.doublet()	Yes	High	Low	P0
Chunked cell processing	Yes	High	Low	P0
Avoid full dense count matrix materialization	Yes, if zero genes still contribute	High	Medium	P0
Prepack RCTD internals into C++ arrays	Yes	High	Medium	P0
C++ implementation of calc_Q_all() / calc_log_l_vec()	Yes	High	Medium	P1
C++ two-type IRWLS scoring	Yes, if exact update path is preserved	Very high	Medium-high	P1
C++ pair-scoring loop	Yes	Very high	Medium	P1
Parallel over cells with RcppParallel/OpenMP	Yes	Very high	Medium	P1
C++ singlet scoring	Yes	Medium	Medium	P2
C++ final selected doublet refit	Yes	Medium	Medium	P2
C++ full all-type fit	Yes, but harder due to QP	Medium-high	High	P3
Optional approximate mode	No	Very high	Medium	Later
Core principle

Exact mode must preserve RCTD’s implementation behavior, not merely the statistical model.

That means preserving:

candidate threshold = all_weights > 0.01
fallback candidate behavior
singlet scoring logic
all pair score matrix semantics
score_mode n.iter = 25
final doublet refit n.iter = 50
MIN_CHANGE behavior
CONFIDENCE_THRESHOLD
DOUBLET_THRESHOLD
spot_class factor levels
first_type / second_type ordering
first_class / second_class behavior
convergence flags
output object schema
2. Core Architecture
2.1 Proposed package/API

Internal package name:

spacexrfast

Main user-facing function:

run.RCTD.fast.doublet(
  rctd,
  max_cores = 8,
  chunk_size = 50000,
  exact = TRUE,
  return_diagnostics = FALSE,
  validate_inputs = TRUE
)

Target usage:

library(spacexr)
library(spacexrfast)

rctd <- create.RCTD(puck, reference, max_cores = 8, test_mode = FALSE)

# Keep upstream RCTD preprocessing:
rctd <- fitBulk(rctd)
rctd <- chooseSigmaC(rctd)

# Replace only the slow doublet pixel/cell fitting:
rctd <- run.RCTD.fast.doublet(
  rctd,
  max_cores = 8,
  chunk_size = 50000,
  exact = TRUE
)

Alternative convenience API:

rctd <- run.RCTD.fast(
  rctd,
  doublet_mode = "doublet",
  max_cores = 8
)
2.2 Architectural philosophy

Keep it simple:

R owns:
  - RCTD object creation
  - gene filtering
  - platform normalization
  - sigma selection
  - likelihood table generation
  - S4 object manipulation
  - output assembly

C++ owns:
  - dense/sparse cell extraction by chunk
  - likelihood interpolation
  - two-type IRWLS
  - singlet score
  - pair score loop
  - final selected doublet refit
  - per-cell classification logic, eventually

No R calls inside parallel workers.
No quadprog calls inside parallel workers.
No GPU in v1.
No approximate pruning in exact mode.

This matters because R is single-threaded at the API level; RcppParallel explicitly warns that worker code should not call the R or Rcpp API because concurrent interaction with R data structures can crash or lead to undefined behavior.

2.3 Execution flow
run.RCTD.fast.doublet(rctd)
    |
    |-- validate RCTD object
    |-- extract immutable RCTD internals
    |-- build FastRCTDContext
    |-- split cells into chunks
    |
    |-- for each chunk:
    |       |
    |       |-- extract counts for chunk
    |       |-- call C++ backend
    |       |       |
    |       |       |-- parallel over cells
    |       |       |-- run exact doublet-mode logic
    |       |       |-- return compact result arrays
    |       |
    |       |-- append results
    |
    |-- assemble results_df
    |-- assemble weights_doublet
    |-- optionally assemble diagnostics
    |-- attach to rctd@results
    |-- return rctd
2.4 Hot path to preserve

Current RCTD doublet mode does this per bead/cell:

1. Build cell_type_profiles[gene_list, ] * UMI_tot
2. Run decompose_full(...)
3. candidates = names(which(all_weights > 0.01))
4. Fallback if candidates length is 0 or 1
5. Compute singlet_scores for candidates
6. For every candidate pair:
       score = decompose_sparse(..., score_mode = TRUE)
7. Pick pair with minimum score
8. Run check_pairs_type(...) for both selected types
9. Decide reject / singlet / doublet_certain / doublet_uncertain
10. Run final decompose_sparse(..., score_mode = FALSE)
11. Return all_weights, spot_class, first_type, second_type, doublet_weights, scores, convergence flags

This is visible in process_bead_doublet(), including the 0.01 candidate threshold, pair loop, classification logic, and final doublet refit.

2.5 Kernel architecture
Kernel A — likelihood interpolation

Port exactly:

calc_Q_all(Y, lambda)
calc_log_l_vec(lambda, Y)
get_d1_d2(Y, lambda)

RCTD clamps Y to K_val, clamps lambda between epsilon = 1e-4 and X_max - epsilon, computes spline indices, fetches from Q_mat / SQ_mat, and performs cubic spline interpolation.

C++ target:

struct LikelihoodTables {
    const double* Q_mat;
    const double* SQ_mat;
    const double* X_vals;
    int K_val;
    int n_x;
    int n_y;
};

inline QDerivatives calc_q_all_one(
    int y,
    double lambda,
    const LikelihoodTables& tables
);
Kernel B — derivatives

Port exactly:

get_der_fast(S, B, gene_list, prediction)

The current function computes d1_vec, d2_vec, gradient, and Hessian using the global S_mat products.

For the two-type case, specialize the Hessian construction without changing math:

S has columns [type_a, type_b]

S_mat has:
  col 0 = S_a * S_a
  col 1 = S_a * S_b
  col 2 = S_b * S_b

hessian:
  h_aa = -sum_g d2_g * S_a_g * S_a_g
  h_ab = -sum_g d2_g * S_a_g * S_b_g
  h_bb = -sum_g d2_g * S_b_g * S_b_g
Kernel C — two-type IRWLS

Port exactly:

solveIRWLS.weights(S, B, nUMI, n.iter = 25 or 50)
solveWLS(S, B, initialSol, nUMI)

RCTD’s solveIRWLS.weights() initializes weights uniformly, clips B > K_val, builds S_mat, iterates until MIN_CHANGE or n.iter, and returns weights plus convergence status.

Each WLS step:

solution = pmax(initialSol, 0)
prediction = abs(S %*% solution)
threshold = max(1e-4, nUMI * 1e-7)
prediction[prediction < threshold] = threshold
derivatives = get_der_fast(...)
d_vec = -grad
D_mat = psd(hess)
norm_factor = norm(D_mat, "2")
D_mat /= norm_factor
d_vec /= norm_factor
D_mat += 1e-7 * I
alpha = 0.3
solution = solution + alpha * QP_step

Those constants and operations are explicit in the R implementation.

For two cell types, replace quadprog::solve.QP with an exact closed-form active-set solver for the same 2D QP subproblem. Do not replace IRWLS with a different optimizer.

Kernel D — pair loop

For each cell:

for candidate i:
    singlet_score[i] = score_singlet(i)

for candidate pair i < j:
    score[i, j] = score_pair(i, j, n_iter = 25)

best_pair = argmin(score)

Parallelization should be across cells, not inside a cell, for simplicity and determinism.

2.6 Parallelism model

Recommended v1:

Parallelization axis: cells
Backend: RcppParallel first, OpenMP optional later
Thread count: max_cores
Scheduling: dynamic chunks of cells
R API in workers: forbidden
Output writes: preallocated per-cell result arrays
Progress: main-thread only, optional

RcppParallel provides high-level parallelFor/parallelReduce abstractions and uses Intel TBB on major platforms, with fallback support elsewhere.

2.7 Exactness contract

Exactness levels:

Level	Definition	Required?
L0	Same spot_class only	Not enough
L1	Same spot_class, first_type, second_type	Minimum
L2	Same labels + weights within tolerance	Required
L3	Same labels + weights + min/singlet scores within tolerance	Required for validation
L4	Bitwise identical floating point	Not realistic

Target:

Categorical outputs: exact equality
Weights: all.equal tolerance <= 1e-10 initially; relax only if needed
Scores: all.equal tolerance <= 1e-8 initially
Convergence flags: exact equality
Tie-breaking: exactly match R candidate order
2.8 Expected speedup logic

A 50× speedup is aggressive but plausible if the pair-scoring loop dominates. With ~20 reference types, the theoretical maximum candidate pair count is:

choose(20, 2) = 190 pair fits / cell

RCTD does not score all 190 pairs unless all 20 pass all_weights > 0.01, but candidate counts of 5–15 still imply 10–105 pair fits per cell. Moving those pair fits from R + quadprog calls to compiled, parallelized C++ is the main opportunity.

Expected staged speedups:

Phase	Expected speedup	Why
Chunking + sparse extraction	1.2–3×	Less memory pressure and less R object churn
C++ likelihood interpolation	2–5×	Removes repeated R vector/matrix overhead
C++ two-type IRWLS	5–20×	Removes repeated quadprog/R optimizer overhead
C++ pair loop + parallel over cells	10–50×	Attacks the dominant doublet-mode cost
C++ full all-type fit	Additional 1.5–5×	Only useful if full fit dominates after pair loop acceleration
3. Data Models
3.1 R-level input model

The R wrapper receives a fitted/prepared RCTD object.

Required fields conceptually:

FastRCTDInput <- list(
  counts = sparse_gene_by_cell_matrix,
  gene_list = character_vector,
  nUMI = numeric_vector_by_cell,
  cell_type_profiles = numeric_matrix_genes_by_celltypes,
  cell_type_names = character_vector,
  class_df = optional_data_frame,
  Q_mat = numeric_matrix,
  SQ_mat = numeric_matrix,
  X_vals = numeric_vector,
  K_val = integer,
  thresholds = list(
    MIN_CHANGE = 0.001,
    CONFIDENCE_THRESHOLD = 10,
    DOUBLET_THRESHOLD = 25
  )
)

The C++ layer should receive explicit copies/views of Q_mat, SQ_mat, X_vals, and K_val. Do not rely on hidden global variables.

3.2 Count matrix model

Target input:

Matrix::dgCMatrix
genes × cells
integer or double counts
filtered to RCTD gene_list before C++ call

Important rule:

Sparse storage is allowed, but zero-count genes cannot be ignored.

RCTD’s likelihood sums over all genes through calc_log_l_vec(lambda, Y), and Y = 0 genes still contribute to the score. The sparse representation is for memory and extraction efficiency, not for dropping zeros.

Implementation choice:

v1 simple path:
  - process one chunk of cells
  - materialize chunk as dense genes × chunk_cells double array
  - run C++ over the chunk
  - free chunk

v2 optimization:
  - sparse-aware likelihood decomposition
  - only if exactness tests pass

For millions of Xenium cells, chunked dense materialization is probably the right simplicity/speed tradeoff. Example:

genes retained: 500–5000
chunk size: 10k–50k cells
dense chunk memory:
  5000 genes × 50000 cells × 8 bytes ≈ 2 GB
  1000 genes × 50000 cells × 8 bytes ≈ 400 MB

So default chunk size should be adaptive.

3.3 Reference profile model

C++ structure:

struct ReferenceProfiles {
    int n_genes;
    int n_types;

    // gene-major: ref[g * n_types + t]
    std::vector<double> ref;

    std::vector<std::string> type_names;

    // Optional class index:
    // class_id[t] = integer class label, or -1 if absent
    std::vector<int> class_id;
};

Reference profiles should be pre-multiplied per cell by UMI_tot exactly as RCTD does:

cell_type_profiles <- cell_type_info[[1]][gene_list, ]
cell_type_profiles <- cell_type_profiles * UMI_tot

That multiplication is currently inside process_bead_doublet().

For performance, avoid rebuilding cell_type_profiles * UMI_tot as a full matrix for every cell. Instead:

lambda_g_for_type_t = UMI_tot_cell * ref_profile[g, t]

This preserves math while avoiding per-cell matrix allocation.

3.4 Likelihood table model
struct LikelihoodTables {
    std::vector<double> Q;      // row-major or column-major, but index must match R
    std::vector<double> SQ;
    std::vector<double> X;
    int K_val;
    double epsilon = 1e-4;
    double delta = 1e-6;
};

RCTD uses cubic spline interpolation from Q_mat, SQ_mat, and X_vals; the index arithmetic must be copied exactly, including the m computation.

3.5 Per-cell model
struct CellInput {
    int cell_index;
    double nUMI;

    // Dense counts over gene_list for this cell.
    // In v1, this is easiest and safest.
    const double* y;  // length n_genes
};
3.6 Candidate set model
struct CandidateSet {
    int n_candidates;
    int candidate_type_ids[MAX_TYPES];  // dynamic vector in practice
};

Candidate generation must match RCTD:

candidates = types where all_weights > 0.01

if no candidates:
    use first min(3, n_types)

if one candidate:
    add first reference type unless candidate is already first;
    otherwise add second reference type

This fallback behavior is explicit in process_bead_doublet().

3.7 Per-cell result model
enum SpotClass {
    REJECT = 0,
    SINGLET = 1,
    DOUBLET_CERTAIN = 2,
    DOUBLET_UNCERTAIN = 3
};

struct CellResult {
    int spot_class;
    int first_type;
    int second_type;

    bool first_class;
    bool second_class;

    double first_weight;
    double second_weight;

    double min_score;
    double singlet_score;

    bool conv_all;
    bool conv_doublet;

    // Optional diagnostics:
    int n_candidates;
    int n_pair_scores;
};
3.8 Batch result model
struct BatchResult {
    std::vector<int> spot_class;
    std::vector<int> first_type;
    std::vector<int> second_type;

    std::vector<uint8_t> first_class;
    std::vector<uint8_t> second_class;

    std::vector<double> first_weight;
    std::vector<double> second_weight;

    std::vector<double> min_score;
    std::vector<double> singlet_score;

    std::vector<uint8_t> conv_all;
    std::vector<uint8_t> conv_doublet;

    // Optional:
    std::vector<double> all_weights;      // n_cells × n_types
    std::vector<double> weights_doublet;  // n_cells × n_types, sparse assembly in R
};

For memory, do not store full score_mat by default. It is useful for debugging but huge at scale.

return_diagnostics = FALSE

should be the production default.

4. Implementation Roadmap: Phased Delivery with ROI Priorities
Phase 0 — Baseline profiler and golden tests

Goal: prove exactly where the 10h runtime goes and establish correctness fixtures before changing anything.

Tasks
Fork or wrap current spacexr.
Add timing instrumentation around:
full all-type fit;
candidate selection;
singlet scoring;
pair scoring loop;
check_pairs_type;
final doublet refit;
output assembly.
Build three benchmark datasets:
tiny: 100 cells × 100 genes × 5 types;
medium: 10k cells × realistic genes × 20 types;
production slice: 100k–500k Xenium cells × 20 types.
Run vanilla RCTD and serialize golden outputs:
results_df;
weights_doublet;
all_weights, if needed;
min_score;
singlet_score;
convergence flags.
Define exactness tests.
Acceptance criteria
expect_equal(fast$results_df$spot_class, vanilla$results_df$spot_class)
expect_equal(fast$results_df$first_type, vanilla$results_df$first_type)
expect_equal(fast$results_df$second_type, vanilla$results_df$second_type)

expect_equal(
  as.matrix(fast$weights_doublet),
  as.matrix(vanilla$weights_doublet),
  tolerance = 1e-10
)
ROI

Very high. Without this, we are guessing.

Deliverable
benchmarks/profile_vanilla_doublet.R
tests/testthat/golden/
docs/profile_report.md
Phase 1 — R wrapper and chunked execution

Goal: create the drop-in API and reduce R memory overhead without touching numerical kernels yet.

Tasks
Implement:
run.RCTD.fast.doublet <- function(
  rctd,
  max_cores = 8,
  chunk_size = NULL,
  exact = TRUE,
  return_diagnostics = FALSE,
  validate_inputs = TRUE
)
Extract required RCTD internals once.
Add adaptive chunk-size selection:
estimate_chunk_size <- function(n_genes, target_gb = 2) {
  floor((target_gb * 1024^3) / (n_genes * 8))
}
Process cells in chunks.
Initially call vanilla RCTD per chunk or per cell to preserve behavior.
Reassemble outputs in original cell order.
Add strict test that chunked vanilla output equals non-chunked vanilla output.
Important design point

Do not optimize too early. First make the drop-in wrapper behave exactly like RCTD.

Acceptance criteria
Same outputs as vanilla RCTD on small and medium fixtures.
No change in algorithm.
Memory peak reduced on large fixtures.
ROI

Medium-high. This creates the production API and prevents all later C++ work from leaking into user-facing R code.

Phase 2 — C++ likelihood kernels

Goal: port the safest hot kernels first.

C++ functions
calc_q_all_cpp()
calc_log_l_vec_cpp()
get_d1_d2_cpp()
score_prediction_cpp()
Tasks
Copy RCTD’s calc_Q_all() logic exactly:
Y[Y > K_val] <- K_val;
epsilon <- 1e-4;
delta <- 1e-6;
clamp lambda;
compute l, m;
fetch from Q_mat, SQ_mat;
cubic spline interpolation.
Implement scalar and vector versions.
Test against R for random:
Y;
lambda;
edge cases near zero;
edge cases near K_val;
edge cases near X_max.
Disable aggressive floating-point optimizations:
no -ffast-math;
avoid changing summation order initially.
Acceptance criteria
all.equal(
  calc_Q_all(Y, lambda),
  calc_q_all_cpp(Y, lambda),
  tolerance = 1e-12
)

all.equal(
  calc_log_l_vec(lambda, Y),
  calc_log_l_vec_cpp(lambda, Y),
  tolerance = 1e-10
)
ROI

High. This removes repeated R vectorization overhead and creates the foundation for the two-type solver.

Phase 3 — C++ exact two-type IRWLS

Goal: replace the repeated two-type decompose_sparse() calls.

Scope

Implement only the two-type case first:

decompose_sparse(
  cell_type_profiles,
  UMI_tot,
  bead,
  type1,
  type2,
  score_mode = TRUE/FALSE,
  constrain = TRUE,
  MIN.CHANGE = 0.001
)

RCTD’s decompose_sparse() uses n.iter = 25 in score mode and n.iter = 50 for final weights.

C++ API
TwoTypeFitResult solve_irwls_two_type_exact(
    const double* y,
    const double* ref_a,
    const double* ref_b,
    double nUMI,
    const LikelihoodTables& tables,
    int n_iter,
    double min_change,
    bool constrain
);
Required behavior

Preserve:

initial weights = [0.5, 0.5]
B clipping to K_val
prediction = abs(S %*% solution)
threshold = max(1e-4, nUMI * 1e-7)
PSD projection of Hessian
spectral norm normalization
epsilon ridge = 1e-7
alpha = 0.3
convergence based on norm(new_solution - solution)

These are present in the current IRWLS/WLS implementation.

QP strategy

For two types, the constrained QP is tiny. Do not call quadprog from workers. Implement an exact 2D active-set solver for the same subproblem.

For constrain = TRUE, RCTD solves for an update u with:

sum(u) = 1 - sum(solution)
u >= -solution

Then:

solution_new = solution + alpha * u

For two types, the equality constraint reduces the update to one scalar. Implement all active-set boundary cases explicitly.

Acceptance criteria

On random fixtures:

for score_mode in c(TRUE, FALSE):
  for constrain in c(TRUE, FALSE):
    compare decompose_sparse_R vs decompose_sparse_cpp

score tolerance <= 1e-8
weights tolerance <= 1e-10
converged exactly equal
ROI

Very high. This is the core acceleration.

Phase 4 — C++ pair-scoring loop

Goal: replace the slowest R-level loop in doublet mode.

Scope

Port this logic:

for(i in 1:(length(candidates)-1)) {
  for(j in (i+1):length(candidates)) {
    score = decompose_sparse(..., score_mode = TRUE)
    score_mat[i,j] = score
    score_mat[j,i] = score

    if(is.null(second_type) || score < min_score) {
      first_type <- type1
      second_type <- type2
      min_score = score
    }
  }
}

This loop is explicit in process_bead_doublet().

C++ API
PairScoreResult score_all_candidate_pairs_exact(
    const double* y,
    const ReferenceProfiles& ref,
    const std::vector<int>& candidates,
    double nUMI,
    const LikelihoodTables& tables,
    double min_change,
    bool constrain,
    bool return_score_mat
);
Rules
Preserve candidate order.
Preserve score < min_score, not <=.
Preserve first minimum tie-breaking.
Compute symmetric score matrix only if diagnostics are requested.
Parallelize across cells, not across pairs, in v1.
Acceptance criteria
Same best pair as vanilla RCTD.
Same min_score within tolerance.
Same score_mat within tolerance when diagnostics are enabled.
ROI

Very high. This should be the first phase where runtime visibly collapses.

Phase 5 — C++ singlet scoring and final doublet refit

Goal: remove remaining repeated R calls inside each cell.

Tasks
Port get_singlet_score().
Port final selected-pair refit:
doublet_results = decompose_sparse(
  ...,
  score_mode = FALSE,
  n.iter = 50
)
Return normalized doublet weights exactly as RCTD does:
results$weights = results$weights / sum(results$weights)

This normalization is in decompose_sparse() when score_mode = FALSE.

Acceptance criteria
Same singlet_score.
Same doublet weights.
Same conv_doublet.
ROI

Medium-high. This finishes the two-type part.

Phase 6 — C++ classification logic

Goal: move all per-cell doublet-mode logic into C++ while keeping S4 assembly in R.

Tasks
Port check_pairs_type().
Preserve class-aware logic.
Preserve class-absent behavior.
Preserve final spot_class factor semantics:
factor(
  spot_class,
  c("reject", "singlet", "doublet_certain", "doublet_uncertain")
)

The factor levels are explicit in the implementation.

Return compact integer-coded outputs to R.
Map integer codes back to R factors and cell type names.
Acceptance criteria
Exact same spot_class.
Exact same first_type / second_type.
Exact same first_class / second_class.
ROI

Medium. This removes R branching overhead and simplifies the C++ batch worker.

Phase 7 — Optional C++ full all-type fit

Goal: accelerate candidate generation if profiling shows it remains dominant.

Why optional

The full all-type fit is harder because it uses a QP with dimension equal to the number of reference cell types. Current solveWLS() calls quadprog::solve.QP() after computing gradients/Hessian and PSD projection.

With ~20 reference types, this is manageable, but reproducing quadprog behavior exactly is more work than the two-type case.

Decision gate

Only start Phase 7 if, after Phases 3–6:

full all-type fit > 30–40% of total runtime
Implementation options
Option	Exactness	Complexity	Recommendation
Keep full fit in R	Exact	Low	Good v1 fallback
Use C++ QP library	Approx/exact depends	Medium-high	Risky
Implement active-set QP	Potentially exact	High	Only if needed
Specialized projected Newton	Not exact	Medium	Not for exact mode
Acceptance criteria
Candidate set exactly matches vanilla RCTD.
all_weights within tolerance.
conv_all exactly matches or is explicitly documented if not.
ROI

Unknown until profiling.

Phase 8 — Production hardening

Goal: make it reliable on millions of Xenium cells.

Tasks
Add adaptive chunk size based on memory.
Add deterministic mode:
deterministic = TRUE
Add optional checkpointing:
checkpoint_dir = NULL
Add failure behavior:
on_cell_error = c("stop", "fallback_to_r", "mark_reject")

Default should be:

on_cell_error = "stop"

Do not silently change outputs in exact mode.

Add telemetry:
cells/sec;
average candidates/cell;
average pair fits/cell;
full fit time;
pair scoring time;
final refit time;
peak memory estimate.
Add production benchmark:
1M cells × realistic Xenium genes × 20 reference cell types
Acceptance criteria
No memory blow-up.
Stable on 8-core machine.
Same outputs as vanilla on validation subsets.
Graceful error messages.
Proposed delivery plan
Milestone 1 — Exact wrapper + profiler

Timeline target: first implementation milestone.

Deliver:

run.RCTD.fast.doublet()
profile_vanilla_doublet()
golden tests
chunked execution

Expected speedup:

1–3×

Main value:

Creates production API and correctness harness.
Milestone 2 — C++ likelihood + two-type solver

Deliver:

calc_Q_all_cpp
calc_log_l_vec_cpp
get_der_fast_two_type_cpp
solve_irwls_two_type_exact_cpp
decompose_sparse_cpp

Expected speedup:

5–20× on pair-scoring-heavy cells

Main value:

Removes repeated R/quadprog overhead for pair scoring.
Milestone 3 — C++ pair loop

Deliver:

score_all_candidate_pairs_exact_cpp
best-pair selection
optional score_mat diagnostics

Expected speedup:

10–50× depending on candidate counts

Main value:

Attacks the dominant doublet-mode loop.
Milestone 4 — Full per-cell C++ doublet path

Deliver:

singlet scoring
final doublet refit
check_pairs_type
classification logic
compact result arrays
R output assembly

Expected speedup:

20–50× realistic target if pair loop dominates

Main value:

Production-ready exact doublet backend.
Milestone 5 — Full all-type fit acceleration, only if needed

Deliver:

C++ full all-type IRWLS/QP or optimized fallback strategy

Expected speedup:

additional 1.5–5× if full fit remains bottleneck

Main value:

Closes remaining gap if 50× is not reached.
Recommended implementation order

Do not start with GPU. Do not start with a full C++ rewrite. Do not start with a new optimizer.

Start here:

1. Build golden tests.
2. Build R wrapper.
3. Port calc_Q_all / calc_log_l_vec.
4. Port exact two-type IRWLS.
5. Replace pair-scoring loop.
6. Validate exactness on thousands of cells.
7. Scale to 100k cells.
8. Scale to millions of cells.
9. Only then decide whether full all-type fit needs C++.

This matches the philosophy: simple, fast, low-overengineering, easy low-hanging fruit first.

Risk register
Risk	Severity	Mitigation
Exactness breaks due to floating-point differences	High	Golden tests, deterministic summation, no fast-math, preserve candidate order
Two-type QP differs from quadprog	High	Unit-test QP step against R on randomized cases; fallback to R during dev
Full all-type fit dominates after pair loop optimization	Medium	Profile after Milestone 3; only then port full fit
Sparse optimization accidentally ignores zero genes	High	v1 uses dense chunk over genes; sparse-aware mode later
R API called from worker threads	High	Use RcppParallel accessors and plain C++ buffers only
Diagnostics too large for millions of cells	Medium	return_diagnostics = FALSE default
GPL-3.0 implications	Medium	Since spacexr is GPL-3.0, be careful if distributing modified/derived code outside internal use.
50× not reached on 8 cores	Medium	Measure candidate counts; consider full fit port or approximate mode later
Final target spec
rctd_fast <- run.RCTD.fast.doublet(
  rctd,
  max_cores = 8,
  chunk_size = "auto",
  exact = TRUE,
  return_diagnostics = FALSE
)

Must satisfy:

Same spot_class as RCTD.
Same first_type / second_type as RCTD.
Same weights_doublet within tight numerical tolerance.
Same confidence logic.
Same factor levels.
Same output schema.
Materially faster on million-cell Xenium slides.
Target: 50× wall-clock speedup from 10h to ~10–15 minutes if pair scoring dominates.

The hard truth: 50× with exact identity is possible only if the current runtime is dominated by repeated two-type pair scoring and R/quadprog overhead. That is likely, but Phase 0 profiling should confirm it before we touch the full all-type solver.

