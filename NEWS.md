# PathMI 0.3.0

* `numerical_tolerance()` replaces `reproducibility_bound()` as the
  tolerance of the protocol: the largest discrepancy over the numerical
  identities (replicate refit, sign-flip refit). A refit after adding a
  null covariate block is a specification check, passed as
  `specification` and reported outside the tolerance.
  `reproducibility_bound()` remains as a deprecated wrapper that treats a
  `null_block` argument this way.
* `select_representation()` now works over any set of candidates (the four
  projections of a pooling rule, or all sixteen pooling x projection
  configurations for the paper's p_max16), stops when the candidates were not
  refitted on the same permuted labels, and returns the null maxima.
* `pathmi_null()` stores its permutation matrix and fit settings.
* New `verify_permutation_alignment()` refits stored permuted fits and checks
  that they reproduce exactly.
* `recoverability_record()` (schema 1.1): `eps_tol` and `incremental`
  replace `eps_rep` and `attributable` (old names still accepted); new
  `selection` argument stores the configuration-selection p value; records
  carry protocol notes, the resolution margin, and the fields
  `supported_interpretations` and `interpretations_not_established`.
* `learning_curve()` and `evidence_status()` take `eps_tol`; documentation
  states that the resolution criterion is operational and that the
  learning-curve and jackknife scalings describe studies of different sizes.
* README: worked example over 16 configurations and a guide that reproduces
  every analysis of the manuscript; `results/` and `study/` synchronised with
  the submitted version (max-16 tests, tolerance, status re-derivation).

# PathMI 0.2.0

* Full measurement toolkit, patient-level results and simulated README.
