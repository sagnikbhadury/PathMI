[strong] 464 slides, 400 patients, prevalence 0.32, D = 64 
[strong] pca_mean         estimate 0.0393  null 0.0014 +/- 0.0029  p 0.010  k* 4 
[strong] cosine_mean      estimate 0.0487  null 0.0014 +/- 0.0024  p 0.010  k* 5 
[strong] rbf_mean         estimate 0.0476  null 0.0016 +/- 0.0025  p 0.010  k* 5 
[strong] poly_mean        estimate 0.0516  null 0.0017 +/- 0.0029  p 0.010  k* 6 
[strong] pca_norm         estimate 0.0304  null 0.0014 +/- 0.0030  p 0.010  k* 3 
[strong] cosine_norm      estimate 0.0434  null 0.0015 +/- 0.0025  p 0.010  k* 5 
[strong] rbf_norm         estimate 0.0414  null 0.0015 +/- 0.0024  p 0.010  k* 5 
[strong] poly_norm        estimate 0.0450  null 0.0015 +/- 0.0025  p 0.010  k* 7 
[strong] pca_moment       estimate 0.5810  null 0.0012 +/- 0.0021  p 0.010  k* 1 
[strong] cosine_moment    estimate 0.5746  null 0.0011 +/- 0.0017  p 0.010  k* 1 
[strong] rbf_moment       estimate 0.5778  null 0.0014 +/- 0.0021  p 0.010  k* 1 
[strong] poly_moment      estimate 0.5811  null 0.0016 +/- 0.0026  p 0.010  k* 1 
[strong] pca_attention    estimate 0.0084  null 0.0013 +/- 0.0025  p 0.050  k* 1 
[strong] cosine_attention estimate 0.0038  null 0.0017 +/- 0.0026  p 0.149  k* 1 
[strong] rbf_attention    estimate 0.0045  null 0.0014 +/- 0.0027  p 0.099  k* 1 
[strong] poly_attention   estimate 0.0160  null 0.0013 +/- 0.0026  p 0.020  k* 2 
[strong] mean pooling: max-T over 4 projections selects poly_mean (0.0516), p_maxT = 0.010 
[strong] norm pooling: max-T over 4 projections selects poly_norm (0.0450), p_maxT = 0.010 
[strong] moment pooling: max-T over 4 projections selects poly_moment (0.5811), p_maxT = 0.010 
[strong] attention pooling: max-T over 4 projections selects poly_attention (0.0160), p_maxT = 0.030 
[strong] max over 16 configurations: poly_moment (0.5811), p_max16 = 0.010 | alignment refits max |diff| = 0 
[weak] 461 slides, 400 patients, prevalence 0.28, D = 64 
[weak] pca_mean         estimate 0.0000  null 0.0015 +/- 0.0024  p 1.000  k* 1 
[weak] cosine_mean      estimate 0.0082  null 0.0014 +/- 0.0023  p 0.040  k* 1 
[weak] rbf_mean         estimate 0.0023  null 0.0015 +/- 0.0024  p 0.228  k* 1 
[weak] poly_mean        estimate 0.0001  null 0.0013 +/- 0.0021  p 0.396  k* 1 
[weak] pca_norm         estimate 0.0000  null 0.0014 +/- 0.0023  p 1.000  k* 1 
[weak] cosine_norm      estimate 0.0079  null 0.0016 +/- 0.0025  p 0.069  k* 1 
[weak] rbf_norm         estimate 0.0015  null 0.0015 +/- 0.0023  p 0.287  k* 1 
[weak] poly_norm        estimate 0.0003  null 0.0011 +/- 0.0018  p 0.446  k* 1 
[weak] pca_moment       estimate 0.2281  null 0.0018 +/- 0.0034  p 0.010  k* 7 
[weak] cosine_moment    estimate 0.2272  null 0.0016 +/- 0.0023  p 0.010  k* 9 
[weak] rbf_moment       estimate 0.2324  null 0.0017 +/- 0.0027  p 0.010  k* 6 
[weak] poly_moment      estimate 0.2039  null 0.0015 +/- 0.0023  p 0.010  k* 9 
[weak] pca_attention    estimate 0.0000  null 0.0016 +/- 0.0024  p 1.000  k* 1 
[weak] cosine_attention estimate 0.0030  null 0.0013 +/- 0.0021  p 0.149  k* 1 
[weak] rbf_attention    estimate 0.0010  null 0.0015 +/- 0.0024  p 0.366  k* 1 
[weak] poly_attention   estimate 0.0000  null 0.0016 +/- 0.0030  p 1.000  k* 1 
[weak] mean pooling: max-T over 4 projections selects cosine_mean (0.0082), p_maxT = 0.059 
[weak] norm pooling: max-T over 4 projections selects cosine_norm (0.0079), p_maxT = 0.079 
[weak] moment pooling: max-T over 4 projections selects rbf_moment (0.2324), p_maxT = 0.010 
[weak] attention pooling: max-T over 4 projections selects cosine_attention (0.0030), p_maxT = 0.366 
[weak] max over 16 configurations: rbf_moment (0.2324), p_max16 = 0.010 | alignment refits max |diff| = 0 
numerical tolerance: replicate max 0.0000, sign-flip max 0.0000 -> eps_tol 0.0000 | specification (null block) max 0.0040, 5 of 8 exact 
[strong] recorded configuration poly_mean | half-sampling s.e. 0.0236, 90% interval 0.0128-0.0903 
[strong] usable ridge 0.0445 (AUC 0.68), usable PathMI readout 0.0303 (AUC 0.64) 
[strong] incremental beyond covariates: covariates 0.0777, combined 0.1382, incremental 0.0605 
[strong] N_evidence 400 
[weak] recorded configuration cosine_mean | half-sampling s.e. 0.0092, 90% interval 0.0000-0.0233 
[weak] usable ridge -0.0020 (AUC 0.51), usable PathMI readout -0.0082 (AUC 0.49) 
[weak] incremental beyond covariates: covariates 0.0615, combined 0.0721, incremental 0.0106 
[weak] N_evidence NA 
tile identity max error 1.33e-15; linear arm estimate 0.0393, k* 4 
burden N90 15, N95 30 tiles per slide 
Recoverability Record
  cohort         simulated strong
  model          simulated encoder (D = 64)
  pooling        mean
  projection     poly
  d              20
  unit           patient
  represented    0.0516 nats (8% of H(Y) = 0.631), k* = 6, n = 400
  null           0.0017 +/- 0.0029 (p = 0.010; max over 4 configs p = 0.010) | eps_tol 0.0000
  90% interval   0.0128 to 0.0903 (s.e. 0.0236)
  incremental    0.0605 (beyond covariates)
  usable         0.0445 (AUC 0.68)
  tiles N90/N95  15 / 30
  status         resolved (N_evidence 400)
Recoverability Record
  cohort         simulated weak
  model          simulated encoder (D = 64)
  pooling        mean
  projection     cosine
  d              20
  unit           patient
  represented    0.0082 nats (1% of H(Y) = 0.588), k* = 1, n = 400
  null           0.0014 +/- 0.0023 (p = 0.040; max over 4 configs p = 0.059) | eps_tol 0.0000
  90% interval   0.0000 to 0.0233 (s.e. 0.0092)
  incremental    0.0106 (beyond covariates)
  usable         -0.0020 (AUC 0.51)
  status         unresolved (N_evidence not reached)
hidden context pi_A 1.00: estimate 0.1713, I(Z;Y) 0.1689, pi_A*I_A 0.1689, I_A 0.1689 
hidden context pi_A 0.75: estimate 0.0776, I(Z;Y) 0.0890, pi_A*I_A 0.1267, I_A 0.1689 
hidden context pi_A 0.50: estimate 0.0217, I(Z;Y) 0.0384, pi_A*I_A 0.0845, I_A 0.1689 
hidden context pi_A 0.25: estimate 0.0028, I(Z;Y) 0.0095, pi_A*I_A 0.0422, I_A 0.1689 
