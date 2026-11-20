[strong] 464 slides, 400 patients, prevalence 0.32, D = 64 
[strong] pca    estimate 0.0393  null 0.0014 +/- 0.0029  p 0.010  k* 4 
[strong] cosine estimate 0.0487  null 0.0014 +/- 0.0024  p 0.010  k* 5 
[strong] rbf    estimate 0.0476  null 0.0016 +/- 0.0025  p 0.010  k* 5 
[strong] poly   estimate 0.0516  null 0.0017 +/- 0.0029  p 0.010  k* 6 
[strong] max-T selects poly: 0.0516, p_maxT = 0.010 
[weak] 461 slides, 400 patients, prevalence 0.28, D = 64 
[weak] pca    estimate 0.0000  null 0.0015 +/- 0.0024  p 1.000  k* 1 
[weak] cosine estimate 0.0082  null 0.0014 +/- 0.0023  p 0.040  k* 1 
[weak] rbf    estimate 0.0023  null 0.0015 +/- 0.0024  p 0.228  k* 1 
[weak] poly   estimate 0.0001  null 0.0013 +/- 0.0021  p 0.396  k* 1 
[weak] max-T selects cosine: 0.0082, p_maxT = 0.059 
eps_rep: replicate max 0.0000, sign-flip max 0.0000, null-block max 0.0040, bound 0.0040 
[strong] half-sampling s.e. 0.0236, 90% interval 0.0128-0.0903 
[strong] usable ridge 0.0445 (AUC 0.68), usable PathMI readout 0.0303 (AUC 0.64) 
[strong] attributable: covariates 0.0777, combined 0.1382, incremental 0.0605 
[strong] N_evidence 400 
[weak] half-sampling s.e. 0.0092, 90% interval 0.0000-0.0233 
[weak] usable ridge -0.0020 (AUC 0.51), usable PathMI readout -0.0082 (AUC 0.49) 
[weak] attributable: covariates 0.0615, combined 0.0721, incremental 0.0106 
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
  null           0.0017 +/- 0.0029 (p = 0.010) | eps_rep 0.0040
  90% interval   0.0128 to 0.0903 (s.e. 0.0236)
  attributable   0.0605
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
  null           0.0014 +/- 0.0023 (p = 0.040) | eps_rep 0.0040
  90% interval   0.0000 to 0.0233 (s.e. 0.0092)
  attributable   0.0106
  usable         -0.0020 (AUC 0.51)
  status         unresolved (N_evidence not reached)
hidden context pi_A 1.00: estimate 0.1713, I(Z;Y) 0.1689, pi_A*I_A 0.1689, I_A 0.1689 
hidden context pi_A 0.75: estimate 0.0776, I(Z;Y) 0.0890, pi_A*I_A 0.1267, I_A 0.1689 
hidden context pi_A 0.50: estimate 0.0217, I(Z;Y) 0.0384, pi_A*I_A 0.0845, I_A 0.1689 
hidden context pi_A 0.25: estimate 0.0028, I(Z;Y) 0.0095, pi_A*I_A 0.0422, I_A 0.1689 
