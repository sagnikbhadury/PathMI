u <- read.csv("C:/Users/bhadury/University of Michigan Dropbox/Sagnik Bhadury/CODES/SBTurboCodes/entropy/Research - DOE Horseshoe/TRIDENT/usable/usable_infold.csv")
u <- u[!u$control, ]; d <- u$MI_refit - u$MI_obs
cat(sprintf("cells %d | exact (<1e-9) %d | |diff| > 0.001: %d | > 0.0095: %d | median |diff| %.4f | max %.4f\n", nrow(u), sum(abs(d) < 1e-9), sum(abs(d) > 1e-3), sum(abs(d) > 0.0095), median(abs(d)), max(abs(d))))
print(table(u$decomp, abs(d) < 1e-9))
print(head(u[order(-abs(d)), c("id", "MI_obs", "MI_refit")], 8), digits = 4)
