pkgload::load_all(".", quiet = TRUE)
ns <- asNamespace("baystability")
data(cultivo2008); data(cultivo2009)
pub <- read.csv("inst/extdata/perez2012.csv"); pl <- pub$mean[grepl("lambda",pub$parameter)]
ex <- list(y = quote(y), gen = quote(entry), env = quote(site), rep = quote(rep))
cd <- ns$.bys_prep("t", cultivo2009, ex$y, ex$gen, ex$env, ex$rep, environment(), allow_missing = TRUE)
K <- 11L
ls <- ns$.bys_ls(cd$Y, K)
tp <- ns$.bys_trial_prior("t", cultivo2008, ex, environment(), cd, K, FALSE, NULL)
pr <- bys_priors(); pr$nu <- tp$a; pr$s2 <- tp$s2
pl0 <- ns$.bys_prior_list(pr, 12, 25, K, tp)
st0 <- utils::modifyList(ls, tp$fit[c("mu","tau","delta","lambda","alpha","gamma")])
go <- function(st, pairs, n = 1500, b = 500) {
  if (!pairs) assignInNamespace(".bys_pair_sweep", function(X, C, one) X, "baystability")
  set.seed(7)
  r <- ns$.bys_chain(cd$Y, 3, cd$ssw, K, pl0, st, list(alpha = tp$fit$alpha), n, (b+1):n)
  colMeans(r$scalar[, 2 + 12 + 25 + 1:K])
}
s1 <- ns$.bys_start(st0, tp$s2, 1)
rnd <- s1; set.seed(99)
rnd$alpha <- ns$.bys_orth1(matrix(rnorm(12*K),12)); rnd$gamma <- ns$.bys_orth1(matrix(rnorm(25*K),25))
a <- go(s1, TRUE); b <- go(s1, FALSE); c <- go(rnd, FALSE)
print(round(cbind(pub = pl, half_l0 = tp$fit$lambda/2, pairs = a, nopairs_lsstart = b, nopairs_random = c), 0))
