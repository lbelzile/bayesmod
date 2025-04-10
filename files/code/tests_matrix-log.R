S <- rWishart(
  n = 1,
  df = 10,
  Sigma = diag(5) + matrix(0.5, 5, 5)
)[,, 1]
# Singular value decomposition
svdS <- La.svd(S)
Sb <- with(svdS, u %*% diag(d) %*% vt)
isTRUE(all.equal(S, Sb))
# Matrix log of matrix S
L <- with(svdS, u %*% diag(0.5 * log(d)) %*% vt)
expL <- with(svdS, u %*% diag(sqrt(d)) %*% vt)
# Construct matrix exponential directly from L
isTRUE(
  all.equal(
    with(La.svd(L), u %*% diag(exp(d)) %*% vt),
    expL
  )
)
isTRUE(all.equal(
  with(La.svd(L), u %*% diag(exp(2 * d)) %*% vt),
  crossprod(expL)
))
isSymmetric(expL)
isTRUE(all.equal(crossprod(expL), S))
# Alternative formulations for log absolute value of determinant
isTRUE(all.equal(
  0.5 * determinant(S)$modulus,
  sum(diag(L)),
  check.attributes = FALSE
))
