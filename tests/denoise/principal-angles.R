principalAngles <- function(A, B) {
  stopifnot(all(dim(A) == dim(B)))
  Qa = qr.Q(qr(A))
  Qb = qr.Q(qr(B))
  C = svd(crossprod(Qa, Qb))$d
  C = vapply(C, function(x) min(x, 1), numeric(1))
  angles = acos(C)
  return(angles)
}
