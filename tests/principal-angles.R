principalAngles <- function(A, B) {
  if (is.vector(A)) A = matrix(A, ncol=1)
  if (is.vector(B)) B = matrix(B, ncol=1)

  Qa = qr.Q(qr(A))
  Qb = qr.Q(qr(B))

  C = svd(crossprod(Qa, Qb))$d
  C = vapply(C, function(x) min(x, 1), numeric(1))
  angles = sort(acos(C), decreasing=TRUE)
  return(angles)
}
