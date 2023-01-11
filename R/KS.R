KS_dot_vec <-
function (A, B, C, inv = FALSE, svd.A = svd(A), svd.C = svd(C),diag.A=FALSE) 
{
    W <- svd.A$d %o% rep(1, ncol(B)) + rep(1, nrow(B)) %o% svd.C$d
    if(inv)
        W <- 1/W
    if( ! diag.A ){
      cur <- t(svd.A$u) %*% B %*% svd.C$u
      cur <- W * cur
      cur <- svd.A$u %*% cur %*% t(svd.C$u)
    } else {
      cur <- B %*% svd.C$u
      cur <- W * cur
      cur <- cur %*% t(svd.C$u)
    }
    return(cur)
}

trp_KS <-
function (keep, kill, svd.keep = svd(keep), lam.kill = svd(kill)$d) 
{
    Lambda <- diag(sapply(svd.keep$d, function(lam_x) sum(1/(lam_x + 
        lam.kill))))
    return(svd.keep$u %*% Lambda %*% t(svd.keep$u))
}

tr <-
function (X){
  sum(diag(X))
}
