# Set example values for A and sigma_eps.
sigma_eps =matrix(c(58282.95, 58876.84, 58876.84, 66030.17), nrow = 2, ncol = 2)
A = array(c(-0.9002804, 1.300596, -0.05519068, -1831.0805486, 5537.174244, 16933.09 ), dim = c(2, 2, 2))


# Set the forecast horizon (H).
H = dim(A)[3]

# When H=1, the above line returns H=na because R thinks there is no third dimension.
# Thus, in that case, we must force H=1.
if (is.na(H) == TRUE) {
  H = 1
}

# Set the number of variables in the VAR (K).
K = dim(sigma_eps)[1]

# Calculate the forecast error covariance matrix (Xi).
if (H > 1) {
  Xi_h = array(0, dim = c(K, K, H))
  for (h in 1:H) {
    Xi_h[,,h] = A[,,h] %*% sigma_eps %*% t(A[,,h])  # Calculated Xi at each h.
  }

  # Sum them along THIRD dimension to form Xi.
  Xi = rowSums(Xi_h, dims = 2)  # Note: dims = 2 sums along the third dimension.
} else {
  Xi = sigma_eps  # When H = 1, Xi is equal to sigma_eps because A_0 = I_K.
}

# Identity matrix.
I_K = diag(1, nrow = K, ncol = K)
# Calculate the gFEVD.
if (H > 1) {
  A_times_sigma_squar_h = array(0, dim = c(K, K, H))
  for (h in 1:H) {
    # Calculate A times Sigma squared (element by element squared) at each h.
    A_times_sigma_squar_h[,,h] = (A[,,h] %*% sigma_eps) * (A[,,h] %*% sigma_eps)
  }

  # Sum them along THIRD dimension
  A_times_sigma_squar_sum = rowSums(A_times_sigma_squar_h, dims = 2)
  # Note that because the above function is a row sum, dims=2 actually sums along...
  # the third dimension.
} else {
  # When H=1, the numerator of the gFEVD is just the element by element squares...
  # of sigma_eps.
  A_times_sigma_squar_sum = sigma_eps * sigma_eps
}

gFEVD = array(0, dim = c(K, K))
for (i in 1:K) {
  for (j in 1:K) {
    gFEVD[i,j] = A_times_sigma_squar_sum[i,j] / (sigma_eps[j,j] * Xi[i,i])
  }
}

# Calculate the gSOT.
gFEVD_row_sum = rowSums(gFEVD)  # Calculate the gFEVD row sums.
gSOT = array(0, dim = c(K, K))
for (i in 1:K) {
  for (j in 1:K) {
    gSOT[i,j] = (gFEVD[i,j]) / gFEVD_row_sum[i]  # Calculate the gSOT.
  }
}

# Calculate the generalized total spillover from all others to variable i...
# (S_gen_from_others) and the generalized total spillover to others from variable...
# i (S_gen_to_others).
S_gen_from_others = array(0, dim = c(K))
S_gen_to_others = array(0, dim = c(K))
for (i in 1:K) {
  # Generalized total spillover from others to variable i.
  S_gen_from_others[i] = sum(gSOT[i,]) - gSOT[i,i]
  # Generalized total spillover to others from variable i.
  S_gen_to_others[i] = sum(gSOT[,i]) - gSOT[i,i]
}

# Calculate the generalized net total spillover (S_gen_net).
S_gen_net = S_gen_to_others - S_gen_from_others
# Calculate the generalized spillover index (gSOI).
gSOI = mean(S_gen_from_others)

# Calculate the elimination matrices. These are usually denoted as a KxK-1 matrix...
# M_i. Here, they are an array where M[,,1]=M_1, and in general M[,,i]=M_i.
M = array(0, dim = c(K, K-1, K))
for (i in 1:K) {
  M[,,i] = I_K[,-i]  # Calculate the elimination matrices.
}
#Calculate    the    joint    total    spillover    from    all    others    to    variable    i...
#(S_jnt_from_others).
#Calculate    the    numerator    of    S_jnt_from_others.
# Calculate the numerator of S_jnt_from_others.
if (H > 1) {
  S_jnt_from_numerator_h = array(0, dim = c(K, H))
  for (i in 1:K) {
    for (h in 1:H) {
      # Calculate the numerator of S_jnt_from_others at each h.
      S_jnt_from_numerator_h[i, h] = I_K[i, ] %*% A[,,h] %*% sigma_eps %*% M[,,i] %*%
                                       (solve(t(M[,,i]) %*% sigma_eps %*% M[,,i])) %*%
                                       t(M[,,i]) %*% sigma_eps %*% t(A[,,h]) %*% I_K[,i]
    }
  }
}

S_jnt_from_numerator = array(0, dim = c(K))
for (i in 1:K) {
  # Calculate the numerator of S_jnt_from_others (sum over h).
  S_jnt_from_numerator[i] = sum(S_jnt_from_numerator_h[i, ])
}

if (H == 1) {
  S_jnt_from_numerator = array(0, dim = c(K))
  for (i in 1:K) {
    S_jnt_from_numerator[i] = I_K[i, ] %*% sigma_eps %*% M[,,i] %*%
                               (solve(t(M[,,i]) %*% sigma_eps %*% M[,,i])) %*%
                               t(M[,,i]) %*% sigma_eps %*% I_K[,i]
    # When H=1, A is just the identity matrix, and thus it cancels out.
  }
}

S_jnt_from_others = array(0, dim = c(K))
for (i in 1:K) {
  S_jnt_from_others[i] = S_jnt_from_numerator[i] / Xi[i, i]
}

# Calculate the joint spillover index (jSOI).
jSOI = mean(S_jnt_from_others)

# Calculate the scaling factor lambda.
lambda = jSOI / gSOI

# Calculate gSOT tilde (gSOT_tilde).
gSOT_tilde = lambda * gSOT

# Calculate the joint total spillover to all others from variable i...
# (S_jnt_to_others).
S_jnt_to_others = array(0, dim = c(K))
for (i in 1:K) {
  # Joint total spillover to others from variable i.
  S_jnt_to_others[i] = sum(gSOT_tilde[,i]) - gSOT_tilde[i,i]
}

# Calculate the joint net total spillover (S_jnt_net).
S_jnt_net = S_jnt_to_others - S_jnt_from_others

##############################################################################
#OUTPUT.
#Print    the    spillover    measures.         All    spillover    measures    (excluding    lambda)    have... #been    multiplied    by    100    and    thus    are    reported    as    a    percent.

#Print    the    generalized    spillover    table.
print(100*gSOT)
#Print    the    generalized    total    spillover    from    all    others    to    the    individual    variables. print(100*S_gen_from_others)
#Print    the    generalized    total    spillover    to    all    others    from    the    individual    variables. print(100*S_gen_to_others)
#Print    the    generalized    net    total    spillover.
print(100*S_gen_net)
#Print    the    generalized    spillover    index.
print(100*gSOI)
#Print    the    joint    total    spillover    from    all    others    to    the    individual    variables. print(100*S_jnt_from_others)
#Print    the    joint    total    spillover    to    all    others    from    the    individual    variables. print(100*S_jnt_to_others)
#Print    the    joint    net    total    spillover. 
print(100*S_jnt_net)
#Print    the    scaling    factor. 
print(lambda)
#Print    the    joint    spillover    index. 
print(100*jSOI)
