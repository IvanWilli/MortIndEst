# function to split child mortality indicators by single age.
# Needs two of the three arguments q0, q1_4 and q0_5.
# Method: use lx at ages 0, 1 and 5 to interpolate using monotonic spline
split_q0_5 <- function(q0 = NULL, q1_4 = NULL, q0_5 = NULL){
  if(!is.null(q0) & !is.null(q1_4)){
    lx <- c(1, 1-q0, (1-q0)*(1-q1_4))
  } 
  if(!is.null(q0) & !is.null(q0_5) & is.null(q1_4)){
    lx <- c(1, 1-q0, (1-q0_5))
  } 
  if(!is.null(q1_4) & !is.null(q0_5) & is.null(q0)){
    lx <- c(1, (1-q0_5)/(1-q1_4), (1-q0_5))
  }
  x <- c(0, 1, 5)
  lx_interp <- splinefun(x = x, y = lx, method = "monoH.FC")(0:5)
  qx_interp <- -diff(lx_interp) / lx_interp[-length(lx_interp)]
  return(data.frame(x = 0:4, lx = lx_interp[1:5], qx = qx_interp))
}

# example: using two Interchangeably gives the same result
# q0 <- .05
# q1_4 <- .005
# q0_5 <- 1-(1-q0)*(1-q1_4)
# split_q0_5(q0 = q0, q1_4 = q1_4)
# split_q0_5(q0 = q0, q0_5 = q0_5)
# split_q0_5(q0 = q0, q1_4 = q1_4)
