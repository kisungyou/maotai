# auxiliary functions to deal with ECDF objects
# (1) elist_check   : list of 'ecdf' objects
# (2) elist_fform   : make a function form in a discrete grid
# (3) elist_epmeans : either a vector or something


# (1) elist_check ---------------------------------------------------------
#' @keywords internal
#' @noRd
elist_check <- function(elist){
  cond1 = (is.list(elist))
  cond2 = all(unlist(lapply(elist, inherits, "ecdf"))==TRUE)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# (2) elist_fform ---------------------------------------------------------
#' @keywords internal
#' @noRd
elist_fform <- function(elist){
  nlist = length(elist)
  # compute knot points
  allknots = array(0,c(nlist,2))
  for (i in 1:nlist){
    tgt = stats::knots(elist[[i]])
    allknots[i,] = c(min(tgt), max(tgt))
  }
  mint = min(allknots[,1]) - 0.01
  maxt = max(allknots[,2]) + 0.01
  ssize = min((maxt-mint)/1000, 0.001)
  tseq  = seq(mint, maxt, by=ssize)
  # return the list of y values
  outY = list()
  for (i in 1:nlist){
    tgt       = elist[[i]]
    outY[[i]] = tgt(tseq)
  }
  # return the result
  output = list()
  output$tseq = tseq
  output$fval = outY # list of function values
  return(output)
}

# (3) elist_epmeans -------------------------------------------------------
#' @keywords internal
#' @noRd
elist_epmeans <- function(elist){
  N      = length(elist)
  output = list()
  for (n in 1:N){
    tgt = elist[[n]]
    if (is.vector(tgt)&&(!any(is.infinite(tgt)))&&(!any(is.na(tgt)))){ # Case 1. just a vector
      output[[n]] = stats::ecdf(tgt)
    } else if (inherits(tgt, "ecdf")){
      output[[n]] = tgt
    } else {
      smsg = paste("* epmeans : ",n,"-th element from 'elist' is neither an 'ecdf' object nor a vector.")
      stop(smsg)
    }
  }
  return(output)
}