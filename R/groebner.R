##' Create  monomial term
##'
##' A monomial term is represented as a coefficient and a vector of
##' variable exponents. This function is used in conjuction with
##' [poly()] to create a polynomial.
##'
##' @title Create a monomial
##' @param coef the coefficient of the term
##' @param expt vector of variable exponents
##' @seealso [poly()]
##' @return a pterm object
##' @export
pterm <- function(coef, expt) {
  r <- list(coef=coef, expt=as.integer(expt))
  class(r) <- "pterm"
  r
}

# As character method for pterm
##' @export
as.character.pterm <- function(x,...) {
  keep <- which(x$expt != 0L)
  if(length(keep)) {
    tm <- paste0("x",keep,ifelse(x$expt[keep]!=1L,paste0("^",x$expt[keep]),""),collapse="*")
    if(abs(x$coef) != 1) paste0(x$coef,tm) else (if(x$coef < 0) paste0("-",tm) else tm)
  } else {
    as.character(x$coef)
  }
}


# Print method for pterm
##' @export
print.pterm <- function(x, ...) {
  cat(as.character(x),"\n")
}


##' Polynomial term orders.
##'
##' Three standard term ordering are available
##'
##' * `lex` lexicographic
##' * `grlex` graded lexicographic
##' * `grevlex` reverse graded lexicographic
##'
##' but any arbitrary term order can be used. The default term order
##' is `grevlex`.
##'
##' The current term can be set and retrieved with the set_term_order
##' and get_term_order functions.
##'
##' If the term order is changed, any existing polynomials are left in
##' an inconsistent state, and it is necessary to reapply [poly()] to
##' recover a valid polynomial.
##'
##' @title Term Orders
##' @param expt1 vector of variable exponents
##' @param expt2 vector of variable exponents
##' @param order a term order function
##' @seealso [reorder_vars()]
##' @return The term order functions return
##'   * -1 if `expt1 < expt2`
##'   * 0 if `expt1 == expt2`
##'   * 1 if `expt1 > expt2`
##' @rdname term_order
##' @export
lex <- function(expt1,expt2) {
  for(i in seq_along(expt1)) {
    if(expt1[i] != expt2[i]) return(if(expt1[i] < expt2[i]) -1L else 1L)
  }
  0L
}

##' @rdname term_order
##' @export
grlex <- function(expt1,expt2) {
  s1 <- sum(expt1)
  s2 <- sum(expt2)
  if(s1 != s2) return(if(s1 < s2) -1L else 1L)
  for(i in seq_along(expt1)) {
    if(expt1[i] != expt2[i]) return(if(expt1[i] < expt2[i]) -1L else 1L)
  }
  0L
}

##' @rdname term_order
##' @export
grevlex <- function(expt1,expt2) {
  s1 <- sum(expt1)
  s2 <- sum(expt2)
  if(s1 != s2) return(if(s1 < s2) -1L else 1L)
  for(i in rev(seq_along(expt1))) {
    if(expt1[i] != expt2[i]) return(if(expt1[i] < expt2[i]) -1L else 1L)
  }
  0L
}


term_order <- new.env(parent = emptyenv())
term_order$order <- grevlex

##' @rdname term_order
##' @export
set_term_order <- function(order = grevlex) {
  term_order$order <- order
}

##' @rdname term_order
##' @export
get_term_order <- function(order = grevlex) {
  term_order$order
}


##' Create a polynomial
##'
##' A polynomial is represented as a list of monomials, sorted
##' according to the term order, where no two terms have the same
##' exponents and no term has coefficient zero.
##'
##' If the term order is changed, any existing polynomials are left in
##' an inconsistent state and [poly()] should be applied to recover a
##' valid polynomial.
##'
##' `poly` reduces the list of terms to normal form, `.poly` assumes
##' the list is already in normal form and simply sets the class.
##
##' @title Polynomials
##' @param tms a list of monomials
##' @return a poly object
##' @seealso [pterm()] [set_term_order()]
##' @export
poly <- function(tms) {
  tms <- sort_terms(tms)

  ## Combine repeated terms
  l <- 0L
  nxt <- 1L
  while(nxt <= length(tms)) {
    tm <- tms[[nxt]]
    nxt <- nxt+1L
    while(nxt <= length(tms) && identical(tm$expt,tms[[nxt]]$expt)) {
      tm$coef <- tm$coef+tms[[nxt]]$coef
      nxt <- nxt+1L
    }
    if(tm$coef != 0L) tms[[l <- l+1L]] <- tm
  }
  if(l < length(tms)) tms <- tms[seq_len(l)]
  .poly(tms)
}

##' @rdname poly
##' @export
.poly <- function(tms) {
  class(tms) <- "poly"
  tms
}


## as.character method for poly
##' @export
as.character.poly <- function(x,...) {
  sign <- function(tm) if(substr(tm,1,1)=="-") tm else paste0("+",tm)
  if(length(x)==0L) return("0")
  x <- vapply(x, as.character,character(1))
  if(length(x) > 1L) x[-1] <- vapply(x[-1],sign,character(1))
  paste0(x,collapse="")
}

## Print method for poly
##' @export
print.poly <- function(x, ...) {
  cat(as.character(x),"\n")
}


##' Compute the degree of a polynomial.
##'
##' @title Degree of a polynomial
##' @param p a polynomial
##' @return the degree of the polynomial
##' @export
poly_degree <- function(p) {
  max(0L,vapply(p,function(tm) sum(tm$expt),integer(1)))
}


##' Parse a string representation of a signed polynomial term.
##'
##' Parses a signed polynomial term, returning the term and the
##' unconsumed component of the string.
##'
##' Note that the variables in `var` are translated to the standard
##' variables `x1,x2,x3,...`.
##'
##' @title Parse term
##' @param str a character vector
##' @param vars a character vector of variable names
##' @return a polynomial term
##' @importFrom gmp as.bigq
##' @export
parse_term <- function(str,vars) {

  isMatch <- function(m) !(length(m[[1L]])==1L && m[[1L]][1L]==-1L)

  str0 <- str
  expt <- integer(length(vars))

  ## Parse sign and the coeff
  if(isMatch(m <- regexec("^([+-])(\\d+(?:/\\d+)?)?",str,perl=TRUE))) {
    cparse <- regmatches(str,m)[[1L]]
    coef <- if(cparse[3L]=="") as.bigq(1L) else as.bigq(cparse[3L])
    if(cparse[2L]=="-") coef <- -coef
    str <- substr(str,attr(m[[1L]],"match.length")[1L]+1L,nchar(str))
    if(cparse[3L]=="") str <- paste0("*",str)
    ## Parse individual variables
    while(isMatch(m <- regexec("^\\*([a-z,A-Z]\\w*)(?:\\^(\\d+))?",str,perl=TRUE))) {
      parse <- regmatches(str,m)[[1L]]
      if(is.na(var <- match(parse[2L],vars))) stop("Unknown variable:",parse[2L])
      pwr <- if(parse[3L]=="") 1L else as.integer(parse[3L])
      expt[var] <- expt[var] + pwr
      str <- substr(str,attr(m[[1L]],"match.length")[1L]+1L,nchar(str))
    }
    ## A valid parse must have either variables or an explicit coef
    if(sum(expt)==0 && cparse[3L]=="")
      list(tm=NULL,str=str0)
    else
      list(tm=pterm(coef,expt),str=str)
  } else {
    list(tm=NULL,str=str0)
  }
}


##' Parse string representations of polynomials
##'
##' The `parse_poly` function converts a string representation of a
##' polynomial in the variables `vars` into a `poly` object. The
##' `parse_polys` function converts a list of strings into a list of
##' `poly` objects.
##'
##' Note that the variables in `var` are translated to the standard
##' variables `x1,x2,x3,...`.
##'
##' @title Polynomial parser
##' @param str a character vector
##' @param vars a character vector of variable names
##' @return a polynomial object
##' @examples
##' parse_poly("x^2*y^3-x*y*z",c("x","y","z"))
##' @export
parse_poly <- function(str,vars) {
   str <- gsub(" ","",str)
   if(!(substr(str,1L,1L)=="+" || substr(str,1L,1L)=="-")) str <- paste0("+",str)

   tms <- list()
   while(nchar(str)>0) {
    parse <- parse_term(str,vars)
    if(is.null(parse$tm)) stop("Parse error at term:",parse$str)
    tms[[length(tms)+1L]] <- parse$tm
    str <- parse$str
  }
  poly(tms)
}

##' @rdname parse_poly
##' @export
parse_polys <- function(str,vars) {
  lapply(str,parse_poly,vars=vars)
}


##' Sort the elements of a list.
##'
##' Sort the elements of a list according to a comparison function
##' `cmp(x,y)` which returns
##'
##' * -1 if `x < y`
##' * 0 if `x==y`
##' * 1 if `x > y`.
##'
##' @title Sort list
##' @param lst a list
##' @param cmp comparison function
##' @return the sorted list
##' @export
sort_list <- function(lst,cmp) {
  # In-place recursive quicksort
  qsort <- function(terms, low, high) {
    if(low < high) {
      pivot <- terms[[high]]  # Use the last element as pivot
      p <- low

      for (j in low:(high-1L)) {
        if (cmp(terms[[j]], pivot) > 0L) {
          terms[c(p, j)] <- terms[c(j, p)]  # Swap elements
          p <- p + 1L
        }
      }
      # Place pivot in the correct position
      terms[c(p, high)] <- terms[c(high, p)]
      terms <- qsort(terms, low, p-1L)
      terms <- qsort(terms, p+1L, high)
    }
    terms
  }

  qsort(lst,1L,length(lst))
}


## @rdname sort_list
## @export
sort_list1 <- function(lst, cmp) {

  # In-place recursive quicksort
  qsort <- function(lst, low, high) {
    if(low < high) {
      p <- low + floor((high - low)/2L)
      pivot <- lst[[p]]
      lst[c(high,p)] <- lst[c(p,high)]

      p <- low
      for (j in low:(high-1L)) {
        if (cmp(lst[[j]], pivot) > 0L) {
          lst[c(p,j)] <- lst[c(j,p)]
          p <- p + 1L
        }
      }
      # Move pivot to its final place
      lst[c(high,p)] <- lst[c(p,high)]
      lst <- qsort(lst, low, p-1L)
      lst <- qsort(lst, p+1L, high)
    }
    lst
  }

  qsort(lst,1L,length(lst))
}


##' Sort a list of polynomial terms.
##'
##'
##' @title Sort terms
##' @param terms a list of monomials
##' @return the sorted list
##' @export
sort_terms <- function(terms) {
  order <- term_order$order
  cmp <- function(x,y) order(x$expt,y$expt)
  sort_list(terms,cmp)
}

## @rdname sort_terms
## @export
sort_terms1 <- function(terms) {
  order <- term_order$order
  cmp <- function(x,y) order(x$expt,y$expt)
  sort_list1(terms,cmp)
}


##' Numerically evaluate a polynomial or a list of polynomials at an 
##' array of points.
##'
##' The function `neval_poly` evaluates a single polynomial at an array
##' of evaluation points returning the result as a numeric or complex 
##' vector.  The function `neval_polys` evaluates a list of polynomials
##' at an array of evaluation points returning the result as a numeric 
##' or complex matrix.
##'
##' @title Evaluate polynomial
##' @param p a polynomial
##' @param ps a list of polynomials
##' @param x a vector or matrix of evaluation points
##' @return the value of the polynomial at each x as complex or double.
##' @export
neval_poly <- function(p,x) {

  ## Coerce coeffs to numeric or complex
  coef <- function(tm) {
    if(is.complex(tm$coef)) tm$coef else as.numeric(tm$coef)
  }

  if(!is.matrix(x)) x <- matrix(x,1L,length(x))
  if(is.complex(x) || any(vapply(p,function(tm) is.complex(tm$coef),logical(1))))
    y <- complex(nrow(x))
  else
    y <- numeric(nrow(x))

  for(k in seq_len(nrow(x)))
    y[k] <- sum(vapply(p,function(tm) coef(tm)*prod(x[k,]^tm$expt),y[1]))

  y
}


##' @rdname neval_poly
##' @export
neval_polys <- function(ps,x) {
  do.call(cbind,lapply(ps,neval_poly,x=x))
}


##' Form the sum a*px + py where px and py are polynomials and a is a
##' monomial term.
##'
##' @title Polynomial axpy
##' @param acoef the coefficient of a
##' @param aexpt the exponent of a
##' @param px a polynomial
##' @param py a polynomial
##' @return the polynomial a*px + py
##' @export
poly_axpy <- function(acoef=1L,aexpt=0L,px,py) {

  lx <- length(px)
  if(lx==0L) return(py)
  ly <- length(py)

  tm_order <- term_order$order

  r <- vector("list",lx+ly)
  l <- 0L
  ix <- 0L
  iy <- 0L

  while(ix < lx || iy < ly) {

    ## py exhausted - copy remaining terms from px
    if(iy == ly) {
      while(ix < lx) {
        ix <- ix+1L
        r[[l<-l+1L]] <- pterm(px[[ix]]$coef*acoef,px[[ix]]$expt+aexpt)
      }
      break
    }

    ## px exhausted - copy remaining terms from py
    if(ix == lx) {
      while(iy < ly) r[[l<-l+1L]] <- py[[iy <- iy+1L]]
      break
    }

    ## Merge terms
    ix <- ix+1L
    axcoef <- px[[ix]]$coef*acoef
    axexpt <- px[[ix]]$expt+aexpt
    ## Insert greater terms from py
    while(iy < ly && (ord <- tm_order(py[[iy+1L]]$expt,axexpt))>0) r[[l<-l+1L]] <- py[[iy<-iy+1L]]
    ## Insert scaled term
    if(iy < ly && ord == 0) {
      coef <- axcoef + py[[iy <- iy+1L]]$coef
      if(coef != 0) r[[l<-l+1L]] <- pterm(coef,axexpt)
    } else {
      r[[l<-l+1L]] <- pterm(axcoef,axexpt)
    }
  }
  if(length(r) != l) r <- r[seq_len(l)]
  class(r) <- "poly"
  r
}


##' Change the variable order of a polynomial
##'
##' Permute the variables in a polynomial according to a permutation
##' vector `order`.  For example, if `order` is `c(3,2,1)`,then the 
##' standard variables `x1` and `x3` are interchanged and the terms 
##' are reordered according to the term order.
##'
##' If `order` specifies more or fewer variables than are present in the
##' current variable set, the exponent vectors are resized to encode the
##' new variable set.
##'
##' @title Reorder variables
##' @param p a polynomial
##' @param order an integer vector defining the new variable order
##' @return a polynomial with the variables in the new order
##' @export
reorder_vars <- function(p,order) {
  if(any(duplicated(order))) warning("Invalid variable permutation")
  l <- length(order)

  reorder <- function(expt) {
    if(length(expt) < l) expt <- c(expt,integer(l-length(expt)))
    expt[order]
  }
  
  poly(lapply(p,function(tm) pterm(tm$coef,reorder(tm$expt))))
}


##' Polynomial arithmetic
##'
##' `poly_add` computes p1+p2, `poly_sub` computes p1-p2, `poly_mul`
##' computes p1*p2, and `poly_scale` computes a*p where a is a term.
##' @title Polynomial arithmetic
##' @param p1,p2 polynomials
##' @param acoef the coefficient of a
##' @param aexpt the exponent of a
##' @return a polynomial
##' @rdname poly_arithmetic
##' @export
poly_add <- function(p1,p2) {
  poly_axpy(1L,0L,p1,p2)
}

##' @rdname poly_arithmetic
##' @export
poly_sub <- function(p1,p2) {
  poly_axpy(-1L,0L,p2,p1)
}

##' @rdname poly_arithmetic
##' @export
poly_mul <- function(p1,p2) {
  r <- list()
  if(length(p1) < length(p2))
    for(tm in p1) r <- poly_axpy(tm$coef,tm$expt,p2,r)
  else
    for(tm in p2) r <- poly_axpy(tm$coef,tm$expt,p1,r)
  .poly(r)
}

##' @rdname poly_arithmetic
##' @export
poly_scale <- function(acoef=1L,aexpt=0L,p1) {
  if(acoef==0) return(.poly(list()))

  for(k in seq_along(p1)) {
    p1[[k]]$coef <- p1[[k]]$coef*acoef
    p1[[k]]$expt <- p1[[k]]$expt+aexpt
  }
  .poly(p1)
}


##' Normalize a polynomial to have leading coefficient 1
##'
##' @title Normalize
##' @param p a polynomial
##' @return the normalized polynomial
##' @export
normalize <- function(p) {
  if((lc <- p[[1L]]$coef) != 1)
    for(k in seq_along(p)) p[[k]]$coef <- p[[k]]$coef/lc
  p
}


##' Compute the derivative of a polynomial.
##'
##' @title Derivative
##' @param p a polynomial
##' @param order vector of non-negative integers specifying the
##' order of the derivative
##' @return the polynomial that is the derivative of `p`.
##' @export
derivative <- function(p,order) {
  order <- as.integer(order)

  dp <- vector("list",length(p))
  l <- 0L
  for(tm in p) {
    if(all(tm$expt>=order)) {
      ord <- order
      while(any(ord>0)) {
        pos <- ord>0
        tm$coef <- tm$coef*prod(tm$expt[pos])
        tm$expt[pos] <- tm$expt[pos]-1L
        ord[pos] <- ord[pos]-1L
      }
      dp[[l<-l+1L]] <- tm
    }
  }
  if(l<length(p)) .poly(dp[seq_len(l)]) else .poly(dp)
}


##' Reduce polynomials.
##'
##' There are four functions for reducing polynomials:
##'
##' * `reduce_pp` reduces a polynomial by another polynomial
##' * `reduce_ps` reduces a polynomial by a list of polynomials
##' * `reduce_sp` reduces a list of polynomials by a polynomial
##' * `reduce_s` reduces a list of polynomials
##'
##' For functions `reduce_pp` and `reduce_sp` the reductions are
##' performed in a single pass, but for `reduce_ps` and `reduce_s` the
##' reductions are iterated until the result cannot be further
##' reduced.
##'
##' @title Reduce
##' @param p1 a polynomial
##' @param p2 a polynomial
##' @param p a polynomial
##' @param ps a list of polynomials
##' @return the reduced polynomial
##' @rdname reduce
##' @export
reduce_pp <- function(p1,p2) {

  l1 <- length(p1)
  l2 <- length(p2)
  if(l1==0 || l2==0) return(p1)
  k <- 1L

  ## Cache active binding
  tm_order <- term_order$order

  ## Scan p1 for terms to reduce
  while(k<=l1) {
    if(tm_order(p2[[1]]$expt,p1[[k]]$expt) > 0) break

    ## If p2[[1]] divides p1[[k]] - merge terms
    if(all((rexpt <- p1[[k]]$expt-p2[[1]]$expt) >= 0)) {
      rcoef <- -p1[[k]]$coef/p2[[1]]$coef

      ## Allocate result
      r <- vector("list",l1+l2-2L)
      r[seq_len(k-1L)] <- p1[seq_len(k-1L)]

      ## Merge terms
      l <- k-1L
      i1 <- k
      i2 <- 1L
      while(i1 < l1 || i2 < l2) {

        ## p1 exhausted - copy remaining terms from p2
        if(i1 == l1) {
          while(i2 < l2) {
            i2 <- i2+1L
            r[[l<-l+1L]] <- pterm(p2[[i2]]$coef*rcoef,p2[[i2]]$expt+rexpt)
          }
          break
        }

        ## p2 exhausted - copy remaining terms from p1
        if(i2 == l2) {
          while(i1<l1) r[[l<-l+1L]] <- p1[[i1<-i1+1L]]
          break
        }

        ## Exponent of reducing term
        expt2 <- p2[[i2 <- i2+1L]]$expt+rexpt
        ## Insert greater terms from p1
        while(i1 < l1 && (ord <- tm_order(p1[[i1+1L]]$expt,expt2))>0) r[[l<-l+1L]] <- p1[[i1<-i1+1L]]
        ## Insert reducing term
        if(i1 < l1 && ord == 0) {
          coef <- p1[[i1 <- i1+1L]]$coef + p2[[i2]]$coef*rcoef
          if(coef != 0) r[[l<-l+1L]] <- pterm(coef,expt2)
        } else {
          r[[l<-l+1L]] <- pterm(p2[[i2]]$coef*rcoef,expt2)
        }
      }

      ## Replace p1 with reduced polynomial
      p1 <- r
      l1 <- l
    } else {
      k <- k+1L
    }
  }
  if(length(p1) !=l1) p1 <-p1[seq_len(l1)]
  class(p1) <- "poly"
  p1
}


##' @rdname reduce
##' @export
reduce_ps <- function(p,ps) {
  if(length(p)==0L) return(p)
  p0 <- NULL
  while(!identical(p,p0)) {
    p0 <- p
    for(q in ps) {
      p <- reduce_pp(p,q)
      if(length(p)==0L) return(p)
    }
  }
  p
}


##' @rdname reduce
##' @export
reduce_sp <- function(ps,p) {
  for(k in seq_along(ps))
    ps[[k]] <- reduce_ps(ps[[k]],ps)
  if(any((ls <- lengths(ps))==0L)) ps[ls>0L] else ps
}


##' @rdname reduce
##' @export
reduce_s <- function(ps) {
  ps0 <- NULL
  while(!identical(ps,ps0)) {
    ps0 <- ps
    for(i in seq_along(ps)) {
      p <- ps[[i]]
      if(length(p)>0L) {
        for(j in seq_along(ps)[-i]) {
          ps[[j]] <- reduce_pp(ps[[j]],p)
        }
      }
    }
    if(any((ls <- lengths(ps))==0L)) ps <- ps[ls>0L]
  }
  ps
}


##' Compute the S-polynomial of two polynomials
##'
##' @title S-polynomial
##' @param p1 a polynomial
##' @param p2 a polynomial
##' @return the S-polynomial
##' @export
spoly <- function(p1,p2) {

  l1 <- length(p1)
  l2 <- length(p2)
  if(l1==1L && l2==1L) return(.poly(list()))

  ## Cache active binding
  tm_order <- term_order$order

  ## Compute leading coefficients
  lc1 <- p1[[1L]]$coef
  lc2 <- p2[[1L]]$coef
  expt <- pmax(p1[[1L]]$expt,p2[[1L]]$expt)
  rexpt1 <- expt - p1[[1L]]$expt
  rexpt2 <- expt - p2[[1L]]$expt

  ## Allocate result
  r <- vector("list",l1+l2-2L)

  ## Merge terms
  l <- 0L
  i1 <- 1L
  i2 <- 1L
  while(i1 < l1 || i2 < l2) {

    ## p1 exhausted - add remaining terms from p2
    if(i1 == l1) {
      while(i2 < l2) {
        i2 <- i2+1L
        r[[l<-l+1L]] <- pterm(-p2[[i2]]$coef/lc2,p2[[i2]]$expt+rexpt2)
      }
      break
    }

    ## p2 exhausted - copy remaining terms from p1
    if(i2 == l2) {
      while(i1 < l1) {
        i1 <- i1+1L
        r[[l<-l+1L]] <- pterm(p1[[i1]]$coef/lc1,p1[[i1]]$expt+rexpt1)
      }
      break
    }

    ## Merge terms
    expt1 <- p1[[i1+1L]]$expt+rexpt1
    expt2 <- p2[[i2+1L]]$expt+rexpt2
    ord <- tm_order(expt1,expt2)
    if(ord==0L) {
      coef <- p1[[i1 <- i1+1L]]$coef/lc1 - p2[[i2 <- i2+1L]]$coef/lc2
      if(coef != 0) r[[l<-l+1L]] <- pterm(coef,expt1)
    } else if(ord > 0L) {
      r[[l<-l+1L]] <- pterm(p1[[i1 <- i1+1L]]$coef/lc1,expt1)
    } else {
      r[[l<-l+1L]] <- pterm(-p2[[i2 <- i2+1L]]$coef/lc2,expt2)
    }
  }
  if(length(r) != l) r <- r[seq_len(l)]
  .poly(r)
}


##' Compute the reduced Groebner basis of a set of polynomials
##'
##' @title Groebner basis
##' @param ps a list of polynomials
##' @return the Groebner basis
##' @export
groebner <- function(ps) {
  ps <- reduce_s(ps)
  while(TRUE) {
    qs <- list()
    for(i in seq_along(ps)) {
      for(j in seq_len(i-1L)) {
        pi <- ps[[i]]
        pj <- ps[[j]]
        if(sum(pi[[1L]]$expt*pj[[1L]]$expt)>0L) {
          p <- spoly(pi,pj)
          p <- reduce_ps(p,ps)
          if(length(p)>0L) qs[[length(qs)+1L]] <- p
        }
      }
    }
    if(length(qs)==0L) break
    ps <- reduce_s(c(ps,qs))
  }
  ps
}


##' Generate a set of monomials upto a given degree.
##'
##' @title Monomial set
##' @param nvar the number of variables
##' @param deg the degree of the monomials
##' @return a list of integer vectors representing the monomial set
##' @export
monomial_set <- function(nvar,deg) {

  ## The monomials
  expts <- list()

  ## The candidate set of monomials
  prv <- list(integer(nvar))
  expts <- prv
  for(k in seq_len(deg)) {
    nxt <- vector("list",nvar*length(prv))
    l <- 0L
    for(i in seq_len(nvar)) {
      for(j in seq_along(prv)) {
        ## Multiply a previous candidate by the i-th variable
        expt <- prv[[j]]
        expt[i] <- expt[i]+1L
        nxt[[l <- l+1L]] <- expt
      }
    }
    ## Add the candidates to the set
    prv <- unique(nxt[seq_len(l)])
    expts <- c(expts,prv)
  }
  expts
}

##' Compute a monomial basis from a Groebner basis or a list of
##' reducible monomials.
##'
##' Given a Groebner basis or a list of reducible monomials terms
##' compute the corresponding monomial basis. The monomial basis is
##' returned as a list of integer vectors of the exponents of the
##' of the standard variables `x1,x2,...` in the monomial term. The
##' `monomial_basis0` function computes a monomial basis given a list
##' of reducible terms, and `monomial_basis` computes a monomial 
##' basis from a Groebner basis.
##'
##' @title Monomial basis
##' @param gb a Groebner basis
##' @param rexpts a list of reducible exponent vectors
##' @return a list of the exponent vectors in the monomial basis
##' @export
monomial_basis <- function(gb) {
  monomial_basis0(lapply(gb,function(p) p[[1L]]$expt))
}


##' @rdname monomial_basis
##' @export
monomial_basis0 <- function(rexpts) {

  nvar <- length(rexpts[[1L]])

  ## The exponent vectors in the monomial basis
  expts <- list()

  ## The candidate set of monomials
  prv <- list(integer(nvar))
  while(length(prv)>0L) {
    nxt <- vector("list",nvar*length(prv))
    l <- 0L
    for(i in seq_len(nvar)) {
      for(j in seq_along(prv)) {
        ## Multiply a previous candidate by the i-th variable
        expt <- prv[[j]]
        expt[i] <- expt[i]+1L
        ## Check if the monomial is reducible
        keep <- TRUE
        for(rexpt in rexpts)
          if(all(rexpt<=expt)) { keep <- FALSE; break }
        if(!keep) next
        nxt[[l <- l+1L]] <- expt
      }
    }
    ## Add the previous candidates to the basis and update
    expts <- c(expts,prv)
    prv <- unique(nxt[seq_len(l)])
  }
  expts
}


##' Compute the multiplication matrix for a variable.
##'
##' Computes the matrix representation in the basis `ms` of
##' multiplication by the `v`-th standard variable.
##'
##' @title Multiplication matrix
##' @param v a variable index
##' @param ms a list of the exponent vectors of the monomial basis
##' @param gb a Groebner basis
##' @return the multiplication matrix
##' @export
multiplication_matrix <- function(v,ms,gb) {

  as.nc<- function(x) {
    if(is.numeric(x) || is.complex(x)) x else as.numeric(x)
  }

  n <- length(ms)
  M <- matrix(0L,n,n)
  for(j in seq_len(n)) {
    m <- ms[[j]]
    m[v] <- m[v]+1L
    p <- .poly(list(pterm(1L,expt=m)))
    p <- reduce_ps(p,gb)
    if(length(p)>0L) {
      is <- match(lapply(p,function(tm) tm$expt),ms)
      cs <- sapply(p,function(tm) as.nc(tm$coef))
      M[is,j] <- cs
    }
  }
  M
}

##' Compute the derivative matrix for a variable.
##'
##' Computes the matrix representation in the basis `ms` of
##' differentiation by the `v`-th standard variable.
##'
##' @title Derivative matrix
##' @param v a variable index
##' @param ms a list of the exponent vectors of the monomial basis
##' @return the derivative matrix
##' @export
derivative_matrix <- function(v,ms) {
  n <- length(ms)
  M <- matrix(0L,n,n)
  for(j in seq_len(n)) {
    m <- ms[[j]]
    if(m[v]>1L) {
      m[v] <- m[v]-1L
      i <- match(list(m),ms)
      if(!is.na(i)) M[i,j] <- m[v]+1
    }
  }
  M
}

##' Compute the Macaulay matrix for a polynomial system.
##'
##' @title Macaulay matrix
##' @param ps a list of polynomials
##' @param tdeg the degree of the Macaulay matrix
##' @return a list with elements
##'   * `M`: the Macaulay matrix
##'   * `ms`: the monomial basis
##' @export
macaulay_matrix <- function(ps,tdeg=sum(degs-1L)+1L) {

  as.nc<- function(x) {
    if(is.numeric(x) || is.complex(x)) x else as.numeric(x)
  }

  degs <- vapply(ps,poly_degree,0L)
  nvar <- length(ps[[1L]][[1L]]$expt)
  ms <- rev(monomial_set(nvar,tdeg))

  ## Scale the polynomials by the monomials
  ps <- unlist(lapply(seq_along(ps),function(k)
    lapply(monomial_set(nvar,tdeg-degs[k]),
           function(m) poly_scale(1L,m,ps[[k]]))),
    recursive=FALSE)

  ## Create the Macaulay matrix
  M <- matrix(0,length(ps),length(ms))
  for(i in seq_along(ps)) {
    p <- ps[[i]]
    js <- match(lapply(p,function(tm) tm$expt),ms)
    cs <- sapply(p,function(tm) as.nc(tm$coef))
    M[i,js] <- cs
  }
  list(M=M,ms=ms)
}


##' reduce a matrix to reduced row echelon form by Gauss-Jordan
##' elimination.
##'
##' @title Gauss-Jordan elimination
##' @param A a matrix
##' @return a list with elements
##'   * `A`: the reduced matrix
##'   * `reduced`: the indices of the reduced columns
##'   * `rank`: the rank of the matrix
##' @export
gauss_jordan <- function(A) {

  ## Rescale rows
  A <- A/apply(abs(A),1L,max)

  reduced <- integer(nrow(A))
  rank <- 0L
  j <- 1L
  ## Loop over rows to reduce
  for(i in seq_len(nrow(A))) {
    ## Find a pivot element
    while(j <= ncol(A) && all((as <- abs(A[seq.int(i,nrow(A)),j]))<1.0E-12)) j <- j+1L
    if(j > ncol(A)) break
    piv <- which.max(as) + i - 1L
    ## Swap rows if necessary and normalize
    if(piv != i) A[c(i,piv),] <- A[c(piv,i),]
    A[i,] <- A[i,]/A[i,j]
    ## Reduce the remaining rows
    A[-i,] <- A[-i,] - outer(A[-i,j],A[i,])
    #for(k in seq_len(nrow(A))[-i])
      #if(A[k,j] != 0) A[k,] <- A[k,]-A[k,j]*A[i,]
    rank <- rank+1L
    reduced[i] <- j
    j <- j+1L
    if(j > ncol(A)) break
  }
  if(rank < length(reduced)) reduced <- reduced[seq_len(rank)]
  list(A=A,reduced=reduced,rank=rank)
}


##' Compute a monomial basis determined from the dependencies of a
##' Macaulay matrix.
##'
##' Computes a monomial basis from the linear dependences of the
##' Macaulay matrix of a polynomial system.  Returns a monomial basis
##' and a reduction matrix that can be used to reduce terms to the
##' monomial basis.
##'
##' @title Macaulay monomial basis
##' @param ps a list of polynomials
##' @return a list with elements
##'   * `A`: reduction matrix
##'   * `term`: reducible terms
##'   * `basis`: the monomial basis
##' @export
macaulay_monomial_basis <- function(ps) {
  mc <- macaulay_matrix(ps)
  gj <- gauss_jordan(mc$M)
  bs <- monomial_basis0(mc$ms[gj$reduced])
  bidx <- match(bs,mc$ms)
  ## Drop terms that reduce to terms not in the basis
  keep <- which(rowSums(abs(gj$A[seq_len(gj$rank),-c(gj$reduced,bidx),drop=FALSE]))==0)
  list(A=rbind(-gj$A[keep,bidx,drop=FALSE],diag(1,length(bs))),
       term=mc$ms[c(gj$reduced[keep],bidx)],
       basis=bs)
}

##' Compute the multiplication matrix for a variable from a Macaulay
##' monomial basis.
##'
##' @title Multiplication matrix from a Macaulay monomial basis
##' @param v a variable index
##' @param mb an object returned by [macaulay_monomial_basis()]
##' @return the multiplication matrix
##' @export
macaulay_multiplication_matrix <- function(v,mb) {

  ms <- mb$basis
  for(k in seq_along(ms)) ms[[k]][v] <- ms[[k]][v]+1L
  is <- match(ms,mb$term)
  if(any(is.na(is))) stop("Irreducible terms found")
  t(mb$A[is,])
}




##' Polish a root of a polynomial system using Newton's method.
##'
##' Given a root of a polynomial system, this function performs `n`
##' iterations of Newton's method to polish the root.
##'
##' @title Newton's method
##' @param ps a list of polynomials
##' @param x a vector or matrix of initial estimates
##' @param n the number of iterations
##' @return a root of the polynomial system
##' @export
newton_polish <- function(ps,x,n) {
  ## Number of variables
  nvars <- length(ps[[1L]][[1L]]$expt)

  ## Compute the Jacobian
  deriv <- function(i) {
    ord <- integer(nvars)
    ord[i] <- 1L
    lapply(ps,function(p) derivative(p,ord))
  }
  dp <- lapply(seq_len(nvars),deriv)

  if(!is.matrix(x)) x <- matrix(x,1,length(x))
  for(k in seq_len(nrow(x))) {
    z <- x[k,]
    ## Newton iterations
    for(i in seq_len(n)) {
      y <- drop(neval_polys(ps,z))
      J <- do.call(cbind,lapply(dp,function(dp) t(neval_polys(dp,z))))
      z <- z - qr.solve(J,y)
    }
    x[k,] <- z
  }
  drop(x)
}


##' Compute an eigenbasis of a matrix.
##'
##' Returns a list of matrices, where the columns of each matrix
##' form an orthonormal basis for an eigenspace of the matrix.
##'
##' If `tol` is too small, equal eigenvalues are erroneously identified
##' as distinct, and the bases for some eigenspaces are split.
##'
##' @title Eigenbasis
##' @param A a matrix
##' @param tol tolerance used to determine if two eigenvalues differ.
##' @return a list of matrices, one for each distinct eigenvalue, the
##'   columns of each are the eigenbasis for that eigenvalue
##' @export
eigenbasis <- function(A, tol=1e-6) {

  ## Compute the independent columns of a matrix
  independent <- function(X) {
    if(ncol(X) == 1L) return(X)
    QR <- qr(X,LAPACK=TRUE)
    rs <- abs(diag(qr.R(QR)))
    rank <- sum(rs > 1.0E-7*max(rs))
    qr.Q(QR)[, seq_len(rank), drop = FALSE]
  }

  ## Compute the sets of distinct eigenvalues
  eig <- eigen(A)
  sets <- cumsum(c(TRUE,abs(diff(eig$values))>tol))

  lapply(unique(sets),
         function(set) independent(eig$vectors[,sets==set,drop=FALSE]))
}


## Alternate implementation of eigenbasis based on SVD
eigenbasis1 <- function(A,tol=1.0E-6) {

  M <- max(abs(A))

  nullspace <- function(A) {
    SVD <- svd(A)
    rank <- sum(SVD$d > 1.0E-12*M)
    if(rank==ncol(A)) return(matrix(0,nrow(A),0))
    if(rank==0L) return(diag(1,nrow(A)))
    SVD$v[,seq.int(rank+1L,ncol(A)),drop=FALSE]
  }

  ## Distinct eigenvalues
  eig <- eigen(A,only.values=TRUE)
  lambda <- eig$values[c(TRUE,abs(diff(eig$values))>tol)]

  lapply(lambda,function(l) nullspace(A - diag(l,nrow(A))))
}


## Alternate implementation of eigenbasis based on QR
eigenbasis2 <- function(A,tol=1.0E-6) {

  M <- max(abs(A))

  nullspace <- function(A) {
    QR <- qr(t(A),LAPACK=TRUE)
    rs <- abs(diag(qr.R(QR)))
    rank <- sum(rs > 1.0E-12*M)
    if(rank==ncol(A)) return(matrix(0,nrow(A),0))
    if(rank==0L) return(diag(1,nrow(A)))
    qr.Q(QR,complete=TRUE)[,-seq_len(rank),drop=FALSE]
  }

  ## Distinct eigenvalues
  eig <- eigen(A,only.values=TRUE)
  lambda <- eig$values[c(TRUE,abs(diff(eig$values))>tol)]

  lapply(lambda,function(l) nullspace(A - diag(l,nrow(A))))
}


##' Compute a common eigenbasis for a list of commuting matrices.
##'
##' Compute the common eigenbasis of a list of commuting matrices. No
##' check is performed to ensure the matrices are commute.
##'
##' `common_eigenbasis` computes a common eigenbasis by computing the
##' eigenbasis of the first matrix and then rotating and splitting the
##' basis to align with the eigenbasis of each subsequent matrix.
##'
##' `common_eigenbasis0` computes a common eigenbasis by computing the
##' eigenbasis of a random linear combination of the matrices.
##'
##'  No check is performed to ensure the matrices are commute.
##'
##' If `tol` is too small, equal eigenvalues are erroneously identified
##' as distinct, and the bases for some eigenspaces are split.
##'
##' @title Common eigenbasis
##' @param As a list of commuting matrices
##' @param tol tolerance used to determine if two eigenvalues differ.
##' @param left if `TRUE`, compute a common left eigenbasis
##' @return a list of matrices, the columns of each matrix are the
##'   common eigenbasis for the corresponding matrix in `As`
##' @export
common_eigenbasis <- function(As,tol=1e-6,left=FALSE) {

  ## Transpose for the left eigenbasis
  if(left) As <- lapply(As,t)

  ## Compute the eigenbasis for the first matrix
  basis <- eigenbasis(As[[1L]])
  n <- sum(vapply(basis,ncol,1L))

  ## Rotate and split the basis for each subsequent matrix
  for(A in As[-1L]) {
    m <- 0L
    basis1 <- vector("list",n)
    for(i in seq_along(basis)) {
      if(ncol(basis[[i]])==1L) {
        basis1[[m <- m+1L]] <- basis[[i]]
      } else {
        Ak <- crossprod(Conj(basis[[i]]),A%*%basis[[i]])
        bs <- eigenbasis(Ak,tol=tol)
        for(k in seq_along(bs))
          basis1[[m <- m+1L]] <- basis[[i]]%*%bs[[k]]
      }
    }
    basis <- basis1[seq_len(m)]
  }
  basis
}


##' @rdname common_eigenbasis
##' @importFrom stats runif
##' @export
common_eigenbasis0 <- function(As,tol=1e-6,left=FALSE) {
  ## Transpose for the left eigenbasis
  if(left) As <- lapply(As,t)
  ## Random linear combination of the matrices
  S <- 0
  for(A in As) S <- S+runif(1)*A
  eigenbasis(S,tol=tol)
}


##' Extract the roots of a polynomial system
##'
##' `roots` computes the roots of a polynomial system from the
##' eigenvalues of its multiplication matrices using a common basis
##' computed by [common_eigenbasis()].  The roots are returned in the
##' same order as the multiplication matrices.
##'
##' @title Roots from eigenvalues and eigenvectors
##' @param Ms a list of multiplication matrices
##' @param Bs an eigenbasis
##' @param unique In degenerate cases, return each root only once
##' @return a matrix where each row is a root of the polynomial system
##' @rdname roots
##' @export
roots <- function(Ms,Bs,unique=TRUE) {
  if(unique) Bs <- lapply(Bs,function(B) B[,1,drop=FALSE])
  do.call(cbind,lapply(Ms,function(M) do.call(c,lapply(Bs,function(B) diag(crossprod(Conj(B),M%*%B))))))
}


##' Extract the roots of a polynomial system from the left eigenvectors
##' of the multiplication matrices.
##'
##' If `Bs` is a common left eigenbasis calculated from the transposed
##' multiplication matrices, then this function will extract the values
##' of the standard variables that appear in the monomials in the basis
##' `ms` at the roots of the polynomial system.
##'
##' @title Roots from eigenvectors
##' @param ms a monomial basis
##' @param Bs a common left eigenbasis
##' @return a matrix where each row is a root of the polynomial system
##' @export
roots_left_eigenbasis <- function(ms,Bs) {

  ## Find constant term
  cnst <- which(vapply(ms,function(m) all(m==0L),logical(1L)))

  ## Find the standard variables in the monomial basis and sort
  idx <- which(vapply(ms,function(m) sum(m)==1L,logical(1L)))
  vr <- vapply(ms[idx],function(m) which(m==1L),integer(1L))
  ord <- order(vr)
  vr <- vr[ord]
  idx <- idx[ord]

  ## Normalize the eigenbasis
  roots <- do.call(rbind,lapply(Bs,function(B) t(B[idx,,drop=FALSE])/B[cnst,]))
  colnames(roots) <- paste0("x",vr)
  roots
}


##' Find the roots of a polynomial system
##'
##' @title Roots of a polynomial system
##' @param ps a list of polynomials
##' @param gb a Groebner basis
##' @param newton the number of Newton iterations used to refine roots
##' @param tol tolerance used to determine if two eigenvalues differ.
##' @return a matrix where each row is a root of the polynomial system
##' @export
solve_polys <- function(ps,gb=groebner(ps),newton=0,tol=1.0E-6) {
  ms <- monomial_basis(gb)
  Ms <- lapply(seq_along(ms[[1]]),function(m) multiplication_matrix(m,ms,gb))
  Bs <- common_eigenbasis(Ms,tol=tol)
  rs <- roots(Ms,Bs)
  if(newton>0) rs <- newton_polish(ps,rs,newton)
  rs
}

