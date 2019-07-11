check_belong                 <- function(value, name, allowed, len) {
  if (is.infinite(len)) {
    for (i in 1:length(value)) {
      if (!(value[i] %in% allowed)){
        stop(name, " must contain values only in ",
             paste(allowed, collapse = ", "))
      }
    }
  } else {
    if (any(length(value) != 1, !(value %in% allowed))) {
      stop(name, " must be set to one of ", paste(allowed, collapse = ", "))
    }
  }
}

check_default                <- function(value, name, default) {
  if (!all(value == default)) {
    stop(name, " has been changed from its default value, but this will have ",
         "no effect given the value of the other arguments")
  }
}

check_des                    <- function(des, type, stages) {
  allowed_classes <- c("ph2rand_des_one_stage", "ph2rand_des_two_stage")
  if (!any(class(des) %in% allowed_classes[stages])) {
    stop("des must be of class ", paste(allowed_classes, sep = " or "))
  }
  if (type == "bernard") {
    if (des$type != "bernard") {
      stop("des must correspond to a design of type ", type)
    }
  } else if (type == "binomial") {
    if (des$type != "binomial") {
      stop("des must correspond to a design of type ", type)
    }
  } else if (type == "fisher") {
    if (des$type != "fisher") {
      stop("des must correspond to a design of type ", type)
    }
  } else if (type == "single_double") {
    if (des$type != "single_double") {
      stop("des must correspond to a design of type ", type)
    }
  }
}

check_des_bernard            <- function(des) {
  if (any(!("ph2rand_bernard" %in% class(des)), !is.list(des))) {
    stop("des must be of class ph2rand_bernard and of class list")
  } else if (any(is.null(des$J), !(des$J %in% c(1, 2)))) {
    stop("des$J must be equal to 1 or 2")
  }
  if (des$J == 1) {
    if (any(is.null(des$n0), length(des$n0) != 1, !is.numeric(des$n0),
            !is.finite(des$n0), des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 1, des$n0 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 1, !is.numeric(des$n1),
            !is.finite(des$n1), des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 1, des$n1 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$e), length(des$e) != 1, !is.numeric(des$e))) {
      stop("For des$J = 1, des$e must be a single numeric")
    } else if (is.infinite(des$e)) {
      warning("For des$J = 1, des$e should typically be finite")
    }
  } else {
    if (any(is.null(des$n0), length(des$n0) != 2, !is.numeric(des$n0),
            des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 2, des$n0 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 2, !is.numeric(des$n1),
            des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 2, des$n1 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$e1), length(des$e1) != 1, !is.numeric(des$e1),
            des$e1 == -Inf)) {
      stop("For des$J = 2, des$e1 must be a single numeric not equal to -Inf")
    }
    if (any(is.null(des$f1), length(des$f1) != 1, !is.numeric(des$f1),
            des$f1 >= des$e1)) {
      stop("For des$J = 2, des$f1 must be a single numeric that is strictly ",
           "less than des$e1")
    }
    if (any(is.null(des$e2), length(des$e2) != 1, !is.numeric(des$e2))) {
      stop("For des$J = 2, des$e2 must be a single numeric")
    } else if (is.infinite(des$e2)) {
      warning("For des$J = 2, des$e2 should typically be finite")
    }
  }
}

check_des_binomial           <- function(des) {
  if (any(!("ph2rand_binomial" %in% class(des)), !is.list(des))) {
    stop("des must be of class ph2rand_binomial and of class list")
  } else if (any(is.null(des$J), !(des$J %in% c(1, 2)))) {
    stop("des$J must be equal to 1 or 2")
  }
  if (des$J == 1) {
    if (any(is.null(des$n0), length(des$n0) != 1, !is.numeric(des$n0),
            !is.finite(des$n0), des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 1, des$n0 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 1, !is.numeric(des$n1),
            !is.finite(des$n1), des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 1, des$n1 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$e), length(des$e) != 1, !is.numeric(des$e))) {
      stop("For des$J = 1, des$e must be a single numeric")
    } else if (any(is.infinite(des$e), des$e%%1 != 0, des$e <= -des$n0,
                   des$e > des$n1)) {
      warning("For des$J = 1, des$e should typically be a single integer in ",
              "(-des$n0, des$n1]")
    }
  } else {
    if (any(is.null(des$n0), length(des$n0) != 2, !is.numeric(des$n0),
            des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 2, des$n0 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 2, !is.numeric(des$n1),
            des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 2, des$n1 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$e1), length(des$e1) != 1, !is.numeric(des$e1),
            des$e1 == -Inf)) {
      stop("For des$J = 2, des$e1 must be a single numeric not equal to -Inf")
    }
    if (any(is.null(des$f1), length(des$f1) != 1, !is.numeric(des$f1),
            des$f1 >= des$e1)) {
      stop("For des$J = 2, des$f1 must be a single numeric that is strictly ",
           "less than des$e1")
    }
    if (any(is.null(des$e2), length(des$e2) != 1, !is.numeric(des$e2))) {
      stop("For des$J = 2, des$e2 must be a single numeric")
    }
    if (all(is.finite(des$e1), is.finite(des$f1))) {
      if (any(des$e1%%1 != 0, des$e1 <= -des$n0[1] + 1, des$e1 > des$n1[1])) {
        warning("For des$J = 2, when des$e1 and des$f1 are finite, des$e1 ",
                "should typically be a single integer in (-des$n0[1] + 1, ",
                "des$n1[1]]")
      }
      if (any(des$f1%%1 != 0, des$f1 < -des$n0[1], des$f1 >= des$n1[1] - 1)) {
        warning("For des$J = 2, when des$e1 and des$f1 are finite, des$f1 ",
                "should typically be a single integer in [-des$n0[1], ",
                "des$n1[1] - 1)")
      }
      if (any(is.infinite(des$e2), des$e2%%1 != 0,
              des$e2 < des$f1 + 2 - des$n0[2],
              des$e2 > des$e1 - 1 + des$n1[2])) {
        warning("For des$J = 2, when des$e1 and des$f1 are finite, des$e2 ",
                "should typically be a single integer in ",
                "[des$f1 + 2 - des$n0[2], des$e1 - 1 + des$n1[2]]")
      }
    } else if (is.finite(des$e1)) {
      if (any(des$e1%%1 != 0, des$e1 <= -des$n0[1], des$e1 > des$n1[1])) {
        warning("For des$J = 2, when des$e1 is finite and des$f1 is not, ",
                "des$e1 should typically be a single integer in (-des$n0[1], ",
                "des$n1[1]]")
      }
      if (any(is.infinite(des$e2), des$e2%%1 != 0,
              des$e2 <= -des$n0[1] - des$n0[2],
              des$e2 > des$e1 - 1 + des$n1[2])) {
        warning("For des$J = 2, when des$e1 is finite and des$f1 is not, ",
                "des$e2 should typically be a single integer in ",
                "(-des$n0[1] - des$n0[2], des$e1 - 1 + des$n1[2]]")
      }
    } else if (is.finite(des$f1)) {
      if (any(des$f1%%1 != 0, des$f1 <= -des$n0[1], des$f1 >= des$n1[1])) {
        warning("For des$J = 2, when des$f1 is finite and des$e1 is not, ",
                "des$f1 should typically be a single integer in [-des$n0[1], ",
                "des$n1[1])")
      }
      if (any(is.infinite(des$e2), des$e2%%1 != 0,
              des$e2 <= des$f1 + 2 - des$n0[2],
              des$e2 > des$n1[1] + des$n1[2])) {
        warning("For des$J = 2, when des$f1 is finite and des$e1 is not, ",
                "des$e2 should typically be a single integer in ",
                "(des$f1 + 2 - des$n0[2], des$n1[1] + des$n1[2])")
      }
    }
  }
}

check_des_fisher             <- function(des) {
  if (any(!("ph2rand_fisher" %in% class(des)), !is.list(des))) {
    stop("des must be of class ph2rand_fisher and of class list")
  } else if (any(is.null(des$J), !(des$J %in% c(1, 2)))) {
    stop("des$J must be equal to 1 or 2")
  }
  if (des$J == 1) {
    if (any(is.null(des$n0), length(des$n0) != 1, !is.numeric(des$n0),
            !is.finite(des$n0), des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 1, des$n0 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 1, !is.numeric(des$n1),
            !is.finite(des$n1), des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 1, des$n1 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$e), length(des$e) != des$n0 + des$n1 + 1,
            !is.numeric(des$e))) {
      stop("For des$J = 1, des$e must be a numeric vector of length des$n0 + ",
           "des$n1 + 1")
    }
  } else {
    if (any(is.null(des$n0), length(des$n0) != 2, !is.numeric(des$n0),
            des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 2, des$n0 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 2, !is.numeric(des$n1),
            des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 2, des$n1 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$e1), length(des$e1) != des$n0[1] + des$n1[1] + 1,
            !is.numeric(des$e1), des$e1 == -Inf)) {
      stop("For des$J = 2, des$e1 must be a numeric vector of length des$n0[1]",
           " + des$n1[1] + 1, with all elements strictly greater than -Inf")
    }
    if (any(is.null(des$f1), length(des$f1) != des$n0[1] + des$n1[1] + 1,
            !is.numeric(des$f1), des$f1 >= des$e1)) {
      stop("For des$J = 2, des$f1 must be a numeric vector of length des$n0[1]",
           " + des$n1[1] + 1, with des$e1 >= des$f1 for all elements")
    }
    if (any(is.null(des$e2), nrow(des$e2) != des$n0[1] + des$n1[1] + 1,
            ncol(des$e2) != des$n0[2] + des$n1[2] + 1, !is.matrix(des$e2),
            !is.numeric(des$e2))) {
      stop("For des$J = 2, des$e2 must be a nnumeric matrix with des$n0[1] + ",
           "des$n1[1] + 1 rows and des$n0[2] + des$n1[2] + 1 columns")
    }
  }
}

check_des_single_double      <- function(des) {
  if (any(!("ph2rand_single_double" %in% class(des)), !is.list(des))) {
    stop("des must be of class ph2rand_single_double and of class list")
  } else if (any(is.null(des$J), !(des$J %in% c(1, 2)))) {
    stop("des$J must be equal to 1 or 2")
  }
  if (des$J == 1) {
    if (any(is.null(des$n0), length(des$n0) != 1, !is.numeric(des$n0),
            !is.finite(des$n0), des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 1, des$n0 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 1, !is.numeric(des$n1),
            !is.finite(des$n1), des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 1, des$n1 must be a single integer in the range ",
           "[1, \u221E)")
    }
    if (any(is.null(des$eS), length(des$eS) != 1, !is.numeric(des$eS))) {
      stop("For des$J = 1, des$eS must be a single numeric")
    } else if (any(is.infinite(des$eS), des$eS <= 0, des$eS > des$n1)) {
      warning("For des$J = 1, des$e should typically be a single integer in ",
              "(0, des$n1]")
    }
    if (any(is.null(des$eT), length(des$eT) != 1, !is.numeric(des$eT))) {
      stop("For des$J = 1, des$eT must be a single numeric")
    } else if (any(is.infinite(des$eT), des$eT%%1 != 0, des$eT <= -des$n0,
                   des$eT > des$n1)) {
      warning("For des$J = 1, des$eT should typically be a single integer in ",
              "(-des$n0, des$n1]")
    }
  } else {
    if (any(is.null(des$n0), length(des$n0) != 2, !is.numeric(des$n0),
            des$n0%%1 != 0, des$n0 < 1)) {
      stop("For des$J = 2, des$n0 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$n1), length(des$n1) != 2, !is.numeric(des$n1),
            des$n1%%1 != 0, des$n1 < 1)) {
      stop("For des$J = 2, des$n1 must be a numeric vector of length two,", 
           " containing integers in the range [1, \u221E)")
    }
    if (any(is.null(des$eS1), length(des$eS1) != 1, !is.numeric(des$eS1),
            des$eS1 == -Inf)) {
      stop("For des$J = 2, des$eS1 must be a single numeric not equal to -Inf")
    }
    if (any(is.null(des$fS1), length(des$fS1) != 1, !is.numeric(des$fS1),
            des$fS1 >= des$eS1)) {
      stop("For des$J = 2, des$fS1 must be a single numeric that is strictly ",
           "less than des$eS1")
    }
    if (any(is.null(des$eT1), length(des$eT1) != 1, !is.numeric(des$eT1),
            des$eT1 == -Inf)) {
      stop("For des$J = 2, des$eT1 must be a single numeric not equal to -Inf")
    }
    if (any(is.null(des$fT1), length(des$fT1) != 1, !is.numeric(des$fT1),
            des$fT1 >= des$eT1)) {
      stop("For des$J = 2, des$fT1 must be a single numeric that is strictly ",
           "less than des$eT1")
    }
    if (any(is.null(des$eS2), length(des$eS2) != 1, !is.numeric(des$eS2))) {
      stop("For des$J = 2, des$eS2 must be a single numeric")
    }
    if (any(is.null(des$eT2), length(des$eT2) != 1, !is.numeric(des$eT2))) {
      stop("For des$J = 2, des$eT2 must be a single numeric")
    }
    if (all(is.finite(des$eS1), is.finite(des$fS1))) {
      if (any(des$eS1%%1 != 0, des$eS1 <= 1, des$eS1 > des$n1[1])) {
        warning("For des$J = 2, when des$eS1 and des$fS1 are finite, des$eS1 ",
                "should typically be a single integer in (1, des$n1[1]]")
      }
      if (any(des$fS1%%1 != 0, des$fS1 < 0, des$fS1 >= des$n1[1] - 1)) {
        warning("For des$J = 2, when des$eS1 and des$fS1 are finite, des$fS1 ",
                "should typically be a single integer in [0, des$n1[1] - 1)")
      }
      if (any(is.infinite(des$eS2), des$eS2%%1 != 0, des$eS2 < des$fS1 + 2,
              des$eS2 > des$eS1 - 1 + des$n1[2])) {
        warning("For des$J = 2, when des$eS1 and des$fS1 are finite, des$eS2 ",
                "should typically be a single integer in ",
                "[des$fS1 + 2, des$eS1 - 1 + des$n1[2]]")
      }
    } else if (is.finite(des$eS1)) {
      if (any(des$eS1%%1 != 0, des$eS1 <= 0, des$e1 > des$n1[1])) {
        warning("For des$J = 2, when des$eS1 is finite and des$fS1 is not, ",
                "des$eS1 should typically be a single integer in (0, ",
                "des$n1[1]]")
      }
      if (any(is.infinite(des$eS2), des$eS2%%1 != 0,
              des$eS2 <= 0,
              des$eS2 > des$e1 - 1 + des$n1[2])) {
        warning("For des$J = 2, when des$eS1 is finite and des$fS1 is not, ",
                "des$eS2 should typically be a single integer in ",
                "(0, des$e1 - 1 + des$n1[2]]")
      }
    } else if (is.finite(des$fS1)) {
      if (any(des$fS1%%1 != 0, des$fS1 < 0, des$fS1 >= des$n1[1])) {
        warning("For des$J = 2, when des$fS1 is finite and des$eS1 is not, ",
                "des$fS1 should typically be a single integer in [0, ",
                "des$n1[1] - 1]")
      }
      if (any(is.infinite(des$eS2), des$eS2%%1 != 0,
              des$eS2 <= des$f1 + 1,
              des$eS2 > des$n1[1] + des$n1[2])) {
        warning("For des$J = 2, when des$fS1 is finite and des$eS1 is not, ",
                "des$eS2 should typically be a single integer in ",
                "(des$fS1 + 1, des$n1[1] + des$n1[2]]")
      }
    }
    if (all(is.finite(des$eT1), is.finite(des$fT1))) {
      if (any(des$eT1%%1 != 0, des$eT1 <= -des$n0[1] + 1,
              des$eT1 > des$n1[1])) {
        warning("For des$J = 2, when des$eT1 and des$fT1 are finite, des$eT1 ",
                "should typically be a single integer in (-des$n0[1] + 1, ",
                "des$n1[1]]")
      }
      if (any(des$fT1%%1 != 0, des$fT1 < -des$n0[1], des$fT1 >= des$n1[1] - 1)) {
        warning("For des$J = 2, when des$eT1 and des$fT1 are finite, des$fT1 ",
                "should typically be a single integer in [-des$n0[1], ",
                "des$n1[1] - 1)")
      }
      if (any(is.infinite(des$eT2), des$eT2%%1 != 0,
              des$eT2 < des$fT1 + 2 - des$n0[2],
              des$eT2 > des$eT1 - 1 + des$n1[2])) {
        warning("For des$J = 2, when des$eT1 and des$fT1 are finite, des$eT2 ",
                "should typically be a single integer in ",
                "[des$fT1 + 2 - des$n0[2], des$eT1 - 1 + des$n1[2]]")
      }
    } else if (is.finite(des$eT1)) {
      if (any(des$eT1%%1 != 0, des$eT1 <= -des$n0[1], des$eT1 > des$n1[1])) {
        warning("For des$J = 2, when des$eT1 is finite and des$fT1 is not, ",
                "des$eT1 should typically be a single integer in (-des$n0[1], ",
                "des$n1[1]]")
      }
      if (any(is.infinite(des$eT2), des$eT2%%1 != 0,
              des$eT2 <= -des$n0[1] - des$n0[2],
              des$eT2 > des$eT1 - 1 + des$n1[2])) {
        warning("For des$J = 2, when des$eT1 is finite and des$fT1 is not, ",
                "des$eT2 should typically be a single integer in ",
                "(-des$n0[1] - des$n0[2], des$eT1 - 1 + des$n1[2]]")
      }
    } else if (is.finite(des$fT1)) {
      if (any(des$fT1%%1 != 0, des$fT1 <= -des$n0[1], des$fT1 >= des$n1[1])) {
        warning("For des$J = 2, when des$fT1 is finite and des$eT1 is not, ",
                "des$fT1 should typically be a single integer in [-des$n0[1], ",
                "des$n1[1])")
      }
      if (any(is.infinite(des$eT2), des$eT2%%1 != 0,
              des$eT2 <= des$fT1 + 2 - des$n0[2],
              des$eT2 > des$n1[1] + des$n1[2])) {
        warning("For des$J = 2, when des$fT1 is finite and des$eT1 is not, ",
                "des$eT2 should typically be a single integer in ",
                "(des$fT1 + 2 - des$n0[2], des$n1[1] + des$n1[2])")
      }
    }
  }
}

check_integer_range          <- function(value, name, range, len) {
  check         <- F
  if (is.finite(len)) {
    if (any(length(value) != len, !is.numeric(value), value%%1 != 0,
            value <= range[1], value >= range[2])) {
      check     <- T
      if (len == 1) {
        segment <- " a single integer that belongs to {"
      } else {
        segment <- paste(" an integer vector of length", len, "whose elements",
                         "all belong to {")
      }
    }
  } else if (any(value%%1 != 0, !is.numeric(value), value <= range[1],
                 value >= range[2])) {
    check       <- T  
    segment     <- " an integer vector whose elements all belong to {"
  }
  if (check) {
    if (range[1] + 2 == range[2]) {
      stop(name, " must be equal to ", range[1] + 1)
    } else if (range[1] + 3 == range[2]) {
      stop(name, segment, range[1] + 1, ", ", range[1] + 2, "}")
    } else if (range[1] + 4 == range[2]) {
      stop(name, segment, range[1] + 1, ", ", range[1] + 2, ", ", range[1] + 3,
           "}")
    } else if (all(is.infinite(range))) {
      stop(name, segment, "..., -1, 0, 1, ...}")
    } else if (is.infinite(range[1])) {
      stop(name, segment, "..., ", range[2] - 2, ", ..., ", range[2] - 1, "}")
    } else if (is.infinite(range[2])) {
      stop(name, segment, range[1] + 1, ", ", range[1] + 2, ", ...}")
    } else {
      stop(name, segment, range[1] + 1, ", ..., ", range[2] - 1, "}")
    }
  }
  return(as.integer(value))
}

check_k                      <- function(k, des) {
  if (missing(k)) {
    if ("ph2rand_des_one_stage" %in% class(des)) {
      return(1)
    } else {
      return(1:2)
    }
  } else {
    if ("ph2rand_des_one_stage" %in% class(des)) {
      allowed <- 1
    } else {
      allowed <- 1:2
    }
    if (any(!(k %in% allowed))) {
      stop("k must contain values in {", paste(allowed, sep = ", "), "}")
    }
  }
}

check_logical                <- function(value, name) {
  if (!is.logical(value)) {
    stop(name, " must be a logical variable")
  }
}

check_fisher_params          <- function(efficacy_type, efficacy_param,
                                         futility_type, futility_param) {
  if (futility_type == 0) {
    if (futility_param != 0) {
      warning("futility_param has been changed from default, but this will ",
              "have no effect given the value of futility_type")
    }
  } else if (futility_type == 1) {
    if (any(length(futility_param) != 1, futility_param%%1 != 0)) {
      stop("For futility_type = 1, futility_param must be a single integer")
    }
  } else {
    if (any(length(futility_param) != 1, futility_param <= 0,
            futility_param >= 1)) {
      stop("For futility_type = 2, futility_param must be a single numeric in ",
           "(0, 1)")
    }
  }
  check              <- F
  if (is.null(efficacy_param)) {
    efficacy_param   <- -0.5
    check            <- T
  }
  if (efficacy_type == 0) {
    if (!check) {
      warning("efficacy_param has been changed from default, but this will ",
              "have no effect given the value of efficacy_type")
    }
  } else if (efficacy_type == 1) {
    if (!check) {
      if (any(length(efficacy_param) != 1, efficacy_param%%1 != 0)) {
        stop("For efficacy_type = 1, efficacy_param must be a single integer")
      }
    }
  } else {
    if (any(length(efficacy_param) != 1, efficacy_param <= 0,
            efficacy_param >= 1)) {
      stop("For efficacy_type = 2, efficacy_param must be a single numeric in ",
           "(0, 1)")
    }
  }
  if (all(efficacy_type == 1, futility_type == 1,
          efficacy_param <= futility_param, efficacy_param != -0.5, !check)) {
    stop("When efficacy_type and futility_type are both equal to 1, futility_",
         "param must be strictly less than efficacy_param")
  }
  return(efficacy_param)
}

check_pi                     <- function(pi, des) {
  if (missing(pi)) {
    pi           <- as.matrix(des$opchar[, 1:2])
    colnames(pi) <- NULL
  } else {
    if (!is.numeric(pi)) {
      stop("pi must be either a numeric vector of length two, or a numeric ",
           "matrix with two columns. In either case all elements must take ",
           "values in [0, 1]")
    } else {
      if (is.matrix(pi)) {
        if (any(ncol(pi) != 2, pi < 0, pi > 1)) {
          stop("pi must be either a numeric vector of length two, or a numeric",
               " matrix with two columns. In either case all elements must ",
               "take values in [0, 1]")
        } else if (sum(duplicated(pi)) > 0) {
          warning("pi contains duplicated rows")
        }
      } else {
        if (any(length(pi) != 2, pi < 0, pi > 1)) {
          stop("pi must be either a numeric vector of length two, or a numeric",
               " matrix with two columns. In either case all elements must ",
               "take values in [0, 1]")
        } else {
          pi     <- matrix(pi, 1)
        }
      }
    }
  }
  pi
}

check_pi0_alt                <- function(pi0_alt, delta) {
  if (!(length(pi0_alt) %in% c(1, 2))) {
    stop("pi0_alt must be a numeric vector of length 1 or 2")
  }
  if (length(pi0_alt) == 1) {
    if (any(pi0_alt < 0, pi0_alt + delta > 1)) {
      stop("If pi0_alt is a numeric vector of length 1, it must belong to ",
           "[0, 1 - delta]")
    }
  } else {
    if (any(pi0_alt < 0, pi0_alt > 1 - delta, pi0_alt[2] <= pi0_alt[1])) {
      stop("If pi0_alt is a numeric vector of length 2, both elements must ",
           "belong to [0, 1 - delta], and the second element must be strictly ",
           "larger than the first")
    }
  }
}

check_pi0_null               <- function(pi0_null, point_null) {
  if (!(length(pi0_null) %in% c(1, 2))) {
    stop("pi0_null must be a numeric vector of length 1 or 2")
  }
  if (length(pi0_null) == 1) {
    if (any(pi0_null < 0, pi0_null > 1)) {
      stop("If pi0_alt is a numeric vector of length 1, it must belong to ",
           "[0, 1]")
    }
  } else {
    if (any(pi0_null < 0, pi0_null > 1, pi0_null[2] <= pi0_null[1])) {
      stop("If pi0_alt is a numeric vector of length 2, both elements must ",
           "belong to [0, 1], and the second element must be strictly larger ",
           "than the first")
    }
  }
}

check_real_range             <- function(value, name, range, len) {
  if (is.finite(len)) {
    if (any(length(value) != len, !is.numeric(value), value < range[1],
            value > range[2])) {
      if (len == 1) {
        stop(name, " must be a single numeric that belongs to [", range[1], ",",
             range[2], "]")
      } else {
        stop(name, " must be a numeric vector of length ", len, ", whose ",
             "elements all belong to [", range[1], ",", range[2], "]")
      }
    }
  } else {
    if (any(!is.numeric(value), value < range[1], value > range[2])) {
      stop(name, " must be a numeric vector whose elements all belong to [",
           range[1], ",", range[2], "]")
    }
  }
}

check_real_range_strict      <- function(value, name, range, len) {
  if (is.finite(len)) {
    if (any(length(value) != len, !is.numeric(value), value <= range[1],
            value >= range[2])) {
      if (len == 1) {
        stop(name, " must be a single numeric that belongs to (", range[1], ",",
             range[2], ")")
      } else {
        stop(name, " must be a numeric vector of length ", len, ", whose ",
             "elements all belong to (", range[1], ",", range[2], ")")
      }
    }
  } else {
    if (any(!is.numeric(value), value <= range[1], value >= range[2])) {
      stop(name, " must be a numeric vector whose elements all belong to (",
           range[1], ",", range[2], ")")
    }
  }
}

check_w                      <- function(w) {
  if (any(length(w) != 5, !is.numeric(w), w < 0, sum(w[-5]) == 0)) {
    stop("w must be a numeric vector of length 5, with all elements greater ",
         "than or equal to 0, and with at least one of the first 4 elements ",
         "strictly positive")
  }
}

check_real_pair_range_strict <- function(value1, value2, name1, name2, range){
  if (any(length(value1) != 1, value1 <= range[1], value1 >= range[2])){
    stop(name1, " must be a single number in (", range[1], ",", range[2], ")")
  }
  if (any(length(value2) != 1, value2 <= range[1], value2 >= range[2])){
    stop(name2, " must be a single number in (", range[1], ",", range[2], ")")
  }
  if (value1 >= value2){
    stop(name1, " must be strictly less than ", name2)
  }
}

check_integer_pair_range     <- function(value1, value2, name1, name2, range){
  if (any(length(value1) != 1, value1 <= range[1], value1 >= range[2])){
    stop(name1, " must be a single integer in (", range[1], ",", range[2], ")")
  }
  if (any(length(value2) != 1, value2 <= range[1], value2 >= range[2])){
    stop(name2, " must be a single integer in (", range[1], ",", range[2], ")")
  }
  if (value1 > value2){
    stop(name1, " must be less than or equal to ", name2)
  }
}

check_stopping               <- function(futility, efficacy){
  if (!is.logical(futility)){
    stop("futility must be a logical variable")
  }
  if (!is.logical(efficacy)){
    stop("efficacy must be a logical variable")
  }
  if (all(!futility, !efficacy)){
    stop("At least one of futility and efficacy must be set to TRUE")
  }
}

check_optimality             <- function(optimality, point_prior, beta_prior,
                                         futility, efficacy){
  if (!(optimality %in% c("minimax", "null_ess", "alt_ess", "null_med",
                          "alt_med", "prior", "admissable"))){
    stop("optimality must be set to one of \"minimax\", \"null_ess\", \"alt_ess\", \"null_med\", \"alt_med\", \"prior\" or \"admissable\"")
  }
  if (optimality == "prior"){
    if (all(!is.null(point_prior), !is.null(beta_prior))){
      stop("point_prior and beta_prior cannot both be specified when optimality is set to \"prior\"")
    }
    if (all(is.null(point_prior), is.null(beta_prior))){
      stop("One of point_prior and beta_prior must be specified when optimality is set to \"prior\"")
    }
    if (!is.null(point_prior)){
      if (any(length(point_prior) != 1, !is.finite(point_prior),
              point_prior < 0, point_prior > 1)){
        stop("point_prior must be a single number in [0,1]")
      }
    }
    if (!is.null(beta_prior)){
      if (any(length(beta_prior) != 2, any(!is.finite(beta_prior)),
              any(beta_prior <= 0))){
        stop("beta_prior must be a vector of length 2, containing values in (0,Inf)")
      }
    }
  } else {
    if (!is.null(point_prior)){
      warning("point_prior is not missing. This will have no effect as optimality is not set to \"prior\"")
    }
    if (!is.null(beta_prior)){
      warning("beta_prior is not missing. This will have no effect as optimality is not set to \"prior\"")
    }
  }
  if (all(any(!futility, !efficacy), optimality == "delta_minimax")){
    warning("If futility or efficacy are set to FALSE the delta-minimax design is equivalent to the minimax design")
  }
}

check_ensign                 <- function(ensign, futility){
  if (!is.logical(ensign)){
    stop("ensign must be a logical variable")
  }
  if (all(!futility, ensign)){
    stop("ensign cannot be set to TRUE and futility to FALSE")
  }
}

check_gs_boundaries          <- function(J, n, a, r) {
  if (any(!is.numeric(n), length(n) != J, n%%1 != 0)) {
    stop("n must be a numeric vector of length ", J,
         ", containing integer elements")
  }
  if (any(!is.numeric(a), length(a) != J)) {
    stop("a must be a numeric vector of length ", J)
  }
  if (any(!is.numeric(r), length(r) != J)) {
    stop("r must be a numeric vector of length ", J)
  }
  if (J > 1) {
    for (j in 1:(J - 1)) {
      if (all(is.finite(a[j]), a[j]%%1 != 0)) {
        stop("Finite elements of a must be integers")
      }
      if (all(is.finite(r[j]), r[j]%%1 != 0)) {
        stop("Finite elements of r must be integers")
      }
    }
  }
  if (!is.finite(a[J])) {
    stop("a[J] must be finite")
  }
  if (!is.finite(r[J])) {
    stop("r[J] must be finite")
  }
  if (a[J] != r[J] - 1) {
    stop("a[J] must be one less than r[J]")
  }
}