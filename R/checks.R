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

check_integer_pair_range     <- function(value1, value2, name1, name2, range) {
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

check_integer_range          <- function(value, name, range, len) {
  check         <- FALSE
  if (is.finite(len)) {
    if (any(length(value) != len, !is.numeric(value), value%%1 != 0,
            value <= range[1], value >= range[2])) {
      check     <- TRUE
      if (len == 1) {
        segment <- " a single integer that belongs to {"
      } else {
        segment <- paste(" an integer vector of length", len, "whose elements",
                         "all belong to {")
      }
    }
  } else if (any(value%%1 != 0, !is.numeric(value), value <= range[1],
                 value >= range[2])) {
    check       <- TRUE  
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
  if (any(!(k %in% 1:des$J))) {
    stop("k must contain values in {", paste(1:des$J, sep = ", "), "}")
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
  check              <- FALSE
  if (is.null(efficacy_param)) {
    efficacy_param   <- -0.5
    check            <- TRUE
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

check_ph2rand_des            <- function(des, type, name = "des") {
  if (!("ph2rand_des" %in% class(des))) {
    stop(name, " must be of class ph2rand_des")
  }
  if (all(type == "any", !(des$type %in% c("barnard", "binomial", "fisher",
                                           "sat")))) {
    stop(name, "$type must be one of \"barnard\", \"binomial\", \"fisher\", or",
         " \"sat\"")
  } else if (all(type == "barnard", des$type != "barnard")) {
    stop(name, "$type must be equal to \"barnard\"")
  } else if (all(type == "binomial", des$type != "binomial")) {
    stop(name, "$type must be equal to \"binomial\"")
  } else if (all(type == "fisher", des$type != "fisher")) {
    stop(name, "$type must be equal to \"fisher\"")
  } else if (all(type == "sat", des$type != "sat")) {
    stop(name, "$type must be equal to \"sat\"")
  }
  if (any(is.null(des$J), !(des$J %in% 1:2))) {
    stop(name, "$J must be equal to 1 or 2")
  } else if (des$J == 1) {
    if (any(is.null(des$nC), length(des$nC) != 1, !is.numeric(des$nC),
            !is.finite(des$nC), des$nC%%1 != 0, des$nC < 1)) {
      stop("For ", name, "$J = 1, ", name, "$nC must be a single integer in ",
           "{1, 2, 3, ...}")
    }
    if (any(is.null(des$nE), length(des$nE) != 1, !is.numeric(des$nE),
            !is.finite(des$nE), des$nE%%1 != 0, des$nE < 1)) {
      stop("For ", name, "$J = 1, ", name, "$nE must be a single integer in ",
           "{1, 2, 3, ...}")
    }
  } else {
    if (any(is.null(des$nC), length(des$nC) != 2, !is.numeric(des$nC),
            des$nC%%1 != 0, des$nC < 1)) {
      stop("For ", name, "$J = 2, ", name, "$nC must be a numeric vector of ", 
           "length two, whose elements are integers in {1, 2, 3, ...}")
    }
    if (any(is.null(des$nE), length(des$nE) != 2, !is.numeric(des$nE),
            des$nE%%1 != 0, des$nE < 1)) {
      stop("For ", name, "$J = 2, ", name, "$nE must be a numeric vector of ", 
           "length two, whose elements are integers in {1, 2, 3, ...}")
    }
  }
  if (des$type == "barnard") {
    if (des$J == 1) {
      if (any(is.null(des$boundaries$e1), length(des$boundaries$e1) != 1,
              !is.numeric(des$boundaries$e1))) {
        stop("For ", name, "$J = 1, ", name, "$boundaries$e1 must be a single ",
             "numeric")
      } else if (is.infinite(des$boundaries$e1)) {
        warning("For ", name, "$J = 1, ", name, "$boundaries$e1 should ",
                "typically be finite")
      }
    } else {
      if (any(is.null(des$boundaries$e1), length(des$boundaries$e1) != 1,
              !is.numeric(des$boundaries$e1), des$boundaries$e1 == -Inf)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$e1 must be a single ",
             "numeric not equal to -Inf")
      }
      if (any(is.null(des$boundaries$f1), length(des$boundaries$f1) != 1,
              !is.numeric(des$boundaries$f1),
              des$boundaries$f1 >= des$boundaries$e1)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$f1 must be a single ",
             "numeric that is strictly less than ", name, "$boundaries$e1")
      }
      if (any(is.null(des$boundaries$e2), length(des$boundaries$e2) != 1,
              !is.numeric(des$boundaries$e2))) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$e2 must be a single ",
             "numeric")
      } else if (is.infinite(des$boundaries$e2)) {
        warning("For ", name, "$J = 2, ", name, "$boundaries$e2 should ",
                "typically be finite")
      }
    }
  } else if (des$type == "binomial") {
    if (des$J == 1) {
      if (any(is.null(des$boundaries$e1), length(des$boundaries$e1) != 1,
              !is.numeric(des$boundaries$e1))) {
        stop("For ", name, "$J = 1, ", name, "$boundaries$e1 must be a single ",
             "numeric")
      } else if (any(is.infinite(des$boundaries$e1), des$boundaries$e1%%1 != 0,
                     des$boundaries$e1 <= -des$nC,
                     des$boundaries$e1 > des$nE)) {
        warning("For ", name, "$J = 1, ", name, "$boundaries$e1 should ",
                "typically be a single integer in (-", name, "$nC, ", name,
                "$nE]")
      }
    } else {
      if (any(is.null(des$boundaries$e1), length(des$boundaries$e1) != 1,
              !is.numeric(des$boundaries$e1), des$boundaries$e1 == -Inf)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$e1 must be a single ",
             "numeric not equal to -Inf")
      }
      if (any(is.null(des$boundaries$f1), length(des$boundaries$f1) != 1,
              !is.numeric(des$boundaries$f1),
              des$boundaries$f1 >= des$boundaries$e1)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$f1 must be a single ",
             "numeric that is strictly less than ", name, "$boundaries$e1")
      }
      if (any(is.null(des$boundaries$e2), length(des$boundaries$e2) != 1,
              !is.numeric(des$boundaries$e2))) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$e2 must be a single ",
             "numeric")
      }
      if (all(is.finite(des$boundaries$e1), is.finite(des$boundaries$f1))) {
        if (any(des$boundaries$e1%%1 != 0, des$boundaries$e1 <= -des$nC[1] + 1,
                des$boundaries$e1 > des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$e1 and ",
                  name, "$boundaries$f1 are finite, ", name, "$boundaries$e1 ",
                  "should typically be a single integer in (-", name,
                  "$nC[1] + 1, ", name, "$nE[1]]")
        }
        if (any(des$boundaries$f1%%1 != 0, des$boundaries$f1 < -des$nC[1],
                des$boundaries$f1 >= des$nE[1] - 1)) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$e1 and ",
                  name, "$boundaries$f1 are finite, ", name, "$boundaries$f1 ",
                  "should typically be a single integer in [-", name,
                  "$nC[1], ", name, "$nE[1] - 1)")
        }
        if (any(is.infinite(des$boundaries$e2), des$boundaries$e2%%1 != 0,
                des$boundaries$e2 < des$boundaries$f1 + 2 - des$nC[2],
                des$boundaries$e2 > des$boundaries$e1 - 1 + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$e1 and ",
                  name, "$boundaries$f1 are finite, ", name, "$boundaries$e2 ",
                  "should typically be a single integer in [", name,
                  "$boundaries$f1 + 2 - ", name, "$nC[2], ", name,
                  "$boundaries$e1 - 1 + ", name, "$nE[2]]")
        }
      } else if (is.finite(des$boundaries$e1)) {
        if (any(des$boundaries$e1%%1 != 0, des$boundaries$e1 <= -des$nC[1],
                des$boundaries$e1 > des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$e1 is ",
                  "finite and ", name, "$boundaries$f1 is not, ", name,
                  "$boundaries$e1 should typically be a single integer in (-",
                  name, "$nC[1], ", name, "$nE[1]]")
        }
        if (any(is.infinite(des$boundaries$e2), des$boundaries$e2%%1 != 0,
                des$boundaries$e2 <= -des$nC[1] - des$nC[2],
                des$boundaries$e2 > des$boundaries$e1 - 1 + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$e1 is ",
                  "finite and ", name, "$boundaries$f1 is not, ", name,
                  "$boundaries$e2 should typically be a single integer in (-",
                  name, "$nC[1] - ", name, "$nC[2], ",
                  name, "$boundaries$e1 - 1 + ", name, "$nE[2]]")
        }
      } else if (is.finite(des$boundaries$f1)) {
        if (any(des$boundaries$f1%%1 != 0, des$boundaries$f1 <= -des$nC[1],
                des$boundaries$f1 >= des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$f1 is ",
                  "finite and ", name, "$boundaries$e1 is not, ", name,
                  "$boundaries$f1 should typically be a single integer in [-",
                  name, "$nC[1], ", name, "$nE[1])")
        }
        if (any(is.infinite(des$boundaries$e2), des$boundaries$e2%%1 != 0,
                des$boundaries$e2 <= des$boundaries$f1 + 2 - des$nC[2],
                des$boundaries$e2 > des$nE[1] + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$f1 is ",
                  "finite and ", name, "$boundaries$e1 is not, ", name,
                  "$boundaries$e2 should typically be a single integer in (",
                  name, "$boundaries$f1 + 2 - ", name, "$nC[2], ", name,
                  "$nE[1] + ", name, "$nE[2])")
        }
      }
    }
  } else if (des$type == "fisher") {
    if (des$J == 1) {
      if (any(is.null(des$boundaries$e1),
              length(des$boundaries$e1) != des$nC + des$nE + 1,
              !is.numeric(des$boundaries$e1))) {
        stop("For ", name, "$J = 1, ", name, "$boundaries$e1 must be a numeric",
             "vector of length ", name, "$boundaries$nC + ", name,
             "$boundaries$nE + 1")
      }
    } else {
      if (any(is.null(des$boundaries$e1),
              length(des$boundaries$e1) != des$nC[1] + des$nE[1] + 1,
              !is.numeric(des$boundaries$e1), des$boundaries$e1 == -Inf)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$e1 must be a numeric",
             " vector of length ", name, "$nC[1] + ", name, "$nE[1] + 1, with ",
             "all elements strictly greater than -Inf")
      }
      if (any(is.null(des$boundaries$f1),
              length(des$boundaries$f1) != des$nC[1] + des$nE[1] + 1,
              !is.numeric(des$boundaries$f1),
              des$boundaries$f1 >= des$boundaries$e1)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$f1 must be a numeric",
             "vector of length ", name, "$nC[1] + ", name, "$nE[1] + 1, with ",
             name, "$boundaries$e1 >= ", name, "$boundaries$f1 for all ",
             "elements")
      }
      if (any(is.null(des$boundaries$e2),
              nrow(des$boundaries$e2) != des$nC[1] + des$nE[1] + 1,
              ncol(des$boundaries$e2) != des$nC[2] + des$nE[2] + 1,
              !is.matrix(des$boundaries$e2), !is.numeric(des$boundaries$e2))) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$e2 must be a numeric",
             " matrix with ", name, "$nC[1] + ", name, "$nE[1] + 1 rows and ",
             name, "$nC[2] + ", name, "$nE[2] + 1 columns")
      }
    }
  } else if (des$type == "sat") {
    if (des$J == 1) {
      if (any(is.null(des$boundaries$eS1), length(des$boundaries$eS1) != 1,
              !is.numeric(des$boundaries$eS1))) {
        stop("For ", name, "$J = 1, ", name, "$boundaries$eS1 must be a single",
             "numeric")
      } else if (any(is.infinite(des$boundaries$eS1), des$boundaries$eS1 <= 0,
                     des$boundaries$eS1 > des$nE)) {
        warning("For ", name, "$J = 1, ", name, "$boundaries$eS1 should ",
                "typically be a single integer in (0, ", name, "$nE]")
      }
      if (any(is.null(des$boundaries$eT1), length(des$boundaries$eT1) != 1,
              !is.numeric(des$boundaries$eT1))) {
        stop("For ", name, "$J = 1, ", name, "$boundaries$eT1 must be a single",
             " numeric")
      } else if (any(is.infinite(des$boundaries$eT1),
                     des$boundaries$eT1%%1 != 0, des$boundaries$eT1 <= -des$nC,
                     des$boundaries$eT1 > des$nE)) {
        warning("For ", name, "$J = 1, ", name, "$boundaries$eT1 should ",
                "typically be a single integer in (-", name, "$nC, ", name,
                "$nE]")
      }
    } else {
      if (any(is.null(des$boundaries$eS1), length(des$boundaries$eS1) != 1,
              !is.numeric(des$boundaries$eS1), des$boundaries$eS1 == -Inf)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$eS1 must be a single",
             " numeric not equal to -Inf")
      }
      if (any(is.null(des$boundaries$fS1), length(des$boundaries$fS1) != 1,
              !is.numeric(des$boundaries$fS1),
              des$boundaries$fS1 >= des$boundaries$eS1)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$fS1 must be a single",
             " numeric that is strictly less than ", name, "$boundaries$eS1")
      }
      if (any(is.null(des$boundaries$eT1), length(des$boundaries$eT1) != 1,
              !is.numeric(des$boundaries$eT1), des$boundaries$eT1 == -Inf)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$eT1 must be a single",
             " numeric not equal to -Inf")
      }
      if (any(is.null(des$boundaries$fT1), length(des$boundaries$fT1) != 1,
              !is.numeric(des$boundaries$fT1),
              des$boundaries$fT1 >= des$boundaries$eT1)) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$fT1 must be a single",
             " numeric that is strictly less than ", name, "$boundaries$eT1")
      }
      if (any(is.null(des$boundaries$eS2), length(des$boundaries$eS2) != 1,
              !is.numeric(des$boundaries$eS2))) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$eS2 must be a single",
             " numeric")
      }
      if (any(is.null(des$boundaries$eT2), length(des$boundaries$eT2) != 1,
              !is.numeric(des$boundaries$eT2))) {
        stop("For ", name, "$J = 2, ", name, "$boundaries$eT2 must be a single",
             " numeric")
      }
      if (all(is.finite(des$boundaries$eS1), is.finite(des$boundaries$fS1))) {
        if (any(des$boundaries$eS1%%1 != 0, des$boundaries$eS1 <= 1,
                des$boundaries$eS1 > des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eS1 and ",
                  name, "$boundaries$fS1 are finite, ", name, "$boundaries$eS1",
                  " should typically be a single integer in (1, ", name,
                  "$nE[1]]")
        }
        if (any(des$boundaries$fS1%%1 != 0, des$boundaries$fS1 < 0,
                des$boundaries$fS1 >= des$nE[1] - 1)) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eS1 and ",
                  name, "$boundaries$fS1 are finite, ", name, "$boundaries$fS1",
                  " should typically be a single integer in [0, ", name,
                  "$nE[1] - 1)")
        }
        if (any(is.infinite(des$boundaries$eS2), des$boundaries$eS2%%1 != 0,
                des$boundaries$eS2 < des$boundaries$fS1 + 2,
                des$boundaries$eS2 > des$boundaries$eS1 - 1 + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eS1 and ",
                  name, "$boundaries$fS1 are finite, ", name, "$boundaries$eS2",
                  " should typically be a single integer in [", name,
                  "$boundaries$fS1 + 2, ", name, "$boundaries$eS1 - 1 + ", name,
                  "$nE[2]]")
        }
      } else if (is.finite(des$boundaries$eS1)) {
        if (any(des$boundaries$eS1%%1 != 0, des$boundaries$eS1 <= 0,
                des$boundaries$eS1 > des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eS1 is ",
                  "finite and ", name, "$boundaries$fS1 is not, ", name,
                  "$boundaries$eS1 should typically be a single integer in ",
                  "(0, ", name, "$nE[1]]")
        }
        if (any(is.infinite(des$boundaries$eS2), des$boundaries$eS2%%1 != 0,
                des$boundaries$eS2 <= 0,
                des$boundaries$eS2 > des$boundaries$eS1 - 1 + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eS1 is ",
                  "finite and ", name, "$boundaries$fS1 is not, ", name,
                  "$boundaries$eS2 should typically be a single integer in ",
                  "(0, ", name, "$boundaries$eS1 - 1 + ", name, "$nE[2]]")
        }
      } else if (is.finite(des$boundaries$fS1)) {
        if (any(des$boundaries$fS1%%1 != 0, des$boundaries$fS1 < 0,
                des$boundaries$fS1 >= des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$fS1 is ",
                  "finite and ", name, "$boundaries$eS1 is not, ", name,
                  "$boundaries$fS1 should typically be a single integer in ",
                  "[0, ", name, "$nE[1] - 1]")
        }
        if (any(is.infinite(des$boundaries$eS2), des$boundaries$eS2%%1 != 0,
                des$boundaries$eS2 <= des$boundaries$fS1 + 1,
                des$boundaries$eS2 > des$nE[1] + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$fS1 is ",
                  "finite and ", name, "$boundaries$eS1 is not, ", name,
                  "$boundaries$eS2 should typically be a single integer in (",
                  name, "$boundaries$fS1 + 1, ", name, "$nE[1] + ", name,
                  "$nE[2]]")
        }
      }
      if (all(is.finite(des$boundaries$eT1), is.finite(des$boundaries$fT1))) {
        if (any(des$boundaries$eT1%%1 != 0,
                des$boundaries$eT1 <= -des$nC[1] + 1,
                des$boundaries$eT1 > des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eT1 and ",
                  name, "$boundaries$fT1 are finite, ", name, "$boundaries$eT1",
                  " should typically be a single integer in (-", name,
                  "$nC[1] + 1, ", name, "$nE[1]]")
        }
        if (any(des$boundaries$fT1%%1 != 0, des$boundaries$fT1 < -des$nC[1],
                des$boundaries$fT1 >= des$nE[1] - 1)) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eT1 and ",
                  name, "$boundaries$fT1 are finite, ", name, "$boundaries$fT1",
                  " should typically be a single integer in [-", name,
                  "$nC[1], ", name, "$nE[1] - 1)")
        }
        if (any(is.infinite(des$boundaries$eT2), des$boundaries$eT2%%1 != 0,
                des$boundaries$eT2 < des$boundaries$fT1 + 2 - des$nC[2],
                des$boundaries$eT2 > des$boundaries$eT1 - 1 + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eT1 and ",
                  name, "$boundaries$fT1 are finite, ", name, "$boundaries$eT2",
                  " should typically be a single integer in [", name,
                  "$boundaries$fT1 + 2 - ", name, "$nC[2], ", name,
                  "$boundaries$eT1 - 1 + ", name, "$nE[2]]")
        }
      } else if (is.finite(des$boundaries$eT1)) {
        if (any(des$boundaries$eT1%%1 != 0, des$boundaries$eT1 <= -des$nC[1],
                des$boundaries$eT1 > des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eT1 is ",
                  "finite and ", name, "$boundaries$fT1 is not, ", name,
                  "$boundaries$eT1 should typically be a single integer in (-",
                  name, "$nC[1], ", name, "$nE[1]]")
        }
        if (any(is.infinite(des$boundaries$eT2), des$boundaries$eT2%%1 != 0,
                des$boundaries$eT2 <= -des$nC[1] - des$nC[2],
                des$boundaries$eT2 > des$boundaries$eT1 - 1 + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$eT1 is ",
                  "finite and ", name, "$boundaries$fT1 is not, ", name,
                  "$boundaries$eT2 should typically be a single integer in (-",
                  name, "$nC[1] - ", name, "$nC[2], ", name,
                  "$boundaries$eT1 - 1 + ", name, "$nE[2]]")
        }
      } else if (is.finite(des$boundaries$fT1)) {
        if (any(des$boundaries$fT1%%1 != 0, des$boundaries$fT1 <= -des$nC[1],
                des$boundaries$fT1 >= des$nE[1])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$fT1 is ",
                  "finite and ", name, "$boundaries$eT1 is not, ", name,
                  "$boundaries$fT1 should typically be a ",
                  "single integer in [-", name, "$nC[1], ", name, "$nE[1])")
        }
        if (any(is.infinite(des$boundaries$eT2), des$boundaries$eT2%%1 != 0,
                des$boundaries$eT2 <= des$boundaries$fT1 + 2 - des$nC[2],
                des$boundaries$eT2 > des$nE[1] + des$nE[2])) {
          warning("For ", name, "$J = 2, when ", name, "$boundaries$fT1 is ",
                  "finite and ", name, "$boundaries$eT1 is not, ", name,
                  "$boundaries$eT2 should typically be a single integer in (",
                  name, "$boundaries$fT1 + 2 - ", name, "$nC[2], ",
                  name, "$nE[1] + ", name, "$nE[2])")
        }
      }
    }
  }
}

check_ph2rand_pmf            <- function(x) {
  check_ph2rand_des(x$des, "any", "x$des")
  if (!tibble::is_tibble(x$pmf)) {
    stop("x$pmf must be a tibble")
  }
  if (x$des$type %in% c("barnard", "binomial")) {
    if (ncol(x$pmf) != 10) {
      stop("For x$des$type = ", x$des$type, ", x$pmf must be a tibble with 10 ",
           "columns")
    }
    if (!all(colnames(x$pmf) == c("piC", "piE", "xC", "xE", "mC", "mE",
                                  "statistic", "decision", "k", "f(x,m|pi)"))) {
      stop("The column names of x$pmf are not as required when x$des$type = ",
           x$des$type)
    }
  } else if (x$des$type == "sat") {
    if (ncol(x$pmf) != 11) {
      stop("For x$des$type = \"sat\", x$pmf must be a tibble with 11",
           " columns")
    }
    if (!all(colnames(x$pmf) == c("piC", "piE", "xC", "xE", "mC", "mE",
                                  "statisticS", "statisticD", "decision", "k",
                                  "f(x,m|pi)"))) {
      stop("The column names of x$pmf are not as required when x$des$type = ",
           "\"sat\"")
    }
  } else if (all(x$des$type == "fisher", x$des$J == 1)) {
    if (ncol(x$pmf) != 11) {
      stop("For x$des$type = \"fisher\", with x$des$J = 1, x$pmf must be a ",
           "tibble with 11 columns")
    }
    if (!all(colnames(x$pmf) == c("piC", "piE", "xC", "xE", "mC", "mE", "z",
                                  "statistic", "decision", "k", "f(x,m|pi)"))) {
      stop("The column names of x$pmf are not as required when x$des$type = ",
           "\"fisher\" and x$des$J = 1")
    }
  } else {
    if (ncol(x$pmf) != 14) {
      stop("For x$des$type = \"fisher\", with x$des$J = 2, x$pmf must be a ",
           "tibble with 14 columns")
    }
    if (!all(colnames(x$pmf) == c("piC", "piE", "xC1", "xE1", "xC2", "xE2",
                                  "mC", "mE", "z1", "z2", "statistic",
                                  "decision", "k", "f(x,m|pi)"))) {
      stop("The column names of x$pmf are not as required when x$des$type = ",
           "\"fisher\" and x$des$J = 2")
    }
  }
}

check_ph2rand_terminal       <- function(x) {
  check_ph2rand_des(x$des, "any", "x$des")
  if (!tibble::is_tibble(x$terminal)) {
    stop("x$terminal must be a tibble")
  }
  if (x$des$type %in% c("barnard", "binomial")) {
    if (ncol(x$terminal) != 7) {
      stop("For x$des$type = ", x$des$type, ", x$terminal must be a tibble ",
           "with 7 columns")
    }
    if (!all(colnames(x$terminal) == c("xC", "xE", "mC", "mE", "statistic",
                                       "decision", "k"))) {
      stop("The column names of x$terminal are not as required when ",
           "x$des$type = ", x$des$type)
    }
  } else if (x$des$type == "sat") {
    if (ncol(x$terminal) != 8) {
      stop("For x$des$type = \"sat\", x$terminal must be a tibble ",
           "with 8 columns")
    }
    if (!all(colnames(x$terminal) == c("xC", "xE", "mC", "mE", "statisticS",
                                       "statisticT", "decision", "k"))) {
      stop("The column names of x$terminal are not as required when ",
           "x$des$type = \"sat\"")
    }
  } else if (all(x$des$type == "fisher", x$des$J == 1)) {
    if (ncol(x$terminal) != 8) {
      stop("For x$des$type = \"fisher\", with x$des$J = 1, x$terminal must be ",
           "a tibble with 8 columns")
    }
    if (!all(colnames(x$terminal) == c("xC", "xE", "mC", "mE", "z", "statistic",
                                       "decision", "k"))) {
      stop("The column names of x$terminal are not as required when ",
           "x$des$type = \"fisher\" and x$des$J = 1")
    }
  } else {
    if (ncol(x$terminal) != 11) {
      stop("For x$des$type = \"fisher\", with x$des$J = 2, x$terminal must be ",
           "a tibble with 11 columns")
    }
    if (!all(colnames(x$terminal) == c("xC1", "xE1", "xC2", "xE2", "mC", "mE",
                                       "z1", "z2", "statistic", "decision",
                                       "k"))) {
      stop("The column names of x$terminal are not as required when ",
           "x$des$type = \"fisher\" and x$des$J = 2")
    }
  }
}

check_pi                     <- function(pi, des) {
  if (all(!is.numeric(pi), !is.data.frame(pi))) {
    stop("pi must be either a numeric vector of length two, or a numeric ",
         "matrix/data frame with two columns. In either case all elements must",
         " take values in [0, 1]")
  } else {
    if (any(is.matrix(pi), is.data.frame(pi))) {
      if (any(ncol(pi) != 2, pi < 0, pi > 1)) {
        stop("pi must be either a numeric vector of length two, or a numeric ",
             "matrix/data frame with two columns. In either case all elements ",
             "must take values in [0, 1]")
      } else if (sum(duplicated(pi)) > 0) {
        warning("pi contains duplicated rows")
      }
      if (is.data.frame(pi)) {
        pi <- as.matrix(pi)
      }
    } else {
      if (any(length(pi) != 2, pi < 0, pi > 1)) {
        stop("pi must be either a numeric vector of length two, or a numeric ",
             "matrix/data frame with two columns. In either case all elements ",
             "must take values in [0, 1]")
      } else {
        pi <- matrix(pi, 1)
      }
    }
  }
  pi
}

check_Pi0                    <- function(Pi0) {
  if (!(length(Pi0) %in% c(1, 2))) {
    stop("Pi0 must be a numeric vector of length 1 or 2")
  }
  if (length(Pi0) == 1) {
    if (any(Pi0 < 0, Pi0 > 1)) {
      stop("If Pi0 is a numeric vector of length 1, it must belong to [0, 1]")
    }
  } else {
    if (any(Pi0 < 0, Pi0 > 1, Pi0[2] <= Pi0[1])) {
      stop("If Pi0 is a numeric vector of length 2, both elements must belong ",
           "to [0, 1], and the second element must be strictly larger than the",
           " first")
    }
  }
}

check_Pi1                    <- function(Pi1, delta) {
  if (!(length(Pi1) %in% c(1, 2))) {
    stop("Pi1 must be a numeric vector of length 1 or 2")
  }
  if (length(Pi1) == 1) {
    if (any(Pi1 < 0, Pi1 + delta > 1)) {
      stop("If Pi1 is a numeric vector of length 1, it must belong to ",
           "[0, 1 - delta]")
    }
  } else {
    if (any(Pi1 < 0, Pi1 > 1 - delta, Pi1[2] <= Pi1[1])) {
      stop("If Pi1 is a numeric vector of length 2, both elements must belong ",
           "to [0, 1 - delta], and the second element must be strictly larger ",
           "than the first")
    }
  }
}

check_real_pair_range_strict <- function(value1, value2, name1, name2, range) {
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

check_w                      <- function(w) {
  if (any(length(w) != 5, !is.numeric(w), w < 0, sum(w[-5]) == 0)) {
    stop("w must be a numeric vector of length 5, with all elements greater ",
         "than or equal to 0, and with at least one of the first 4 elements ",
         "strictly positive")
  }
}