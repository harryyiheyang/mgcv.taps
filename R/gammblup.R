#' Selected Individual BLUPs from a gammfast Fit
#'
#' Returns random-trajectory BLUPs and row-level random-effect predictions for
#' explicitly requested individuals. `ids` is mandatory: the function never
#' constructs predictions for every fitted individual by default. Stored BLUPs
#' are used for fitted IDs. For new IDs in `newdata`, `new.level` explicitly
#' chooses an error, a zero population prediction, or conditional BLUP
#' estimation from the supplied outcomes while all fitted mean and covariance
#' parameters remain fixed.
#'
#' Duplicate values in `ids` are retained as separate requests in input order.
#' Within each request, rows retain their original fitted-data or `newdata`
#' order. A requested ID absent from the selected data is an error. New-data
#' custom-fs covariates outside their fitted ranges are clamped to the fitted
#' boundaries and marked in the returned row table.
#'
#' @param object A fitted `gammfast` object.
#' @param ids Mandatory vector of individual IDs to return, in requested order.
#' @param newdata Optional data frame. If omitted, fitted rows and stored BLUPs
#'   are used. If supplied, it must contain the fitted subject ID and all
#'   random-effect covariates.
#' @param new.level Handling of IDs not present in the fitted object:
#'   `"error"`, `"zero"`, or `"estimate"`. Estimation additionally requires
#'   the response and all global mean predictors in `newdata`.
#' @param weights Optional prior weights for `newdata`; used only when
#'   `new.level = "estimate"`.
#' @param max.iter Maximum frozen-parameter PIRLS/BLUP iterations for new IDs.
#' @param tol Relative convergence tolerance for new-ID BLUP estimation.
#' @param n_threads Number of threads for the subject-block BLUP kernel.
#'
#' @return A `gammblup` list. `subjects` has one row per ID request and its
#'   BLUP coefficients. `rows` maps each request to original row numbers and
#'   random-effect predictions on the link scale. `blup` is the coefficient
#'   matrix in request order.
#' @export
gammblup <- function(object, ids, newdata = NULL,
                     new.level = c("error", "zero", "estimate"),
                     weights = NULL, max.iter = 50L, tol = 1e-8,
                     n_threads = 1L) {
  if (!inherits(object, "gammfast")) {
    stop("object must be a 'gammfast' fit.")
  }
  if (missing(ids)) {
    stop("ids is required; gammblup never predicts all individuals by default.")
  }
  if (!is.atomic(ids) || is.matrix(ids) || !length(ids)) {
    stop("ids must be a non-empty atomic vector.")
  }
  requested <- as.character(ids)
  if (anyNA(requested) || any(!nzchar(requested))) {
    stop("ids cannot contain missing or empty values.")
  }
  new.level <- match.arg(new.level)
  if (length(max.iter) != 1L || !is.finite(max.iter) || max.iter < 1 ||
      max.iter != as.integer(max.iter)) {
    stop("max.iter must be a positive integer.")
  }
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be one positive finite value.")
  }
  if (length(n_threads) != 1L || !is.finite(n_threads) || n_threads < 1 ||
      n_threads != as.integer(n_threads)) {
    stop("n_threads must be a positive integer.")
  }
  max.iter <- as.integer(max.iter)
  n_threads <- as.integer(n_threads)

  random <- gammfast_random_info(object)
  fitted_levels <- as.character(random$id.levels)
  fitted_index <- match(requested, fitted_levels)
  k <- ncol(random$B)
  coefficient_names <- colnames(object$random.effects)
  if (is.null(coefficient_names)) {
    coefficient_names <- paste0("cos", seq_len(k) - 1L)
  }

  if (is.null(newdata)) {
    if (!is.null(weights)) {
      stop("weights can be supplied only with newdata.")
    }
    absent <- unique(requested[is.na(fitted_index)])
    if (length(absent)) {
      stop(
        "Requested fitted-data IDs are unknown: ",
        paste(absent, collapse = ", "), "."
      )
    }
    blup <- as.matrix(object$random.effects[fitted_index, , drop = FALSE])
    rownames(blup) <- paste0("request", seq_along(requested))
    rows <- vector("list", length(requested))
    n_rows <- integer(length(requested))
    for (i in seq_along(requested)) {
      ii <- which(random$id.index == fitted_index[i])
      n_rows[i] <- length(ii)
      random_link <- rowSums(
        random$B[ii, , drop = FALSE] *
          matrix(blup[i, ], nrow = length(ii), ncol = k, byrow = TRUE)
      )
      rows[[i]] <- data.frame(
        request.index = i,
        individual = requested[i],
        data.row = ii,
        row.within.individual = seq_along(ii),
        known = TRUE,
        source = "stored",
        input.clamped = FALSE,
        time.clamped = FALSE,
        random.link = random_link,
        stringsAsFactors = FALSE
      )
    }
    subjects <- data.frame(
      request.index = seq_along(requested),
      individual = requested,
      known = TRUE,
      source = "stored",
      n.rows = n_rows,
      iterations = 0L,
      stringsAsFactors = FALSE
    )
    subjects <- cbind(subjects, as.data.frame(blup, check.names = FALSE))
    ans <- list(
      subjects = subjects,
      rows = do.call(rbind, rows),
      blup = blup,
      requested.ids = requested,
      data.source = "fitted",
      new.level = NA_character_,
      groups = random$groups,
      full.random.design = FALSE,
      call = match.call()
    )
    class(ans) <- "gammblup"
    return(ans)
  }

  if (!is.data.frame(newdata)) stop("newdata must be a data frame.")
  id_name <- random$id
  if (!id_name %in% names(newdata)) {
    stop(
      "newdata is missing the required subject ID variable: ", id_name, "."
    )
  }
  row_id <- as.character(newdata[[id_name]])
  selected <- !is.na(row_id) & row_id %in% unique(requested)
  absent <- requested[!requested %in% row_id[selected]]
  if (length(absent)) {
    stop(
      "Requested IDs have no rows in newdata: ",
      paste(unique(absent), collapse = ", "), "."
    )
  }
  selected_rows <- which(selected)
  if (!is.null(weights)) {
    if (!is.numeric(weights) || length(weights) != nrow(newdata) ||
        any(!is.finite(weights)) || any(weights < 0)) {
      stop("weights must be finite, non-negative, and have one value per newdata row.")
    }
  } else {
    weights <- rep(1, nrow(newdata))
  }

  selected_design <- gammfast_random_design(
    random, newdata[selected_rows, , drop = FALSE]
  )
  B_selected <- selected_design$B
  input_clamped <- selected_design$clamped
  selected_id <- row_id[selected_rows]
  unique_requested <- unique(requested)
  unique_fitted_index <- match(unique_requested, fitted_levels)
  unknown_ids <- unique_requested[is.na(unique_fitted_index)]
  if (length(unknown_ids) && new.level == "error") {
    stop(
      "newdata contains unknown fitted IDs: ",
      paste(unknown_ids, collapse = ", "),
      ". Set new.level to 'zero' or 'estimate' explicitly."
    )
  }

  unique_blup <- matrix(
    0, nrow = length(unique_requested), ncol = k,
    dimnames = list(unique_requested, coefficient_names)
  )
  known_unique <- !is.na(unique_fitted_index)
  unique_blup[known_unique, ] <- object$random.effects[
    unique_fitted_index[known_unique], , drop = FALSE
  ]
  iteration_count <- stats::setNames(
    integer(length(unique_requested)), unique_requested
  )
  if (length(unknown_ids) && new.level == "estimate") {
    estimate_rows_local <- which(selected_id %in% unknown_ids)
    estimate_rows <- selected_rows[estimate_rows_local]
    required <- unique(c(
      names(object$global$var.summary),
      all.vars(object$global.formula[[2L]])
    ))
    missing_required <- setdiff(required, names(newdata))
    if (length(missing_required)) {
      stop(
        "newdata is missing variables required to estimate new-ID BLUPs: ",
        paste(missing_required, collapse = ", "), "."
      )
    }
    if (any(vapply(newdata[estimate_rows, required, drop = FALSE],
                   anyNA, logical(1)))) {
      stop("Variables used to estimate new-ID BLUPs cannot contain missing values.")
    }

    g <- object$global
    g$coefficients <- object$coefficients
    g$sp <- object$sp
    nd <- newdata[estimate_rows, , drop = FALSE]
    eta_global <- as.numeric(stats::predict(g, newdata = nd, type = "link"))
    if (length(eta_global) != nrow(nd) || any(!is.finite(eta_global))) {
      stop("The global mean model produced invalid newdata predictions.")
    }
    response <- eval(
      object$global.formula[[2L]], envir = nd,
      enclos = environment(object$global.formula)
    )
    if (is.factor(response)) response <- as.integer(response)
    response <- as.numeric(response)
    if (length(response) != nrow(nd) || any(!is.finite(response))) {
      stop("The newdata response must evaluate to one finite value per selected row.")
    }
    estimate_id <- match(selected_id[estimate_rows_local], unknown_ids)
    B_estimate <- B_selected[estimate_rows_local, , drop = FALSE]
    u <- matrix(0, nrow = length(unknown_ids), ncol = k)
    converged <- FALSE
    for (iteration in seq_len(max.iter)) {
      eta <- eta_global + rowSums(B_estimate * u[estimate_id, , drop = FALSE])
      work <- gammfast_working(
        object$family, response, eta, weights[estimate_rows],
        nthreads = n_threads
      )
      sw <- sqrt(work$w / object$sigma2)
      moments <- gammfast_gaussian_moments(
        residual = (work$z - eta_global) * sw,
        B = B_estimate * sw,
        id = estimate_id,
        G = object$G,
        sigma2 = 1,
        n_threads = n_threads
      )
      u_new <- moments$u
      change <- max(abs(u_new - u)) / (1 + max(abs(u)))
      u <- u_new
      if (change < tol) {
        converged <- TRUE
        break
      }
    }
    if (!converged) {
      stop("New-ID BLUP estimation reached max.iter before convergence.")
    }
    unique_blup[unknown_ids, ] <- u
    iteration_count[unknown_ids] <- iteration
  }

  request_map <- match(requested, unique_requested)
  blup <- unique_blup[request_map, , drop = FALSE]
  rownames(blup) <- paste0("request", seq_along(requested))
  request_known <- !is.na(fitted_index)
  request_source <- ifelse(
    request_known, "stored",
    if (new.level == "estimate") "estimated" else "zero"
  )
  rows <- vector("list", length(requested))
  n_rows <- integer(length(requested))
  for (i in seq_along(requested)) {
    ii_local <- which(selected_id == requested[i])
    n_rows[i] <- length(ii_local)
    random_link <- rowSums(
      B_selected[ii_local, , drop = FALSE] *
        matrix(blup[i, ], nrow = length(ii_local), ncol = k, byrow = TRUE)
    )
    rows[[i]] <- data.frame(
      request.index = i,
      individual = requested[i],
      data.row = selected_rows[ii_local],
      row.within.individual = seq_along(ii_local),
      known = request_known[i],
      source = request_source[i],
      input.clamped = input_clamped[ii_local],
      time.clamped = input_clamped[ii_local],
      random.link = random_link,
      stringsAsFactors = FALSE
    )
  }
  subjects <- data.frame(
    request.index = seq_along(requested),
    individual = requested,
    known = request_known,
    source = request_source,
    n.rows = n_rows,
    iterations = unname(iteration_count[requested]),
    stringsAsFactors = FALSE
  )
  subjects <- cbind(subjects, as.data.frame(blup, check.names = FALSE))
  ans <- list(
    subjects = subjects,
    rows = do.call(rbind, rows),
    blup = blup,
    requested.ids = requested,
    data.source = "newdata",
    new.level = new.level,
    groups = random$groups,
    full.random.design = FALSE,
    call = match.call()
  )
  class(ans) <- "gammblup"
  ans
}
