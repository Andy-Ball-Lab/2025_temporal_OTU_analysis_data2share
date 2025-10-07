# --- Patch for plsmod mixOmics bug ("incorrect number of dimensions") ---
fix_multi_numeric_preds <- function(object, new_data, comps = NULL) {
  tmp_pred <- predict(object$fit, new_data, dist = "mahalanobis.dist")
  tmp_pred <- tmp_pred$predict
  n <- dim(tmp_pred)[1]
  p <- dim(tmp_pred)[2]
  q <- dim(tmp_pred)[3]

  if (is.null(comps)) comps <- q
  comps <- comps[comps <= q]
  tmp_grid <- tibble::tibble(num_comp = comps)
  tmp_pred <- tmp_pred[, , comps, drop = FALSE]

  if (p > 1) {
    nms <- dimnames(tmp_pred)[[2]]
    new_nms <- paste0(".pred_", nms)
    tmp_pred <- purrr::map(1:n, function(.x)
      t(as.matrix(tmp_pred[.x, , , drop = FALSE])))
  } else {
    new_nms <- ".pred"
    tmp_pred <- tmp_pred[, 1, , drop = FALSE]  # keep 3D structure
    tmp_pred <- purrr::map(
      1:n,
      function(.x) data.frame(.pred = as.numeric(tmp_pred[.x, , , drop = FALSE]))
    )
  }

  tmp_pred <- purrr::map(
    tmp_pred,
    function(.x)
      dplyr::bind_cols(tmp_grid, stats::setNames(tibble::as_tibble(.x), new_nms))
  )

  tibble::tibble(.pred = tmp_pred)
}
