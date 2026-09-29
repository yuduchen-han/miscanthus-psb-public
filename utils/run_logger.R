# ==============================================================================
# run_logger.R
# 简单科研分析运行日志
# ==============================================================================

start_log <- function(script_name, log_dir = "logs") {
  dir.create(
    log_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  log_file <- file.path(
    log_dir,
    paste0(
      tools::file_path_sans_ext(basename(script_name)),
      "_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      ".log"
    )
  )

  log_con <- file(log_file, open = "wt")

  sink(log_con, split = TRUE)
  sink(log_con, type = "message")

  cat("============================================================\n")
  cat("SCRIPT:", script_name, "\n")
  cat("START :", format(Sys.time()), "\n")
  cat("============================================================\n\n")

  invisible(log_con)
}


end_log <- function(log_con) {
  cat("\n============================================================\n")
  cat("FINISH:", format(Sys.time()), "\n")
  cat("============================================================\n")

  cat("\nSESSION INFO\n")
  cat("------------------------------------------------------------\n")
  print(sessionInfo())

  sink(type = "message")
  sink()

  close(log_con)
}
