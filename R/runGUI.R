#' A function to launch Repliscope in interactive mode (Shiny app).
#' @keywords shiny genomics bioinformatics replication
#' @import ggplot2 shiny
#' @importFrom colourpicker colourInput colourPicker
#' @export

runGUI <- function() {
  appDir <- system.file("shinyApp", package = "Repliscope")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing Repliscope.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
