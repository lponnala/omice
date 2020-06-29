
# Use this .Rprofile on Windows, by copying it to *one* of the following locations:
# - R.home()/etc/Rprofile.site
# - Sys.getenv("R_USER")/.Rprofile

cat("Getting ready ... ", "\n", sep = "")

# set default CRAN mirror
options(repos=structure(c(CRAN="http://cloud.r-project.org/")))

# set library path
.libPaths("C:/Users/lalit_ponnala/Documents/R")

# do not ask "save workspace image?" when quitting via q()
utils::assignInNamespace(
  "q", function(save = "no", status = 0, runLast = TRUE) {
    .Internal(quit(save, status, runLast))
  }, "base")

# do not convert strings to factors
options("stringsAsFactors" = FALSE)

# set timezone
Sys.setenv(TZ = "America/New_York") # local time
# Sys.setenv(TZ = "UTC") # UTC makes xts work smoother

# !!BEGIN!! start commenting here if you are updating/installing packages (be sure to re-start R)
options(
  radian.color_scheme = "native",
  radian.editing_mode = "vi",
  radian.indent_lines = FALSE,
  radian.auto_match = FALSE, # FALSE enables you to send lines to terminal using VScode
  radian.auto_indentation = TRUE,
  radian.tab_size = 4,
  radian.complete_while_typing = TRUE,
  radian.completion_timeout = 0.05,
  radian.auto_width = TRUE,
  radian.insert_new_line = FALSE
)

# httr::set_config(httr::use_proxy(url = "http://web.url.com", port = 8080, username = "lponnala", password = rawToChar(base64enc::base64decode(Sys.getenv("USERPASSWORD")))))
options("pillar.subtle" = FALSE)
options("pillar.neg" = FALSE)
# options("shiny.launch.browser" = TRUE)
# alphavantager::av_api_key(Sys.getenv("ALPHAVANTAGE_API_KEY"))
# quantmod::setDefaults("getSymbols.av", api.key = Sys.getenv("ALPHAVANTAGE_API_KEY"))
# tidyquant::quandl_api_key(Sys.getenv("QUANDL_API_KEY"))

# -- reasons for loading various packages --
# dplyr : for many reasons (avoid printing large datasets, access to tibble, etc)
# ggplot2 : for plots (otherwise typing ggplot2:: for every aes, geom, etc becomes cumbersome)
# shiny : for apps (otherwise typing shiny:: for every tabPanel, fluidRow, etc becomes cumbersome)
# magrittr : for additional pipe operators
# lpu : miscellaneous functions
invisible(
  # %>% and %<>% are exported from lpu, so no need to load magrittr
  sapply(intersect(c("ggplot2","shiny","magrittr"),rownames(utils::installed.packages())), function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, logical.return = TRUE))
  })
)
# !!END!! stop commenting here if you are updating/installing packages (be sure to re-start R)

cat("done!", "\n", "Attached packages: ", paste(sort(.packages()), collapse = ","), "\n", sep = "")
