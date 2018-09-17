.onAttach <- function(...) { # nocov start
    if(interactive()){
        msg <- c("Thank you for using clustRviz!",
                 "The current logging level is",
                 sQuote(paste0(clustRviz_logger_level(), ".")),
                 "To change this, see ?clustRviz_logging.")

        packageStartupMessage(paste(msg, collapse=" "))
    }
}                            # nocov end
