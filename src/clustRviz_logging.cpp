#include "clustRviz.h"

// [[Rcpp::export]]
void clustRviz_set_logger_level_cpp(int level){
    auto logger_level = static_cast<ClustRVizLoggerLevel>(level);
    ClustRVizLogger::set_level(logger_level);
}

// [[Rcpp::export]]
int clustRviz_get_logger_level_cpp(){
    auto logger_level = static_cast<int>(ClustRVizLogger::get_level());
    return logger_level;
}

// [[Rcpp::export]]
void clustRviz_log_cpp(int level, Rcpp::StringVector x){
    auto msg_level = static_cast<ClustRVizLoggerLevel>(level);
    std::string msg = Rcpp::as<std::string>(x[0]);
    if(msg_level >= ClustRVizLoggerLevel::ERRORS){
        ClustRVizLogger::error(msg);
    } else if(msg_level >= ClustRVizLoggerLevel::WARNING){
        ClustRVizLogger::warning(msg);
    } else if(msg_level >= ClustRVizLoggerLevel::MESSAGES){
        ClustRVizLogger::message(msg);
    } else if(msg_level >= ClustRVizLoggerLevel::INFO){
        ClustRVizLogger::info(msg);
    } else if(msg_level >= ClustRVizLoggerLevel::DEBUG){
        ClustRVizLogger::debug(msg);
    }
}
