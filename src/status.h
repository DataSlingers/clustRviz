// This file is a heavily adapted and simplified version of RProgress.h from
// the `progress` package
//
// The original is located at
// https://github.com/r-lib/progress/blob/master/inst/include/RProgress.h
//
// Original authors are Gábor Csárdi & Rich FitzJohn
// Original copyright holders are Gábor Csárdi & RStudio Inc.

#ifndef CLUSTRVIZ_STATUS_H
#define CLUSTRVIZ_STATUS_H 1

#include <unistd.h>
#include <sys/time.h>

#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <ctime>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

#include "clustRviz_base.h"

class StatusPrinter {

public:
  StatusPrinter(bool show_progress,
                uint total_fusions,
                std::string format = "",
                int width = Rf_GetOptionWidth() - 2) :

  show_progress(show_progress),
  count(0),
  output_supported(is_supported()),
  format(format),
  width(width),
  last_tick(std::clock()),
  total_fusions(total_fusions),
  fusions(0),
  v_norm_init(0),
  v_norm(0),
  iter(0),
  gamma(0){}

  ~StatusPrinter() {
    this->update_format();
  }

  void set_width(int width){
    this->width = width;
  }

  void set_v_norm_init(double v_norm_init){
    this->v_norm_init = v_norm_init;
  }

  void update(int fusions_, double v_norm_, uint iter_, double gamma_) {
    double time_since_last_tick = (std::clock() - last_tick) / ((double) CLOCKS_PER_SEC);

    if(time_since_last_tick < CLUSTRVIZ_STATUS_UPDATE_TIME_SECS){
      return; // Don't need to print too frequently
    }

    force_update(fusions_, v_norm_, iter_, gamma_);
  }

  void force_update(int fusions_, double v_norm_, uint iter_, double gamma_){
    Rcpp::checkUserInterrupt(); // Check essentially every time we do the progress bar

    if(!show_progress){
      return; // Not printing, so return after checking for CTRL-C
    }

    // Update internal counters
    fusions = fusions_;
    v_norm  = v_norm_;
    iter    = iter_;
    gamma   = gamma_;

    if( (count % CLUSTRVIZ_STATUS_WIDTH_CHECK) == 0){
      set_width(Rf_GetOptionWidth() - 2);
      this->update_format();
    }

    count++;

    if (fusions == total_fusions){
      this->terminate();
    } else {
      this->render();
    }
  }

private:
  const bool show_progress;    // Do we print output?
  uint count;                  // Number of times we have updated the printer
  const bool output_supported; // Do we support the desired print location?
  std::string format;          // Bar template -- hard coded above
  int width;                   // Width of progress bar (used to set template and control print length)
  std::string last_draw;       // Last progress bar drawn
  std::clock_t last_tick;      // Time of last PB tick

  // ClustRviz-specific measures of progress
  const uint total_fusions; // Number of edges that need to be fused before termination
  uint fusions;             // Number of edges fused so far
  double v_norm_init;       // Initial squared Frobenius norm of V
  double v_norm;           // Current squared Frobenius norm of V
  uint iter;               // What iteration we are on
  double gamma;            // Current regularization level

  void update_format(){
    if(this->width >= 120){
      this->format = ":spin Iteration :iter. Current Gamma :gamma. Fused :edges / :total edges. Percent fusion achieved :pct%";
    } else if(this->width >= 80){
      this->format = ":spin Iteration :iter. Current Gamma :gamma. Fused :edges / :total edges";
    } else if(this->width >= 40){
      this->format = ":spin Current Iteration :iter. Current Gamma :gamma";
    } else if(this->width >= 20){
      this->format = ":spin Current Iteration :iter";
    } else if (this->width >= 5){
      this->format = ":spin";
    } else {
      this->format = "";
    }
  }

  void render() {
    // Update the progress output
    if (!output_supported){
      // If we are printing to some weird output location that we don't support, skip
      return;
    }

    std::string str = format;
    std::stringstream buffer;

    // Iteration
    buffer << iter;
    replace_all(str, ":iter", buffer.str());
    buffer.str(""); buffer.clear();

    // Current edges
    buffer << fusions;
    replace_all(str, ":edges", buffer.str());
    buffer.str(""); buffer.clear();

    // Max edges
    buffer << total_fusions;
    replace_all(str, ":total", buffer.str());
    buffer.str(""); buffer.clear();

    // Percent fusion achieved
    double pct_fusion = std::fmax(0, 1.0 - v_norm / v_norm_init) * 100;
    buffer << std::setw(5) << pct_fusion;
    replace_all(str, ":pct", buffer.str());
    buffer.str(""); buffer.clear();

    // Current gamma
    buffer << gamma;
    replace_all(str, ":gamma", buffer.str());
    buffer.str(""); buffer.clear();

    // spin
    replace_all(str, ":spin", spin_symbol());

    if (last_draw != str) {
      if (last_draw.length() > str.length()) { clear_line(width); }
      cursor_to_start();
      Rprintf(str.c_str());
      last_draw = str;
    }
  }

  void terminate() {
    if (!output_supported){
      // Never printed -- ok to return w/o doing anything
      return;
    }

    // Clean up and exit
    clear_line(width);
    cursor_to_start();

    // RStudio appears to not implement \r quite right,
    // so let's add a new line to be safe
    if(is_r_studio()){
      Rprintf("\n");
    }
  }

  std::string spin_symbol() {
    const char symbols[4] = {'-', '\\', '|', '/'};
    return std::string(1, symbols[(count - 1) % 4]);
  }

  void clear_line(int width) {

    char *spaces = (char*) calloc(width + 2, sizeof(char));
    if (!spaces) Rf_error("Progress bar: out of memory");
    for (int i = 1; i <= width; i++) spaces[i] = ' ';
    spaces[0] = '\r';
    spaces[width + 1] = '\0';

    Rprintf(spaces);
    free(spaces);
  }

  void cursor_to_start() {
    Rprintf("\r");
  }

  bool is_r_app() const {
    char *v = std::getenv("R_GUI_APP_VERSION");
    return v != 0;
  }

  bool is_r_studio() const {
    char *v = std::getenv("RSTUDIO");
    return v != 0 && v[0] == '1' && v[1] == '\0';
  }

  // If stdout is a terminal, or R Studio or macOS R.app
  // On windows, stdout is a terminal, apparently
  bool is_supported() {
    return (isatty(1) || is_r_studio() || is_r_app());
  }

  static void replace_all(std::string& str, const std::string& from,
                          const std::string& to) {
    if (from.empty()) return;

    size_t start_pos = 0;

    while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
      str.replace(start_pos, from.length(), to);
      start_pos += to.length();
    }
  }

};

#endif
