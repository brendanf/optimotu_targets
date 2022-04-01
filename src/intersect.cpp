#include <Rcpp.h>

// [[Rcpp::export]]
int intersect_length_(const std::vector<int> &c, const std::vector<int> &k) {
  int l = 0;
  auto ci = c.begin();
  auto ki = k.begin();
  while (ci != c.end() && ki != k.end()) {
    if (*ci < *ki) {
      ++ci;
    } else if (*ki < *ci) {
      ++ki;
    } else {
      ++l;
      ++ci;
      ++ki;
    }
  }
  return l;
}
