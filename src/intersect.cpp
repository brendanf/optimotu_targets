#include <Rcpp.h>

// [[Rcpp::export]]
int intersect_length_(const std::vector<int> &c, const std::vector<int> &k) {
  int l = 0;
  auto ci = c.begin();
  auto ki = k.begin();
  auto cend = c.end();
  auto kend = k.end();
  while (true) {
    if (*ci < *ki) {
      if (++ci == cend) return l;
    } else if (*ki < *ci) {
      if (++ki == kend) return l;
    } else {
      ++l;
      if (++ci == cend || ++ki == kend) return l;
    }
  }
  return l;
}
