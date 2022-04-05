#include <Rcpp.h>
#include <fstream>

// [[Rcpp::export]]
Rcpp::NumericVector single_linkage(std::string file) {
  Rcpp::NumericVector out = Rcpp::NumericVector();
  std::ifstream infile(file);
  int seq1, seq2;
  float dist;
  while(infile >> seq1 >> seq2 >> dist) {
    std::cout << "seq1: " << seq1 << "; seq2: " << seq2 << "; dist: " << dist << "\n";
    out.push_back(dist);
  }
  return out;
}
