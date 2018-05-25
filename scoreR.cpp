#include <Rcpp.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
DoubleVector update_mass(DoubleVector peak_masses, IntegerMatrix rule) {

  int id = floor(unif_rand() * rule.nrow());
  int beg = rule(id, 0);
  int end = rule(id, 1);

  double beg_mass = peak_masses[beg], end_mass = peak_masses[end];
  double delta = R::runif(-beg_mass, end_mass);

  DoubleVector res = clone(peak_masses);

  res[beg] = beg_mass + delta;
  res[end] = end_mass - delta;

  return res;
}


// [[Rcpp::export]]
double score_peak(DoubleVector spectrum, DoubleVector peak_masses, double nlp_mass, double MASS_PROTON, bool keep_both) {
  std::sort(peak_masses.begin(), peak_masses.end(), [](double a, double b) { return a < b; });

  double score = 0;
  double product_ion_thresh = 0.02;
  auto pmb = peak_masses.begin();
  auto pme = std::reverse_iterator<decltype(peak_masses.end())>(peak_masses.end());
  auto rpme = std::reverse_iterator<decltype(peak_masses.begin())>(peak_masses.begin());

  for (const auto& rp: spectrum) {
    double thr = rp - product_ion_thresh - MASS_PROTON;

    while (pmb != peak_masses.end() && *pmb <= thr) {
      pmb++;
      // std::cout << *pmb << std::endl;
    }

    if (pmb != peak_masses.end()) {
      if (std::abs(*pmb + MASS_PROTON - rp) < product_ion_thresh) {
        score += 1;
        // std::cout << *pmb << std::endl;
        continue;
      }
    }

    if (keep_both) {
      if (pmb == peak_masses.end())
        break;
      continue;
    }

    while (pme != rpme && nlp_mass - *pme <= thr) {
      pme++;
      // std::cout << nlp_mass - *pme + MASS_PROTON << std::endl;
    }

    if (pme != rpme) {
      if (std::abs(nlp_mass - *pme + MASS_PROTON - rp) < product_ion_thresh) {
        score += 1;
        // std::cout << nlp_mass - *pme + MASS_PROTON << std::endl;
        continue;
      }
    }

    if (pmb == peak_masses.end() && pme == rpme) {
      break;
    }
  }

  return score;
}
