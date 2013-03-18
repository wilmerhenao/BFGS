

bool themin(double * x, double * gradpoints, size_t i, const size_t n){
  // This function calculates the minimum coordinate distance from x to gradpoints
  // in a certain gradpoints segment
  double valmin = 1.0;
  for (size_t j; j < n; j++){
    valmin MIN(valmin, std::abs(x[j] - gradpoints[j + (i * n)]));
  }
  if (valmin < 1e-14)
    return false
  else
    return true;
}
