#ifndef UTILS_HPP
#define UTILS_HPP

#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

#include <algorithm>
#include <charconv>
#include <chrono>
#include <iostream>
#include <string>
#include <sundials/sundials_core.hpp>
#include <vector>

static const double h = 6.62607015e-34;
static const double c = 299792458;
static const double frq = c / 532e-9;

static int check_flags(const int flag, const std::string function) {
  if (flag < 0) {
    std::cerr << "Error: " << function << " returned error code " << flag
              << "\n";
    return 1;
  }
  return 0;
}

static int check_func(const void *ptr, const std::string function) {
  if (ptr == nullptr) {
    std::cerr << "Error: " << function << " function doesn't exit\n";
    return 1;
  }
  return 0;
}

template <typename T>
inline void find_arg(std::vector<std::string> &args, const std::string key,
                     T &dest) {
  auto it = std::find(args.begin(), args.end(), key);
  std::string str = *(it + 1);
  if (it != args.end()) {
    std::from_chars(str.data(), str.data() + str.size(), dest);
    args.erase(it, it + 2);
  }
}
#endif // !UTILS_HPP
