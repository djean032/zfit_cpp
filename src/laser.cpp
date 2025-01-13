#include "laser.hpp"
#include <cmath>
#include <numbers>

double waist(double z, double zr, double w0) { 
  double waist = w0 * sqrt(1 + std::pow((z/zr), 2));
  return waist; 
}

double power(double t, double ep, double tau, double wid) {
  double power = (0.94 * ep / tau) * std::exp(std::pow(-t / wid, 2));
  return power;
}

double irradiance(double r, double t, double z, double zr, double w0, double ep,
                  double tau, double wid) { 
  double w_z = waist(z, zr, w0);
  double phi = power(t, ep, tau, wid);

  double irr = 0.0001 * (2 * phi) / std::numbers::pi * std::pow(w_z, 2) *
               std::exp(-2 * std::pow(r / w_z, 2));

  return irr; 
}
