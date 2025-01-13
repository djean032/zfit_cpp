#ifndef LASER_HPP
#define LASER_HPP

double waist(double z, double zr, double w0);

double power(double t, double ep, double tau, double wid);

double irradiance(double r, double t, double z, double zr, double w0, double ep,
                  double tau, double wid);
#endif // !LASER_HPP
