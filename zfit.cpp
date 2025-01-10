#include "zfit.hpp"
#include "src/chem.hpp"


int main(int argc, char *argv[]) {
  auto t1 = std::chrono::high_resolution_clock::now();

// Initialize needed values
  int t_slices = 10;
    
  int z_slices = 10;

  int z_pos_slices = 101;

  std::vector<double> z_pos(z_pos_slices);

  std::vector<double> t(t_slices);

  std::vector<double> z(z_slices);

  std::vector<std::vector<double>> intensities(z_pos_slices,
                                               std::vector<double>(t_slices));

  std::vector<std::vector<double>> final_intensities(
      z_pos_slices, std::vector<double>(t_slices));

  std::vector<std::vector<double>> populations(z_pos_slices,
                                               std::vector<double>(5));

  std::vector<double> total_intensities(z_pos_slices);

  std::vector<double> initial_intensities(z_pos_slices);

  sunrealtype tret = 0.0;
  sunrealtype tout = tret + 1e-8;

  sunrealtype zret = 0.0;
  sunrealtype zout = zret;

  std::vector<double> initial_conditions = {2.75e18, 0, 0, 0, 0};
  sunrealtype params[10] = {1.0e-9,  frq,     5.48e-18, 16.0e-18, 28.1e-18,
                            1.0e-12, 1.0e-12, 1.0e-12,  1.05e-7,  1.0e12};
  for (int i = 0; i < 10; i++) {
    solve_rates(chem_func, tret, tout, initial_conditions, params);
  }

  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout
        << "Populations: "
        << initial_conditions[0] << " "
        << initial_conditions[1] << " "
        << initial_conditions[2] << " "
        << initial_conditions[3] << " "
        << initial_conditions[4]
        << "\n";
  std::cout
      << "Time: "
      << std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() /
             1e6
      << " ms\n";
  return 0;
}
