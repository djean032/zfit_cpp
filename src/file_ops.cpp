#include <fstream>
#include <sstream>

struct experimental_data {
  std::vector<double> position;
  std::vector<double> ai0;
  std::vector<double> ai1;
  std::vector<double> ai2;
};

experimental_data read_data(const std::string &filename) {
  experimental_data data;
  std::ifstream file(filename);
  std::string line;
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file");
  }
  while (!file.eof()) {
    std::getline(file, line);
    if (line.empty()) {
      continue;
    }
    else if (std::strcmp(line.substr(0, 1).c_str(), "z") == 0) {
      std::istringstream iss(line);
      double z, a0, a1, a2;
      iss >> z >> a0 >> a1 >> a2;
      data.position.push_back(z);
      data.ai0.push_back(a0);
      data.ai1.push_back(a1);
      data.ai2.push_back(a2);
    }
    else {
      continue;
    }
  }
  return data;
};
