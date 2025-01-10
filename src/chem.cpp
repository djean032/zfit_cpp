#include "chem.hpp"
int solve_rates(CVRhsFn rhs, double tin, double tout,
               std::vector<double> &initial_y,
               sunrealtype params[10]) {
  sundials::Context sunctx;

  // Initial Condition Vec
  N_Vector y = N_VNew_Serial(5, sunctx);
  sunrealtype *ydata = N_VGetArrayPointer(y);
  ydata[0] = initial_y[0];
  ydata[1] = initial_y[1];
  ydata[2] = initial_y[2];
  ydata[3] = initial_y[3];
  ydata[4] = initial_y[4];

  // Create matrix and solver
  SUNMatrix A = SUNDenseMatrix(5, 5, sunctx);

  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  // Create CVODE mem
  void *cvode_mem = CVodeCreate(CV_BDF, sunctx);

  // Initialize CVODE
  CVodeInit(cvode_mem, rhs, 0.0, y);

  // tolerances
  CVodeSStolerances(cvode_mem, 1e-6, 1e-10);

  // Attach matrix and solver
  CVodeSetLinearSolver(cvode_mem, LS, A);

  CVodeSetUserData(cvode_mem, params);

  CVode(cvode_mem, tout, y, &tin, CV_NORMAL);
 
  initial_y[0] = ydata[0];
  initial_y[1] = ydata[1];
  initial_y[2] = ydata[2];
  initial_y[3] = ydata[3];
  initial_y[4] = ydata[4];
  

  N_VDestroy(y);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  CVodeFree(&cvode_mem);

  return 0;
}

int solve_intensity(CVRhsFn rhs, double tin, double tout,
               std::vector<double> &initial_y,
               sunrealtype params[10]) {
  sundials::Context sunctx;
  N_Vector y = N_VNew_Serial(1, sunctx);
  sunrealtype *ydata = N_VGetArrayPointer(y);
  ydata[0] = initial_y[0];

  // Create matrix and solver
  SUNMatrix A = SUNDenseMatrix(1, 1, sunctx);
  SUNLinearSolver LS = SUNLinSol_Dense(y, A, sunctx);

  // Create CVODE mem
  void *cvode_mem = CVodeCreate(CV_BDF, sunctx);

  // Initialize CVODE
  CVodeInit(cvode_mem, rhs, 0.0, y);

  // tolerances
  CVodeSStolerances(cvode_mem, 1e-6, 1e-10);

  // Attach matrix and solver
  CVodeSetLinearSolver(cvode_mem, LS, A);

  CVodeSetUserData(cvode_mem, params);

  CVode(cvode_mem, tout, y, &tin, CV_NORMAL);

  initial_y[0] = ydata[0];

  N_VDestroy(y);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  CVodeFree(&cvode_mem);

  return 0;
}

int chem_func(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data) {
  sunrealtype *udata = static_cast<sunrealtype *>(user_data);
  const sunrealtype I = udata[0];
  const sunrealtype frq = udata[1];
  const sunrealtype SigmaG = udata[2];
  const sunrealtype SigmaS = udata[3];
  const sunrealtype SigmaT = udata[4];
  const sunrealtype t10 = udata[5];
  const sunrealtype t13 = udata[6];
  const sunrealtype t21 = udata[7];
  const sunrealtype t30 = udata[8];
  const sunrealtype t43 = udata[9];

  sunrealtype *ydata = N_VGetArrayPointer(y);
  sunrealtype *ydotdata = N_VGetArrayPointer(ydot);

  const sunrealtype nS0 = ydata[0];
  const sunrealtype nS1 = ydata[1];
  const sunrealtype nS2 = ydata[2];
  const sunrealtype nT1 = ydata[3];
  const sunrealtype nT2 = ydata[4];

  ydotdata[0] = -(SigmaG * I * nS0) / (h * frq) + nS1 / t10 + nT1 / t30;
  ydotdata[1] = (SigmaG * I * nS0) / (h * frq) - nS1 / t10 -
                (SigmaS * I * nS1) / (h * frq) + nS2 / t21 - nS1 / t13;
  ydotdata[2] = (SigmaS * I * nS1) / (h * frq) - nS2 / t21;
  ydotdata[3] =
      -(SigmaT * I * nT1) / (h * frq) + nT2 / t43 + nS1 / t13 - nT1 / t30;
  ydotdata[4] = (SigmaT * I * nT1) / (h * frq) - nT2 / t43;

  return 0;
}

int intensity_func(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data) {
  sunrealtype *udata = static_cast<sunrealtype *>(user_data);
  //const sunrealtype I = udata[0];
  const sunrealtype frq = udata[1];
  const sunrealtype SigmaG = udata[2];
  const sunrealtype SigmaS = udata[3];
  const sunrealtype SigmaT = udata[4];
  const sunrealtype t10 = udata[5];
  const sunrealtype t13 = udata[6];
  const sunrealtype t21 = udata[7];
  const sunrealtype t30 = udata[8];
  const sunrealtype t43 = udata[9];

  sunrealtype *ydata = N_VGetArrayPointer(y);
  sunrealtype *ydotdata = N_VGetArrayPointer(ydot);

  const sunrealtype I = ydata[0];

  y

  return 0;
};

