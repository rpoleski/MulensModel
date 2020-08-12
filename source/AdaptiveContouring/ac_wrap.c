#define PY_SSIZE_T_CLEAN
#include <Python.h>

#ifndef FLOAT
#define FLOAT double
#endif


FLOAT ld_linear(int n, FLOAT gam[], FLOAT rho);

FLOAT mag_binext(FLOAT y1, FLOAT y2, FLOAT rho, FLOAT d, FLOAT q,
        FLOAT (*ld_func)(int,FLOAT*,FLOAT), int n, FLOAT gam[],
        FLOAT acc, FLOAT ld_acc);

FLOAT mag_binpt(FLOAT y1, FLOAT y2, FLOAT d, FLOAT q);


static PyObject * Adaptive_Contouring_Linear_wrapper(PyObject *self, PyObject *args) {
  double d, q, y1, y2, rho, gamma, acc, ld_acc;
  double mag, gam[1];

  if (!PyArg_ParseTuple(args, "dddddddd", &d, &q, &y1, &y2, &rho, &gamma, &acc, &ld_acc)) return NULL;

  if (gamma == 0.0) {
    mag = mag_binext(y1, y2, rho, d, q, NULL, -1, NULL, acc, ld_acc);
  }
  else {
    gam[0] = gamma;
    mag = mag_binext(y1, y2, rho, d, q, ld_linear, 1, gam, acc, ld_acc);
  }

  return Py_BuildValue("d", mag);
}

static PyMethodDef ACMethods[] = {
    {"Adaptive_Contouring_Linear", Adaptive_Contouring_Linear_wrapper, METH_VARARGS, "some notes here"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef ACmodule = {
    PyModuleDef_HEAD_INIT,
    "AdaptiveContouring",
    NULL,
    -1,
    ACMethods
};

PyMODINIT_FUNC PyInit_AdaptiveContouring(void) {
  return PyModule_Create(&ACmodule);
}
