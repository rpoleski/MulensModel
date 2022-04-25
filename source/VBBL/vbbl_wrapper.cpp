#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "VBBinaryLensingLibrary.h"

/*
Based on code written by Przemek Mroz
*/

static PyObject * 
VBBinaryLensing_BinaryMagDark_wrapper(PyObject *self, PyObject *args) {
  double a, q, y1, y2, RSv, a1, Tol, mag;
  static VBBinaryLensing VBBL;

  if (!PyArg_ParseTuple(args, "ddddddd", &a, &q, &y1, &y2, &RSv, &a1, &Tol)) return NULL;

  mag = VBBL.BinaryMagDark(a, q, y1, y2, RSv, a1, Tol);

  return Py_BuildValue("d", mag);
}

static PyObject * 
VBBinaryLensing_BinaryMag0_wrapper(PyObject *self, PyObject *args) {
  double a, q, y1, y2, K, G, Gi, mag;
  static VBBinaryLensing VBBL;

  if (!PyArg_ParseTuple(args, "ddddddd", &a, &q, &y1, &y2, &K, &G, &Gi)) return NULL;

  mag = VBBL.BinaryMag0_shear(a, q, y1, y2, K, G, Gi);

  return Py_BuildValue("d", mag);
}

PyObject * makelist(double *array, size_t size) {
    PyObject *l = PyList_New(size);
    for (size_t i = 0; i != size; ++i) {
        PyList_SET_ITEM(l,i,PyFloat_FromDouble(array[i]));
    }
    return l;
}

static PyObject *
VBBL_SG12_5_wrapper(PyObject *self, PyObject *args) {
  static VBBinaryLensing VBBL;
  complex complex_poly[6], complex_roots[5];
  double roots[10];
  double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11;
  int i;

  if (!PyArg_ParseTuple(args, "dddddddddddd", &p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8, &p9, &p10, &p11)) return NULL;

    complex_poly[0] = complex(p0, p6);
    complex_poly[1] = complex(p1, p7);
    complex_poly[2] = complex(p2, p8);
    complex_poly[3] = complex(p3, p9);
    complex_poly[4] = complex(p4, p10);
    complex_poly[5] = complex(p5, p11);

    VBBL.cmplx_roots_gen(complex_roots, complex_poly, 5, true, true);

    for (i=0; i<5; i++) {
        roots[i] = complex_roots[i].re;
        roots[i+5] = complex_roots[i].im;
    }

  return makelist(roots, 10);
}

static PyObject * 
VBBinaryLensing_BinaryMag_wrapper(PyObject *self, PyObject *args) {
  double a, q, y1, y2, mag;
  static VBBinaryLensing VBBL;

  if (!PyArg_ParseTuple(args, "dddd", &a, &q, &y1, &y2)) return NULL;

  mag = VBBL.BinaryMag0(a, q, y1, y2);

  return Py_BuildValue("d", mag);
}

static PyObject *
VBBL_SG12_9_wrapper(PyObject *self, PyObject *args) {
  static VBBinaryLensing VBBL;
  complex complex_poly[10], complex_roots[9];
  double roots[18];
  double p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19;
  int i;

  if (!PyArg_ParseTuple(args, "dddddddddddddddddddd", &p0, &p1, &p2, &p3, &p4, &p5, &p6, &p7, &p8, &p9, &p10, &p11,
                                            &p12, &p13, &p14, &p15, &p16, &p17, &p18, &p19)) return NULL;
    complex_poly[0] = complex(p0, p10);
    complex_poly[1] = complex(p1, p11);
    complex_poly[2] = complex(p2, p12);
    complex_poly[3] = complex(p3, p13);
    complex_poly[4] = complex(p4, p14);
    complex_poly[5] = complex(p5, p15);
    complex_poly[6] = complex(p6, p16);
    complex_poly[7] = complex(p7, p17);
    complex_poly[8] = complex(p8, p18);
    complex_poly[9] = complex(p9, p19);

    VBBL.cmplx_roots_gen(complex_roots, complex_poly, 9, true, true);

    for (i=0; i<9; i++) {
        roots[i] = complex_roots[i].re;
        roots[i+9] = complex_roots[i].im;
    }

  return makelist(roots, 18);
}

static PyMethodDef VBBLMethods[] = {
    {"VBBinaryLensing_BinaryMagDark", VBBinaryLensing_BinaryMagDark_wrapper, METH_VARARGS, "some notes here"},
    {"VBBinaryLensing_BinaryMag0", VBBinaryLensing_BinaryMag0_wrapper, METH_VARARGS, "some notes here"},
    {"VBBL_SG12_5", VBBL_SG12_5_wrapper, METH_VARARGS, "some notes here"},
    {"VBBL_BinaryMag", VBBinaryLensing_BinaryMag_wrapper, METH_VARARGS, "some notes here"},
    {"VBBL_SG12_9", VBBL_SG12_9_wrapper, METH_VARARGS, "some notes here"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef VBBLmodule = {
    PyModuleDef_HEAD_INIT,
    "VBBL",
    NULL,
    -1,
    VBBLMethods
};

extern "C" {
PyMODINIT_FUNC PyInit_VBBL(void) {
  return PyModule_Create(&VBBLmodule);
}
}

