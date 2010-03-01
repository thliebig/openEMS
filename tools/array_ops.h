#ifndef ARRAY_OPS_H
#define ARRAY_OPS_H

#include "../FDTD/operator.h"

FDTD_FLOAT*** Create3DArray(unsigned int* numLines);
void Delete3DArray(FDTD_FLOAT*** array, unsigned int* numLines);

FDTD_FLOAT**** Create_N_3DArray(unsigned int* numLines);
void Delete_N_3DArray(FDTD_FLOAT**** array, unsigned int* numLines);

#endif // ARRAY_OPS_H
