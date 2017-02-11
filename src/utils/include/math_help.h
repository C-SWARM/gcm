#ifndef H__H__MATH_HELP__H__H
#define H__H__MATH_HELP__H__H

#define DEBUG_PRINT_ERROR 0

#ifdef __cplusplus
extern "C" {
#endif /* #ifdef __cplusplus */

int inverse(double const* A, const int M, double *A_I);
void symmetric_part(double *sym, const double *mat, const int dim);

#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif

