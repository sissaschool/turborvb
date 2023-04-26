/*
Copyright (C) 2022 TurboRVB group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdlib.h>
#include <stddef.h>
#include <ctype.h>
typedef size_t devptr_t;
struct complex
{
    double re;
    double im;
};
#ifdef _CUBLAS
#include <cublas.h>
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_dger_offload_(const int * M, const int * N, const double * alpha, const devptr_t * devPtrX, const int * incx, const devptr_t * devPtrY, const int * incy, const devptr_t * devPtrA, const int * lda)
{
    #pragma omp target data use_device_ptr(devPtrX, devPtrY, devPtrA)
    {
    double * devPtrX_ = (double *) devPtrX;
    double * devPtrY_ = (double *) devPtrY;
    double * devPtrA_ = (double *) devPtrA;
    cublasDger(* M, * N, * alpha, devPtrX_, * incx, devPtrY_, * incy, devPtrA_, * lda);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_zgeru_offload_(const int * M, const int * N, const cuDoubleComplex * alpha, const devptr_t * devPtrX, const int * incx, const devptr_t * devPtrY, const int * incy, const devptr_t * devPtrA, const int * lda)
{
    #pragma omp target data use_device_ptr(devPtrX, devPtrY, devPtrA)
    {
    cuDoubleComplex * devPtrX_ = (cuDoubleComplex *) devPtrX;
    cuDoubleComplex * devPtrY_ = (cuDoubleComplex *) devPtrY;
    cuDoubleComplex * devPtrA_ = (cuDoubleComplex *) devPtrA;
    cublasZgeru(* M, * N, * alpha, devPtrX_, * incx, devPtrY_, * incy, devPtrA_, * lda);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_sgemm_offload_(const char * transa, const char * transb, const int * M, const int * N, const int * K, const float * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrB, const int * ldb, const float * beta, const devptr_t * devPtrC, const int * ldc)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrB, devPtrC)
    {
    float * devPtrA_ = (float *) devPtrA;
    float * devPtrB_ = (float *) devPtrB;
    float * devPtrC_ = (float *) devPtrC;
    cublasSgemm(* transa, * transb, * M, * N, * K, * alpha, devPtrA_, * lda, devPtrB_, * ldb, * beta, devPtrC_, * ldc);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_dgemm_offload_(const char * transa, const char * transb, const int * M, const int * N, const int * K, const double * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrB, const int * ldb, const double * beta, const devptr_t * devPtrC, const int * ldc)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrB, devPtrC)
    {
    double * devPtrA_ = (double *) devPtrA;
    double * devPtrB_ = (double *) devPtrB;
    double * devPtrC_ = (double *) devPtrC;
    cublasDgemm(* transa, * transb, * M, * N, * K, * alpha, devPtrA_, * lda, devPtrB_, * ldb, * beta, devPtrC_, * ldc);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_zgemm_offload_(const char * transa, const char * transb, const int * M, const int * N, const int * K, const cuDoubleComplex * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrB, const int * ldb, const cuDoubleComplex * beta, const devptr_t * devPtrC, const int * ldc)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrB, devPtrC)
    {
    cuDoubleComplex * devPtrA_ = (cuDoubleComplex *) devPtrA;
    cuDoubleComplex * devPtrB_ = (cuDoubleComplex *) devPtrB;
    cuDoubleComplex * devPtrC_ = (cuDoubleComplex *) devPtrC;
    cublasZgemm(* transa, * transb, * M, * N, * K, * alpha, devPtrA_, * lda, devPtrB_, * ldb, * beta, devPtrC_, * ldc);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_sgemv_offload_(const char * trans, const int * M, const int * N, const float * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrX, const int * incx, const float * beta, const devptr_t * devPtrY, const int * incy)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrX, devPtrY)
    {
    float * devPtrA_ = (float *) devPtrA;
    float * devPtrX_ = (float *) devPtrX;
    float * devPtrY_ = (float *) devPtrY;
    cublasSgemv(* trans, * M, * N, * alpha, devPtrA_, * lda, devPtrX_, * incx, * beta, devPtrY_, * incy);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_dgemv_offload_(const char * trans, const int * M, const int * N, const double * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrX, const int * incx, const double * beta, const devptr_t * devPtrY, const int * incy)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrX, devPtrY)
    {
    double * devPtrA_ = (double *) devPtrA;
    double * devPtrX_ = (double *) devPtrX;
    double * devPtrY_ = (double *) devPtrY;
    cublasDgemv(* trans, * M, * N, * alpha, devPtrA_, * lda, devPtrX_, * incx, * beta, devPtrY_, * incy);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_zgemv_offload_(const char * trans, const int * M, const int * N, const cuDoubleComplex * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrX, const int * incx, const cuDoubleComplex * beta, const devptr_t * devPtrY, const int * incy)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrX, devPtrY)
    {
    cuDoubleComplex * devPtrA_ = (cuDoubleComplex *) devPtrA;
    cuDoubleComplex * devPtrX_ = (cuDoubleComplex *) devPtrX;
    cuDoubleComplex * devPtrY_ = (cuDoubleComplex *) devPtrY;
    cublasZgemv(* trans, * M, * N, * alpha, devPtrA_, * lda, devPtrX_, * incx, * beta, devPtrY_, * incy);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_dtrsm_offload_(const char * side, const char * uplo, const char * transa, const char * diag, const int * M, const int * N, const double * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrB, const int * ldb)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrB)
    {
    double * devPtrA_ = (double *) devPtrA;
    double * devPtrB_ = (double *) devPtrB;
    cublasDtrsm(* side, * uplo, * transa, * diag, * M, * N, * alpha, devPtrA_, * lda, devPtrB_, * ldb);
    }
}
#endif /* _CUBLAS */
#ifdef _CUBLAS
void cublas_ztrsm_offload_(const char * side, const char * uplo, const char * transa, const char * diag, const int * M, const int * N, const cuDoubleComplex * alpha, const devptr_t * devPtrA, const int * lda, const devptr_t * devPtrB, const int * ldb)
{
    #pragma omp target data use_device_ptr(devPtrA, devPtrB)
    {
    cuDoubleComplex * devPtrA_ = (cuDoubleComplex *) devPtrA;
    cuDoubleComplex * devPtrB_ = (cuDoubleComplex *) devPtrB;
    cublasZtrsm(* side, * uplo, * transa, * diag, * M, * N, * alpha, devPtrA_, * lda, devPtrB_, * ldb);
    }
}
#endif /* _CUBLAS */
#ifdef _CUSOLVER
#include "cusolverDn.h"
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_handle_init_(long int * handle)
{
    cusolverDnHandle_t *h;
    h = (cusolverDnHandle_t*) malloc(sizeof(cusolverDnHandle_t));
    cusolverDnCreate(h);
    *handle = (long int) h;
}
void cusolver_handle_destroy_(long int * handle)
{
    cusolverDnHandle_t *h = (cusolverDnHandle_t *) *handle;
    cusolverDnDestroy(*h);
    free(h);
}
void cudasync_(void)
{
    cudaDeviceSynchronize();
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dgetrf_buffersize_(const long int * handle,  int * status,  int * M,  int * N,  void * A,  int * lda,  int * workspace)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A)
    {
    void * A_ = (void *) A;
    *status = cusolverDnDgetrf_bufferSize(* handle__, * M, * N, A_, * lda, workspace);
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zgetrf_buffersize_(const long int * handle,  int * status,  int * M,  int * N,  void * A,  int * lda,  int * workspace)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A)
    {
    void * A_ = (void *) A;
    *status = cusolverDnZgetrf_bufferSize(* handle__, * M, * N, A_, * lda, workspace);
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dpotrf_buffersize_(const long int * handle,  int * status,  char * fillm,  int * M,  void * A,  int * lda,  int * workspace)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A)
    {
    void * A_ = (void *) A;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnDpotrf_bufferSize(* handle__, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, workspace);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnDpotrf_bufferSize(* handle__, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, workspace);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zpotrf_buffersize_(const long int * handle,  int * status,  char * fillm,  int * M,  void * A,  int * lda,  int * workspace)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A)
    {
    void * A_ = (void *) A;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnZpotrf_bufferSize(* handle__, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, workspace);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnZpotrf_bufferSize(* handle__, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, workspace);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dpotri_buffersize_(const long int * handle,  int * status,  char * fillm,  int * M,  void * A,  int * lda,  int * workspace)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A)
    {
    void * A_ = (void *) A;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnDpotri_bufferSize(* handle__, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, workspace);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnDpotri_bufferSize(* handle__, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, workspace);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zpotri_buffersize_(const long int * handle,  int * status,  char * fillm,  int * M,  void * A,  int * lda,  int * workspace)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A)
    {
    void * A_ = (void *) A;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnZpotri_bufferSize(* handle__, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, workspace);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnZpotri_bufferSize(* handle__, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, workspace);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dsyevd_buffersize_(const long int * handle,  int * status,  char * jobmode,  char * fillm,  int * M,  void * A,  int * lda,  void * W,  int * lwork)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, W)
    {
    void * A_ = (void *) A;
    void * W_ = (void *) W;
    if ((toupper(*fillm) == *"U") && (toupper(*jobmode) == *"N")) {
    *status = cusolverDnDsyevd_bufferSize(* handle__, CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, W_, lwork);
    }
    if ((toupper(*fillm) == *"L") && (toupper(*jobmode) == *"N")) {
    *status = cusolverDnDsyevd_bufferSize(* handle__, CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, W_, lwork);
    }
    if ((toupper(*fillm) == *"U") && (toupper(*jobmode) == *"V")) {
    *status = cusolverDnDsyevd_bufferSize(* handle__, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, W_, lwork);
    }
    if ((toupper(*fillm) == *"L") && (toupper(*jobmode) == *"V")) {
    *status = cusolverDnDsyevd_bufferSize(* handle__, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, W_, lwork);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dgetrf_(const long int * handle,  int * status,  int * M,  int * N,  void * A,  int * lda,  void * workspace,  void * devIpiv,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, workspace, devIpiv, devInfo)
    {
    double * A_ = (double *) A;
    double * workspace_ = (double *) workspace;
    int * devIpiv_ = (int *) devIpiv;
    int * devInfo_ = (int *) devInfo;
    *status = cusolverDnDgetrf(* handle__, * M, * N, A_, * lda, workspace_, devIpiv_, devInfo_);
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zgetrf_(const long int * handle,  int * status,  int * M,  int * N,  void * A,  int * lda,  void * workspace,  void * devIpiv,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, workspace, devIpiv, devInfo)
    {
    cuDoubleComplex * A_ = (cuDoubleComplex *) A;
    cuDoubleComplex * workspace_ = (cuDoubleComplex *) workspace;
    int * devIpiv_ = (int *) devIpiv;
    int * devInfo_ = (int *) devInfo;
    *status = cusolverDnZgetrf(* handle__, * M, * N, A_, * lda, workspace_, devIpiv_, devInfo_);
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dgetrs_(const long int * handle,  int * status,  char * trans,  int * M,  int * nrhs,  void * A,  int * lda,  void * devIpiv,  void * B,  int * ldb,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, devIpiv, B, devInfo)
    {
    double * A_ = (double *) A;
    int * devIpiv_ = (int *) devIpiv;
    double * B_ = (double *) B;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*trans) == *"N")) {
    *status = cusolverDnDgetrs(* handle__, CUBLAS_OP_N, * M, * nrhs, A_, * lda, devIpiv_, B_, * ldb, devInfo_);
    }
    if ((toupper(*trans) == *"L")) {
    *status = cusolverDnDgetrs(* handle__, CUBLAS_OP_T, * M, * nrhs, A_, * lda, devIpiv_, B_, * ldb, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zgetrs_(const long int * handle,  int * status,  char * trans,  int * M,  int * nrhs,  void * A,  int * lda,  void * devIpiv,  void * B,  int * ldb,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, devIpiv, B, devInfo)
    {
    cuDoubleComplex * A_ = (cuDoubleComplex *) A;
    int * devIpiv_ = (int *) devIpiv;
    cuDoubleComplex * B_ = (cuDoubleComplex *) B;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*trans) == *"N")) {
    *status = cusolverDnZgetrs(* handle__, CUBLAS_OP_N, * M, * nrhs, A_, * lda, devIpiv_, B_, * ldb, devInfo_);
    }
    if ((toupper(*trans) == *"L")) {
    *status = cusolverDnZgetrs(* handle__, CUBLAS_OP_T, * M, * nrhs, A_, * lda, devIpiv_, B_, * ldb, devInfo_);
    }
    if ((toupper(*trans) == *"C")) {
    *status = cusolverDnZgetrs(* handle__, CUBLAS_OP_C, * M, * nrhs, A_, * lda, devIpiv_, B_, * ldb, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dpotrs_(const long int * handle,  int * status,  char * fillm,  int * M,  int * nrhs,  void * A,  int * lda,  void * B,  int * ldb,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, B, devInfo)
    {
    double * A_ = (double *) A;
    double * B_ = (double *) B;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnDpotrs(* handle__, CUBLAS_FILL_MODE_UPPER, * M, * nrhs, A_, * lda, B_, * ldb, devInfo_);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnDpotrs(* handle__, CUBLAS_FILL_MODE_LOWER, * M, * nrhs, A_, * lda, B_, * ldb, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zpotrs_(const long int * handle,  int * status,  char * fillm,  int * M,  int * nrhs,  void * A,  int * lda,  void * B,  int * ldb,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, B, devInfo)
    {
    cuDoubleComplex * A_ = (cuDoubleComplex *) A;
    cuDoubleComplex * B_ = (cuDoubleComplex *) B;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnZpotrs(* handle__, CUBLAS_FILL_MODE_UPPER, * M, * nrhs, A_, * lda, B_, * ldb, devInfo_);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnZpotrs(* handle__, CUBLAS_FILL_MODE_LOWER, * M, * nrhs, A_, * lda, B_, * ldb, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dpotri_(const long int * handle,  int * status,  char * fillm,  int * M,  void * A,  int * lda,  void * workspace,  int * lwork,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, workspace, devInfo)
    {
    double * A_ = (double *) A;
    double * workspace_ = (double *) workspace;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnDpotri(* handle__, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, workspace_, * lwork, devInfo_);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnDpotri(* handle__, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, workspace_, * lwork, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_zpotri_(const long int * handle,  int * status,  char * fillm,  int * M,  void * A,  int * lda,  void * workspace,  int * lwork,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, workspace, devInfo)
    {
    cuDoubleComplex * A_ = (cuDoubleComplex *) A;
    cuDoubleComplex * workspace_ = (cuDoubleComplex *) workspace;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*fillm) == *"U")) {
    *status = cusolverDnZpotri(* handle__, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, workspace_, * lwork, devInfo_);
    }
    if ((toupper(*fillm) == *"L")) {
    *status = cusolverDnZpotri(* handle__, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, workspace_, * lwork, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */
#ifdef _CUSOLVER
void cusolver_dsyevd_(const long int * handle,  int * status,  char * jobmode,  char * fillm,  int * M,  void * A,  int * lda,  void * W,  void * workspace,  int * lwork,  void * devInfo)
{
    cusolverDnHandle_t * handle__ = (cusolverDnHandle_t *) *handle;
    #pragma omp target data use_device_ptr(A, W, workspace, devInfo)
    {
    void * A_ = (void *) A;
    double * W_ = (double *) W;
    double * workspace_ = (double *) workspace;
    int * devInfo_ = (int *) devInfo;
    if ((toupper(*fillm) == *"U") && (toupper(*jobmode) == *"N")) {
    *status = cusolverDnDsyevd(* handle__, CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, W_, workspace_, * lwork, devInfo_);
    }
    if ((toupper(*fillm) == *"L") && (toupper(*jobmode) == *"N")) {
    *status = cusolverDnDsyevd(* handle__, CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, W_, workspace_, * lwork, devInfo_);
    }
    if ((toupper(*fillm) == *"U") && (toupper(*jobmode) == *"V")) {
    *status = cusolverDnDsyevd(* handle__, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_UPPER, * M, A_, * lda, W_, workspace_, * lwork, devInfo_);
    }
    if ((toupper(*fillm) == *"L") && (toupper(*jobmode) == *"V")) {
    *status = cusolverDnDsyevd(* handle__, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER, * M, A_, * lda, W_, workspace_, * lwork, devInfo_);
    }
    }
}
#endif /* _CUSOLVER */