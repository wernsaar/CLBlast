
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This file contains the common kernels shared among different BLAS routines. This file contains
// kernels to copy matrices.
//
// =================================================================================================

// Enables loading of this file using the C++ pre-processor's #include (C++11 standard raw string
// literal). Comment-out this line for syntax-highlighting when developing.
R"(

// =================================================================================================

// Data-widths
#if COPY_VW == 1
  #define COPY_VW_SHIFT 0
  typedef real realC;
#elif COPY_VW == 2
  #define COPY_VW_SHIFT 1
  typedef real2 realC;
#elif COPY_VW == 4
  #define COPY_VW_SHIFT 2
  typedef real4 realC;
#elif COPY_VW == 8
  #define COPY_VW_SHIFT 3
  typedef real8 realC;
#elif COPY_VW == 16
  #define COPY_VW_SHIFT 4
  typedef real16 realC;
#endif


#if COPY_WPT == 1
  #define COPY_WPT_SHIFT 0
#elif COPY_WPT == 2
  #define COPY_WPT_SHIFT 1
#elif COPY_WPT == 4
  #define COPY_WPT_SHIFT 2
#elif COPY_WPT == 8
  #define COPY_WPT_SHIFT 3
#elif COPY_WPT == 16
  #define COPY_WPT_SHIFT 4
#endif

#if COPY_DIMY == 1
  #define COPY_DIMY_SHIFT 0
#elif COPY_DIMY == 2
  #define COPY_DIMY_SHIFT 1
#elif COPY_DIMY == 4
  #define COPY_DIMY_SHIFT 2
#elif COPY_DIMY == 8
  #define COPY_DIMY_SHIFT 3
#elif COPY_DIMY == 16
  #define COPY_DIMY_SHIFT 4
#elif COPY_DIMY == 32
  #define COPY_DIMY_SHIFT 5
#elif COPY_DIMY == 64
  #define COPY_DIMY_SHIFT 6
#endif



// =================================================================================================

// Fast copy kernel. Requires 'ld' and the number of threads in dimension 0 to be a multiple of
// COPY_VW. Also requires both matrices to be of the same dimensions and without offset.
__attribute__((reqd_work_group_size(COPY_DIMX, COPY_DIMY, 1)))
__kernel void CopyMatrixFast(const int ld,
                             __global const realC* restrict src,
                             __global realC* dest,
                             const __constant real* restrict arg_alpha) {


  #if USE_MAD24 == 1
    uint GroupID1_M2 = mad24(((int) get_group_id(1) << (COPY_WPT_SHIFT + COPY_DIMY_SHIFT)) + (int) get_local_id(1), (ld >> COPY_VW_SHIFT) , (int) get_global_id(0)) ;
  #else
    uint GroupID1_M2 = (((uint) get_group_id(1) << (COPY_WPT_SHIFT + COPY_DIMY_SHIFT)) + (uint) get_local_id(1)) * (ld >> COPY_VW_SHIFT) + (uint) get_global_id(0) ;
  #endif

  uint ld_D1 = (ld >> COPY_VW_SHIFT);
   
  uint ld_D1P = 0;

  #pragma unroll COPY_WPT 
  for (uint w_one=0; w_one<COPY_WPT; w_one++) {

    uint id = ld_D1P + GroupID1_M2;
    #if COPY_WPT > 1
      ld_D1P += ld_D1 << COPY_DIMY_SHIFT;
    #endif

    dest[id] = src[id];

    }

}

// =================================================================================================

// End of the C++11 raw string literal
)"

// =================================================================================================
