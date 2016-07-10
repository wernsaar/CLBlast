
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This file contains the common kernels shared among different BLAS functions. This file contains
// kernels to copy and pad matrices in various ways, including:
// 1) copying into a larger matrix by adding padding
// 2) copying into a smaller matrix by optionally removing padding. This is the general version
//    without restrictions, see the 'copy.opencl' file for a faster but more restricted copy kernel.
//
// =================================================================================================

// Enables loading of this file using the C++ pre-processor's #include (C++11 standard raw string
// literal). Comment-out this line for syntax-highlighting when developing.
R"(

#if PAD_DIMY == 1
  #define PAD_DIMY_SHIFT 0
#elif PAD_DIMY == 2
  #define PAD_DIMY_SHIFT 1
#elif PAD_DIMY == 4
  #define PAD_DIMY_SHIFT 2
#elif PAD_DIMY == 8
  #define PAD_DIMY_SHIFT 3
#elif PAD_DIMY == 16
  #define PAD_DIMY_SHIFT 4
#elif PAD_DIMY == 32
  #define PAD_DIMY_SHIFT 5
#elif PAD_DIMY == 64
  #define PAD_DIMY_SHIFT 6
#elif PAD_DIMY == 128
  #define PAD_DIMY_SHIFT 7
#elif PAD_DIMY == 256
  #define PAD_DIMY_SHIFT 8
#endif

#if PAD_DIMX == 1
  #define PAD_DIMX_SHIFT 0
#elif PAD_DIMX == 2
  #define PAD_DIMX_SHIFT 1
#elif PAD_DIMX == 4
  #define PAD_DIMX_SHIFT 2
#elif PAD_DIMX == 8
  #define PAD_DIMX_SHIFT 3
#elif PAD_DIMX == 16
  #define PAD_DIMX_SHIFT 4
#elif PAD_DIMX == 32
  #define PAD_DIMX_SHIFT 5
#elif PAD_DIMX == 64
  #define PAD_DIMX_SHIFT 6
#elif PAD_DIMX == 128
  #define PAD_DIMX_SHIFT 7
#elif PAD_DIMX == 256
  #define PAD_DIMX_SHIFT 8
#endif


#if PAD_WPTX == 1
  #define PAD_WPTX_SHIFT 0
#elif PAD_WPTX == 2
  #define PAD_WPTX_SHIFT 1
#elif PAD_WPTX == 4
  #define PAD_WPTX_SHIFT 2
#elif PAD_WPTX == 8
  #define PAD_WPTX_SHIFT 3
#elif PAD_WPTX == 16
  #define PAD_WPTX_SHIFT 4
#elif PAD_WPTX == 32
  #define PAD_WPTX_SHIFT 5
#elif PAD_WPTX == 64
  #define PAD_WPTX_SHIFT 6
#elif PAD_WPTX == 128
  #define PAD_WPTX_SHIFT 7
#elif PAD_WPTX == 256
  #define PAD_WPTX_SHIFT 8
#endif

#if PAD_WPTY == 1
  #define PAD_WPTY_SHIFT 0
#elif PAD_WPTY == 2
  #define PAD_WPTY_SHIFT 1
#elif PAD_WPTY == 4
  #define PAD_WPTY_SHIFT 2
#elif PAD_WPTY == 8
  #define PAD_WPTY_SHIFT 3
#elif PAD_WPTY == 16
  #define PAD_WPTY_SHIFT 4
#elif PAD_WPTY == 32
  #define PAD_WPTY_SHIFT 5
#elif PAD_WPTY == 64
  #define PAD_WPTY_SHIFT 6
#elif PAD_WPTY == 128
  #define PAD_WPTY_SHIFT 7
#elif PAD_WPTY == 256
  #define PAD_WPTY_SHIFT 8
#endif




// =================================================================================================

// Copies a matrix from source to destination. The output is padded with zero values in case the
// destination matrix dimensions are larger than the source matrix dimensions. Additionally, the ld
// value and offset can be different.
__attribute__((reqd_work_group_size(PAD_DIMX, PAD_DIMY, 1)))
__kernel void CopyPadMatrix(const int src_one, const int src_two,
                            const int src_ld, const int src_offset,
                            __global const real* restrict src,
                            const int dest_one, const int dest_two,
                            const int dest_ld, const int dest_offset,
                            __global real* dest,
                            const __constant real* restrict arg_alpha,
                            const int do_conjugate) {
  const real alpha = arg_alpha[0];

  int GroupID0xPAD_WPTX = (int) get_group_id(0) << PAD_WPTX_SHIFT;
  int GroupID1xPAD_WPTY = (int) get_group_id(1) << PAD_WPTY_SHIFT;

  // Loops over the work per thread in both dimensions

  #pragma unroll
  for (int w_one=0; w_one<PAD_WPTX; ++w_one) {

    const int id_one = ((GroupID0xPAD_WPTX + w_one) << PAD_DIMX_SHIFT) + get_local_id(0);

    if ( id_one < dest_one ) {

      #pragma unroll
      for (int w_two=0; w_two<PAD_WPTY; ++w_two) {
  
        const int id_two = ((GroupID1xPAD_WPTY + w_two) << PAD_DIMY_SHIFT) + get_local_id(1);
  
        if (id_two < dest_two) {
  
          // Loads data if the thread IDs are within bounds of the source matrix. Otherwise, set the
          // value to be written to zero.
          real value;
          SetToZero(value);
          bool condition = id_two < src_two && id_one < src_one;
          if (condition) {
  
              value = src[id_two*src_ld + id_one + src_offset];
  
          }
  
          // Stores the value in the destination matrix
  
          if (do_conjugate == 1) { COMPLEX_CONJUGATE(value); }
  
          #if PRECISION == 32 
            dest[id_two*dest_ld + id_one + dest_offset] = value;
          #else
            Multiply(dest[id_two*dest_ld + id_one + dest_offset], alpha, value);
          #endif
        }
      }
    }
  }
}

// =================================================================================================

// Same as above, but now un-pads a matrix. This kernel reads data from a padded source matrix, but
// writes only the actual data back to the destination matrix. Again, the ld value and offset can
// be different.
__attribute__((reqd_work_group_size(PAD_DIMX, PAD_DIMY, 1)))
__kernel void CopyMatrix(const int src_one, const int src_two,
                         const int src_ld, const int src_offset,
                         __global const real* restrict src,
                         const int dest_one, const int dest_two,
                         const int dest_ld, const int dest_offset,
                         __global real* dest,
                         const __constant real* restrict arg_alpha,
                         const int upper, const int lower,
                         const int diagonal_imag_zero) {
  const real alpha = arg_alpha[0];

  int GroupID0xPAD_WPTX = get_group_id(0) << PAD_WPTX_SHIFT;
  int GroupID1xPAD_WPTY = get_group_id(1) << PAD_WPTY_SHIFT;

  // Loops over the work per thread in both dimensions
  #pragma unroll
  for (int w_one=0; w_one<PAD_WPTX; ++w_one) {

    const int id_one = ((GroupID0xPAD_WPTX + w_one) << PAD_DIMX_SHIFT) + get_local_id(0);

    if ( id_one < dest_one ) {

      #pragma unroll
      for (int w_two=0; w_two<PAD_WPTY; ++w_two) {
  
        const int id_two = ((GroupID1xPAD_WPTY + w_two) << PAD_DIMY_SHIFT) + get_local_id(1);
  
        // Masking in case of triangular matrices: updates only the upper or lower part
        bool condition = true;
  
        #if defined(ROUTINE_SYRK) || defined(ROUTINE_HERK) || defined(ROUTINE_SYR2K) || defined(ROUTINE_HER2K)
          if (upper == 1) { condition = (id_two >= id_one); }
          else if (lower == 1) { condition = (id_two <= id_one); }
        #endif
  
        if (condition) {
  
          // Copies the value into the destination matrix. This is always within bounds of the source
          // matrix, as we know that the destination matrix is smaller or equal to the source.
  
          if ( id_two < dest_two ) {
  
            #if USE_MAD24 == 1
              real value = src[mad24(id_two,src_ld , id_one + src_offset)];
            #else 
              real value = src[id_two*src_ld + id_one + src_offset];
            #endif

            bool condition2 = diagonal_imag_zero == 1 && id_one == id_two;
            if (condition2) { ImagToZero(value); }

            #if (PRECISION == 32) || (PRECISION == 3232) || (PRECISION == 64) || (PRECISION == 6464)
              #if USE_MAD24 == 1
                dest[mad24(id_two,dest_ld , id_one + dest_offset)] = value;
              #else
                dest[id_two*dest_ld + id_one + dest_offset] = value;
              #endif
            #else
              #if USE_MAD24 == 1
                Multiply(dest[mad24(id_two,dest_ld , id_one + dest_offset)], alpha, value);
              #else
                Multiply(dest[id_two*dest_ld + id_one + dest_offset], alpha, value);
              #endif
            #endif
  
          }
        }
      }
    }
  }
}

// =================================================================================================

// End of the C++11 raw string literal
)"

// =================================================================================================
