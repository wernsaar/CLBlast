
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

typedef union PAD_Ptr {
    __global singlereal        *f;
    __global const singlereal  *cf;
    __global singlereal2       *f2;
    __global const singlereal2 *cf2;
    __global real              *s;
    __global const real        *cs;
} PAD_Ptr;



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
  // const real alpha = arg_alpha[0];

  const uint GroupID1xPAD_WPTY = ((uint) get_group_id(1) << (PAD_WPTY_SHIFT+PAD_DIMY_SHIFT))+get_local_id(1);
  const uint GroupID0xPAD_WPTX = ((uint) get_group_id(0) << (PAD_WPTX_SHIFT+PAD_DIMX_SHIFT))+get_local_id(0);
  const uint GroupID1xPAD_WPTYxSRC_LD = GroupID1xPAD_WPTY*src_ld;
  const uint GroupID1xPAD_WPTYxDEST_LD = GroupID1xPAD_WPTY*dest_ld;

  // Loops over the work per thread in both dimensions

  #pragma unroll 
  for (uint w_one=0; w_one<PAD_WPTX; ++w_one) {

    const uint id_one = GroupID0xPAD_WPTX + (w_one << PAD_DIMX_SHIFT) ;
    const uint id_onePLUSdest_offset = id_one + dest_offset+GroupID1xPAD_WPTYxDEST_LD;
    const uint id_onePLUSsrc_offset = id_one + src_offset+GroupID1xPAD_WPTYxSRC_LD;

    if ( id_one < dest_one ) {

        #if PAD_WPTY > 1
          uint src_ldp = 0;
          uint dest_ldp = 0;
        #endif

        #pragma unroll
        for (uint w_two=0; w_two<PAD_WPTY; ++w_two) {
  
          const uint id_two = GroupID1xPAD_WPTY + (w_two << PAD_DIMY_SHIFT);

          #if (USE_VLOAD == 1) && ((PRECISION == 3232) || (PRECISION == 6464))

            if (id_two < dest_two) {
  
              singlereal2 value;
	      PAD_Ptr SRC;
	      PAD_Ptr DEST;

              bool condition = (id_one < src_one) && (id_two < src_two) ;
              if (condition) {
                #if PAD_WPTY > 1
                  SRC.cs =  &src[(src_ldp<<PAD_DIMY_SHIFT) + id_onePLUSsrc_offset];
                  value = vload2(0,SRC.cf);
                #else
                  SRC.cs =  &src[id_onePLUSsrc_offset];
                  value = vload2(0,SRC.cf);
                #endif
                #if (PRECISION == 3232) || (PRECISION == 6464) 
                  if (do_conjugate == 1) { value.y = -value.y; }
                #endif

              } else {
                value = (singlereal2) ZERO ;
              }
  
              // Stores the value in the destination matrix
              #if PAD_WPTY > 1  
                DEST.s = &dest[(dest_ldp<<PAD_DIMY_SHIFT) + id_onePLUSdest_offset];
                vstore2(value,0,DEST.f);
                src_ldp  += src_ld;
                dest_ldp += dest_ld;
              #else
                DEST.s = &dest[id_onePLUSdest_offset];
                vstore2(value,0,DEST.f);
              #endif

          #else

            if (id_two < dest_two) {
  
              // Loads data if the thread IDs are within bounds of the source matrix. Otherwise, set the
              // value to be written to zero.
              real value;
              bool condition = (id_one < src_one) && (id_two < src_two) ;
              if (condition) {
                #if PAD_WPTY > 1
                  value = src[(src_ldp<<PAD_DIMY_SHIFT) + id_onePLUSsrc_offset];
                #else
                  value = src[id_onePLUSsrc_offset];
                #endif
                #if (PRECISION == 3232) || (PRECISION == 6464) 
                  if (do_conjugate == 1) { COMPLEX_CONJUGATE(value); }
                #endif

              } else {
                SetToZero(value);
              }
  
              // Stores the value in the destination matrix
              #if PAD_WPTY > 1  
                dest[(dest_ldp<<PAD_DIMY_SHIFT) + id_onePLUSdest_offset] = value;
                src_ldp  += src_ld;
                dest_ldp += dest_ld;
              #else
                dest[id_onePLUSdest_offset] = value;
              #endif

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
