
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This file contains the common kernels shared among different BLAS functions. This file contains
// kernels to transpose matrices in various ways, including:
// 1) transposing into a larger matrix by adding padding
// 2) transposing into a smaller matrix by optionally removing padding. This is the general version
//    without restrictions, see the 'transpose.opencl' file for a faster but more restricted
//    transpose kernel.
//
// =================================================================================================

// Enables loading of this file using the C++ pre-processor's #include (C++11 standard raw string
// literal). Comment-out this line for syntax-highlighting when developing.
R"(

#if PADTRA_TILE == 1
  #define PADTRA_TILE_SHIFT 0
#elif PADTRA_TILE == 2
  #define PADTRA_TILE_SHIFT 1
#elif PADTRA_TILE == 4
  #define PADTRA_TILE_SHIFT 2
#elif PADTRA_TILE == 8
  #define PADTRA_TILE_SHIFT 3
#elif PADTRA_TILE == 16
  #define PADTRA_TILE_SHIFT 4
#elif PADTRA_TILE == 32
  #define PADTRA_TILE_SHIFT 5
#elif PADTRA_TILE == 64
  #define PADTRA_TILE_SHIFT 6
#elif PADTRA_TILE == 128
  #define PADTRA_TILE_SHIFT 7
#elif PADTRA_TILE == 256
  #define PADTRA_TILE_SHIFT 8
#endif

#if PADTRA_WPT == 1
  #define PADTRA_WPT_SHIFT 0
#elif PADTRA_WPT == 2
  #define PADTRA_WPT_SHIFT 1
#elif PADTRA_WPT == 4
  #define PADTRA_WPT_SHIFT 2
#elif PADTRA_WPT == 8
  #define PADTRA_WPT_SHIFT 3
#elif PADTRA_WPT == 16
  #define PADTRA_WPT_SHIFT 4
#elif PADTRA_WPT == 32
  #define PADTRA_WPT_SHIFT 5
#elif PADTRA_WPT == 64
  #define PADTRA_WPT_SHIFT 6
#elif PADTRA_WPT == 128
  #define PADTRA_WPT_SHIFT 7
#elif PADTRA_WPT == 256
  #define PADTRA_WPT_SHIFT 8
#endif

typedef union PADTRA_Ptr {
    __global singlereal        *f;
    __local  singlereal        *lf;
    __global const singlereal  *cf;
    __global singlereal2       *f2;
    __global const singlereal2 *cf2;
    __global real              *s;
    __local  real              *ls;
    __global const real        *cs;
} PADTRA_Ptr;



// =================================================================================================

// Transposes a matrix from source to destination. The output is padded with zero values in case the
// destination matrix dimensions are larger than the transposed source matrix dimensions.
__attribute__((reqd_work_group_size(PADTRA_TILE, PADTRA_TILE, 1)))
__kernel void TransposePadMatrix(const int src_one, const int src_two,
                                 const int src_ld, const int src_offset,
                                 __global const real* restrict src,
                                 const int dest_one, const int dest_two,
                                 const int dest_ld, const int dest_offset,
                                 __global real* dest,
                                 const __constant real* restrict arg_alpha,
                                 const int do_conjugate) {
  // const real alpha = arg_alpha[0];


  const uint LocalID0xPADTRA_WPT = (uint) get_local_id(0) << PADTRA_WPT_SHIFT;
  const uint LocalID1xPADTRA_WPT = (uint) get_local_id(1) << PADTRA_WPT_SHIFT;


  // Local memory to store a tile of the matrix (for coalescing)
  __local real tile[PADTRA_WPT<<PADTRA_TILE_SHIFT][(PADTRA_WPT<<PADTRA_TILE_SHIFT) + PADTRA_PAD];


  // Loop over the work per thread
  const uint GroupID0xPADTRA_WPT = ((uint) get_group_id(0) << (PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT)) + (uint) get_local_id(1);
  const uint GroupID1xPADTRA_WPT = ((uint) get_group_id(1) << (PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT)) + (uint) get_local_id(0);
  uint src_ldp = 0;
  const uint GroupID0xPADTRA_WPTxSRC_LD = GroupID0xPADTRA_WPT*src_ld + src_offset;

  #pragma unroll 
  for (uint w_one=0; w_one<PADTRA_WPT; ++w_one) {

    const uint id_src_two  = GroupID0xPADTRA_WPT + (w_one << PADTRA_TILE_SHIFT);
    const uint w_one_1    = LocalID0xPADTRA_WPT + w_one;
    const uint id_src_twoXsrc_ld = GroupID0xPADTRA_WPTxSRC_LD + (src_ldp << PADTRA_TILE_SHIFT);

    #pragma unroll 
    for (uint w_two=0; w_two<PADTRA_WPT; ++w_two) {

      const uint id_src_one = GroupID1xPADTRA_WPT + (w_two << PADTRA_TILE_SHIFT);

      bool condition = (id_src_two < src_two) && (id_src_one < src_one);
      #if (PRECISION == 32) || (PRECISION == 64)
        tile[LocalID1xPADTRA_WPT + w_two][w_one_1] = condition ? src[id_src_twoXsrc_ld + id_src_one] : (real) ZERO ;
      #else
        if (condition) {
          #if (USE_VLOAD == 1) 

            PADTRA_Ptr S;
            PADTRA_Ptr L;
            S.cs = &src[id_src_twoXsrc_ld + id_src_one];
            L.ls = &tile[LocalID1xPADTRA_WPT + w_two][w_one_1];
            vstore2(vload2(0,S.cf),0,L.lf);

          #else
            tile[LocalID1xPADTRA_WPT + w_two][w_one_1] = src[id_src_twoXsrc_ld + id_src_one];
          #endif

        } else {
             real value;
             SetToZero(value);
             tile[LocalID1xPADTRA_WPT + w_two][w_one_1] = value;
        }
      #endif

    }
    src_ldp += src_ld;

  }

  // Synchronizes all threads in a workgroup
  barrier(CLK_LOCAL_MEM_FENCE);

  const uint GroupID1xPADTRA_WPT_1 = ((uint) get_group_id(1) << (PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT)) + (uint) get_local_id(1);
  const uint GroupID0xPADTRA_WPT_1 = ((uint) get_group_id(0) << (PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT)) + (uint) get_local_id(0);
  const uint GroupID1xPADTRA_WPT_1xDEST_LD = GroupID1xPADTRA_WPT_1*dest_ld + dest_offset;

  // Loop over the work per thread
  #pragma unroll
  for (uint w_one=0; w_one<PADTRA_WPT; ++w_one) {

    const uint id_dest_one = GroupID0xPADTRA_WPT_1 + (w_one << PADTRA_TILE_SHIFT);
    const uint id_dest_onep = GroupID1xPADTRA_WPT_1xDEST_LD+id_dest_one;

    if ( id_dest_one < dest_one ) {

      const uint w_one_1 = LocalID1xPADTRA_WPT + w_one; 
      
      uint dest_ldp = 0;

      #pragma unroll
      for (uint w_two=0; w_two<PADTRA_WPT; ++w_two) {
  
        // Computes the identifiers for the destination matrix
  
        const int id_dest_two = GroupID1xPADTRA_WPT_1 + (w_two << PADTRA_TILE_SHIFT);
  
        // Stores the transposed value in the destination matrix
        if ( id_dest_two < dest_two ) {
  
          #if (PRECISION == 32) || (PRECISION == 64)

            dest[(dest_ldp << PADTRA_TILE_SHIFT) + id_dest_onep] = tile[LocalID0xPADTRA_WPT + w_two][w_one_1];

          #else  

            if (do_conjugate == 1) {
              real value = tile[LocalID0xPADTRA_WPT + w_two][w_one_1];
              COMPLEX_CONJUGATE(value);
              dest[(dest_ldp << PADTRA_TILE_SHIFT) + id_dest_onep] = value;
            } else {

              #if USE_VLOAD == 1

                PADTRA_Ptr L;
                PADTRA_Ptr D;
                L.ls = &tile[LocalID0xPADTRA_WPT + w_two][w_one_1];
                D.s  = &dest[(dest_ldp << PADTRA_TILE_SHIFT) + id_dest_onep];
                vstore2(vload2(0,L.lf),0,D.f);

              #else
                dest[(dest_ldp << PADTRA_TILE_SHIFT) + id_dest_onep] = tile[LocalID0xPADTRA_WPT + w_two][w_one_1];
              #endif
            }

          #endif
        }
        dest_ldp += dest_ld;
      }
    }
  }
}

// =================================================================================================

// Transposes a matrix, while considering possible padding in the source matrix. Data is read from a
// padded source matrix, but only the actual data is written back to the transposed destination
// matrix. This kernel optionally checks for upper/lower triangular matrices.
__attribute__((reqd_work_group_size(PADTRA_TILE, PADTRA_TILE, 1)))
__kernel void TransposeMatrix(const int src_one, const int src_two,
                              const int src_ld, const int src_offset,
                              __global const real* restrict src,
                              const int dest_one, const int dest_two,
                              const int dest_ld, const int dest_offset,
                              __global real* dest,
                              const __constant real* restrict arg_alpha,
                              const int upper, const int lower,
                              const int diagonal_imag_zero) {
  const real alpha = arg_alpha[0];


  const int GroupID1xPADTRA_WPT = (get_group_id(1) << ( PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT )) + get_local_id(0);
  const int GroupID0xPADTRA_WPT = (get_group_id(0) << ( PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT )) + get_local_id(1);

  const int GroupID1xPADTRA_WPT_1 = (get_group_id(1) << ( PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT )) + get_local_id(1);
  const int GroupID0xPADTRA_WPT_1 = (get_group_id(0) << ( PADTRA_WPT_SHIFT + PADTRA_TILE_SHIFT )) + get_local_id(0);

  const int LocalID0xPADTRA_WPT = get_local_id(0) << PADTRA_WPT_SHIFT;
  const int LocalID1xPADTRA_WPT = get_local_id(1) << PADTRA_WPT_SHIFT;

  // Local memory to store a tile of the matrix (for coalescing)
  __local real tile[PADTRA_WPT*PADTRA_TILE][PADTRA_WPT*PADTRA_TILE + PADTRA_PAD];

    // Loop over the work per thread
  #pragma unroll PADTRA_WPT
  for (int w_one=0; w_one<PADTRA_WPT; ++w_one) {
  
    const int id_src_two = GroupID0xPADTRA_WPT + ( w_one << PADTRA_TILE_SHIFT ) ;
  
    if ( id_src_two < src_two ) {

      const int w_one_1 = LocalID0xPADTRA_WPT + w_one;

      #pragma unroll PADTRA_WPT
      for (int w_two=0; w_two<PADTRA_WPT; ++w_two) {
  
        // Computes the identifiers for the source matrix. Note that the local and global dimensions
        // do not correspond to each other!
  
        const int id_src_one = GroupID1xPADTRA_WPT + (w_two << PADTRA_TILE_SHIFT);
  
        // Loads data into the local memory if the thread IDs are within bounds of the source matrix.
        if (id_src_one < src_one) {
  
            #if USE_MAD24 == 1
              real value = src[mad24(id_src_two,src_ld , id_src_one + src_offset)];
            #else
              real value = src[id_src_two*src_ld + id_src_one + src_offset];
            #endif

            tile[LocalID1xPADTRA_WPT + w_two][w_one_1] = value;
  
  
        }
      }
    }
  }

  // Synchronizes all threads in a workgroup
  barrier(CLK_LOCAL_MEM_FENCE);

  // Loop over the work per thread
  #pragma unroll PADTRA_WPT
  for (int w_one=0; w_one<PADTRA_WPT; ++w_one) {

    const int id_dest_one = GroupID0xPADTRA_WPT_1 + (w_one << PADTRA_TILE_SHIFT);

    if ( id_dest_one < dest_one ) {

      const int w_one_1 = LocalID1xPADTRA_WPT + w_one; 

      #pragma unroll PADTRA_WPT
      for (int w_two=0; w_two<PADTRA_WPT; ++w_two) {
  
        // Computes the identifiers for the destination matrix
  
        const int id_dest_two = GroupID1xPADTRA_WPT_1 + (w_two << PADTRA_TILE_SHIFT);
  
  
        // Masking in case of triangular matrices: updates only the upper or lower part
        bool condition = true;
  
        #if defined(ROUTINE_SYRK) || defined(ROUTINE_HERK) || defined(ROUTINE_SYR2K) || defined(ROUTINE_HER2K)
          if (upper == 1) { condition = (id_dest_one >= id_dest_two); }
          else if (lower == 1) { condition = (id_dest_one <= id_dest_two); }
        #endif
  
        if (condition) {
  
          // Stores the transposed value in the destination matrix
          if ( id_dest_two < dest_two ) {
  
              real value = tile[LocalID0xPADTRA_WPT + w_two][w_one_1];
  
              #if defined(ROUTINE_SYRK) || defined(ROUTINE_HERK) || defined(ROUTINE_SYR2K) || defined(ROUTINE_HER2K)
                if (diagonal_imag_zero == 1 && id_dest_one == id_dest_two) { ImagToZero(value); }
              #endif
  
              #if (PRECISION == 32) || (PRECISION == 3232) || (PRECISION == 64) || (PRECISION == 6464)
                #if USE_MAD24 == 1
                  dest[mad24(id_dest_two,dest_ld , id_dest_one + dest_offset)] = value;
                #else
                  dest[id_dest_two*dest_ld + id_dest_one + dest_offset] = value;
                #endif
              #else
                #if USE_MAD24 == 1
                  Multiply(dest[mad24(id_dest_two,dest_ld , id_dest_one + dest_offset)], alpha, value);
                #else
                  Multiply(dest[id_dest_two*dest_ld + id_dest_one + dest_offset], alpha, value);
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
