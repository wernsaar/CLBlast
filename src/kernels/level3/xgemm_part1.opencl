
// =================================================================================================
// This file is part of the CLBlast project. The project is licensed under Apache Version 2.0. This
// project loosely follows the Google C++ styleguide and uses a tab-size of two spaces and a max-
// width of 100 characters per line.
//
// Author(s):
//   Cedric Nugteren <www.cedricnugteren.nl>
//
// This file contains an optimized matrix-multiplication kernel inspired by the paper by Matsumoto
// et al. and the tutorial on http://www.cedricnugteren.nl/tutorial.php. It is fully configurable
// (and tunable!) using more or less the same parameters/naming conventions as in the paper. It
// supports different data-types (SGEMM/DGEMM/CGEMM/ZGEMM/HGEMM) through a pre-processor define.
//
// Matrices are accessed as follows:
// A: [k*M + m], with 'k' ranging from 0:K and 'm' from 0:M (m,k,m)
// B: [k*N + n], with 'k' ranging from 0:K and 'n' from 0:N (n,k,n)
// C: [n*M + m], with 'n' ranging from 0:N and 'm' from 0:M (m,n,m)
//
// Or as an image (assuming column-major)
//       K                      
//    o-------o                 
//    |       |                 
//  N | [B^T] |                 
//    |       |                 
//    o-------o                 
//        K               N     
//    o-------o        o-----o  
//  M |  [A]  |      M | [C] |  
//    |       |        |     |  
//    o-------o        o-----o  
//                              
//
// This kernel is seperated into two files. This is part 1 out of 2.
//
// =================================================================================================

// Enables loading of this file using the C++ pre-processor's #include (C++11 standard raw string
// literal). Comment-out this line for syntax-highlighting when developing.
R"(

// =================================================================================================

// Parameters set by the tuner or by the database. Here they are given a basic default value in case
// this kernel file is used outside of the CLBlast library.
#ifndef MWG
  #define MWG 8      // Tile-size in dimension M (e.g. 64, 128)
#endif
#ifndef NWG
  #define NWG 8      // Tile-size in dimension N (e.g. 64, 128)
#endif
#ifndef KWG
  #define KWG 8      // Tile-size in dimension K (e.g. 8, 16)
#endif
#ifndef MDIMC
  #define MDIMC 8    // Threads per workgroup in M-dimension (e.g. 8, 16, 32)
#endif
#ifndef NDIMC
  #define NDIMC 8    // Threads per workgroup in N-dimension (e.g. 8, 16, 32)
#endif
#ifndef MDIMA
  #define MDIMA 8    // Re-shaped tile dimension of matrix A: KDIMA * MDIMA
#endif
#ifndef NDIMB
  #define NDIMB 8    // Re-shaped tile dimension of matrix B: KDIMB * NDIMB
#endif
#ifndef KWI
  #define KWI 1      // Unroll factor of the KWG loop (smaller or equal than KWG)
#endif
#ifndef VWM
  #define VWM 1      // Vector width of matrices A and C
#endif
#ifndef VWN
  #define VWN 1      // Vector width of matrix B
#endif
#ifndef STRM
  #define STRM 0     // Use strided access within a thread in the M-dimension (1) or not (0)
#endif
#ifndef STRN
  #define STRN 0     // Use strided access within a thread in the N-dimension (1) or not (0)
#endif
#ifndef SA
  #define SA 0       // Use local/shared memory to cache matrix A (1) or not (0)
#endif
#ifndef SB
  #define SB 0       // Use local/shared memory to cache matrix B (1) or not (0)
#endif



// Helper parameters based on the above tuning parameters
#define MWI (MWG/MDIMC)               // Work per work-item (M-dimension)
#define NWI (NWG/NDIMC)               // Work per work-item (N-dimension)
#define KDIMA ((MDIMC*NDIMC)/(MDIMA)) // Re-shaped tile dimension of matrix A: KDIMA * MDIMA
#define KDIMB ((MDIMC*NDIMC)/(NDIMB)) // Re-shaped tile dimension of matrix B: KDIMB * NDIMB
#define MWA (MWG/MDIMA)               // Amount of loads-per-thread for matrix A (M-dimension)
#define KWA (KWG/KDIMA)               // Amount of loads-per-thread for matrix A (K-dimension)
#define KWB (KWG/KDIMB)               // Amount of loads-per-thread for matrix B (K-dimension)
#define NWB (NWG/NDIMB)               // Amount of loads-per-thread for matrix B (N-dimension)

// Settings
#ifndef USE_VECTOR_MAD
  #define USE_VECTOR_MAD 0      // Unroll (0) or don't (1) unroll the vector MAD manually
#endif
#ifndef GLOBAL_MEM_FENCE
  #define GLOBAL_MEM_FENCE 0    // Global synchronisation barrier for potential better performance
#endif

// =================================================================================================

// Data-widths in dimension M
#if VWM == 1
    #define VWM_SHIFT 0
    typedef real realM;
#elif VWM == 2
    #define VWM_SHIFT 1
    typedef real2 realM;
#elif VWM == 4
    #define VWM_SHIFT 2
    typedef real4 realM;
#elif VWM == 8
    #define VWM_SHIFT 3
    typedef real8 realM;
#elif VWM == 16
    #define VWM_SHIFT 4
    typedef real16 realM;
#endif

// Data-widths in dimension N
#if VWN == 1
    #define VWN_SHIFT 0
    typedef real realN;
#elif VWN == 2
    #define VWN_SHIFT 1
    typedef real2 realN;
#elif VWN == 4
    #define VWN_SHIFT 2
    typedef real4 realN;
#elif VWN == 8
    #define VWN_SHIFT 3
    typedef real8 realN;
#elif VWN == 16
    #define VWN_SHIFT 4
    typedef real16 realN;
#endif

#if MDIMA == 1
    #define MDIMA_SHIFT 0
#elif MDIMA == 2
    #define MDIMA_SHIFT 1
#elif MDIMA == 4
    #define MDIMA_SHIFT 2
#elif MDIMA == 8
    #define MDIMA_SHIFT 3
#elif MDIMA == 16
    #define MDIMA_SHIFT 4
#elif MDIMA == 32
    #define MDIMA_SHIFT 5
#elif MDIMA == 64
    #define MDIMA_SHIFT 6
#endif

#if NDIMB == 1
    #define NDIMB_SHIFT 0
#elif NDIMB == 2
    #define NDIMB_SHIFT 1
#elif NDIMB == 4
    #define NDIMB_SHIFT 2
#elif NDIMB == 8
    #define NDIMB_SHIFT 3
#elif NDIMB == 16
    #define NDIMB_SHIFT 4
#elif NDIMB == 32
    #define NDIMB_SHIFT 5
#elif NDIMB == 64
    #define NDIMB_SHIFT 6
#endif



// =================================================================================================

// Initializes the accumulation registers to zero
inline void InitAccRegisters(realM cpm[NWI][MWI >> VWM_SHIFT]) {
  #pragma unroll
  for (int mi=0; mi<(MWI>>VWM_SHIFT); ++mi) {
    #pragma unroll
    for (int ni=0; ni<NWI; ++ni) {
      #if VWM == 1
        SetToZero(cpm[ni][mi]);
      #elif VWM == 2
        SetToZero(cpm[ni][mi].x);
        SetToZero(cpm[ni][mi].y);
      #elif VWM == 4
        SetToZero(cpm[ni][mi].x);
        SetToZero(cpm[ni][mi].y);
        SetToZero(cpm[ni][mi].z);
        SetToZero(cpm[ni][mi].w);
      #elif VWM == 8
        SetToZero(cpm[ni][mi].s0);
        SetToZero(cpm[ni][mi].s1);
        SetToZero(cpm[ni][mi].s2);
        SetToZero(cpm[ni][mi].s3);
        SetToZero(cpm[ni][mi].s4);
        SetToZero(cpm[ni][mi].s5);
        SetToZero(cpm[ni][mi].s6);
        SetToZero(cpm[ni][mi].s7);
      #elif VWM == 16
        SetToZero(cpm[ni][mi].s0);
        SetToZero(cpm[ni][mi].s1);
        SetToZero(cpm[ni][mi].s2);
        SetToZero(cpm[ni][mi].s3);
        SetToZero(cpm[ni][mi].s4);
        SetToZero(cpm[ni][mi].s5);
        SetToZero(cpm[ni][mi].s6);
        SetToZero(cpm[ni][mi].s7);
        SetToZero(cpm[ni][mi].s8);
        SetToZero(cpm[ni][mi].s9);
        SetToZero(cpm[ni][mi].sA);
        SetToZero(cpm[ni][mi].sB);
        SetToZero(cpm[ni][mi].sC);
        SetToZero(cpm[ni][mi].sD);
        SetToZero(cpm[ni][mi].sE);
        SetToZero(cpm[ni][mi].sF);
      #endif
    }
  }
}

// =================================================================================================

// Caches global off-chip memory into local (shared) memory on-chip. This function is specific for
// caching the A input matrix.
#if SA == 1
inline void GlobalToLocalA(const __global realM* restrict agm, __local realM* alm,
                           const int kSizeM, const int tid, const int kwg) {
  const int la0 = tid - (tid & -MDIMA);
  const int la1 = tid >> MDIMA_SHIFT;

  #pragma unroll
  for (int mia=0; mia<(MWA >> VWM_SHIFT); ++mia) {

    #if USE_CL_MAD == 1

        // Computes the indices based on strided/non-strided access

        #if STRM == 0
          int mg = mad24(la0,(MWA >> VWM_SHIFT),mia);
        #elif STRM == 1
          int mg = mad24(mia,MDIMA,la0);
        #endif

    #else

        #if STRM == 0
          int mg = mia + la0*(MWA >> VWM_SHIFT);
        #elif STRM == 1
          int mg = la0 + mia*MDIMA;
        #endif

    #endif

    #pragma unroll
    for (int kia=0; kia<KWA; ++kia) {

      #if USE_CL_MAD == 1
        // Computes the indices for the global memory
        int kg =  mad24(la1,KWA,kia);
        int idm = mad24((int) GetGroupID0() , (MWG >> VWM_SHIFT),mg);
        int idk = kg + kwg;

        // Loads the data from global memory (not transposed) into the local memory
        alm[mad24(kg,(MWG/VWM) , mg)] = agm[mad24(idk,(kSizeM >> VWM_SHIFT) , idm)];

     #else
        // Computes the indices for the global memory
        int kg = kia + la1*KWA;
        int idm = mg + GetGroupID0() * (MWG >> VWM_SHIFT);
        int idk = kg + kwg;

        // Loads the data from global memory (not transposed) into the local memory
        alm[kg*(MWG/VWM) + mg] = agm[idk*(kSizeM >> VWM_SHIFT) + idm];


     #endif
    }
  }
}
#endif

// Same as above, but now for the B input matrix
#if SB == 1
inline void GlobalToLocalB(const __global realN* restrict bgm, __local realN* blm,
                           const int kSizeN, const int tid, const int kwg) {
  const int lb0 = tid - (tid & -NDIMB);
  const int lb1 = tid >> NDIMB_SHIFT;

  #pragma unroll
  for (int kib=0; kib<KWB; ++kib) {

    #if USE_CL_MAD == 1
      int kg = mad24(lb1,KWB,kib);
      int idk = kg + kwg;
    #else
      int kg = kib + lb1*KWB;
      int idk = kg + kwg;
    #endif

    #pragma unroll
    for (int nib=0; nib<(NWB >> VWN_SHIFT); ++nib) {

      #if USE_CL_MAD == 1
        // Computes the indices based on strided/non-strided access
        #if STRN == 0
          int ng = mad24(lb0,(NWB >> VWN_SHIFT),nib);
        #elif STRN == 1
          int ng = mad24(nib,NDIMB,lb0);
        #endif

        // Computes the indices for the global memory
        int idn = mad24((int) GetGroupID1() , (NWG >> VWN_SHIFT),ng);

        // Loads the data from global memory (transposed) into the local memory
        blm[mad24(kg,(NWG >> VWN_SHIFT) , ng)] = bgm[mad24(idk,(kSizeN >> VWN_SHIFT) , idn)];

     #else

        // Computes the indices based on strided/non-strided access
        #if STRN == 0
          int ng = nib + lb0*(NWB >> VWN_SHIFT);
        #elif STRN == 1
          int ng = lb0 + nib*NDIMB;
        #endif

        // Computes the indices for the global memory
        int idn = ng + GetGroupID1() * (NWG >> VWN_SHIFT);

        // Loads the data from global memory (transposed) into the local memory
        blm[kg*(NWG >> VWN_SHIFT) + ng] = bgm[idk*(kSizeN >> VWN_SHIFT) + idn];

     #endif


    }
  }
}
#endif

// =================================================================================================

// Caches global off-chip memory directly into per-thread private memory (registers). This function
// is specific for caching the A input matrix.
#if SA == 0
inline void GlobalToPrivateA(const __global realM* restrict agm, realM apm[MWI>>VWM_SHIFT],
                             const int kSizeM, const int idk, const int kwg) {
  #pragma unroll
  for (int mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {

    #if USE_CL_MAD == 1

      // Computes the indices based on strided/non-strided access
      #if STRM == 0
        int mg = mad24((int) get_local_id(0),(MWI >> VWM_SHIFT),mi);
      #elif STRM == 1
        int mg = mad24(mi,MDIMC,(int) get_local_id(0));
      #endif

      // Computes the indices for the global memory
      int idm = mad24((int) GetGroupID0() , (MWG >> VWM_SHIFT), mg);

      // Loads the data from global memory (not transposed) and stores into registers
      apm[mi] = agm[mad24(idk,(kSizeM >> VWM_SHIFT) , idm)];

    #else

      // Computes the indices based on strided/non-strided access
      #if STRM == 0
        int mg = mi + get_local_id(0)*(MWI >> VWM_SHIFT);
      #elif STRM == 1
        int mg = get_local_id(0) + mi*MDIMC;
      #endif

      // Computes the indices for the global memory
      int idm = mg + GetGroupID0() * (MWG >> VWM_SHIFT);

      // Loads the data from global memory (not transposed) and stores into registers
      apm[mi] = agm[idk*(kSizeM >> VWM_SHIFT) + idm];

    #endif

  }
}
#endif

// Same as above, but now for the B input matrix
#if SB == 0
inline void GlobalToPrivateB(const __global realN* restrict bgm, realN bpm[NWI >> VWN_SHIFT],
                             const int kSizeN, const int idk) {
  #pragma unroll
  for (int ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {

    #if USE_CL_MAD == 1

      // Computes the indices based on strided/non-strided access
      #if STRN == 0
        int ng = mad24((int) get_local_id(1),(NWI >> VWN_SHIFT),ni);
      #elif STRN == 1
        int ng = mad24(ni,NDIMC,(int) get_local_id(1));
      #endif

      // Computes the indices for the global memory
      int idn = mad24((int) GetGroupID1() , (NWG >> VWN_SHIFT), ng);

      // Loads the data from global memory (transposed) and stores into registers
      bpm[ni] = bgm[mad24(idk,(kSizeN >> VWN_SHIFT) , idn)];

    #else

      // Computes the indices based on strided/non-strided access
      #if STRN == 0
        int ng = ni + get_local_id(1)*(NWI >> VWN_SHIFT);
      #elif STRN == 1
        int ng = get_local_id(1) + ni*NDIMC;
      #endif

      // Computes the indices for the global memory
      int idn = ng + GetGroupID1() * (NWG >> VWN_SHIFT);

      // Loads the data from global memory (transposed) and stores into registers
      bpm[ni] = bgm[idk*(kSizeN >> VWN_SHIFT) + idn];


    #endif

  }
}
#endif

// =================================================================================================

// Caches on-chip local memory into per-thread private memory (registers). This function is specific
// for caching the A input matrix.
#if SA == 1
inline void LocalToPrivateA(__local realM* alm, realM apm[MWI >> VWM_SHIFT], const int kg) {

  #pragma unroll
  for (int mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {

    #if USE_CL_MAD == 2		// bad performance of mad24

      #if STRM == 0
        int mg = mad24((int) get_local_id(0),(MWI >> VWM_SHIFT),mi);
      #elif STRM == 1
        int mg = mad24(mi,MDIMC,(int) get_local_id(0));
      #endif
      apm[mi] = alm[mad24(kg,(MWG >> VWM_SHIFT) , mg)];

    #else

      #if STRM == 0
        int mg = mi + get_local_id(0)*(MWI >> VWM_SHIFT);
      #elif STRM == 1
        int mg = get_local_id(0) + mi*MDIMC;
      #endif
      apm[mi] = alm[kg*(MWG >> VWM_SHIFT) + mg];

    #endif


  }
}
#endif

// Same as above, but now for the B input matrix
#if SB == 1
inline void LocalToPrivateB(__local realN* blm, realN bpm[NWI >> VWN_SHIFT], const int kg) {
  #pragma unroll
  for (int ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {
    #if STRN == 0
      int ng = ni + get_local_id(1)*(NWI >> VWN_SHIFT);
    #elif STRN == 1
      int ng = get_local_id(1) + ni*NDIMC;
    #endif
    bpm[ni] = blm[kg*(NWG >> VWN_SHIFT) + ng];
  }
}
#endif

// =================================================================================================

// End of the C++11 raw string literal
)"

// =================================================================================================
