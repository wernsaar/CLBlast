
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


// Settings
#ifndef USE_VECTOR_MAD
  #define USE_VECTOR_MAD 0      // Unroll (0) or don't (1) unroll the vector MAD manually
#endif
#ifndef GLOBAL_MEM_FENCE
  #define GLOBAL_MEM_FENCE 0    // Global synchronisation barrier for potential better performance
#endif
#ifndef USE_INITIALIZED_ARRAYS
  #define USE_INITIALIZED_ARRAYS 1
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
#elif MDIMA == 128
    #define MDIMA_SHIFT 7
#elif MDIMA == 256
    #define MDIMA_SHIFT 8
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
#elif NDIMB == 128
    #define NDIMB_SHIFT 7
#elif NDIMB == 256
    #define NDIMB_SHIFT 8
#endif

#if MDIMC == 1
    #define MDIMC_SHIFT 0
#elif MDIMC == 2
    #define MDIMC_SHIFT 1
#elif MDIMC == 4
    #define MDIMC_SHIFT 2
#elif MDIMC == 8
    #define MDIMC_SHIFT 3
#elif MDIMC == 16
    #define MDIMC_SHIFT 4
#elif MDIMC == 32
    #define MDIMC_SHIFT 5
#elif MDIMC == 64
    #define MDIMC_SHIFT 6
#elif MDIMC == 128
    #define MDIMC_SHIFT 7
#elif MDIMC == 256
    #define MDIMC_SHIFT 8
#endif

#if NDIMC == 1
    #define NDIMC_SHIFT 0
#elif NDIMC == 2
    #define NDIMC_SHIFT 1
#elif NDIMC == 4
    #define NDIMC_SHIFT 2
#elif NDIMC == 8
    #define NDIMC_SHIFT 3
#elif NDIMC == 16
    #define NDIMC_SHIFT 4
#elif NDIMC == 32
    #define NDIMC_SHIFT 5
#elif NDIMC == 64
    #define NDIMC_SHIFT 6
#elif NDIMC == 128
    #define NDIMC_SHIFT 7
#elif NDIMC == 256
    #define NDIMC_SHIFT 8
#endif


#if NWG == 1
    #define NWG_SHIFT 0
#elif NWG == 2
    #define NWG_SHIFT 1
#elif NWG == 4
    #define NWG_SHIFT 2
#elif NWG == 8
    #define NWG_SHIFT 3
#elif NWG == 16
    #define NWG_SHIFT 4
#elif NWG == 32
    #define NWG_SHIFT 5
#elif NWG == 64
    #define NWG_SHIFT 6
#elif NWG == 128
    #define NWG_SHIFT 7
#elif NWG == 256
    #define NWG_SHIFT 8
#endif

#if MWG == 1
    #define MWG_SHIFT 0
#elif MWG == 2
    #define MWG_SHIFT 1
#elif MWG == 4
    #define MWG_SHIFT 2
#elif MWG == 8
    #define MWG_SHIFT 3
#elif MWG == 16
    #define MWG_SHIFT 4
#elif MWG == 32
    #define MWG_SHIFT 5
#elif MWG == 64
    #define MWG_SHIFT 6
#elif MWG == 128
    #define MWG_SHIFT 7
#elif MWG == 256
    #define MWG_SHIFT 8
#endif

#if KWG == 1
    #define KWG_SHIFT 0
#elif KWG == 2
    #define KWG_SHIFT 1
#elif KWG == 4
    #define KWG_SHIFT 2
#elif KWG == 8
    #define KWG_SHIFT 3
#elif KWG == 16
    #define KWG_SHIFT 4
#elif KWG == 32
    #define KWG_SHIFT 5
#elif KWG == 64
    #define KWG_SHIFT 6
#elif KWG == 128
    #define KWG_SHIFT 7
#elif KWG == 256
    #define KWG_SHIFT 8
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

#define MWI_SHIFT   ((MWG_SHIFT - MDIMC_SHIFT) > 0 ? (MWG_SHIFT - MDIMC_SHIFT) : 0 )
#define NWI_SHIFT   ((NWG_SHIFT - NDIMC_SHIFT) > 0 ? (NWG_SHIFT - NDIMC_SHIFT) : 0 )
#define KDIMA_SHIFT ((MDIMC_SHIFT + NDIMC_SHIFT - MDIMA_SHIFT) > 0 ? (MDIMC_SHIFT + NDIMC_SHIFT - MDIMA_SHIFT) : 0 )
#define KDIMB_SHIFT ((MDIMC_SHIFT + NDIMC_SHIFT - NDIMB_SHIFT) > 0 ? (MDIMC_SHIFT + NDIMC_SHIFT - NDIMB_SHIFT) : 0 )
#define MWA_SHIFT   ((MWG_SHIFT - MDIMA_SHIFT) > 0 ? (MWG_SHIFT - MDIMA_SHIFT) : 0 )
#define KWA_SHIFT   ((KWG_SHIFT - KDIMA_SHIFT) > 0 ? (KWG_SHIFT - KDIMA_SHIFT) : 0 )
#define KWB_SHIFT   ((KWG_SHIFT - KDIMB_SHIFT) > 0 ? (KWG_SHIFT - KDIMB_SHIFT) : 0 )
#define NWB_SHIFT   ((NWG_SHIFT - NDIMB_SHIFT) > 0 ? (NWG_SHIFT - NDIMB_SHIFT) : 0 )

// =================================================================================================

// Initializes the accumulation registers to zero
inline void InitAccRegisters(realM cpm[NWI][MWI >> VWM_SHIFT]) {

  #if (USE_VECTOR_MAD == 1) && ((PRECISION == 32) || (PRECISION == 64))
    const realM z = (realM) ZERO;
  #endif

  // #pragma unroll 
  for (uint mi=0; mi<(MWI>>VWM_SHIFT); ++mi) {
    // #pragma unroll
    for (uint ni=0; ni<NWI; ++ni) {
      #if (USE_VECTOR_MAD == 1) && ((PRECISION == 32) || (PRECISION == 64))
        cpm[ni][mi] = z;
      #else
  
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

  const uint  kSizeMxVWM = kSizeM >> VWM_SHIFT;
  const uint  la1_M1 = (tid >> MDIMA_SHIFT) << KWA_SHIFT;

  #if STRM == 0
    const uint la0_M1 = (tid - (tid & -MDIMA)) << ((MWA_SHIFT - VWM_SHIFT) > 0 ? (MWA_SHIFT - VWM_SHIFT) : 0);
    #if USE_MAD24 == 1
      const uint  kSizeMxVWMM = mad24((uint)(la1_M1+kwg),kSizeMxVWM ,(uint) ((GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) > 0 ? (MWG_SHIFT - VWM_SHIFT) : 0)) + la0_M1));
    #else
      const uint  kSizeMxVWMM = (la1_M1+kwg)*kSizeMxVWM + (GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) > 0 ? (MWG_SHIFT - VWM_SHIFT) : 0)) + la0_M1;
    #endif
  #else
    const uint la0 = tid - (tid & -MDIMA);
    #if USE_MAD24 == 1
      const uint  kSizeMxVWMM = mad24((uint)(la1_M1+kwg),kSizeMxVWM ,(uint) ((GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) > 0 ? (MWG_SHIFT - VWM_SHIFT) : 0)) + la0));
    #else
      const uint  kSizeMxVWMM = (la1_M1+kwg)*kSizeMxVWM + (GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) > 0 ? (MWG_SHIFT - VWM_SHIFT) : 0)) + la0;
    #endif
  #endif

  #if (MWA >> VWM_SHIFT) == 1 

        uint kSizeMxVWMP = 0;
    
        #pragma unroll 
        for (uint kia=0; kia<KWA; ++kia) {
    
            // Loads the data from global memory (not transposed) into the local memory
            #if STRM == 0
              alm[((kia + la1_M1) << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0 )) + la0_M1] = agm[kSizeMxVWMP + kSizeMxVWMM];
            #else
              alm[((kia + la1_M1) << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0 )) + la0] = agm[kSizeMxVWMP + kSizeMxVWMM];
            #endif
            kSizeMxVWMP+=kSizeMxVWM;
    
        }

  #else

      #pragma unroll
      for (uint mia=0; mia<(MWA >> VWM_SHIFT); ++mia) {
    
        #if STRM == 0
          const uint mg = mia + la0_M1;
          const uint idm = mia + kSizeMxVWMM;
        #else
          const uint mg = la0 + (mia << MDIMA_SHIFT);
          const uint idm = (mia << MDIMA_SHIFT) + kSizeMxVWMM;
        #endif
    
        uint kSizeMxVWMP = 0;
    
        #pragma unroll 
        for (uint kia=0; kia<KWA; ++kia) {
    
            // Loads the data from global memory (not transposed) into the local memory
            alm[((kia + la1_M1) << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0 )) + mg] = agm[kSizeMxVWMP + idm];
            kSizeMxVWMP+=kSizeMxVWM;
    
        }
      }

   #endif 

}
#endif

// Same as above, but now for the B input matrix
#if SB == 1
inline void GlobalToLocalB(const __global realN* restrict bgm, __local realN* blm,
                           const int kSizeN, const int tid, const int kwg) {

  const uint lb1_M1 = (tid >> NDIMB_SHIFT) << KWB_SHIFT;
  const uint kSizeNxVWN = kSizeN >> VWN_SHIFT;

  #if STRN == 0
    const uint lb0_M1 = (tid - (tid & -NDIMB)) << ((NWB_SHIFT - VWN_SHIFT) >0 ? (NWB_SHIFT - VWN_SHIFT) : 0);
  #else
    const uint lb0 = tid - (tid & -NDIMB);
  #endif

  #if STRN == 0
    #if USE_MAD24 == 1
      const uint kSizeNxVWNM = mad24((uint) (lb1_M1 + kwg),kSizeNxVWN ,(uint) ((GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0_M1));
    #else
      const uint kSizeNxVWNM = (lb1_M1 + kwg) * kSizeNxVWN + (GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0_M1;
    #endif
  #else
    #if USE_MAD24 == 1
      const uint kSizeNxVWNM = mad24((uint)(lb1_M1 + kwg),kSizeNxVWN,(uint) ((GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0));
    #else
      const uint kSizeNxVWNM = (lb1_M1 + kwg) * kSizeNxVWN + (GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0;
    #endif
  #endif

  uint kSizeNxVWNP = 0;

  #if (NWB >> VWN_SHIFT) == 1

      #pragma unroll
      for (uint kib=0; kib<KWB; ++kib) {
    
        #if STRN == 0
          blm[((kib + lb1_M1) << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0_M1] = bgm[kSizeNxVWNP + kSizeNxVWNM];
        #else
          blm[((kib + lb1_M1) << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0] = bgm[kSizeNxVWNP + kSizeNxVWNM];
        #endif
    
        kSizeNxVWNP+=kSizeNxVWN;
      }

  #else

      #pragma unroll
      for (uint kib=0; kib<KWB; ++kib) {
    
        const uint idk_M1 = kSizeNxVWNP + kSizeNxVWNM;

        #if STRN == 0
          const uint kg_M1 = ((kib + lb1_M1) << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0_M1;
        #else
          const uint kg_M1 = ((kib + lb1_M1) << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + lb0;
        #endif
    
        #pragma unroll
        for (uint nib=0; nib<(NWB >> VWN_SHIFT); ++nib) {
    
            #if STRN == 0
              blm[kg_M1 + nib] = bgm[idk_M1 + nib];
            #else
              blm[kg_M1 + (nib << NDIMB_SHIFT)] = bgm[idk_M1 + (nib << NDIMB_SHIFT)];
            #endif
    
        }
        kSizeNxVWNP+=kSizeNxVWN;
      }

  #endif

}
#endif

// =================================================================================================

// Caches global off-chip memory directly into per-thread private memory (registers). This function
// is specific for caching the A input matrix.
#if SA == 0
inline void GlobalToPrivateA(const __global realM* restrict agm, realM apm[MWI>>VWM_SHIFT],
                             const int kSizeM, const int idk, const int kwg) {

  #if (MWI >> VWM_SHIFT) == 1
 
      apm[0] = agm[idk*(kSizeM >> VWM_SHIFT) + (GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0)) + get_local_id(0)];

  #else

      #if STRM == 0
        uint LocalID0_M1 = get_local_id(0) << ((MWI_SHIFT - VWM_SHIFT) > 0 ? (MWI_SHIFT - VWM_SHIFT) : 0);
      #else
        uint local_id0 = get_local_id(0);
      #endif 
    
      #if USE_MAD24 == 1
        uint idk_M1 = mad24((uint) idk,(uint) (kSizeM >> VWM_SHIFT) , (uint) GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
      #else
        uint idk_M1 = idk*(kSizeM >> VWM_SHIFT) + ((uint) GetGroupID0() << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
      #endif
    
      #pragma unroll
      for (uint mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {
    
          #if STRM == 0
            apm[mi] = agm[idk_M1 + mi + LocalID0_M1];
          #else
            apm[mi] = agm[idk_M1 + (local_id0 + (mi << MDIMC_SHIFT))];
          #endif
    
      }

   #endif

}
#endif

// Same as above, but now for the B input matrix
#if SB == 0
inline void GlobalToPrivateB(const __global realN* restrict bgm, realN bpm[NWI >> VWN_SHIFT],
                             const int kSizeN, const int idk) {

  #if (NWI >> VWN_SHIFT) == 1

    #if USE_MAD24 == 1
      bpm[0] = bgm[mad24((uint) idk,(uint) (kSizeN >> VWN_SHIFT) ,(uint) ((GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)) + get_local_id(1)))];
    #else

      bpm[0] = bgm[idk*(kSizeN >> VWN_SHIFT) + (GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0))
             + get_local_id(1)];

    #endif

  #else

       #if STRN == 0
         uint LocalID1_N1 = get_local_id(1) << ((NWI_SHIFT - VWN_SHIFT) >0 ? (NWI_SHIFT - VWN_SHIFT) : 0);
       #else
         uint local_id1 = get_local_id(1); 
       #endif
     
       #if USE_MAD24 == 1
         uint idk_N1 = mad24((uint) idk,(uint) (kSizeN >> VWN_SHIFT) , ((uint) GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0)));
       #else
         uint idk_N1 = idk*(kSizeN >> VWN_SHIFT) + ((uint) GetGroupID1() << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0));
       #endif
     
       #pragma unroll
       for (uint ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {
     
           #if STRN == 0
             bpm[ni] = bgm[idk_N1 + ni + LocalID1_N1];
           #else
             bpm[ni] = bgm[idk_N1 + (local_id1 + (ni << NDIMC_SHIFT))];
           #endif
     
       }
  #endif
}
#endif

// =================================================================================================

// Caches on-chip local memory into per-thread private memory (registers). This function is specific
// for caching the A input matrix.
#if SA == 1
inline void LocalToPrivateA(__local realM* alm, realM apm[MWI >> VWM_SHIFT], const int kg) {

  #if (MWI >> VWM_SHIFT) == 1

      #if STRM == 0
        apm[0] = alm[(get_local_id(0) << ((MWI_SHIFT - VWM_SHIFT) >0 ? (MWI_SHIFT - VWM_SHIFT) : 0)) + (kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0))];
      #else
        apm[0] = alm[get_local_id(0) + (kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0))];
      #endif

  #else

    #if STRM == 0
      const uint LocalID0_M1 = (get_local_id(0) << ((MWI_SHIFT - VWM_SHIFT) >0 ? (MWI_SHIFT - VWM_SHIFT) : 0)) + (kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
    #else
      const uint local_id0 = get_local_id(0) + (kg << ((MWG_SHIFT - VWM_SHIFT) >0 ? (MWG_SHIFT - VWM_SHIFT) : 0));
    #endif

    #pragma unroll
    for (uint mi=0; mi<(MWI >> VWM_SHIFT); ++mi) {

      #if STRM == 0
        apm[mi] = alm[mi + LocalID0_M1];
      #else
        apm[mi] = alm[local_id0 + (mi << MDIMC_SHIFT)];
      #endif

    }
  #endif


}
#endif

// Same as above, but now for the B input matrix
#if SB == 1
inline void LocalToPrivateB(__local realN* blm, realN bpm[NWI >> VWN_SHIFT], const int kg) {
    
   #if (NWI >> VWN_SHIFT) == 1

      #if STRN == 0
          bpm[0] = blm[(get_local_id(1) << ((NWI_SHIFT - VWN_SHIFT) >0 ? (NWI_SHIFT - VWN_SHIFT) : 0)) +
              (kg << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0))];
      #else
          bpm[0] = blm[get_local_id(1) + 
              (kg << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0))];
      #endif
    
   #else

      #if STRN == 0
        const uint local_id1_M1 = (get_local_id(1) << ((NWI_SHIFT - VWN_SHIFT) >0 ? (NWI_SHIFT - VWN_SHIFT) : 0)) +  
              (kg << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0));
      #else
        const uint local_id1 = get_local_id(1) + 
              (kg << ((NWG_SHIFT - VWN_SHIFT) >0 ? (NWG_SHIFT - VWN_SHIFT) : 0));
      #endif
    
      #pragma unroll
      for (uint ni=0; ni<(NWI >> VWN_SHIFT); ++ni) {
    
        #if STRN == 0
          bpm[ni] = blm[ni + local_id1_M1];
        #else
          bpm[ni] = blm[local_id1 + (ni << NDIMC_SHIFT)];
        #endif
    
      }

   #endif
}
#endif

// =================================================================================================

// End of the C++11 raw string literal
)"

// =================================================================================================
