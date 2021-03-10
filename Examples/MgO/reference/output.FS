  Using      4 MPI processes

      
     For 2n+1 calculations of force constants please cite:
      1. Lorenzo Paulatto, Francesco Mauri, and Michele Lazzeri
         Phys. Rev. B 87, 214303 (2013)
      
     For thermal transport calculations please cite:
      2. Giorgia Fugallo, Michele Lazzeri, Lorenzo Paulatto, and Francesco Mauri
         Phys. Rev. B 88, 045430 (2013)
      3. A. Cepellotti, G. Fugallo, L. Paulatto, M. Lazzeri, F. Mauri, N. Marzari,
         Nature communications 6 (2015)
      4. G. Fugallo, A. Cepellotti, L. Paulatto, M. Lazzeri, N. Marzari, F. Mauri,
         Nano letters 14 (11), 6109-6114 (2014)
      
     For spectral function calculations also cite:
      5. Lorenzo Paulatto, Ion Errea, Matteo Calandra, and Francesco Mauri,
         Phys. Rev. B 91, 054304 (2015)

  Reading input file 'input.FS'
  merging with command line arguments
&LWINPUT
 CALCULATION="final           ",
 OUTDIR="./LW/                                                                                                                                                                                                                                                           ",
 PREFIX="final_q0                                                                                                                                                                                                                                                        ",
 FILE_MAT2="mat2R                                                                                                                                                                                                                                                           ",
 FILE_MAT3="mat3R.asr.sparse                                                                                                                                                                                                                                                ",
 ASR2="simple  ",
 NCONF=-1         ,
 SKIP_Q=0          ,
 NQ=-1         ,
 NK= 3*17         ,
 GRID_TYPE="simple",
 XK0= 3*1.0000000000000000       ,
 OPTIMIZE_GRID=F,
 OPTIMIZE_GRID_THR=  1.0000000000000000E-002,
 NE=2000       ,
 DE= 0.50000000000000000     ,
 E0=  0.0000000000000000     ,
 SIGMA_E= -1.0000000000000000     ,
 NU_INITIAL=6          ,
 E_INITIAL=  682.00000000000000     ,
 Q_INITIAL=  1.0000000000000000E-003, 2*0.0000000000000000       ,
 Q_RESOLVED=T,
 Q_SUMMED=F,
 SIGMAQ= 0.10000000000000001     ,
 EXP_T_FACTOR=F,
 SORT_FREQ="overlap",
 ISOTOPIC_DISORDER=F,
 CASIMIR_SCATTERING=F,
 SAMPLE_LENGTH_AU= -1.0000000000000000     ,
 SAMPLE_LENGTH_MU= -1.0000000000000000     ,
 SAMPLE_LENGTH_MM= -1.0000000000000000     ,
 SAMPLE_DIR= 3*0.0000000000000000       ,
 MAX_SECONDS=-1         ,
 MAX_TIME= -1.0000000000000000     ,
 /
 ** Reading FC2 file: mat2R
     alat (bohr)    7.984849
     tpiba (1/bohr)    0.786888
     at_1 = (   -0.500000    0.000000    0.500000 )
     at_2 = (    0.000000    0.500000    0.500000 )
     at_3 = (   -0.500000    0.500000    0.000000 )
     bg_1 = (   -1.000000   -1.000000    1.000000 )
     bg_2 = (    1.000000    1.000000    1.000000 )
     bg_3 = (   -1.000000    1.000000   -1.000000 )
     volume (bohr^3)  127.274119
 Original FC2 grid:           4           4           4
 ** Reading sparse FC3 file mat3R.asr.sparse
     alat (bohr)    7.984849
     tpiba (1/bohr)    0.786888
     at_1 = (   -0.500000    0.000000    0.500000 )
     at_2 = (    0.000000    0.500000    0.500000 )
     at_3 = (   -0.500000    0.500000    0.000000 )
     bg_1 = (   -1.000000   -1.000000    1.000000 )
     bg_2 = (    1.000000    1.000000    1.000000 )
     bg_3 = (   -1.000000    1.000000   -1.000000 )
     volume (bohr^3)  127.274117
    Original FC3 grid:           2           2           2
    Number of R:               523
 WARNING! FC2 and FC3 systems DO NOT MATCH !!!
 Reading : done.
 Memory used :            0 Mb
 calculation: final           full            
 Reading CONFIGS           1
  Temperatures:
    300.0
  Smearings:   
    5.000

 Reading QPOINTS
  Path turning point:     1.000000
  Path turning point:     2.000000
  Path turning point:     3.414214
  Read   5 lines, set-up    111 q-points, cartesian coordinates

 --> Setting up inner grid
  Applying grid shift   -0.029412    0.029412    0.029412
  Setup a simple grid of     4913 q-points
 opening ./LW//final_q0_T300_s5.out
  Going to compute final state decomposition for    0.001000    0.000000    0.000000  mode:   6  energy:  682.000000cm^-1
         0.0% | STEP TIME:         0.0s | STEP END:          -.-s | WALL:         0.6s 
         0.1% | STEP TIME:         0.0s | STEP END:          3.7s | WALL:         0.6s 
        10.0% | STEP TIME:         0.4s | STEP END:          3.9s | WALL:         1.0s 
        20.0% | STEP TIME:         0.8s | STEP END:          4.1s | WALL:         1.4s 
        30.0% | STEP TIME:         1.2s | STEP END:          4.0s | WALL:         1.8s 
        40.1% | STEP TIME:         1.6s | STEP END:          4.0s | WALL:         2.2s 
        50.0% | STEP TIME:         2.0s | STEP END:          4.0s | WALL:         2.6s 
        60.0% | STEP TIME:         2.4s | STEP END:          4.1s | WALL:         3.0s 
        70.0% | STEP TIME:         2.8s | STEP END:          4.0s | WALL:         3.4s 
        80.0% | STEP TIME:         3.2s | STEP END:          4.0s | WALL:         3.8s 
        90.1% | STEP TIME:         3.6s | STEP END:          4.0s | WALL:         4.2s 
       100.0% | STEP TIME:         4.0s | STEP END:          4.0s | WALL:         4.6s 
   * WALL :       6.9307 s
   *                          *              ms *         ms/call *          ms*cpu *     ms*cpu/call *  % wtime *        calls *
   *        spectral function *     3902.031207 *        0.794226 *       15608.125 *        3.176903 *   56.296 *        4913. *
   *    q-resolved i/o & comm *      652.296033 *       81.537004 *        2609.184 *      326.148017 *    9.411 *           8. *
   *             read fc data *       64.450785 *       16.112696 *         257.803 *       64.450785 *    0.930 *           4. *
   *                read data *       67.112056 *       16.778014 *         268.448 *       67.112056 *    0.968 *           4. *
*** * Contributions to spectral function time:
   *           q-resolved spf *     2622.427510 *        0.533773 *       10489.710 *        2.135093 *   37.835 *        4913. *
   *             bose distrib *        2.366378 *        0.000482 *           9.466 *        0.001927 *    0.034 *        4913. *
   *                sum modes *      466.279971 *        0.094907 *        1865.120 *        0.379630 *    6.727 *        4913. *
   *          fc3 interpolate *      340.217440 *        0.069248 *        1360.870 *        0.276994 *    4.908 *        4913. *
   *           fc3 modulus sq *        1.712756 *        0.000349 *           6.851 *        0.001394 *    0.025 *        4913. *
   *               fc3 rotate *      216.919575 *        0.044152 *         867.678 *        0.176609 *    3.130 *        4913. *

      
     For 2n+1 calculations of force constants please cite:
      1. Lorenzo Paulatto, Francesco Mauri, and Michele Lazzeri
         Phys. Rev. B 87, 214303 (2013)
      
     For thermal transport calculations please cite:
      2. Giorgia Fugallo, Michele Lazzeri, Lorenzo Paulatto, and Francesco Mauri
         Phys. Rev. B 88, 045430 (2013)
      3. A. Cepellotti, G. Fugallo, L. Paulatto, M. Lazzeri, F. Mauri, N. Marzari,
         Nature communications 6 (2015)
      4. G. Fugallo, A. Cepellotti, L. Paulatto, M. Lazzeri, N. Marzari, F. Mauri,
         Nano letters 14 (11), 6109-6114 (2014)
      
     For spectral function calculations also cite:
      5. Lorenzo Paulatto, Ion Errea, Matteo Calandra, and Francesco Mauri,
         Phys. Rev. B 91, 054304 (2015)

