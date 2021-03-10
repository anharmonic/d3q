&lwinput
  calculation = 'final'
  prefix="final_q0"
  file_mat2 = 'mat2R'
  file_mat3 = 'mat3R.asr.sparse'
  outdir    = './LW/'  ! where to save the data
  asr2 = 'simple'   ! second order sum rule, analogous to "crystal"
  nk =  17, 17, 17     ! size of the integration grid
  xk0 = 1,1,1
  sort_freq = 'overlap'
  ne = 2000
  de = 0.5

   nu_initial = 6
   q_initial = 0.001,    0.,    0.
  e_initial = 682.  ! leave undefined to use the energy of phonon band nu_initial at this q

   q_resolved = .true.
   sigmaq= 0.1

/
  CONFIGS
   1
  5.0 300
  QPOINTS
5
0.0000000000    0.0000000000    0.0000000000 1 Γ
1.0000000000    0.0000000000    0.0000000000 20 X
1.0000000000    1.0000000000    0.0000000000 20 X
0.0000000000    0.0000000000    0.0000000000 50 Γ
0.5000000000    0.5000000000    0.5000000000 20 L
