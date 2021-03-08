#!/usr/bin/octave-cli -qf
# boltzmann constant
global kB_ry = 1.3806504E-23 / (4.35974394E-18/2)
# energy conversion
global ry_to_cmm1 = 13.6058 * 8065.5;
global cmm1_to_ry = 1 / ry_to_cmm1
global ry_to_joule =  0.5* 4.35974394e-18
# time in rydberg is funny:
global ry_to_second = 2* 2.418884326505e-17
global ry_to_meter = 5.2917721092e-11
global ry_to_wattmm1km1 = ry_to_joule / (ry_to_second * ry_to_meter)

function n = bose(freq, T)
  global kB_ry;
  beta = 1./(kB_ry*T);
  n = 1./(exp(beta.*freq)-1);
endfunction

for T = [1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 80 90 100 120 140 160 180 200 220 240 260 280 300 320 340 360 380 400 450 500 550 600];

#T=__T #300
volume=951.106863;
# when reading files produced by d3_tk.x, with option store_lw=.true.

  f_vel=load("vel.Def.0005.11x11x11.out");
  f_lw=load( strrep("lw.11x11x11_T__T_s2.out", "__T", mat2str(T)));
  #f_lwiso=load( strrep("lwiso.Natural.11x11x11_T__T_s2.out", "__T", mat2str(T)));
  f_lwiso=load( strrep("lwiso.Def.0005.11x11x11_T__T_s2.out", "__T", mat2str(T)));
  f_lwcas=load("lwcas.250mu.11x11x11.out");
  f_freq=load("freq.Def.0005.11x11x11.out");

  nat3=(size(f_lw,2));
  nq = size(f_lw,1);

  vel =f_vel; #already in Ry (bohr length/Ry time)
  freq=f_freq*cmm1_to_ry;
  lw  =2*(f_lw+f_lwiso+f_lwcas)*cmm1_to_ry;
  #lw  =2*(f_lw+f_lwiso)*cmm1_to_ry;
  clear f_vel f_lw f_freq f_lwiso f_lwcas

n   = bose(freq, T);

pref=1/(volume*T^2*kB_ry*nq) * ry_to_wattmm1km1;
tk = zeros(3,3);
for iq = [1:nq];
  #
  for nu  = [1:nat3]; 
    #
    aux = pref * freq(iq,nu)^2 * n(iq,nu)*(n(iq,nu)+1) / lw(iq,nu);
    #
    for alpha = [1:3];
    for beta  = [1:3];
      tk(alpha,beta) = tk(alpha,beta) + aux* vel(iq,3*(nu-1)+alpha)*vel(iq,3*(nu-1)+beta);
    end
    end
  end 
  #
end 
printf("0 0  %.0f   %.5f   %.5f  %.5f\n", T, tk(1,1), tk(2,2), tk(3,3))

end # Temperatures

#tk_base=[0.918628E+00   -0.353574E-01   0.491818E-02
#         -0.353574E-01  0.889474E+00 -0.244528E-02
#         0.491818E-02 -0.244528E-02  0.261615E+00]
#
#tk ./ tk_base




