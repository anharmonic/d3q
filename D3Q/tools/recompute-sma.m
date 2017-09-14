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

f_vel=load("vel.test.5x5x5.out");
f_lw=load("lw.test.5x5x5_T297.1_s5.out");
f_freq=load("freq.test.5x5x5.out");
# these two have to be set by hand:
T=297.1
volume=126.0312

nat3=(size(f_lw,2))
nat=nat3/3
nq = size(f_vel,1)

vel =f_vel; #already in Ry (bohr length/Ry time)
freq=f_freq*cmm1_to_ry;
lw  =2*f_lw*cmm1_to_ry;
clear f_vel f_lw f_freq

n   = bose(freq, T);


save vel

#pref=1/(T^2 * kB_ry * nq* volume)

pref=1/(volume*T^2*kB_ry*nq) * ry_to_wattmm1km1
tk = zeros(3,3);
for iq = [1:nq];
  #
  for band  = [1:nat3]; 
    #
    aux = pref * freq(iq,band)^2 * n(iq,band)*(n(iq,band)+1) / lw(iq,band);
    #
    for alpha = [1:3];
    for beta  = [1:3];
      tk(alpha,beta) = tk(alpha,beta) + aux* vel(iq,3*(band-1)+alpha)*vel(iq,3*(band-1)+beta);
    end
    end
  end 
  #
end 
tk 




