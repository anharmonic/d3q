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

T=297.1
volume=126.0312
reading_from_store_lw=0
readin_from_lw_and_r2q=0
# when reading files produced by d3_tk.x, with option store_lw=.true.
if(reading_from_store_lw);
  f_vel=load("vel.test.5x5x5.out");
  f_lw=load("lw.test.5x5x5_T297.1_s5.out");
  f_freq=load("freq.test.5x5x5.out");

  nat3=(size(f_lw,2))
  nq = size(f_lw,1)

  vel =f_vel; #already in Ry (bohr length/Ry time)
  freq=f_freq*cmm1_to_ry;
  lw  =2*f_lw*cmm1_to_ry;
  clear f_vel f_lw f_freq
# when reading from list of q-points produced by d3_lw.x and d3_r2q.x (with print_vel=.true.)
elseif (readin_from_lw_and_r2q);
  f_vel=load("test_vel.out");
  f_lw=load("test_T297.1_s25.out");

  nat3=(size(f_vel,2)-5)/3
  nq = size(f_vel,1)

  vel =f_vel(:,6:5+3*nat3)*cmm1_to_ry;
  freq=f_lw(:,6:5+nat3)*cmm1_to_ry;
  lw  =2*f_lw(:,6+nat3:5+nat3*2)*cmm1_to_ry;
else
  error "\n*********************\nThis is a matlab/octave script. Open it with a text editor! You need to specify the input format of data, the cell volume and the temperature!! \n*********************\n "
end


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
tk 




