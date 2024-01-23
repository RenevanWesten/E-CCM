To run the TAMS procedure on the various model settings, use the following commands:

%CCM, on-to-off transitions
python TAMS_Excim.py --N_Ea 101 --N_noise 101 --start "on"

%CCM, off-to-on transitions
python TAMS_Excim.py --N_Ea 101 --N_noise 101 --start "off"

%E-CCM (temperature only), on-to-off transitions
python TAMS_Excim.py --N_Ea 101 --N_noise 101 --start "on" --temp_params 0.0000035
 
%E-CCM (temperature only), off-to-on transitions
python TAMS_Excim.py --N_Ea 101 --N_noise 101 --start "off" --temp_params 0.0000035

%E-CCM, on-to-off transitions
python TAMS_Excim.py --N_Ea 101 --N_noise 101 --start "on" --temp_params 0.0000035 --ice_params 20 0.020 --min_level 1e6

%E-CCM, off-to-on transitions
python TAMS_Excim.py --N_Ea 101 --N_noise 101 --start "off" --temp_params 0.0000035 --ice_params 20 0.020 --min_level 1e6