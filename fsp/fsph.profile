Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls   s/call   s/call  name    
 30.46     25.11    25.11       19     1.32     1.65  .__simulation_box_NMOD_minimum_image
 19.95     41.56    16.45       17     0.97     0.97  .__neighbour_list_NMOD_compress_nlist
 11.28     50.86     9.30       17     0.55     2.19  .__neighbour_list_NMOD_find_pair_separations
  9.08     58.35     7.49       23     0.33     0.33  .__simulation_box_NMOD_rboxv
  5.18     62.62     4.27       16     0.27     0.90  .__sphforce_NMOD_calc_sphforce
  4.60     66.41     3.79       16     0.24     0.24  .__sphforce_NMOD_calc_capillary_pressure
  3.17     69.02     2.61       16     0.16     0.16  .__sphforce_NMOD_calc_grad_v
  3.14     71.61     2.59   961096     0.00     0.00  .__kernel_NMOD_calckernel
  2.28     73.49     1.88       16     0.12     0.12  .__sphforce_NMOD_calc_viscous_entropy_os
  1.49     74.72     1.23       16     0.08     0.08  .__sphforce_NMOD_calc_heat_flux
  1.30     75.79     1.07        1     1.07     2.72  .__g_r_NMOD_accum_g_r
  1.30     76.86     1.07        1     1.07     3.33  .__neighbour_list_NMOD_form_nlist
  1.20     77.85     0.99       17     0.06     0.06  .__density_NMOD_sum_density
  0.93     78.62     0.77       33     0.02     0.10  .__kernel_NMOD_calc_kernels_all
  0.75     79.24     0.62        4     0.16    18.18  .__sphstep_NMOD_sph_step_runge_kutta
  0.74     79.85     0.61        1     0.61     0.61  .__neighbour_list_NMOD_all_pairs_method
  0.63     80.37     0.52        4     0.13     0.13  .__neighbour_list_NMOD_increment_nlist_drtot
  0.55     80.82     0.45       16     0.03     0.03  .__neighbour_list_NMOD_calc_dv
  0.53     81.26     0.44   455624     0.00     0.00  .__core_potential_NMOD_apply_core_repulsion
  0.53     81.70     0.44        4     0.11     0.11  .__particle_NMOD_calc_smoothed_properties
  0.27     81.92     0.22        4     0.06     0.06  .__neighbour_list_NMOD_reform_nlist_now
  0.25     82.13     0.21                             .fsph
  0.17     82.27     0.14        1     0.14     0.14  .__neighbour_list_NMOD_init_nlist
  0.07     82.33     0.06   455624     0.00     0.00  .__art_viscosity_NMOD_calc_art_viscosity
  0.02     82.35     0.02       16     0.00     0.00  .__sphforce_NMOD_calc_grad_v_os
  0.02     82.37     0.02       16     0.00     0.00  .__sphforce_NMOD_calc_pi_os
  0.02     82.39     0.02        4     0.01     0.33  .__simulation_box_NMOD_ac_apply_pbcs
  0.01     82.40     0.01    18000     0.00     0.00  .__eos_NMOD_update_temperature
  0.01     82.41     0.01    15300     0.00     0.00  .__eos_NMOD_vdweos_attractive
  0.01     82.42     0.01     4500     0.00     0.00  .__eos_NMOD_vdwenergy
  0.01     82.43     0.01        1     0.01     0.01  .__particle_NMOD_initialise_particles
  0.01     82.44     0.01                             .__art_viscosity_NMOD__&&_art_viscosity
  0.01     82.45     0.01                             .__core_potential_NMOD__&&_core_potential
  0.00     82.45     0.00    18000     0.00     0.00  .__eos_NMOD_vdwtemp
  0.00     82.45     0.00    15300     0.00     0.00  .__eos_NMOD_vdweos_repulsive
  0.00     82.45     0.00       20     0.00     0.00  .__reader_NMOD_read_dbl
  0.00     82.45     0.00       20     0.00     0.00  .__reader_NMOD_read_int
  0.00     82.45     0.00       16     0.00     0.00  .__sphforce_NMOD_calc_div_v
  0.00     82.45     0.00       16     0.00     0.00  .__sphforce_NMOD_calc_iso_pressure
  0.00     82.45     0.00       16     0.00     0.00  .__sphforce_NMOD_calc_pi_one
  0.00     82.45     0.00        8     0.00     0.00  .__reader_NMOD_read_bool
  0.00     82.45     0.00        8     0.00     0.00  .__system_properties_NMOD_average_scalar_property
  0.00     82.45     0.00        4     0.00     0.00  .__particle_NMOD_check_velocity
  0.00     82.45     0.00        4     0.00     0.00  .__system_properties_NMOD_average_temperature
  0.00     82.45     0.00        4     0.00     0.00  .__system_properties_NMOD_isolated_internal_energy
  0.00     82.45     0.00        4     0.00     0.00  .__system_properties_NMOD_total_internal_energy
  0.00     82.45     0.00        4     0.00     0.00  .__system_properties_NMOD_total_kinetic_energy
  0.00     82.45     0.00        4     0.00     0.00  .__thermostat_NMOD_apply_thermostat
  0.00     82.45     0.00        4     0.00     0.00  .__writer_NMOD_write_properties
  0.00     82.45     0.00        1     0.00     0.00  .__boundary_NMOD_rd_boundary_input
  0.00     82.45     0.00        1     0.00     0.00  .__eos_NMOD_rd_eos_input
  0.00     82.45     0.00        1     0.00     0.00  .__g_r_NMOD_create_g_r
  0.00     82.45     0.00        1     0.00     0.00  .__g_r_NMOD_destroy_g_r
  0.00     82.45     0.00        1     0.00     0.00  .__g_r_NMOD_rdinput_g_r
  0.00     82.45     0.00        1     0.00     0.00  .__g_r_NMOD_write_g_r
  0.00     82.45     0.00        1     0.00     0.00  .__kernel_NMOD_rd_kernel_input
  0.00     82.45     0.00        1     0.00     0.00  .__neighbour_list_NMOD_rd_nlist_input
  0.00     82.45     0.00        1     0.00     0.00  .__particle_NMOD_create_particles
  0.00     82.45     0.00        1     0.00     0.00  .__particle_NMOD_destroy_particles
  0.00     82.45     0.00        1     0.00     0.00  .__particle_NMOD_rd_particle_input
  0.00     82.45     0.00        1     0.00     0.00  .__reader_NMOD_read_input
  0.00     82.45     0.00        1     0.00     0.00  .__simulation_box_NMOD_create_box
  0.00     82.45     0.00        1     0.00     0.00  .__simulation_box_NMOD_destroy_box
  0.00     82.45     0.00        1     0.00     0.00  .__simulation_box_NMOD_init_box
  0.00     82.45     0.00        1     0.00     0.00  .__simulation_box_NMOD_rd_box_input
  0.00     82.45     0.00        1     0.00     0.00  .__simulation_box_NMOD_volume
  0.00     82.45     0.00        1     0.00     0.00  .__sllod_NMOD_rd_sllod_input
  0.00     82.45     0.00        1     0.00     0.00  .__thermostat_NMOD_rd_thermostat_input
  0.00     82.45     0.00        1     0.00     0.00  .__writer_NMOD_clearfile
  0.00     82.45     0.00        1     0.00     0.00  .__writer_NMOD_write_state

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.
 
 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this 
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

		     Call graph (explanation follows)


granularity: each sample hit covers 2 byte(s) for 0.01% of 82.45 seconds

index % time    self  children    called     name
                                   1             .fsph [1]
[1]    100.0    0.21   82.22       0+1       .fsph [1]
                0.62   72.08       4/4           .__sphstep_NMOD_sph_step_runge_kutta [2]
                1.07    2.26       1/1           .__neighbour_list_NMOD_form_nlist [9]
                1.07    1.65       1/1           .__g_r_NMOD_accum_g_r [11]
                0.55    1.65       1/17          .__neighbour_list_NMOD_find_pair_separations [3]
                0.97    0.00       1/17          .__neighbour_list_NMOD_compress_nlist [5]
                0.14    0.00       1/1           .__neighbour_list_NMOD_init_nlist [24]
                0.02    0.08       1/33          .__kernel_NMOD_calc_kernels_all [10]
                0.06    0.00       1/17          .__density_NMOD_sum_density [17]
                0.01    0.00       1/1           .__particle_NMOD_initialise_particles [31]
                0.00    0.00     900/4500        .__eos_NMOD_vdwenergy [30]
                0.00    0.00     900/15300       .__eos_NMOD_vdweos_attractive [29]
                0.00    0.00     900/15300       .__eos_NMOD_vdweos_repulsive [37]
                0.00    0.00       4/4           .__writer_NMOD_write_properties [49]
                0.00    0.00       4/4           .__particle_NMOD_check_velocity [44]
                0.00    0.00       1/1           .__reader_NMOD_read_input [61]
                0.00    0.00       1/1           .__eos_NMOD_rd_eos_input [51]
                0.00    0.00       1/1           .__kernel_NMOD_rd_kernel_input [56]
                0.00    0.00       1/1           .__boundary_NMOD_rd_boundary_input [50]
                0.00    0.00       1/1           .__particle_NMOD_create_particles [58]
                0.00    0.00       1/1           .__particle_NMOD_rd_particle_input [60]
                0.00    0.00       1/1           .__thermostat_NMOD_rd_thermostat_input [68]
                0.00    0.00       1/1           .__simulation_box_NMOD_rd_box_input [65]
                0.00    0.00       1/1           .__simulation_box_NMOD_create_box [62]
                0.00    0.00       1/1           .__neighbour_list_NMOD_rd_nlist_input [57]
                0.00    0.00       1/1           .__sllod_NMOD_rd_sllod_input [67]
                0.00    0.00       1/1           .__g_r_NMOD_rdinput_g_r [54]
                0.00    0.00       1/1           .__simulation_box_NMOD_init_box [64]
                0.00    0.00       1/1           .__g_r_NMOD_create_g_r [52]
                0.00    0.00       1/1           .__writer_NMOD_write_state [70]
                0.00    0.00       1/1           .__g_r_NMOD_write_g_r [55]
                0.00    0.00       1/1           .__g_r_NMOD_destroy_g_r [53]
                0.00    0.00       1/1           .__particle_NMOD_destroy_particles [59]
                0.00    0.00       1/1           .__simulation_box_NMOD_destroy_box [63]
                0.00    0.00       1/1           .__writer_NMOD_clearfile [69]
                                   1             .fsph [1]
-----------------------------------------------
                0.62   72.08       4/4           .fsph [1]
[2]     88.2    0.62   72.08       4         .__sphstep_NMOD_sph_step_runge_kutta [2]
                8.75   26.36      16/17          .__neighbour_list_NMOD_find_pair_separations [3]
               15.48    0.00      16/17          .__neighbour_list_NMOD_compress_nlist [5]
                4.27   10.06      16/16          .__sphforce_NMOD_calc_sphforce [6]
                0.75    2.46      32/33          .__kernel_NMOD_calc_kernels_all [10]
                0.02    1.30       4/4           .__simulation_box_NMOD_ac_apply_pbcs [15]
                0.93    0.04      16/17          .__density_NMOD_sum_density [17]
                0.52    0.00       4/4           .__neighbour_list_NMOD_increment_nlist_drtot [19]
                0.45    0.00      16/16          .__neighbour_list_NMOD_calc_dv [20]
                0.44    0.01       4/4           .__particle_NMOD_calc_smoothed_properties [21]
                0.22    0.00       4/4           .__neighbour_list_NMOD_reform_nlist_now [23]
                0.01    0.00   18000/18000       .__eos_NMOD_update_temperature [28]
                0.00    0.01       4/4           .__thermostat_NMOD_apply_thermostat [35]
-----------------------------------------------
                0.55    1.65       1/17          .fsph [1]
                8.75   26.36      16/17          .__sphstep_NMOD_sph_step_runge_kutta [2]
[3]     45.2    9.30   28.00      17         .__neighbour_list_NMOD_find_pair_separations [3]
               22.47    5.54      17/19          .__simulation_box_NMOD_minimum_image [4]
-----------------------------------------------
                1.32    0.33       1/19          .__neighbour_list_NMOD_form_nlist [9]
                1.32    0.33       1/19          .__g_r_NMOD_accum_g_r [11]
               22.47    5.54      17/19          .__neighbour_list_NMOD_find_pair_separations [3]
[4]     38.0   25.11    6.19      19         .__simulation_box_NMOD_minimum_image [4]
                6.19    0.00      19/23          .__simulation_box_NMOD_rboxv [7]
-----------------------------------------------
                0.97    0.00       1/17          .fsph [1]
               15.48    0.00      16/17          .__sphstep_NMOD_sph_step_runge_kutta [2]
[5]     20.0   16.45    0.00      17         .__neighbour_list_NMOD_compress_nlist [5]
-----------------------------------------------
                4.27   10.06      16/16          .__sphstep_NMOD_sph_step_runge_kutta [2]
[6]     17.4    4.27   10.06      16         .__sphforce_NMOD_calc_sphforce [6]
                3.79    0.00      16/16          .__sphforce_NMOD_calc_capillary_pressure [8]
                2.61    0.00      16/16          .__sphforce_NMOD_calc_grad_v [12]
                1.88    0.00      16/16          .__sphforce_NMOD_calc_viscous_entropy_os [14]
                1.23    0.00      16/16          .__sphforce_NMOD_calc_heat_flux [16]
                0.44    0.00  455624/455624      .__core_potential_NMOD_apply_core_repulsion [22]
                0.06    0.00  455624/455624      .__art_viscosity_NMOD_calc_art_viscosity [25]
                0.02    0.00      16/16          .__sphforce_NMOD_calc_grad_v_os [26]
                0.02    0.00      16/16          .__sphforce_NMOD_calc_pi_os [27]
                0.00    0.01      16/16          .__sphforce_NMOD_calc_iso_pressure [34]
                0.00    0.00      16/16          .__sphforce_NMOD_calc_div_v [40]
                0.00    0.00      16/16          .__sphforce_NMOD_calc_pi_one [41]
-----------------------------------------------
                1.30    0.00       4/23          .__simulation_box_NMOD_ac_apply_pbcs [15]
                6.19    0.00      19/23          .__simulation_box_NMOD_minimum_image [4]
[7]      9.1    7.49    0.00      23         .__simulation_box_NMOD_rboxv [7]
-----------------------------------------------
                3.79    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[8]      4.6    3.79    0.00      16         .__sphforce_NMOD_calc_capillary_pressure [8]
-----------------------------------------------
                1.07    2.26       1/1           .fsph [1]
[9]      4.0    1.07    2.26       1         .__neighbour_list_NMOD_form_nlist [9]
                1.32    0.33       1/19          .__simulation_box_NMOD_minimum_image [4]
                0.61    0.00       1/1           .__neighbour_list_NMOD_all_pairs_method [18]
-----------------------------------------------
                0.02    0.08       1/33          .fsph [1]
                0.75    2.46      32/33          .__sphstep_NMOD_sph_step_runge_kutta [2]
[10]     4.0    0.77    2.54      33         .__kernel_NMOD_calc_kernels_all [10]
                2.54    0.00  942196/961096      .__kernel_NMOD_calckernel [13]
-----------------------------------------------
                1.07    1.65       1/1           .fsph [1]
[11]     3.3    1.07    1.65       1         .__g_r_NMOD_accum_g_r [11]
                1.32    0.33       1/19          .__simulation_box_NMOD_minimum_image [4]
-----------------------------------------------
                2.61    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[12]     3.2    2.61    0.00      16         .__sphforce_NMOD_calc_grad_v [12]
-----------------------------------------------
                0.01    0.00    3600/961096      .__particle_NMOD_calc_smoothed_properties [21]
                0.04    0.00   15300/961096      .__density_NMOD_sum_density [17]
                2.54    0.00  942196/961096      .__kernel_NMOD_calc_kernels_all [10]
[13]     3.1    2.59    0.00  961096         .__kernel_NMOD_calckernel [13]
-----------------------------------------------
                1.88    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[14]     2.3    1.88    0.00      16         .__sphforce_NMOD_calc_viscous_entropy_os [14]
-----------------------------------------------
                0.02    1.30       4/4           .__sphstep_NMOD_sph_step_runge_kutta [2]
[15]     1.6    0.02    1.30       4         .__simulation_box_NMOD_ac_apply_pbcs [15]
                1.30    0.00       4/23          .__simulation_box_NMOD_rboxv [7]
-----------------------------------------------
                1.23    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[16]     1.5    1.23    0.00      16         .__sphforce_NMOD_calc_heat_flux [16]
-----------------------------------------------
                0.06    0.00       1/17          .fsph [1]
                0.93    0.04      16/17          .__sphstep_NMOD_sph_step_runge_kutta [2]
[17]     1.3    0.99    0.04      17         .__density_NMOD_sum_density [17]
                0.04    0.00   15300/961096      .__kernel_NMOD_calckernel [13]
-----------------------------------------------
                0.61    0.00       1/1           .__neighbour_list_NMOD_form_nlist [9]
[18]     0.7    0.61    0.00       1         .__neighbour_list_NMOD_all_pairs_method [18]
-----------------------------------------------
                0.52    0.00       4/4           .__sphstep_NMOD_sph_step_runge_kutta [2]
[19]     0.6    0.52    0.00       4         .__neighbour_list_NMOD_increment_nlist_drtot [19]
-----------------------------------------------
                0.45    0.00      16/16          .__sphstep_NMOD_sph_step_runge_kutta [2]
[20]     0.5    0.45    0.00      16         .__neighbour_list_NMOD_calc_dv [20]
-----------------------------------------------
                0.44    0.01       4/4           .__sphstep_NMOD_sph_step_runge_kutta [2]
[21]     0.5    0.44    0.01       4         .__particle_NMOD_calc_smoothed_properties [21]
                0.01    0.00    3600/961096      .__kernel_NMOD_calckernel [13]
-----------------------------------------------
                0.44    0.00  455624/455624      .__sphforce_NMOD_calc_sphforce [6]
[22]     0.5    0.44    0.00  455624         .__core_potential_NMOD_apply_core_repulsion [22]
-----------------------------------------------
                0.22    0.00       4/4           .__sphstep_NMOD_sph_step_runge_kutta [2]
[23]     0.3    0.22    0.00       4         .__neighbour_list_NMOD_reform_nlist_now [23]
-----------------------------------------------
                0.14    0.00       1/1           .fsph [1]
[24]     0.2    0.14    0.00       1         .__neighbour_list_NMOD_init_nlist [24]
-----------------------------------------------
                0.06    0.00  455624/455624      .__sphforce_NMOD_calc_sphforce [6]
[25]     0.1    0.06    0.00  455624         .__art_viscosity_NMOD_calc_art_viscosity [25]
-----------------------------------------------
                0.02    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[26]     0.0    0.02    0.00      16         .__sphforce_NMOD_calc_grad_v_os [26]
-----------------------------------------------
                0.02    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[27]     0.0    0.02    0.00      16         .__sphforce_NMOD_calc_pi_os [27]
-----------------------------------------------
                0.01    0.00   18000/18000       .__sphstep_NMOD_sph_step_runge_kutta [2]
[28]     0.0    0.01    0.00   18000         .__eos_NMOD_update_temperature [28]
                0.00    0.00   18000/18000       .__eos_NMOD_vdwtemp [36]
-----------------------------------------------
                0.00    0.00     900/15300       .fsph [1]
                0.01    0.00   14400/15300       .__sphforce_NMOD_calc_iso_pressure [34]
[29]     0.0    0.01    0.00   15300         .__eos_NMOD_vdweos_attractive [29]
-----------------------------------------------
                0.00    0.00     900/4500        .fsph [1]
                0.01    0.00    3600/4500        .__thermostat_NMOD_apply_thermostat [35]
[30]     0.0    0.01    0.00    4500         .__eos_NMOD_vdwenergy [30]
-----------------------------------------------
                0.01    0.00       1/1           .fsph [1]
[31]     0.0    0.01    0.00       1         .__particle_NMOD_initialise_particles [31]
-----------------------------------------------
                                                 <spontaneous>
[32]     0.0    0.01    0.00                 .__art_viscosity_NMOD__&&_art_viscosity [32]
-----------------------------------------------
                                                 <spontaneous>
[33]     0.0    0.01    0.00                 .__core_potential_NMOD__&&_core_potential [33]
-----------------------------------------------
                0.00    0.01      16/16          .__sphforce_NMOD_calc_sphforce [6]
[34]     0.0    0.00    0.01      16         .__sphforce_NMOD_calc_iso_pressure [34]
                0.01    0.00   14400/15300       .__eos_NMOD_vdweos_attractive [29]
                0.00    0.00   14400/15300       .__eos_NMOD_vdweos_repulsive [37]
-----------------------------------------------
                0.00    0.01       4/4           .__sphstep_NMOD_sph_step_runge_kutta [2]
[35]     0.0    0.00    0.01       4         .__thermostat_NMOD_apply_thermostat [35]
                0.01    0.00    3600/4500        .__eos_NMOD_vdwenergy [30]
                0.00    0.00       4/4           .__system_properties_NMOD_average_temperature [45]
-----------------------------------------------
                0.00    0.00   18000/18000       .__eos_NMOD_update_temperature [28]
[36]     0.0    0.00    0.00   18000         .__eos_NMOD_vdwtemp [36]
-----------------------------------------------
                0.00    0.00     900/15300       .fsph [1]
                0.00    0.00   14400/15300       .__sphforce_NMOD_calc_iso_pressure [34]
[37]     0.0    0.00    0.00   15300         .__eos_NMOD_vdweos_repulsive [37]
-----------------------------------------------
                0.00    0.00       1/20          .__thermostat_NMOD_rd_thermostat_input [68]
                0.00    0.00       2/20          .__kernel_NMOD_rd_kernel_input [56]
                0.00    0.00       3/20          .__particle_NMOD_rd_particle_input [60]
                0.00    0.00       4/20          .__eos_NMOD_rd_eos_input [51]
                0.00    0.00      10/20          .__reader_NMOD_read_input [61]
[38]     0.0    0.00    0.00      20         .__reader_NMOD_read_dbl [38]
-----------------------------------------------
                0.00    0.00       1/20          .__eos_NMOD_rd_eos_input [51]
                0.00    0.00       1/20          .__thermostat_NMOD_rd_thermostat_input [68]
                0.00    0.00       2/20          .__simulation_box_NMOD_rd_box_input [65]
                0.00    0.00       2/20          .__kernel_NMOD_rd_kernel_input [56]
                0.00    0.00      14/20          .__reader_NMOD_read_input [61]
[39]     0.0    0.00    0.00      20         .__reader_NMOD_read_int [39]
-----------------------------------------------
                0.00    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[40]     0.0    0.00    0.00      16         .__sphforce_NMOD_calc_div_v [40]
-----------------------------------------------
                0.00    0.00      16/16          .__sphforce_NMOD_calc_sphforce [6]
[41]     0.0    0.00    0.00      16         .__sphforce_NMOD_calc_pi_one [41]
-----------------------------------------------
                0.00    0.00       4/8           .__reader_NMOD_read_input [61]
                0.00    0.00       4/8           .__boundary_NMOD_rd_boundary_input [50]
[42]     0.0    0.00    0.00       8         .__reader_NMOD_read_bool [42]
-----------------------------------------------
                0.00    0.00       8/8           .__writer_NMOD_write_properties [49]
[43]     0.0    0.00    0.00       8         .__system_properties_NMOD_average_scalar_property [43]
-----------------------------------------------
                0.00    0.00       4/4           .fsph [1]
[44]     0.0    0.00    0.00       4         .__particle_NMOD_check_velocity [44]
-----------------------------------------------
                0.00    0.00       4/4           .__thermostat_NMOD_apply_thermostat [35]
[45]     0.0    0.00    0.00       4         .__system_properties_NMOD_average_temperature [45]
-----------------------------------------------
                0.00    0.00       4/4           .__writer_NMOD_write_properties [49]
[46]     0.0    0.00    0.00       4         .__system_properties_NMOD_isolated_internal_energy [46]
-----------------------------------------------
                0.00    0.00       4/4           .__writer_NMOD_write_properties [49]
[47]     0.0    0.00    0.00       4         .__system_properties_NMOD_total_internal_energy [47]
-----------------------------------------------
                0.00    0.00       4/4           .__writer_NMOD_write_properties [49]
[48]     0.0    0.00    0.00       4         .__system_properties_NMOD_total_kinetic_energy [48]
-----------------------------------------------
                0.00    0.00       4/4           .fsph [1]
[49]     0.0    0.00    0.00       4         .__writer_NMOD_write_properties [49]
                0.00    0.00       8/8           .__system_properties_NMOD_average_scalar_property [43]
                0.00    0.00       4/4           .__system_properties_NMOD_total_kinetic_energy [48]
                0.00    0.00       4/4           .__system_properties_NMOD_total_internal_energy [47]
                0.00    0.00       4/4           .__system_properties_NMOD_isolated_internal_energy [46]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[50]     0.0    0.00    0.00       1         .__boundary_NMOD_rd_boundary_input [50]
                0.00    0.00       4/8           .__reader_NMOD_read_bool [42]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[51]     0.0    0.00    0.00       1         .__eos_NMOD_rd_eos_input [51]
                0.00    0.00       4/20          .__reader_NMOD_read_dbl [38]
                0.00    0.00       1/20          .__reader_NMOD_read_int [39]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[52]     0.0    0.00    0.00       1         .__g_r_NMOD_create_g_r [52]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[53]     0.0    0.00    0.00       1         .__g_r_NMOD_destroy_g_r [53]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[54]     0.0    0.00    0.00       1         .__g_r_NMOD_rdinput_g_r [54]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[55]     0.0    0.00    0.00       1         .__g_r_NMOD_write_g_r [55]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[56]     0.0    0.00    0.00       1         .__kernel_NMOD_rd_kernel_input [56]
                0.00    0.00       2/20          .__reader_NMOD_read_int [39]
                0.00    0.00       2/20          .__reader_NMOD_read_dbl [38]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[57]     0.0    0.00    0.00       1         .__neighbour_list_NMOD_rd_nlist_input [57]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[58]     0.0    0.00    0.00       1         .__particle_NMOD_create_particles [58]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[59]     0.0    0.00    0.00       1         .__particle_NMOD_destroy_particles [59]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[60]     0.0    0.00    0.00       1         .__particle_NMOD_rd_particle_input [60]
                0.00    0.00       3/20          .__reader_NMOD_read_dbl [38]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[61]     0.0    0.00    0.00       1         .__reader_NMOD_read_input [61]
                0.00    0.00      14/20          .__reader_NMOD_read_int [39]
                0.00    0.00      10/20          .__reader_NMOD_read_dbl [38]
                0.00    0.00       4/8           .__reader_NMOD_read_bool [42]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[62]     0.0    0.00    0.00       1         .__simulation_box_NMOD_create_box [62]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[63]     0.0    0.00    0.00       1         .__simulation_box_NMOD_destroy_box [63]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[64]     0.0    0.00    0.00       1         .__simulation_box_NMOD_init_box [64]
                0.00    0.00       1/1           .__simulation_box_NMOD_volume [66]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[65]     0.0    0.00    0.00       1         .__simulation_box_NMOD_rd_box_input [65]
                0.00    0.00       2/20          .__reader_NMOD_read_int [39]
-----------------------------------------------
                0.00    0.00       1/1           .__simulation_box_NMOD_init_box [64]
[66]     0.0    0.00    0.00       1         .__simulation_box_NMOD_volume [66]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[67]     0.0    0.00    0.00       1         .__sllod_NMOD_rd_sllod_input [67]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[68]     0.0    0.00    0.00       1         .__thermostat_NMOD_rd_thermostat_input [68]
                0.00    0.00       1/20          .__reader_NMOD_read_int [39]
                0.00    0.00       1/20          .__reader_NMOD_read_dbl [38]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[69]     0.0    0.00    0.00       1         .__writer_NMOD_clearfile [69]
-----------------------------------------------
                0.00    0.00       1/1           .fsph [1]
[70]     0.0    0.00    0.00       1         .__writer_NMOD_write_state [70]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.


Index by function name

  [32] .__art_viscosity_NMOD__&&_art_viscosity [9] .__neighbour_list_NMOD_form_nlist [8] .__sphforce_NMOD_calc_capillary_pressure
  [25] .__art_viscosity_NMOD_calc_art_viscosity [19] .__neighbour_list_NMOD_increment_nlist_drtot [40] .__sphforce_NMOD_calc_div_v
  [50] .__boundary_NMOD_rd_boundary_input [24] .__neighbour_list_NMOD_init_nlist [12] .__sphforce_NMOD_calc_grad_v
  [33] .__core_potential_NMOD__&&_core_potential [57] .__neighbour_list_NMOD_rd_nlist_input [26] .__sphforce_NMOD_calc_grad_v_os
  [22] .__core_potential_NMOD_apply_core_repulsion [23] .__neighbour_list_NMOD_reform_nlist_now [16] .__sphforce_NMOD_calc_heat_flux
  [17] .__density_NMOD_sum_density [21] .__particle_NMOD_calc_smoothed_properties [34] .__sphforce_NMOD_calc_iso_pressure
  [51] .__eos_NMOD_rd_eos_input [44] .__particle_NMOD_check_velocity [41] .__sphforce_NMOD_calc_pi_one
  [28] .__eos_NMOD_update_temperature [58] .__particle_NMOD_create_particles [27] .__sphforce_NMOD_calc_pi_os
  [30] .__eos_NMOD_vdwenergy  [59] .__particle_NMOD_destroy_particles [6] .__sphforce_NMOD_calc_sphforce
  [29] .__eos_NMOD_vdweos_attractive [31] .__particle_NMOD_initialise_particles [14] .__sphforce_NMOD_calc_viscous_entropy_os
  [37] .__eos_NMOD_vdweos_repulsive [60] .__particle_NMOD_rd_particle_input [2] .__sphstep_NMOD_sph_step_runge_kutta
  [36] .__eos_NMOD_vdwtemp    [42] .__reader_NMOD_read_bool [43] .__system_properties_NMOD_average_scalar_property
  [11] .__g_r_NMOD_accum_g_r  [38] .__reader_NMOD_read_dbl [45] .__system_properties_NMOD_average_temperature
  [52] .__g_r_NMOD_create_g_r [61] .__reader_NMOD_read_input [46] .__system_properties_NMOD_isolated_internal_energy
  [53] .__g_r_NMOD_destroy_g_r [39] .__reader_NMOD_read_int [47] .__system_properties_NMOD_total_internal_energy
  [54] .__g_r_NMOD_rdinput_g_r [15] .__simulation_box_NMOD_ac_apply_pbcs [48] .__system_properties_NMOD_total_kinetic_energy
  [55] .__g_r_NMOD_write_g_r  [62] .__simulation_box_NMOD_create_box [35] .__thermostat_NMOD_apply_thermostat
  [10] .__kernel_NMOD_calc_kernels_all [63] .__simulation_box_NMOD_destroy_box [68] .__thermostat_NMOD_rd_thermostat_input
  [13] .__kernel_NMOD_calckernel [64] .__simulation_box_NMOD_init_box [69] .__writer_NMOD_clearfile
  [56] .__kernel_NMOD_rd_kernel_input [4] .__simulation_box_NMOD_minimum_image [49] .__writer_NMOD_write_properties
  [18] .__neighbour_list_NMOD_all_pairs_method [7] .__simulation_box_NMOD_rboxv [70] .__writer_NMOD_write_state
  [20] .__neighbour_list_NMOD_calc_dv [65] .__simulation_box_NMOD_rd_box_input [1] .fsph
   [5] .__neighbour_list_NMOD_compress_nlist [66] .__simulation_box_NMOD_volume
   [3] .__neighbour_list_NMOD_find_pair_separations [67] .__sllod_NMOD_rd_sllod_input
