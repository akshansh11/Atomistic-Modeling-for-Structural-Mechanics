# CREATED BY Akshansh Mishra on 16/11/2025 at 10.48 pm Central European Time
# Modified for laser sintering of two aluminum particles - OPTIMIZED I/O
# This work is licensed under Creative Commons Attribution 4.0 International 
# LAMMPS script for laser sintering/welding of two Al particles on substrate

# Initialization
units           metal
dimension       3
atom_style      atomic
boundary        p p p

# Define lattice for Aluminum
lattice         fcc 4.05

# Create simulation box
region          box block -40 40 -40 40 0 50
create_box      3 box

# Create substrate/surface region FIRST
region          substrate block -40 40 -40 40 0 8
create_atoms    3 region substrate

# Create two spherical particles in contact (slight overlap for neck formation)
# Particle 1 - Left sphere (radius 12)
region          particle1 sphere -10 0 20 12
create_atoms    1 region particle1

# Particle 2 - Right sphere (radius 12, touching particle 1)
region          particle2 sphere 10 0 20 12
create_atoms    2 region particle2

# Define atomic masses (all aluminum)
mass            1 26.9815
mass            2 26.9815
mass            3 26.9815

# Define interatomic potential for Aluminum
pair_style      eam/alloy
pair_coeff      * * Al99.eam.alloy Al Al Al

# Delete overlapping atoms at particle-particle and particle-substrate interfaces
delete_atoms    overlap 0.5 all all

# Neighbor list settings
neighbor        2.0 bin
neigh_modify    delay 0 every 1 check yes

# Minimize initial structure
minimize        1.0e-4 1.0e-6 1000 10000

# Define groups
group           particle1 type 1
group           particle2 type 2
group           particles union particle1 particle2
group           substrate type 3

# KEEP SUBSTRATE COMPLETELY FROZEN (acts as heat sink/build platform)
fix             freeze substrate setforce 0.0 0.0 0.0
velocity        substrate set 0.0 0.0 0.0

# Define laser irradiation zones
# Laser focused on the neck region between particles
region          neck_zone block -15 15 -8 8 14 26
region          particle1_region sphere -10 0 20 12 side in
region          particle2_region sphere 10 0 20 12 side in

# Intersection: laser hits both particles in neck area
region          laser_p1 intersect 2 neck_zone particle1_region
region          laser_p2 intersect 2 neck_zone particle2_region

group           laser_zone1 region laser_p1
group           laser_zone2 region laser_p2
group           laser_zones union laser_zone1 laser_zone2

# Heat affected zone (HAZ) - neck region and surrounding area
region          haz_region block -18 18 -10 10 12 28
group           haz region haz_region

# Interface atoms (for neck formation analysis)
region          interface_region block -5 5 -6 6 16 24
group           interface_atoms region interface_region

# Computes for analysis
compute         p1_temp particle1 temp
compute         p2_temp particle2 temp
compute         laser_temp laser_zones temp
compute         interface_temp interface_atoms temp

compute         ke_per_atom all ke/atom
compute         pe_per_atom all pe/atom
compute         coord all coord/atom cutoff 3.5
compute         csym all centro/atom fcc

# Temperature per atom for visualization
variable        temp_atom atom c_ke_per_atom/(1.5*8.617333e-5)

# Stress analysis
compute         stress_tensor all stress/atom NULL
variable        s11 atom c_stress_tensor[1]
variable        s22 atom c_stress_tensor[2]
variable        s33 atom c_stress_tensor[3]
variable        von_mises atom sqrt(0.5*((v_s11-v_s22)^2+(v_s22-v_s33)^2+(v_s33-v_s11)^2+6*(c_stress_tensor[4]^2+c_stress_tensor[5]^2+c_stress_tensor[6]^2)))

# Thermodynamic output
thermo          100
thermo_style    custom step temp c_p1_temp c_p2_temp c_laser_temp c_interface_temp ke pe etotal press atoms vol

# ============================================
# OPTIMIZED DUMPS - REDUCED FREQUENCY
# ============================================

# Main dump - REDUCED to every 1000 steps (was 200)
dump            1 all custom 1000 laser_sintering.lammpstrj id type x y z &
                v_temp_atom c_coord c_csym
dump_modify     1 sort id element Al Al Al

# Neck formation tracking - REDUCED to every 2000 steps (was 100)
dump            neck all custom 2000 neck_formation.lammpstrj id type x y z &
                v_temp_atom c_coord c_csym
dump_modify     neck element Al Al Al

# Temperature field - REDUCED to every 2000 steps (was 200)
dump            temp all custom 2000 temperature_field.lammpstrj id type x y z v_temp_atom
dump_modify     temp element Al Al Al

# Stress field - REDUCED to every 5000 steps (was 500)
dump            stress all custom 5000 stress_field.lammpstrj id type x y z &
                v_von_mises c_csym
dump_modify     stress element Al Al Al

# Interface tracking - REDUCED to every 1000 steps (was 100)
dump            interface interface_atoms custom 1000 interface_atoms.lammpstrj id type x y z &
                v_temp_atom c_coord
dump_modify     interface element Al Al Al

# Restart files - REDUCED to every 50000 steps (was 10000)
restart         50000 restart.laser_sinter

# Smaller timestep for stability
timestep        0.0005

# Initial equilibration - particles at room temperature
print "Initial equilibration at 300K..."
velocity        particles create 300.0 12345 dist gaussian
fix             equil particles nvt temp 300.0 300.0 0.1
run             15000
unfix           equil

# LASER SINTERING PROCESS
print "=========================================="
print "LASER SINTERING - Neck Formation Between Particles"
print "=========================================="

# Phase 1: Laser heating - Focus on neck region
print "Phase 1: Laser heating neck region (0-40 ps)..."

# Apply laser heat source to neck area (both particles)
fix             laser_heat1 laser_zone1 heat 1 150.0 region laser_p1
fix             laser_heat2 laser_zone2 heat 1 150.0 region laser_p2

# Heat neck region to near melting point
fix             laser_temp1 laser_zone1 langevin 300.0 1800.0 0.05 48279
fix             laser_temp2 laser_zone2 langevin 300.0 1800.0 0.05 48279

# Rest of particles follow heat conduction
fix             particles_dynamics particles nve

# Temperature control for non-laser regions
group           particles_no_laser subtract particles laser_zones
fix             particles_temp particles_no_laser temp/rescale 10 300.0 900.0 0.02 1.0

run             80000

print "Phase 1 completed - Neck region heated"

# Phase 2: Hold at sintering temperature - neck formation
print "Phase 2: Steady-state sintering (40-70 ps)..."
unfix           laser_temp1
unfix           laser_temp2
fix             laser_temp1 laser_zone1 langevin 1800.0 1800.0 0.05 48279
fix             laser_temp2 laser_zone2 langevin 1800.0 1800.0 0.05 48279
run             60000

print "Phase 2 completed - Neck formed and growing"

# Phase 3: Extended sintering - particle bonding
print "Phase 3: Extended sintering for stronger bonding..."
run             40000

print "Phase 3 completed - Strong neck bond established"

# Phase 4: Laser shutdown and cooling
print "Phase 4: Laser off - Cooling begins..."
unfix           laser_heat1
unfix           laser_heat2
unfix           laser_temp1
unfix           laser_temp2
unfix           particles_temp

# Natural cooling through conduction
fix             cooling particles langevin 1800.0 1000.0 0.1 48279
run             40000

print "Phase 5: Continued cooling..."
unfix           cooling
fix             cooling particles langevin 1000.0 600.0 0.1 48279
run             40000

print "Phase 6: Final cooling to room temperature..."
unfix           cooling
fix             cooling particles langevin 600.0 300.0 0.1 48279
run             50000

print "Phase 7: Final equilibration..."
unfix           cooling
fix             final_equil particles nvt temp 300.0 300.0 0.1
run             30000

print "=========================================="
print "Laser sintering simulation completed"
print "=========================================="

# Final structure analysis
print "Analyzing sintered structure and neck formation..."
unfix           particles_dynamics
unfix           final_equil
fix             final_nve particles nve
run             10000

# Final detailed dump - only 50 snapshots
dump            final all custom 1000 final_sintered_structure.lammpstrj id type x y z &
                v_temp_atom c_coord c_csym v_von_mises
dump_modify     final element Al Al Al

run             5000

# Write final configuration
write_data      final_sintered_particles.data

# Analysis output
variable        avg_coord_all equal ave(c_coord)
variable        avg_csym_all equal ave(c_csym)

print           "=========================================="
print           "Final Sintering Analysis:"
print           "Average coordination number: ${avg_coord_all}"
print           "Average centro-symmetry parameter: ${avg_csym_all}"
print           "=========================================="
print           "Visualization tips:"
print           "- Color by type to distinguish particles: Particle1=Red, Particle2=Blue"
print           "- Color by v_temp_atom to see temperature distribution"
print           "- Color by c_coord to identify neck formation (coord 10-12 = bonded)"
print           "- Color by c_csym to see grain boundaries at interface"
print           "- Check interface_atoms.lammpstrj for detailed neck evolution"
print           "=========================================="
