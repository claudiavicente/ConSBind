# PyMOL script for visualizing predicted binding sites
load 4dfr_predicted.pdb, main_obj
hide everything
bg_color white
set antialias, 2
set light_count, 2
set specular, 0.1
set sphere_quality, 2
set cartoon_fancy_helices, 1
set depth_cue, 0
set label_size, 10
set label_position, (0,0,0)
set label_color, black
set label_depth_mask, 0
set float_labels, on
show cartoon, main_obj and not chain X
color gray80, main_obj and not chain X
set cartoon_transparency, 0.5

# Show ligands in magenta
select ligands, main_obj and hetatm and not resn HOH and not chain X
show sticks, ligands
color magenta, ligands

# Binding site 1
select site_1_center, (main_obj and chain X and resi 1 and name O)
select site_1_points, (main_obj and chain X and resi 1 and name H)
show spheres, site_1_center
color black, site_1_center
set sphere_scale, 1.0, site_1_center
show spheres, site_1_points
color lime, site_1_points
set sphere_scale, 0.6, site_1_points
select site_1_res, ((main_obj and chain A and resi 2) or (main_obj and chain A and resi 3) or (main_obj and chain A and resi 4) or (main_obj and chain A and resi 5) or (main_obj and chain A and resi 6) or (main_obj and chain A and resi 30) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 33) or (main_obj and chain A and resi 34) or (main_obj and chain A and resi 35) or (main_obj and chain A and resi 38) or (main_obj and chain A and resi 90) or (main_obj and chain A and resi 91) or (main_obj and chain A and resi 92) or (main_obj and chain A and resi 93) or (main_obj and chain A and resi 94) or (main_obj and chain A and resi 109) or (main_obj and chain A and resi 110) or (main_obj and chain A and resi 111) or (main_obj and chain A and resi 112) or (main_obj and chain A and resi 154) or (main_obj and chain A and resi 155))
show sticks, site_1_res
color lime, site_1_res
create site_1_surface, site_1_res
show surface, site_1_surface
set surface_color, lime, site_1_surface
set transparency, 0.3, site_1_surface

# Add label sphere at binding site center
pseudoatom label_obj_1, pos=[26.73563610423695, 60.04827187278054, 15.190273024819113]
set sphere_scale, 0.01, label_obj_1
color black, label_obj_1
label label_obj_1, "Site 1 (Score: 15.3)"

# Site 1 details:
# Consensus Score: 3.90
# Binding Potential Score: 15.34
# Size: 22
# Detection Methods: geometric

# Binding site 2
select site_2_center, (main_obj and chain X and resi 2 and name O)
select site_2_points, (main_obj and chain X and resi 2 and name H)
show spheres, site_2_center
color black, site_2_center
set sphere_scale, 1.0, site_2_center
show spheres, site_2_points
color limon, site_2_points
set sphere_scale, 0.6, site_2_points
select site_2_res, ((main_obj and chain B and resi 4) or (main_obj and chain B and resi 5) or (main_obj and chain B and resi 6) or (main_obj and chain B and resi 7) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 35) or (main_obj and chain B and resi 40) or (main_obj and chain B and resi 41) or (main_obj and chain B and resi 42) or (main_obj and chain B and resi 43) or (main_obj and chain B and resi 46) or (main_obj and chain B and resi 50) or (main_obj and chain B and resi 54) or (main_obj and chain B and resi 92) or (main_obj and chain B and resi 93) or (main_obj and chain B and resi 94) or (main_obj and chain B and resi 95) or (main_obj and chain B and resi 96) or (main_obj and chain B and resi 100))
show sticks, site_2_res
color limon, site_2_res
create site_2_surface, site_2_res
show surface, site_2_surface
set surface_color, limon, site_2_surface
set transparency, 0.3, site_2_surface

# Add label sphere at binding site center
pseudoatom label_obj_2, pos=[19.729142597743444, 68.32099914550781, 49.39157172611782]
set sphere_scale, 0.01, label_obj_2
color black, label_obj_2
label label_obj_2, "Site 2 (Score: 10.8)"

# Site 2 details:
# Consensus Score: 3.00
# Binding Potential Score: 10.80
# Size: 7
# Detection Methods: geometric, energy

# Set up view
orient
zoom main_obj, 20
center main_obj

# Ray trace settings for better quality
set ray_shadows, 1
set ray_shadow_decay_factor, 0.1
set ray_trace_mode, 1
set ray_trace_color, black
set ray_trace_gain, 0.8
set ray_opaque_background, on

# Group all objects for easier manipulation
group binding_sites_centers, site_*_center
group binding_sites_points, site_*_points
group binding_sites_residues, site_*_res
group binding_sites, binding_sites_*
group labels, label_obj_*
group surfaces, *_surface

# Save session
save 4dfr_predicted_binding_sites.pse
