# PyMOL script for visualizing predicted binding sites
load pdb4dfr_predicted.pdb, main_obj
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
color forest, site_1_points
set sphere_scale, 0.6, site_1_points
select site_1_res, ((main_obj and chain B and resi 4) or (main_obj and chain B and resi 5) or (main_obj and chain B and resi 6) or (main_obj and chain B and resi 7) or (main_obj and chain B and resi 8) or (main_obj and chain B and resi 14) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 27) or (main_obj and chain B and resi 28) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 46) or (main_obj and chain B and resi 93) or (main_obj and chain B and resi 94) or (main_obj and chain B and resi 95) or (main_obj and chain B and resi 96) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 111) or (main_obj and chain B and resi 113) or (main_obj and chain B and resi 125))
show sticks, site_1_res
color forest, site_1_res
create site_1_surface, site_1_res
show surface, site_1_surface
set surface_color, forest, site_1_surface
set transparency, 0.3, site_1_surface

# Add label sphere at binding site center
pseudoatom label_obj_1, pos=[19.480666796366375, 65.48733266194661, 48.62966696421305]
set sphere_scale, 0.01, label_obj_1
color black, label_obj_1
label label_obj_1, "Site 1 (Score: 16.7)"

# Site 1 details:
# Consensus Score: 5.00
# Binding Potential Score: 16.66
# Size: 3
# Detection Methods: energy, geometric

# Binding site 2
select site_2_center, (main_obj and chain X and resi 2 and name O)
select site_2_points, (main_obj and chain X and resi 2 and name H)
show spheres, site_2_center
color black, site_2_center
set sphere_scale, 1.0, site_2_center
show spheres, site_2_points
color limon, site_2_points
set sphere_scale, 0.6, site_2_points
select site_2_res, ((main_obj and chain B and resi 8) or (main_obj and chain B and resi 24) or (main_obj and chain B and resi 113) or (main_obj and chain B and resi 114) or (main_obj and chain B and resi 115) or (main_obj and chain B and resi 116) or (main_obj and chain B and resi 117) or (main_obj and chain B and resi 139) or (main_obj and chain B and resi 140) or (main_obj and chain B and resi 141) or (main_obj and chain B and resi 142) or (main_obj and chain B and resi 143) or (main_obj and chain B and resi 147) or (main_obj and chain B and resi 148) or (main_obj and chain B and resi 149) or (main_obj and chain B and resi 150) or (main_obj and chain B and resi 151) or (main_obj and chain B and resi 152) or (main_obj and chain B and resi 153))
show sticks, site_2_res
color limon, site_2_res
create site_2_surface, site_2_res
show surface, site_2_surface
set surface_color, limon, site_2_surface
set transparency, 0.3, site_2_surface

# Add label sphere at binding site center
pseudoatom label_obj_2, pos=[11.240420793232165, 54.689420198139395, 44.06826345544113]
set sphere_scale, 0.01, label_obj_2
color black, label_obj_2
label label_obj_2, "Site 2 (Score: 14.9)"

# Site 2 details:
# Consensus Score: 3.00
# Binding Potential Score: 14.91
# Size: 19
# Detection Methods: geometric

# Binding site 3
select site_3_center, (main_obj and chain X and resi 3 and name O)
select site_3_points, (main_obj and chain X and resi 3 and name H)
show spheres, site_3_center
color black, site_3_center
set sphere_scale, 1.0, site_3_center
show spheres, site_3_points
color limon, site_3_points
set sphere_scale, 0.6, site_3_points
select site_3_res, ((main_obj and chain A and resi 5) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 32) or (main_obj and chain A and resi 34) or (main_obj and chain A and resi 35) or (main_obj and chain A and resi 36) or (main_obj and chain A and resi 37) or (main_obj and chain A and resi 38) or (main_obj and chain A and resi 39) or (main_obj and chain A and resi 40) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 50) or (main_obj and chain A and resi 52) or (main_obj and chain A and resi 53) or (main_obj and chain A and resi 54) or (main_obj and chain A and resi 55) or (main_obj and chain A and resi 56) or (main_obj and chain A and resi 57) or (main_obj and chain A and resi 58) or (main_obj and chain A and resi 59) or (main_obj and chain A and resi 60) or (main_obj and chain A and resi 92) or (main_obj and chain A and resi 94))
show sticks, site_3_res
color limon, site_3_res
create site_3_surface, site_3_res
show surface, site_3_surface
set surface_color, limon, site_3_surface
set transparency, 0.3, site_3_surface

# Add label sphere at binding site center
pseudoatom label_obj_3, pos=[17.871999740600586, 64.73276385139016, 17.72770617989933]
set sphere_scale, 0.01, label_obj_3
color black, label_obj_3
label label_obj_3, "Site 3 (Score: 12.1)"

# Site 3 details:
# Consensus Score: 3.00
# Binding Potential Score: 12.14
# Size: 17
# Detection Methods: geometric

# Binding site 4
select site_4_center, (main_obj and chain X and resi 4 and name O)
select site_4_points, (main_obj and chain X and resi 4 and name H)
show spheres, site_4_center
color black, site_4_center
set sphere_scale, 1.0, site_4_center
show spheres, site_4_points
color limon, site_4_points
set sphere_scale, 0.6, site_4_points
select site_4_res, ((main_obj and chain A and resi 7) or (main_obj and chain A and resi 14) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 18) or (main_obj and chain A and resi 19) or (main_obj and chain A and resi 20) or (main_obj and chain A and resi 21) or (main_obj and chain A and resi 22) or (main_obj and chain A and resi 23) or (main_obj and chain A and resi 24) or (main_obj and chain A and resi 49) or (main_obj and chain A and resi 115) or (main_obj and chain A and resi 119) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 23) or (main_obj and chain B and resi 145) or (main_obj and chain B and resi 146) or (main_obj and chain B and resi 147) or (main_obj and chain B and resi 148))
show sticks, site_4_res
color limon, site_4_res
create site_4_surface, site_4_res
show surface, site_4_surface
set surface_color, limon, site_4_surface
set transparency, 0.3, site_4_surface

# Add label sphere at binding site center
pseudoatom label_obj_4, pos=[18.948922817523663, 55.32099914550781, 31.039923374469463]
set sphere_scale, 0.01, label_obj_4
color black, label_obj_4
label label_obj_4, "Site 4 (Score: 11.4)"

# Site 4 details:
# Consensus Score: 3.00
# Binding Potential Score: 11.44
# Size: 13
# Detection Methods: geometric

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
save pdb4dfr_predicted_binding_sites.pse
