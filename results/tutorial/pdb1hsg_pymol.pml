# PyMOL script for visualizing predicted binding sites
load pdb1hsg_predicted.pdb, main_obj
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
select site_1_res, ((main_obj and chain A and resi 25) or (main_obj and chain A and resi 26) or (main_obj and chain A and resi 27) or (main_obj and chain A and resi 87) or (main_obj and chain A and resi 96) or (main_obj and chain A and resi 97) or (main_obj and chain A and resi 98) or (main_obj and chain A and resi 99) or (main_obj and chain B and resi 3) or (main_obj and chain B and resi 4) or (main_obj and chain B and resi 5) or (main_obj and chain B and resi 7) or (main_obj and chain B and resi 8) or (main_obj and chain B and resi 9) or (main_obj and chain B and resi 10) or (main_obj and chain B and resi 11) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 23) or (main_obj and chain B and resi 24) or (main_obj and chain B and resi 25) or (main_obj and chain B and resi 26) or (main_obj and chain B and resi 85) or (main_obj and chain B and resi 90) or (main_obj and chain B and resi 95) or (main_obj and chain B and resi 97))
show sticks, site_1_res
color lime, site_1_res
create site_1_surface, site_1_res
show surface, site_1_surface
set surface_color, lime, site_1_surface
set transparency, 0.3, site_1_surface

# Add label sphere at binding site center
pseudoatom label_obj_1, pos=[16.829999923706055, 32.00600051879883, -2.7510025024414064]
set sphere_scale, 0.01, label_obj_1
color black, label_obj_1
label label_obj_1, "Site 1 (Score: 13.9)"

# Site 1 details:
# Consensus Score: 3.90
# Binding Potential Score: 13.90
# Size: 5
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
select site_2_res, ((main_obj and chain A and resi 13) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 30) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 32) or (main_obj and chain A and resi 33) or (main_obj and chain A and resi 38) or (main_obj and chain A and resi 59) or (main_obj and chain A and resi 62) or (main_obj and chain A and resi 63) or (main_obj and chain A and resi 64) or (main_obj and chain A and resi 66) or (main_obj and chain A and resi 71) or (main_obj and chain A and resi 72) or (main_obj and chain A and resi 73) or (main_obj and chain A and resi 74) or (main_obj and chain A and resi 75) or (main_obj and chain A and resi 76) or (main_obj and chain A and resi 84) or (main_obj and chain A and resi 85) or (main_obj and chain A and resi 86) or (main_obj and chain A and resi 88) or (main_obj and chain A and resi 89) or (main_obj and chain A and resi 90))
show sticks, site_2_res
color limon, site_2_res
create site_2_surface, site_2_res
show surface, site_2_surface
set surface_color, limon, site_2_surface
set transparency, 0.3, site_2_surface

# Add label sphere at binding site center
pseudoatom label_obj_2, pos=[12.969534807426985, 35.51762842577557, 14.425741683605105]
set sphere_scale, 0.01, label_obj_2
color black, label_obj_2
label label_obj_2, "Site 2 (Score: 14.2)"

# Site 2 details:
# Consensus Score: 3.00
# Binding Potential Score: 14.20
# Size: 43
# Detection Methods: geometric, energy

# Binding site 3
select site_3_center, (main_obj and chain X and resi 3 and name O)
select site_3_points, (main_obj and chain X and resi 3 and name H)
show spheres, site_3_center
color black, site_3_center
set sphere_scale, 1.0, site_3_center
show spheres, site_3_points
color limon, site_3_points
set sphere_scale, 0.6, site_3_points
select site_3_res, ((main_obj and chain B and resi 11) or (main_obj and chain B and resi 12) or (main_obj and chain B and resi 13) or (main_obj and chain B and resi 15) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 21) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 23) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 32) or (main_obj and chain B and resi 33) or (main_obj and chain B and resi 34) or (main_obj and chain B and resi 36) or (main_obj and chain B and resi 64) or (main_obj and chain B and resi 75) or (main_obj and chain B and resi 76) or (main_obj and chain B and resi 80) or (main_obj and chain B and resi 82) or (main_obj and chain B and resi 83) or (main_obj and chain B and resi 84) or (main_obj and chain B and resi 85))
show sticks, site_3_res
color limon, site_3_res
create site_3_surface, site_3_res
show surface, site_3_surface
set surface_color, limon, site_3_surface
set transparency, 0.3, site_3_surface

# Add label sphere at binding site center
pseudoatom label_obj_3, pos=[12.829999923706055, 21.35266653696696, -8.05566660563151]
set sphere_scale, 0.01, label_obj_3
color black, label_obj_3
label label_obj_3, "Site 3 (Score: 10.7)"

# Site 3 details:
# Consensus Score: 3.00
# Binding Potential Score: 10.72
# Size: 3
# Detection Methods: energy, geometric

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
save pdb1hsg_predicted_binding_sites.pse
