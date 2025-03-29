# PyMOL script for visualizing predicted binding sites
load 1stp_predicted.pdb, main_obj
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
color salmon, site_1_points
set sphere_scale, 0.6, site_1_points
select site_1_res, ((main_obj and chain A and resi 23) or (main_obj and chain A and resi 27) or (main_obj and chain A and resi 28) or (main_obj and chain A and resi 29) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 44) or (main_obj and chain A and resi 45) or (main_obj and chain A and resi 56) or (main_obj and chain A and resi 75) or (main_obj and chain A and resi 76) or (main_obj and chain A and resi 77) or (main_obj and chain A and resi 79) or (main_obj and chain A and resi 90) or (main_obj and chain A and resi 91) or (main_obj and chain A and resi 92) or (main_obj and chain A and resi 93) or (main_obj and chain A and resi 106) or (main_obj and chain A and resi 107) or (main_obj and chain A and resi 108) or (main_obj and chain A and resi 127) or (main_obj and chain A and resi 128) or (main_obj and chain A and resi 130))
show sticks, site_1_res
color salmon, site_1_res
create site_1_surface, site_1_res
show surface, site_1_surface
set surface_color, salmon, site_1_surface
set transparency, 0.3, site_1_surface

# Add label sphere at binding site center
pseudoatom label_obj_1, pos=[14.101984508453853, 2.3535873473636686, -4.588349175831628]
set sphere_scale, 0.01, label_obj_1
color black, label_obj_1
label label_obj_1, "Site 1 (Score: 18.7)"

# Site 1 details:
# Consensus Score: 1.00
# Binding Potential Score: 18.70
# Size: 126
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
save 1stp_predicted_binding_sites.pse
