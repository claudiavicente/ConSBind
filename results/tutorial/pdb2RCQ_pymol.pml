# PyMOL script for visualizing predicted binding sites
load pdb2RCQ_predicted.pdb, main_obj
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
select site_1_res, ((main_obj and chain A and resi 1) or (main_obj and chain A and resi 3) or (main_obj and chain A and resi 4) or (main_obj and chain A and resi 8) or (main_obj and chain A and resi 40) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 91) or (main_obj and chain A and resi 107) or (main_obj and chain A and resi 108) or (main_obj and chain A and resi 109) or (main_obj and chain A and resi 110) or (main_obj and chain A and resi 111) or (main_obj and chain A and resi 112) or (main_obj and chain A and resi 113) or (main_obj and chain A and resi 114) or (main_obj and chain A and resi 115) or (main_obj and chain A and resi 116) or (main_obj and chain A and resi 117) or (main_obj and chain A and resi 128) or (main_obj and chain A and resi 129) or (main_obj and chain A and resi 130) or (main_obj and chain A and resi 131))
show sticks, site_1_res
color salmon, site_1_res
create site_1_surface, site_1_res
show surface, site_1_surface
set surface_color, salmon, site_1_surface
set transparency, 0.3, site_1_surface

# Add label sphere at binding site center
pseudoatom label_obj_1, pos=[15.903151132727182, -1.4346873785859795, -3.571956470448484]
set sphere_scale, 0.01, label_obj_1
color black, label_obj_1
label label_obj_1, "Site 1 (Score: 23.5)"

# Site 1 details:
# Consensus Score: 1.00
# Binding Potential Score: 23.50
# Size: 186
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
save pdb2RCQ_predicted_binding_sites.pse
