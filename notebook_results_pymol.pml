# PyMOL script for visualizing predicted binding sites
load notebook_results_predicted.pdb, main_obj
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
select site_1_res, ((main_obj and chain B and resi 32) or (main_obj and chain B and resi 35) or (main_obj and chain B and resi 36) or (main_obj and chain B and resi 37) or (main_obj and chain B and resi 38) or (main_obj and chain B and resi 39) or (main_obj and chain B and resi 40) or (main_obj and chain B and resi 54) or (main_obj and chain B and resi 55) or (main_obj and chain B and resi 56) or (main_obj and chain B and resi 57) or (main_obj and chain B and resi 58) or (main_obj and chain B and resi 59) or (main_obj and chain B and resi 90) or (main_obj and chain B and resi 92))
show sticks, site_1_res
color lime, site_1_res
create site_1_surface, site_1_res
show surface, site_1_surface
set surface_color, lime, site_1_surface
set transparency, 0.3, site_1_surface

# Add label sphere at binding site center
pseudoatom label_obj_1, pos=[14.971999740600586, 80.32099914550781, 47.46300029754639]
set sphere_scale, 0.01, label_obj_1
color black, label_obj_1
label label_obj_1, "Site 1 (Score: 14.3)"

# Site 1 details:
# Consensus Score: 3.90
# Binding Potential Score: 14.33
# Size: 10
# Detection Methods: geometric

# Binding site 2
select site_2_center, (main_obj and chain X and resi 2 and name O)
select site_2_points, (main_obj and chain X and resi 2 and name H)
show spheres, site_2_center
color black, site_2_center
set sphere_scale, 1.0, site_2_center
show spheres, site_2_points
color lime, site_2_points
set sphere_scale, 0.6, site_2_points
select site_2_res, ((main_obj and chain B and resi 2) or (main_obj and chain B and resi 4) or (main_obj and chain B and resi 41) or (main_obj and chain B and resi 62) or (main_obj and chain B and resi 77) or (main_obj and chain B and resi 78) or (main_obj and chain B and resi 79) or (main_obj and chain B and resi 80) or (main_obj and chain B and resi 81) or (main_obj and chain B and resi 82) or (main_obj and chain B and resi 91) or (main_obj and chain B and resi 93) or (main_obj and chain B and resi 99) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 101) or (main_obj and chain B and resi 102) or (main_obj and chain B and resi 103) or (main_obj and chain B and resi 104) or (main_obj and chain B and resi 105) or (main_obj and chain B and resi 106) or (main_obj and chain B and resi 107))
show sticks, site_2_res
color lime, site_2_res
create site_2_surface, site_2_res
show surface, site_2_surface
set surface_color, lime, site_2_surface
set transparency, 0.3, site_2_surface

# Add label sphere at binding site center
pseudoatom label_obj_2, pos=[21.471999740600587, 74.32099914550781, 61.763000297546384]
set sphere_scale, 0.01, label_obj_2
color black, label_obj_2
label label_obj_2, "Site 2 (Score: 14.2)"

# Site 2 details:
# Consensus Score: 3.90
# Binding Potential Score: 14.21
# Size: 5
# Detection Methods: geometric, energy

# Binding site 3
select site_3_center, (main_obj and chain X and resi 3 and name O)
select site_3_points, (main_obj and chain X and resi 3 and name H)
show spheres, site_3_center
color black, site_3_center
set sphere_scale, 1.0, site_3_center
show spheres, site_3_points
color lime, site_3_points
set sphere_scale, 0.6, site_3_points
select site_3_res, ((main_obj and chain A and resi 29) or (main_obj and chain A and resi 30) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 32) or (main_obj and chain A and resi 33) or (main_obj and chain A and resi 34) or (main_obj and chain A and resi 35) or (main_obj and chain A and resi 36) or (main_obj and chain A and resi 37) or (main_obj and chain A and resi 38) or (main_obj and chain A and resi 56) or (main_obj and chain A and resi 57) or (main_obj and chain A and resi 90) or (main_obj and chain A and resi 92))
show sticks, site_3_res
color lime, site_3_res
create site_3_surface, site_3_res
show surface, site_3_surface
set surface_color, lime, site_3_surface
set transparency, 0.3, site_3_surface

# Add label sphere at binding site center
pseudoatom label_obj_3, pos=[26.538666407267254, 68.32099914550781, 15.963000297546387]
set sphere_scale, 0.01, label_obj_3
color black, label_obj_3
label label_obj_3, "Site 3 (Score: 14.2)"

# Site 3 details:
# Consensus Score: 3.90
# Binding Potential Score: 14.20
# Size: 6
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
select site_4_res, ((main_obj and chain A and resi 2) or (main_obj and chain A and resi 3) or (main_obj and chain A and resi 4) or (main_obj and chain A and resi 5) or (main_obj and chain A and resi 39) or (main_obj and chain A and resi 40) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 60) or (main_obj and chain A and resi 91) or (main_obj and chain A and resi 92) or (main_obj and chain A and resi 93) or (main_obj and chain A and resi 94) or (main_obj and chain A and resi 95) or (main_obj and chain A and resi 96) or (main_obj and chain A and resi 99) or (main_obj and chain A and resi 100) or (main_obj and chain A and resi 103))
show sticks, site_4_res
color limon, site_4_res
create site_4_surface, site_4_res
show surface, site_4_surface
set surface_color, limon, site_4_surface
set transparency, 0.3, site_4_surface

# Add label sphere at binding site center
pseudoatom label_obj_4, pos=[18.071999740600585, 57.170999145507814, 14.163000297546386]
set sphere_scale, 0.01, label_obj_4
color black, label_obj_4
label label_obj_4, "Site 4 (Score: 12.1)"

# Site 4 details:
# Consensus Score: 3.00
# Binding Potential Score: 12.06
# Size: 20
# Detection Methods: geometric

# Binding site 5
select site_5_center, (main_obj and chain X and resi 5 and name O)
select site_5_points, (main_obj and chain X and resi 5 and name H)
show spheres, site_5_center
color black, site_5_center
set sphere_scale, 1.0, site_5_center
show spheres, site_5_points
color salmon, site_5_points
set sphere_scale, 0.6, site_5_points
select site_5_res, ((main_obj and chain A and resi 19) or (main_obj and chain A and resi 20) or (main_obj and chain A and resi 21) or (main_obj and chain A and resi 22) or (main_obj and chain A and resi 23) or (main_obj and chain A and resi 24) or (main_obj and chain A and resi 115) or (main_obj and chain A and resi 117) or (main_obj and chain A and resi 119) or (main_obj and chain A and resi 148) or (main_obj and chain A and resi 149) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 21) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 23) or (main_obj and chain B and resi 24) or (main_obj and chain B and resi 148) or (main_obj and chain B and resi 149))
show sticks, site_5_res
color salmon, site_5_res
create site_5_surface, site_5_res
show surface, site_5_surface
set surface_color, salmon, site_5_surface
set transparency, 0.3, site_5_surface

# Add label sphere at binding site center
pseudoatom label_obj_5, pos=[20.90836337696422, 54.957362781871446, 35.09027302481911]
set sphere_scale, 0.01, label_obj_5
color black, label_obj_5
label label_obj_5, "Site 5 (Score: 11.6)"

# Site 5 details:
# Consensus Score: 1.30
# Binding Potential Score: 11.61
# Size: 55
# Detection Methods: geometric

# Binding site 6
select site_6_center, (main_obj and chain X and resi 6 and name O)
select site_6_points, (main_obj and chain X and resi 6 and name H)
show spheres, site_6_center
color black, site_6_center
set sphere_scale, 1.0, site_6_center
show spheres, site_6_points
color salmon, site_6_points
set sphere_scale, 0.6, site_6_points
select site_6_res, ((main_obj and chain A and resi 22) or (main_obj and chain A and resi 23) or (main_obj and chain A and resi 24) or (main_obj and chain A and resi 25) or (main_obj and chain A and resi 26) or (main_obj and chain A and resi 115) or (main_obj and chain A and resi 116) or (main_obj and chain A and resi 141) or (main_obj and chain A and resi 142) or (main_obj and chain A and resi 143) or (main_obj and chain A and resi 144) or (main_obj and chain A and resi 145) or (main_obj and chain A and resi 146) or (main_obj and chain A and resi 147) or (main_obj and chain A and resi 148) or (main_obj and chain A and resi 149) or (main_obj and chain A and resi 150) or (main_obj and chain A and resi 151) or (main_obj and chain B and resi 18) or (main_obj and chain B and resi 19) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 21))
show sticks, site_6_res
color salmon, site_6_res
create site_6_surface, site_6_res
show surface, site_6_surface
set surface_color, salmon, site_6_surface
set transparency, 0.3, site_6_surface

# Add label sphere at binding site center
pseudoatom label_obj_6, pos=[29.621999740600586, 57.84877692328559, 35.21300029754639]
set sphere_scale, 0.01, label_obj_6
color black, label_obj_6
label label_obj_6, "Site 6 (Score: 11.3)"

# Site 6 details:
# Consensus Score: 1.30
# Binding Potential Score: 11.33
# Size: 36
# Detection Methods: geometric

# Binding site 7
select site_7_center, (main_obj and chain X and resi 7 and name O)
select site_7_points, (main_obj and chain X and resi 7 and name H)
show spheres, site_7_center
color black, site_7_center
set sphere_scale, 1.0, site_7_center
show spheres, site_7_points
color salmon, site_7_points
set sphere_scale, 0.6, site_7_points
select site_7_res, ((main_obj and chain B and resi 7) or (main_obj and chain B and resi 24) or (main_obj and chain B and resi 25) or (main_obj and chain B and resi 26) or (main_obj and chain B and resi 27) or (main_obj and chain B and resi 28) or (main_obj and chain B and resi 30) or (main_obj and chain B and resi 112) or (main_obj and chain B and resi 113) or (main_obj and chain B and resi 114) or (main_obj and chain B and resi 115) or (main_obj and chain B and resi 116) or (main_obj and chain B and resi 139) or (main_obj and chain B and resi 140) or (main_obj and chain B and resi 141) or (main_obj and chain B and resi 142) or (main_obj and chain B and resi 147) or (main_obj and chain B and resi 149) or (main_obj and chain B and resi 150) or (main_obj and chain B and resi 151) or (main_obj and chain B and resi 152) or (main_obj and chain B and resi 153) or (main_obj and chain B and resi 154))
show sticks, site_7_res
color salmon, site_7_res
create site_7_surface, site_7_res
show surface, site_7_surface
set surface_color, salmon, site_7_surface
set transparency, 0.3, site_7_surface

# Add label sphere at binding site center
pseudoatom label_obj_7, pos=[10.594221962822807, 59.82099914550781, 43.90744474199083]
set sphere_scale, 0.01, label_obj_7
color black, label_obj_7
label label_obj_7, "Site 7 (Score: 10.0)"

# Site 7 details:
# Consensus Score: 1.30
# Binding Potential Score: 10.05
# Size: 18
# Detection Methods: geometric

# Binding site 8
select site_8_center, (main_obj and chain X and resi 8 and name O)
select site_8_points, (main_obj and chain X and resi 8 and name H)
show spheres, site_8_center
color black, site_8_center
set sphere_scale, 1.0, site_8_center
show spheres, site_8_points
color salmon, site_8_points
set sphere_scale, 0.6, site_8_points
select site_8_res, ((main_obj and chain B and resi 26) or (main_obj and chain B and resi 27) or (main_obj and chain B and resi 30) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 33) or (main_obj and chain B and resi 34) or (main_obj and chain B and resi 111) or (main_obj and chain B and resi 112) or (main_obj and chain B and resi 113) or (main_obj and chain B and resi 114) or (main_obj and chain B and resi 137) or (main_obj and chain B and resi 138) or (main_obj and chain B and resi 139) or (main_obj and chain B and resi 140) or (main_obj and chain B and resi 141) or (main_obj and chain B and resi 151) or (main_obj and chain B and resi 152) or (main_obj and chain B and resi 153) or (main_obj and chain B and resi 154) or (main_obj and chain B and resi 155))
show sticks, site_8_res
color salmon, site_8_res
create site_8_surface, site_8_res
show surface, site_8_surface
set surface_color, salmon, site_8_surface
set transparency, 0.3, site_8_surface

# Add label sphere at binding site center
pseudoatom label_obj_8, pos=[6.781090649691495, 63.50281732732599, 47.41754575209184]
set sphere_scale, 0.01, label_obj_8
color black, label_obj_8
label label_obj_8, "Site 8 (Score: 9.9)"

# Site 8 details:
# Consensus Score: 1.30
# Binding Potential Score: 9.93
# Size: 11
# Detection Methods: geometric

# Binding site 9
select site_9_center, (main_obj and chain X and resi 9 and name O)
select site_9_points, (main_obj and chain X and resi 9 and name H)
show spheres, site_9_center
color black, site_9_center
set sphere_scale, 1.0, site_9_center
show spheres, site_9_points
color salmon, site_9_points
set sphere_scale, 0.6, site_9_points
select site_9_res, ((main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 17) or (main_obj and chain A and resi 18) or (main_obj and chain A and resi 19) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 44) or (main_obj and chain A and resi 45) or (main_obj and chain A and resi 46) or (main_obj and chain A and resi 47) or (main_obj and chain A and resi 48) or (main_obj and chain A and resi 49) or (main_obj and chain A and resi 50) or (main_obj and chain A and resi 51) or (main_obj and chain A and resi 96) or (main_obj and chain B and resi 144) or (main_obj and chain B and resi 145) or (main_obj and chain B and resi 146))
show sticks, site_9_res
color salmon, site_9_res
create site_9_surface, site_9_res
show surface, site_9_surface
set surface_color, salmon, site_9_surface
set transparency, 0.3, site_9_surface

# Add label sphere at binding site center
pseudoatom label_obj_9, pos=[10.417454286055131, 56.04827187278054, 25.14481847936457]
set sphere_scale, 0.01, label_obj_9
color black, label_obj_9
label label_obj_9, "Site 9 (Score: 9.9)"

# Site 9 details:
# Consensus Score: 1.30
# Binding Potential Score: 9.87
# Size: 11
# Detection Methods: geometric

# Binding site 10
select site_10_center, (main_obj and chain X and resi 10 and name O)
select site_10_points, (main_obj and chain X and resi 10 and name H)
show spheres, site_10_center
color black, site_10_center
set sphere_scale, 1.0, site_10_center
show spheres, site_10_points
color salmon, site_10_points
set sphere_scale, 0.6, site_10_points
select site_10_res, ((main_obj and chain B and resi 1) or (main_obj and chain B and resi 2) or (main_obj and chain B and resi 3) or (main_obj and chain B and resi 4) or (main_obj and chain B and resi 5) or (main_obj and chain B and resi 89) or (main_obj and chain B and resi 90) or (main_obj and chain B and resi 91) or (main_obj and chain B and resi 92) or (main_obj and chain B and resi 93) or (main_obj and chain B and resi 103) or (main_obj and chain B and resi 104) or (main_obj and chain B and resi 106) or (main_obj and chain B and resi 107) or (main_obj and chain B and resi 108) or (main_obj and chain B and resi 109) or (main_obj and chain B and resi 110) or (main_obj and chain B and resi 111) or (main_obj and chain B and resi 133) or (main_obj and chain B and resi 156) or (main_obj and chain B and resi 157) or (main_obj and chain B and resi 158))
show sticks, site_10_res
color salmon, site_10_res
create site_10_surface, site_10_res
show surface, site_10_surface
set surface_color, salmon, site_10_surface
set transparency, 0.3, site_10_surface

# Add label sphere at binding site center
pseudoatom label_obj_10, pos=[12.271999740600586, 70.97814200265067, 59.22014315468925]
set sphere_scale, 0.01, label_obj_10
color black, label_obj_10
label label_obj_10, "Site 10 (Score: 8.9)"

# Site 10 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.89
# Size: 35
# Detection Methods: geometric

# Binding site 11
select site_11_center, (main_obj and chain X and resi 11 and name O)
select site_11_points, (main_obj and chain X and resi 11 and name H)
show spheres, site_11_center
color black, site_11_center
set sphere_scale, 1.0, site_11_center
show spheres, site_11_points
color salmon, site_11_points
set sphere_scale, 0.6, site_11_points
select site_11_res, ((main_obj and chain B and resi 8) or (main_obj and chain B and resi 9) or (main_obj and chain B and resi 10) or (main_obj and chain B and resi 114) or (main_obj and chain B and resi 115) or (main_obj and chain B and resi 116) or (main_obj and chain B and resi 117) or (main_obj and chain B and resi 118) or (main_obj and chain B and resi 140) or (main_obj and chain B and resi 141) or (main_obj and chain B and resi 142) or (main_obj and chain B and resi 143) or (main_obj and chain B and resi 149) or (main_obj and chain B and resi 150) or (main_obj and chain B and resi 151) or (main_obj and chain B and resi 152))
show sticks, site_11_res
color salmon, site_11_res
create site_11_surface, site_11_res
show surface, site_11_surface
set surface_color, salmon, site_11_surface
set transparency, 0.3, site_11_surface

# Add label sphere at binding site center
pseudoatom label_obj_11, pos=[12.621999740600586, 51.44599914550781, 44.83800029754639]
set sphere_scale, 0.01, label_obj_11
color black, label_obj_11
label label_obj_11, "Site 11 (Score: 8.7)"

# Site 11 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.75
# Size: 8
# Detection Methods: geometric

# Binding site 12
select site_12_center, (main_obj and chain X and resi 12 and name O)
select site_12_points, (main_obj and chain X and resi 12 and name H)
show spheres, site_12_center
color black, site_12_center
set sphere_scale, 1.0, site_12_center
show spheres, site_12_points
color salmon, site_12_points
set sphere_scale, 0.6, site_12_points
select site_12_res, ((main_obj and chain A and resi 4) or (main_obj and chain A and resi 5) or (main_obj and chain A and resi 6) or (main_obj and chain A and resi 7) or (main_obj and chain A and resi 26) or (main_obj and chain A and resi 27) or (main_obj and chain A and resi 29) or (main_obj and chain A and resi 30) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 32) or (main_obj and chain A and resi 34) or (main_obj and chain A and resi 92) or (main_obj and chain A and resi 110) or (main_obj and chain A and resi 111) or (main_obj and chain A and resi 112) or (main_obj and chain A and resi 113) or (main_obj and chain A and resi 114) or (main_obj and chain A and resi 137) or (main_obj and chain A and resi 152) or (main_obj and chain A and resi 153) or (main_obj and chain A and resi 154) or (main_obj and chain A and resi 155))
show sticks, site_12_res
color salmon, site_12_res
create site_12_surface, site_12_res
show surface, site_12_surface
set surface_color, salmon, site_12_surface
set transparency, 0.3, site_12_surface

# Add label sphere at binding site center
pseudoatom label_obj_12, pos=[28.071999740600585, 58.72099914550781, 21.163000297546386]
set sphere_scale, 0.01, label_obj_12
color black, label_obj_12
label label_obj_12, "Site 12 (Score: 8.5)"

# Site 12 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.53
# Size: 5
# Detection Methods: geometric

# Binding site 13
select site_13_center, (main_obj and chain X and resi 13 and name O)
select site_13_points, (main_obj and chain X and resi 13 and name H)
show spheres, site_13_center
color black, site_13_center
set sphere_scale, 1.0, site_13_center
show spheres, site_13_points
color salmon, site_13_points
set sphere_scale, 0.6, site_13_points
select site_13_res, ((main_obj and chain A and resi 7) or (main_obj and chain A and resi 8) or (main_obj and chain A and resi 9) or (main_obj and chain A and resi 10) or (main_obj and chain A and resi 11) or (main_obj and chain A and resi 12) or (main_obj and chain A and resi 13) or (main_obj and chain A and resi 14) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 20) or (main_obj and chain A and resi 22) or (main_obj and chain A and resi 97) or (main_obj and chain A and resi 100) or (main_obj and chain A and resi 119) or (main_obj and chain A and resi 120) or (main_obj and chain A and resi 121) or (main_obj and chain A and resi 122) or (main_obj and chain A and resi 123) or (main_obj and chain A and resi 124) or (main_obj and chain A and resi 125) or (main_obj and chain A and resi 126))
show sticks, site_13_res
color salmon, site_13_res
create site_13_surface, site_13_res
show surface, site_13_surface
set surface_color, salmon, site_13_surface
set transparency, 0.3, site_13_surface

# Add label sphere at binding site center
pseudoatom label_obj_13, pos=[20.871999740600586, 47.19599914550781, 25.588000297546387]
set sphere_scale, 0.01, label_obj_13
color black, label_obj_13
label label_obj_13, "Site 13 (Score: 8.5)"

# Site 13 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.52
# Size: 8
# Detection Methods: geometric

# Binding site 14
select site_14_center, (main_obj and chain X and resi 14 and name O)
select site_14_points, (main_obj and chain X and resi 14 and name H)
show spheres, site_14_center
color black, site_14_center
set sphere_scale, 1.0, site_14_center
show spheres, site_14_points
color salmon, site_14_points
set sphere_scale, 0.6, site_14_points
select site_14_res, ((main_obj and chain A and resi 12) or (main_obj and chain A and resi 13) or (main_obj and chain A and resi 14) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 17) or (main_obj and chain A and resi 97) or (main_obj and chain A and resi 98) or (main_obj and chain A and resi 101) or (main_obj and chain A and resi 120) or (main_obj and chain A and resi 121) or (main_obj and chain A and resi 122) or (main_obj and chain A and resi 123) or (main_obj and chain A and resi 124) or (main_obj and chain A and resi 125))
show sticks, site_14_res
color salmon, site_14_res
create site_14_surface, site_14_res
show surface, site_14_surface
set surface_color, salmon, site_14_surface
set transparency, 0.3, site_14_surface

# Add label sphere at binding site center
pseudoatom label_obj_14, pos=[16.371999740600586, 44.94599914550781, 24.713000297546387]
set sphere_scale, 0.01, label_obj_14
color black, label_obj_14
label label_obj_14, "Site 14 (Score: 8.3)"

# Site 14 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.34
# Size: 8
# Detection Methods: geometric

# Binding site 15
select site_15_center, (main_obj and chain X and resi 15 and name O)
select site_15_points, (main_obj and chain X and resi 15 and name H)
show spheres, site_15_center
color black, site_15_center
set sphere_scale, 1.0, site_15_center
show spheres, site_15_points
color salmon, site_15_points
set sphere_scale, 0.6, site_15_points
select site_15_res, ((main_obj and chain B and resi 7) or (main_obj and chain B and resi 8) or (main_obj and chain B and resi 9) or (main_obj and chain B and resi 10) or (main_obj and chain B and resi 11) or (main_obj and chain B and resi 12) or (main_obj and chain B and resi 13) or (main_obj and chain B and resi 14) or (main_obj and chain B and resi 15) or (main_obj and chain B and resi 16) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 96) or (main_obj and chain B and resi 97) or (main_obj and chain B and resi 98) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 121) or (main_obj and chain B and resi 122) or (main_obj and chain B and resi 123) or (main_obj and chain B and resi 124) or (main_obj and chain B and resi 125) or (main_obj and chain B and resi 126))
show sticks, site_15_res
color salmon, site_15_res
create site_15_surface, site_15_res
show surface, site_15_surface
set surface_color, salmon, site_15_surface
set transparency, 0.3, site_15_surface

# Add label sphere at binding site center
pseudoatom label_obj_15, pos=[24.621999740600586, 58.57099914550781, 52.71300029754639]
set sphere_scale, 0.01, label_obj_15
color black, label_obj_15
label label_obj_15, "Site 15 (Score: 8.2)"

# Site 15 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.17
# Size: 8
# Detection Methods: geometric

# Binding site 16
select site_16_center, (main_obj and chain X and resi 16 and name O)
select site_16_points, (main_obj and chain X and resi 16 and name H)
show spheres, site_16_center
color black, site_16_center
set sphere_scale, 1.0, site_16_center
show spheres, site_16_points
color salmon, site_16_points
set sphere_scale, 0.6, site_16_points
select site_16_res, ((main_obj and chain A and resi 14) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 17) or (main_obj and chain A and resi 18) or (main_obj and chain A and resi 19) or (main_obj and chain A and resi 20) or (main_obj and chain A and resi 21) or (main_obj and chain A and resi 22) or (main_obj and chain A and resi 119) or (main_obj and chain A and resi 120) or (main_obj and chain B and resi 23) or (main_obj and chain B and resi 143) or (main_obj and chain B and resi 144) or (main_obj and chain B and resi 145) or (main_obj and chain B and resi 146) or (main_obj and chain B and resi 147) or (main_obj and chain B and resi 148) or (main_obj and chain B and resi 149))
show sticks, site_16_res
color salmon, site_16_res
create site_16_surface, site_16_res
show surface, site_16_surface
set surface_color, salmon, site_16_surface
set transparency, 0.3, site_16_surface

# Add label sphere at binding site center
pseudoatom label_obj_16, pos=[16.271999740600585, 52.12099914550781, 33.16300029754639]
set sphere_scale, 0.01, label_obj_16
color black, label_obj_16
label label_obj_16, "Site 16 (Score: 8.2)"

# Site 16 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.17
# Size: 5
# Detection Methods: geometric

# Binding site 17
select site_17_center, (main_obj and chain X and resi 17 and name O)
select site_17_points, (main_obj and chain X and resi 17 and name H)
show spheres, site_17_center
color black, site_17_center
set sphere_scale, 1.0, site_17_center
show spheres, site_17_points
color salmon, site_17_points
set sphere_scale, 0.6, site_17_points
select site_17_res, ((main_obj and chain A and resi 23) or (main_obj and chain A and resi 143) or (main_obj and chain A and resi 144) or (main_obj and chain A and resi 145) or (main_obj and chain A and resi 146) or (main_obj and chain A and resi 147) or (main_obj and chain A and resi 148) or (main_obj and chain A and resi 149) or (main_obj and chain A and resi 150) or (main_obj and chain B and resi 15) or (main_obj and chain B and resi 16) or (main_obj and chain B and resi 17) or (main_obj and chain B and resi 18) or (main_obj and chain B and resi 19) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 21) or (main_obj and chain B and resi 22))
show sticks, site_17_res
color salmon, site_17_res
create site_17_surface, site_17_res
show surface, site_17_surface
set surface_color, salmon, site_17_surface
set transparency, 0.3, site_17_surface

# Add label sphere at binding site center
pseudoatom label_obj_17, pos=[28.471999740600587, 57.72099914550781, 41.16300029754639]
set sphere_scale, 0.01, label_obj_17
color black, label_obj_17
label label_obj_17, "Site 17 (Score: 8.1)"

# Site 17 details:
# Consensus Score: 1.30
# Binding Potential Score: 8.13
# Size: 5
# Detection Methods: geometric

# Binding site 18
select site_18_center, (main_obj and chain X and resi 18 and name O)
select site_18_points, (main_obj and chain X and resi 18 and name H)
show spheres, site_18_center
color black, site_18_center
set sphere_scale, 1.0, site_18_center
show spheres, site_18_points
color salmon, site_18_points
set sphere_scale, 0.6, site_18_points
select site_18_res, ((main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 17) or (main_obj and chain A and resi 18) or (main_obj and chain A and resi 19) or (main_obj and chain A and resi 20) or (main_obj and chain A and resi 45) or (main_obj and chain A and resi 46) or (main_obj and chain A and resi 48) or (main_obj and chain A and resi 49) or (main_obj and chain B and resi 143) or (main_obj and chain B and resi 144) or (main_obj and chain B and resi 145) or (main_obj and chain B and resi 146) or (main_obj and chain B and resi 147))
show sticks, site_18_res
color salmon, site_18_res
create site_18_surface, site_18_res
show surface, site_18_surface
set surface_color, salmon, site_18_surface
set transparency, 0.3, site_18_surface

# Add label sphere at binding site center
pseudoatom label_obj_18, pos=[11.157714026314872, 52.32099914550781, 29.10585744040353]
set sphere_scale, 0.01, label_obj_18
color black, label_obj_18
label label_obj_18, "Site 18 (Score: 8.0)"

# Site 18 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.98
# Size: 7
# Detection Methods: geometric

# Binding site 19
select site_19_center, (main_obj and chain X and resi 19 and name O)
select site_19_points, (main_obj and chain X and resi 19 and name H)
show spheres, site_19_center
color black, site_19_center
set sphere_scale, 1.0, site_19_center
show spheres, site_19_points
color salmon, site_19_points
set sphere_scale, 0.6, site_19_points
select site_19_res, ((main_obj and chain B and resi 40) or (main_obj and chain B and resi 41) or (main_obj and chain B and resi 42) or (main_obj and chain B and resi 43) or (main_obj and chain B and resi 44) or (main_obj and chain B and resi 45) or (main_obj and chain B and resi 46) or (main_obj and chain B and resi 47) or (main_obj and chain B and resi 59) or (main_obj and chain B and resi 60) or (main_obj and chain B and resi 61) or (main_obj and chain B and resi 62) or (main_obj and chain B and resi 63) or (main_obj and chain B and resi 65) or (main_obj and chain B and resi 72) or (main_obj and chain B and resi 73) or (main_obj and chain B and resi 74) or (main_obj and chain B and resi 75))
show sticks, site_19_res
color salmon, site_19_res
create site_19_surface, site_19_res
show surface, site_19_surface
set surface_color, salmon, site_19_surface
set transparency, 0.3, site_19_surface

# Add label sphere at binding site center
pseudoatom label_obj_19, pos=[27.721999740600587, 77.1709991455078, 50.71300029754639]
set sphere_scale, 0.01, label_obj_19
color black, label_obj_19
label label_obj_19, "Site 19 (Score: 7.9)"

# Site 19 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.92
# Size: 20
# Detection Methods: geometric

# Binding site 20
select site_20_center, (main_obj and chain X and resi 20 and name O)
select site_20_points, (main_obj and chain X and resi 20 and name H)
show spheres, site_20_center
color black, site_20_center
set sphere_scale, 1.0, site_20_center
show spheres, site_20_points
color salmon, site_20_points
set sphere_scale, 0.6, site_20_points
select site_20_res, ((main_obj and chain A and resi 145) or (main_obj and chain B and resi 14) or (main_obj and chain B and resi 15) or (main_obj and chain B and resi 16) or (main_obj and chain B and resi 17) or (main_obj and chain B and resi 18) or (main_obj and chain B and resi 19) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 45) or (main_obj and chain B and resi 46) or (main_obj and chain B and resi 49) or (main_obj and chain B and resi 96) or (main_obj and chain B and resi 97) or (main_obj and chain B and resi 122) or (main_obj and chain B and resi 123))
show sticks, site_20_res
color salmon, site_20_res
create site_20_surface, site_20_res
show surface, site_20_surface
set surface_color, salmon, site_20_surface
set transparency, 0.3, site_20_surface

# Add label sphere at binding site center
pseudoatom label_obj_20, pos=[27.705333073933918, 62.65433247884115, 47.12966696421305]
set sphere_scale, 0.01, label_obj_20
color black, label_obj_20
label label_obj_20, "Site 20 (Score: 7.9)"

# Site 20 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.88
# Size: 6
# Detection Methods: geometric

# Binding site 21
select site_21_center, (main_obj and chain X and resi 21 and name O)
select site_21_points, (main_obj and chain X and resi 21 and name H)
show spheres, site_21_center
color black, site_21_center
set sphere_scale, 1.0, site_21_center
show spheres, site_21_points
color salmon, site_21_points
set sphere_scale, 0.6, site_21_points
select site_21_res, ((main_obj and chain A and resi 14) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 17) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 45) or (main_obj and chain A and resi 46) or (main_obj and chain A and resi 94) or (main_obj and chain A and resi 95) or (main_obj and chain A and resi 96) or (main_obj and chain A and resi 97) or (main_obj and chain A and resi 98) or (main_obj and chain A and resi 99) or (main_obj and chain A and resi 100) or (main_obj and chain A and resi 101) or (main_obj and chain A and resi 123) or (main_obj and chain A and resi 124) or (main_obj and chain A and resi 125) or (main_obj and chain A and resi 126))
show sticks, site_21_res
color salmon, site_21_res
create site_21_surface, site_21_res
show surface, site_21_surface
set surface_color, salmon, site_21_surface
set transparency, 0.3, site_21_surface

# Add label sphere at binding site center
pseudoatom label_obj_21, pos=[15.094221962822807, 49.43211025661893, 20.963000297546387]
set sphere_scale, 0.01, label_obj_21
color black, label_obj_21
label label_obj_21, "Site 21 (Score: 7.9)"

# Site 21 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.87
# Size: 9
# Detection Methods: geometric

# Binding site 22
select site_22_center, (main_obj and chain X and resi 22 and name O)
select site_22_points, (main_obj and chain X and resi 22 and name H)
show spheres, site_22_center
color black, site_22_center
set sphere_scale, 1.0, site_22_center
show spheres, site_22_points
color salmon, site_22_points
set sphere_scale, 0.6, site_22_points
select site_22_res, ((main_obj and chain B and resi 42) or (main_obj and chain B and resi 43) or (main_obj and chain B and resi 44) or (main_obj and chain B and resi 45) or (main_obj and chain B and resi 46) or (main_obj and chain B and resi 47) or (main_obj and chain B and resi 48) or (main_obj and chain B and resi 49) or (main_obj and chain B and resi 50) or (main_obj and chain B and resi 61) or (main_obj and chain B and resi 62) or (main_obj and chain B and resi 94) or (main_obj and chain B and resi 95) or (main_obj and chain B and resi 96) or (main_obj and chain B and resi 97) or (main_obj and chain B and resi 99))
show sticks, site_22_res
color salmon, site_22_res
create site_22_surface, site_22_res
show surface, site_22_surface
set surface_color, salmon, site_22_surface
set transparency, 0.3, site_22_surface

# Add label sphere at binding site center
pseudoatom label_obj_22, pos=[28.23563610423695, 68.86645369096236, 48.78118211572821]
set sphere_scale, 0.01, label_obj_22
color black, label_obj_22
label label_obj_22, "Site 22 (Score: 7.8)"

# Site 22 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.85
# Size: 11
# Detection Methods: geometric

# Binding site 23
select site_23_center, (main_obj and chain X and resi 23 and name O)
select site_23_points, (main_obj and chain X and resi 23 and name H)
show spheres, site_23_center
color black, site_23_center
set sphere_scale, 1.0, site_23_center
show spheres, site_23_points
color salmon, site_23_points
set sphere_scale, 0.6, site_23_points
select site_23_res, ((main_obj and chain A and resi 7) or (main_obj and chain A and resi 12) or (main_obj and chain A and resi 13) or (main_obj and chain A and resi 14) or (main_obj and chain A and resi 15) or (main_obj and chain A and resi 16) or (main_obj and chain A and resi 17) or (main_obj and chain A and resi 18) or (main_obj and chain A and resi 19) or (main_obj and chain A and resi 20) or (main_obj and chain A and resi 22) or (main_obj and chain A and resi 45) or (main_obj and chain A and resi 46) or (main_obj and chain A and resi 49) or (main_obj and chain A and resi 96) or (main_obj and chain A and resi 97) or (main_obj and chain A and resi 100) or (main_obj and chain A and resi 119) or (main_obj and chain A and resi 123) or (main_obj and chain A and resi 124))
show sticks, site_23_res
color salmon, site_23_res
create site_23_surface, site_23_res
show surface, site_23_surface
set surface_color, salmon, site_23_surface
set transparency, 0.3, site_23_surface

# Add label sphere at binding site center
pseudoatom label_obj_23, pos=[16.871999740600586, 50.920999145507814, 26.163000297546386]
set sphere_scale, 0.01, label_obj_23
color black, label_obj_23
label label_obj_23, "Site 23 (Score: 7.8)"

# Site 23 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.83
# Size: 5
# Detection Methods: geometric

# Binding site 24
select site_24_center, (main_obj and chain X and resi 24 and name O)
select site_24_points, (main_obj and chain X and resi 24 and name H)
show spheres, site_24_center
color black, site_24_center
set sphere_scale, 1.0, site_24_center
show spheres, site_24_points
color salmon, site_24_points
set sphere_scale, 0.6, site_24_points
select site_24_res, ((main_obj and chain B and resi 6) or (main_obj and chain B and resi 7) or (main_obj and chain B and resi 8) or (main_obj and chain B and resi 9) or (main_obj and chain B and resi 10) or (main_obj and chain B and resi 12) or (main_obj and chain B and resi 13) or (main_obj and chain B and resi 14) or (main_obj and chain B and resi 15) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 24) or (main_obj and chain B and resi 27) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 112) or (main_obj and chain B and resi 113) or (main_obj and chain B and resi 114) or (main_obj and chain B and resi 115) or (main_obj and chain B and resi 119) or (main_obj and chain B and resi 125))
show sticks, site_24_res
color salmon, site_24_res
create site_24_surface, site_24_res
show surface, site_24_surface
set surface_color, salmon, site_24_surface
set transparency, 0.3, site_24_surface

# Add label sphere at binding site center
pseudoatom label_obj_24, pos=[18.60533307393392, 59.520999145507815, 48.363000297546385]
set sphere_scale, 0.01, label_obj_24
color black, label_obj_24
label label_obj_24, "Site 24 (Score: 7.7)"

# Site 24 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.69
# Size: 15
# Detection Methods: geometric

# Binding site 25
select site_25_center, (main_obj and chain X and resi 25 and name O)
select site_25_points, (main_obj and chain X and resi 25 and name H)
show spheres, site_25_center
color black, site_25_center
set sphere_scale, 1.0, site_25_center
show spheres, site_25_points
color salmon, site_25_points
set sphere_scale, 0.6, site_25_points
select site_25_res, ((main_obj and chain B and resi 8) or (main_obj and chain B and resi 9) or (main_obj and chain B and resi 10) or (main_obj and chain B and resi 11) or (main_obj and chain B and resi 13) or (main_obj and chain B and resi 16) or (main_obj and chain B and resi 20) or (main_obj and chain B and resi 21) or (main_obj and chain B and resi 22) or (main_obj and chain B and resi 115) or (main_obj and chain B and resi 116) or (main_obj and chain B and resi 117) or (main_obj and chain B and resi 118) or (main_obj and chain B and resi 119) or (main_obj and chain B and resi 120) or (main_obj and chain B and resi 121) or (main_obj and chain B and resi 149))
show sticks, site_25_res
color salmon, site_25_res
create site_25_surface, site_25_res
show surface, site_25_surface
set surface_color, salmon, site_25_surface
set transparency, 0.3, site_25_surface

# Add label sphere at binding site center
pseudoatom label_obj_25, pos=[21.171999740600587, 51.62099914550781, 45.96300029754639]
set sphere_scale, 0.01, label_obj_25
color black, label_obj_25
label label_obj_25, "Site 25 (Score: 7.5)"

# Site 25 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.52
# Size: 10
# Detection Methods: geometric

# Binding site 26
select site_26_center, (main_obj and chain X and resi 26 and name O)
select site_26_points, (main_obj and chain X and resi 26 and name H)
show spheres, site_26_center
color black, site_26_center
set sphere_scale, 1.0, site_26_center
show spheres, site_26_points
color salmon, site_26_points
set sphere_scale, 0.6, site_26_points
select site_26_res, ((main_obj and chain A and resi 7) or (main_obj and chain A and resi 22) or (main_obj and chain A and resi 24) or (main_obj and chain A and resi 27) or (main_obj and chain A and resi 28) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 32) or (main_obj and chain A and resi 46) or (main_obj and chain A and resi 49) or (main_obj and chain A and resi 50) or (main_obj and chain A and resi 54) or (main_obj and chain A and resi 94) or (main_obj and chain B and resi 23))
show sticks, site_26_res
color salmon, site_26_res
create site_26_surface, site_26_res
show surface, site_26_surface
set surface_color, salmon, site_26_surface
set transparency, 0.3, site_26_surface

# Add label sphere at binding site center
pseudoatom label_obj_26, pos=[19.187789214284795, 60.95257809287623, 25.699842402809544]
set sphere_scale, 0.01, label_obj_26
color black, label_obj_26
label label_obj_26, "Site 26 (Score: 7.2)"

# Site 26 details:
# Consensus Score: 1.30
# Binding Potential Score: 7.24
# Size: 19
# Detection Methods: geometric

# Binding site 27
select site_27_center, (main_obj and chain X and resi 27 and name O)
select site_27_points, (main_obj and chain X and resi 27 and name H)
show spheres, site_27_center
color black, site_27_center
set sphere_scale, 1.0, site_27_center
show spheres, site_27_points
color salmon, site_27_points
set sphere_scale, 0.6, site_27_points
select site_27_res, ((main_obj and chain B and resi 6) or (main_obj and chain B and resi 8) or (main_obj and chain B and resi 9) or (main_obj and chain B and resi 11) or (main_obj and chain B and resi 12) or (main_obj and chain B and resi 13) or (main_obj and chain B and resi 14) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 104) or (main_obj and chain B and resi 110) or (main_obj and chain B and resi 112) or (main_obj and chain B and resi 114) or (main_obj and chain B and resi 124) or (main_obj and chain B and resi 125) or (main_obj and chain B and resi 126) or (main_obj and chain B and resi 127) or (main_obj and chain B and resi 128))
show sticks, site_27_res
color salmon, site_27_res
create site_27_surface, site_27_res
show surface, site_27_surface
set surface_color, salmon, site_27_surface
set transparency, 0.3, site_27_surface

# Add label sphere at binding site center
pseudoatom label_obj_27, pos=[18.300571169172013, 58.03528485979353, 57.10585744040353]
set sphere_scale, 0.01, label_obj_27
color black, label_obj_27
label label_obj_27, "Site 27 (Score: 6.9)"

# Site 27 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.94
# Size: 7
# Detection Methods: geometric

# Binding site 28
select site_28_center, (main_obj and chain X and resi 28 and name O)
select site_28_points, (main_obj and chain X and resi 28 and name H)
show spheres, site_28_center
color black, site_28_center
set sphere_scale, 1.0, site_28_center
show spheres, site_28_points
color salmon, site_28_points
set sphere_scale, 0.6, site_28_points
select site_28_res, ((main_obj and chain A and resi 39) or (main_obj and chain A and resi 40) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 58) or (main_obj and chain A and resi 60) or (main_obj and chain A and resi 73) or (main_obj and chain A and resi 75) or (main_obj and chain A and resi 77) or (main_obj and chain A and resi 78) or (main_obj and chain A and resi 79) or (main_obj and chain A and resi 80) or (main_obj and chain A and resi 81) or (main_obj and chain A and resi 82) or (main_obj and chain A and resi 83) or (main_obj and chain A and resi 84) or (main_obj and chain A and resi 85) or (main_obj and chain A and resi 86) or (main_obj and chain A and resi 87) or (main_obj and chain A and resi 88) or (main_obj and chain A and resi 91) or (main_obj and chain A and resi 103))
show sticks, site_28_res
color salmon, site_28_res
create site_28_surface, site_28_res
show surface, site_28_surface
set surface_color, salmon, site_28_surface
set transparency, 0.3, site_28_surface

# Add label sphere at binding site center
pseudoatom label_obj_28, pos=[14.121999740600586, 61.23766581217448, 6.29633363087972]
set sphere_scale, 0.01, label_obj_28
color black, label_obj_28
label label_obj_28, "Site 28 (Score: 6.8)"

# Site 28 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.84
# Size: 12
# Detection Methods: geometric

# Binding site 29
select site_29_center, (main_obj and chain X and resi 29 and name O)
select site_29_points, (main_obj and chain X and resi 29 and name H)
show spheres, site_29_center
color black, site_29_center
set sphere_scale, 1.0, site_29_center
show spheres, site_29_points
color salmon, site_29_points
set sphere_scale, 0.6, site_29_points
select site_29_res, ((main_obj and chain A and resi 2) or (main_obj and chain A and resi 3) or (main_obj and chain A and resi 4) or (main_obj and chain A and resi 5) or (main_obj and chain A and resi 6) or (main_obj and chain A and resi 7) or (main_obj and chain A and resi 8) or (main_obj and chain A and resi 14) or (main_obj and chain A and resi 93) or (main_obj and chain A and resi 94) or (main_obj and chain A and resi 95) or (main_obj and chain A and resi 96) or (main_obj and chain A and resi 97) or (main_obj and chain A and resi 99) or (main_obj and chain A and resi 100) or (main_obj and chain A and resi 101) or (main_obj and chain A and resi 103) or (main_obj and chain A and resi 104) or (main_obj and chain A and resi 107) or (main_obj and chain A and resi 109) or (main_obj and chain A and resi 110) or (main_obj and chain A and resi 111) or (main_obj and chain A and resi 112) or (main_obj and chain A and resi 125) or (main_obj and chain A and resi 126) or (main_obj and chain A and resi 128) or (main_obj and chain A and resi 133))
show sticks, site_29_res
color salmon, site_29_res
create site_29_surface, site_29_res
show surface, site_29_surface
set surface_color, salmon, site_29_surface
set transparency, 0.3, site_29_surface

# Add label sphere at binding site center
pseudoatom label_obj_29, pos=[22.871999740600586, 51.90433247884115, 15.879666964213053]
set sphere_scale, 0.01, label_obj_29
color black, label_obj_29
label label_obj_29, "Site 29 (Score: 6.8)"

# Site 29 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.81
# Size: 12
# Detection Methods: geometric

# Binding site 30
select site_30_center, (main_obj and chain X and resi 30 and name O)
select site_30_points, (main_obj and chain X and resi 30 and name H)
show spheres, site_30_center
color black, site_30_center
set sphere_scale, 1.0, site_30_center
show spheres, site_30_points
color salmon, site_30_points
set sphere_scale, 0.6, site_30_points
select site_30_res, ((main_obj and chain B and resi 9) or (main_obj and chain B and resi 10) or (main_obj and chain B and resi 11) or (main_obj and chain B and resi 12) or (main_obj and chain B and resi 13) or (main_obj and chain B and resi 117) or (main_obj and chain B and resi 118) or (main_obj and chain B and resi 119) or (main_obj and chain B and resi 120) or (main_obj and chain B and resi 121) or (main_obj and chain B and resi 124))
show sticks, site_30_res
color salmon, site_30_res
create site_30_surface, site_30_res
show surface, site_30_surface
set surface_color, salmon, site_30_surface
set transparency, 0.3, site_30_surface

# Add label sphere at binding site center
pseudoatom label_obj_30, pos=[20.671999740600587, 49.920999145507814, 52.363000297546385]
set sphere_scale, 0.01, label_obj_30
color black, label_obj_30
label label_obj_30, "Site 30 (Score: 6.7)"

# Site 30 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.66
# Size: 5
# Detection Methods: geometric

# Binding site 31
select site_31_center, (main_obj and chain X and resi 31 and name O)
select site_31_points, (main_obj and chain X and resi 31 and name H)
show spheres, site_31_center
color black, site_31_center
set sphere_scale, 1.0, site_31_center
show spheres, site_31_points
color salmon, site_31_points
set sphere_scale, 0.6, site_31_points
select site_31_res, ((main_obj and chain B and resi 110) or (main_obj and chain B and resi 128) or (main_obj and chain B and resi 129) or (main_obj and chain B and resi 130) or (main_obj and chain B and resi 131) or (main_obj and chain B and resi 132) or (main_obj and chain B and resi 133) or (main_obj and chain B and resi 134) or (main_obj and chain B and resi 135) or (main_obj and chain B and resi 156) or (main_obj and chain B and resi 157) or (main_obj and chain B and resi 158) or (main_obj and chain B and resi 159))
show sticks, site_31_res
color salmon, site_31_res
create site_31_surface, site_31_res
show surface, site_31_surface
set surface_color, salmon, site_31_surface
set transparency, 0.3, site_31_surface

# Add label sphere at binding site center
pseudoatom label_obj_31, pos=[9.871999740600586, 60.32099914550781, 64.36300029754639]
set sphere_scale, 0.01, label_obj_31
color black, label_obj_31
label label_obj_31, "Site 31 (Score: 6.6)"

# Site 31 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.57
# Size: 5
# Detection Methods: geometric

# Binding site 32
select site_32_center, (main_obj and chain X and resi 32 and name O)
select site_32_points, (main_obj and chain X and resi 32 and name H)
show spheres, site_32_center
color black, site_32_center
set sphere_scale, 1.0, site_32_center
show spheres, site_32_points
color salmon, site_32_points
set sphere_scale, 0.6, site_32_points
select site_32_res, ((main_obj and chain B and resi 5) or (main_obj and chain B and resi 28) or (main_obj and chain B and resi 29) or (main_obj and chain B and resi 30) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 32) or (main_obj and chain B and resi 33) or (main_obj and chain B and resi 34) or (main_obj and chain B and resi 35) or (main_obj and chain B and resi 36) or (main_obj and chain B and resi 37) or (main_obj and chain B and resi 38) or (main_obj and chain B and resi 40) or (main_obj and chain B and resi 57) or (main_obj and chain B and resi 90) or (main_obj and chain B and resi 92) or (main_obj and chain B and resi 111))
show sticks, site_32_res
color salmon, site_32_res
create site_32_surface, site_32_res
show surface, site_32_surface
set surface_color, salmon, site_32_surface
set transparency, 0.3, site_32_surface

# Add label sphere at binding site center
pseudoatom label_obj_32, pos=[10.538666407267252, 73.82099914550781, 47.29633363087972]
set sphere_scale, 0.01, label_obj_32
color black, label_obj_32
label label_obj_32, "Site 32 (Score: 6.5)"

# Site 32 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.47
# Size: 6
# Detection Methods: geometric

# Binding site 33
select site_33_center, (main_obj and chain X and resi 33 and name O)
select site_33_points, (main_obj and chain X and resi 33 and name H)
show spheres, site_33_center
color black, site_33_center
set sphere_scale, 1.0, site_33_center
show spheres, site_33_points
color salmon, site_33_points
set sphere_scale, 0.6, site_33_points
select site_33_res, ((main_obj and chain B and resi 5) or (main_obj and chain B and resi 6) or (main_obj and chain B and resi 27) or (main_obj and chain B and resi 28) or (main_obj and chain B and resi 29) or (main_obj and chain B and resi 30) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 32) or (main_obj and chain B and resi 33) or (main_obj and chain B and resi 34) or (main_obj and chain B and resi 35) or (main_obj and chain B and resi 36) or (main_obj and chain B and resi 40) or (main_obj and chain B and resi 50) or (main_obj and chain B and resi 54) or (main_obj and chain B and resi 57) or (main_obj and chain B and resi 92) or (main_obj and chain B and resi 94))
show sticks, site_33_res
color salmon, site_33_res
create site_33_surface, site_33_res
show surface, site_33_surface
set surface_color, salmon, site_33_surface
set transparency, 0.3, site_33_surface

# Add label sphere at binding site center
pseudoatom label_obj_33, pos=[16.371999740600586, 70.19599914550781, 45.71300029754639]
set sphere_scale, 0.01, label_obj_33
color black, label_obj_33
label label_obj_33, "Site 33 (Score: 6.4)"

# Site 33 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.41
# Size: 8
# Detection Methods: geometric

# Binding site 34
select site_34_center, (main_obj and chain X and resi 34 and name O)
select site_34_points, (main_obj and chain X and resi 34 and name H)
show spheres, site_34_center
color black, site_34_center
set sphere_scale, 1.0, site_34_center
show spheres, site_34_points
color salmon, site_34_points
set sphere_scale, 0.6, site_34_points
select site_34_res, ((main_obj and chain B and resi 3) or (main_obj and chain B and resi 5) or (main_obj and chain B and resi 31) or (main_obj and chain B and resi 39) or (main_obj and chain B and resi 40) or (main_obj and chain B and resi 41) or (main_obj and chain B and resi 42) or (main_obj and chain B and resi 43) or (main_obj and chain B and resi 44) or (main_obj and chain B and resi 46) or (main_obj and chain B and resi 47) or (main_obj and chain B and resi 50) or (main_obj and chain B and resi 59) or (main_obj and chain B and resi 60) or (main_obj and chain B and resi 61) or (main_obj and chain B and resi 62) or (main_obj and chain B and resi 91) or (main_obj and chain B and resi 92) or (main_obj and chain B and resi 93) or (main_obj and chain B and resi 94) or (main_obj and chain B and resi 95) or (main_obj and chain B and resi 96) or (main_obj and chain B and resi 99) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 103))
show sticks, site_34_res
color salmon, site_34_res
create site_34_surface, site_34_res
show surface, site_34_surface
set surface_color, salmon, site_34_surface
set transparency, 0.3, site_34_surface

# Add label sphere at binding site center
pseudoatom label_obj_34, pos=[22.246999740600586, 72.44599914550781, 52.58800029754639]
set sphere_scale, 0.01, label_obj_34
color black, label_obj_34
label label_obj_34, "Site 34 (Score: 6.4)"

# Site 34 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.40
# Size: 8
# Detection Methods: geometric

# Binding site 35
select site_35_center, (main_obj and chain X and resi 35 and name O)
select site_35_points, (main_obj and chain X and resi 35 and name H)
show spheres, site_35_center
color black, site_35_center
set sphere_scale, 1.0, site_35_center
show spheres, site_35_points
color salmon, site_35_points
set sphere_scale, 0.6, site_35_points
select site_35_res, ((main_obj and chain A and resi 2) or (main_obj and chain A and resi 3) or (main_obj and chain A and resi 4) or (main_obj and chain A and resi 99) or (main_obj and chain A and resi 100) or (main_obj and chain A and resi 101) or (main_obj and chain A and resi 102) or (main_obj and chain A and resi 103) or (main_obj and chain A and resi 104) or (main_obj and chain A and resi 105) or (main_obj and chain A and resi 106) or (main_obj and chain A and resi 107) or (main_obj and chain A and resi 108) or (main_obj and chain A and resi 109) or (main_obj and chain A and resi 110) or (main_obj and chain A and resi 133) or (main_obj and chain A and resi 158))
show sticks, site_35_res
color salmon, site_35_res
create site_35_surface, site_35_res
show surface, site_35_surface
set surface_color, salmon, site_35_surface
set transparency, 0.3, site_35_surface

# Add label sphere at binding site center
pseudoatom label_obj_35, pos=[22.094221962822807, 49.87655470106337, 7.963000297546387]
set sphere_scale, 0.01, label_obj_35
color black, label_obj_35
label label_obj_35, "Site 35 (Score: 6.2)"

# Site 35 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.21
# Size: 9
# Detection Methods: geometric

# Binding site 36
select site_36_center, (main_obj and chain X and resi 36 and name O)
select site_36_points, (main_obj and chain X and resi 36 and name H)
show spheres, site_36_center
color black, site_36_center
set sphere_scale, 1.0, site_36_center
show spheres, site_36_points
color salmon, site_36_points
set sphere_scale, 0.6, site_36_points
select site_36_res, ((main_obj and chain A and resi 24) or (main_obj and chain A and resi 25) or (main_obj and chain A and resi 26) or (main_obj and chain A and resi 27) or (main_obj and chain A and resi 28) or (main_obj and chain A and resi 29) or (main_obj and chain A and resi 30) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 32) or (main_obj and chain A and resi 33) or (main_obj and chain A and resi 34) or (main_obj and chain A and resi 57) or (main_obj and chain A and resi 153))
show sticks, site_36_res
color salmon, site_36_res
create site_36_surface, site_36_res
show surface, site_36_surface
set surface_color, salmon, site_36_surface
set transparency, 0.3, site_36_surface

# Add label sphere at binding site center
pseudoatom label_obj_36, pos=[25.671999740600587, 66.72099914550782, 25.763000297546387]
set sphere_scale, 0.01, label_obj_36
color black, label_obj_36
label label_obj_36, "Site 36 (Score: 6.1)"

# Site 36 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.14
# Size: 5
# Detection Methods: geometric

# Binding site 37
select site_37_center, (main_obj and chain X and resi 37 and name O)
select site_37_points, (main_obj and chain X and resi 37 and name H)
show spheres, site_37_center
color black, site_37_center
set sphere_scale, 1.0, site_37_center
show spheres, site_37_points
color salmon, site_37_points
set sphere_scale, 0.6, site_37_points
select site_37_res, ((main_obj and chain B and resi 2) or (main_obj and chain B and resi 3) or (main_obj and chain B and resi 4) or (main_obj and chain B and resi 5) or (main_obj and chain B and resi 93) or (main_obj and chain B and resi 99) or (main_obj and chain B and resi 100) or (main_obj and chain B and resi 101) or (main_obj and chain B and resi 102) or (main_obj and chain B and resi 103) or (main_obj and chain B and resi 104) or (main_obj and chain B and resi 105) or (main_obj and chain B and resi 106) or (main_obj and chain B and resi 107) or (main_obj and chain B and resi 108) or (main_obj and chain B and resi 109) or (main_obj and chain B and resi 110) or (main_obj and chain B and resi 125) or (main_obj and chain B and resi 126) or (main_obj and chain B and resi 133) or (main_obj and chain B and resi 158))
show sticks, site_37_res
color salmon, site_37_res
create site_37_surface, site_37_res
show surface, site_37_surface
set surface_color, salmon, site_37_surface
set transparency, 0.3, site_37_surface

# Add label sphere at binding site center
pseudoatom label_obj_37, pos=[18.371999740600586, 68.57099914550781, 60.71300029754639]
set sphere_scale, 0.01, label_obj_37
color black, label_obj_37
label label_obj_37, "Site 37 (Score: 6.1)"

# Site 37 details:
# Consensus Score: 1.30
# Binding Potential Score: 6.14
# Size: 8
# Detection Methods: geometric

# Binding site 38
select site_38_center, (main_obj and chain X and resi 38 and name O)
select site_38_points, (main_obj and chain X and resi 38 and name H)
show spheres, site_38_center
color black, site_38_center
set sphere_scale, 1.0, site_38_center
show spheres, site_38_points
color salmon, site_38_points
set sphere_scale, 0.6, site_38_points
select site_38_res, ((main_obj and chain A and resi 60) or (main_obj and chain A and resi 61) or (main_obj and chain A and resi 62) or (main_obj and chain A and resi 63) or (main_obj and chain A and resi 73) or (main_obj and chain A and resi 74) or (main_obj and chain A and resi 75) or (main_obj and chain A and resi 76) or (main_obj and chain A and resi 77) or (main_obj and chain A and resi 78) or (main_obj and chain A and resi 79) or (main_obj and chain A and resi 80) or (main_obj and chain A and resi 81) or (main_obj and chain A and resi 82) or (main_obj and chain A and resi 83) or (main_obj and chain A and resi 84) or (main_obj and chain A and resi 103))
show sticks, site_38_res
color salmon, site_38_res
create site_38_surface, site_38_res
show surface, site_38_surface
set surface_color, salmon, site_38_surface
set transparency, 0.3, site_38_surface

# Add label sphere at binding site center
pseudoatom label_obj_38, pos=[8.20533307393392, 55.98766581217448, 6.79633363087972]
set sphere_scale, 0.01, label_obj_38
color black, label_obj_38
label label_obj_38, "Site 38 (Score: 5.9)"

# Site 38 details:
# Consensus Score: 1.30
# Binding Potential Score: 5.92
# Size: 6
# Detection Methods: geometric

# Binding site 39
select site_39_center, (main_obj and chain X and resi 39 and name O)
select site_39_points, (main_obj and chain X and resi 39 and name H)
show spheres, site_39_center
color black, site_39_center
set sphere_scale, 1.0, site_39_center
show spheres, site_39_points
color salmon, site_39_points
set sphere_scale, 0.6, site_39_points
select site_39_res, ((main_obj and chain A and resi 42) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 44) or (main_obj and chain A and resi 60) or (main_obj and chain A and resi 61) or (main_obj and chain A and resi 62) or (main_obj and chain A and resi 63) or (main_obj and chain A and resi 64) or (main_obj and chain A and resi 65) or (main_obj and chain A and resi 66) or (main_obj and chain A and resi 67) or (main_obj and chain A and resi 73) or (main_obj and chain A and resi 74) or (main_obj and chain A and resi 75) or (main_obj and chain A and resi 76) or (main_obj and chain A and resi 77) or (main_obj and chain A and resi 78) or (main_obj and chain A and resi 80))
show sticks, site_39_res
color salmon, site_39_res
create site_39_surface, site_39_res
show surface, site_39_surface
set surface_color, salmon, site_39_surface
set transparency, 0.3, site_39_surface

# Add label sphere at binding site center
pseudoatom label_obj_39, pos=[5.871999740600586, 55.57099914550781, 12.713000297546387]
set sphere_scale, 0.01, label_obj_39
color black, label_obj_39
label label_obj_39, "Site 39 (Score: 5.3)"

# Site 39 details:
# Consensus Score: 1.00
# Binding Potential Score: 5.25
# Size: 12
# Detection Methods: geometric

# Binding site 40
select site_40_center, (main_obj and chain X and resi 40 and name O)
select site_40_points, (main_obj and chain X and resi 40 and name H)
show spheres, site_40_center
color black, site_40_center
set sphere_scale, 1.0, site_40_center
show spheres, site_40_points
color salmon, site_40_points
set sphere_scale, 0.6, site_40_points
select site_40_res, ((main_obj and chain A and resi 5) or (main_obj and chain A and resi 6) or (main_obj and chain A and resi 7) or (main_obj and chain A and resi 27) or (main_obj and chain A and resi 31) or (main_obj and chain A and resi 35) or (main_obj and chain A and resi 40) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 46) or (main_obj and chain A and resi 50) or (main_obj and chain A and resi 54) or (main_obj and chain A and resi 93) or (main_obj and chain A and resi 94) or (main_obj and chain A and resi 95) or (main_obj and chain A and resi 96) or (main_obj and chain A and resi 100))
show sticks, site_40_res
color salmon, site_40_res
create site_40_surface, site_40_res
show surface, site_40_surface
set surface_color, salmon, site_40_surface
set transparency, 0.3, site_40_surface

# Add label sphere at binding site center
pseudoatom label_obj_40, pos=[18.814000129699707, 57.48733266194662, 21.29633363087972]
set sphere_scale, 0.01, label_obj_40
color black, label_obj_40
label label_obj_40, "Site 40 (Score: 4.8)"

# Site 40 details:
# Consensus Score: 1.00
# Binding Potential Score: 4.82
# Size: 3
# Detection Methods: energy

# Binding site 41
select site_41_center, (main_obj and chain X and resi 41 and name O)
select site_41_points, (main_obj and chain X and resi 41 and name H)
show spheres, site_41_center
color black, site_41_center
set sphere_scale, 1.0, site_41_center
show spheres, site_41_points
color salmon, site_41_points
set sphere_scale, 0.6, site_41_points
select site_41_res, ((main_obj and chain A and resi 2) or (main_obj and chain A and resi 4) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 62) or (main_obj and chain A and resi 77) or (main_obj and chain A and resi 78) or (main_obj and chain A and resi 79) or (main_obj and chain A and resi 80) or (main_obj and chain A and resi 81) or (main_obj and chain A and resi 82) or (main_obj and chain A and resi 83) or (main_obj and chain A and resi 85) or (main_obj and chain A and resi 91) or (main_obj and chain A and resi 93) or (main_obj and chain A and resi 99) or (main_obj and chain A and resi 102) or (main_obj and chain A and resi 103) or (main_obj and chain A and resi 104) or (main_obj and chain A and resi 106) or (main_obj and chain A and resi 107))
show sticks, site_41_res
color salmon, site_41_res
create site_41_surface, site_41_res
show surface, site_41_surface
set surface_color, salmon, site_41_surface
set transparency, 0.3, site_41_surface

# Add label sphere at binding site center
pseudoatom label_obj_41, pos=[16.480666796366375, 54.820665995279946, 7.29633363087972]
set sphere_scale, 0.01, label_obj_41
color black, label_obj_41
label label_obj_41, "Site 41 (Score: 4.8)"

# Site 41 details:
# Consensus Score: 1.00
# Binding Potential Score: 4.80
# Size: 6
# Detection Methods: energy

# Binding site 42
select site_42_center, (main_obj and chain X and resi 42 and name O)
select site_42_points, (main_obj and chain X and resi 42 and name H)
show spheres, site_42_center
color black, site_42_center
set sphere_scale, 1.0, site_42_center
show spheres, site_42_points
color salmon, site_42_points
set sphere_scale, 0.6, site_42_points
select site_42_res, ((main_obj and chain A and resi 38) or (main_obj and chain A and resi 39) or (main_obj and chain A and resi 40) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 47) or (main_obj and chain A and resi 52) or (main_obj and chain A and resi 53) or (main_obj and chain A and resi 54) or (main_obj and chain A and resi 55) or (main_obj and chain A and resi 56) or (main_obj and chain A and resi 57) or (main_obj and chain A and resi 58) or (main_obj and chain A and resi 59) or (main_obj and chain A and resi 60) or (main_obj and chain A and resi 61) or (main_obj and chain A and resi 71) or (main_obj and chain A and resi 72) or (main_obj and chain A and resi 73))
show sticks, site_42_res
color salmon, site_42_res
create site_42_surface, site_42_res
show surface, site_42_surface
set surface_color, salmon, site_42_surface
set transparency, 0.3, site_42_surface

# Add label sphere at binding site center
pseudoatom label_obj_42, pos=[13.70533307393392, 66.48766581217448, 14.963000297546387]
set sphere_scale, 0.01, label_obj_42
color black, label_obj_42
label label_obj_42, "Site 42 (Score: 4.4)"

# Site 42 details:
# Consensus Score: 1.00
# Binding Potential Score: 4.44
# Size: 6
# Detection Methods: geometric

# Binding site 43
select site_43_center, (main_obj and chain X and resi 43 and name O)
select site_43_points, (main_obj and chain X and resi 43 and name H)
show spheres, site_43_center
color black, site_43_center
set sphere_scale, 1.0, site_43_center
show spheres, site_43_points
color salmon, site_43_points
set sphere_scale, 0.6, site_43_points
select site_43_res, ((main_obj and chain A and resi 40) or (main_obj and chain A and resi 41) or (main_obj and chain A and resi 42) or (main_obj and chain A and resi 43) or (main_obj and chain A and resi 58) or (main_obj and chain A and resi 59) or (main_obj and chain A and resi 60) or (main_obj and chain A and resi 61) or (main_obj and chain A and resi 62) or (main_obj and chain A and resi 63) or (main_obj and chain A and resi 72) or (main_obj and chain A and resi 73) or (main_obj and chain A and resi 74) or (main_obj and chain A and resi 75) or (main_obj and chain A and resi 76) or (main_obj and chain A and resi 77) or (main_obj and chain A and resi 80) or (main_obj and chain A and resi 81) or (main_obj and chain A and resi 84))
show sticks, site_43_res
color salmon, site_43_res
create site_43_surface, site_43_res
show surface, site_43_surface
set surface_color, salmon, site_43_surface
set transparency, 0.3, site_43_surface

# Add label sphere at binding site center
pseudoatom label_obj_43, pos=[9.071999740600585, 60.12099914550781, 11.963000297546387]
set sphere_scale, 0.01, label_obj_43
color black, label_obj_43
label label_obj_43, "Site 43 (Score: 4.3)"

# Site 43 details:
# Consensus Score: 1.00
# Binding Potential Score: 4.34
# Size: 5
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
save notebook_results_predicted_binding_sites.pse
