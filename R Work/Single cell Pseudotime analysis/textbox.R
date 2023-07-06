# Empty vectors for text box coords
x_L_coord_list <- vector("numeric", length = 4)
x_R_coord_list <- vector("numeric", length = 4)
y_B_coord_list <- vector("numeric", length = 4)
y_T_coord_list <- vector("numeric", length = 4)

# For loop calculating the box coords based on box center coords
for (clus in 1:4){
  # define coordinates for rectangle
  x_center = c(-0.015,0.025,0.02,-0.015)[clus]
  y_center = c(-0.015,-0.02,0.02,0.03)[clus]
  
  # X coords
  x_len = str_length(cluster_names[clus])/650
  x_L_coord = x_center - (x_len*0.5)
  x_R_coord = x_center + (x_len*0.5)
  
  # Y coords
  y_len = str_width(cluster_names[clus])/1400
  y_B_coord = y_center - (y_len*0.5)
  y_T_coord = y_center + (y_len*0.5)
  
  # appending calculated coords to list
  x_L_coord_list[clus] <- x_L_coord
  x_R_coord_list[clus] <- x_R_coord
  y_B_coord_list[clus] <- y_B_coord
  y_T_coord_list[clus] <- y_T_coord
}