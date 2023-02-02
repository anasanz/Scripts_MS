buffer_rectangle <- function(point, x_length, y_length, degree = 0){
  
  #Error handlers for input type
  if(!("sf" %in% class(point))) stop("Input one point of class sf")
  if(nrow(point) != 1) stop("Input one point of class sf")
  if(!sf::st_is(point, "POINT")) stop("Input one point of class sf")
  
  if(!is.numeric(degree) || degree < 0 || degree > 360){
    stop("Input degree as numeric between 0 and 360")
  }
  
  #function starts
  point_coordinates <- sf::st_coordinates(point)
  
  radians <- degree * 0.0174532925
  
  #create an empty matrix
  pts_df = as.data.frame(matrix(nrow = 4, ncol = 2))
  colnames(pts_df) <- c("x", "y")
  
  #set the points
  pts_df[1,1] <- (point_coordinates[1] - x_length / 2) 
  pts_df[1,2] <- (point_coordinates[2] + y_length / 2)
  
  pts_df[2,1] <- (point_coordinates[1] + x_length / 2) 
  pts_df[2,2] <- (point_coordinates[2] + y_length / 2)
  
  pts_df[3,1] <- (point_coordinates[1] - x_length / 2) 
  pts_df[3,2] <- (point_coordinates[2] - y_length / 2) 
  
  pts_df[4,1] <- (point_coordinates[1] + x_length / 2) 
  pts_df[4,2] <- (point_coordinates[2] - y_length / 2) 
  
  #convert to sf
  pts_sf <- sf::st_as_sf(pts_df, coords = c("x", "y"), crs = sf::st_crs(point))
  
  ##create the convex hull
  rectangular_sf <- sf::st_convex_hull(sf::st_union(pts_sf))
  
  (rectangular_sf - point_coordinates) * rotation_f(radians) + point_coordinates -> rectangular_sf
  
  #return
  return(rectangular_sf)  
}

buffer_square <- function(point, length, degree = 0){
  return(buffer_rectangle(point, length, length, degree))  
}



mcp <- function(x, percentile=95){
  
  centroid <- sf::st_centroid(sf::st_union(x))
  dist <- as.numeric(sf::st_distance(x, centroid))
  within_percentile_range <- dist <= quantile(dist, percentile/100)
  x_filter <- st_union(x[within_percentile_range,])
  st_convex_hull(x_filter)
  
}
