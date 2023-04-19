#' Create Trapezoid Shape
#'
#' This function creates a trapezoid shape with given coordinates.
#'
#' @param xPos A numeric vector of length 4 representing x-coordinates of the
#' trapezoid in order.
#' @param yPos A numeric vector of length 4 representing y-coordinates of the
#' trapezoid in order.
#' @param nDivisions An integer value indicating the number of divisions to
#' create the trapezoid shape. Default is 2000.
#' @param reverse A boolean value indicating whether to create the trapezoid in
#' reverse order or not. Default is FALSE.
#'
#' @return A data frame with id, x and y coordinates.
#'
#' @importFrom plyr ldply
#' @import tidyr
#' @export
createTrapezoid <- function(xPos = NULL,
                            yPos = NULL,
                            nDivisions = 2000,
                            reverse = FALSE){
  n = nDivisions
  # j = 1
  plyr::ldply(1:n,function(j){
    # whether reverse Trapezoid graph
    if(reverse == TRUE){
      x1 = xPos[2];x2 = xPos[1];x3 = xPos[4];x4 = xPos[3]
    }else{
      x1 = xPos[1];x2 = xPos[2];x3 = xPos[3];x4 = xPos[4]
    }
    y1 = yPos[1];y2 = yPos[2];y3 = yPos[3];y4 = yPos[4]

    # calculate positions
    xp <- c(x1 + (x2-x1)/n*(j-1),x1 + (x2-x1)/n*j,x4 - (x4-x3)/n*j,x4 - (x4-x3)/n*(j-1))
    yp <- c(y1 + (y2-y1)/n*(j-1),y1 + (y2-y1)/n*j,y1 + (y2-y1)/n*j,y4 + (y3-y4)/n*(j-1))

    # add id
    res <- data.frame(id = j,x = xp,y = yp)
    return(res)
  }) -> trapezoid.df
  return(trapezoid.df)
}


#' Create segment coordinates
#'
#' This function takes starting and ending coordinates and creates segments in
#' equal divisions or given relative length. It returns a data frame with the x
#' and y coordinates of each segment.
#'
#' @param xPos A vector of length 2 containing the x-coordinates of the starting
#' and ending positions for the segment.
#' @param yPos A vector of length 2 containing the y-coordinates of the starting
#' and ending positions for the segment.
#' @param n_division An integer value indicating the number of segments the line
#' should be divided into.
#' @param rel_len A decimal value indicating the relative length of each segment.
#' @param abs_len A absolute length for generating segments.
#'
#' @return A data frame with the x and y-coordinates of each segment.
#'
#' @examples
#' createSegment(xPos = c(0.5, 0.8), yPos = c(0.2, 0.5), n_division = 4)
#' createSegment(xPos = c(0.5, 0.8), yPos = c(0.2, 0.5), rel_len = 0.2)
#' createSegment(xPos = c(0.5, 0.8), yPos = c(0.2, 0.5), abs_len = 0.1)
#'
#' @export
createSegment <- function(xPos = NULL,
                          yPos = NULL,
                          n_division = NULL,
                          rel_len = NULL,
                          abs_len = NULL){
  if(xPos[1] > xPos[2]){
    x = xPos[2];xend = xPos[1]
  }else{
    x = xPos[1];xend = xPos[2]
  }

  if(yPos[1] > yPos[2]){
    y = yPos[2];yend = yPos[1]
  }else{
    y = yPos[1];yend = yPos[2]
  }

  # get seqs
  if(y == yend){
    seg_len = abs(x - xend)
    if(!is.null(n_division)){
      seqs <- seq(from = x,to = xend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*seg_len
      }else{
        by = abs_len
      }
      seqs <- seq(from = x,to = xend,by = by)
      if(!(xend %in% seqs)) seqs <- c(seqs,xend)
    }
    res <- data.frame(id = 1:(length(seqs) - 1),
                      x = seqs[1:(length(seqs) - 1)],
                      xend = seqs[2:length(seqs)],
                      y = y,
                      yend = y)
  }else if(x == xend){
    seg_len = abs(y - yend)
    if(!is.null(n_division)){
      seqs <- seq(from = y,to = yend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*seg_len
      }else{
        by = abs_len
      }
      seqs <- seq(from = y,to = yend,by = by)
      if(!(yend %in% seqs)) seqs <- c(seqs,yend)
    }
    res <- data.frame(id = 1:(length(seqs) - 1),
                      x = x,
                      xend = x,
                      y = seqs[1:(length(seqs) - 1)],
                      yend = seqs[2:length(seqs)])
  }else{
    xseg_len = abs(xend - x)
    yseg_len = abs(yend - y)
    if(!is.null(n_division)){
      xseqs <- seq(from = x,to = xend,length.out = n_division)
      yseqs <- seq(from = y,to = yend,length.out = n_division)
    }else{
      # check type
      if(is.null(abs_len)){
        by = rel_len*xseg_len
      }else{
        by = abs_len
      }
      xseqs <- seq(from = x,to = xend,by = by)
      if(!(xend %in% xseqs)) xseqs <- c(xseqs,xend)
      yseqs <- seq(from = y,to = yend,length.out = length(xseqs))
    }
    res <- data.frame(id = 1:(length(xseqs) - 1),
                      x = xseqs[1:(length(xseqs) - 1)],
                      xend = xseqs[2:length(xseqs)],
                      y = yseqs[1:(length(yseqs) - 1)],
                      yend = yseqs[2:length(yseqs)])
  }
  return(res)
}


#' Create a Pair of Polygons
#'
#' This function creates a pair of polygons, each consisting of two line segments,
#' with a shared mid-point and a user-defined start and end point.
#'
#' @param start A numeric value indicating the start point of the polygon.
#' Defaults to NULL.
#' @param end A numeric value indicating the end point of the polygon. Defaults
#' to NULL.
#' @param open A numeric value indicating the openness of the polygon. Defaults
#' to 1.
#' @param npoints An integer value indicating the number of points to be used to
#' construct the polygon. Defaults to 1000.
#' @param reverse A logical value indicating whether to reverse the order of the
#' polygons.
#'
#' @return A data.frame object containing the x and y coordinates of the two
#' polygons and the type of polygon.
#'
#' @examples
#' createPairPolygon(start = 0, end = 10, open = 1, npoints = 1000, reverse = FALSE)
#'
#' @export
createPairPolygon <- function(start = NULL,
                              end = NULL,
                              open = 1,
                              npoints = 1000,
                              reverse = FALSE){
  mid = (start + end)/2
  interval_1 = seq(start,mid,length.out = npoints)

  # check start+end = 0
  if(mid == 0){
    y_interval = abs(seq(0,sum(abs(start) + abs(end))/2,
                         length.out = npoints))
  }else{
    y_interval = seq(0,mid,length.out = npoints)
  }

  ratio = sqrt(abs(max(y_interval))/2)
  x = c(interval_1,rev(interval_1))
  y = c(-(sqrt(abs(rev(y_interval)/2))*open),
        sqrt(abs(y_interval)/2)*open)/ratio


  interval_2 = seq(mid,end,length.out = npoints)
  x1 = c(rev(interval_2),interval_2)
  y1 = y

  if(reverse == TRUE){
    ypos = c(x1,x);xpos = c(y,y1)
  }else{
    xpos = c(x1,x);ypos = c(y,y1)
  }

  res <- data.frame(xpos = xpos,ypos = ypos,
                    type = rep(c(1,2),each = length(x)))
}



#' Add a custom annotation layer to a ggplot
#'
#' This function adds a custom annotation layer to a ggplot using the provided grob object.
#' The annotation can be placed within specified x and y limits.
#'
#' @param grob A grob object that will be used as the custom annotation.
#' @param xmin The minimum x value where the custom annotation will be placed (default is -Inf).
#' @param xmax The maximum x value where the custom annotation will be placed (default is Inf).
#' @param ymin The minimum y value where the custom annotation will be placed (default is -Inf).
#' @param ymax The maximum y value where the custom annotation will be placed (default is Inf).
#' @param data The data that the custom annotation layer will be applied to.
#' @return A ggplot layer with the custom annotation.
#' @export
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data){
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob,
                                          xmin = xmin, xmax = xmax,
                                          ymin = ymin, ymax = ymax))
}


#' getRotatedPolygon
#'
#' This function rotates a polygon based on a given angle and returns the rotated
#' coordinates along with a polygon of a given window size around the rotated polygon.
#'
#' @param data A data frame containing the original polygon coordinates.
#' @param rx The column name of the x-coordinates in the \code{data} parameter.
#' @param ry The column name of the y-coordinates in the \code{data} parameter.
#' @param value The column name of the values associated with each point in the polygon.
#' @param theta The angle in degrees to rotate the polygon.
#' @param workers The number of cores to use for parallel processing. Default is 1.
#' @param window The size of the window around the rotated polygon. Can be a numeric
#' value or the name of a column in the \code{data} parameter. Default is 1.
#'
#' @return A list containing the rotated coordinates and the polygon coordinates.
#'
#' @export
getRotatedPolygon <- function(data = NULL, rx = NULL, ry = NULL,
                              value = NULL,theta = 45, workers = 1,
                              window = 1) {
  # Convert theta to radians
  theta_rad <- pi * (theta / 180)

  # Set a "plan" for how the code should run.
  # future::plan(future::multisession, workers = workers)

  # Vectorize the operations
  data$x <- data[[rx]]
  data$y <- data[[ry]]
  data$xr <- data$x * cos(theta_rad) + data$y * sin(theta_rad)
  data$yr <- data$y * cos(theta_rad) - data$x * sin(theta_rad)
  data$xr <- data$xr * cos(theta_rad)
  data$yr <- data$yr * sin(theta_rad)

  data$yr <- round(data$yr,digits = 1)
  # Combine the results
  rotated_coods <- furrr::future_map_dfc(data, identity)

  # ============================================================================
  # check window size
  if(is.numeric(window)){
    window_size = window*0.5
  }else if(is.character(window)){
    window_size = rotated_coods[[window]]*0.5
  }else{
    message("please supply column name or numeric value!")
  }

  # get cordinates
  polygon_x <- cbind(rotated_coods$xr - window_size, rotated_coods$xr,
                     rotated_coods$xr + window_size, rotated_coods$xr)
  polygon_y <- cbind(rotated_coods$yr, rotated_coods$yr + window_size,
                     rotated_coods$yr, rotated_coods$yr - window_size)
  polygon_id <- rep(1:nrow(rotated_coods), 4)
  polygon_value <- rep(rotated_coods[[value]], 4)

  # polygon coords
  polygon_coods <- data.frame(xp = as.vector(polygon_x),
                              yp = as.vector(polygon_y),
                              id = polygon_id,
                              value = polygon_value)

  # ============================================================================
  # output
  return(list(rotated_coods = rotated_coods,
              polygon_coods = polygon_coods))
}


#' Create Triangle Coordinates
#'
#' This function creates a data frame of x and y coordinates for a triangle
#' defined by three vertices.
#' The user can specify the number of divisions and which vertex to use as the base.
#'
#' @param x numeric vector of length 3, containing the x coordinates of the
#' triangle vertices.
#' @param y numeric vector of length 3, containing the y coordinates of the
#' triangle vertices.
#' @param nDivisions integer specifying the number of divisions between each
#' vertex. Default is 1.
#' @param vertice integer specifying which vertex to use as the base. Must be 1,
#' 2 or 3. Default is 2.
#'
#' @return A data frame with columns "x", "y" and "id". The "id" column specifies
#' which division the coordinate belongs to.
#'
#' @examples
#' \dontrun{
#' createTriangle(x = c(0, 5, 3), y = c(0, 0, 4))
#' createTriangle(x = c(0, 5, 3), y = c(0, 0, 4), nDivisions = 3, vertice = 1)
#' }
#'
#' @export
createTriangle <- function(x = NULL,y = NULL,
                           nDivisions = 1,vertice = 2){
  if(vertice == 2){
    x1 = x[1];x2 = x[2];x3 = x[3]
    y1 = y[1];y2 = y[2];y3 = y[3]
  }else if(vertice == 1){
    x1 = x[3];x2 = x[1];x3 = x[2]
    y1 = y[3];y2 = y[1];y3 = y[2]
  }else if(vertice == 3){
    x1 = x[2];x2 = x[3];x3 = x[1]
    y1 = y[2];y2 = y[3];y3 = y[1]
  }else{
    message("please define 1/2/3 for vertice!")
  }

  # iter = 1
  id_num <- nDivisions:1
  plyr::ldply(1:nDivisions,function(iter){
    x_coord <- c(x2 - (iter/nDivisions)*(x2 - x1),
                 x2,
                 x2 + (iter/nDivisions)*(x3 - x2))
    y_coord <- c(y2 - (iter/nDivisions)*(y2 - y1),
                 y2,
                 y2 - (iter/nDivisions)*(y2 - y3))

    res <- data.frame(x = x_coord,y = y_coord,id = id_num[iter])
    return(res)
  }) -> coord_res

  return(coord_res)
}
