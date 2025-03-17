border_check = function(lscape_defs, xmax, ymax) {

  border_defs <- tibble::tibble(
    X = c(list(1:xmax), xmax, list(1:xmax), 1),
    Y = c(1, list(1:ymax), ymax, list(1:ymax)),
    directions = c(3 * pi / 2, 0, pi / 2, pi)
  )
  
  border_errors <- c()
  for (ii in 1:4) {
    lscape_border <- lscape_defs %>% 
      dplyr::filter(
        X %in% border_defs$X[[ii]],
        Y %in% border_defs$Y[[ii]],
        Road == "On Trail"
      )
    
    if (any(abs(lscape_border$Direction %% (2 * pi) - 
                border_defs$directions[ii] %% (2 * pi)) < 0.001)) {
      border_errors <- border_errors %>% 
        dplyr::bind_rows(
          lscape_border %>% 
            dplyr::filter(abs(Direction %% (2 * pi) - 
                                border_defs$directions[ii] %% (2 * pi)) < 0.001)
          )
    }
  }
  
  return(border_errors)
}
