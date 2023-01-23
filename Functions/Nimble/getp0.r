
getp0 <- nimbleFunction(
  run = function(
    p0 = double(1), country = double(1)) {
    nsites <- length(country)
    out <- country
    
    for(i in 1:nsites){
      out[i] <- p0[country[i]]
    }
    return(out)
    returnType(double(1))
  })