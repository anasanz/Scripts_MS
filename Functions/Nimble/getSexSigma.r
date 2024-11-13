getSexSigma <- nimbleFunction(
  run=function(
    age.cat         = double(0),
    sex             = double(0)
  ){
    returnType(double(0))
    if( (age.cat == 1) | (age.cat == 2) ) {sex.age.idx <-0 } else {
      sex.age.idx <- sex}
    return(sex.age.idx)
  })
