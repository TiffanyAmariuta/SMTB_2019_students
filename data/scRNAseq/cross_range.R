cross = function(vec1, vec2, characters=1000) {
  out = c()
  for (i in 1:length(vec1)) {
    out = c(out, vec1[i], vec2[i])
  }
  
  return(unique(out)[1:characters])
}

cross(c(1, 2, 3), c(1, 20, 30), characters=3)
  