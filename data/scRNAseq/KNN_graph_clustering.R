a = readRDS("~/R/SMTB_2019_students/data/scRNAseq/shrinked_table.rds")

count_dists = function(expression) {
  count_dist = function(another_expression) {
    return(sqrt(sum((expression - another_expression)^2)))
  }
  return(apply(a, 2, count_dist))
}

dist_matrix = dist(a)
