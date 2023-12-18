# The "match" variable identifies the SOILWAT2 layer above (col 1) and 
# below (col 2) the associated observed soil depth (rows). If the match
# variable is the same (col 1 = col 2 value) for a particular row (obs
# layer depth), then the observed depth occurs "inside" the corresponding
# SOILWAT2 layer.

match.ses = matrix(data = c(1,1,
                        3,3,
                        5,5,
                        7,7,
                        9,9), nrow=5, ncol=2,byrow=TRUE)

match.seg = matrix(data = c(1,1,
                            3,3,
                            5,5,
                            7,7,
                            9,9), nrow=5, ncol=2,byrow=TRUE)

match.wjs = matrix(data = c(1,2,
                            2,3,
                            6,7), nrow=3, ncol=2,byrow=TRUE)

match.mpj = matrix(data = c(1,2,
                            2,3,
                            6,7), nrow=3, ncol=2,byrow=TRUE)

match.vcp = matrix(data = c(1, 2,
                            2, 3,
                            4, 5,
                            6, 7,
                            8, 9,
                            9, 9), nrow=6, ncol=2,byrow=TRUE)

match.vcm = matrix(data = c(1, 2,
                            2, 3,
                            4, 5,
                            6, 7,
                            8, 9,
                            9, 9), nrow=6, ncol=2,byrow=TRUE)

match.vcs = matrix(data = c(1, 2,
                            2, 3,
                            6, 7,
                            9, 9), nrow=4, ncol=2,byrow=TRUE)