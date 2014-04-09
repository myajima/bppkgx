# vector will be row or column vector so the dimension is at least 2 
read.arma.ascii <- function( filename ){
    dm <- as.integer( read.table( filename, skip = 1, nrows = 1 ) ) # dim
    dt <- as.matrix( read.table( filename, skip = 2 ) ) # data
    if( length( dm ) == 2 ){
        rs = dt
    } else if ( length( dm ) == 3 ){
        rs <- array( NA, dm )
        for( i in 1:dm[3] ) {
            rs[,,i] <- dt[ ( 1:dm[1] + dm[1] * ( i - 1 ) ),]
        }
    }
    return( rs )
}


if(0){
ddd =read.arma.ascii("PostPkPg/Y.txt")
}