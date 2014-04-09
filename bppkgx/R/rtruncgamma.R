rtruncgamma.gibbs <- function( n, a, b, shape, scale, 
                         x0 = runif( 1, a, b ), y0 = runif( 1, 0, exp( -x0 * scale ) ) ){
    y    <- rep( NA, n )
    x    <- rep( NA, n )
    y[1] <- y0
    x[1] <- x0
    for( i in 2:n ){
        y[i] <- runif( 1, 0, exp( -x[i-1] * scale ) )
        x[i] <- ( shape * runif( 1, a^shape /shape, min( b, -log( y[i] )/scale )^shape/shape ) )^( 1/shape )
    }
    return( x )
}

# rtruncbeta.gibbs <- function( n, a, b, shape1, shape2, 
#                          x0 = runif( 1, a, b ), y0 = runif( 1, 0, (1-x0)^(shape2-1) ) ){
#     y    <- rep( NA, n )
#     x    <- rep( NA, n )
#     y[1] <- y0
#     x[1] <- x0
#     for( i in 2:n ){
#         y[i] <- runif( 1, 0, ( 1 - x[i-1])^(shape2-1) )
#         if(shape2==1)
#         #x[i] <- ( shape * runif( 1, a^shape /shape, min( b, -log( y[i] )/scale )^shape/shape ) )^( 1/shape )
#     }
#     return( x )
# }

rtruncnorm <- function( n, lower = -Inf, upper = Inf, mean = 0, sd = 1 ){

    return( qnorm( runif( n, pnorm( lower, mean = mean, sd = sd ),
                             pnorm( upper, mean = mean, sd = sd  ) ), 
                    mean = mean, sd = sd ) )
}


# Sample from truncated normal N(mu,sig^2)1_{lower<r<upper}
rtruncnorm <- function( mu, sig, lower, upper){

    x1 =  ( lower-mu ) / sig ;
    x2 =  ( upper-mu ) / sig ;
    p1 = pnorm( x1, 0.0, 1.0 );
    p2 = pnorm( x2, 0.0, 1.0 );
    if( p2 < 1e-30 ){
        x = upper;
    } else {
        u = ( p2 - p1 ) * runif( 1, 0, 1 ) + p1;
        u = min( max( 1e-30, u ), 1-1e-16 );
        x4 = qnorm( u, 0.0, 1.0 );
        x = x4 * sig + mu;
    } 
    return(x);
}
rtgamma <- function( n, lower = 0, upper = Inf, shape, scale ){
    return( qgamma( runif( n, pgamma( lower, shape = shape, scale = scale ),
                              pgamma( upper, shape = shape, scale = scale ) ), 
                    shape = shape, scale = scale ) )
}

rtbeta <- function( n, lower = 0, upper = 1, shape1, shape2, ncp = 0 ){
    return( qbeta( runif( n, pbeta( lower, shape1=shape1, shape2=shape2, ncp = ncp ),
                              pbeta( upper, shape1=shape1, shape2=shape2, ncp = ncp ) ), 
                    shape1=shape1, shape2=shape2, ncp = ncp ) )
}

#test
test=0
if(test==1){
    x=rtruncgamma(10000,1,3,2,2)
    hist(x)
    hist( rtgamma(10000,1,3,2,2))
}