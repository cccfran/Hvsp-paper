### find a beta that gives a predefined out-in ratio in balanced binary string SBM
FindBeta <- function(DD,rho){
    test.func <- function(beta,DD,rho){
        if(beta==0.5) return(1)
        print((beta*((2*beta)^DD-1))/((2*beta-1)*(2^DD-1)) - rho)
        return((beta*((2*beta)^DD-1))/((2*beta-1)*(2^DD-1)) - rho)
    }
    beta.left <- 0
    beta.right <- 0.8
    f1 <- test.func(beta.left,DD,rho)
    f2 <- test.func(beta.right,DD,rho)

    for(k in 1:15){
        beta.new <- 0.5*beta.left+0.5*beta.right
        f.new <- test.func(beta.new,DD,rho)
        if(f1*f.new < 0){
            f2 <- f.new
            beta.right <- beta.new
        }else{
            f1 <- f.new
            beta.left <- beta.new
        }
    }
    return(beta.new)
}
