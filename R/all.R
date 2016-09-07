if(FALSE){
  setClass("tensor", representation(A="array",alongA="integer",
                                    alongB="integer"), prototype=
             list(A=array(), 
                  alongA=integer(),
                  alongB=integer()),
           validity=function(object){
             if(is.list(object)){
               if(length(object)==3){
                 if(is.array(object)[[1]]&is.integer(object[[2]])&is.integer(object[[3]])){
                   if(is.array(object[[1]])){if(identical(names(object),c("A","alongA","alongB"))){
                     identical(dim(object),c(object[[2]],object[[3]]))}
                     else{FALSE}}
                   else{FALSE}}
                 else{FALSE}}
               else{FALSE}}
             else{FALSE}
           })
  as.tensor<-function(A){list(A=array(A,dim(A)),alongA=dim(A),alongB=integer(0))}
  as.tensor(matrix(1:20,2))
}
#' Convert an array to an array of dimension 2 or 1
#' 
#' @param A An array of dimensions (a_1 x ... x a_{dim(A)})
#' @param n An integer between 0 and length(dim(A))
#' @return M the matrix of dimension (a_1 x ... x a_n) x (a_{n+1} x ... x a_{dim(A)})
#' @examples
#' A=array(rnorm(prod(2:5)),2:5);M=a2m(A,2);dim(A);dim(M);dim(a2m(A))
a2m<-function(A,n=length(dim(A))){
  p=length(dim(A))-n
  array(A,c(if(n==0){integer(0)}else{prod(dim(A)[1:n])},if(p==0){ingeger(0)}else{prod(dim(A)[(n+1):(n+p)])}))}

#' Convert  a matrix to an array 
#' 
#' @param M A matrix of dimensions (a_1 x ... x a_{dim(A)})
#' @param a A vector of integers (numeric(0) is accepted)
#' @return b A vector of integers (numeric(0) is accepted), such that prod(a)*prod(b)=prod(dim(M))
#' @examples
#' M<-matrix(1:(prod(2:5)),prod(2:3),prod(4:5));m2a(M,2:3,4:5);identical(m2a(M),M)
m2a<-function(M,a=dim(M),b=numeric(0)){
  if(prod(a)*prod(b)!=prod(dim(M))){stop("m2a: prod(a)*prod(b)!=prod(dim(M))")}
  array(M,c(a,b))}

#' Computes a matrix that is a linear combinaison of the rotation group mis estimates
#' 
#' @param X An array of dimension b_1 x ... x b_p
#' @param W An array of dimension a_1 x ... x a_n x b_1 x ... x b_p
#' @return Y the array of dimension a_1 x ... x a_n Y[i_1,...,i_n]=sum(W[i_1,...,i_n,,...,] )
#' @examples
#' W=array(1:(prod(2:5)),2:5);X=array(1:(prod(4:5)),4:5); W%.%.;try(W[,,,-1]%.%.);X%.%X; sum(X*X);X%.%t(X);sum(c(X)*c(t(X)))
#' X=array(1:(prod(4:6)),4:6); "%.%"(W,X,j=2);
#' W%.%.%.%t(X);

"%.%" <-
  function(W,X,j=NULL,k=0){
    p<-dim(X)[0:j]
    if(is.null(j)&is.null(k)){k=0;j=length(dim(X))}
    if(is.null(j)){j<-length(dim(X))-k}
    q<-dim(X)[-(0:j)]
    i<-length(dim(W))-j
    if(i<0|k<0|j<0){stop(paste0("non-conformable arguments"))}
    n<-dim(W)[0:i]
    if(prod(dim(W))!=prod(n)*prod(p)){stop(paste0("non-conformable arguments"))}
    Y<-m2a(a2m(W,i)%*%a2m(X,j),c(n,q))
    dimnames(Y)<-c(dimnames(W)[1:i],dimnames(X)[-(0:j)])
    names(dimnames(Y))<-names(dimnames(W)[1:i])
    return(Y)}


#' Computes a matrix that is a linear combinaison of the rotation group mis estimates
#' 
#' @param X An array of dimension b_1 x ... x b_p
#' @param W An array of dimension a_1 x ... x a_n x b_1 x ... x b_p
#' @return Y the array of dimension a_1 x ... x a_n Y[i_1,...,i_n]=sum(W[i_1,...,i_n,,...,] )
#' @examples
#' W=array(1:(prod(2:5)),2:5);X=array(1:(prod(4:5)),4:5); W%.%.;try(W[,,,-1]%.%.);X%.%.; sum(X*X);X%.%t(X);sum(c(X)*c(t(X)))
#' X=array(1:(prod(4:6)),4:6); "%.%"(W,X,j=2);
#' W%.%.%.%t(X);

array.prod<-function(...,j=NULL,k=0){
  "%.k%"<-function(W,X){"%.%"(W,X,j=j)}
  apply(..., 1, function(x) Reduce("%.k%", x,accumulate = FALSE))
}

"%.k1%"<-function(W,X){"%.%"(W,X,k=1)}
"%.k2%"<-function(W,X){"%.%"(W,X,k=2)}
"%.k3%"<-function(W,X){"%.%"(W,X,k=3)}