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
  array(A,c(if(n==0){integer(0)}else{prod(dim(A)[1:n])},if(p==0){integer(0)}else{prod(dim(A)[(n+1):(n+p)])}))}


A2M<-function(A,n,p=if(is.character(n)){setdiff(names(dimnames(A),n))}else{setdiff(1:length(dim(A)),n)}){
  if(is.character(n)){n<-match(n,names(dimnames(A)))}
  if(is.character(p)){p<-match(p,names(dimnames(A)))}
  array(aperm(A,c(n,p)),c(prod(dim(A)[n]),prod(dim(A)[p])))}





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


#' @examples
#' A=array(1:(prod(2:4)),2:4);
#' dimnames(A)<-sapply(dim(A),seq_len)
#' names(dimnames(A))<-paste0("x",2:4)
#' extractA(A,integer(0),integer(0));
#' extractA(A,"x3",2);
#' extractA(A,c("x4","x3"),1,2);

extractA<-function(A,a,...){
  i=lapply(dimnames(A),function(x){TRUE})
  x=list(...)
  names(x)<-NULL
  i[a]<-x
  do.call("[",c(list(A),i,drop=FALSE))}


#' Define a tensor product
#' 
#' @param A An array of dimension $\prod_{i\in I_A} a_i$
#' @param B An array of dimension $\prod_{i\in I_B} b_i$
#' @param I_A, I_B list of subvectors from names(dimnames(A)) and names(dimnames(B))
#' @return C=AB the array of dimension prod_{i\in Ia\cup Ib\cup I}a_1 x ... x a_q x b_1 x ... x b_n  Y[i_1,...,i_q,j_1,...,j_n]=sum_{k_1,...,k_p}(W[i_1,...,i_q,j_1,...,j_n,k_1,...,k_p]*X[j_1,...,i_q,k_1,...,k_p])
#' @examples
#' A=array(1:(prod(2:6)),2:6);
#' dimnames(A)<-sapply(dim(A),seq_len)
#' names(dimnames(A))<-paste0("x",2:6)
#' B=array(1:(prod(3:7)),3:7);
#' dimnames(B)<-sapply(dim(B),seq_len)
#' names(dimnames(B))<-paste0("y",3:7)
#' I_A=list(c=c("x3","x5"),n=c("x2","x4"),p="x6")
#' I_B=list(c=c("y3","y5"),q=c("y4","y7"),p="y6")
#' "%.%"(A,B)
#' "%.%"(A,B,I_A,I_B)
#' W%.%.%.%t(X);


"%.%" <-
  function(A,B,I_A=list(c=integer(0),n=1:length(dim(A)),p=integer(0)),
           I_B=list(c=integer(0),p=integer(0),q=1:length(dim(B)))){
    dime<-c(dim(A)[match(I_A$n,names(dimnames(A)))],dim(B)[match(I_B$q,names(dimnames(B)))])
    D=if(length(I_A$c)>0){expand.grid(dimnames(A)[I_A$c])}else{data.frame(row.names = 1)}
    
    C<-plyr::maply(.data=D,.fun=function(...){
      plyr::aaply(extractA(A,I_A$c,...),I_A$a,function(a){
        plyr::aaply(extractA(B,I_B$c,...),I_B$b,function(b){
          array(A2M(a,c(I_A$c,I_A$n),c(I_A$p))%*%A2M(b,c(I_B$p),c(I_B$c,I_B$q)),
                dime)
        })})})
    dimnames(C)<-c(dimnames(A)[c(I_A$c,I_A$n)],dimnames(B)[I_B$q])
  C
    }




"%.old%" <-
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