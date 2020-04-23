#' Converts Array to matrix
#' @param A An (eventually named) array
#' @param n a subvector of 1:length(dim(A)) or names(dimnames(A)).
#' @param p a subvector of 1:length(dim(A)) or names(dimnames(A)).
A2M<-function(A,n,p=if(is.character(n)){setdiff(names(dimnames(A),n))}else{setdiff(1:length(dim(A)),n)}){
  if(is.character(n)){n<-match(n,names(dimnames(A)))}
  if(is.character(p)){p<-match(p,names(dimnames(A)))}
  array(aperm(A,c(n,p)),c(prod(dim(A)[n]),prod(dim(A)[p])))}





#' Extracts dimension from array
#' @param A a named array
#' @param a a list of dimensions of A ( a subvector of 1:length(dim(A)) or a subvector of dimnames(A))
#' @param ... a vector the same length of a of integers. necessarily ...[i]<=dim(A)[a[i]]
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
#' @param A An (eventually named) array of dimension \eqn{\mathrm{dim}(A)=(a_i)_{i\in I_A}}
#' @param B An (eventually named) array of dimension \eqn{\mathrm{dim}(B)=(b_j)_{j\in I_B}}
#' @param I_A a named list of subvectors from names(dimnames(A)) or from 1:length(dim(A)). \eqn{I_A=(I_A^{(c)},I_A^{(n)},I_A^{(p)})}. 
#' @param I_B a named list of subvectors from names(dimnames(B)) or from 1:length(dim(B)). \eqn{I_B=(I_B^{(c)},I_B^{(p)},I_B^{(q)})}. 
#' Necessarily, \eqn{(\mathrm{dim}(A))_{I_A^{(p)}}=(\mathrm{dim}(B))_{I_B^{(p)}}} (e.g dim(A)[I_A$p]==dim(B)[I_B$p]) and \eqn{(\mathrm{dim}(A))_{i\in I_A^{(c)}}=(\mathrm{dim}(B))_{i\in I_B^{(c)}}} e.g dim(A)[I_A$c]==dim(B)[I_B$c])
#' @return C=AB the array of dimension \eqn{\left((a_\ell)_{\ell\in I_A^{(c)}},(a_i)_{i\in I_A^{(n)}}(b_j)_{j\in I_B^{(q)}}\right)} 
#'  defined by 
#'  \deqn{\forall (\ell_1,\ldots,\ell_C)\in\prod_{i\in I_A^{(c)}}\{1,\ldots, a_i\},}
#'  \deqn{\forall (i_1,\ldots,i_N)\in\prod_{i\in I_A^{(n)}}\{1,\ldots ,a_i\},} 
#'  \deqn{\forall (j_1,\ldots,j_Q)\in\prod_{j\in I_B^{(q)}}\{1,\ldots ,b_j\},} 
#' \deqn{C[\ell_1,\ldots,\ell_C,i_1,\ldots,i_N,j_1,\ldots,j_Q]=\sum_{k_1=1}^{K_1}\!\!\ldots\!\!\sum_{k_P=1}^{K_P}A^\star[\ell_1,\ldots,\ell_C,i_1,...,i_q,k_1,...,k_p]\times B^\star[\ell_1,\ldots,\ell_C,k_1,...,k_p,j_1,...,j_n])}
#' where \eqn{A^\star} and \eqn{B^\star} are multidimensional transposition of \eqn{A} and \eqn{B} and \eqn{K_1,\ldots,K_P=\mathrm{dim}(A)_{I_A^{(p)}}}.
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
           I_B=list(c=integer(0),p=integer(0),q=1:length(dim(B))),
           requiresameindices=F){
    I_A=lapply(I_A,function(x){if(is.character(x)){
      match(x,names(dimnames(A)))}else{x}})
    I_B=lapply(I_B,function(x){if(is.character(x)){
      match(x,names(dimnames(B)))}else{x}})
    dime<-c(dim(A)[I_A$n],dim(B)[I_B$q])
    if(length(I_A$c)>0){D=expand.grid(dimnames(A)[I_A$c],stringsAsFactors = FALSE)
    C<-plyr::maply(.data=D,.fun=function(...){
      a=extractA(A,I_A$c,...)
      b=extractA(B,I_B$c,...)
      if(requiresameindices){b<-extractA(b,I_B$p,dimnames(A)[I_A$p])}
      array(A2M(a,c(I_A$c,I_A$n),c(I_A$p))%*%(A2M(b,c(I_B$p),c(I_B$c,I_B$q))),
            dime)
    })
      
      }else{
        C<-array(A2M(A,c(I_A$n),c(I_A$p))%*%A2M(B,I_B$p,I_B$q),
                    dime)
            }
    dimnames(C)<-c(dimnames(A)[c(I_A$c,I_A$n)],dimnames(B)[I_B$q])
    C
  }

