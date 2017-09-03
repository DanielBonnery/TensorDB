
A2M<-function(A,n,p=if(is.character(n)){setdiff(names(dimnames(A),n))}else{setdiff(1:length(dim(A)),n)}){
  if(is.character(n)){n<-match(n,names(dimnames(A)))}
  if(is.character(p)){p<-match(p,names(dimnames(A)))}
  array(aperm(A,c(n,p)),c(prod(dim(A)[n]),prod(dim(A)[p])))}






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
    D=if(length(I_A$c)>0){expand.grid(dimnames(A)[I_A$c])
      C<-plyr::maply(.data=D,.fun=function(...){
        plyr::aaply(extractA(A,I_A$n,...),I_A$c,function(a){
          plyr::aaply(extractA(B,I_B$q,...),I_B$c,function(b){
            array(A2M(a,c(I_A$c,I_A$n),c(I_A$p))%*%A2M(b,c(I_B$p),c(I_B$c,I_B$q)),
                  dime)
          })})})
      
      }else{
        C<-array(A2M(A,c(I_A$n),c(I_A$p))%*%A2M(B,I_B$p,I_B$q),
                    dime)
            }
    dimnames(C)<-c(dimnames(A)[c(I_A$c,I_A$n)],dimnames(B)[I_B$q])
    C
  }

