#include <RcppArmadillo.h> 
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma; 
using namespace Rcpp;


List VecSort(NumericVector IDdata,IntegerVector Auxdata) {
   NumericVector sortX = clone(IDdata);
   IntegerVector sortY = clone(Auxdata);
   
   int n = sortX.size();
   
   for(int i=0; i<(n-1); i++){
      for(int j=(i+1); j<n; j++){
          if(sortX(j) < sortX(i)){
            double temp = sortX(i);
            sortX(i) = sortX(j);
            sortX(j)=temp;
            int tempA = sortY(i);
            sortY(i) = sortY(j);
            sortY(j) = tempA;
          }
      }
   }
   
   return List::create(_["sortID"]=sortX,_["sortAux"]=sortY);
}

// [[Rcpp::export]]
NumericMatrix NormalizeD(NumericMatrix rawdata){

   int n=rawdata.nrow(),p=rawdata.ncol();
   
   NumericMatrix normdata(n,p);
   NumericVector vmean(p),vs(p);
   
   for(int k=0;k<p;k++){
      for(int i=0; i<n; i++){
        vmean(k) += rawdata(i,k);
        vs(k) += pow(rawdata(i,k),2);
     }
     vmean(k) /=n;
     vs(k) = pow((vs(k) - n*pow(vmean(k),2))/(n-1),0.5);
   }
   for(int k=0;k<p;k++){
      for(int i=0; i<n; i++){
        normdata(i,k) = (rawdata(i,k)-vmean(k))/vs(k);
     }
   }
   return normdata;
 
}


NumericMatrix NormalizeProj(NumericMatrix proj){

   int p=proj.nrow(),q=proj.ncol();
   
   NumericMatrix normproj(p,q);
   
   for(int j=0;j<q;j++){
      double ss=0;
      for(int i=0; i<p; i++){
         ss += proj(i,j)*proj(i,j);
      }
      for(int i=0; i<p; i++){
         normproj(i,j) = proj(i,j)/sqrt(ss);
      }   
   }
   return normproj;
 
}
// [[Rcpp::export]]
NumericVector NormalizeProjV(NumericVector proj){

   int p=proj.size();
   
   NumericVector normproj(p);
   
      double ss=0;
      for(int i=0; i<p; i++){
         ss += proj(i)*proj(i);
      }
      for(int i=0; i<p; i++){
         normproj(i) = proj(i)/sqrt(ss);
      }   
   return normproj;
 
}
// [[Rcpp::export]]
int FindMaxID(NumericVector X){

   int p=X.size();
   
      double maxtemp = X(0); int maxID=0;
      for(int i=1; i<p; i++){
         if(maxtemp < X(i)){
            maxtemp = X(i);
            maxID = i;
           }
      }
  
   return maxID;
 
}


// [[Rcpp::export]]

double LDAindex1(IntegerVector projclass, NumericMatrix projdata){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();  
   NumericMatrix W(p,p),WB(p,p),gsum(p,g);
   NumericVector allmean(p);
   
/*   projdata = NormalizeD(projdata);*/
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) += projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }
   
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (projdata(i,j1)-gsum(j1,l)/gn(l))*(projdata(i,j2)-gsum(j2,l)/gn(l));
                W(j2,j1) = W(j1,j2);
                double temp = (projdata(i,j1)-gsum(j1,l)/gn(l))*(projdata(i,j2)-gsum(j2,l)/gn(l))+
                               (gsum(j1,l)/gn(l)-allmean(j1))*(gsum(j2,l)/gn(l)-allmean(j2));
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
   return index;
}


// [[Rcpp::export]]
double LDAindex2(IntegerVector projclass, NumericMatrix projdata){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gnt=table(projclass);
   int g=gnt.size();   
   NumericMatrix W(p,p),WB(p,p),gsum(p,g);
   NumericVector allmean(p);
   
   projdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) += projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }
   double gn = n/g;
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (projdata(i,j1)-gsum(j1,l)/gn)*(projdata(i,j2)-gsum(j2,l)/gn);
                W(j2,j1) = W(j1,j2);
                double temp = (projdata(i,j1)-gsum(j1,l)/gn)*(projdata(i,j2)-gsum(j2,l)/gn)+
                               (gsum(j1,l)/gn-allmean(j1))*(gsum(j2,l)/gn-allmean(j2));
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
    
   return index;
}


// [[Rcpp::export]]
double Lpindex(IntegerVector projclass, NumericMatrix projdata,int r){
   printf("In Lpindex");
   double index,B=0,W=0;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix gsum(p,g);
   NumericVector allmean(p);
   
   projdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        allmean(k) += projdata(i,k)/n;
        gsum(k,(projclass(i)-1)) += projdata(i,k);
     }
   }

   for(int i=0; i<n; i++){
         int l = projclass(i)-1;
         for(int j=0;j<p;j++){
                W += pow(abs(projdata(i,j)-gsum(j,l)/gn(l)),r);
                B += pow(abs(gsum(j,l)/gn(l)-allmean(j)),r);
        }
   }
   

   index = pow(B/W,1/r);
      printf("Out Lpindex");
   return index;
}

// [[Rcpp::export]]
double PDAindex1(IntegerVector projclass, NumericMatrix projdata,double lambda){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size();   
   NumericMatrix W(p,p),WB(p,p),gsum(p,g),normdata(n,p);
   
   normdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        gsum(k,(projclass(i)-1)) += normdata(i,k);
     }
   }
   
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*(normdata(i,j2)-gsum(j2,l)/gn(l));
                
                W(j2,j1) = W(j1,j2);
                double temp = (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn(l))*(normdata(i,j2)-gsum(j2,l)/gn(l))+
                               (1-lambda)*(gsum(j1,l)/gn(l))*(gsum(j2,l)/gn(l));
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   for(int j=0; j<p; j++){
      W(j,j) = W(j,j)+ n*lambda;
      WB(j,j) = WB(j,j)+ n*lambda;  
   }
   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
    
   return index;
}


// [[Rcpp::export]]
double PDAindex2(IntegerVector projclass, NumericMatrix projdata,double lambda){

   double index;
   int n=projdata.nrow(),p=projdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gnt=table(projclass);
   int g=gnt.size();   
   NumericMatrix W(p,p),WB(p,p),gsum(p,g),normdata(n,p);

   
   normdata = NormalizeD(projdata);
   for(int i=0; i<n; i++){
      for(int k=0;k<p;k++){
        gsum(k,(projclass(i)-1)) += normdata(i,k);
     }
   }
      double gn = n/g;
   for(int i=0; i<n; i++){
         int l = projclass[i]-1;
         for(int j1=0;j1<p;j1++){
            for(int j2=0; j2<=j1; j2++) {
                W(j1,j2) += (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn)*(normdata(i,j2)-gsum(j2,l)/gn);
                W(j2,j1) = W(j1,j2);
                double temp = (1-lambda)*(normdata(i,j1)-gsum(j1,l)/gn)*(normdata(i,j2)-gsum(j2,l)/gn)+
                               (1-lambda)*(gsum(j1,l)/gn)*(gsum(j2,l)/gn);
                WB(j1,j2) += temp;
                WB(j2,j1) = WB(j1,j2);
            }
         }
   }
   for(int j=0; j<p; j++){
      W(j,j) = W(j,j)+ n*lambda;
      WB(j,j) = WB(j,j)+ n*lambda;  
   }

   Function det = base["det"];
   index =1.0-as<double>(det(wrap(W)))/as<double>(det(wrap(WB)));
    
   return index;
}



// [[Rcpp::export]]
double GINIindex1D(IntegerVector projclass, NumericVector projdata){
 
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size(), n = projclass.size();
   double n1,n2;
   double index=0,tempindex;
   List VecSortdata =  VecSort(projdata,projclass);
   NumericVector sortdata = as<NumericVector>(VecSortdata["sortID"]);
   IntegerVector sortclass = as<IntegerVector>(VecSortdata["sortAux"]);
   IntegerVector part1,part2,temptable1,temptable2;

   for(int j=0; j<g; j++){
      index += (gn(j)/n)*(1-gn(j)/n);
   }  
   
   for(int i=1; i<(n-1); i++){  
       part1 = sortclass[sortdata<=sortdata[i]];
       part2 = sortclass[sortdata>sortdata[i]];
       n1 = part1.size();
       n2 = part2.size();
       temptable1 = table(part1); int g1 = temptable1.size();  tempindex=0;
       temptable2 = table(part2); int g2 = temptable2.size(); 
       for(int j=0; j<g1; j++){
         tempindex += ((n1)/n)*(temptable1(j)/n1)*(1-temptable1(j)/n1);
       }
       for(int j=0; j<g2; j++){         
         tempindex += ((n2)/n)*(temptable2(j)/n2)*(1-temptable2(j)/n2);        
       } 
      if (tempindex < index) index = tempindex; 
            
   }   
   return 1-index;
}


// [[Rcpp::export]]

double ENTROPYindex1D(IntegerVector projclass, NumericVector projdata){
 
   Environment base("package:base");
   Function table = base["table"];
   NumericVector gn=table(projclass);
   int g=gn.size(), n = projclass.size();
   double n1,n2;
   double index=0,tempindex;
   List VecSortdata =  VecSort(projdata,projclass);
   NumericVector sortdata = as<NumericVector>(VecSortdata["sortID"]);
   IntegerVector sortclass = as<IntegerVector>(VecSortdata["sortAux"]);
   IntegerVector part1,part2,temptable1,temptable2;

   for(int j=0; j<g; j++){
      index -= (gn(j)/n)*log(gn(j)/n);
   }  
   
   for(int i=1; i<(n-1); i++){  
       part1 = sortclass[sortdata<=sortdata[i]];
       part2 = sortclass[sortdata>sortdata[i]];
       n1 = part1.size();
       n2 = part2.size();
       temptable1 = table(part1); int g1 = temptable1.size();  tempindex=0;
       temptable2 = table(part2); int g2 = temptable2.size(); 
       for(int j=0; j<g1; j++){
         tempindex -= ((n1)/n)*(temptable1(j)/n1)*log(temptable1(j)/n1);
       }
       for(int j=0; j<g2; j++){         
         tempindex -= ((n2)/n)*(temptable2(j)/n2)*log(temptable2(j)/n2);        
       } 
      if (tempindex < index) index = tempindex; 
            
   }   
   return 1-index;
}

// [[Rcpp::export]]
List PPoptimize1D(IntegerVector origclass, NumericMatrix origdata,std::string method,
             int r=1,double lambda=0,double TOL=0.00001, int maxiter=5000){

   int n=origdata.nrow(),p=origdata.ncol(),q=1;
   Environment base("package:base");
   Function table = base["table"];
   
   NumericMatrix projbest(p,q);
   for(int k=0; k<q; k++){
         projbest(_,k)= rnorm(p);
      }
   NumericMatrix projdata(n,q);
   projbest = NormalizeProj(projbest);  
   for(int i=0; i<n; i++){
      for(int k=0; k<q; k++){
         projdata(i,k)=0;
         for(int j=0; j<p; j++){
            projdata(i,k) += origdata(i,j)*projbest(j,k);
         }
      }
   }
  
  double indexbest,newindex;
  if(method=="LDA1"){
   indexbest = LDAindex1(origclass,projdata);
  } else if(method=="LDA2"){
   indexbest = LDAindex2(origclass,projdata);
  } else if(method=="Lp"){
   indexbest = Lpindex(origclass,projdata,r);
  } else if(method=="PDA1"){
   indexbest = PDAindex1(origclass,projdata,lambda);
  }  else if(method=="PDA2"){
   indexbest = PDAindex2(origclass,projdata,lambda);
  }  else if(method=="GINI"){
   indexbest = GINIindex1D(origclass,projdata);
  }  else if(method=="ENTROPY"){
   indexbest = ENTROPYindex1D(origclass,projdata);
  } 
   double temp=1;
   int kk=0;
   
   NumericVector indexkeep(maxiter);
   while(temp > TOL && kk < maxiter){
      NumericMatrix projnew(p,q);
      for(int k=0; k<q; k++){
         projnew(_,k)= rnorm(p);
      }
      projnew = NormalizeProj(projnew);  
        
      for(int i=0; i<n; i++){
         for(int k=0; k<q; k++){
            projdata(i,k)=0;
            for(int j=0; j<p; j++){
               projdata(i,k) += origdata(i,j)*projnew(j,k);
            }
         }
      }

      if(method=="LDA1"){
         newindex = LDAindex1(origclass,projdata);
      } else if(method=="LDA2"){
         newindex = LDAindex2(origclass,projdata);
      } else if(method=="Lp"){
         newindex = Lpindex(origclass,projdata,r);
      } else if(method=="PDA1"){
         newindex = PDAindex1(origclass,projdata,lambda);
      } else if(method=="PDA2"){
         newindex = PDAindex2(origclass,projdata,lambda);
      } else if(method=="GINI"){
         newindex = GINIindex1D(origclass,projdata);
      } else if(method=="ENTROPY"){
         newindex = ENTROPYindex1D(origclass,projdata);
      } 
 
      if(newindex > indexbest){
          temp = newindex - indexbest;
          indexbest = newindex; 
          projbest = projnew;     

      }
          indexkeep[kk]=newindex;      
      kk++;    
   }
         indexkeep[kk-1]=100;
   return Rcpp::List::create(Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj") = projbest,
                          Rcpp::Named("optindex") = indexbest);
}


// [[Rcpp::export]]
List PPoptimizeqD(IntegerVector origclass, NumericMatrix origdata,std::string method,int q,
             int r=1,double lambda=0,double TOL=0.00001, int maxiter=5000){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   
   NumericMatrix projbest(p,q);
   for(int k=0; k<q; k++){
         projbest(_,k)= rnorm(p);
      }
   NumericMatrix projdata(n,q);
   projbest = NormalizeProj(projbest);  
   for(int i=0; i<n; i++){
      for(int k=0; k<q; k++){
         projdata(i,k)=0;
         for(int j=0; j<p; j++){
            projdata(i,k) += origdata(i,j)*projbest(j,k);
         }
      }
   }
  
  double indexbest,newindex;
  if(method=="LDA1"){
   indexbest = LDAindex1(origclass,projdata);
  } else if(method=="LDA2"){
   indexbest = LDAindex2(origclass,projdata);
  } else if(method=="Lp"){
   indexbest = Lpindex(origclass,projdata,r);
  } else if(method=="PDA1"){
   indexbest = PDAindex1(origclass,projdata,lambda);
  }  else if(method=="PDA2"){
   indexbest = PDAindex2(origclass,projdata,lambda);
  }  else if(method=="GINI"){
   indexbest = GINIindex1D(origclass,projdata);
  }  else if(method=="ENTROPY"){
   indexbest = ENTROPYindex1D(origclass,projdata);
  } 
   double temp=1;
   int kk=0;
   
   NumericVector indexkeep(maxiter);
   while(temp > TOL && kk < maxiter){
      NumericMatrix projnew(p,q);
      for(int k=0; k<q; k++){
         projnew(_,k)= rnorm(p);
      }
      projnew = NormalizeProj(projnew);  
        
      for(int i=0; i<n; i++){
         for(int k=0; k<q; k++){
            projdata(i,k)=0;
            for(int j=0; j<p; j++){
               projdata(i,k) += origdata(i,j)*projnew(j,k);
            }
         }
      }

      if(method=="LDA1"){
         newindex = LDAindex1(origclass,projdata);
      } else if(method=="LDA2"){
         newindex = LDAindex2(origclass,projdata);
      } else if(method=="Lp"){
         newindex = Lpindex(origclass,projdata,r);
      } else if(method=="PDA1"){
         newindex = PDAindex1(origclass,projdata,lambda);
      } else if(method=="PDA2"){
         newindex = PDAindex2(origclass,projdata,lambda);
      } else if(method=="GINI"){
         newindex = GINIindex1D(origclass,projdata);
      } else if(method=="ENTROPY"){
         newindex = ENTROPYindex1D(origclass,projdata);
      } 
 
      if(newindex > indexbest){
          temp = newindex - indexbest;
          indexbest = newindex; 
          projbest = projnew;     

      }
          indexkeep[kk]=newindex;      
      kk++;    
   }
         indexkeep[kk-1]=100;
   return Rcpp::List::create(Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj") = projbest,
                          Rcpp::Named("optindex") = indexbest);
}


// [[Rcpp::export]]
List PPoptimizeAnnealqD(IntegerVector origclass, NumericMatrix origdata,std::string method,int q,
                        int r=1,double lambda=0, double TOL=0.001, int maxiter=500000,double energy=0.01, double cooling=0.999){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
   Function table = base["table"];
   double cool=cooling;
   
   NumericMatrix projbest(p,q);
   for(int k=0; k<q; k++){
         projbest(_,k)= rnorm(p);
      }
   projbest = NormalizeProj(projbest);  
   NumericMatrix projdata(n,q);  
  
   for(int i=0; i<n; i++){
      for(int k=0; k<q; k++){
         projdata(i,k)=0;
         for(int j=0; j<p; j++){
            projdata(i,k) += origdata(i,j)*projbest(j,k);
         }
      }
   }
  
   double indexbest,newindex;
   if(method=="LDA1"){
      indexbest = LDAindex1(origclass,projdata);
   } else if(method=="LDA2"){
      indexbest = LDAindex2(origclass,projdata);
   } else if(method=="Lp"){
      indexbest = Lpindex(origclass,projdata,r);
   } else if(method=="PDA1"){
      indexbest = PDAindex1(origclass,projdata,lambda);
   } else if(method=="PDA2"){
      indexbest = PDAindex2(origclass,projdata,lambda);
   } 
   
   double temp=1;
   int kk=0;
   
   NumericVector indexkeep(maxiter);
   double diff = 100;
   while((temp > TOL || abs(diff)>(energy/1000000)) && kk < maxiter){
      double tempp = energy/log(kk+1)/10000;
      temp = temp*cool;
      NumericMatrix projnew(p,q);
      for(int k=0; k<q; k++){
         projnew(_,k)= rnorm(p);
      }
      projnew = NormalizeProj(projnew);  
        
      for(int i=0; i<n; i++){
         for(int k=0; k<q; k++){
            projdata(i,k)=0;
            for(int j=0; j<p; j++){
               projdata(i,k) += origdata(i,j)*projnew(j,k);
            }
         }
      }
      if(method=="LDA1"){
         newindex = LDAindex1(origclass,projdata);
      } else if(method=="LDA2"){
         newindex = LDAindex2(origclass,projdata);
      } else if(method=="Lp"){
         newindex = Lpindex(origclass,projdata,r);
      } else if(method=="PDA1"){
         newindex = PDAindex1(origclass,projdata,lambda);
      } else if(method=="PDA2"){
         newindex = PDAindex2(origclass,projdata,lambda);
      } 
      
      NumericVector prob = runif(1);
      diff = newindex - indexbest;
      double e = exp(diff/tempp);
      if(prob[1] < e){
          for(int i=0; i<p; i++){
             for(int j=0; j<q; j++){
                projbest(i,j) = projnew(i,j);
             }
          }
          indexbest = newindex;    

      }
           indexkeep[kk]=newindex;     
      kk++;    
   }
         indexkeep[kk-1]=100;
   return Rcpp::List::create(Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj") = projbest,
                          Rcpp::Named("optindex") = indexbest);
}


// [[Rcpp::export]]
List SPSoptimize1D(IntegerVector origclass, NumericMatrix origdata,std::string method,
             int r=1, double lambda=0.7,int maxiter=100,int maxN=50){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
/*   int S=10+round(2*sqrt(p)); */
   int S = maxN; 
/*   origdata = NormalizeD(origdata); */

  NumericMatrix X(p,S),V(p,S),P(p,S),L(p,S);
  NumericVector Xindex(S),NXindex(S+2),seqID(S+2),Xindexnew(S);
  seqID(0) = S-1; seqID(S+1) = 0;
  double optindex=-1.0;
  NumericVector optproj(p);
  
   for(int k=0; k<S; k++){
         X(_,k)= NormalizeProjV(rnorm(p));
         V(_,k) = NormalizeProjV((NormalizeProjV(rnorm(p))-X(_,k))/2.0);
         P(_,k) = X(_,k);
         NumericMatrix projdata(n,1);
         for(int i=0; i<n; i++){
            projdata(i,0)=0.0;
            for(int j=0; j<p; j++){
               projdata(i,0) = projdata(i,0)+origdata(i,j)*X(j,k);
            }
         }      
         if(method=="LDA1"){
            Xindex(k) = LDAindex1(origclass,projdata);
         } else if(method=="LDA2"){
            Xindex(k) = LDAindex2(origclass,projdata);
         } else if(method=="Lp"){
            Xindex(k) = Lpindex(origclass,projdata,r);
         } else if(method=="PDA1"){
            Xindex(k) = PDAindex1(origclass,projdata,lambda);
         } else if(method=="PDA2"){
            Xindex(k) = PDAindex2(origclass,projdata,lambda);
         } else if(method=="GINI"){
            Xindex(k) = GINIindex1D(origclass,projdata);
         } else if(method=="ENTROPY"){
            Xindex(k) = ENTROPYindex1D(origclass,projdata);
         } 
     
         NXindex(k+1) = Xindex(k);
         seqID(k+1) = k; 
      }
 
      NXindex(0) = NXindex(S);NXindex(S+1) = NXindex(1); 
      
   int maxid;
   for(int i=0; i<S; i++){
      NumericVector tempindex(3); 
      for(int j=0; j<3; j++)
          tempindex(j)=NXindex(i+j);
      maxid = FindMaxID(tempindex);
      L(_,i) = P(_,seqID(i+maxid));  
      if(optindex<Xindex(seqID(i+maxid))){
         optindex = Xindex(seqID(i+maxid));
         optproj = L(_,i);
      }
   }   
   double c=(0.5+log(2));
   double w=1.0/(2*log(2));
   NumericMatrix projdata(n,1); 
   NumericVector indexkeep(maxiter);

   for(int t=0; t<maxiter;t++){
      for(int k=0; k<S; k++){  
       for(int d=0; d<p; d++){
            V(d,k)= w*V(d,k)+rnorm(1)[0]*c*(P(d,k)-X(d,k))+rnorm(1)[0]*c*(L(d,k)-X(d,k));
        } 

        X(_,k) = X(_,k)+V(_,k);
        X(_,k) = NormalizeProjV(X(_,k));
        for(int i=0; i<n; i++){
            projdata(i,0)=0;
            for(int j=0; j<p; j++){
               projdata(i,0) = projdata(i,0)+origdata(i,j)*X(j,k);
            }
        } 
        if(method=="LDA1"){
           Xindexnew(k) = LDAindex1(origclass,projdata);
        } else if(method=="LDA2"){
           Xindexnew(k) = LDAindex2(origclass,projdata);
        } else if(method=="Lp"){
           Xindexnew(k) = Lpindex(origclass,projdata,r);
        } else if(method=="PDA1"){
           Xindexnew(k) = PDAindex1(origclass,projdata,lambda);
        } else if(method=="PDA2"){
           Xindexnew(k) = PDAindex2(origclass,projdata,lambda);
        } else if(method=="GINI"){
           Xindexnew(k) = GINIindex1D(origclass,projdata);
        } else if(method=="ENTROPY"){
           Xindexnew(k) = ENTROPYindex1D(origclass,projdata);
        } 
        if(optindex<Xindexnew(k)){
           optindex = Xindexnew(k);
           optproj = X(_,k);
        }
      }
      int maxid;
      for(int i=0; i<S; i++){
         NumericVector tempindex(3);
         for(int j=0; j<3; j++)
             tempindex(j)=NXindex(i+j);
         maxid = FindMaxID(tempindex);
         L(_,i) = P(_,seqID(i+maxid));  
         if(Xindex(i)<Xindexnew(i)){
           P(_,i) = X(_,i);           
         }

      }
      Xindex = Xindexnew;
      indexkeep[t] = optindex;
   }   
 
   return Rcpp::List::create(Rcpp::Named("optindex") = optindex,
                          Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj")=optproj);
} 


// [[Rcpp::export]]
List SPSoptimizeqD(IntegerVector origclass, NumericMatrix origdata,std::string method,int q,
             int r=1, double lambda=0.7,int maxiter=100,int maxN=50){

   int n=origdata.nrow(),p=origdata.ncol();
   Environment base("package:base");
/*   int S=10+round(2*sqrt(p)); */
   int S = maxN; 
   origdata = NormalizeD(origdata);

  NumericMatrix X(p*q,S),V(p*q,S),P(p*q,S),L(p*q,S);
  NumericVector Xindex(S),NXindex(S+2),seqID(S+2),Xindexnew(S);
  seqID(0) = S-1; seqID(S+1) = 0;
  double optindex=-1.0;
  NumericMatrix optproj(p,q);
  
   for(int k=0; k<S; k++){
     NumericMatrix Xtemp(p,q),Vtemp(p,q),Ptemp(p,q);
     
     for(int l=0; l<q; l++){
      Xtemp(_,l) = NormalizeProjV(rnorm(p));
      Vtemp(_,l) = NormalizeProjV((NormalizeProjV(rnorm(p))-Xtemp(_,l))/2.0);
     }
         Ptemp = Xtemp;
         NumericMatrix projdata(n,q);
         for(int i=0; i<n; i++){
            for(int j=0; j<q; j++){
              projdata(i,j)=0.0;
              for(int jj=0; jj<p;jj++){
               projdata(i,j) = projdata(i,j)+origdata(i,jj)*Xtemp(jj,j);
              }
            }
         }  

         if(method=="LDA1"){
            Xindex(k) = LDAindex1(origclass,projdata);
         } else if(method=="LDA2"){
            Xindex(k) = LDAindex2(origclass,projdata);
         } else if(method=="Lp"){
            Xindex(k) = Lpindex(origclass,projdata,r);
         } else if(method=="PDA1"){
            Xindex(k) = PDAindex1(origclass,projdata,lambda);
         } else if(method=="PDA2"){
            Xindex(k) = PDAindex2(origclass,projdata,lambda);
         } 
         NXindex(k+1) = Xindex(k);
         seqID(k+1) = k; 
         for(int jj=0; jj<q;jj++){
           for(int j=0; j<p; j++){
              X(p*jj+j,k)=Xtemp(j,jj); 
              V(p*jj+j,k)=Vtemp(j,jj); 
              P(p*jj+j,k)=Ptemp(j,jj); 
           }
         }
      }
 
      NXindex(0) = NXindex(S);NXindex(S+1) = NXindex(1); 
    
   int maxid;
   for(int i=0; i<S; i++){
      NumericVector tempindex(3); 
      for(int j=0; j<3; j++)
          tempindex(j)=NXindex(i+j);
      maxid = FindMaxID(tempindex);
      L(_,i) = P(_,seqID(i+maxid));  
      if(optindex<Xindex(seqID(i+maxid))){
         optindex = Xindex(seqID(i+maxid));
         for(int d=0;d<p;d++)
           for(int dd=0;dd<q;dd++)
               optproj(d,dd) = L(p*dd+d,i);
      }
   }   
   double c=(0.5+log(2));
   double w=1.0/(2*log(2));
   NumericMatrix projdata(n,q); 
   NumericVector indexkeep(maxiter);

   for(int t=0; t<maxiter;t++){
       NumericMatrix Xtemp(p,q),Vtemp(p,q),Ptemp(p,q),Ltemp(p,q);
      for(int k=0; k<S; k++){  
 
         for(int jj=0; jj<q;jj++){
           for(int j=0; j<p; j++){
              Xtemp(j,jj)=X(p*jj+j,k); 
              Vtemp(j,jj)=V(p*jj+j,k); 
              Ptemp(j,jj)=P(p*jj+j,k); 
              Ltemp(j,jj)=L(p*jj+j,k); 
           }
         }
        for(int d=0; d<p; d++){
          for(int dd=0; dd<q; dd++){
            Vtemp(d,dd)= w*Vtemp(d,dd)+rnorm(1)[0]*c*(P(d,dd)-X(d,dd))+rnorm(1)[0]*c*(L(d,dd)-X(d,dd));
          }
        } 
        Vtemp = NormalizeProj(Vtemp);
        for(int dd=0; dd<q; dd++)
          Xtemp(_,dd) = Xtemp(_,dd)+Vtemp(_,dd);
        Xtemp = NormalizeProj(Xtemp);
         for(int i=0; i<n; i++){
            for(int j=0; j<q; j++){
              projdata(i,j)=0.0;
              for(int jj=0; jj<p;jj++){
               projdata(i,j) = projdata(i,j)+origdata(i,jj)*Xtemp(jj,j);
              }
            }
         }         

        if(method=="LDA1"){
           Xindexnew(k) = LDAindex1(origclass,projdata);
        } else if(method=="LDA2"){
           Xindexnew(k) = LDAindex2(origclass,projdata);
        } else if(method=="Lp"){
           Xindexnew(k) = Lpindex(origclass,projdata,r);
        } else if(method=="PDA1"){
           Xindexnew(k) = PDAindex1(origclass,projdata,lambda);
        } else if(method=="PDA2"){
           Xindexnew(k) = PDAindex2(origclass,projdata,lambda);
        } 
        if(optindex<Xindexnew(k)){
           optindex = Xindexnew(k);
           for(int dd=0; dd<q; dd++)
              optproj(_,dd) = Xtemp(_,dd);
        }
        for(int jj=0; jj<q;jj++){
            for(int j=0; j<p; j++){
              X(p*jj+j,k)=Xtemp(j,jj); 
              V(p*jj+j,k)=Vtemp(j,jj); 
           }
        }        
      }
      int maxid;
      for(int i=0; i<S; i++){
         if(Xindex(i)<Xindexnew(i)){ 
            for(int dd=0; dd<q; dd++)
              Ptemp(_,dd) = Xtemp(_,dd);
         }
         NumericVector tempindex(3);
         for(int j=0; j<3; j++)
             tempindex(j)=NXindex(i+j);
         maxid = FindMaxID(tempindex);
         L(_,i) = X(_,seqID(i+maxid));  
         for(int jj=0; jj<q;jj++){
            for(int j=0; j<p; j++){
              P(p*jj+j,i)=Ptemp(j,jj); 
           }
         } 
      }
      Xindex = Xindexnew;
      indexkeep[t] = optindex;
    }   
   return Rcpp::List::create(Rcpp::Named("optindex") = optindex,
                          Rcpp::Named("indexkeep") = indexkeep,
                          Rcpp::Named("optproj")=optproj);
}  


// [[Rcpp::export]]
List TRIBESoptimize1D(IntegerVector origclass, NumericMatrix origdata,std::string method,
             int r, double lambda,int maxiter,double TOL){

   int n=origdata.nrow(),p=origdata.ncol(); 
   Environment base("package:base");
   Function Concat = base["c"];

   List tribeALL(maxiter);
   List A;
   int L=0;
   int  nTRIBE=0; 
   int NB2 = 0, NL=0;
   NumericVector G(p),tempRN(p);
   double Gindex;    
   NumericVector Gindexkeep(maxiter);
   double diff;

   for(int ITER=0; ITER<maxiter;ITER++){       
  /*   printf("ITER=%d\n",ITER);*/
      if(ITER==0){
   /* first TRIBE */ 
         List tribeTEMP;
         NumericMatrix temp(p,1);   
         tempRN = rnorm(p);  
         tempRN = NormalizeProjV(tempRN); 
         temp(_,0) =  tempRN;

         tribeTEMP["particle"] = temp;
         tribeTEMP["P"] = temp;
   
         NumericMatrix projdata(n,1);
         for(int i=0; i<n; i++){
            projdata(i,0)=0.0;
            for(int j=0; j<p; j++){
               projdata(i,0) = projdata(i,0)+origdata(i,j)*temp(j,0);
            }
         }  
         NumericVector tempIndex(1);
         CharacterVector tempPstatus(1);
         tempPstatus = "-";
         tempIndex(0)=LDAindex1(origclass,projdata);
         if(method=="LDA1"){
            tempIndex(0) = LDAindex1(origclass,projdata);
         } else if(method=="LDA2"){
            tempIndex(0) = LDAindex2(origclass,projdata);
         } else if(method=="Lp"){
            tempIndex(0) = Lpindex(origclass,projdata,r);
         } else if(method=="PDA1"){
            tempIndex(0) = PDAindex1(origclass,projdata,lambda);
         } else if(method=="PDA2"){
            tempIndex(0) = PDAindex2(origclass,projdata,lambda);
         } else if(method=="GINI"){
            tempIndex(0) = GINIindex1D(origclass,projdata);
         } else if(method=="ENTROPY"){
            tempIndex(0) = ENTROPYindex1D(origclass,projdata);
         }          
         tribeTEMP["Pindex"]=tempIndex;
         tribeTEMP["Pstatus"]=tempPstatus;

         tribeTEMP["X"] = temp(_,0);tribeTEMP["Xindex"] = tempIndex;
         G = temp(_,0); Gindex = tempIndex(0);
         tribeTEMP["TribeStatus"]="Bad";
         CharacterVector dispStatus(1);
         dispStatus(0) = "pivot";
         tribeTEMP["DispStatus"] = dispStatus;
         tribeALL[nTRIBE] = tribeTEMP;
         NB2 = NB2+1;
      }   else if(ITER==1){
   /*    2nd TRIBE */
         nTRIBE++;
         List tribeTEMP;
         NumericMatrix tempM(p,2);   
         List tempL;
         tempL = tribeALL[0];
         A=tempL;
         tempM(_,0) = NormalizeProjV(runif(p)); 
         tempRN = rnorm(p); tempM(_,1) = NormalizeProjV(tempRN);

         NumericVector IndexK(2),XbestP(p);
         XbestP = tempL["X"]; 
         double r=0;
         for(int i=0; i<p; i++){
            r = r + pow(abs(XbestP(i)-G(i)),2);
         }
         r = pow(r,0.5);
         for(int i=0; i<p; i++){
            tempM(i,1) = tempM(i,1)*r+G(i);         
         }

         tempM(_,0) = NormalizeProjV(tempM(_,0)); 
         tempM(_,1) = NormalizeProjV(tempM(_,1));
 
         tribeTEMP["particle"]= tempM;
         tribeTEMP["P"] = tempM;
 
         NumericMatrix projdata(n,1);

         for(int k=0; k<2; k++){
            for(int i=0; i<n; i++){
               projdata(i,0)=0.0;
               for(int j=0; j<p; j++){
                  projdata(i,0) = projdata(i,0)+origdata(i,j)*tempM(j,k);
               }
            }  

            if(method=="LDA1"){
               IndexK(k) = LDAindex1(origclass,projdata);
            } else if(method=="LDA2"){
               IndexK(k) = LDAindex2(origclass,projdata);
            } else if(method=="Lp"){
               IndexK(k) = Lpindex(origclass,projdata,r);
            } else if(method=="PDA1"){
               IndexK(k) = PDAindex1(origclass,projdata,lambda);
            } else if(method=="PDA2"){
               IndexK(k) = PDAindex2(origclass,projdata,lambda);
            } else if(method=="GINI"){
               IndexK(k) = GINIindex1D(origclass,projdata);
            } else if(method=="ENTROPY"){
               IndexK(k) = ENTROPYindex1D(origclass,projdata);
            } 
         }   
         tribeTEMP["Pindex"]=IndexK;  
         tribeTEMP["X"]=tempM(_,0); 
         NumericVector Xbest(1); Xbest(0) = IndexK(0);       
         
         CharacterVector Pstatus(2),dispStatus(2);
    
         for(int k=0; k<2; k++){
            if(IndexK(k)>Gindex){
               Gindex = IndexK(k);
               G=tempM(_,k);           
            }  
            if(IndexK(k)>Xbest(0)){
               Xbest(0) = IndexK(k);
               tribeTEMP["X"]=tempM(_,k);           
            }          
            Pstatus(k) = "-";     
            dispStatus(k) = "pivot";
         }

         NumericVector Pbest(2);
         Pbest = IndexK;
         tribeTEMP["Pstatus"]=Pstatus;
         tribeTEMP["Xindex"] = Xbest;
         tribeTEMP["TribeStatus"]="Bad"; 
         tribeTEMP["DispStatus"] = dispStatus;

         tribeALL[nTRIBE]= tribeTEMP;
         NB2 = NB2+2*2;
         NL = NB2+nTRIBE*(nTRIBE+1);
        
      }   else {
         /* displacement */

         for(int kk=0; kk<=nTRIBE; kk++) {
            List tempTRIBE;
            tempTRIBE=tribeALL[kk];

            CharacterVector oldStatus = tempTRIBE["Pstatus"];
            CharacterVector Dstatus = tempTRIBE["DispStatus"];

            NumericMatrix oldParticle = tempTRIBE["particle"],oldP = tempTRIBE["P"],newP = oldP;
            int pOP = oldParticle.ncol(), nOP = oldParticle.nrow();
            CharacterVector tempPstatus(pOP),tempDstatus(pOP),oldPstatus = tempTRIBE["Pstatus"];            
            NumericVector tempIndex(pOP),oldPindex = tempTRIBE["Pindex"];
            NumericVector  XbestIndexT = tempTRIBE["Xindex"];
            double XbestIndex = XbestIndexT(0);
            NumericMatrix newParticle(nOP,pOP);
            NumericVector Xbest(p);
            double tempXindex;
            int nGood=0;
            for(int i=0; i<pOP; i++){
               tempDstatus(i) = "pivot";
               if(Dstatus(i) == "pivot"){
                  NumericVector temp1(p),temp2(p);
                  double P1 = XbestIndex/(Gindex+XbestIndex),P2 = Gindex/(Gindex+XbestIndex);
                  double r=0;
                  for(int j=0; j<p; j++){
                     r = r+pow((oldP(j,i)-G(j)),2);
                  }
                  r = pow(r,0.5);
                 
                  temp1 = NormalizeProjV(rnorm(p)); temp2 = NormalizeProjV(rnorm(p));
                  for(int j=0;j<p; j++){
                     newParticle(j,i) = P1*(temp1(j)*r+ oldP(j,i))+P2*(temp2(j)*r+ G(j));
                  }
                  newParticle(_,i) = NormalizeProjV(newParticle(_,i));    

               } else if (Dstatus(i) == "Dpivot"){
                  NumericVector temp1(p),temp2(p);
                  double P1 = XbestIndex/(Gindex+XbestIndex),P2 = Gindex/(Gindex+XbestIndex);
                  double r=0;
                  for(int j=0; j<p; j++){
                     r = r+pow((oldP(j,i)-G(j)),2);
                  }
                  r = pow(r,0.5);
                 
                  temp1 = NormalizeProjV(rnorm(p)); temp2 = NormalizeProjV(rnorm(p));
                  NumericVector b= rnorm(1);
                  b(0) = b(0) * abs(P1-P2);
                  
                  for(int j=0;j<p; j++){
                     newParticle(j,i) = (P1*(temp1(j)*r+ oldP(j,i))+P2*(temp2(j)*r+ G(j)))*(1+b(0));
                  }                               
                  newParticle(_,i) = NormalizeProjV(newParticle(_,i));  

               } else if(Dstatus(i) == "indepG"){
                  NumericVector temp1=rnorm(p);              
                  for(int j=0 ; j<p; j++){
                     newParticle(j,i) = 2*G(j)+temp1(j)*pow(abs(oldParticle(j,i)-G(j)),0.5)+oldParticle(j,i);
                  }                     
                  newParticle(_,i) = NormalizeProjV(newParticle(_,i)); 

               }
               NumericMatrix projdata(n,1);
               for(int ii=0; ii<n; ii++){
                  projdata(ii,0)=0.0;   
                  for(int j=0; j<p; j++){
                      projdata(ii,0) = projdata(ii,0)+origdata(ii,j)*newParticle(j,i);                
                  }
               }    
 
               if(method=="LDA1"){
                  tempIndex(i) = LDAindex1(origclass,projdata);
               } else if(method=="LDA2"){
                  tempIndex(i) = LDAindex2(origclass,projdata);
               } else if(method=="Lp"){
                  tempIndex(i) = Lpindex(origclass,projdata,r);
               } else if(method=="PDA1"){
                  tempIndex(i) = PDAindex1(origclass,projdata,lambda);
               } else if(method=="PDA2"){
                  tempIndex(i) = PDAindex2(origclass,projdata,lambda);
               } else if(method=="GINI"){
                  tempIndex(i) = GINIindex1D(origclass,projdata);
               } else if(method=="ENTROPY"){
                  tempIndex(i) = ENTROPYindex1D(origclass,projdata);
               } 
               if(oldPindex(i)>tempIndex(i)){
                 tempPstatus(i) = "-"; 
                 tempDstatus(i) = "pivot";                
               } else if (oldPindex(i)==tempIndex(i)){
                 tempPstatus(i) ="=";
                 if(oldPstatus(i) =="+"){
                    tempDstatus(i) =  "Dpivot";
                 } else{
                    tempDstatus(i) =  "pivot";                   
                 }
               } else{
                  tempPstatus(i) = "+";
                  newP(_,i) = newParticle(_,i);
                  if(oldPstatus(i) =="-"){
                     tempDstatus(i) =  "Dpivot";
                  } else{
                     tempDstatus(i) =  "indepG";                   
                  }
               }
               if(i==0){
                  XbestIndex = tempIndex(i);
                  Xbest = newParticle(_,i);
               } else if(XbestIndex < tempIndex(i)){
                  XbestIndex = tempIndex(i);
                  Xbest = newParticle(_,i);                                
               } 
               if(Gindex < tempIndex(i)){
                  Gindex = tempIndex(i);
                  G = newParticle(_,i);
               }
               if(tempPstatus(i)=="+")
                  nGood++;
               
            }   
            NumericVector tempU=runif(1);
            CharacterVector TribeStatus(1);
            if(nGood ==0){
               tempTRIBE["TribeStatus"] = "Bad";
            } else if(tempU(0)<0.5){
               tempTRIBE["TribeStatus"] = "Good"; 
            } else{
               tempTRIBE["TribeStatus"] = "Bad"; 
            }
            tempTRIBE["Pindex"]=tempIndex;
            tempTRIBE["particle"] =newParticle;
            tempTRIBE["Pstatus"] = tempPstatus;
            tempTRIBE["X"]=Xbest;
            tempTRIBE["Xindex"]=XbestIndex;
            tempTRIBE["DispStatus"]=tempDstatus;
            tribeALL[kk] = tempTRIBE;
         }  
         if(ITER < (NL*0.5)) {
          /* evolution of Tribe */
            List tempP;
            NumericVector Gkeep(p),Gindexkeep(1);
            
            for(int kk=0; kk <nTRIBE; kk++){
               List tempL = tribeALL[kk];
               CharacterVector tempTstatus = tempL["TribeStatus"];
               NumericMatrix oldParticle=tempL["particle"];
               NumericMatrix oldP = tempL["P"];
               NumericVector oldPindex = tempL["Pindex"];
               CharacterVector oldPstatus = tempL["Pstatus"];
               CharacterVector oldDstatus = tempL["DispStatus"];
               NumericVector oldX = tempL["X"];
               NumericVector oldXindex = tempL["Xindex"];
               int nP = oldPindex.size();
               if(tempTstatus(0)=="Good"){          
                 if(nP >1){
                 /* Removal from good Tribe */
                  int worstID=0;
                  double tempMin=oldPindex(0);
                  for(int i=0; i<nP; i++){
                      if(tempMin >oldPindex(i)){
                        worstID = i; 
                      }
                  }

                  int newnp = nP-1;
                  NumericMatrix tempParticle(p,newnp);
                  NumericMatrix tempP(p,newnp);
                  NumericVector tempPindex(newnp);
                  CharacterVector tempPstatus(newnp);
                  CharacterVector tempDstatus(newnp);
                  for(int i=0;i<nP; i++){
                     if(i < worstID){
                        tempParticle(_,i) = oldParticle(_,i); 
                        tempP(_,i) = oldP(_,i); 
                        tempPindex(i) = oldPindex(i); 
                        tempPstatus(i) = oldPstatus(i); 
                        tempDstatus(i) = oldDstatus(i);                         
                     } else if(i > worstID) {
                        tempParticle(_,(i-1)) = oldParticle(_,i); 
                        tempP(_,(i-1)) = oldP(_,i); 
                        tempPindex(i-1) = oldPindex(i); 
                        tempPstatus(i-1) = oldPstatus(i); 
                        tempDstatus(i-1) = oldDstatus(i);                                          
                     }
                  }   
                  tempL["particle"] = tempParticle;
                  tempL["P"] = tempP;
                  tempL["Pindex"] = tempPindex;
                  tempL["Pstatus"] = tempPstatus;
                  tempL["DispStatus"] = tempDstatus;
                  tribeALL[kk] = tempL;   
                  NB2 = NB2 - nP*nP + newnp*newnp;
                 }
               } else if(tempTstatus(0)=="Bad"){
                   /* Addition for bad tribes */
                  NumericVector temp1(p),temp2(p); 
                  temp1 = NormalizeProjV(runif(p)); 
                  temp2 = NormalizeProjV(rnorm(p));
                  NumericVector XbestP(p);
                  XbestP = tempL["X"];                
                  double r=0;
                  for(int i=0; i<p; i++){
                     r = r + pow((XbestP(i)-G(i)),2);
                  }
                  r = pow(r,0.5);
                  for(int i=0; i<p; i++){
                     temp2(i) = temp2(i)*r+G(i);         
                  }                  
                  temp1 = NormalizeProjV(temp1); 
                  temp2 = NormalizeProjV(temp2);                 
                  List T;
                  T["1"] = temp1;T["2"]=temp2;                  
                  tempP = Concat(tempP,T);
               } 
            }
            int newPn = tempP.size();
            if(newPn!=0){
               List tempLL;
               NumericVector tempPart(p);
               NumericMatrix tParticle(p,newPn),tP(p,newPn);
               NumericVector tPindex(newPn),tXindex(1),tX(p),tGindex(1);
               CharacterVector tPstatus(newPn),tDstatus(newPn),tTstatus(1);
                    
               for(int i=0; i<newPn; i++){
                  tempPart = tempP[i];
                  tParticle(_,i) = tempPart;
                  tDstatus(i) = "pivot";
                  tP(_,i) = tempPart;
                  NumericMatrix projdata(n,1);
                  for(int i=0; i<n; i++){
                     projdata(i,0)=0.0;
                     for(int j=0; j<p; j++){
                        projdata(i,0) = projdata(i,0)+origdata(i,j)*tempPart(j);
                     }
                  }  
                  double tempI;                       
                  if(method=="LDA1"){
                     tempI = LDAindex1(origclass,projdata);
                  } else if(method=="LDA2"){
                     tempI = LDAindex2(origclass,projdata);
                  } else if(method=="Lp"){
                     tempI = Lpindex(origclass,projdata,r);
                  } else if(method=="PDA1"){
                     tempI = PDAindex1(origclass,projdata,lambda);
                  } else if(method=="PDA2"){
                     tempI = PDAindex2(origclass,projdata,lambda);
                  } else if(method=="GINI"){
                     tempI = GINIindex1D(origclass,projdata);
                  } else if(method=="ENTROPY"){
                     tempI = ENTROPYindex1D(origclass,projdata);
                  }  
                  tPindex(i) = tempI;
                  tPstatus(i) = "-";
                                          
                  if(i==0){
                     tXindex(0) = tempI; 
                     tX = tempPart;
                  } else if(tXindex(0) < tempI){
                     tXindex(0) = tempI;
                     tX = tempPart;
                  }
                  if(Gindex < tempI){
                     Gindex = tempI;
                     G = tempPart;
                  }                                                 
               }

               tTstatus(0) = "Bad";

               tempLL["particle"] = tParticle;
               tempLL["P"] = tP;
               tempLL["Pindex"] = tPindex;
               tempLL["Pstatus"] = tPstatus;
               tempLL["DispStatus"] = tDstatus;   
               tempLL["X"] = tX;
               tempLL["Xindex"] = tXindex;
               tempLL["TribeStatus"]=tTstatus;
               nTRIBE ++;
               tribeALL[nTRIBE] = tempLL;
               NB2 = NB2+newPn*newPn;
            }
            NL = NB2+nTRIBE*(nTRIBE+1);
         }
      }
      Gindexkeep(ITER) = Gindex;

      if(ITER>30){
        diff = Gindex-Gindexkeep(ITER-10);
      /*  printf("diff=%f\n",diff); */
        if(diff <TOL){
           return Rcpp::List::create(Rcpp::Named("optindex") = Gindex,
                          Rcpp::Named("indexkeep") = Gindexkeep,
                          Rcpp::Named("optproj")=G);
         
        }
        
      } 
   }
  
   return Rcpp::List::create(Rcpp::Named("optindex") = Gindex,
                          Rcpp::Named("indexkeep") = Gindexkeep,
                          Rcpp::Named("optproj")=G);
} 

