---
title: ""
author: "Polymorphism calculator for molecular markers"
date: "Apr 14, 2017"
output: html_document
runtime: shiny
---
```{r echo=FALSE}
options(width = 1000)
fluidPage(br(),
  HTML('<!DOCTYPE html>
       <html>
       <head>
       <style type="text/css">
       h1 {color:#951405;font-size: 40px;}
       h4 {color:black;font-size: 19px;}
       h3 {color:#550404;}
       p {color:#683203;}
       </style>
       </head>
       <body>
       <div style="padding:20px;position:relative;border-radius:300px;border:20px solid #951405";">
       <div style="opacity:0.1;position:absolute;bottom:0px;right:0px;left:0px;height:144px;border-radius:300px;background-color:#F3EC0D"></div>
       <h1 style="font-family:Copperplate Gothic Bold;letter-spacing:8px;text-align:center;"><i>i</i>MEC</h1>
       
       <h4 style = "font-family:Charlesworth;text-align:center;color:navy;"><em> Marker Efficieny Calculator</em> </h4>        
       
       </div>
       </body>
       </html>
       
       '))
```
<style type="text/css">

h1.title {
  font-family: Copperplate Gothic Bold;
  letter-spacing:12px;
  font-size: 38px;
  color: DarkRed;
  text-align: center;
}
h4.author { /* Header 4 - and the author and data headers use this too  */
    font-size: 27px;
  font-family: "Times New Roman", Times, serif;
  color: DarkRed;
  text-align: center;
}
h4.date { /* Header 4 - and the author and data headers use this too  */
  font-size: 18px;
  font-family: "Times New Roman", Times, serif;
  color: DarkBlue;
  text-align: center;
}
</style>


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r, echo=FALSE}
library(WriteXLS)
```

```{r, include=FALSE}

#Defining the functions for reading the nexus file

read.phy<- function(file){#reading the phylip file and return a list of species name and their corresponding characters
  p<- readLines(file)
  s<-p[1]
  d<-as.numeric(unlist(strsplit(s, " ")))
  i<-1
  while(is.na(d[i])){#eliminating the empty spaces in reading the nomber of species
    i<- i+1  
  }
  no.sp<- d[i]
  while(is.na(d[i+1])){#eliminating the empty spaces in reading the number of characters
    i<-i+1
  }
  no.ch<- d[i+1]
  
  if(length(p)-1 != no.sp){
    stop("Check the Phylip file: The number of species specified in the header is not equal to the one presented in the file")
  }
  sp.name<- character(no.sp)
  data<- list()
  for (i in 2:length(p)){
    ss<- unlist(strsplit(p[i], " "))
    sp.name[(i-1)]<- ss[1]
    if(length(ss)==2){#check if the characters are binary and not seperated
      ss<- c(ss[1], unlist(strsplit(ss[2], "")))
    }
    if(length(ss)-1 != no.ch){
      stop("Check the Phylip file: The number of characters specified in the header is not equal to the one presented in the file")
    }
    data[[(i-1)]]<- as.vector(ss[2:length(ss)])
  }
  names(data)<- sp.name
  return(data)
}

is.binary<- function(vec){#for a vector checks if the elements are mostly 0 and 1 (>80%)
  l<- length(vec)
  if (is.na(vec[1])){return(FALSE)}
  else if((sum(vec=="1")+sum(vec=="0"))>l*0.8){return(TRUE)}
  else {return(FALSE)}
}

is.input<- function(line){#on the top of the previous funciton for each line of the nexus file return a true or false if that is the input file
  if(is.binary(strsplit(strsplit(line, " ")[[1]][2], "")[[1]])){return(TRUE)}
  else{return(FALSE)}
}

input.cut<- function(nexus){#give a nexus input file and cut the input section of it
  out<- list()
  for (i in 1:length(nexus)){
    if(is.input(nexus[i])){
      out[[i]]<- nexus[i]
    } 
  }
  return(out[which(out!="NULL")])
}


read.nex<- function(file){#making the phylip block area without the header
  n<- readLines(file)
  out<- input.cut(n)
  for (i in 1:length(out)){
    names(out)[i]<- strsplit(out[[i]], " ")[[1]][1]
    out[[i]]<- unlist(strsplit(strsplit(out[[i]], " ")[[1]][2], ""))
  }
  return(out)
}


ls2li<- function(list){
  out<- character(length(list)+1)
  out[1]<- paste(length(list),length(list[[1]]))#length(strsplit(strsplit(list[[1]], " ")[[1]], "")[[1]]), sep=" ")
  for(i in 1:length(list)){
    out[(i+1)]<- paste(names(list)[i], paste(list[[i]], collapse = "") , sep=" ")#strsplit(list[[i]], " ")[[1]], sep=" ")
  }
  return(out)
}


```


<br></br>
<center><p> `r icon("home", lib = "glyphicon")` **iMEC** calculates the polymorphism indices as defined in *iMEC: Online Marker Efficiency Calculator.* Applications in Plant Science (*submitted*).</p></center>


<center> <h3>**INPUT**</h3> </center>

##### To start calculating, please provide your [`Phylip`](http://evolution.genetics.washington.edu/phylip.html) or [`Nexus`](https://www.ncbi.nlm.nih.gov/pubmed/11975335) file. This should result in updated **data** similar to the template. 

`r icon("exclamation-sign", lib = "glyphicon")` If the **data** did not appear or does not look like exisiting example, it means that you file not compatible with a standard input files. Please check your file and upload again. Here are `r downloadLink('biT', 'binary')` or `r downloadLink('coT', 'codominant')` `Phylip` and `r downloadLink('neX', 'nexus')` file as guidlines and test.
```{r, echo=FALSE}
binary <- readLines("Marker_Phylip_binary.txt")
codominant <- readLines("Marker_Phylip_codominant.txt")
nexxxus<- readLines("nexus.txt")
output$coT <- downloadHandler(
    filename = "Marker_Phylip_codominant.txt",
    content = function(file) {
      writeLines(codominant, file)
      #write.csv(out, file, row.names = FALSE,fileEncoding = "macroman", sep="\t")
    }
  )
output$biT <- downloadHandler(
    filename = "Marker_Phylip_binary.txt",
    content = function(file) {
      writeLines(binary, file)
      #write.csv(out, file, row.names = FALSE,fileEncoding = "macroman", sep="\t")
    }
  )
output$neX <- downloadHandler(
    filename = "nexus.txt",
    content = function(file) {
      writeLines(nexxxus, file)
      #write.csv(out, file, row.names = FALSE,fileEncoding = "macroman", sep="\t")
    }
  )
```







```{r, echo=FALSE}
inputPanel(
    radioButtons("inputType", "Input format", c("Phylip" = "phy", "Nexus"= "nex")),
    fileInput("inpuut", label = "Upload file"))

renderTable({
  
  data<-input$inpuut
  po<-data$datapath
  if (is.null(po)) {print(c("3 10", "SpeciesA 0111011001", "SpeciesB 1010100110", "SpeciesC 0001011011"))}
  else { if (input$inputType=="phy"){
    print(head(readLines(po)))
  }
    else if (input$inputType=="nex"){
      print(head(ls2li(read.nex(po))))
    }
    else {print("no File")}
  }
})
```


##### If the heading of your **data** above seems satisfactory, please provide the number of primers followed by the number of bands for each primer as comma separated numbers and click ***e*PIC**.

`r icon("info-sign", lib = "glyphicon")` You can skip providing input for this section. This will return the *SSR* values. 





```{r, echo=FALSE}
inputPanel(
  sliderInput("n_bins", label = "Number of primers:",
              min= 1, max=40, value=1, step= 1),
  
  textInput("band_sizes", "Number of bands with each primer", placeholder = "10, 20, 15, ...")
)
actionButton("Go", "iMEC")
```

<center> <h3>**OUTPUT** </h3> </center>
##### If you have followed the upper steps succesfully, below you find the result of the calculations. 




```{r, echo=FALSE}

renderDataTable({
  if(input$Go){
    data<-input$inpuut
    po<-data$datapath
    if (is.null(po)) {print("No File input")}
    else { if(input$inputType=="phy"){
    phylip<- read.phy(po)}
      else if(input$inputType=="nex") {
        phylip<- read.nex(po)}
    if (input$band_sizes==""){vec<-1}
    else {
    vec<- as.numeric(unlist(strsplit(input$band_sizes, split = ",")))
    }
    out<<-Calculator(phylip, vec)
    print(out)
    }
  }
})


```


*If your output seems satisfactory, you may download that via tab below.* 

```{r, echo=FALSE}

inputPanel(
downloadButton("Download", "Download Table")
)

output$Download <- downloadHandler(
    filename = function() {
      paste("Table-", Sys.Date(), ".xls", sep="")
    },
    content = function(file) {
      WriteXLS(out, file)
      #write.csv(out, file, row.names = FALSE,fileEncoding = "macroman", sep="\t")
    }
  )
```



`r icon("exclamation-sign", lib = "glyphicon")` If you did not got the output or prompted with an error, it means that either your data is not compatible or you have provided a number of bands for each primer which does not perfectly partition the number of characters. Amend your inputs and rerun.

`r icon("info-sign", lib = "glyphicon")` Indices

Acronyms| Definition | Details                                   
------ | ------------------------------  |  --------------------------------------------------------------------------------------------------------
**H** | *Heterozygocity index* | Also called expected heterozygosity *HE* of diversity index *DI* in some cases. It is defined as the probability that an individual is heterozygous for the locus in the population.
**PIC** | *Polymorphic Information Content* | Probability that the marker genotype of a given offspring will allow deduction, in the absence of crossing over, of which of the two marker alleles of the affected parents it received.  
**E** | *Effective multiplex ratio*  | Also abbreviated as *EMR*, it is the product of the fraction of polymorphic loci for an individual assay. In other words the number of loci polymorphic in the germplasm set of interest analyzed per experiment. Except for *SSR*s (or other co-dominant markers) where *E* is 1 based on the assumption that each assay reveals a single locus.  
**H.av**| *Arithmatic mean of H* | Also referred to as average heterozygosity for all markers, equals the average heterozygosity for polymorphic markers multiplied by the fraction of markers that are polymorphic 
**MI**| *Marker Index* | This is defined as the product of the effective multiplex ratio *E* and the arithmetic mean heterozygosity *H.av* for the polymorphic bands in an assay.
**D**| *Discriminating power* |  Represents the probability that two randomly chosen individuals have different patterns, and thus are distinguishable from one another. 
**R**| *Resolving power* | Also abbreviated as *Rp*, is the ability of a primer or technique to distinguish between large number of genotypes. It is calculated as the sum of band informativeness *Ib* calculated for the banding pattern generated by the assay, where *Ib* measures the similarity of a band to the optimal condition, that *50%* (presence or absence) of genotypes contain the band.    


`r icon("info-sign", lib = "glyphicon")` Each of the values are defined precisely in the related paper mentioned above. You can restrict your view with searching a specific index or a primer. ’0’ refers to the *SSR* values, ’1’ to first primer, ’2’ to second primer and so forth. Note that *resolving power* is not defined for codominant markers.

`r icon("info-sign", lib = "glyphicon")` Note that the calculator reads the numbers of primers from the vector of numbers you provided for the number of band section. 

******

`r icon("book", lib = "glyphicon")` The formulae related to these calculations can be found in the original paper mentioned above. Please cite that upon using this program.

`r icon("cog", lib = "glyphicon")` The codes written in R can be downloaded `r downloadLink('Downloadcode', 'here')`.


`r icon("envelope", lib = "glyphicon")` ali.amiryousefi@helsinki.fi or peter.poczai@helsinki.fi.




```{r, echo=FALSE}
#inputPanel(
#downloadButton("Downloadcode", "Download Code")
#)
Rcodes <- readLines("PIC_v0.04.R")
output$Downloadcode <- downloadHandler(
    filename = function() {
      paste("iMEC", Sys.Date(), ".R", sep="")
    },
    content = function(cong) {
      writeLines(Rcodes, cong)
      #write.csv(out, file, row.names = FALSE,fileEncoding = "macroman", sep="\t")
    }
  )
```




```{r, echo=FALSE}
Calculator<- function(phy, Bl=c(1)){#the first output is H and the second is PIC
  
  #Defining the functions
  
  DisPow<- function(m){
    if (max(m)==1){
      n<- length(m)
      p1<- sum(m==1)/n
      D<- 1-(p1*(sum(m==1)-1)/(n-1))
    }
    else {
      coln<- dim(m)[2]
      Ds<- numeric(coln)
      for (i in 1:coln){
        col<- m[,i]
        uni<-unique(col)
        p<- numeric(length(uni))
        c<- numeric(length(uni))
        for (j in 1:length(uni)){
          p[j]<- sum(uni[j]==col)/length(col)
          c[j]<- p[j]*(sum(uni[j]==col)-1)/(length(col)-1)
        }
        C<- sum(c)
        Ds[i]<-1-C 
      }
      D<-mean(Ds)
    }
    return(D)
  }
  
  ResPow<- function(m){
    if(max(m)==1){
      coln<- dim(m)[2]
      Ib<- numeric(coln)
      for (i in 1:coln){
        col<- m[,i]
        p<- sum(col==1)/length(col)
        Ib[i]<- 1-(2*abs(0.5-p))
      }
      Rp<- sum(Ib)
      return(Rp)
    }
  }
  
  m<- matrix(0, length(phy), length(phy[[1]]))
  for (i in 1:length(phy)){
    m[i,]<- (phy[[i]])
  }
  alfbt<-unique(as.vector(m))
  n<- dim(m)[1]*dim(m)[2]
  prop<- numeric(length(alfbt))
  for (i in 1:length(prop)){
    prop[i]<- sum(alfbt[i]==m)/n
  }
  h<-1-sum(prop^2)
  p2<-(prop^2 %o% prop^2)
  ex<- 2*sum(p2[lower.tri(p2)])
  pic<- h-ex
  if(max(m)==1){
    n1<- sum(m==1)
    beta<- n1/n
    e<-dim(m)[2]*beta
    n2<- sum(m==0)#2
    h1<-1-(n1/n)^2
    h2<-1-(n2/n)^2
    h_av<- (h)/n#(h1+h2)/n
    mi<-  h1*e#dim(m)[1]*dim(m)[2]*h_av
  }
  else {
    mi<-h
    h_av<-h
    e<-1
  }
  r<-NULL
  d<-DisPow(m)
  r<-ResPow(m)
  hs<-NULL
  pics<-NULL
  ds<- NULL
  rs<-NULL
  es<-NULL
  h_avs<-NULL
  mis<- NULL
  lBl<- length(Bl)
  if(lBl>1){
    BB<- c(0, Bl)
    hs<- numeric(lBl)
    pics<- numeric(lBl)
    ds<- numeric(lBl)
    es<- numeric(lBl)
    h_avs<- numeric(lBl)
    mis<- numeric(lBl)
    if (max(m)==1){rs<- numeric(lBl)}
    M<-m
    for (j in 1:lBl){
      m<- as.matrix(M[,(sum(BB[1:j])+1):(sum(BB[1:j])+BB[(j+1)])])
      alfbt<-unique(as.vector(m))
      n<- dim(m)[1]*dim(m)[2]
      prop<- numeric(length(alfbt))
      for (i in 1:length(prop)){
        prop[i]<- sum(alfbt[i]==m)/n
      }
      hs[j]<-1-sum(prop^2)
      p2<-(prop^2 %o% prop^2)
      ex<- 2*sum(p2[lower.tri(p2)])
      pics[j]<- h-ex
      ds[j]<- DisPow(m)
      if (max(M)==1){
        rs[j]<-ResPow(m)
        n1<- sum(m==1)
        beta<- n1/n
        es[j]<-dim(m)[2]*beta
        n2<- sum(m==0)#2
        h1<-1-(n1/n)^2
        h2<-1-(n2/n)^2
        h_avs[j]<- hs[j]/n#(h1+h2)/n
        mis[j]<- h1*es[j]#dim(m)[1]*dim(m)[2]*h_avs[j]#dim(m)[2]*h_av
      }
      else {
        mis[j]<-hs[j]
        h_avs[j]<-hs[j]
        es[j]<-1
      }
    }
  }
  H<- c(h, hs)
  PIC<- c(pic, pics)
  E<- c(e, es)
  H_av<- c(h_av, h_avs)
  MI<- c(mi, mis)
  D<- c(d, ds)
  R<- c(r, rs)
  H.name<-numeric(length(H))
  PIC.name<- numeric(length(PIC))
  E.name<- numeric(length(E))
  H_av.name<- numeric(length(H_av))
  MI.name<- numeric(length(MI))
  D.name<- numeric(length(D))
  if (is.null(R)){R.name<- NULL}
  else {R.name<- numeric(length(R))}
  for (i in 1:length(H)){
    H.name[i]<- paste("H_", i-1, sep="")
  }
  for (i in 1:length(PIC)){
    PIC.name[i]<- paste("PIC_", i-1, sep="")
  }
  for (i in 1:length(E)){
    E.name[i]<- paste("E_", i-1, sep="")
  }
  for (i in 1:length(H_av)){
    H_av.name[i]<- paste("H.av_", i-1, sep="")
  }
  for (i in 1:length(MI)){
    MI.name[i]<- paste("MI_", i-1, sep="")
  }
  for (i in 1:length(D)){
    D.name[i]<- paste("D_", i-1, sep="")
  }
  if (!is.null(R.name)){for (i in 1:length(R)){
    R.name[i]<- paste("R_", i-1, sep="")
  }}
  Value<- c(H, PIC, E, H_av, MI, D, R)
  Index<- c(H.name, PIC.name, E.name, H_av.name, MI.name, D.name, R.name)
  return(data.frame(Index, Value))
}
```