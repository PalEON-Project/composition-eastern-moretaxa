# auxiliary code for creating the graph structure for binary (standard CAR-style) or thin-plate spline neighborhoods

graphCreate <- function(m1, m2, type = "bin", dir = getwd(), fn = NULL, fn_cats = NULL) {
  # this creates a file with nbhd info; cells are numbered from NW to NE and then in latitude/y rows
  # cells are numbered starting with 0

  m=m1*m2
  ids=0:(m-1)
  xgr=0:(m1-1)
  ygr=0:(m2-1)
  
  gr=expand.grid(xgr,rev(ygr))
  
  grid=0:(m-1)

  if(is.null(fn)) {
    typeUC = type
    substring(typeUC, 1 ,1) = toupper(substring(typeUC, 1, 1))
    fn = paste('graph', typeUC, '-',  m1, 'x', m2, '.csv', sep='')
  }

  if(is.null(fn_cats)) {
    typeUC = type
    substring(typeUC, 1 ,1) = toupper(substring(typeUC, 1, 1))
    fn_cats <- paste('graphCats', typeUC, '-', m1, 'x', m2, '.csv', sep='')
  }
  
  if(type == 'bin') {

    write(m,file.path(dir, fn))

    nbsx=c(0,-1,1,0)
    nbsy=c(-1,0,0,1)
    for(i in 1:m){
      x=gr[i,1]
      y=gr[i,2]
      id=i-1
      tmp=cbind(nbsx+x,nbsy+y)
      tmp=tmp[tmp[,1]>=0&tmp[,1]<m1&tmp[,2]>=0&tmp[,2]<m2,]
      nbs=sort((tmp[,1])+(m2-1 - tmp[,2])*m1)
      write(c(grid[i],length(nbs),nbs),file.path(dir, fn),ncol=length(nbs)+2,sep=' ',append=TRUE)
      fn_cats <- ""
    }
  } else { # tps
    write(m,file.path(dir, fn))

    nbsx=c(0,-1:1,-2,-1,1,2,-1:1,0)
    nbsy=c(-2,rep(-1,3),rep(0,4),rep(1,3),2)
    for(i in 1:m){
      x=gr[i,1]
      y=gr[i,2]
      id=i-1
      tmp=cbind(nbsx+x,nbsy+y)
      tmp=tmp[tmp[,1]>=0&tmp[,1]<m1&tmp[,2]>=0&tmp[,2]<m2,]
      nbs=(tmp[,1])+(m2-1 - tmp[,2])*m1
      write(c(grid[i],length(nbs),nbs),file.path(dir, fn),ncol=length(nbs)+2,sep=' ',append=TRUE)
    }

    category <- matrix(6,nr=m2,nc=m1)
    category[,2] <- 5
    category[,m1-1] <- 5
    category[2,] <- 5
    category[m2-1,] <- 5
    category[,1] <- 2
    category[,m1] <- 2
    category[1,] <- 2
    category[m2,] <- 2
    category[1,1] <- category[1,m1] <- category[m2,1] <- category[m2,m1] <- 1
    category[2,2] <- category[2,m1-1] <- category[m2-1,2] <- category[m2-1,m1-1] <- 4
    category[1,2] <- category[2,1] <- category[m2-1,1] <- category[m2,2] <- category[1,m1-1] <- category[2,m1] <- category[m2-1,m1] <- category[m2,m1-1] <- 3
    write(c(t(category)), file.path(dir, fn_cats) , ncol = 1)
  }
  return(c(fn, fn_cats))
}

graphRead <- function(fileName, catsFileName = "", m1, m2, type = 'bin', dir = getwd()) {
  # reads graph structure from file created by graphCreate() and produces sparse matrix precision matrix representing neighborhood structure
  require(spam)
  
  if(type %in% c("tps", "lindgren_nu1")){
    map=c(rep(0,4),4,0,10,11,0,18,19,20)
  } else{
    map=1:4
  }
  m = m1*m2
  
  Q=as.spam(diag(rep(1,2)))
  graphVals=as.integer(scan(file.path(dir, fileName)))
  if(type %in% c("tps", "lindgren_nu1"))
    cats=scan(file.path(dir, catsFileName))
  
  Q@dimension=c(graphVals[1],graphVals[1])
  dim=Q@dimension[1]
  rowpointers=rep(0,dim+1)
  
  graph=list()
  graph$n=graphVals[1]
  graph$cells=NULL
  graph$nnbs=NULL
  graph$nbs=NULL
  graph$ptr=NULL
  counter=1
  nbCounter=nbCounterGraph=1
  entries=colindices=NULL
  
  for(node in 0:(m-1)){
    nodex=node%%m1
    nodey=floor(node/m1)
    counter=counter+1
    graph$cells[node+1]=graphVals[counter]
    counter=counter+1
    nnbs=graph$nnbs[node+1]=graphVals[counter]
    rowpointers[node+1]=nbCounter
    graph$nbPtr[node+1]=nbCounterGraph
    
    if(nnbs){
      graph$nbs[(nbCounterGraph):(nbCounterGraph+graph$nnbs[node+1]-1)]=graphVals[(counter+1):(counter+graph$nnbs[node+1])]
      nbCounterGraph=nbCounterGraph+graph$nnbs[node+1]
      nbs=graphVals[(counter+1):(counter+nnbs)]
      nbs=sort(c(nbs,node))
      nbCounter=nbCounter+nnbs+1
    } else{
      nbs=node
    }
    counter=counter+nnbs
    
    nnodex=nbs%%m1
    nnodey=floor(nbs/m1)
    tmpent=NULL
    for(nnodeIndex in 1:length(nbs)){
      nnode=nbs[nnodeIndex]
      if(node==nnode){
        value=map[nnbs]
      } else{
        if(type %in% c('tps', 'lindgren_nu1')){
          if(abs(nodex-nnodex[nnodeIndex])==2 | abs(nodey-nnodey[nnodeIndex])==2){
            value=1
          } else{
            if(abs(nodex-nnodex[nnodeIndex])==1 & abs(nodey-nnodey[nnodeIndex])==1){
              value=2
            }
            else{
              value=-8
              if(cats[node+1]==2 | cats[node+1]==3){
                value=-6
              }
              if(cats[node+1]==1){
                value=-4
              }
              if((cats[node+1]==3 & cats[nnode+1]==1) | (cats[node+1]==4 & cats[nnode+1]==3) |(cats[node+1]==5 & cats[nnode+1]==2)){
                value=value+2
              }
            }
          }
        } else{
          if(abs(nodex-nnodex[nnodeIndex])+abs(nodey-nnodey[nnodeIndex])==1){
            value=-1
          }
        }
      }
      tmpent=c(tmpent,value)
    }
    entries=c(entries,tmpent)
    colindices=c(colindices,nbs+1)
    if((node+1)%%1000==0){print(node+1)}
  }
  rowpointers[dim+1]=length(entries)+1
  Q@rowpointers=as.integer(rowpointers)
  Q@colindices=as.integer(colindices)
  Q@entries=entries
  return(Q)
}
