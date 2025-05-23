context("Typical session")

test_that("Typical session can be run.", {
  set.seed(0)
  neg=100/365
  off.r=5
  w.shape=10
  w.scale=0.1
  pi=0.25
  dateT=2008
    expect_silent(simu <- simulateOutbreak(pi=pi,off.r=off.r,w.shape=w.shape,
                             w.scale=w.scale,dateStartOutbreak=2005,dateT=dateT))
   expect_is(simu,'ctree')
   expect_silent(plotCTree(simu))
   expect_is(capture_output(print(simu)),'character')
   expect_silent(plot(simu))
   expect_silent(ttree<-extractTTree(simu))
   expect_is(ttree,'ttree')
   expect_is(capture_output(print(ttree)),'character')
   expect_silent(plot(ttree,w.shape,w.scale))
   expect_silent(ptree<-extractPTree(simu))
   expect_silent(plot(ptree))
   expect_is(ptree,'ptree')
   expect_is(capture_output(print(ptree)),'character')
   expect_silent(p<-phyloFromPTree(ptree))
   expect_is(p,'phylo')
   expect_silent(ptree<-ptreeFromPhylo(p,dateLastSample=max(simu$ctree[,1])))
   expect_equal(dateLastSample(simu),dateLastSample(ttree))
   expect_equal(dateLastSample(ptree),dateLastSample(ttree))
   expect_is(ptree,'ptree')
   expect_is(capture_output(res<-inferTTree(ptree,mcmcIterations=100,w.shape=w.shape,w.scale=w.scale,dateT=dateT,verbose=T)),'character')
   expect_is(res,'resTransPhylo')
   expect_is(capture_output(print(res)),'character')
   expect_is(capture_output(summary(res)),'character')
   expect_silent(last<-extractCTree(res,100))
   expect_is(last,'ctree')
   expect_silent(plot(last))
   expect_silent(plot(res))
   expect_silent(plotTraces(res,burnin=0.5))
   expect_silent(mcmc<-convertToCoda(res))
   expect_is(mcmc,'mcmc')
   expect_is(capture_output(cons<-consTTree(res)),'character')
   expect_is(cons,'ctree')
   expect_silent(plot(cons))
   expect_silent(ttree<-extractTTree(cons))
   expect_is(ttree,'ttree')
   expect_silent(plot(ttree))
   expect_silent(plot(ttree,type='detailed',w.shape,w.scale))
   expect_silent(mat<-computeMatWIW(res))
   expect_is(mat,'matrix')
   expect_silent(mat<-computeMatTDist(res))
   expect_is(mat,'matrix')
   expect_silent(a<-getIncidentCases(res,show.plot = T))
   expect_is(a,'list')
   expect_silent(a<-getGenerationTimeDist(res,show.plot = T))
   expect_is(a,'matrix')
   expect_silent(a<-getSamplingTimeDist(res,show.plot = T))
   expect_is(a,'matrix')
   expect_silent(a<-getInfectionTimeDist(res,k='1.1',show.plot = T))
   expect_is(a,'numeric')
   expect_silent(a<-getOffspringDist(res,k='1.1',show.plot = T))
   expect_is(a,'numeric')
})




