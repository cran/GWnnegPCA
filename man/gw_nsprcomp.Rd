\name{gw_nsprcomp}
\alias{gw_nsprcomp}
\title{Geographically Weighted non-negative Principal Component Analysis}
\usage{
gw_nsprcomp(data, elocat, vars, bw, k = 2, kernel, adaptive = TRUE,
          p = 2, theta = 0, longlat = FALSE, dMat = NULL, n.obs = NA,
          n.iter = 1, ncomp = k, nneg = TRUE, localcenter = TRUE, localscale = FALSE,...)
}
\description{
Implementation of geographically weighted non-negative principal component analysis, whcih consists of the fusion of GWPCA and sparse non-negative PCA.
}

\arguments{
  \item{data}{
a Spatial*DataFrame either SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp.

}
  \item{elocat}{
Same as GWmodel::gwpca. two-column numeric array or Spatial*DataFrame object for providing evaluation locations, i.e. SpatialPointsDataFrame or SpatialPolygonsDataFrame as defined in package sp.

}
  \item{vars}{
the number of retained components; k must be less than the number of variables.

}

  \item{bw}{
bandwidth used in the weighting function, possibly calculated by bw.gwpca;fixed (distance) or adaptive bandwidth(number of nearest neighbours).
}


  \item{k}{
the number of retained components; k must be less than the number of variables
}

  \item{kernel}{
Same as GWmodel::gwpca. Function chosen as follows:
gaussian: wgt = exp(-.5*(vdist/bw)^2);
exponential: wgt = exp(-vdist/bw);
bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
boxcar: wgt=1 if dist < bw, wgt=0 otherwise

see help(GWmodel::gw.weight) more detail.
}
 \item{adaptive}{
 if TRUE calculate an adaptive kernel where the bandwidth corresponds to the number of nearest neighbours (i.e. adaptive distance); default is FALSE, where a fixed kernel is found (bandwidth is a fixed distance).

}


  \item{p}{
the power of the Minkowski distance, default is 2, i.e. the Euclidean distance.

}
  \item{theta}{
an angle in radians to rotate the coordinate system, default is 0.

}
  \item{longlat}{
if TRUE, great circle distances will be calculated.

}
  \item{dMat}{
a pre-specified distance matrix, it can be calculated by the function gw.dist .

}
  \item{n.obs}{
Number of observations used to find the correlation matrix if using a correlation matrix. Used for finding the goodness of fit statistics. Must be specified if using a correlaton matrix and finding confidence intervals.

}

  \item{n.iter}{
Same as psych::fa. Number of bootstrap interations to do in fa or fa.poly

}
  \item{ncomp}{
the number of principal components (PCs) to be computed. With the default setting, PCs are computed until x is fully deflated. ncomp can be specified implicitly if k is given as a vector.

}

  \item{nneg}{

a logical value indicating whether the loadings should be non-negative, i.e. the PAs should be constrained to the non-negative orthant.

}

  \item{localcenter}{

If TRUE, local weighted x is centered. The default is TRUE.

}

  \item{localscale}{

If TRUE, local weighted x is scaled. The default is FALSE.

}


  \item{...}{
arguments passed to or from other methods.

}

}

\value{
  \item{loadings }{The localized loadings}
  \item{score }{The PC score by the localized non-negative PCA.}
  \item{sdev }{The localized standard deviation of the principal components.}
}
\references{
Tsutsumida N., Murakami D.,Yoshida T., Nakaya T. Exploring geographically weighted non negative principal component analysis for producing index. The 27th meeting of GIS association of Japan, Tokyo, 20-21 October, 2018 (Japanese), http://www.gisa-japan.org/conferences/proceedings/2018/papers/C52.pdf

Tsutsumida N., Murakami D., Yoshida T., Nakaya T., Lu B., and P. Harris. Geographically Weighted Non-negative Principal Component Analysis for Exploring Spatial Variation in Multidimensional Composite Index, Geocomputation 2019, https://doi.org/10.17608/k6.auckland.9850826.v1


}
\author{
N. Tsutsumida
}

\examples{

### This example implements the GWnnegPCA for the same data used in the GWmodel example.
\donttest{
library(GWmodel)
data(DubVoter)
Data.scaled <- scale(as.matrix(Dub.voter@data[, 4:11]))
Coords <- as.matrix(cbind(Dub.voter$X, Dub.voter$Y))
Data.scaled.spdf <- SpatialPointsDataFrame(Coords, as.data.frame(Data.scaled))

bw_gwnnegpca_ans <- 200

gwnnegpca_ans <- gw_nsprcomp(
  data = Data.scaled.spdf,
  vars = colnames(Data.scaled.spdf@data),
  bw = bw_gwnnegpca_ans,
  k = 8,
  kernel="bisquare",
  adaptive = TRUE,
  nneg=TRUE,
  center=FALSE
)

Dub.voter$gwnnegpcaload_DiffAdd <- gwnnegpca_ans$loadings[,"DiffAdd","PC1"]
Dub.voter$gwnnegpcaload_LARent <- gwnnegpca_ans$loadings[,"LARent","PC1"]
Dub.voter$gwnnegpcaload_SC1 <- gwnnegpca_ans$loadings[,"SC1","PC1"]

spplot(Dub.voter[c("gwnnegpcaload_DiffAdd","gwnnegpcaload_LARent","gwnnegpcaload_SC1")])

### Compare with standard gwpca
bw_gwpca_ans <- bw.gwpca(
  Data.scaled.spdf,
  vars = colnames(Data.scaled.spdf@data),
  k = 3,
  robust = FALSE,
  adaptive = TRUE
)

gwpca_ans <- gwpca(
  Data.scaled.spdf,
  vars = colnames(Data.scaled.spdf@data),
  bw = bw.gwpca.basic,
  k = 8,
  robust = FALSE,
  adaptive = TRUE
)

Dub.voter$gwpcaload_DiffAdd <- gwpca_ans$loadings[,"DiffAdd","PC1"]
Dub.voter$gwpcaload_LARent <- gwpca_ans$loadings[,"LARent","PC1"]
Dub.voter$gwpcaload_SC1 <- gwpca_ans$loadings[,"SC1","PC1"]

spplot(Dub.voter[c("gwpcaload_DiffAdd","gwpcaload_LARent","gwpcaload_SC1")])

}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "SPATIAL" }% use one of  RShowDoc("KEYWORDS")
\keyword{ "MULTIVALIATE" }% __ONLY ONE__ keyword per line
