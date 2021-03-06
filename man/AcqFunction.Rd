% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AcqFunction.R
\name{AcqFunction}
\alias{AcqFunction}
\title{Acquisition Function Base Class}
\description{
Based on a surrogate model, the acquisition function encodes the preference to evaluate
a new point for evaluation.
}
\seealso{
Other Acquisition Function: 
\code{\link{AcqFunctionAEI}},
\code{\link{AcqFunctionCB}},
\code{\link{AcqFunctionEIPS}},
\code{\link{AcqFunctionEI}},
\code{\link{AcqFunctionMESN}},
\code{\link{AcqFunctionSmsEgo}}
}
\concept{Acquisition Function}
\section{Super classes}{
\code{\link[bbotk:Objective]{bbotk::Objective}} -> \code{\link[bbotk:ObjectiveRFunDt]{bbotk::ObjectiveRFunDt}} -> \code{AcqFunction}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{surrogate}}{\link{Surrogate}.}

\item{\code{direction}}{(\code{character(1)})\cr
Must be \code{"same"}, \code{"minimize"}, or \code{"maximize"}.}

\item{\code{surrogate_max_to_min}}{(\code{numeric(1)}).
Optimization direction of the objective function: 1 for minimization, -1 for maximization.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{AcqFunction$new()}}
\item \href{#method-setup}{\code{AcqFunction$setup()}}
\item \href{#method-update}{\code{AcqFunction$update()}}
\item \href{#method-clone}{\code{AcqFunction$clone()}}
}
}
\if{html}{
\out{<details open ><summary>Inherited methods</summary>}
\itemize{
\item \out{<span class="pkg-link" data-pkg="bbotk" data-topic="Objective" data-id="eval">}\href{../../bbotk/html/Objective.html#method-eval}{\code{bbotk::Objective$eval()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="bbotk" data-topic="Objective" data-id="format">}\href{../../bbotk/html/Objective.html#method-format}{\code{bbotk::Objective$format()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="bbotk" data-topic="Objective" data-id="print">}\href{../../bbotk/html/Objective.html#method-print}{\code{bbotk::Objective$print()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="bbotk" data-topic="ObjectiveRFunDt" data-id="eval_dt">}\href{../../bbotk/html/ObjectiveRFunDt.html#method-eval_dt}{\code{bbotk::ObjectiveRFunDt$eval_dt()}}\out{</span>}
\item \out{<span class="pkg-link" data-pkg="bbotk" data-topic="ObjectiveRFunDt" data-id="eval_many">}\href{../../bbotk/html/ObjectiveRFunDt.html#method-eval_many}{\code{bbotk::ObjectiveRFunDt$eval_many()}}\out{</span>}
}
\out{</details>}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of this \link[R6:R6Class]{R6} class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AcqFunction$new(id, constants = ParamSet$new(), surrogate, direction, fun)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{id}}{(\code{character(1)}).}

\item{\code{constants}}{(\link[paradox:ParamSet]{paradox::ParamSet}).}

\item{\code{surrogate}}{\link{Surrogate}.}

\item{\code{direction}}{(\code{character(1)}).}

\item{\code{fun}}{(\verb{function(xdt)}).
Must be \code{"same"}, \code{"minimize"}, or \code{"maximize"}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setup"></a>}}
\if{latex}{\out{\hypertarget{method-setup}{}}}
\subsection{Method \code{setup()}}{
Sets up the acquisition function.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AcqFunction$setup(archive)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{archive}}{\link[bbotk:Archive]{bbotk::Archive}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-update"></a>}}
\if{latex}{\out{\hypertarget{method-update}{}}}
\subsection{Method \code{update()}}{
Update the acquisition function given an \link[bbotk:Archive]{bbotk::Archive}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AcqFunction$update(archive)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{archive}}{\link[bbotk:Archive]{bbotk::Archive}.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{AcqFunction$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
