---
title: "Regression modelling using I-priors"
subtitle: NUS Department of Statistics & Data Science Seminar
author: "Haziq Jamil"
date: "Wednesday, 16 November 2022"
institute: |
  | Mathematical Sciences, Faculty of Science, UBD
  | \url{https://haziqj.ml}
output: 
  beamer_presentation:
    template: ubd_beamer_rmd.tex
    latex_engine: xelatex
    slide_level: 3
    keep_tex: false
    citation_package: biblatex
    pandoc_args: ["--lua-filter=/Library/Frameworks/R.framework/Versions/4.2/Resources/library/bookdown/rmarkdown/lua/custom-environment.lua"]    
    # includes:
    #   after_body: afterbody.txt
toc: false
toctitle: "Overview"
banner: true
logo: true
progressdots: true
transitions: true
handout: false
bibliography: refs.bib
nocite: '@*'
refslide: false
aspectratio: 43
editor_options: 
  markdown: 
    wrap: 72
# header-includes:
#   - \usetikzlibrary{backgrounds,calc,intersections}
#   - \usepackage{pgfplots}
---




### Abstract

Regression analysis is undoubtedly an important tool to understand the relationship between one or more explanatory and independent variables of interest. 
The problem of estimating a generic regression function in a model with normal errors is considered.
For this purpose, a novel objective prior for the regression function is proposed, defined as the distribution maximizing entropy (subject to a suitable constraint) based on the Fisher information on the regression function.
This prior is called the I-prior.
The regression function is then estimated by its posterior mean under the I-prior, and accompanying hyperparameters are estimated via maximum marginal likelihood. 
Estimation of I-prior models is simple and inference straightforward, while predictive performances are comparative, and often better, to similar leading state-of-the-art models--as will be illustrated by several data examples.
Further plans for research in this area are also presented, including variable selection for interaction effects and extending the I-prior methodology to non-Gaussian errors.
Please visit the project website for further details: https://phd.haziqj.ml/

\vspace{1em}
Keywords: Bayes, Fisher information, RKHS, Gaussian process, EM algorithm

### Plan

- Introduction
- Some basic functional analysis (?)
- The I-prior 
- Estimation
- Inference
- Examples
- Further work (variable selection, interaction effects, non-gaussian errors)

# Introduction


For $i=1,\dots,n$, consider the regression model

\begin{equation}\label{mod1}
\begin{gathered}
y_i = f(x_i) + \epsilon_i \\
(\epsilon_1,\dots,\epsilon_n)^\top \sim \N_n(0, \Psi^{-1})
\end{gathered}
\end{equation}

where each $y_i\in\bbR$, $x_i\in \cX$ (some set of covariates), and $f$ is a regression function.
This forms the basis for a multitude of statistical models:

1. Ordinary linear regression when $f$ is parameterised linearly.
<!-- $f(x_i) = x_i^\top\beta$, with $\cX,\beta\in\bbR^p$. -->

2. Varying intercepts/slopes model when $\cX$ is grouped.
<!-- $f(x_{ij},j) = f_1(x_i) + f_2(j) + f_{12}(x_{ij},j)$. -->

3. Smoothing models when $f$ is a smooth function.

4. Functional regression when $\cX$ is functional.



::: {.block}
#### Goal
To estimate the regression function $f$ given the observations $\{(y_i,x_i)\}_{i=1}^n$.
:::

### Ordinary linear regression

Suppose $f(x_i) = x_i^\top \beta$ for $i=1,\dots,n$, where $x_i,\beta \in \bbR^p$.

\vspace{1em}

![](figure/linear-reg-1.pdf)<!-- --> 

### Varying intercepts/slopes model

Suppose each unit $i=1,\dots,n$ relates to the $k$th observation in group $j\in\{1,\dots,m\}$. 
Model the function $f$ additively:
\only<1>{
$$
f(x_{kj}, j) = f_1(x_{kj}) + f_2(j) + f_{12}(x_{kj},j).
\phantom{\myunderbrace{x_{kj}^\top\beta_1}{f_1}}
$$
}
\only<2>{
$$
f(x_{kj}, j) = 
\myunderbrace{x_{kj}^\top\beta_1}{f_1} + 
\myunderbrace{\beta_{0j}}{f_2} + 
\myunderbrace{x_{kj}^\top\beta_{1j}}{f_{12}}
\phantom{f_1(x_{kj})}
$$
}
\vspace{-0.5em}

::: {.onlyenv latex=<1>}

![](figure/var-int-slope1-1.pdf)<!-- --> 

:::

::: {.onlyenv latex=<2>}

![](figure/var-int-slope2-1.pdf)<!-- --> 

:::

### Smoothing models

Suppose $f\in\cF$ where $\cF$ is a space of "smoothing functions" (models like LOESS, kernel regression, smoothing splines, etc.).

\vspace{1em}

![](figure/smooth1-1.pdf)<!-- --> 

### Functional regression

Suppose the input set $\cX$ is functional. 
The (linear) regression aims to estimate a coefficient function $\beta:\cT\to\bbR$
$$
y_i = \myunderbrace{\int_\cT x_i(t)\beta(t) \dint t}{f(x_i)} + \epsilon_i
$$

![](figure/functionalx-1.pdf)<!-- --> 

### The I-prior 

For the regression model stated in \eqref{mod1}, we assume that $f$ lies in some RKHS of functions $\cF$, with reproducing kernel $h$ over $\cX$.

::: {.definition name="I-prior"} 
\label{def:iprior}
The entropy maximising prior distribution for $f$, subject to constraints, is \vspace{-0.7em}
\begin{equation}
\begin{gathered}\label{iprior}
f(x) = \sum_{i=1}^n h(x,x_i)w_i \\
(w_1,\dots,w_n)^\top \sim \N_n(0, \Psi) \vspace{-0.7em}
\end{gathered}
\end{equation}
:::

Therefore, the covariance kernel of $\mathbf f = \big(f(x_1),\dots, f(x_n) \big)^\top$ is determined by the function \vspace{-0.4em}
$$
k(x,x') = \sum_{i=1}^n\sum_{j=1}^n \Psi_{i,j} h(x,x_i )h(x',x_j),
$$
which happens to be **Fisher information** between two linear forms of $f$.

### The I-prior (cont.)

Interpretation:

<!-- > The more information about $f$, the larger its prior variance, and hence the smaller the influence of the prior mean (and vice versa). -->

\vspace{-1em}

\begin{empheq}[box=\tcbhighmath]{align*}
\begin{split}
&\color{navyblue}\text{The more information about } f, \text{ the larger its prior variance,} \\[-0.2em]
&\color{navyblue}\text{and hence the smaller the influence of the prior mean (and} \\[-0.2em]
&\color{navyblue}\text{vice versa).}
\end{split}
\end{empheq}

\pause

\vspace{0.2em}

Of interest then are

1. Posterior distribution for the regression function,
$$
p(\mathbf f \,|\, \mathbf y) =
\frac{p(\mathbf y \,|\, \mathbf f)p(\mathbf f)}
{\int p(\mathbf y \,|\, \mathbf f)p(\mathbf f) \dint \mathbf f}.
$$


2. Posterior predictive distribution (given a new data point $x_{new}$)
$$
p(y_{new} \,|\, \mathbf y) = \int p( y_{new} \,|\, f_{new}) p( f_{new} \,|\, \mathbf y) \dint f_{new},
$$
where $f_{new} = f(x_{new})$.

### Introduction (cont.)

Observations $\{(y_i,x_i) \mid y_i,x_i\in\bbR \ \forall i=1,\dots,n\}$.

\vspace{1em }

![](figure/datapoints-1.pdf)<!-- --> 

### Introduction (cont.)

Choose $h(x,x') = e^{-\frac{\lVert x - x' \rVert^2}{2l^2}}$ (Gaussian kernel).
Sample paths from the I-prior:

\vspace{1em}

![](figure/priorsamp-1.pdf)<!-- --> 

### Introduction (cont.)

Sample paths from the posterior of $f$:

\vspace{1em}

![](figure/postsamp-1.pdf)<!-- --> 

### Introduction (cont.)

Posterior mean estimate for $y=f(x)$ and its 95% credibility interval.

\vspace{1em}

![](figure/ipriorplot-1.pdf)<!-- --> 

### Why I-priors?




::: {.columns}

::: {.column width=48%}


Advantages

- Provides a unifying methodology for regression.

- Simple and parsimonious model specification and estimation.

- Often yield comparable (or better) predictions than competing ML algorithms.

<!-- Competitors -->

<!-- - Tikhonov regulariser (e.g. cubic spline smoother) -->

<!-- - Gaussian process regression -->

:::

::: {.column width=48%}

\includegraphics[width=0.9\linewidth]{figure/wordcloud} 
:::

:::

Competitors:

- Tikhonov regulariser (e.g. cubic spline smoother)
$$
\hat f = \argmin_f \sum_{i=1}^n (y_i - f(x_i))^2 + \lambda \int f''(x)^2 \dint x
$$

- Gaussian process regression

# Regression using I-priors

## Reproducing kernel Hilbert spaces



## The Fisher information

Suppose further that $f\in\cF$ where $\cF$ is a reproducing kernel Hilbert space (RKHS) with reproducing kernel $h:\cX\times\cX\to\bbR$.
Then \eqref{mod1} can be expressed as 
\begin{equation}\label{mod2}
\begin{gathered}
y_i = \big\langle f, h(\cdot,x_i)\big\rangle_\cF + \epsilon_i \\
(\epsilon_1,\dots,\epsilon_n)^\top \sim \N(\bzero, \bPsi^{-1})
\end{gathered}
\end{equation}

The Fisher information for $f$ is given by
$$
\cI_f = \sum_{i=1}^n\sum_{j=1}^n \psi_{ij}h(\cdot,x_i) \otimes h(\cdot,x_j)
$$


It's helpful to think of $\cI_f$ as a bilinear form $\cI_f:\cF \times \cF \to \bbR$ defined by
$$
\cI_f = -\E\nabla^2 L(f|y)
$$
so between two linear functionals of f....

###

where each $y_i\in\bbR$,, and $f\in\cF$ a reproducing kernel Hilbert space (RKHS) with kernel $h:\cX\times\cX\to\bbR$.
The I-prior [@bergsma2019] for the regression function $f$ is the random function defined
\begin{equation}
\begin{gathered}
f(x_i) = f_0(x_i) + \sum_{k=1}^n h(x_i,x_k)w_k \\
(w_1,\dots,w_n)^\top \sim \N(\bzero, \bPsi)
\end{gathered}
\end{equation}
where $f_0$ is some prior mean for the regression function.


<!-- ### Advantages and current state -->

<!-- Advantages -->



<!-- \pause -->

<!-- Current state of research -->

<!-- - **1 x PhD Thesis**: @jamil2018phdthesis  -->
<!--     - Theoretical foundations -->
<!--     - Computational methods -->
<!--     - Bayesian variable selection -->
<!--     - \includegraphics[height=.75em]{1F3C6} 2020 Zellner Thesis Award (honourable mention) -->
<!-- - **3 x manuscripts**: @bergsma2019, @bergsma2020regression, and @jamil4bayesian. -->
<!-- - **1 x `R` package**: @jamil2017iprior -->



# Estimation

# Examples

# Further research

Hello







