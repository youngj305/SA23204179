% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/homework.R
\name{solve_game}
\alias{solve_game}
\title{The Morra game's solution}
\usage{
solve_game(A)
}
\arguments{
\item{A}{a matrix.}
}
\value{
a list to store the solution of the Morra game.
}
\description{
In the Morra game, the set of optimal strategies are not changed if a constant is subtracted from every entry of the payoff matrix, or a positive constant is multiplied times every entry of the payoff matrix. However, the simplex algorithm may terminate at a different basic feasible point (also optimal).
}
\examples{
\dontrun{
    A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
    2,0,0,0,-3,-3,4,0,0,
    2,0,0,3,0,0,0,-4,-4,
    -3,0,-3,0,4,0,0,5,0,
    0,3,0,-4,0,-4,0,5,0,
    0,3,0,0,4,0,-5,0,-5,
    -4,-4,0,0,0,5,0,0,6,
    0,0,4,-5,-5,0,0,0,6,
    0,0,4,0,0,5,-6,-6,0), 9, 9)
    result <- solve_game(A)
}
}
