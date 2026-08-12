module Probability.Error
(ProbabilityError(..)) where

data ProbabilityError
    = ProbabilityOutOfRange Double
    | NegativeValue Double
    | InvalidBounds Double Double
    | InvalidSampleSize Int
    | InvalidSuccessCount Int
    | InvalidFailureCount Int
    | InvalidRate Double
    | InvalidScale Double
    | InvalidShape Double
    | InvalidPopulation Int
    |InvalidProbabilityVector [Double]
    |InvalidCountVector [Int]
    |DimensionMismatch Int Int
    |InvalidStandardDeviation Double
    deriving(Show, Eq)
