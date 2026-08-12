module Probability.Continuous
(
    normalPDF
    ,normalLogPDF
    ,normalCDF
    ,normalMean
    ,normalVar
    ,exponentialLogPDF
    ,exponentialPDF
    ,exponentialCDF
    ,exponentialMean
    ,exponentialVar
    ,gammaLogPDF
    ,gammaPDF
    ,gammaMean
    ,gammaVar
    ,uniformLogPDF
    ,uniformPDF
    ,uniformCDF
    ,uniformMean
    ,uniformVar
)where

import Probability.Error (ProbabilityError(..))
import Probability.Math
    (logFactorial
    , logCombination
    , logGamma
    , logBeta
    ,erf
    ,expm1)
normalLogPDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
normalLogPDF x mu sigma
    | sigma <= 0 =
        Left (InvalidStandardDeviation sigma)

    | otherwise =
        let z = (x - mu) / sigma

            logP =
                - log sigma
                - 0.5 * log (2.0 * pi)
                - 0.5 * z * z

        in Right logP


normalPDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
normalPDF x mu sigma
    | sigma <= 0 =
        Left (InvalidStandardDeviation sigma)

    | otherwise = do
        logP <- normalLogPDF x mu sigma
        pure (exp logP)


normalCDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
normalCDF x mu sigma
    | sigma <= 0 =
        Left (InvalidStandardDeviation sigma)

    | otherwise =
        let z = (x - mu) / (sigma * sqrt 2.0)

        in Right
            (0.5 * (1.0 + erf z))


normalMean
    :: Double
    -> Double
    -> Either ProbabilityError Double
normalMean mu sigma
    | sigma <= 0 =
        Left (InvalidStandardDeviation sigma)

    | otherwise =
        Right mu


normalVar
    :: Double
    -> Double
    -> Either ProbabilityError Double
normalVar mu sigma
    | sigma <= 0 =
        Left (InvalidStandardDeviation sigma)

    | otherwise =
        Right (sigma * sigma)


exponentialLogPDF
    :: Double
    -> Double
    -> Either ProbabilityError Double
exponentialLogPDF x lambda
    | lambda <= 0 =
        Left (InvalidRate lambda)

    | x < 0 =
        Right (-1.0 / 0.0)

    | otherwise =
        Right (log lambda - lambda * x)


exponentialPDF
    :: Double
    -> Double
    -> Either ProbabilityError Double
exponentialPDF x lambda
    | lambda <= 0 =
        Left (InvalidRate lambda)

    | x < 0 =
        Right 0.0

    | otherwise = do
        logP <- exponentialLogPDF x lambda
        pure (exp logP)


exponentialCDF
    :: Double
    -> Double
    -> Either ProbabilityError Double
exponentialCDF x lambda
    | lambda <= 0 =
        Left (InvalidRate lambda)

    | x < 0 =
        Right 0.0

    | otherwise =
        Right (- expm1 (-lambda * x))


exponentialMean
    :: Double
    -> Either ProbabilityError Double
exponentialMean lambda
    | lambda <= 0 =
        Left (InvalidRate lambda)

    | otherwise =
        Right (1.0 / lambda)


exponentialVar
    :: Double
    -> Either ProbabilityError Double
exponentialVar lambda
    | lambda <= 0 =
        Left (InvalidRate lambda)

    | otherwise =
        Right (1.0 / (lambda * lambda))

gammaLogPDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
gammaLogPDF x alpha lambda
    | alpha <= 0 =
        Left (InvalidShape alpha)

    | lambda <= 0 =
        Left (InvalidRate lambda)

    | x < 0 =
        Right (-1.0 / 0.0)

    | x == 0 && alpha < 1 =
        Right (1.0 / 0.0)

    | x == 0 && alpha == 1 =
        Right (log lambda)

    | x == 0 =
        Right (-1.0 / 0.0)

    | otherwise = do
        logG <- logGamma alpha

        let logP =
                alpha * log lambda
                - logG
                + (alpha - 1.0) * log x
                - lambda * x

        pure logP


gammaPDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
gammaPDF x alpha lambda
    | alpha <= 0 =
        Left (InvalidShape alpha)

    | lambda <= 0 =
        Left (InvalidRate lambda)

    | x < 0 =
        Right 0.0

    | x == 0 && alpha < 1 =
        Right (1.0 / 0.0)

    | x == 0 && alpha == 1 =
        Right lambda

    | x == 0 =
        Right 0.0

    | otherwise = do
        logP <- gammaLogPDF x alpha lambda
        pure (exp logP)


gammaMean
    :: Double
    -> Double
    -> Either ProbabilityError Double
gammaMean alpha lambda
    | alpha <= 0 =
        Left (InvalidShape alpha)

    | lambda <= 0 =
        Left (InvalidRate lambda)

    | otherwise =
        Right (alpha / lambda)


gammaVar
    :: Double
    -> Double
    -> Either ProbabilityError Double
gammaVar alpha lambda
    | alpha <= 0 =
        Left (InvalidShape alpha)

    | lambda <= 0 =
        Left (InvalidRate lambda)

    | otherwise =
        Right (alpha / (lambda * lambda))

uniformLogPDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
uniformLogPDF x a b
    | b <= a =
        Left (InvalidBounds a b)

    | x < a || x > b =
        Right (-1.0 / 0.0)

    | otherwise =
        Right (- log (b - a))


uniformPDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
uniformPDF x a b
    | b <= a =
        Left (InvalidBounds a b)

    | x < a || x > b =
        Right 0.0

    | otherwise =
        Right (1.0 / (b - a))


uniformCDF
    :: Double
    -> Double
    -> Double
    -> Either ProbabilityError Double
uniformCDF x a b
    | b <= a =
        Left (InvalidBounds a b)

    | x < a =
        Right 0.0

    | x > b =
        Right 1.0

    | otherwise =
        Right ((x - a) / (b - a))


uniformMean
    :: Double
    -> Double
    -> Either ProbabilityError Double
uniformMean a b
    | b <= a =
        Left (InvalidBounds a b)

    | otherwise =
        Right ((a + b) / 2.0)


uniformVar
    :: Double
    -> Double
    -> Either ProbabilityError Double
uniformVar a b
    | b <= a =
        Left (InvalidBounds a b)

    | otherwise =
        let width = b - a
        in Right (width * width / 12.0)
