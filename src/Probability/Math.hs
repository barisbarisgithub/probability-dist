module Probability.Math
    ( logFactorial
    , logCombination
    , logGamma
    , logBeta
    ,erf
    ,expm1
    ) where

import Probability.Error (ProbabilityError(..))


logFactorial :: Int -> Either ProbabilityError Double
logFactorial n
    | n < 0  = Left (NegativeValue (fromIntegral n))
    | n < 35 = Right (sum [log (fromIntegral i) | i <- [1 .. n]])
    | otherwise = logGamma (fromIntegral n + 1)


logCombination :: Int -> Int -> Either ProbabilityError Double
logCombination n r
    | n < 0 = Left (NegativeValue (fromIntegral n))
    | r < 0 = Left (NegativeValue (fromIntegral r))
    | r > n = Left (InvalidSuccessCount r)
    | otherwise = do
        lnN  <- logFactorial n
        lnR  <- logFactorial r
        lnNR <- logFactorial (n - r)

        pure (lnN - lnR - lnNR)


logGamma :: Double -> Either ProbabilityError Double
logGamma z
    | z <= 0 = Left (NegativeValue z)
    | otherwise = Right result
  where
    g :: Double
    g = 7.0

    coefficients :: [Double]
    coefficients =
        [ 0.99999999999980993
        , 676.5203681218851
        , -1259.1392167224028
        , 771.32342877765313
        , -176.61502916214059
        , 12.507343278686905
        , -0.13857109526572012
        , 9.9843695780195716e-6
        , 1.5056327351493116e-7
        ]

    x :: Double
    x = z - 1.0

    a :: Double
    a =
        foldl
            (\acc (i, c) ->
                acc + c / (x + fromIntegral i))
            (head coefficients)
            (zip [1 ..] (tail coefficients))

    t :: Double
    t = x + g + 0.5

    result :: Double
    result =
        0.5 * log (2.0 * pi)
            + (x + 0.5) * log t
            - t
            + log a


logBeta :: Double -> Double -> Either ProbabilityError Double
logBeta a b
    | a <= 0 = Left (InvalidShape a)
    | b <= 0 = Left (InvalidShape b)
    | otherwise = do
        logA   <- logGamma a
        logB   <- logGamma b
        logAB  <- logGamma (a + b)

        pure (logA + logB - logAB)

erf
    :: Double
    -> Double
erf x
    | x == 0.0 = 0.0
    | otherwise =
        let sign = if x < 0.0 then -1.0 else 1.0
            ax   = abs x

            p    = 0.3275911
            t    = 1.0 / (1.0 + p * ax)

            poly =
                t * exp
                    ( -ax * ax
                        -1.26551223
                        + t * ( 1.00002368
                        + t * ( 0.37409196
                        + t * ( 0.09678418
                        + t * (-0.18628806
                        + t * (0.27886807
                        + t * (-1.13520398
                        + t * (1.48851587
                        + t * (-0.82215223
                        + t * 0.17087277)))))))))

        in sign * (1.0 - poly)

expm1
    :: Double
    -> Double
expm1 x
    | abs x < 1.0e-5 =
        x
        + x^2 / 2.0
        + x^3 / 6.0
        + x^4 / 24.0
        + x^5 / 120.0
        + x^6 / 720.0

    | otherwise =
        exp x - 1.0
