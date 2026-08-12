module Probability.Discrete
    ( bernoulliPMF
    , bernoulliMean
    , bernoulliVar
    ,binomialPMF
    ,binomialLogPMF
    ,binomialMean
    ,binomialVar
    ,poissonPMF
    ,poissonLogPMF
    ,poissonMean
    ,poissonVar
    ,negativeBinomialPMF
    ,negativeBinomialLogPMF
    ,negativeBinomialMean
    ,negativeBinomialVar
    ,geometricLogPMF
    ,geometricPMF
    ,geometricMean
    ,geometricVar
    ,hypergeometricPMF
    ,hypergeometricLogPMF
    ,hypergeometricMean
    ,hypergeometricVar
    ,multinomialLogPMF
    ,multinomialPMF
    ) where

import Probability.Error (ProbabilityError(..))
import Probability.Math
    (logCombination
    ,logFactorial)

bernoulliPMF :: Int -> Double -> Either ProbabilityError Double
bernoulliPMF x p
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | x == 0         = Right (1 - p)
    | x == 1         = Right p
    | otherwise      = Right 0.0


bernoulliMean :: Double -> Either ProbabilityError Double
bernoulliMean p
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | otherwise      = Right p


bernoulliVar :: Double -> Either ProbabilityError Double
bernoulliVar p
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | otherwise      = Right (p * (1 - p))



binomialLogPMF
    :: Int
    -> Int
    -> Double
    -> Either ProbabilityError Double
binomialLogPMF k n p
    | n < 0 = Left (InvalidSampleSize n)
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | k < 0 || k > n = Right 0.0
    | otherwise = do
        logC <- logCombination n k

        let termK =
                if k == 0
                    then 0.0
                    else fromIntegral k * log p

            termNMinusK =
                if k == n
                    then 0.0
                    else fromIntegral (n - k) * log (1 - p)

        return (logC + termK + termNMinusK)


binomialPMF
    :: Int
    -> Int
    -> Double
    -> Either ProbabilityError Double
binomialPMF k n p
    | n < 0 = Left (InvalidSampleSize n)
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | k < 0 || k > n = Left(InvalidSuccessCount k)
    | otherwise = do
        logP <- binomialLogPMF k n p
        pure (exp logP)

binomialMean
    :: Int
    -> Double
    -> Either ProbabilityError Double
binomialMean n p
    | n < 0 = Left (InvalidSampleSize n)
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | otherwise = Right (fromIntegral n * p)

binomialVar
    :: Int
    -> Double
    -> Either ProbabilityError Double
binomialVar n p
    | n < 0 = Left (InvalidSampleSize n)
    | p < 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | otherwise = Right (fromIntegral n * p * (1 - p))


poissonLogPMF
    :: Int
    -> Double
    -> Either ProbabilityError Double
poissonLogPMF k lambda
    | lambda <= 0 = Left (InvalidRate lambda)
    | k < 0       = Right 0.0
    | otherwise   = do
        logFact <- logFactorial k

        let logP =
                -lambda
                + if k == 0
                    then 0.0
                    else fromIntegral k * log lambda
                - logFact

        pure logP


poissonPMF
    :: Int
    -> Double
    -> Either ProbabilityError Double
poissonPMF k lambda
    | lambda <= 0 = Left (InvalidRate lambda)
    | k < 0       = Right 0.0
    | otherwise   = do
        logP <- poissonLogPMF k lambda
        pure (exp logP)


poissonMean
    :: Double
    -> Either ProbabilityError Double
poissonMean lambda
    | lambda <= 0 = Left (InvalidRate lambda)
    | otherwise   = Right lambda


poissonVar
    :: Double
    -> Either ProbabilityError Double
poissonVar lambda
    | lambda <= 0 = Left (InvalidRate lambda)
    | otherwise   = Right lambda


negativeBinomialLogPMF
    :: Int
    -> Int
    -> Double
    -> Either ProbabilityError Double
negativeBinomialLogPMF k r p
    | r <= 0 = Left (InvalidSuccessCount r)
    | p <= 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | k < 0 = Left (InvalidFailureCount k)
    | otherwise = do
        logC <- logCombination (r + k - 1) k

        let logP =
                logC
                + fromIntegral r * log p
                + if k == 0
                    then 0.0
                    else fromIntegral k * log (1 - p)

        pure logP


negativeBinomialPMF
    :: Int
    -> Int
    -> Double
    -> Either ProbabilityError Double
negativeBinomialPMF k r p
    | r <= 0 = Left (InvalidSuccessCount r)
    | p <= 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | k < 0 = Left (InvalidFailureCount k)
    | otherwise = do
        logP <- negativeBinomialLogPMF k r p
        pure (exp logP)


negativeBinomialMean
    :: Int
    -> Double
    -> Either ProbabilityError Double
negativeBinomialMean r p
    | r <= 0 = Left (InvalidSuccessCount r)
    | p <= 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | otherwise =
        Right (fromIntegral r * (1 - p) / p)


negativeBinomialVar
    :: Int
    -> Double
    -> Either ProbabilityError Double
negativeBinomialVar r p
    | r <= 0 = Left (InvalidSuccessCount r)
    | p <= 0 || p > 1 = Left (ProbabilityOutOfRange p)
    | otherwise =
        Right (fromIntegral r * (1 - p) / (p * p))


geometricLogPMF
    :: Int
    -> Double
    -> Either ProbabilityError Double
geometricLogPMF k p =
    negativeBinomialLogPMF k 1 p


geometricPMF
    :: Int
    -> Double
    -> Either ProbabilityError Double
geometricPMF k p =
    negativeBinomialPMF k 1 p


geometricMean
    :: Double
    -> Either ProbabilityError Double
geometricMean p =
    negativeBinomialMean 1 p


geometricVar
    :: Double
    -> Either ProbabilityError Double
geometricVar p =
    negativeBinomialVar 1 p


hypergeometricLogPMF
    :: Int
    -> Int
    -> Int
    -> Int
    -> Either ProbabilityError Double
hypergeometricLogPMF k sampleSize successes populationSize
    | populationSize <= 0 =
        Left (InvalidPopulation populationSize) --size

    | successes < 0 || successes > populationSize =
        Left (InvalidSuccessCount successes)

    | sampleSize < 0 || sampleSize > populationSize =
        Left (InvalidSampleSize sampleSize)

    | k < 0
        || k > sampleSize
        || k > successes
        || sampleSize - k > populationSize - successes =
        Right 0.0

    | otherwise = do
        logC1 <- logCombination successes k
        logC2 <- logCombination
                    (populationSize - successes)
                    (sampleSize - k)
        logC3 <- logCombination populationSize sampleSize

        pure (logC1 + logC2 - logC3)


hypergeometricPMF
    :: Int
    -> Int
    -> Int
    -> Int
    -> Either ProbabilityError Double
hypergeometricPMF k sampleSize successes populationSize
    | populationSize <= 0 =
        Left (InvalidPopulation populationSize) --size

    | successes < 0 || successes > populationSize =
        Left (InvalidSuccessCount successes)

    | sampleSize < 0 || sampleSize > populationSize =
        Left (InvalidSampleSize sampleSize)

    | k < 0
        || k > sampleSize
        || k > successes
        || sampleSize - k > populationSize - successes =
        Right 0.0

    | otherwise = do
        logP <- hypergeometricLogPMF
                    k
                    sampleSize
                    successes
                    populationSize

        pure (exp logP)


hypergeometricMean
    :: Int
    -> Int
    -> Int
    -> Either ProbabilityError Double
hypergeometricMean sampleSize successes populationSize
    | populationSize <= 0 =
        Left (InvalidPopulation populationSize) -- size kısımını düzelttik

    | successes < 0 || successes > populationSize =
        Left (InvalidSuccessCount successes)

    | sampleSize < 0 || sampleSize > populationSize =
        Left (InvalidSampleSize sampleSize)

    | otherwise =
        Right
            ( fromIntegral sampleSize
            * fromIntegral successes
            / fromIntegral populationSize
            )


hypergeometricVar
    :: Int
    -> Int
    -> Int
    -> Either ProbabilityError Double
hypergeometricVar sampleSize successes populationSize
    | populationSize <= 1 =
        Left (InvalidPopulation populationSize)   -- düzeltme yaptık

    | successes < 0 || successes > populationSize =
        Left (InvalidSuccessCount successes)

    | sampleSize < 0 || sampleSize > populationSize =
        Left (InvalidSampleSize sampleSize)

    | otherwise =
        let n = fromIntegral sampleSize
            k = fromIntegral successes
            nPop = fromIntegral populationSize

            p = k / nPop

        in Right
            ( n
            * p
            * (1 - p)
            * ((nPop - n) / (nPop - 1))
            )

multinomialLogPMF
    :: Int
    -> [Int]
    -> [Double]
    -> Either ProbabilityError Double
multinomialLogPMF n counts probabilities
    | n < 0 =
        Left (InvalidSampleSize n)

    | null probabilities =
        Left (InvalidProbabilityVector probabilities)

    | length counts /= length probabilities =
        Left (DimensionMismatch
                (length counts)
                (length probabilities))

    | any (< 0) counts =
        Left (InvalidCountVector counts)

    | sum counts /= n =
        Left (InvalidSampleSize n)

    | any (< 0) probabilities =
        Left (InvalidProbabilityVector probabilities)

    | abs (sum probabilities - 1.0) > 1.0e-12 =
        Left (InvalidProbabilityVector probabilities)

    | otherwise = do
        logNFact <- logFactorial n

        logCountsFact <- mapM logFactorial counts

        let logCountTerm =
                logNFact - sum logCountsFact

            logProbabilityTerm =
                sum
                    [ if k == 0
                        then 0.0
                        else fromIntegral k * log p
                    | (k, p) <- zip counts probabilities
                    ]

        pure (logCountTerm + logProbabilityTerm)


multinomialPMF
    :: Int
    -> [Int]
    -> [Double]
    -> Either ProbabilityError Double
multinomialPMF n counts probabilities
    | otherwise = do
        logP <- multinomialLogPMF n counts probabilities
        pure (exp logP)
