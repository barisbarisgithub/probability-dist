module Main where

import Control.Monad (unless)
import Data.Either (isLeft)
import System.Exit (exitFailure)

--import Probability.Math
import Probability.Discrete
import Probability.Continuous


-- ============================================================
-- Test helpers
-- ============================================================

epsilon :: Double
epsilon = 1.0e-12


approxEqual :: Double -> Double -> Bool
approxEqual x y =
    abs (x - y) < epsilon


assert :: Bool -> String -> IO ()
assert condition message =
    unless condition $ do
        putStrLn ("FAIL: " ++ message)
        exitFailure


assertApprox :: Double -> Double -> String -> IO ()
assertApprox actual expected message =
    assert (approxEqual actual expected) message


assertRightApprox
    :: Either e Double
    -> Double
    -> String
    -> IO ()
assertRightApprox result expected message =
    case result of
        Right value ->
            assertApprox value expected message
        Left _ ->
            assert False (message ++ " returned Left")


assertLeft
    :: Either e a
    -> String
    -> IO ()
assertLeft result message =
    assert (isLeft result) message

{-
-- ============================================================
-- Math
-- ============================================================

testMath :: IO ()
testMath = do
    putStrLn "Testing Math..."

    -- logGamma
    assertApprox
        (logGamma 4)
        1.791759469228055
        "logGamma 4"

    assertApprox
        (logGamma 5)
        3.178053830347944
        "logGamma 5"

    assertApprox
        (logGamma 10)
        12.801827480081474
        "logGamma 10"

    -- logFactorial
    assertApprox
        (logFactorial 0)
        0.0
        "logFactorial 0"

    assertApprox
        (logFactorial 10)
        15.104412573075516
        "logFactorial 10"

    assertApprox
        (logFactorial 40)
        110.32063971475739
        "logFactorial 40"

    -- logCombination
    assertRightApprox
        (logCombination 10 0)
        0.0
        "logCombination 10 0"

    assertRightApprox
        (logCombination 10 10)
        0.0
        "logCombination 10 10"

    assertLeft
        (logCombination 10 11)
        "logCombination 10 11"

    assertLeft
        (logCombination (-1) 2)
        "logCombination -1 2"

    -- logBeta
    assertApprox
        (logBeta 2 2)
        (-1.791759469228055)
        "logBeta 2 2"

    assertApprox
        (logBeta 1 2)
        (-0.6931471805599472)
        "logBeta 1 2"

    -- erf
    assertApprox
        (erf 0.0)
        0.0
        "erf 0"

    -- expm1
    assertApprox
        (expm1 0.0)
        0.0
        "expm1 0"

    assertApprox
        (expm1 1.0)
        (exp 1.0 - 1.0)
        "expm1 1"

    putStrLn "Math: OK"

-}
-- ============================================================
-- Discrete
-- ============================================================

testDiscrete :: IO ()
testDiscrete = do
    putStrLn "Testing Discrete..."

    -- --------------------------------------------------------
    -- Binomial
    -- --------------------------------------------------------

    assertRightApprox
        (binomialPMF 0 10 0.5)
        0.0009765625
        "binomialPMF 10 0 0.5"

    assertRightApprox
        (binomialPMF 10 10 0.5)
        0.0009765625
        "binomialPMF 10 10 0.5"

    assertLeft
        (binomialPMF 11 10 0.5)
        "binomialPMF invalid success count"

    assertLeft
        (binomialPMF (-1) 2 0.5)
        "binomialPMF negative population"

    -- --------------------------------------------------------
    -- Geometric
    -- --------------------------------------------------------

    assertRightApprox
        (geometricPMF 1 0.5)
        0.25
        "geometricPMF 1 0.5"

    -- --------------------------------------------------------
    -- Negative Binomial
    -- --------------------------------------------------------

    assertRightApprox
        (negativeBinomialPMF 0 1 0.5)
        0.5
        "negativeBinomialPMF r=1"

    -- --------------------------------------------------------
    -- Multinomial
    -- --------------------------------------------------------

    assertRightApprox
        (multinomialPMF
            10
            [10, 0]
            [1.0, 0.0])
        1.0
        "multinomialPMF degenerate case"

    -- --------------------------------------------------------
    -- Poisson
    -- --------------------------------------------------------

    assertRightApprox
        (poissonPMF 0 2.0)
        (exp (-2.0))
        "poissonPMF 0 2"

    -- --------------------------------------------------------
    -- Bernoulli
    -- --------------------------------------------------------

    assertRightApprox
        (bernoulliPMF 1 0.7)
        0.7
        "bernoulliPMF success"

    assertRightApprox
        (bernoulliPMF 0 0.7)
        0.3
        "bernoulliPMF failure"

    putStrLn "Discrete: OK"


-- ============================================================
-- Continuous
-- ============================================================

testContinuous :: IO ()
testContinuous = do
    putStrLn "Testing Continuous..."

    -- --------------------------------------------------------
    -- Normal
    -- --------------------------------------------------------

    assertRightApprox
        (normalPDF 0.0 0.0 1.0)
        0.3989422804014327
        "normalPDF standard normal at 0"

    assertRightApprox
        (normalCDF 0.0 0.0 1.0)
        0.5
        "normalCDF standard normal at 0"

    -- Symmetry
    let leftCDF  = normalCDF (-1.0) 0.0 1.0
        rightCDF = normalCDF 1.0 0.0 1.0

    case (leftCDF, rightCDF) of
        (Right l, Right r) ->
            assertApprox
                (l + r)
                1.0
                "normalCDF symmetry"
        _ ->
            assert False "normalCDF symmetry returned Left"

    -- Invalid standard deviation
    assertLeft
        (normalPDF 0.0 0.0 (-1.0))
        "normalPDF invalid standard deviation"

    -- --------------------------------------------------------
    -- Exponential
    -- --------------------------------------------------------

    assertRightApprox
        (exponentialPDF 0.0 2.0)
        2.0
        "exponentialPDF x=0"

    assertRightApprox
        (exponentialCDF 0.0 2.0)
        0.0
        "exponentialCDF x=0"

    assertRightApprox
        (exponentialCDF 1.0 2.0)
        (1.0 - exp (-2.0))
        "exponentialCDF x=1 lambda=2"

    assertLeft
        (exponentialPDF 1.0 (-1.0))
        "exponentialPDF invalid lambda"

    -- --------------------------------------------------------
    -- Gamma
    -- --------------------------------------------------------

    -- Gamma(1, lambda) = Exponential(lambda)
    assertRightApprox
        (gammaPDF 0.0 1.0 2.0)
        2.0
        "gammaPDF alpha=1 x=0"

    assertRightApprox
        (gammaPDF 1.0 2.0 1.0)
        0.36787944117144233
        "gammaPDF alpha=2 lambda=1 x=1"

    assertRightApprox
        (gammaMean 2.0 1.0)
        2.0
        "gammaMean"

    assertRightApprox
        (gammaVar 2.0 1.0)
        2.0
        "gammaVar"

    assertLeft
        (gammaPDF 1.0 (-1.0) 1.0)
        "gammaPDF invalid shape"

    assertLeft
        (gammaPDF 1.0 1.0 (-1.0))
        "gammaPDF invalid rate"

    -- --------------------------------------------------------
    -- Uniform
    -- --------------------------------------------------------

    assertRightApprox
        (uniformPDF 0.5 0.0 1.0)
        1.0
        "uniformPDF U(0,1)"

    assertRightApprox
        (uniformCDF 0.5 0.0 1.0)
        0.5
        "uniformCDF U(0,1)"

    assertRightApprox
        (uniformMean 0.0 1.0)
        0.5
        "uniformMean U(0,1)"

    assertRightApprox
        (uniformVar 0.0 1.0)
        (1.0 / 12.0)
        "uniformVar U(0,1)"

    -- CDF boundaries
    assertRightApprox
        (uniformCDF (-1.0) 0.0 1.0)
        0.0
        "uniformCDF below lower bound"

    assertRightApprox
        (uniformCDF 2.0 0.0 1.0)
        1.0
        "uniformCDF above upper bound"

    -- Invalid bounds
    assertLeft
        (uniformPDF 0.0 1.0 1.0)
        "uniformPDF invalid bounds"

    assertLeft
        (uniformCDF 0.0 2.0 1.0)
        "uniformCDF invalid bounds"

    putStrLn "Continuous: OK"


-- ============================================================
-- Main
-- ============================================================

main :: IO ()
main = do
    putStrLn "========================================"
    putStrLn " probability-dist test suite"
    putStrLn "========================================"

    --testMath
    testDiscrete
    testContinuous

    putStrLn "========================================"
    putStrLn " All tests passed."
    putStrLn "========================================"
