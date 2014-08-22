module HaskSee.ImageLibrary
(
    applyFilter, 
    colorToGray,
    medianFilter,
    Roberts(AbsoluteValue, SquareRoot),
    convolutionFilter,
    sobelRowFilter,
    sobelColFilter,
    laplaceFilter,
    meanFilter,
    distanceBetween,
    thresholdFilter,
    writeImage,
    readImage,
    Feature,
    FeatureValue,
    rowFeature,
    colFeature,
    rectFeature,
    circularFeature,
    perimeterFeature,
    centerSquareFeature,
) where

import Data.Image.Internal hiding (medianFilter)
import Data.Image.Boxed hiding (medianFilter)
import Data.List (sort)

applyFilter :: (Image i) => ([[Pixel i]] -> Pixel i) -> Int -> Int -> i -> i
applyFilter filter width height image 
    --First we need to make sure that the user is not wanting a filter that is of even size, cause that doesn't really make sense
    | (width `mod` 2 == 0) || (height `mod` 2 == 0) = error "You cannot have a filter that is an even size"
    --So we make an image of the same size that consists of the filter function supplied the 2D array of pixels that surround the row and column that we are currently looking at
    | otherwise = makeImage (rows image) (cols image) (\r c -> filter (pixels r c))
    where pixels row col = map (map (\(c,r) -> safeRef image r c)) $ filterPixelLocations row col
          filterPixelLocations row col = 
                [[(x,y) | 
                    x <- [(col - (div (width+1) 2) + 1)..(col + (div (width+1) 2) - 1)]] | 
                    y <- [(row - (div (height+1) 2) + 1)..(row + (div (height+1) 2) - 1)]]

safeRef :: (Image i) => i -> Int -> Int -> Pixel i
safeRef image r c = ref image row col
    where row = if r < 0 then 0 else (if r >= imageRows then imageRows - 1 else r)
          col = if c < 0 then 0 else (if c >= imageCols then imageCols - 1 else c)
          imageRows = rows image
          imageCols = cols image


colorToGray :: ColorImage -> GrayImage
colorToGray image = makeImage (rows image) (cols image) (\r c -> grayValue (ref image r c))
    where grayValue (RGB (r, g, b)) = (r+g+b)/3

medianFilter :: [[Pixel GrayImage]] -> Pixel GrayImage
--First we flatten out the pixels and the sort it and take the middle element
medianFilter array = (!! ((length flattened) `div` 2)) . sort $ flattened
    where flattened = foldr (++) [] array

data Roberts = SquareRoot | AbsoluteValue

aMinusD :: [[Pixel GrayImage]] -> Pixel GrayImage
aMinusD array = ((array !! 1) !! 1) - ((array !! 0) !! 0)

bMinusC :: [[Pixel GrayImage]] -> Pixel GrayImage
bMinusC array = ((array !! 0) !! 1) - ((array !! 1) !! 0)

robertsFilter :: Roberts -> [[Pixel GrayImage]] -> Pixel GrayImage
robertsFilter SquareRoot array = sqrt (((aMinusD array)*(aMinusD array)) + ((bMinusC array)*(bMinusC array)))
robertsFilter AbsoluteValue array = abs (aMinusD array) + abs (bMinusC array)


flatten :: [[a]] -> [a]
flatten = foldr (++) []

convolutionFilter :: [[Pixel GrayImage]] -> [[Pixel GrayImage]] -> Pixel GrayImage
convolutionFilter filter array = sum . map (\(f, i) -> f*i) . zip (flatten filter) . flatten $ array

sobelRowFilter :: [[Pixel GrayImage]]
sobelRowFilter = [[-1,0,1],
                  [-2,0,2],
                  [-1,0,1]]

sobelColFilter :: [[Pixel GrayImage]]
sobelColFilter = [[-1,-2,-1],
                  [0,0,0],
                  [1,2,1]]

laplaceFilter :: [[Pixel GrayImage]]
laplaceFilter = [[1,-2,1],
                 [-2,5,-2],
                 [1,-2,1]]

numLength :: (Num b) => [a] -> b
numLength [] = 0
numLength (_:tail) = 1 + (numLength tail)

meanFilter3x3 :: [[Pixel GrayImage]]
meanFilter3x3 = [repeat (1/9)]

meanFilter :: [[Pixel GrayImage]] -> Pixel GrayImage
meanFilter pixels = fromIntegral . flip div (numLength allPix) . sum . map round $ allPix
    where allPix = flatten pixels 

thresholdFilter :: Pixel GrayImage -> [[Pixel GrayImage]] -> Pixel GrayImage
thresholdFilter threshold ((pixel:_):_) = if pixel > threshold then 255 else 0 

distanceBetween :: GrayImage -> GrayImage -> GrayImage
distanceBetween rowImage colImage = makeImage (rows rowImage) (cols colImage) (\r c -> distance (ref rowImage r c) (ref colImage r c))
    where distance x y = sqrt ((x*x)+(y*y))

type Feature i = i -> FeatureValue i
type FeatureValue i = [Pixel i]

rowFeature :: (Image i) => Int -> Int -> Int -> Feature i
rowFeature row col len image = map (\index -> ref image index col) [row..(row+len-1)]

colFeature :: (Image i) => Int -> Int -> Int -> Feature i
colFeature row col len image = map (\index -> ref image row index) [col..(col+len-1)]

rectFeature :: (Image i) => Int -> Int -> Int -> Int -> Feature i
rectFeature row col width height image = flatten . map (\colIndex -> rowFeature row colIndex width image) [col..(col+height-1)]

circularFeature :: (Image i) => Int -> Int -> Int -> Feature i
circularFeature row col radius image = map (\(r, c) -> ref image r c) . generateCirclePoints (row, col) $ radius

perimeterFeature :: (Image i) => Int -> Int -> Int -> Int -> Feature i
perimeterFeature row col width height image = topLine ++ leftLine ++ botLine ++ rightLine
    where topLine = rowFeature row col width image 
          botLine = rowFeature row (col+height-1) width image
          leftLine = colFeature row (col + 1) (height - 2) image
          rightLine = colFeature (row + width - 1) (col + 1) (height - 2) image

centerSquareFeature :: (Image i) => Int -> Feature i
centerSquareFeature size image = rectFeature row col size size image 
    where row = ((rows image) `div` 2) - (size `div` 2)
          col = ((cols image) `div` 2) - (size `div` 2)

allPixelsFeature :: (Image i) => Feature i
allPixelsFeature = pixelList

allRowsFeature :: (Image i) => Feature i
allRowsFeature = allPixelsFeature

--This is, in essence, the same as the row feature although the pixels are in a different order
allColsFeature :: (Image i) => Feature i
allColsFeature = allPixelsFeature . transpose

{-
    This is not my code, this is just a haskell implementation taken from the rosetta code
    This is an implementation of the Midpoint Circle algorithm, that is used for finding
    all points in a circle given a center point and a radius.
-}

type Point = (Int, Int)

-- Takes the center of the circle and radius, and returns the circle points
generateCirclePoints :: Point -> Int -> [Point]
  -- Four initial points, plus the generated points
generateCirclePoints (x0, y0) radius = (x0, y0 + radius) : (x0, y0 - radius) : (x0 + radius, y0) : (x0 - radius, y0) : points
    where
      -- Creates the (x, y) octet offsets, then maps them to absolute points in all octets.
      points = concatMap generatePoints $ unfoldr step initialValues
      generatePoints (x, y) = [(xop x0 x', yop y0 y') | (x', y') <- [(x, y), (y, x)], xop <- [(+), (-)], yop <- [(+), (-)]]

      -- The initial values for the loop
      initialValues = (1 - radius, 1, (-2) * radius, 0, radius)

      -- One step of the loop. The loop itself stops at Nothing.
      step (f, ddf_x, ddf_y, x, y) | x >= y = Nothing
                                   | otherwise = Just ((x', y'), (f', ddf_x', ddf_y', x', y'))
                                     where
                                       (f', ddf_y', y') | f >= 0 = (f + ddf_y' + ddf_x', ddf_y + 2, y - 1)
                                                        | otherwise = (f + ddf_x, ddf_y, y)
                                       ddf_x' = ddf_x + 2
                                       x' = x + 1
