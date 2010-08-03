module Data.Weighting.GSC where

import Data.List (foldl')
import Data.Clustering.Hierarchical (Dendrogram(..))

-- | /O(n^2)/ Calculates the Gerstein\/Sonnhammer\/Chothia
-- weights for all elements of a dendrogram.  Weights are
-- annotated to the leafs of the dendrogram while distances in
-- branches are kept unchanged.
--
-- Distances @d@ in branches should be non-increasing and between
-- @0@ (in the leafs) and @1@.  The final weights are normalized
-- to average to @1@ (i.e. sum to the number of sequences, the
-- same they would sum if all weights were @1@).
gsc :: Fractional d => Dendrogram d a -> Dendrogram d (a, d)
gsc (Leaf x) = Leaf (x,1)
gsc dendro   = let (wsum, nsum, r) = go wfinal undefined [] dendro
                   wfinal = (wsum / fromIntegral nsum)
               in r
    where
      position (Leaf _)       = 0 -- no difference from itself
      position (Branch d _ _) = d

      go wfinal _ cs (Branch d l r) =
          let (wl, nl, l') = go wfinal d ((el / wl) : cs) l
              (wr, nr, r') = go wfinal d ((er / wr) : cs) r

              el = d - position l -- edge length to left branch
              er = d - position r --          ...to right branch

              wsum = wl + wr + el + er
              nsum = nl + nr
          in wsum `seq` nsum `seq` (wsum, nsum, Branch d l' r')
      go wfinal d cs (Leaf x) =
          -- O(n) worst case, O(log n) best case (balanced dendrogram)
          let w = foldl' (\curw c -> curw + curw * c) d (tail cs)
          in (0, 1 :: Int, Leaf (x, w / wfinal))
