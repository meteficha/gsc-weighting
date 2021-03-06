module Data.Weighting.GSC where

import Data.List (foldl')
import Data.Clustering.Hierarchical (Dendrogram(..), Distance)

-- | /O(n^2)/ Calculates the Gerstein\/Sonnhammer\/Chothia
-- weights for all elements of a dendrogram.  Weights are
-- annotated to the leafs of the dendrogram while distances in
-- branches are kept unchanged.
--
-- Distances @d@ in branches should be non-increasing and between
-- @0@ (in the leafs) and @1@.  The final weights are normalized
-- to average to @1@ (i.e. sum to the number of sequences, the
-- same they would sum if all weights were @1@).
--
-- For example, suppose we have
--
-- @
-- dendro = Branch 0.8
--            (Branch 0.5
--              (Branch 0.2
--                (Leaf 'A')
--                (Leaf 'B'))
--              (Leaf 'C'))
--            (Leaf 'D')
-- @
--
-- This is the same as GSC paper's example, however they used
-- similarities while we are using distances (i.e. applying
-- @(1-)@ to the distances would give exactly their example).
-- Then @gsc dendro@ is
--
-- @
-- gsc dendro == Branch 0.8
--                 (Branch 0.5
--                   (Branch 0.2
--                     (Leaf ('A',0.7608695652173914))
--                     (Leaf ('B',0.7608695652173914)))
--                   (Leaf ('C',1.0869565217391306)))
--                 (Leaf ('D',1.3913043478260871))
-- @
--
-- which is exactly what they calculated.
gsc :: Dendrogram a -> Dendrogram (a, Distance)
gsc (Leaf x)                  = Leaf (x,1)
gsc b@(Branch d _ _) | d == 0 = fmap (flip (,) 1) b
gsc dendrogram                = ret
    where
      (wsumF, nsumF, ret) = go undefined [] dendrogram
      wfinal              = (wsumF / fromIntegral nsumF)

      position (Leaf _)       = 0 -- no difference from itself
      position (Branch d _ _) = d

      go d' cs b@(Branch d _ _) | d == 0 = zeroBranch d' cs b
      go _  cs   (Branch d l r) =
          let (wl, nl, l') = go d ((el / wl) : cs) l
              (wr, nr, r') = go d ((er / wr) : cs) r

              el = d - position l -- edge length to left branch
              er = d - position r --          ...to right branch

              wsum = wl + wr + el + er
              nsum = nl + nr
          in wsum `seq` nsum `seq` (wsum, nsum, Branch d l' r')
      go d cs (Leaf x) =
          -- O(n) worst case, O(log n) best case (balanced dendrogram)
          let w = foldl' (\curw c -> curw + curw * c) d (tail cs)
          in (0, 1 :: Int, Leaf (x, w / wfinal))

      -- special case for branches where the distance is zero.
      zeroBranch d' cs b =
          let w  = foldl' (\curw c -> curw + curw * c) d' (tail cs)
              wf = w / (wfinal * fromIntegral nsum)
              setWeights (Leaf x) = (1, Leaf (x, wf))
              setWeights (Branch d l r) =
                  let (nl, l') = setWeights l
                      (nr, r') = setWeights r
                      ns = nl + nr
                  in ns `seq` (ns, Branch d l' r')
              (nsum, retz) = setWeights b
          in (0, nsum, retz)