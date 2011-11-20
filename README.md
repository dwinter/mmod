*Modern Measures of Differentiation

This is an R-module providing functions to calculate measures of 
genetic differentiation. At the moment this README is designed to help
me lay out the work that needs to be done. 

At the moment, there are top-level function for three stats

+D.Jost()
+Gst.Nei()
+Gst.Hedrick()

These functions replicate each others work quite a lot so
there is another to do all three at once:

+diff.stats()

Finally, there are a couple of ways to check how robust the results are.
Pairwise examples of each measure:

+pairwise.D()
+pairwise.Gst.Net()
+pairwise.Gst.Hedrick()

And a way of taking jack-knife samples across populations:

+jacknife.pop()

**TODO
DOCUMENTATION!!
Possibly bring the pairwise functions together into one in the style of diff.stats()












