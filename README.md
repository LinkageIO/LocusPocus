# LocusPocus - Genetic coordinates so simple, it feels like MAGIC!


![locus pocus use cases](img/UseCases.png)

A sequenced reference genome establishes a genetic coordinate system in which
sequencing technologies can relate various genomic features. However, the
genetic coordinate unit, or locus, is not well represented by common data
structures available in contemporary programming languages, often leading to
error-prone, ad-hoc implementations. Here, we implement the Locus as a
fundamental, object-oriented datatype in pure-Python, which enables convenient
relational algebra including powerful expressions to be built which make sense
in the scope of loci. 

Operations such as calculating distance between two
loci or whether they overlap are achieved by overloading the “in” (x in y) and
“subtraction” (x - y) Python mathematical operators making comparisons
intuitive; relative genomic positions and equality can be tested using built-in
Python “greater-than” (x > y), “less-than” (x < y), and “equals” (x == y)
operators. Extending these basic comparison methods, a nested sub_locus feature
is implemented by overloading the “addition” operator (x + y) allowing the
construction of more complicated genomic features such as combining SNPs which
have overlapping, user-defined windows. More specific features such as genes
and entire chromosomes can be easily implemented by extending the base Locus
class and adding feature-specific logic accordingly. Implementation of a
genetic coordinate class in Python allows for quick and easy execution of
common locus related tasks. The extensive overloading of built-in Python
operators makes code more readable and less prone to bugs. The latest version
of LocusPocus can be found on PyPi. Upstream development can be accessed via
GitHub: https://github.com/schae234/LocusPocus.

LocusPocus Model
----------------
![locuspocus model](img/model.png)


License
-------
LocusPocus is freely available under the MIT license, see LICENSE for more info
