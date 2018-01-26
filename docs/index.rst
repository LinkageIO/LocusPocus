.. locuspocus documentation master file, created by
   sphinx-quickstart on Tue Nov  7 15:11:02 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _overview:
.. currentmodule:: locuspocus

#########################################################
LocusPocus: genetic coordinates so easy, it's like magic!
#########################################################


**locuspocus** is a `Python <http://www.python.org>`__ package for creating, 
storing, manipulating, and comparing genetic coordinates, i.e. Loci. Locuspocus
subscribes to the *build once, use many times* philosophy. The LocusPocus 
library provies an way to build single use loci objects, *or* persistant, named
sets of Loci which are *freezalbe* and backed by the *minus80* python library.


To get started, import the locuspocus module

.. ipython:: python
    
    import locuspocus as lp


Since LocusPocus is backed by *minus80*, it has a similar persistance model. 
Before diving into the locuspocus documentation, it might be useful to refresh
on how *minus80* models data.


LocuPocus offers access to the following objects:

- Locus (a genetic coordinate)
- RefLoci (a *freezable* set of Reference Loci)
- Term (a subset of related loci, e.g. GO/KEGG/MapMan Terms)
- Ontology (a *freezable* set of Reference Terms)

The differneces between these objects is important but subtle. As always, perhaps the best
way to understand how these work together is through an example.  






Two main objects
are available: a :class:`Locus` class that represents a genetic coordinate, and
a :class:`RefLoci` object that represents a set of Loci which have a collective
meaning. 





Table of Contents
=================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Overview <index>


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
