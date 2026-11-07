"""crosslabel: exact statistics for stochastically multiplexed localization microscopy.

Each molecule carries one of K labels drawn at random and kept for the whole
acquisition (for example DNA-PAINT docking sequences imaged in separate rounds,
as in RESI). Localizations or molecule positions with different labels then
come from different molecules, which turns cross-label pair and triplet
statistics into exact, clustering-free measures of molecular organization.

Main entry points
-----------------
cross_label_g           cross-label pair correlation with a shifted-label null
excess_neighbours       labelled partners within a radius (model-free)
connected_triplets_m2   labelled partner pairs per molecule; zero for any monomer/dimer mixture
triplet_angle_profile   angles between a molecule's two partners (assembly shape)
label_offsets           fiducial-free registration offsets between labels
assembly_stoichiometry  copy number and efficiency of discrete assemblies
design                  distinct-label probabilities and round-number planning
simulate                forward model for testing
"""
from . import design, labels, model, oracle, simulate
from .correlation import (ANGLE_EDGES, EDGES, TRIPLET_ANNULUS, TRIPLET_RADIUS, excess_neighbours,
                          excess_with_ci, label_translations, plateau_density, point_pair_g,
                          registration_offsets, triplet_angles, triplet_m2)
from .labels import label_moments
from . import stoichiometry
from .stoichiometry import group_counts, group_neighbour_statistics
from .stoichiometry import stoichiometry as assembly_stoichiometry

__version__ = '0.1.0'

cross_label_g = point_pair_g
connected_triplets_m2 = triplet_m2
triplet_angle_profile = triplet_angles
label_offsets = registration_offsets

__all__ = ['cross_label_g', 'excess_neighbours', 'excess_with_ci', 'connected_triplets_m2',
           'triplet_angle_profile', 'label_offsets', 'label_translations', 'plateau_density',
           'assembly_stoichiometry', 'stoichiometry', 'group_counts', 'group_neighbour_statistics', 'label_moments',
           'design', 'simulate', 'model', 'labels', 'oracle', 'EDGES', 'ANGLE_EDGES',
           'TRIPLET_RADIUS', 'TRIPLET_ANNULUS']
