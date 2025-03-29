#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Scoring and Filtering Module
============================
This module contains functions for scoring and filtering binding site predictions.
"""

import numpy as np
import logging
from scipy.stats import zscore
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist

# Use the same logger as main to prevent duplicate messages
logger = logging.getLogger('ConSBind')

def final_scoring(consensus_pockets, protein_type):
    """
    Adjust pocket scores based on known protein function and automatically filter
    significant pockets without requiring a threshold parameter.
    
    Parameters:
    -----------
    consensus_pockets : list
        List of pocket dictionaries
    protein_type : str
        Type of protein ('enzyme', 'transporter', 'receptor', or 'unknown')
        
    Returns:
    --------
    list
        List of pocket dictionaries with adjusted scores and filtered to significant pockets
    """
    if not consensus_pockets:
        logger.warning("No binding sites found")
        return []
        
    # Apply protein-type specific adjustments
    for pocket in consensus_pockets:
        # For enzymes, prioritize pockets with catalytic residue patterns 
        if protein_type == 'enzyme':
            if 'knowledge_score' in pocket and pocket['knowledge_score'] > 1.5:
                pocket['consensus_score'] *= 1.3    # Boost enzymatic sites 
        
        # For transporters, prioritize channel-like cavities 
        elif protein_type == 'transporter':
            if pocket['size'] > 300:    # Larger pockets for transporters 
                pocket['consensus_score'] *= 1.2
        
        # For receptors, prioritize larger, moderately hydrophobic pockets 
        elif protein_type == 'receptor':
            if 'druggability' in pocket and pocket['druggability'] > 0.6:
                pocket['consensus_score'] *= 1.15
        
        # Ensure all pockets have the methods list if it doesn't exist
        if 'methods' not in pocket:
            pocket['methods'] = [pocket['method']] if 'method' in pocket else []
    
    # Ensure all pockets have final score attribute
    for pocket in consensus_pockets:
        pocket['final_score'] = (
            pocket['consensus_score'] * 3.0 +   # Primary criterion: consensus
            pocket.get('score', 0) * 0.5        # Secondary criterion: knowledge/druggability score
        )
    
    # Sort by consensus score and then by final score
    consensus_pockets = sorted(consensus_pockets, key=lambda x: (x['consensus_score'], x['final_score']), reverse=True)
    
    # Automatic filtering based on statistical properties of the scores
    if len(consensus_pockets) > 1:
        # Extract consensus scores
        scores = np.array([p['consensus_score'] for p in consensus_pockets])
        
        # Method 1: Use natural breaks in the data (Jenks Natural Breaks optimization)
        # Simplified version using hierarchical clustering
        if len(scores) >= 3:
            try:
                # Reshape for clustering
                score_matrix = scores.reshape(-1, 1)
                # Calculate distance matrix
                distances = pdist(score_matrix)
                # Perform hierarchical clustering
                Z = linkage(distances, method='ward')
                # Find optimal number of clusters (2-4 clusters)
                max_clusters = min(4, len(scores))
                best_clusters = 2
                best_score = -np.inf
                
                for n_clusters in range(2, max_clusters + 1):
                    clusters = fcluster(Z, n_clusters, criterion='maxclust')
                    # Calculate cluster separation metric
                    cluster_means = []
                    for i in range(1, n_clusters + 1):
                        cluster_scores = scores[clusters == i]
                        if len(cluster_scores) > 0:  # Check if cluster has elements
                            cluster_means.append(np.mean(cluster_scores))
                        else:
                            cluster_means.append(0)  # Assign 0 to empty clusters
                            
                    if len(cluster_means) > 1:
                        sorted_means = sorted(cluster_means, reverse=True)
                        separation = sorted_means[0] - sorted_means[1]  # Explicit calculation instead of np.diff
                        if separation > best_score:
                            best_score = separation
                            best_clusters = n_clusters
                
                # Get the clusters
                clusters = fcluster(Z, best_clusters, criterion='maxclust')
                # Find the cluster with the highest mean score
                cluster_means = []
                for i in range(1, best_clusters + 1):
                    cluster_scores = scores[clusters == i]
                    if len(cluster_scores) > 0:  # Check if cluster has elements
                        cluster_means.append((i, np.mean(cluster_scores)))
                    else:
                        cluster_means.append((i, 0))  # Assign 0 to empty clusters
                
                if cluster_means:  # Make sure we have clusters
                    top_cluster = max(cluster_means, key=lambda x: x[1])[0]
                    # Keep only pockets in the top cluster
                    filtered_pockets = [p for i, p in enumerate(consensus_pockets) if clusters[i] == top_cluster]
                    
                    if filtered_pockets:  # Make sure we have pockets after filtering
                        logger.info(f"Automatic filtering identified {len(filtered_pockets)} significant pockets out of {len(consensus_pockets)} total")
                        consensus_pockets = filtered_pockets
                    else:
                        logger.warning("Clustering resulted in empty top cluster, keeping all pockets")
                else:
                    logger.warning("No valid clusters found, keeping all pockets")
            except Exception as e:
                logger.warning(f"Clustering-based filtering failed: {str(e)}. Using all pockets.")
        
        # Method 2: Fallback - Use statistical outlier detection if clustering fails or for small datasets
        if len(consensus_pockets) <= 2 or len(consensus_pockets) == len(scores):
            # Skip statistical filtering if there's only one pocket
            if len(consensus_pockets) == 1:
                logger.info("Only one pocket found, keeping it without statistical filtering")
            else:
                try:
                    # Only calculate z-scores if we have more than one score
                    if len(scores) > 1:
                        # Check if all scores are nearly identical
                        score_std = np.std(scores)
                        if score_std < 1e-6:  # Very small standard deviation indicates nearly identical values
                            # If scores are nearly identical, just keep the top one
                            consensus_pockets = [consensus_pockets[0]]
                            logger.info("Scores are nearly identical. Keeping only the top-scoring pocket.")
                        else:
                            # Suppress warnings during z-score calculation
                            import warnings
                            with warnings.catch_warnings():
                                warnings.filterwarnings('ignore', category=RuntimeWarning)
                                # Calculate Z-scores of the consensus scores
                                z_scores = zscore(scores)
                            
                            # Keep pockets with z-score > 0 (above average) or top pocket if all are below average
                            significant_indices = np.where(z_scores > 0)[0]
                            
                            if len(significant_indices) > 0:
                                filtered_pockets = [consensus_pockets[i] for i in significant_indices]
                                logger.info(f"Statistical filtering identified {len(filtered_pockets)} significant pockets out of {len(consensus_pockets)} total")
                                consensus_pockets = filtered_pockets
                            else:
                                # If all scores are below average, just keep the top one
                                consensus_pockets = [consensus_pockets[0]]
                                logger.info("Keeping only the top-scoring pocket")
                    else:
                        # If we only have one score, keep it
                        logger.info("Only one pocket found, keeping it")
                except Exception as e:
                    logger.warning(f"Statistical filtering failed: {str(e)}. Using all pockets.")
    
    # Re-sort the filtered pockets
    consensus_pockets.sort(key=lambda x: (x['consensus_score'], x['final_score']), reverse=True)
    
    # Log the results
    if consensus_pockets:
        logger.info(f"Final selection: {len(consensus_pockets)} pockets with scores ranging from {consensus_pockets[0]['final_score']:.2f} to {consensus_pockets[-1]['final_score']:.2f}")
    else:
        logger.warning("No binding sites saved: None were found or all were filtered out")
        
    return consensus_pockets
