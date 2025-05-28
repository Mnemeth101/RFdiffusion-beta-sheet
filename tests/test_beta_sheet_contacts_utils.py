import unittest
import numpy as np
import torch
import random

from rfdiffusion.inference.utils import encode_beta_strand_adjacency, contig_indexed_residues_to_idx0

class TestBetaSheetContactsUtils(unittest.TestCase):
    def test_encode_beta_strand_adjacency_basic(self):
        # Set up a small adjacency matrix and parameters
        # Initialize full_adj as a 3D tensor of shape [L, L, 3] with mask values [0, 0, 1]
        # Make sure matrix is large enough for binderlen + target positions
        matrix_size = 15  # Large enough for our test case
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Set mask value to 1
        
        binderlen = 5
        # These are target-relative indices (0-indexed in the target protein)
        # They should be small enough that when converted to global indices (by adding binderlen)
        # they still fit within our matrix
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]  
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Not flexible (deterministic)
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=False
        )
        
        # Check that the correct contacts are set
        # For each residue in the beta sheet, there should be a contact between binder and target
        for i, t_idx in enumerate(target_beta_sheet_idx0):
            b_idx = binder_beta_sheet_init_position + i
            
            # Convert target index to global index (target index + binderlen)
            t_idx_global = t_idx + binderlen
            
            # Check direct matching contacts - should be [0, 1, 0] (contact, not mask)
            self.assertEqual(adj_out[t_idx_global, b_idx, 1].item(), 1.0)  # contact = [0,1,0]
            self.assertEqual(adj_out[b_idx, t_idx_global, 1].item(), 1.0)
            
            # Check adjacent contacts
            # For each position, verify contacts with adjacent positions (±1)
            for offset in [-1, 1]:  # Check adjacent positions
                target_idx = i + offset
                if 0 <= target_idx < len(target_beta_sheet_idx0):
                    adj_target_idx = target_beta_sheet_idx0[target_idx]
                    adj_target_global = adj_target_idx + binderlen
                    self.assertEqual(adj_out[b_idx, adj_target_global, 1].item(), 1.0)
                    self.assertEqual(adj_out[adj_target_global, b_idx, 1].item(), 1.0)
                    
                binder_idx = i + offset
                if 0 <= binder_idx < binder_beta_sheet_length:
                    adj_binder_idx = binder_beta_sheet_init_position + binder_idx
                    self.assertEqual(adj_out[t_idx_global, adj_binder_idx, 1].item(), 1.0)
                    self.assertEqual(adj_out[adj_binder_idx, t_idx_global, 1].item(), 1.0)

    def test_encode_beta_strand_adjacency_flexible(self):
        # Same as above but with flexible=True
        # Initialize full_adj as a 3D tensor of shape [L, L, 3] with mask values [0, 0, 1]
        # Make sure matrix is large enough for binderlen + target positions
        matrix_size = 15  # Large enough for our test case
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Set mask value to 1
        
        binderlen = 5
        # These are target-relative indices (0-indexed in the target protein)
        # They should be small enough that when converted to global indices (by adding binderlen)
        # they still fit within our matrix
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]  
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Set random seed for reproducibility since flexible mode uses randomness
        torch.manual_seed(42)
        random.seed(42)
        
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=True
        )
        
        # Should still have some contacts, but not necessarily all due to randomness
        # Convert target indices to global indices for checking contacts
        contacts = [adj_out[t_idx + binderlen, binder_beta_sheet_init_position + i, 1].item() == 1.0 for i, t_idx in enumerate(target_beta_sheet_idx0)]
        self.assertTrue(any(contacts))

    def test_contacts_exist(self):
        """Test that contacts are created between binder and target strands."""
        # Setup a standard test case
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value
        
        binderlen = 5
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Execute function
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=False
        )
        
        # Verify contacts exist
        target_global = [idx + binderlen for idx in target_beta_sheet_idx0]
        binder_positions = list(range(binder_beta_sheet_init_position, 
                               binder_beta_sheet_init_position + binder_beta_sheet_length))
        
        contact_count = sum(adj_out[t, b, 1].item() > 0 
                            for t in target_global for b in binder_positions)
        self.assertGreater(contact_count, 0, "Should create some contacts")

    def test_diagonal_contacts(self):
        """Test that matching diagonal positions have contacts."""
        # Setup
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value
        
        binderlen = 5
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Execute
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=False
        )
        
        # Verify diagonal contacts (matching positions)
        for i in range(len(target_beta_sheet_idx0)):
            t_global = target_beta_sheet_idx0[i] + binderlen
            b = binder_beta_sheet_init_position + i
            self.assertTrue(
                adj_out[t_global, b, 1].item() > 0, 
                f"Main diagonal contact missing at {t_global},{b}"
            )
    
    def test_adjacent_contacts(self):
        """Test that adjacent positions have contacts."""
        # Setup
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value
        
        binderlen = 5
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Execute
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=False
        )
        
        # Verify adjacent contacts (not just matching positions)
        # Check the middle position which should have connections to both sides
        middle_idx = len(target_beta_sheet_idx0) // 2
        t_global_middle = target_beta_sheet_idx0[middle_idx] + binderlen
        b_middle = binder_beta_sheet_init_position + middle_idx
        
        # Check if middle target residue has contact with adjacent binder residues
        b_prev = b_middle - 1
        b_next = b_middle + 1
        if b_prev >= binder_beta_sheet_init_position:
            self.assertTrue(adj_out[t_global_middle, b_prev, 1].item() > 0, 
                          "Missing contact with previous binder residue")
        if b_next < binder_beta_sheet_init_position + binder_beta_sheet_length:
            self.assertTrue(adj_out[t_global_middle, b_next, 1].item() > 0, 
                          "Missing contact with next binder residue")

        # Check if middle binder residue has contact with adjacent target residues
        t_global_prev = target_beta_sheet_idx0[middle_idx-1] + binderlen if middle_idx > 0 else None
        t_global_next = target_beta_sheet_idx0[middle_idx+1] + binderlen if middle_idx < len(target_beta_sheet_idx0)-1 else None
        
        if t_global_prev is not None:
            self.assertTrue(adj_out[b_middle, t_global_prev, 1].item() > 0, 
                          "Missing contact with previous target residue")
        if t_global_next is not None:
            self.assertTrue(adj_out[b_middle, t_global_next, 1].item() > 0, 
                          "Missing contact with next target residue")
    
    def test_contact_symmetry(self):
        """Test that contacts are symmetric between binder and target."""
        # Setup
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value
        
        binderlen = 5
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Execute
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=False
        )
        
        # Verify contacts are symmetric
        target_global = [idx + binderlen for idx in target_beta_sheet_idx0]
        binder_positions = list(range(binder_beta_sheet_init_position, 
                               binder_beta_sheet_init_position + binder_beta_sheet_length))
        
        for t in target_global:
            for b in binder_positions:
                self.assertEqual(
                    adj_out[t, b, 1].item(),
                    adj_out[b, t, 1].item(),
                    "Contacts should be symmetric"
                )

    def test_different_length_beta_sheets(self):
        """Test that the function handles different lengths of binder and target beta sheets correctly."""
        # Setup
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value

        binderlen = 5
        # Case where binder is longer than target
        target_beta_sheet_shorter = [0, 1, 2]  # Shorter target sheet
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_longer = 7  # Longer binder sheet

        # Execute
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_shorter,
            binder_beta_sheet_init_position, binder_beta_sheet_longer, flexible=False
        )

        # Verify the effective length is correctly calculated as the minimum
        effective_length = min(len(target_beta_sheet_shorter), binder_beta_sheet_longer)
        self.assertEqual(effective_length, 3)

        # Check contacts exist for positions within the effective length
        for i in range(effective_length):
            t_global = target_beta_sheet_shorter[i] + binderlen
            b = binder_beta_sheet_init_position + i
            self.assertTrue(adj_out[t_global, b, 1].item() > 0, 
                          f"Contact missing at effective position {i}")

        # Check no contacts exist beyond the effective length
        for i in range(effective_length, binder_beta_sheet_longer):
            b = binder_beta_sheet_init_position + i
            contact_beyond_effective = False
            for t in [idx + binderlen for idx in target_beta_sheet_shorter]:
                if adj_out[t, b, 1].item() > 0:
                    contact_beyond_effective = True
                    break
            self.assertFalse(contact_beyond_effective, 
                           f"Found unexpected contact beyond effective length at position {i}")

    def test_non_contiguous_target_indices(self):
        """Test that the function correctly handles non-contiguous target indices."""
        # Setup with non-contiguous target indices
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value
        
        binderlen = 5
        # Non-contiguous target indices
        target_beta_sheet_idx0 = [0, 2, 4]  # Gaps between indices
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 3
        
        # Execute
        adj_out = encode_beta_strand_adjacency(
            full_adj.clone(), binderlen, target_beta_sheet_idx0,
            binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=False
        )
        
        # Verify that correct positions have contacts
        for i, t_idx in enumerate(target_beta_sheet_idx0):
            t_global = t_idx + binderlen
            b = binder_beta_sheet_init_position + i
            self.assertTrue(adj_out[t_global, b, 1].item() > 0, 
                          f"Contact missing at position {i}")
            
    def test_flexible_mode_randomness(self):
        """Test that flexible mode introduces randomness in contacts."""
        # Setup
        matrix_size = 15
        full_adj = torch.zeros((matrix_size, matrix_size, 3), dtype=torch.float)
        full_adj[:, :, 2] = 1.0  # Mask value
        
        binderlen = 5
        target_beta_sheet_idx0 = [0, 1, 2, 3, 4]
        binder_beta_sheet_init_position = 0
        binder_beta_sheet_length = 5
        
        # Run multiple iterations with different seeds to verify randomness
        results = []
        for seed in range(5):
            torch.manual_seed(seed)
            random.seed(seed)
            
            adj_out = encode_beta_strand_adjacency(
                full_adj.clone(), binderlen, target_beta_sheet_idx0,
                binder_beta_sheet_init_position, binder_beta_sheet_length, flexible=True
            )
            
            target_global = [idx + binderlen for idx in target_beta_sheet_idx0]
            binder_positions = list(range(binder_beta_sheet_init_position, 
                                   binder_beta_sheet_init_position + binder_beta_sheet_length))
            
            # Count contacts
            contact_count = sum(adj_out[t, b, 1].item() > 0 
                                for t in target_global for b in binder_positions)
            results.append(contact_count)
        
        # At least some of the runs should have different contact counts
        # If all runs have exactly the same count, flexibility isn't working
        self.assertTrue(len(set(results)) > 1, 
                       "Flexible mode doesn't appear to introduce randomness")

if __name__ == "__main__":
    unittest.main()
