"""
Tests for the Plotly-based root locus GUI.

These tests verify the functionality of the interactive root locus plotting
using Plotly.
"""

import pytest
import numpy as np
import control as ct

# Try to import plotly, skip tests if not available
try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Try to import the GUI module
try:
    from control.interactive.rlocus_gui import root_locus_gui, rlocus_gui
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False


@pytest.mark.skipif(not PLOTLY_AVAILABLE, reason="Plotly not available")
@pytest.mark.skipif(not GUI_AVAILABLE, reason="GUI module not available")
class TestRootLocusGUI:
    """Test cases for the root locus GUI."""
    
    def setup_method(self):
        """Set up test systems."""
        s = ct.tf('s')
        self.sys1 = 1 / (s**2 + 2*s + 1)  # Simple second-order system
        self.sys2 = (s + 1) / (s**3 + 3*s**2 + 2*s)  # Third-order with zero
        self.sys3 = 1 / (s**3 + 4*s**2 + 5*s + 2)  # Third-order system
    
    def test_basic_functionality(self):
        """Test basic root locus GUI creation."""
        fig = root_locus_gui(self.sys1)
        
        assert isinstance(fig, go.Figure)
        assert len(fig.data) > 0  # Should have at least one trace
        
        # Check that the figure has the expected layout
        assert 'xaxis' in fig.layout
        assert 'yaxis' in fig.layout
        assert fig.layout.title.text == "Root Locus"
    
    def test_siso_requirement(self):
        """Test that non-SISO systems raise an error."""
        # Create a MIMO system
        s = ct.tf('s')
        mimo_sys = ct.tf([[1, 1], [0, 1]], [[s+1, 0], [0, s+2]])
        
        with pytest.raises(ValueError, match="System must be single-input single-output"):
            root_locus_gui(mimo_sys)
    
    def test_hover_info_options(self):
        """Test different hover information options."""
        hover_options = ['all', 'gain', 'damping', 'frequency']
        
        for option in hover_options:
            fig = root_locus_gui(self.sys1, hover_info=option)
            assert isinstance(fig, go.Figure)
    
    def test_grid_options(self):
        """Test grid display options."""
        # Test with grid
        fig_with_grid = root_locus_gui(self.sys1, grid=True, show_grid_lines=True)
        assert isinstance(fig_with_grid, go.Figure)
        
        # Test without grid
        fig_no_grid = root_locus_gui(self.sys1, grid=False, show_grid_lines=False)
        assert isinstance(fig_no_grid, go.Figure)
    
    def test_poles_zeros_display(self):
        """Test poles and zeros display options."""
        # Test with poles and zeros
        fig_with_pz = root_locus_gui(self.sys2, show_poles_zeros=True)
        assert isinstance(fig_with_pz, go.Figure)
        
        # Test without poles and zeros
        fig_no_pz = root_locus_gui(self.sys2, show_poles_zeros=False)
        assert isinstance(fig_no_pz, go.Figure)
    
    def test_custom_gains(self):
        """Test custom gain ranges."""
        custom_gains = np.logspace(-1, 2, 50)
        fig = root_locus_gui(self.sys1, gains=custom_gains)
        assert isinstance(fig, go.Figure)
    
    def test_custom_limits(self):
        """Test custom axis limits."""
        fig = root_locus_gui(self.sys1, xlim=(-3, 1), ylim=(-2, 2))
        assert isinstance(fig, go.Figure)
        
        # Check that limits are set correctly
        assert fig.layout.xaxis.range == [-3, 1]
        assert fig.layout.yaxis.range == [-2, 2]
    
    def test_custom_title(self):
        """Test custom title."""
        custom_title = "My Custom Root Locus"
        fig = root_locus_gui(self.sys1, title=custom_title)
        assert fig.layout.title.text == custom_title
    
    def test_custom_size(self):
        """Test custom figure size."""
        height, width = 700, 900
        fig = root_locus_gui(self.sys1, height=height, width=width)
        assert fig.layout.height == height
        assert fig.layout.width == width
    
    def test_convenience_function(self):
        """Test the convenience function rlocus_gui."""
        fig = rlocus_gui(self.sys1)
        assert isinstance(fig, go.Figure)
    
    def test_complex_system(self):
        """Test with a more complex system."""
        s = ct.tf('s')
        complex_sys = (s**2 + 2*s + 2) / (s**4 + 5*s**3 + 8*s**2 + 6*s + 2)
        
        fig = root_locus_gui(complex_sys, hover_info='all')
        assert isinstance(fig, go.Figure)
    
    def test_damping_frequency_lines(self):
        """Test damping and frequency line options."""
        # Test damping lines only
        fig_damping = root_locus_gui(self.sys1, damping_lines=True, frequency_lines=False)
        assert isinstance(fig_damping, go.Figure)
        
        # Test frequency lines only
        fig_freq = root_locus_gui(self.sys1, damping_lines=False, frequency_lines=True)
        assert isinstance(fig_freq, go.Figure)
        
        # Test both
        fig_both = root_locus_gui(self.sys1, damping_lines=True, frequency_lines=True)
        assert isinstance(fig_both, go.Figure)
    
    def test_data_consistency(self):
        """Test that the GUI data is consistent with the original root locus."""
        # Get data from the GUI
        fig = root_locus_gui(self.sys1)
        
        # Get data from the original root locus function
        rl_data = ct.root_locus_map(self.sys1)
        
        # Check that we have the same number of loci
        if rl_data.loci is not None:
            num_loci = rl_data.loci.shape[1]
            # The GUI should have traces for the loci plus poles/zeros
            assert len(fig.data) >= num_loci


@pytest.mark.skipif(not PLOTLY_AVAILABLE, reason="Plotly not available")
@pytest.mark.skipif(not GUI_AVAILABLE, reason="GUI module not available")
def test_import_functionality():
    """Test that the GUI module can be imported and used."""
    from control.interactive.rlocus_gui import root_locus_gui, rlocus_gui
    
    # Create a simple system
    s = ct.tf('s')
    sys = 1 / (s**2 + 2*s + 1)
    
    # Test both functions
    fig1 = root_locus_gui(sys)
    fig2 = rlocus_gui(sys)
    
    assert isinstance(fig1, go.Figure)
    assert isinstance(fig2, go.Figure)


if __name__ == "__main__":
    # Run a simple test if executed directly
    if PLOTLY_AVAILABLE and GUI_AVAILABLE:
        s = ct.tf('s')
        sys = 1 / (s**2 + 2*s + 1)
        fig = root_locus_gui(sys, title="Test Plot")
        print("Test successful! Created root locus GUI.")
        # Uncomment the next line to show the plot
        # fig.show()
    else:
        print("Plotly or GUI module not available for testing.") 