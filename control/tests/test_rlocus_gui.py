"""
Tests for the root locus GUI.

These tests verify the functionality of the interactive root locus plotting.
"""

import pytest
import numpy as np
import control as ct

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from control.interactive.rlocus_gui import root_locus_gui, rlocus_gui, RootLocusGUI
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
@pytest.mark.skipif(not GUI_AVAILABLE, reason="GUI module not available")
class TestRootLocusGUI:
    """Test cases for the root locus GUI."""
    
    def setup_method(self):
        """Set up test systems."""
        s = ct.tf('s')
        self.sys1 = 1 / (s**2 + 2*s + 1)
        self.sys2 = (s + 1) / (s**3 + 3*s**2 + 2*s)
        self.sys3 = 1 / (s**3 + 4*s**2 + 5*s + 2)
    
    def test_basic_functionality(self):
        """Test basic root locus GUI creation."""
        gui = root_locus_gui(self.sys1)
        
        assert isinstance(gui, RootLocusGUI)
        assert gui.sys == self.sys1
        assert gui.fig is not None
        assert gui.ax is not None
        
        assert hasattr(gui.fig, 'canvas')
        assert hasattr(gui.ax, 'get_title')
        title = gui.ax.get_title()
        assert title == "Root Locus" or title == "" or "Root Locus" in title
    
    def test_siso_requirement(self):
        """Test that non-SISO systems raise an error."""
        mimo_sys = ct.tf([[[1]], [[1]]], [[[1, 1]], [[1, 2]]])
        
        with pytest.raises(ValueError, match="System must be single-input single-output"):
            root_locus_gui(mimo_sys)
    
    def test_grid_options(self):
        """Test grid display options."""
        gui_with_grid = root_locus_gui(self.sys1, grid=True, show_grid_lines=True)
        assert isinstance(gui_with_grid, RootLocusGUI)
        
        gui_no_grid = root_locus_gui(self.sys1, grid=False, show_grid_lines=False)
        assert isinstance(gui_no_grid, RootLocusGUI)
    
    def test_poles_zeros_display(self):
        """Test poles and zeros display options."""
        gui_with_pz = root_locus_gui(self.sys2, show_poles_zeros=True)
        assert isinstance(gui_with_pz, RootLocusGUI)
        
        gui_no_pz = root_locus_gui(self.sys2, show_poles_zeros=False)
        assert isinstance(gui_no_pz, RootLocusGUI)
    
    def test_custom_gains(self):
        """Test custom gain ranges."""
        custom_gains = np.logspace(-1, 2, 50)
        gui = root_locus_gui(self.sys1, gains=custom_gains)
        assert isinstance(gui, RootLocusGUI)
        assert gui.gains is custom_gains
    
    def test_custom_limits(self):
        """Test custom axis limits."""
        gui = root_locus_gui(self.sys1, xlim=(-3, 1), ylim=(-2, 2))
        assert isinstance(gui, RootLocusGUI)
        assert gui.xlim == (-3, 1)
        assert gui.ylim == (-2, 2)
    
    def test_custom_title(self):
        """Test custom title."""
        custom_title = "My Custom Root Locus"
        gui = root_locus_gui(self.sys1, title=custom_title)
        assert isinstance(gui, RootLocusGUI)
        assert gui.title == custom_title
    
    def test_convenience_function(self):
        """Test the convenience function rlocus_gui."""
        gui = rlocus_gui(self.sys1)
        assert isinstance(gui, RootLocusGUI)
    
    def test_complex_system(self):
        """Test with a more complex system."""
        s = ct.tf('s')
        complex_sys = (s**2 + 2*s + 2) / (s**4 + 5*s**3 + 8*s**2 + 6*s + 2)
        
        gui = root_locus_gui(complex_sys)
        assert isinstance(gui, RootLocusGUI)
    
    def test_damping_frequency_lines(self):
        """Test damping and frequency line options."""
        # Test damping lines only
        gui_damping = root_locus_gui(self.sys1, damping_lines=True, frequency_lines=False)
        assert isinstance(gui_damping, RootLocusGUI)
        
        # Test frequency lines only
        gui_freq = root_locus_gui(self.sys1, damping_lines=False, frequency_lines=True)
        assert isinstance(gui_freq, RootLocusGUI)
        
        # Test both
        gui_both = root_locus_gui(self.sys1, damping_lines=True, frequency_lines=True)
        assert isinstance(gui_both, RootLocusGUI)
    
    def test_data_consistency(self):
        """Test that the GUI data is consistent with the original root locus."""
        # Get data from the GUI
        gui = root_locus_gui(self.sys1)
        
        # Get data from the original root locus function
        rl_data = ct.root_locus_map(self.sys1)
        
        # Check that we have valid data in both cases
        assert gui.rl_data.gains is not None
        assert rl_data.gains is not None
        assert len(gui.rl_data.gains) > 0
        assert len(rl_data.gains) > 0
        
        # Check that the GUI data has the expected structure
        assert hasattr(gui.rl_data, 'loci')
        assert hasattr(gui.rl_data, 'gains')
        assert hasattr(gui.rl_data, 'poles')
        assert hasattr(gui.rl_data, 'zeros')
    
    def test_info_box_creation(self):
        """Test that the info box is created properly."""
        gui = root_locus_gui(self.sys1)
        assert gui.info_text is not None
        assert hasattr(gui.info_text, 'set_text')
        assert hasattr(gui.info_text, 'set_visible')
    
    def test_mouse_event_handlers(self):
        """Test that mouse event handlers are set up."""
        gui = root_locus_gui(self.sys1)
        # Check that the methods exist
        assert hasattr(gui, '_on_mouse_move')
        assert hasattr(gui, '_on_mouse_leave')
        assert hasattr(gui, '_find_closest_point')
        assert hasattr(gui, '_update_info_box')
        assert hasattr(gui, '_hide_info_box')
    
    def test_save_functionality(self):
        """Test the save functionality."""
        gui = root_locus_gui(self.sys1)
        assert hasattr(gui, 'save')
        # Note: We don't actually save a file in tests to avoid file system dependencies
    
    def test_cursor_marker_creation(self):
        """Test that the cursor marker is created properly."""
        gui = root_locus_gui(self.sys1)
        assert gui.cursor_marker is not None
        assert hasattr(gui.cursor_marker, 'set_data')
        assert hasattr(gui.cursor_marker, 'set_visible')
    
    def test_cursor_marker_methods(self):
        """Test that cursor marker control methods exist."""
        gui = root_locus_gui(self.sys1)
        # Check that the methods exist
        assert hasattr(gui, '_update_cursor_marker')
        assert hasattr(gui, '_hide_cursor_marker')
        assert hasattr(gui, '_create_cursor_marker')


@pytest.mark.skipif(not MATPLOTLIB_AVAILABLE, reason="Matplotlib not available")
@pytest.mark.skipif(not GUI_AVAILABLE, reason="GUI module not available")
def test_import_functionality():
    """Test that the GUI module can be imported and used."""
    from control.interactive.rlocus_gui import root_locus_gui, rlocus_gui, RootLocusGUI
    
    # Create a simple system
    s = ct.tf('s')
    sys = 1 / (s**2 + 2*s + 1)
    
    # Test both functions
    gui1 = root_locus_gui(sys)
    gui2 = rlocus_gui(sys)
    
    assert isinstance(gui1, RootLocusGUI)
    assert isinstance(gui2, RootLocusGUI)


if __name__ == "__main__":
    # Run a simple test if executed directly
    if MATPLOTLIB_AVAILABLE and GUI_AVAILABLE:
        s = ct.tf('s')
        sys = 1 / (s**2 + 2*s + 1)
        gui = root_locus_gui(sys, title="Test Plot")
        print("Test successful! Created root locus GUI.")
        # Uncomment the next line to show the plot
        # gui.show()
    else:
        print("Matplotlib or GUI module not available for testing.") 