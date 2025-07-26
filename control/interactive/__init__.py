"""
Interactive plotting tools for the Python Control Systems Library.

This module provides interactive plotting capabilities using matplotlib,
including root locus analysis with hover functionality.
"""

# Import matplotlib-based functions
try:
    from .rlocus_gui import root_locus_gui, rlocus_gui, root_locus_gui_advanced
    __all__ = ['root_locus_gui', 'rlocus_gui', 'root_locus_gui_advanced']
except ImportError as e:
    # If there's an import error, provide informative error messages
    def root_locus_gui(*args, **kwargs):
        raise ImportError(
            f"root_locus_gui could not be imported: {e}. "
            "Make sure matplotlib is installed: pip install matplotlib"
        )
    
    def rlocus_gui(*args, **kwargs):
        raise ImportError(
            f"rlocus_gui could not be imported: {e}. "
            "Make sure matplotlib is installed: pip install matplotlib"
        )
    
    def root_locus_gui_advanced(*args, **kwargs):
        raise ImportError(
            f"root_locus_gui_advanced could not be imported: {e}. "
            "Make sure matplotlib is installed: pip install matplotlib"
        )
    
    __all__ = ['root_locus_gui', 'rlocus_gui', 'root_locus_gui_advanced'] 