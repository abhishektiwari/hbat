"""GUI components for HBAT."""


def main():
    """Launch the HBAT GUI application."""
    import sys

    try:
        from .main_window import MainWindow

        app = MainWindow()
        app.run()

    except ImportError as e:
        print(f"Error: Missing required dependencies: {e}")
        print("Please ensure all required packages are installed.")
        sys.exit(1)
    except Exception as e:
        print(f"Error starting HBAT GUI: {e}")
        sys.exit(1)
